/*!
 * \file CSysMatrixGPU.cu
 * \brief Implementations of Kernels and Functions for Matrix Operations on the GPU
 * \author A. Raj, Jesse Li, P. Gomes
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>

#include "../../include/linear_algebra/CMatrixInverse.hpp"
#include "../../include/linear_algebra/CSysMatrix.inl"
#include "../../include/linear_algebra/GPUComms.cuh"

namespace {

template <class ScalarType>
__global__ void ApplyJacobiPreconditionerKernel(const ScalarType* invM, const ScalarType* vec, ScalarType* prod,
                                                unsigned long nPointDomain, unsigned long nVar) {
  const auto iPoint = static_cast<unsigned long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (iPoint >= nPointDomain) return;

  const auto block = &invM[iPoint * nVar * nVar];
  const auto rhs = &vec[iPoint * nVar];
  auto out = &prod[iPoint * nVar];

  for (auto iVar = 0ul; iVar < nVar; ++iVar) {
    auto sum = ScalarType(0);
    for (auto jVar = 0ul; jVar < nVar; ++jVar) {
      sum += block[iVar * nVar + jVar] * rhs[jVar];
    }
    out[iVar] = sum;
  }
}

/*--- ILU. The factorization and both triangular solves are level scheduled: the rows of a
 * level are independent of each other and only depend on rows of previous levels, so each
 * level is one kernel launch and the launch boundaries provide the synchronization. The rows
 * of a level are scattered through the matrix, hence the indirection through the level table.
 * Throughout, one CUDA block works on one row. ---*/

/*!
 * \brief The pointers of an LDU-partitioned matrix, all in device memory. This mirrors the
 *        private CSysMatrix::LDU, which the kernels cannot name.
 */
template <class ScalarType>
struct DeviceLDU {
  ScalarType* d;
  ScalarType* l;
  ScalarType* u;
  const su2uint* row_ptr_l;
  const su2uint* col_ind_l;
  const su2uint* row_ptr_u;
  const su2uint* col_ind_u;
};

/*!
 * \brief Start of block (i,j), or nullptr if it is not a nonzero of the pattern.
 */
template <class ScalarType>
__device__ FORCEINLINE ScalarType* GetBlockILU(const DeviceLDU<ScalarType>& M, unsigned long nVar,
                                               unsigned long block_i, unsigned long block_j) {
  const auto blockSize = nVar * nVar;
  if (block_i == block_j) return M.d + block_i * blockSize;

  const bool lower = block_j < block_i;
  const auto* row_ptr = lower ? M.row_ptr_l : M.row_ptr_u;
  const auto* col_ind = lower ? M.col_ind_l : M.col_ind_u;
  auto* vals = lower ? M.l : M.u;

  for (auto k = row_ptr[block_i]; k < row_ptr[block_i + 1]; ++k) {
    if (col_ind[k] == block_j) return vals + k * blockSize;
  }
  return nullptr;
}

/*!
 * \brief Invert the diagonal blocks of the matrix, device version of InverseDiagonalBlock.
 * \note Grid: one block per row, blockDim.x == nVar*nVar (one thread per block entry, they
 *       stage the input in shared memory). Dynamic shared memory: nVar*nVar scalars.
 */
template <class ScalarType>
__global__ void InvertDiagonalBlocksKernel(unsigned long nRows, unsigned long nVar,
                                           const ScalarType* __restrict__ mat_d, ScalarType* __restrict__ invM) {
  const unsigned long iRow = blockIdx.x;
  if (iRow >= nRows) return;

  const auto blockSize = nVar * nVar;
  const unsigned long tid = threadIdx.x;

  /*--- The inversion destroys its input, so it cannot work on the matrix itself. ---*/
  extern __shared__ __align__(sizeof(double)) char smem[];
  auto* work = reinterpret_cast<ScalarType*>(smem);

  work[tid] = mat_d[iRow * blockSize + tid];
  __syncthreads();

  if (tid == 0) SU2_LinAlg::MatrixInverse(nVar, work, invM + iRow * blockSize);
}

/*!
 * \brief Factorize the rows of one color, one sweep of an iterative (colored Gauss-Seidel)
 *        ILU factorization: same order/pattern as the exact level-scheduled algorithm (this
 *        does not change L/U membership, so it converges to the exact same fixed point), but
 *        colors are true independent sets (zero dependency between same-colored rows in either
 *        direction), so a color can be processed with zero races in far fewer, wider launches
 *        than the number of levels — at the cost of needing several sweeps (repeated passes
 *        over all colors) instead of one exact pass, because for a fixed order the level count
 *        is already the minimum number of race-free single-pass groups (Mirsky's theorem).
 * \note Every visit of a row (there is one per sweep) resets it from the original matrix first
 *       (folding in the device version of the InitIluRow helper of BuildILUPreconditioner),
 *       because the elimination below is a re-evaluation of the row's defining equation using
 *       the current (possibly stale) values of other rows, not an incremental accumulation.
 *       Grid: one block per row of the color, blockDim.x == nVar*nVar (one thread per block
 *       entry, so that the small matrix products are one dot product per thread). Dynamic
 *       shared memory: 2*nVar*nVar scalars.
 */
template <class ScalarType>
__global__ void IluFactorColorKernel(const su2uint* __restrict__ color_idx, unsigned long color_begin,
                                     unsigned long color_size, unsigned long nRows, unsigned long nVar,
                                     DeviceLDU<ScalarType> A, DeviceLDU<ScalarType> M) {
  if (blockIdx.x >= color_size) return;

  const unsigned long iRow = color_idx[color_begin + blockIdx.x];
  const auto blockSize = nVar * nVar;
  const unsigned long tid = threadIdx.x;
  const auto iVar = tid / nVar, jVar = tid % nVar;

  extern __shared__ __align__(sizeof(double)) char smem[];
  auto* Lij = reinterpret_cast<ScalarType*>(smem);
  auto* work = Lij + blockSize;

  /*--- Reset this row to the raw matrix entries (device version of InitIluRow, but for one
   * row instead of the whole matrix, since here it runs once per row per sweep). ---*/
  M.d[iRow * blockSize + tid] = A.d[iRow * blockSize + tid];
  auto scatter = [&](const su2uint* a_row_ptr, const su2uint* a_col_ind, const ScalarType* a_vals,
                     const su2uint* m_row_ptr, const su2uint* m_col_ind, ScalarType* m_vals) {
    auto ka = a_row_ptr[iRow];
    const auto ka_end = a_row_ptr[iRow + 1];
    for (auto k = m_row_ptr[iRow]; k < m_row_ptr[iRow + 1]; ++k) {
      const auto jPoint = m_col_ind[k];
      while (ka < ka_end && a_col_ind[ka] < jPoint) ++ka;
      if (ka < ka_end && a_col_ind[ka] == jPoint) {
        m_vals[k * blockSize + tid] = a_vals[ka * blockSize + tid];
      } else {
        m_vals[k * blockSize + tid] = ScalarType(0);
      }
    }
  };
  scatter(A.row_ptr_l, A.col_ind_l, A.l, M.row_ptr_l, M.col_ind_l, M.l);
  scatter(A.row_ptr_u, A.col_ind_u, A.u, M.row_ptr_u, M.col_ind_u, M.u);
  __syncthreads();

  /*--- For this row (unknown), loop over its lower diagonal entries. ---*/
  for (auto kl = M.row_ptr_l[iRow]; kl < M.row_ptr_l[iRow + 1]; ++kl) {
    /*--- All threads must be done with the previous entry: Lij is about to be overwritten,
     * and the blocks of this row updated below are read here across threads. ---*/
    __syncthreads();

    /*--- jPoint is the column index (jPoint < iRow). ---*/
    const unsigned long jPoint = M.col_ind_l[kl];

    /*--- Multiply the block by the inverse of the corresponding diagonal block. ---*/
    auto* Block_ij = M.l + kl * blockSize;
    const auto* invUjj = M.d + jPoint * blockSize;

    ScalarType sum = 0;
    for (auto k = 0ul; k < nVar; ++k) sum += Block_ij[iVar * nVar + k] * invUjj[k * nVar + jVar];
    Lij[tid] = sum;
    __syncthreads();

    /*--- Lij holds Aij*inv(Ujj). Jump to the upper part of the jPoint row. ---*/
    for (auto ku = M.row_ptr_u[jPoint]; ku < M.row_ptr_u[jPoint + 1]; ++ku) {
      /*--- Get the column index (kPoint > jPoint), halo columns are not factorized. ---*/
      const unsigned long kPoint = M.col_ind_u[ku];
      if (kPoint >= nRows) break;

      /*--- If Aik exists, update it: Aik -= Lij * Ujk ---*/
      auto* Block_ik = GetBlockILU(M, nVar, iRow, kPoint);
      if (Block_ik == nullptr) continue;

      /*--- Block_ik cannot alias Block_ij because kPoint > jPoint. ---*/
      const auto* Ujk = M.u + ku * blockSize;
      ScalarType prod = 0;
      for (auto k = 0ul; k < nVar; ++k) prod += Lij[iVar * nVar + k] * Ujk[k * nVar + jVar];
      Block_ik[tid] -= prod;
    }

    /*--- Store Lij in the lower triangular part, each thread only writes its own entry. ---*/
    Block_ij[tid] = Lij[tid];
  }

  /*--- Invert the diagonal entry, Uii, for the rows that depend on it. The loop above may have
   * updated it (when kPoint == iRow), so the whole block has to be done first. ---*/
  __syncthreads();
  work[tid] = M.d[iRow * blockSize + tid];
  __syncthreads();
  if (tid == 0) SU2_LinAlg::MatrixInverse(nVar, work, M.d + iRow * blockSize);
}

/*!
 * \brief Forward substitution for the rows of one level, (L+I).prod = vec.
 * \note One thread per block entry (like the factorization kernel), so the inner dot product
 *       over a neighbor block is spread across nVar threads instead of done serially by one.
 *       Unlike the factorization, forward/backward only ever *read* already-finalized prod
 *       values (no row-to-row write-through during the loop), so each thread can accumulate
 *       its own (iVar,jVar) partial product across every neighbor with no synchronization at
 *       all, and only the final nVar-way reduction (summing over jVar for each iVar) needs one
 *       __syncthreads() — not one per neighbor. Grid: one block per row of the level,
 *       blockDim.x == nVar*nVar. Dynamic shared memory: nVar*nVar scalars.
 */
template <class ScalarType>
__global__ void IluForwardLevelKernel(const su2uint* __restrict__ level_idx, unsigned long level_begin,
                                      unsigned long level_size, unsigned long nVar, DeviceLDU<ScalarType> M,
                                      const ScalarType* __restrict__ vec, ScalarType* __restrict__ prod) {
  if (blockIdx.x >= level_size) return;

  const unsigned long iRow = level_idx[level_begin + blockIdx.x];
  const unsigned long tid = threadIdx.x;
  const auto iVar = tid / nVar, jVar = tid % nVar;

  extern __shared__ __align__(sizeof(double)) char smem[];
  auto* partial = reinterpret_cast<ScalarType*>(smem);

  ScalarType acc = 0;
  /*--- The columns of L are rows of previous levels, so prod is final for all of them. ---*/
  for (auto kl = M.row_ptr_l[iRow]; kl < M.row_ptr_l[iRow + 1]; ++kl) {
    const unsigned long jPoint = M.col_ind_l[kl];
    const auto* blk = M.l + kl * nVar * nVar;
    acc += blk[iVar * nVar + jVar] * prod[jPoint * nVar + jVar];
  }
  partial[tid] = acc;
  __syncthreads();

  if (jVar == 0) {
    ScalarType sum = vec[iRow * nVar + iVar];
    for (auto j = 0ul; j < nVar; ++j) sum -= partial[iVar * nVar + j];
    prod[iRow * nVar + iVar] = sum;
  }
}

/*!
 * \brief Backward substitution for the rows of one level, U.prod = prod.
 * \note Same idea as IluForwardLevelKernel: each thread accumulates its own (iVar,jVar)
 *       partial product across every neighbor with no synchronization, then one
 *       __syncthreads() to reduce over jVar and get the elimination result per iVar, then a
 *       second __syncthreads() before the diagonal multiply (which needs every iVar's result).
 *       Grid: one block per row of the level, blockDim.x == nVar*nVar. Dynamic shared memory:
 *       nVar*nVar + nVar scalars.
 */
template <class ScalarType>
__global__ void IluBackwardLevelKernel(const su2uint* __restrict__ level_idx, unsigned long level_begin,
                                       unsigned long level_size, unsigned long nRows, unsigned long nVar,
                                       DeviceLDU<ScalarType> M, ScalarType* __restrict__ prod) {
  if (blockIdx.x >= level_size) return;

  const unsigned long iRow = level_idx[level_begin + blockIdx.x];
  const auto blockSize = nVar * nVar;
  const unsigned long tid = threadIdx.x;
  const auto iVar = tid / nVar, jVar = tid % nVar;

  extern __shared__ __align__(sizeof(double)) char smem[];
  auto* partial = reinterpret_cast<ScalarType*>(smem);
  auto* aux = partial + blockSize;

  ScalarType acc = 0;
  /*--- The columns of U are rows of later levels, already updated by this sweep. ---*/
  for (auto ku = M.row_ptr_u[iRow]; ku < M.row_ptr_u[iRow + 1]; ++ku) {
    const unsigned long jPoint = M.col_ind_u[ku];
    if (jPoint >= nRows) break;
    const auto* blk = M.u + ku * blockSize;
    acc += blk[iVar * nVar + jVar] * prod[jPoint * nVar + jVar];
  }
  partial[tid] = acc;
  __syncthreads();

  if (jVar == 0) {
    ScalarType sum = prod[iRow * nVar + iVar];
    for (auto j = 0ul; j < nVar; ++j) sum -= partial[iVar * nVar + j];
    /*--- The diagonal blocks are stored inverted by the factorization. ---*/
    aux[iVar] = sum;
  }
  __syncthreads();

  if (jVar == 0) {
    const auto* invUii = M.d + iRow * blockSize;
    ScalarType out = 0;
    for (auto k = 0ul; k < nVar; ++k) out += invUii[iVar * nVar + k] * aux[k];
    prod[iRow * nVar + iVar] = out;
  }
}

/*!
 * \brief Block-LDU SpMV kernel: y[iRow] = (L + D + U) * x per block-row.
 *        One CUDA block per block-row; threadIdx.x indexes output variable (0..nVar-1).
 */
template <class ScalarType>
__global__ void BlockLDU_SpMV_kernel(unsigned long nRows, unsigned long nVar,
                                     const su2uint* __restrict__ row_ptr_l,
                                     const su2uint* __restrict__ col_ind_l,
                                     const ScalarType* __restrict__ mat_l,
                                     const ScalarType* __restrict__ mat_d,
                                     const su2uint* __restrict__ row_ptr_u,
                                     const su2uint* __restrict__ col_ind_u,
                                     const ScalarType* __restrict__ mat_u,
                                     const ScalarType* __restrict__ x, ScalarType* __restrict__ y) {
  const unsigned long iRow = blockIdx.x;
  const unsigned long iVar = threadIdx.x;
  if (iRow >= nRows || iVar >= nVar) return;

  ScalarType sum = 0;
  /* Lower */
  for (auto k = row_ptr_l[iRow]; k < row_ptr_l[iRow + 1]; ++k) {
    const auto col = col_ind_l[k];
    const ScalarType* blk = mat_l + k * nVar * nVar + iVar * nVar;
    for (unsigned long jVar = 0; jVar < nVar; ++jVar) sum += blk[jVar] * x[col * nVar + jVar];
  }
  /* Diagonal */
  {
    const ScalarType* blk = mat_d + iRow * nVar * nVar + iVar * nVar;
    for (unsigned long jVar = 0; jVar < nVar; ++jVar) sum += blk[jVar] * x[iRow * nVar + jVar];
  }
  /* Upper */
  for (auto k = row_ptr_u[iRow]; k < row_ptr_u[iRow + 1]; ++k) {
    const auto col = col_ind_u[k];
    const ScalarType* blk = mat_u + k * nVar * nVar + iVar * nVar;
    for (unsigned long jVar = 0; jVar < nVar; ++jVar) sum += blk[jVar] * x[col * nVar + jVar];
  }
  y[iRow * nVar + iVar] = sum;
}

}  // namespace

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeJacobiPreconditionerGPU(const CSysVector<ScalarType>& vec,
                                                            CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                            const CConfig* config) const {
  (void)geometry;
  (void)config;

  SU2_ZONE_SCOPED

  if (d_invM == nullptr) {
    SU2_MPI::Error("CUDA Jacobi preconditioner used before BuildJacobiPreconditionerGPU.", CURRENT_FUNCTION);
  }

  constexpr unsigned threadsPerBlock = 128;
  const auto blocks = static_cast<unsigned>((nPointDomain + threadsPerBlock - 1) / threadsPerBlock);
  ApplyJacobiPreconditionerKernel<<<blocks, threadsPerBlock>>>(d_invM, vec.GetDevicePointer(), prod.GetDevicePointer(),
                                                               nPointDomain, nVar);
  /*--- Sync so the zone above actually times the kernel, not just the (async) launch call. ---*/
  gpuErrChk(cudaStreamSynchronize(nullptr));
  gpuErrChk(cudaGetLastError());
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildJacobiPreconditionerGPU() {
  SU2_ZONE_SCOPED

  if (d_invM == nullptr) {
    SU2_MPI::Error("CUDA Jacobi preconditioner used without device storage.", CURRENT_FUNCTION);
  }
  if (nVar != nEqn) {
    SU2_MPI::Error("CUDA Jacobi preconditioner requires square blocks.", CURRENT_FUNCTION);
  }
  if (nVar * nVar > 1024) {
    SU2_MPI::Error("CUDA Jacobi preconditioner uses one thread per block entry, nVar is too large.", CURRENT_FUNCTION);
  }
  if (nPointDomain == 0) return;

  /*--- The matrix is expected to be on the device already, it is uploaded once per solve by
   * CSysMatrixVectorProduct, which is created before the preconditioner is built. ---*/
  const auto blockSize = static_cast<unsigned>(nVar * nVar);
  InvertDiagonalBlocksKernel<ScalarType>
      <<<static_cast<unsigned>(nPointDomain), blockSize, blockSize * sizeof(ScalarType)>>>(nPointDomain, nVar, gpu.d,
                                                                                           d_invM);
  /*--- Sync so the zone above actually times the kernel, not just the (async) launch call. ---*/
  gpuErrChk(cudaStreamSynchronize(nullptr));
  gpuErrChk(cudaGetLastError());
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildILUPreconditionerGPU() {
  SU2_ZONE_SCOPED

  if (gpu_ilu.d == nullptr) {
    SU2_MPI::Error("CUDA ILU preconditioner used without device storage.", CURRENT_FUNCTION);
  }
  if (nVar != nEqn) {
    SU2_MPI::Error("CUDA ILU factorization requires square blocks.", CURRENT_FUNCTION);
  }
  if (nVar * nVar > 1024) {
    SU2_MPI::Error("CUDA ILU factorization uses one thread per block entry, nVar is too large.", CURRENT_FUNCTION);
  }
  if (nPointDomain == 0) return;

  /*--- The matrix is expected to be on the device already, it is uploaded once per solve by
   * CSysMatrixVectorProduct, which is created before the preconditioner is built. ---*/
  const DeviceLDU<ScalarType> A{gpu.d,         gpu.l,         gpu.u,        gpu.row_ptr_l,
                                gpu.col_ind_l, gpu.row_ptr_u, gpu.col_ind_u};
  const DeviceLDU<ScalarType> M{gpu_ilu.d,         gpu_ilu.l,         gpu_ilu.u,        gpu_ilu.row_ptr_l,
                                gpu_ilu.col_ind_l, gpu_ilu.row_ptr_u, gpu_ilu.col_ind_u};

  const auto blockSize = static_cast<unsigned>(nVar * nVar);
  const auto shared = 2 * blockSize * sizeof(ScalarType);

  /*--- The legacy default stream cannot be captured, so the graph lives on its own stream,
   * created once. Every launch below is followed by a sync back to the host, so this does not
   * change execution order relative to the rest of the (single-stream) solver. ---*/
  if (ilu_stream == nullptr) gpuErrChk(cudaStreamCreate(&ilu_stream));

  /*--- The launch sequence (ilu_gpu_color_sweeps passes over all colors) is identical on every
   * call: the grid and block sizes only depend on the (fixed) sparsity pattern/coloring and the
   * device pointers are fixed members, allocated once. Capture it into a CUDA graph the first
   * time and replay that from then on, which removes the per-launch host-side overhead without
   * touching the parallelization of any individual kernel (unlike a persistent cooperative-
   * groups kernel, this does not cap per-color parallelism to an occupancy-resident block
   * count). See IluFactorColorKernel for why several sweeps over the colors are needed.
   * Note that factors are not reset between calls to BuildILUPreconditionerGPU, so with
   * LINEAR_SOLVER_ILU_GPU_SWEEPS set low (even 1), each call refines the previous one's
   * result rather than reconverging from scratch, relying on the matrix changing little
   * between outer/pseudo-time iterations. ---*/
  if (ilu_build_graph_exec == nullptr) {
    cudaGraph_t graph;
    gpuErrChk(cudaStreamBeginCapture(ilu_stream, cudaStreamCaptureModeThreadLocal));

    for (unsigned short sweep = 0; sweep < ilu_gpu_color_sweeps; ++sweep) {
      for (auto color = 0ul; color + 1 < ilu_color_ptr.size(); ++color) {
        const auto begin = ilu_color_ptr[color];
        const auto size = ilu_color_ptr[color + 1] - begin;
        if (size == 0) continue;
        IluFactorColorKernel<ScalarType>
            <<<size, blockSize, shared, ilu_stream>>>(d_ilu_color_idx, begin, size, nPointDomain, nVar, A, M);
      }
    }

    gpuErrChk(cudaStreamEndCapture(ilu_stream, &graph));
    gpuErrChk(cudaGraphInstantiate(&ilu_build_graph_exec, graph, nullptr, nullptr, 0));
    gpuErrChk(cudaGraphDestroy(graph));
  }

  gpuErrChk(cudaGraphLaunch(ilu_build_graph_exec, ilu_stream));
  gpuErrChk(cudaStreamSynchronize(ilu_stream));
  gpuErrChk(cudaGetLastError());
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeILUPreconditionerGPU(const CSysVector<ScalarType>& vec,
                                                         CSysVector<ScalarType>& prod) const {
  SU2_ZONE_SCOPED

  if (gpu_ilu.d == nullptr) {
    SU2_MPI::Error("CUDA ILU preconditioner used before BuildILUPreconditionerGPU.", CURRENT_FUNCTION);
  }
  if (nPointDomain == 0) return;

  const DeviceLDU<ScalarType> M{gpu_ilu.d,         gpu_ilu.l,         gpu_ilu.u,        gpu_ilu.row_ptr_l,
                                gpu_ilu.col_ind_l, gpu_ilu.row_ptr_u, gpu_ilu.col_ind_u};

  auto* d_vec = vec.GetDevicePointer();
  auto* d_prod = prod.GetDevicePointer();

  const auto nLevels = ilu_level_ptr.size() - 1;

  /*--- One thread per block entry, like the factorization kernel: spreads each row's neighbor
   * dot products over nVar*nVar threads instead of doing them serially in nVar threads, without
   * changing the number of blocks (still one per row), so this does not trade away SM coverage
   * the way batching several rows into a block did. ---*/
  const auto threads = static_cast<unsigned>(nVar * nVar);
  const auto sharedForward = threads * sizeof(ScalarType);
  const auto sharedBackward = (threads + nVar) * sizeof(ScalarType);

  if (ilu_stream == nullptr) gpuErrChk(cudaStreamCreate(&ilu_stream));

  /*--- Same idea as BuildILUPreconditionerGPU: the launch sequence only depends on the (fixed)
   * level structure, plus the vec/prod device pointers. Those normally are the same temporary
   * buffers on every call (owned by CSysSolve / CSysVector, allocated once), so the graph is
   * captured once and replayed; if the pointers ever do change the graph is recaptured, which
   * is no worse than the un-graphed loop, just not free. ---*/
  if (ilu_apply_graph_exec == nullptr || ilu_apply_graph_vec != d_vec || ilu_apply_graph_prod != d_prod) {
    if (ilu_apply_graph_exec != nullptr) {
      gpuErrChk(cudaGraphExecDestroy(ilu_apply_graph_exec));
      ilu_apply_graph_exec = nullptr;
    }

    cudaGraph_t graph;
    gpuErrChk(cudaStreamBeginCapture(ilu_stream, cudaStreamCaptureModeThreadLocal));

    /*--- Forward substitution: colored-iterative, not level-scheduled. IluForwardLevelKernel
     * only ever reads already-finalized prod values within a single exact level-scheduled pass;
     * across colors that guarantee is gone (a color is wider/less ordered than a level), so an
     * L-neighbor in a not-yet-processed color this sweep is read at whatever value it currently
     * holds. Repeated full passes over all colors still converge to the same forward-solve
     * result (validated on the host: same fixed point, ~4.7x error reduction per sweep, colors
     * are true independent sets so this is race-free). Zero prod first: on the very first call
     * this buffer may hold unrelated leftover data, and reading it would otherwise make the
     * first sweep ill-defined. Unlike the level-scheduled solve this replaces, the sweep count is
     * paid on every apply call, not amortized — see ilu_gpu_fwd_sweeps.
     * The backward solve stays exact/level-scheduled below: a host experiment showed the
     * colored-iterative version diverges to NaN for this class of matrix, so it is not used. ---*/
    gpuErrChk(cudaMemsetAsync(d_prod, 0, nPointDomain * nVar * sizeof(ScalarType), ilu_stream));
    for (unsigned short sweep = 0; sweep < ilu_gpu_fwd_sweeps; ++sweep) {
      for (auto color = 0ul; color < ilu_color_ptr.size() - 1; ++color) {
        const auto begin = ilu_color_ptr[color];
        const auto size = ilu_color_ptr[color + 1] - begin;
        if (size == 0) continue;
        IluForwardLevelKernel<ScalarType>
            <<<size, threads, sharedForward, ilu_stream>>>(d_ilu_color_idx, begin, size, nVar, M, d_vec, d_prod);
      }
    }

    /*--- Backward substitution, levels in decreasing order. ---*/
    for (auto level = nLevels; level > 0;) {
      --level;  // unsigned type
      const auto begin = ilu_level_ptr[level];
      const auto size = ilu_level_ptr[level + 1] - begin;
      if (size == 0) continue;
      IluBackwardLevelKernel<ScalarType><<<size, threads, sharedBackward, ilu_stream>>>(
          d_ilu_level_idx, begin, size, nPointDomain, nVar, M, d_prod);
    }

    gpuErrChk(cudaStreamEndCapture(ilu_stream, &graph));
    gpuErrChk(cudaGraphInstantiate(&ilu_apply_graph_exec, graph, nullptr, nullptr, 0));
    gpuErrChk(cudaGraphDestroy(graph));
    ilu_apply_graph_vec = d_vec;
    ilu_apply_graph_prod = d_prod;
  }

  gpuErrChk(cudaGraphLaunch(ilu_apply_graph_exec, ilu_stream));
  gpuErrChk(cudaStreamSynchronize(ilu_stream));
  gpuErrChk(cudaGetLastError());
}

template <class ScalarType>
void CSysMatrix<ScalarType>::HtDTransfer(bool trigger) const {
  SU2_ZONE_SCOPED
  if (!trigger) return;
  gpuErrChk(cudaMemcpy(gpu.d, mat.d, sizeof(ScalarType) * nPoint * nVar * nEqn, cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy(gpu.l, mat.l, sizeof(ScalarType) * mat.nnz_l * nVar * nEqn, cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy(gpu.u, mat.u, sizeof(ScalarType) * mat.nnz_u * nVar * nEqn, cudaMemcpyHostToDevice));
}

template <class ScalarType>
void CSysMatrix<ScalarType>::MatrixVectorProductGPU(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                    CGeometry* geometry, const CConfig* config) const {
  SU2_ZONE_SCOPED
  if (nVar != nEqn) {
    SU2_MPI::Error("CUDA CSysMatrix block-LDU SpMV requires square blocks.", CURRENT_FUNCTION);
  }

  ScalarType* d_vec = vec.GetDevicePointer();
  ScalarType* d_prod = prod.GetDevicePointer();

  dim3 blockDim(static_cast<unsigned>(nVar), 1, 1);
  dim3 gridDim(static_cast<unsigned>(nPointDomain), 1, 1);
  BlockLDU_SpMV_kernel<ScalarType><<<gridDim, blockDim>>>(
      nPointDomain, nVar, gpu.row_ptr_l, gpu.col_ind_l, gpu.l, gpu.d,
      gpu.row_ptr_u, gpu.col_ind_u, gpu.u, d_vec, d_prod);
  /*--- Sync so the zone above actually times the kernel, not just the (async) launch call. ---*/
  gpuErrChk(cudaStreamSynchronize(nullptr));
  gpuErrChk(cudaGetLastError());
}

#define INSTANTIATE_MATRIX(TYPE)                                                            \
template void CSysMatrix<TYPE>::HtDTransfer(bool trigger) const;                            \
template void CSysMatrix<TYPE>::MatrixVectorProductGPU(const CSysVector<TYPE>& vec,         \
                                                       CSysVector<TYPE>& prod,              \
                                                       CGeometry* geometry,                 \
                                                       const CConfig* config) const;        \
template void CSysMatrix<TYPE>::BuildJacobiPreconditionerGPU();                             \
template void CSysMatrix<TYPE>::BuildILUPreconditionerGPU();                                \
template void CSysMatrix<TYPE>::ComputeILUPreconditionerGPU(const CSysVector<TYPE>& vec,    \
                                                            CSysVector<TYPE>& prod) const;  \
template void CSysMatrix<TYPE>::ComputeJacobiPreconditionerGPU(const CSysVector<TYPE>& vec, \
                                                               CSysVector<TYPE>& prod,      \
                                                               CGeometry* geometry,         \
                                                               const CConfig* config) const;
INSTANTIATE_MATRIX(su2mixedfloat)

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
INSTANTIATE_MATRIX(passivedouble)
#endif
