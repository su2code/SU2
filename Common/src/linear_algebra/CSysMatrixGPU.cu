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
#include <cstring>

#include "../../include/linear_algebra/CMatrixInverse.hpp"
#include "../../include/linear_algebra/CSysMatrix.inl"
#include "../../include/linear_algebra/GPUComms.cuh"

namespace {

/*!
 * \brief Apply the Jacobi preconditioner: prod = invM * vec (block-diagonal). Points are batched
 *        into blocks of ~128 threads (see ComputeJacobiPreconditionerGPU), threadIdx.x mapping
 *        to (point-within-block, output variable) via divmod by nVar - the same layout as
 *        BlockLDU_SpMV_kernel. Occupancy was already fine here (one thread per point, full
 *        warps), but consecutive threads used to land nVar^2 elements apart in invM (one point's
 *        whole dense block per thread); with this mapping they land nVar elements apart instead,
 *        a real (if partial, since it is still not stride-1) coalescing win that needs no shared
 *        memory or synchronization, unlike a fully-coalesced one-thread-per-block-entry version
 *        would.
 */
template <class ScalarType>
__global__ void ApplyJacobiPreconditionerKernel(const ScalarType* __restrict__ invM,
                                                const ScalarType* __restrict__ vec, ScalarType* __restrict__ prod,
                                                unsigned long nPointDomain, unsigned long nVar) {
  const unsigned long rowsPerBlock = blockDim.x / nVar;
  const unsigned long iVar = threadIdx.x % nVar;
  const unsigned long iPoint = static_cast<unsigned long>(blockIdx.x) * rowsPerBlock + threadIdx.x / nVar;
  if (iPoint >= nPointDomain) return;

  const auto* block = &invM[iPoint * nVar * nVar + iVar * nVar];
  const auto* rhs = &vec[iPoint * nVar];

  auto sum = ScalarType(0);
  for (unsigned long jVar = 0; jVar < nVar; ++jVar) sum += block[jVar] * rhs[jVar];
  prod[iPoint * nVar + iVar] = sum;
}

/*--- ILU. The factorization is scheduled by coloring: colors are true independent sets of the
 * ILU dependency graph (no two same-colored rows depend on each other in either direction), so a
 * color's rows can be processed with zero races in one kernel launch, but since a color is wider
 * and less ordered than a level, one pass over all colors is only an approximation, not an exact
 * result — several sweeps (repeated passes) are needed to converge. The triangular solves are
 * scheduled by level instead: a level's rows only depend on earlier levels (already finalized),
 * so one pass over the levels, in order, is exact — no sweeping needed there. Both the color and
 * level rows are scattered through the matrix, hence the indirection through their tables.
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
 * \brief Invert the diagonal blocks of the matrix, parallel-Gauss-Jordan device version of
 *        InverseDiagonalBlock. Unlike SU2_LinAlg::MatrixInverse (host, and still used by
 *        IluFactorColorKernel below) - which does forward elimination on one thread followed by
 *        a serial back-substitution, an inherently sequential ~nVar^3 chain - this eliminates
 *        each pivot column from every OTHER row simultaneously (both above and below the pivot,
 *        no back-substitution needed), so every entry update within one pivot step is independent
 *        of every other entry in that step: nVar serial steps (each with a couple of barriers),
 *        instead of one thread serially doing all the work. One thread per block entry (i,j);
 *        __syncthreads() (not __syncwarp()) is used so this stays correct even when nVar*nVar
 *        exceeds one warp (nVar > ~5), at the cost of a full block-wide barrier even for the
 *        common case where the whole block already fits in one warp.
 * \note Grid: one block per row, blockDim.x == nVar*nVar. Dynamic shared memory: 2*nVar*nVar
 *       scalars (the working copy of the block, and the inverse being accumulated in place of
 *       the old identity-then-eliminate scheme).
 */
template <class ScalarType>
__global__ void InvertDiagonalBlocksKernel(unsigned long nRows, unsigned long nVar,
                                           const ScalarType* __restrict__ mat_d, ScalarType* __restrict__ invM) {
  const unsigned long iRow = blockIdx.x;
  if (iRow >= nRows) return;

  const auto blockSize = nVar * nVar;
  const unsigned long tid = threadIdx.x;
  const unsigned long i = tid / nVar;
  const unsigned long j = tid % nVar;

  /*--- The inversion destroys its input, so it cannot work on the matrix itself. ---*/
  extern __shared__ __align__(sizeof(double)) char smem[];
  auto* A = reinterpret_cast<ScalarType*>(smem);
  auto* Inv = A + blockSize;

  A[tid] = mat_d[iRow * blockSize + tid];
  Inv[tid] = ScalarType(i == j);
  __syncthreads();

  for (auto k = 0ul; k < nVar; ++k) {
    /*--- Regularize the pivot (shared with the host path, same clamp value). ---*/
    if (i == k && j == k) SU2_LinAlg::RegularizePivot(A[k * nVar + k]);
    __syncthreads();

    /*--- Normalize the pivot row. ---*/
    const ScalarType pivot = A[k * nVar + k];
    if (i == k) {
      A[tid] /= pivot;
      Inv[tid] /= pivot;
    }
    __syncthreads();

    /*--- Eliminate column k from every other row; A(k,*) and Inv(k,*) are already finalized for
     * this step (previous barrier), and each thread only ever writes its own (i,j), so this
     * needs no further synchronization until the next pivot's regularization reads A(k+1,k+1). ---*/
    if (i != k) {
      const ScalarType factor = A[i * nVar + k];
      A[tid] -= factor * A[k * nVar + j];
      Inv[tid] -= factor * Inv[k * nVar + j];
    }
    __syncthreads();
  }

  invM[iRow * blockSize + tid] = Inv[tid];
}

/*!
 * \brief Quantize the diagonal blocks straight from the device diagonal (gpu.d), device
 *        counterpart of CSysMatrix::QuantizeBlock (CSysMatrix.cpp), applied row by row. Calls the
 *        exact same EncodeQuantBlock encoding routine (CSysMatrix.hpp, SU2_CUDA_HOST_DEVICE) that
 *        the host path uses, rather than a separate device copy of the encoding logic. One
 *        thread per row.
 */
template <class ScalarType>
__global__ void QuantizeDiagonalBlocksKernel(unsigned long nRows, unsigned long nVar,
                                             const ScalarType* __restrict__ mat_d, int8_t* __restrict__ q_scale_d,
                                             int8_t* __restrict__ q_blocks_d) {
  const auto iRow = static_cast<unsigned long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (iRow >= nRows) return;

  const auto* blk = mat_d + iRow * nVar * nVar;
  EncodeQuantBlock([&](unsigned long r, unsigned long c) { return blk[r * nVar + c]; }, q_scale_d + iRow * nVar,
                   q_blocks_d + iRow * nVar * nVar, nVar);
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
 * \brief Exact forward substitution for the rows of one level, (L+I).prod = vec.
 * \note Every row in a level only depends on rows in earlier levels, which are already
 *       finalized (see CSysMatrix::levels_ilu), so one pass over the levels in increasing order
 *       gives the exact result, unlike the colored factorization above. One thread per block
 *       entry, so the inner dot product over a neighbor block is spread across nVar threads
 *       instead of done serially by one; each thread accumulates its own (iVar,jVar) partial
 *       product across every neighbor with no synchronization at all, and only the final
 *       nVar-way reduction (summing over jVar for each iVar) needs one __syncthreads(). Grid:
 *       one block per row of the level, blockDim.x == nVar*nVar. Dynamic shared memory:
 *       nVar*nVar scalars.
 */
template <class ScalarType>
__global__ void IluForwardKernel(const su2uint* __restrict__ level_idx, unsigned long level_begin,
                                 unsigned long level_size, unsigned long nVar, DeviceLDU<ScalarType> M,
                                 const ScalarType* __restrict__ vec, ScalarType* __restrict__ prod) {
  if (blockIdx.x >= level_size) return;

  const unsigned long iRow = level_idx[level_begin + blockIdx.x];
  const unsigned long tid = threadIdx.x;
  const auto iVar = tid / nVar, jVar = tid % nVar;

  extern __shared__ __align__(sizeof(double)) char smem[];
  auto* partial = reinterpret_cast<ScalarType*>(smem);

  ScalarType acc = 0;
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
 * \brief Exact backward substitution for the rows of one level, U.prod = prod.
 * \note The right-hand side is read directly from \p prod: a level is visited exactly once, so
 *       prod[iRow] is still the untouched forward-solve result when its row is processed (unlike
 *       a colored sweep, which revisits every row and would need a separate fixed buffer to tell
 *       the right-hand side apart from a solution estimate). Levels are processed in decreasing
 *       order so every U-neighbor (a higher row index) is already finalized. Same thread layout
 *       as IluForwardKernel, plus one extra __syncthreads() before the diagonal multiply (which
 *       needs every iVar's reduced sum). Grid: one block per row of the level,
 *       blockDim.x == nVar*nVar. Dynamic shared memory: nVar*nVar + nVar scalars.
 */
template <class ScalarType>
__global__ void IluBackwardKernel(const su2uint* __restrict__ level_idx, unsigned long level_begin,
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
    aux[iVar] = sum;
  }
  __syncthreads();

  if (jVar == 0) {
    /*--- The diagonal blocks are stored inverted by the factorization. ---*/
    const auto* invUii = M.d + iRow * blockSize;
    ScalarType out = 0;
    for (auto k = 0ul; k < nVar; ++k) out += invUii[iVar * nVar + k] * aux[k];
    prod[iRow * nVar + iVar] = out;
  }
}

/*!
 * \brief Block-LDU SpMV kernel: y[iRow] = (L + D + U) * x per block-row. Several rows are
 *        batched into one CUDA block (blockDim.x / nVar of them, see MatrixVectorProductGPU)
 *        instead of one row per block: nVar is typically ~4-6, so one-row-per-block leaves most
 *        of a warp's lanes permanently idle and caps occupancy at a few resident (mostly-empty)
 *        warps per SM, well before DRAM bandwidth is the limit. threadIdx.x indexes
 *        (row-within-block, output variable) as (threadIdx.x / nVar, threadIdx.x % nVar).
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
  const unsigned long rowsPerBlock = blockDim.x / nVar;
  const unsigned long iVar = threadIdx.x % nVar;
  const unsigned long iRow = static_cast<unsigned long>(blockIdx.x) * rowsPerBlock + threadIdx.x / nVar;
  if (iRow >= nRows) return;

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

/*!
 * \brief Quantized block-LDU SpMV kernel: y[iRow] = (L + D + U) * x per block-row, reading int8
 *        row-scaled quantized blocks instead of full precision ones. Device version of
 *        QuantizedRowProduct/QuantizedMatVecAdd (CSysMatrix.inl). Rows are batched per block the
 *        same way as BlockLDU_SpMV_kernel (see its comment). Each row of quantized mantissas is
 *        also read 4 bytes at a time (one 32-bit load instead of four 8-bit ones): dp4a does not
 *        apply here since only the matrix is quantized, not x, so the accumulation itself is
 *        still one scalar FMA per element in the exact same order as the scalar loop (bit-
 *        identical result) - the win is fewer, wider load instructions, not fewer FMAs.
 */
template <class ScalarType>
__global__ void QuantizedBlockLDU_SpMV_kernel(
    unsigned long nRows, unsigned long nVar, const su2uint* __restrict__ row_ptr_l,
    const su2uint* __restrict__ col_ind_l, const int8_t* __restrict__ q_scale_l,
    const int8_t* __restrict__ q_blocks_l, const int8_t* __restrict__ q_scale_d,
    const int8_t* __restrict__ q_blocks_d, const su2uint* __restrict__ row_ptr_u,
    const su2uint* __restrict__ col_ind_u, const int8_t* __restrict__ q_scale_u,
    const int8_t* __restrict__ q_blocks_u, const ScalarType* __restrict__ x, ScalarType* __restrict__ y) {
  const unsigned long rowsPerBlock = blockDim.x / nVar;
  const unsigned long iVar = threadIdx.x % nVar;
  const unsigned long iRow = static_cast<unsigned long>(blockIdx.x) * rowsPerBlock + threadIdx.x / nVar;
  if (iRow >= nRows) return;

  auto addBlock = [&](const int8_t* __restrict__ qs, const int8_t* __restrict__ qv, const ScalarType* __restrict__ xk) {
    const float row_scale = DecodeQuantScale(qs[iVar]);
    const int8_t* __restrict__ row = qv + iVar * nVar;
    ScalarType partial = 0;
    unsigned long jVar = 0;
    for (; jVar + 4 <= nVar; jVar += 4) {
      /*--- Row bytes are not generally 4-byte aligned (nVar*nVar need not be a multiple of 4),
       * so this must go through memcpy rather than a reinterpret_cast<const uint32_t*> deref. ---*/
      uint32_t packed;
      memcpy(&packed, row + jVar, sizeof(packed));
      partial += static_cast<int8_t>(packed) * xk[jVar];
      partial += static_cast<int8_t>(packed >> 8) * xk[jVar + 1];
      partial += static_cast<int8_t>(packed >> 16) * xk[jVar + 2];
      partial += static_cast<int8_t>(packed >> 24) * xk[jVar + 3];
    }
    for (; jVar < nVar; ++jVar) partial += row[jVar] * xk[jVar];
    return static_cast<ScalarType>(row_scale) * partial;
  };

  ScalarType sum = 0;
  /* Lower */
  for (auto k = row_ptr_l[iRow]; k < row_ptr_l[iRow + 1]; ++k) {
    const auto col = col_ind_l[k];
    sum += addBlock(q_scale_l + k * nVar, q_blocks_l + k * nVar * nVar, x + col * nVar);
  }
  /* Diagonal */
  sum += addBlock(q_scale_d + iRow * nVar, q_blocks_d + iRow * nVar * nVar, x + iRow * nVar);
  /* Upper */
  for (auto k = row_ptr_u[iRow]; k < row_ptr_u[iRow + 1]; ++k) {
    const auto col = col_ind_u[k];
    sum += addBlock(q_scale_u + k * nVar, q_blocks_u + k * nVar * nVar, x + col * nVar);
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

  constexpr unsigned long targetThreadsPerBlock = 128;
  const auto rowsPerBlock = std::max<unsigned long>(1, targetThreadsPerBlock / nVar);
  const auto threadsPerBlock = static_cast<unsigned>(rowsPerBlock * nVar);
  const auto blocks = static_cast<unsigned>((nPointDomain + rowsPerBlock - 1) / rowsPerBlock);
  ApplyJacobiPreconditionerKernel<<<blocks, threadsPerBlock>>>(d_invM, vec.GetDevicePointer(), prod.GetDevicePointer(),
                                                               nPointDomain, nVar);
  /*--- Sync so the zone above actually times the kernel, not just the (async) launch call. ---*/
  gpuErrChk(cudaStreamSynchronize(nullptr));
  gpuErrChk(cudaGetLastError());
}

template <class ScalarType>
void CSysMatrix<ScalarType>::QuantizeDiagonalBlocksGPU() {
  SU2_ZONE_SCOPED

  if (nPointDomain == 0) return;

  /*--- The matrix is expected to be on the device already, it is uploaded once per solve by
   * CSysMatrixVectorProduct, which is created before the preconditioner is built. ---*/
  constexpr unsigned threadsPerBlock = 128;
  const auto blocks = static_cast<unsigned>((nPointDomain + threadsPerBlock - 1) / threadsPerBlock);
  QuantizeDiagonalBlocksKernel<ScalarType>
      <<<blocks, threadsPerBlock>>>(nPointDomain, nVar, gpu.d, d_q_scale.d, d_q_blocks.d);
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
  if (nPointDomain == 0) return;

  /*--- The matrix is expected to be on the device already, it is uploaded once per solve by
   * CSysMatrixVectorProduct, which is created before the preconditioner is built. ---*/
  const auto blockSize = static_cast<unsigned>(nVar * nVar);
  InvertDiagonalBlocksKernel<ScalarType><<<static_cast<unsigned>(nPointDomain), blockSize,
                                           2 * blockSize * sizeof(ScalarType)>>>(nPointDomain, nVar, gpu.d, d_invM);
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

  /*--- The launch sequence (ilu_gpu_sweeps passes over all colors) is identical on every call:
   * the grid and block sizes only depend on the (fixed) sparsity pattern/coloring and the device
   * pointers are fixed members, allocated once. Capture it into a CUDA graph the first time and
   * replay that from then on, which removes the per-launch host-side overhead without touching
   * the parallelization of any individual kernel (unlike a persistent cooperative-groups kernel,
   * this does not cap per-color parallelism to an occupancy-resident block count). See
   * IluFactorColorKernel for why several sweeps over the colors are needed. Note that factors
   * are not reset between calls to BuildILUPreconditionerGPU, so with LINEAR_SOLVER_ILU_GPU_SWEEPS
   * set low (even 1), each call refines the previous one's result rather than reconverging from
   * scratch, relying on the matrix changing little between outer/pseudo-time iterations. ---*/
  if (ilu_build_graph_exec == nullptr) {
    cudaGraph_t graph;
    gpuErrChk(cudaStreamBeginCapture(ilu_stream, cudaStreamCaptureModeThreadLocal));

    for (unsigned short sweep = 0; sweep < ilu_gpu_sweeps; ++sweep) {
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
   * captured once and replayed; if the pointers ever do change the graph is recaptured, which is
   * no worse than the un-graphed loop, just not free. ---*/
  if (ilu_apply_graph_exec == nullptr || ilu_apply_graph_vec != d_vec || ilu_apply_graph_prod != d_prod) {
    if (ilu_apply_graph_exec != nullptr) {
      gpuErrChk(cudaGraphExecDestroy(ilu_apply_graph_exec));
      ilu_apply_graph_exec = nullptr;
    }

    cudaGraph_t graph;
    gpuErrChk(cudaStreamBeginCapture(ilu_stream, cudaStreamCaptureModeThreadLocal));

    /*--- Forward substitution: one exact pass over the levels in increasing order,
     * (L+I).prod = vec, see IluForwardKernel. ---*/
    for (auto level = 0ul; level < nLevels; ++level) {
      const auto begin = ilu_level_ptr[level];
      const auto size = ilu_level_ptr[level + 1] - begin;
      if (size == 0) continue;
      IluForwardKernel<ScalarType>
          <<<size, threads, sharedForward, ilu_stream>>>(d_ilu_level_idx, begin, size, nVar, M, d_vec, d_prod);
    }

    /*--- Backward substitution: one exact pass over the levels in decreasing order,
     * U.prod = prod, see IluBackwardKernel. ---*/
    for (auto level = nLevels; level > 0;) {
      --level;
      const auto begin = ilu_level_ptr[level];
      const auto size = ilu_level_ptr[level + 1] - begin;
      if (size == 0) continue;
      IluBackwardKernel<ScalarType>
          <<<size, threads, sharedBackward, ilu_stream>>>(d_ilu_level_idx, begin, size, nPointDomain, nVar, M, d_prod);
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
  if (quantized_mode) {
    /*--- No gpu.l/gpu.u to transfer (never allocated); mirror the host quantized off-diagonal
     * storage instead (the diagonal mirrors, d_q_scale.d/d_q_blocks.d, are not touched here at
     * all - QuantizeDiagonalBlocksGPU() populates them straight from gpu.d, just uploaded above,
     * with no host round trip). Issued as async copies on the default stream, then left in
     * flight: the caller (CSysMatrixVectorProduct's constructor) returns right after this, and
     * whatever the preconditioner's Build() does next runs concurrently with the transfer still
     * draining. This is only genuinely asynchronous (i.e. the host thread does not block here
     * waiting for the copy) because q_scale.l/q_blocks.l/q_scale.u/q_blocks.u are pinned host
     * memory, see the comment on those members / Initialize() - cudaMemcpyAsync silently
     * degrades to a blocking copy from regular pageable memory. Any later kernel that reads
     * d_q_scale.l/.u/d_q_blocks.l/.u is issued on the same default stream too, so it correctly
     * waits for these by stream ordering without an explicit sync here. ---*/
    gpuErrChk(cudaMemcpyAsync(d_q_scale.l, q_scale.l, sizeof(QuantType) * mat.nnz_l * nVar, cudaMemcpyHostToDevice));
    gpuErrChk(cudaMemcpyAsync(d_q_blocks.l, q_blocks.l, sizeof(QuantType) * mat.nnz_l * nVar * nEqn,
                              cudaMemcpyHostToDevice));
    gpuErrChk(cudaMemcpyAsync(d_q_scale.u, q_scale.u, sizeof(QuantType) * mat.nnz_u * nVar, cudaMemcpyHostToDevice));
    gpuErrChk(cudaMemcpyAsync(d_q_blocks.u, q_blocks.u, sizeof(QuantType) * mat.nnz_u * nVar * nEqn,
                              cudaMemcpyHostToDevice));
  } else {
    gpuErrChk(cudaMemcpy(gpu.l, mat.l, sizeof(ScalarType) * mat.nnz_l * nVar * nEqn, cudaMemcpyHostToDevice));
    gpuErrChk(cudaMemcpy(gpu.u, mat.u, sizeof(ScalarType) * mat.nnz_u * nVar * nEqn, cudaMemcpyHostToDevice));
  }
}

template <class ScalarType>
void CSysMatrix<ScalarType>::MatrixVectorProductGPU(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                    CGeometry* geometry, const CConfig* config) const {
  SU2_ZONE_SCOPED

  ScalarType* d_vec = vec.GetDevicePointer();
  ScalarType* d_prod = prod.GetDevicePointer();

  /*--- Batch several rows per block (see BlockLDU_SpMV_kernel's comment): nVar is small
   * (typically ~4-6), so one row per block would leave most of a warp idle. Aim for ~128
   * threads/block, the largest whole number of rows that fits. ---*/
  constexpr unsigned long targetThreadsPerBlock = 128;
  const auto rowsPerBlock = std::max<unsigned long>(1, targetThreadsPerBlock / nVar);
  dim3 blockDim(static_cast<unsigned>(rowsPerBlock * nVar), 1, 1);
  dim3 gridDim(static_cast<unsigned>((nPointDomain + rowsPerBlock - 1) / rowsPerBlock), 1, 1);
  if (quantized_mode) {
    QuantizedBlockLDU_SpMV_kernel<ScalarType><<<gridDim, blockDim>>>(
        nPointDomain, nVar, gpu.row_ptr_l, gpu.col_ind_l, d_q_scale.l, d_q_blocks.l, d_q_scale.d, d_q_blocks.d,
        gpu.row_ptr_u, gpu.col_ind_u, d_q_scale.u, d_q_blocks.u, d_vec, d_prod);
  } else {
    BlockLDU_SpMV_kernel<ScalarType><<<gridDim, blockDim>>>(
        nPointDomain, nVar, gpu.row_ptr_l, gpu.col_ind_l, gpu.l, gpu.d,
        gpu.row_ptr_u, gpu.col_ind_u, gpu.u, d_vec, d_prod);
  }
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
template void CSysMatrix<TYPE>::QuantizeDiagonalBlocksGPU();                                \
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
