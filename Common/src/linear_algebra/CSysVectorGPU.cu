/*!
 * \file CSysVectorGPU.cu
 * \brief Implementations of Kernels and Functions for Vector Operations on the GPU
 * \author A. Raj, D. Di giusto
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

#include "../../include/linear_algebra/CSysVector.hpp"
#include "../../include/linear_algebra/GPUComms.cuh"
#include <cublas_v2.h>
#include <algorithm>
#include <type_traits>
#include <vector>

namespace {

/*--- cuBLAS handle for the reductions. Created on first use and kept for the lifetime of
 * the program, matching the fact that CUDA is either on or off for the whole run. ---*/
cublasHandle_t solver_blas_handle = nullptr;

cublasHandle_t GetBlasHandle() {
  if (solver_blas_handle == nullptr) {
    if (cublasCreate(&solver_blas_handle) != CUBLAS_STATUS_SUCCESS) {
      SU2_MPI::Error("cuBLAS handle creation failed for the GPU linear algebra.", CURRENT_FUNCTION);
    }
  }
  return solver_blas_handle;
}

/*!
 * \brief Cap on n*m (the number of dot products one MultiDotKernel launch computes). Bounds the
 *        per-thread accumulator array below (necessarily a compile-time size - CUDA has no
 *        runtime-sized register/local arrays) and the block's shared-memory reduction buffer.
 *        FGCRODR's Ritz-value path (CSysSolve.cpp) calls this with n = m+1 where m is
 *        LINEAR_SOLVER_RESTART_DEFLATION (user-configurable, default 4 but not uncommon to raise
 *        into the tens), so n*m grows quadratically with that setting - e.g. m=10 already needs
 *        110. Sized generously above that; raise it further if a caller ever needs more (a clear
 *        SU2_MPI::Error fires instead of silently producing a wrong/truncated result). This is a
 *        different constraint from MULTIDOT_MAX_VEC below (n and m individually can each reach
 *        that cap, but not both at once - n*m up to MULTIDOT_MAX_VEC^2 would need an
 *        impractically large per-thread array).
 */
constexpr unsigned int MULTIDOT_MAX_NM = 1024;

/*!
 * \brief Cap on the individual vector counts n and m. Bounds the fixed-size pointer arrays
 *        passed into MultiDotKernel by value as ordinary launch parameters, rather than as a
 *        device-side array of pointers uploaded via a separate cudaMemcpy - CUDA already
 *        marshals launch parameters for you as part of the launch itself, regardless of host
 *        memory pinning, which sidesteps needing pinned host buffers here the way HtDTransfer's
 *        L/U copy does (see aux_stream/htd_event in CSysMatrix, CUDA silently downgrades
 *        cudaMemcpyAsync to a blocking copy from ordinary pageable host memory). n or m can
 *        independently reach the Krylov restart length in this codebase's two call sites
 *        (ModGramSchmidt in CSysSolve.cpp grows m up to the restart length with n=1; FGCRODR's
 *        Ritz-value path grows n up to the restart length with m bounded by
 *        LINEAR_SOLVER_RESTART_DEFLATION), so this has to cover restart length, not just the
 *        (usually much smaller) deflation count - 128 is generous headroom over realistic
 *        restart lengths (typically 10-50).
 */
constexpr unsigned int MULTIDOT_MAX_VEC = 128;

/*!
 * \brief Fixed-size array of device pointers, passed to MultiDotKernel by value (see
 *        MULTIDOT_MAX_VEC for why).
 */
template <class ScalarType>
struct MultiDotPointers {
  const ScalarType* ptr[MULTIDOT_MAX_VEC];
};

/*!
 * \brief Compute the n*m matrix of dot products C(i,j) = <V[i], W[j]> in a single pass over the
 *        data, replacing a previous cublas<t>gemmBatched implementation that treated each (i,j)
 *        pair as its own independent (M=1, N=1, K=size) GEMM. That was a poor fit twice over:
 *        (a) it re-read every input vector once per pair instead of once total (n*m reads of
 *        length size instead of n+m), and (b) GEMM kernels parallelize across output (M,N)
 *        tiles and iterate the reduction (K) dimension within one thread block's kernel
 *        invocation - for an M=N=1 tile with K in the tens of millions, that leaves most of the
 *        GPU idle, unlike a dedicated reduction primitive (cublas<t>dot, or this kernel's
 *        grid-stride loop across many blocks) which parallelizes across all of K.
 * \note Each thread accumulates its own private running n*m sums while striding over k (kept in
 *       thread-local storage, capped at MULTIDOT_MAX_NM - small enough that it should stay
 *       resident in registers or L1 for realistic n*m, cheap either way next to the K-length
 *       main loop's DRAM traffic), so the single pass over V/W stays memory-bound. Only after
 *       that loop do threads combine: a warp-shuffle reduction, then one shared-memory combine
 *       and one atomicAdd per (i,j) per block, so total atomics are O(nm * numBlocks), not
 *       O(nm * size).
 */
template <class ScalarType>
__global__ void MultiDotKernel(MultiDotPointers<ScalarType> V, unsigned int n, MultiDotPointers<ScalarType> W,
                               unsigned int m, unsigned long size, ScalarType* __restrict__ C) {
  const unsigned int nm = n * m;

  ScalarType local[MULTIDOT_MAX_NM];
  for (unsigned int t = 0; t < nm; ++t) local[t] = ScalarType(0);

  const auto stride = static_cast<unsigned long>(blockDim.x) * gridDim.x;
  for (auto k = static_cast<unsigned long>(blockIdx.x) * blockDim.x + threadIdx.x; k < size; k += stride) {
    for (unsigned int i = 0; i < n; ++i) {
      const ScalarType v = V.ptr[i][k];
      for (unsigned int j = 0; j < m; ++j) local[i * m + j] += v * W.ptr[j][k];
    }
  }

  /*--- Warp-level reduction: after this, lane 0 of every warp holds that warp's true sum. ---*/
  for (unsigned int t = 0; t < nm; ++t) {
    ScalarType val = local[t];
    for (int offset = 16; offset > 0; offset >>= 1) val += __shfl_down_sync(0xFFFFFFFFu, val, offset);
    local[t] = val;
  }

  /*--- Block-level combine: warp leaders stage their partials in shared memory, thread 0 sums
   * them and issues the one atomicAdd per (i,j) this block contributes to the global result
   * (pre-zeroed by the caller). ---*/
  extern __shared__ __align__(sizeof(double)) char smem[];
  auto* warpPartials = reinterpret_cast<ScalarType*>(smem);
  const unsigned int lane = threadIdx.x % 32u;
  const unsigned int warpId = threadIdx.x / 32u;
  const unsigned int warpsPerBlock = blockDim.x / 32u;

  if (lane == 0) {
    for (unsigned int t = 0; t < nm; ++t) warpPartials[warpId * nm + t] = local[t];
  }
  __syncthreads();

  if (threadIdx.x == 0) {
    for (unsigned int t = 0; t < nm; ++t) {
      ScalarType sum = 0;
      for (unsigned int w = 0; w < warpsPerBlock; ++w) sum += warpPartials[w * nm + t];
      atomicAdd(&C[t], sum);
    }
  }
}

}  // namespace

namespace VecExpr {

namespace {
/*--- Deliberately not thread local: every thread of an OpenMP team has to agree on this,
 * otherwise the team splits over the worksharing constructs in CSysVector. It is only
 * written by CSysSolve, outside any parallel region over the linear system. ---*/
bool use_device_expressions = false;
}  // namespace

bool UseDeviceExpressions() { return use_device_expressions; }

void SetUseDeviceExpressions(bool use) { use_device_expressions = use; }

}  // namespace VecExpr

template <class ScalarType>
void CSysVector<ScalarType>::HtDTransfer(bool trigger) const {
  SU2_ZONE_SCOPED
  if (trigger)
    gpuErrChk(cudaMemcpy((void*)(d_vec_val), (void*)&vec_val[0], (sizeof(ScalarType) * nElm), cudaMemcpyHostToDevice));
}

template <class ScalarType>
void CSysVector<ScalarType>::DtHTransfer(bool trigger) const {
  SU2_ZONE_SCOPED
  if (trigger)
    gpuErrChk(cudaMemcpy((void*)(&vec_val[0]), (void*)d_vec_val, (sizeof(ScalarType) * nElm), cudaMemcpyDeviceToHost));
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::dotGPU(const CSysVector& other) const {
  SU2_ZONE_SCOPED
  /*--- Both operands are already on the device, the caller owns the transfers. This
   * reduces over MPI, so it must be called by a single thread (see SU2_DEVICE_REGION). ---*/
  cublasHandle_t handle = GetBlasHandle();
  cublasStatus_t status = CUBLAS_STATUS_SUCCESS;

  ScalarType local_dot = ScalarType(0);

  if constexpr (std::is_same_v<ScalarType, float>) {
    status = cublasSdot(handle, static_cast<int>(nElmDomain), GetDevicePointer(), 1, other.GetDevicePointer(), 1,
                        &local_dot);
  } else if constexpr (std::is_same_v<ScalarType, double>) {
    status = cublasDdot(handle, static_cast<int>(nElmDomain), GetDevicePointer(), 1, other.GetDevicePointer(), 1,
                        &local_dot);
  } else {
    SU2_MPI::Error("Unsupported ScalarType in CSysVector::dotGPU.", CURRENT_FUNCTION);
    return ScalarType(0);
  }

  if (status != CUBLAS_STATUS_SUCCESS) {
    SU2_MPI::Error("cuBLAS dot failed in CSysVector::dotGPU.", CURRENT_FUNCTION);
    return ScalarType(0);
  }

  ScalarType global_dot = ScalarType(0);
  const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
  SelectMPIWrapper<ScalarType>::W::Allreduce(&local_dot, &global_dot, 1, mpi_type, MPI_SUM, SU2_MPI::GetComm());

  return global_dot;
}

/*!
 * \brief Multi vector dot product, C(i,j) = <V[i0+i], W[j]>, via MultiDotKernel (see its comment
 *        for why that beats the batched-GEMM approach this replaced).
 */
template <class ScalarType>
su2matrix<ScalarType> CSysVector<ScalarType>::multiDotGPU(const std::vector<CSysVector<ScalarType>>& V, const size_t i0,
                                                          const size_t n, const std::vector<CSysVector<ScalarType>>& W,
                                                          const size_t m) {
  if (n > MULTIDOT_MAX_VEC || m > MULTIDOT_MAX_VEC) {
    SU2_MPI::Error("CSysVector::multiDotGPU: n or m exceeds MULTIDOT_MAX_VEC, raise that constant in "
                   "CSysVectorGPU.cu.",
                   CURRENT_FUNCTION);
  }
  const size_t nm = n * m;
  if (nm > MULTIDOT_MAX_NM) {
    SU2_MPI::Error("CSysVector::multiDotGPU: n*m exceeds MULTIDOT_MAX_NM, raise that constant in "
                   "CSysVectorGPU.cu.",
                   CURRENT_FUNCTION);
  }

  const size_t size = V[0].nElmDomain;

  su2matrix<ScalarType> local;
  local.resize(n, m);
  if (nm == 0) return local;

  /*--- Persistent device workspace for the output only, cached across calls and freed
   * automatically when the program exits (static local destruction), instead of leaking. The
   * V/W pointers themselves need no device buffer at all: they go to MultiDotKernel as ordinary
   * by-value launch parameters (see MultiDotPointers/MULTIDOT_MAX_VEC), which CUDA marshals
   * itself as part of the launch - unlike a separate cudaMemcpy, that needs no pinned host
   * buffer to actually be asynchronous. ---*/
  struct Workspace {
    ScalarType* d_C = nullptr;
    cudaStream_t stream = nullptr;
    size_t cCapacity = 0;

    void EnsureCapacity(size_t nm) {
      if (stream == nullptr) gpuErrChk(cudaStreamCreate(&stream));
      if (nm > cCapacity) {
        cudaFree(d_C);
        gpuErrChk(cudaMalloc(&d_C, nm * sizeof(ScalarType)));
        cCapacity = nm;
      }
    }

    ~Workspace() {
      cudaFree(d_C);
      if (stream != nullptr) cudaStreamDestroy(stream);
    }
  };
  static Workspace ws;
  ws.EnsureCapacity(nm);

  MultiDotPointers<ScalarType> vPtrs{}, wPtrs{};
  for (size_t i = 0; i < n; ++i) vPtrs.ptr[i] = V[i0 + i].GetDevicePointer();
  for (size_t j = 0; j < m; ++j) wPtrs.ptr[j] = W[j].GetDevicePointer();

  gpuErrChk(cudaMemsetAsync(ws.d_C, 0, nm * sizeof(ScalarType), ws.stream));

  constexpr unsigned int threadsPerBlock = 256;
  const auto blocks = static_cast<unsigned int>(std::min<size_t>((size + threadsPerBlock - 1) / threadsPerBlock, 1024));
  const auto sharedBytes = static_cast<size_t>(threadsPerBlock / 32u) * nm * sizeof(ScalarType);

  /*--- sharedBytes can exceed the default 48KB static shared-memory limit for large nm (the
   * FGCRODR deflation matrix in particular, see MULTIDOT_MAX_NM); opt in to the device's larger
   * "dynamic" limit once, the first time it is actually needed, rather than always paying for
   * the query. ---*/
  static size_t optedInSharedBytes = 0;
  if (sharedBytes > optedInSharedBytes) {
    gpuErrChk(cudaFuncSetAttribute(MultiDotKernel<ScalarType>, cudaFuncAttributeMaxDynamicSharedMemorySize,
                                   static_cast<int>(sharedBytes)));
    optedInSharedBytes = sharedBytes;
  }

  MultiDotKernel<ScalarType><<<blocks, threadsPerBlock, sharedBytes, ws.stream>>>(
      vPtrs, static_cast<unsigned int>(n), wPtrs, static_cast<unsigned int>(m), size, ws.d_C);

  gpuErrChk(cudaMemcpyAsync(local.data(), ws.d_C, nm * sizeof(ScalarType), cudaMemcpyDeviceToHost, ws.stream));
  gpuErrChk(cudaStreamSynchronize(ws.stream));
  gpuErrChk(cudaGetLastError());

  return local;
}

/*--- Every expression the solvers assign to a CSysVector needs its assignment kernel
 * instantiated here; the host compiler cannot emit one. A shape that is missing shows up
 * as an undefined reference to VecExpr::AssignDeviceExpression at link time, and is fixed
 * by adding a line to DEVICE_EXPRESSION_SHAPES below. The aliases use CSysVector (not
 * CVectorView) because that is how the operator overloads name their operands; store_t
 * turns it into a view when the node is built. ---*/
namespace {

template <class S>
using Vec = CSysVector<S>;
template <class S>
using Sca = VecExpr::Bcast<S>;

/*--- Leaves and the shapes of the FGMRES/GMRES basis updates. ---*/
template <class S>
using DeviceBcast = Sca<S>;
template <class S>
using DeviceView = VecExpr::CVectorView<S>;
template <class S>
using DeviceNeg = VecExpr::minus_<Vec<S>, S>;

/*--- vector * scalar and scalar * vector are distinct types, both are used. ---*/
template <class S>
using DeviceScale = VecExpr::mul_<Vec<S>, Sca<S>, S>;
template <class S>
using DeviceLScale = VecExpr::mul_<Sca<S>, Vec<S>, S>;
template <class S>
using DeviceDivScale = VecExpr::div_<Vec<S>, Sca<S>, S>;

/*--- Linear combinations, CSysSolve unrolls them up to four terms. ---*/
template <class S>
using DeviceScale2 = VecExpr::add_<DeviceScale<S>, DeviceScale<S>, S>;
template <class S>
using DeviceScale3 = VecExpr::add_<DeviceScale2<S>, DeviceScale<S>, S>;
template <class S>
using DeviceScale4 = VecExpr::add_<DeviceScale3<S>, DeviceScale<S>, S>;

/*--- r = b - A_x, in CG, BCGSTAB, Smoother and FGCRODR. ---*/
template <class S>
using DeviceSub = VecExpr::sub_<Vec<S>, Vec<S>, S>;

/*--- p = beta * p + z, in CG. ---*/
template <class S>
using DeviceLScalePlus = VecExpr::add_<DeviceLScale<S>, Vec<S>, S>;

/*--- p = beta * (p - omega * v) + r, in BCGSTAB. ---*/
template <class S>
using DeviceSubLScale = VecExpr::sub_<Vec<S>, DeviceLScale<S>, S>;
template <class S>
using DeviceLScaleSub = VecExpr::mul_<Sca<S>, DeviceSubLScale<S>, S>;
template <class S>
using DeviceBcgsDir = VecExpr::add_<DeviceLScaleSub<S>, Vec<S>, S>;

}  // namespace

#define DEVICE_EXPRESSION_SHAPES(SCALAR)                  \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceBcast);    \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceView);     \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceNeg);      \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceScale);    \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceLScale);   \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceDivScale); \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceScale2);   \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceScale3);   \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceScale4);   \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceSub);      \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceLScalePlus); \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceSubLScale);  \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceLScaleSub);  \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceBcgsDir)

#define INSTANTIATE_DEVICE_ASSIGN(SCALAR, OP, EXPR)                                          \
  template void VecExpr::AssignDeviceExpression<VecExpr::DeviceAssignOp::OP, SCALAR, EXPR<SCALAR>>( \
      SCALAR*, unsigned long, const VecExpr::CVecExpr<EXPR<SCALAR>, SCALAR>&)

#define INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, EXPR) \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Assign, EXPR);   \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Add, EXPR);      \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Subtract, EXPR); \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Multiply, EXPR); \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Divide, EXPR)

DEVICE_EXPRESSION_SHAPES(su2mixedfloat);

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
DEVICE_EXPRESSION_SHAPES(passivedouble);
#endif

#undef DEVICE_EXPRESSION_SHAPES
#undef INSTANTIATE_DEVICE_ASSIGN_EXPR
#undef INSTANTIATE_DEVICE_ASSIGN


template void CSysVector<su2mixedfloat>::HtDTransfer(bool trigger) const;
template void CSysVector<su2mixedfloat>::DtHTransfer(bool trigger) const;
template su2mixedfloat CSysVector<su2mixedfloat>::dotGPU(const CSysVector<su2mixedfloat>& other) const;
template su2matrix<su2mixedfloat> CSysVector<su2mixedfloat>::multiDotGPU(
    const std::vector<CSysVector<su2mixedfloat>>& V, size_t i0, size_t n,
    const std::vector<CSysVector<su2mixedfloat>>& W, size_t m);

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template void CSysVector<passivedouble>::HtDTransfer(bool trigger) const;
template void CSysVector<passivedouble>::DtHTransfer(bool trigger) const;
template passivedouble CSysVector<passivedouble>::dotGPU(const CSysVector<passivedouble>& other) const;
template su2matrix<passivedouble> CSysVector<passivedouble>::multiDotGPU(
    const std::vector<CSysVector<passivedouble>>& V, size_t i0, size_t n,
    const std::vector<CSysVector<passivedouble>>& W, size_t m);
#endif
