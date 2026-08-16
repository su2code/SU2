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
#include <algorithm>
#include <vector>

namespace {

/*!
 * \brief Fixed launch shape for DotKernel, so the number of warps per block (and thus the static
 *        shared-memory reduction buffer's size) is a compile-time constant - no dynamic shared
 *        memory needed for a single running sum, unlike MultiDotKernel's per-(i,j) buffer.
 */
constexpr unsigned int DOT_THREADS_PER_BLOCK = 256;
constexpr unsigned int DOT_WARPS_PER_BLOCK = DOT_THREADS_PER_BLOCK / 32u;

/*!
 * \brief Single dot product result[0] = <x, y>, same single-pass grid-stride + warp-shuffle +
 *        block-combine + atomicAdd reduction as MultiDotKernel (see its comment), specialized for
 *        the one-running-sum case: no per-thread array, no dynamic shared memory, just a plain
 *        register accumulator and a small static warpSums buffer.
 */
template <class ScalarType>
__global__ void DotKernel(const ScalarType* __restrict__ x, const ScalarType* __restrict__ y, unsigned long size,
                          ScalarType* __restrict__ result) {
  ScalarType sum = 0;
  const auto stride = static_cast<unsigned long>(blockDim.x) * gridDim.x;
  for (auto k = static_cast<unsigned long>(blockIdx.x) * blockDim.x + threadIdx.x; k < size; k += stride) {
    sum += x[k] * y[k];
  }

  /*--- Warp-level reduction: after this, lane 0 of every warp holds that warp's true sum. ---*/
  for (int offset = 16; offset > 0; offset >>= 1) sum += __shfl_down_sync(0xFFFFFFFFu, sum, offset);

  /*--- Block-level combine: warp leaders stage their partials in static shared memory, thread 0
   * sums them and issues the one atomicAdd this block contributes to the global result
   * (pre-zeroed by the caller). ---*/
  __shared__ ScalarType warpSums[DOT_WARPS_PER_BLOCK];
  const unsigned int lane = threadIdx.x % 32u;
  const unsigned int warpId = threadIdx.x / 32u;

  if (lane == 0) warpSums[warpId] = sum;
  __syncthreads();

  if (threadIdx.x == 0) {
    ScalarType blockSum = 0;
    for (unsigned int w = 0; w < DOT_WARPS_PER_BLOCK; ++w) blockSum += warpSums[w];
    atomicAdd(result, blockSum);
  }
}

/*!
 * \brief Cap on n*m (the number of dot products one MultiDotKernel launch computes), sizing the
 *        per-thread accumulator array below and the block's shared-memory reduction buffer -
 *        both compile-time sized, since CUDA has no runtime-sized register/local arrays.
 *        FGCRODR's Ritz-value path (CSysSolve.cpp) can reach n = m+1 with m =
 *        LINEAR_SOLVER_RESTART_DEFLATION (user-configurable, default 4 but not uncommon to raise
 *        into the tens), so n*m grows quadratically with that setting - e.g. m=10 already needs
 *        110. Sized generously above realistic usage; multiDotGPU raises a clear SU2_MPI::Error
 *        rather than silently truncating if a caller ever needs more.
 */
constexpr unsigned int MULTIDOT_MAX_NM = 1024;

/*!
 * \brief Caps on the individual vector counts n and m, bounding the fixed-size pointer arrays
 *        passed into MultiDotKernel by value as ordinary launch parameters (see
 *        MultiDotPointers) - CUDA marshals those for you as part of the launch itself, regardless
 *        of host memory pinning.
 *        Asymmetric because the two call sites in this codebase produce genuinely asymmetric
 *        shapes: ModGramSchmidt (CSysSolve.cpp) grows m up to the Krylov restart length with n=1
 *        (a row vector, not really a rectangle), while FGCRODR's Ritz-value path grows n up to
 *        the restart length with m bounded by the usually much smaller
 *        LINEAR_SOLVER_RESTART_DEFLATION - a genuine, independently-sized rectangle. multiDotGPU
 *        always puts whichever of n, m is larger into the "large" kernel argument (swapping V/W
 *        and transposing the result back if needed), so one axis only needs to cover the
 *        deflation-count case while the other covers the restart-length case.
 * \note Unrelated to MULTIDOT_MAX_NM above: n and m can each independently reach this cap without
 *       n*m approaching MULTIDOT_MAX_VEC_LARGE * MULTIDOT_MAX_VEC_SMALL, since multiDotGPU checks
 *       n*m against MULTIDOT_MAX_NM separately (that one is a genuine per-thread storage cost
 *       paid regardless of the runtime n*m; these two only bound a launch parameter's marshaling
 *       size).
 */
constexpr unsigned int MULTIDOT_MAX_VEC_LARGE = 256;
constexpr unsigned int MULTIDOT_MAX_VEC_SMALL = 64;

/*!
 * \brief Fixed-size array of device pointers, passed to MultiDotKernel by value (see
 *        MULTIDOT_MAX_VEC_LARGE/SMALL for why).
 */
template <class ScalarType, unsigned int MaxCount>
struct MultiDotPointers {
  const ScalarType* ptr[MaxCount];
};

/*!
 * \brief Compute the aCount*bCount matrix of dot products D(a,b) = <A[a], B[b]> in a single pass
 *        over the data: every thread reads each of the aCount+bCount vectors once per k and
 *        forms all aCount*bCount products from that same read, so the vectors are only ever read
 *        once total (aCount+bCount reads of length size, not aCount*bCount), and the reduction
 *        stays memory-bound.
 * \note Each thread accumulates its own private running aCount*bCount sums while striding over k
 *       (kept in thread-local storage, capped at MULTIDOT_MAX_NM - small enough that it should
 *       stay resident in registers or L1 for realistic shapes, cheap either way next to the
 *       K-length main loop's DRAM traffic). Only after that loop do threads combine: a
 *       warp-shuffle reduction, then one shared-memory combine and one atomicAdd per (a,b) per
 *       block, so total atomics are O(aCount*bCount * numBlocks), not O(aCount*bCount * size).
 * \note A is the "large" argument (up to MULTIDOT_MAX_VEC_LARGE), B the "small" one (up to
 *       MULTIDOT_MAX_VEC_SMALL) - the caller (multiDotGPU) is responsible for putting the larger
 *       of its two vector counts into A, and transposing the result back if it had to swap.
 */
template <class ScalarType>
__global__ void MultiDotKernel(MultiDotPointers<ScalarType, MULTIDOT_MAX_VEC_LARGE> A, unsigned int aCount,
                               MultiDotPointers<ScalarType, MULTIDOT_MAX_VEC_SMALL> B, unsigned int bCount,
                               unsigned long size, ScalarType* __restrict__ D) {
  const unsigned int nm = aCount * bCount;

  ScalarType local[MULTIDOT_MAX_NM];
  for (unsigned int t = 0; t < nm; ++t) local[t] = ScalarType(0);

  const auto stride = static_cast<unsigned long>(blockDim.x) * gridDim.x;
  for (auto k = static_cast<unsigned long>(blockIdx.x) * blockDim.x + threadIdx.x; k < size; k += stride) {
    for (unsigned int a = 0; a < aCount; ++a) {
      const ScalarType v = A.ptr[a][k];
      for (unsigned int b = 0; b < bCount; ++b) local[a * bCount + b] += v * B.ptr[b][k];
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
      atomicAdd(&D[t], sum);
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
  /*--- Both operands are already on the device, the caller owns the transfers. This reduces over
   * MPI, so it must be called by a single thread (see SU2_DEVICE_REGION).
   * \note Runs on the default stream: x/y may have just been written by
   *       VecExpr::AssignDeviceExpression, which launches on the default stream with no
   *       synchronization of its own, so a dedicated stream here would have no guaranteed
   *       ordering against that write - CUDA only orders operations within the same stream. This
   *       function also always synchronizes before returning (it needs a real value for the MPI
   *       reduction), so a dedicated stream would not buy any overlap either. ---*/
  static ScalarType* d_result = nullptr;
  if (d_result == nullptr) gpuErrChk(cudaMalloc(&d_result, sizeof(ScalarType)));

  gpuErrChk(cudaMemsetAsync(d_result, 0, sizeof(ScalarType)));

  const auto blocks = static_cast<unsigned int>(
      std::min<unsigned long>((nElmDomain + DOT_THREADS_PER_BLOCK - 1) / DOT_THREADS_PER_BLOCK, 1024));
  DotKernel<ScalarType><<<blocks, DOT_THREADS_PER_BLOCK>>>(GetDevicePointer(), other.GetDevicePointer(), nElmDomain,
                                                           d_result);

  ScalarType local_dot = ScalarType(0);
  gpuErrChk(cudaMemcpyAsync(&local_dot, d_result, sizeof(ScalarType), cudaMemcpyDeviceToHost));
  gpuErrChk(cudaStreamSynchronize(nullptr));
  gpuErrChk(cudaGetLastError());

  ScalarType global_dot = ScalarType(0);
  const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
  SelectMPIWrapper<ScalarType>::W::Allreduce(&local_dot, &global_dot, 1, mpi_type, MPI_SUM, SU2_MPI::GetComm());

  return global_dot;
}

/*!
 * \brief Multi vector dot product, local(i,j) = <V[i0+i], W[j]>, via MultiDotKernel.
 * \note Whichever of n, m is larger is passed as MultiDotKernel's "A" (large-capacity) argument
 *       and the other as "B" (small-capacity) - see MULTIDOT_MAX_VEC_LARGE/SMALL for why - so if
 *       m > n, V and W (and their counts) are swapped for the call, and the resulting m*n matrix
 *       is transposed back into the n*m shape the caller expects (cheap: this matrix is at most
 *       MULTIDOT_MAX_NM elements, nowhere near the size of the reduction itself).
 */
template <class ScalarType>
su2matrix<ScalarType> CSysVector<ScalarType>::multiDotGPU(const std::vector<CSysVector<ScalarType>>& V, const size_t i0,
                                                          const size_t n, const std::vector<CSysVector<ScalarType>>& W,
                                                          const size_t m) {
  const size_t nm = n * m;
  if (nm > MULTIDOT_MAX_NM) {
    SU2_MPI::Error("CSysVector::multiDotGPU: n*m exceeds MULTIDOT_MAX_NM, raise that constant in "
                   "CSysVectorGPU.cu.",
                   CURRENT_FUNCTION);
  }

  su2matrix<ScalarType> local;
  local.resize(n, m);
  if (nm == 0) return local;

  const bool swap = m > n;
  const size_t aCount = swap ? m : n;
  const size_t bCount = swap ? n : m;
  if (aCount > MULTIDOT_MAX_VEC_LARGE || bCount > MULTIDOT_MAX_VEC_SMALL) {
    SU2_MPI::Error("CSysVector::multiDotGPU: n or m exceeds MULTIDOT_MAX_VEC_LARGE/SMALL, raise "
                   "those constants in CSysVectorGPU.cu.",
                   CURRENT_FUNCTION);
  }

  const size_t size = V[0].nElmDomain;

  /*--- Persistent device workspace for the output only, cached across calls and freed
   * automatically when the program exits (static local destruction), instead of leaking. The
   * V/W pointers themselves need no device buffer at all: they go to MultiDotKernel as ordinary
   * by-value launch parameters (see MultiDotPointers/MULTIDOT_MAX_VEC_LARGE/SMALL).
   * \note Runs on the default stream, same reasoning as dotGPU (see its comment): V/W may have
   *       just been written by VecExpr::AssignDeviceExpression (e.g. LinearCombination building
   *       w[i+1] right before ModGramSchmidt's multiDot call on it), which also runs on the
   *       default stream with no synchronization of its own. ---*/
  struct Workspace {
    ScalarType* d_D = nullptr;
    size_t capacity = 0;

    void EnsureCapacity(size_t nm) {
      if (nm > capacity) {
        cudaFree(d_D);
        gpuErrChk(cudaMalloc(&d_D, nm * sizeof(ScalarType)));
        capacity = nm;
      }
    }

    ~Workspace() { cudaFree(d_D); }
  };
  static Workspace ws;
  ws.EnsureCapacity(nm);

  MultiDotPointers<ScalarType, MULTIDOT_MAX_VEC_LARGE> aPtrs{};
  MultiDotPointers<ScalarType, MULTIDOT_MAX_VEC_SMALL> bPtrs{};
  if (!swap) {
    for (size_t i = 0; i < n; ++i) aPtrs.ptr[i] = V[i0 + i].GetDevicePointer();
    for (size_t j = 0; j < m; ++j) bPtrs.ptr[j] = W[j].GetDevicePointer();
  } else {
    for (size_t i = 0; i < m; ++i) aPtrs.ptr[i] = W[i].GetDevicePointer();
    for (size_t j = 0; j < n; ++j) bPtrs.ptr[j] = V[i0 + j].GetDevicePointer();
  }

  gpuErrChk(cudaMemsetAsync(ws.d_D, 0, nm * sizeof(ScalarType)));

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

  MultiDotKernel<ScalarType><<<blocks, threadsPerBlock, sharedBytes>>>(
      aPtrs, static_cast<unsigned int>(aCount), bPtrs, static_cast<unsigned int>(bCount), size, ws.d_D);

  if (!swap) {
    gpuErrChk(cudaMemcpyAsync(local.data(), ws.d_D, nm * sizeof(ScalarType), cudaMemcpyDeviceToHost));
    gpuErrChk(cudaStreamSynchronize(nullptr));
  } else {
    /*--- D is aCount*bCount = m*n row-major (D(a,b) = <W[a],V[b]>); local is n*m with
     * local(i,j) = <V[i],W[j]> = D(j,i), i.e. local is D transposed. ---*/
    static std::vector<ScalarType> D;
    D.resize(nm);
    gpuErrChk(cudaMemcpyAsync(D.data(), ws.d_D, nm * sizeof(ScalarType), cudaMemcpyDeviceToHost));
    gpuErrChk(cudaStreamSynchronize(nullptr));
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < m; ++j) local(i, j) = D[j * n + i];
  }
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
