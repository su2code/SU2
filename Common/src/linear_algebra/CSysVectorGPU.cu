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
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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
  if (trigger)
    gpuErrChk(cudaMemcpy((void*)(d_vec_val), (void*)&vec_val[0], (sizeof(ScalarType) * nElm), cudaMemcpyHostToDevice));
}

template <class ScalarType>
void CSysVector<ScalarType>::DtHTransfer(bool trigger) const {
  if (trigger)
    gpuErrChk(cudaMemcpy((void*)(&vec_val[0]), (void*)d_vec_val, (sizeof(ScalarType) * nElm), cudaMemcpyDeviceToHost));
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::GPUDot(const CSysVector& other) const {
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
    SU2_MPI::Error("Unsupported ScalarType in CSysVector::GPUDot.", CURRENT_FUNCTION);
    return ScalarType(0);
  }

  if (status != CUBLAS_STATUS_SUCCESS) {
    SU2_MPI::Error("cuBLAS dot failed in CSysVector::GPUDot.", CURRENT_FUNCTION);
    return ScalarType(0);
  }

  ScalarType global_dot = ScalarType(0);
  const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
  SelectMPIWrapper<ScalarType>::W::Allreduce(&local_dot, &global_dot, 1, mpi_type, MPI_SUM, SU2_MPI::GetComm());

  return global_dot;
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::GPUNorm() const {
  return sqrt(GPUDot(*this));
}

/*!
 * \brief multi vector product CUDA kernel one line of blocks per pair V[i0+i],W[j];
 *        Configurable multiple blocks reducing over the size of the vectors
 */
template <class ScalarType>
__global__ void GPUmultiDot(const ScalarType* const* __restrict__ d_V, const size_t n,
                            const ScalarType* const* __restrict__ d_W, const size_t m, const size_t size,
                            ScalarType* __restrict__ d_local)
{
  // Map each x,y block to the specific (i,j) dot product
  const size_t pair_idx = blockIdx.y;
  if (pair_idx >= n * m) return;

  const size_t i = pair_idx / m;
  const size_t j = pair_idx % m;

  //get the corresponding vectors
  const ScalarType* __restrict__ vi = d_V[i];
  const ScalarType* __restrict__ wj = d_W[j];

  // grid strided loop over the vector elements
  ScalarType local_sum = 0.0;
  const size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  const size_t stride = gridDim.x * blockDim.x;

  for (size_t k = tid; k < size; k += stride)
  {
      local_sum += vi[k] * wj[k];
  }

  // shared memory reduction within the block
  extern __shared__ char shared_mem[];
  ScalarType* sdata = reinterpret_cast<ScalarType*>(shared_mem);

  sdata[threadIdx.x] = local_sum;
  __syncthreads();

  // parallel reduction on the block
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
  {
      if (threadIdx.x < s)
      {
        sdata[threadIdx.x] += sdata[threadIdx.x + s];
      }
      __syncthreads();
  }

  // atomic add of each block partial sum to the output matrix, operated by thread 0 of each block
  if (threadIdx.x == 0)
  {
      atomicAdd(&d_local[i * m + j], sdata[0]);
  }
}

/*!
 * \brief multi vector dot produt method for GPU
 * \note this is a vectors-read only method that returns an array of scalars
 */
template <class ScalarType>
const su2matrix<ScalarType>& CSysVector<ScalarType>::multiDotGPU(const std::vector<CSysVector<ScalarType>>& V,
                                                                 const size_t i0, const size_t n,
                                                                 const std::vector<CSysVector<ScalarType>>& W,
                                                                 const size_t m) {

  static su2matrix<ScalarType> shared;
  if (n == 0 || m == 0) return shared;

  const size_t size = V[0].nElmDomain;

  // get all the device pointers for V and W in one array, resize if needed
  static std::vector<const ScalarType*> h_V_W_ptrs;
  h_V_W_ptrs.resize(n + m);

  for (size_t i = 0; i < n; ++i){
    h_V_W_ptrs[i] = V[i0 + i].GetDevicePointer();
  }
  for (size_t j = 0; j < m; ++j){
    h_V_W_ptrs[j + n] = W[j].GetDevicePointer();
  }

  // persistent device pointer storing all vectors, resizes when needed
  static const ScalarType** d_V_W_ptrs = nullptr;
  static size_t ptrs_capacity = 0; // current capacity
  const size_t ptrs_needed = n + m; // needed capacity for both V and W

  if (ptrs_needed > ptrs_capacity) { // if not enough capacity, enlarge by re-allocation on device
    if (d_V_W_ptrs) gpuErrChk(cudaFree(d_V_W_ptrs));
    gpuErrChk(cudaMalloc(&d_V_W_ptrs, ptrs_needed * sizeof(ScalarType*)));
    ptrs_capacity = ptrs_needed; // update current capacity
  }
  // copy pointers to device
  gpuErrChk(cudaMemcpy(d_V_W_ptrs, h_V_W_ptrs.data(), ptrs_needed * sizeof(ScalarType*), cudaMemcpyHostToDevice));

  // allocate persisten result buffer that grows if needed
  static ScalarType* d_local = nullptr;
  static size_t local_capacity = 0;
  const size_t local_needed = n * m;

  if (local_needed > local_capacity) { // if not enough capacity, enlarge by re-allocation on device
    if (d_local) gpuErrChk(cudaFree(d_local));
    gpuErrChk(cudaMalloc(&d_local, local_needed * sizeof(ScalarType)));
    local_capacity = local_needed;
  }
  // zero out the result buffer
  gpuErrChk(cudaMemset(d_local, 0, local_needed * sizeof(ScalarType)));


  dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE,1,1);
  int numBlocksPerPair = KernelParameters::round_up_division(KernelParameters::MVP_BLOCK_SIZE, size);
  dim3 gridDim(numBlocksPerPair, n * m, 1);

  GPUmultiDot<<<gridDim, blockDim, KernelParameters::MVP_BLOCK_SIZE * sizeof(ScalarType)>>>(&d_V_W_ptrs[0], n, &d_V_W_ptrs[n], m, size, d_local);
  gpuErrChk(cudaGetLastError());

  // copy result to host for MPI reduce
  su2matrix<ScalarType> local(n,m);
  gpuErrChk(cudaMemcpy(local.data(), d_local, n * m * sizeof(ScalarType), cudaMemcpyDeviceToHost));

  /*--- Single AllReduce of the result, only the master thread communicates. ---*/
  // this is a duplicate. Ideally the cuda section should return local but that depends on the intended OpenMP/CUDA combined usage
  SU2_OMP_MASTER {
    shared.resize(n, m);

    const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
    SelectMPIWrapper<ScalarType>::W::Allreduce(local.data(), shared.data(), n * m, mpi_type, MPI_SUM,
                                               SU2_MPI::GetComm());
  }
  END_SU2_OMP_MASTER

  /*--- All threads have the same view of the result. ---*/
  SU2_OMP_BARRIER

  return shared;
}

template<class ScalarType, int N>
struct WeightedVecs {
  const ScalarType* ptrs[N];
  ScalarType weights[N];
};

/*!
 * \brief linear combination kernel to calculate the next vector v from weights and vectors
 */
template<class ScalarType, int N>
__global__ void LinearCombinationKernel(ScalarType* __restrict__ v, WeightedVecs<ScalarType, N> wv,
                                        int n, unsigned long nElm, bool inc)
{
  const unsigned long k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= nElm) return;

  //handle overwriting or combination with existing
  ScalarType result = inc ? v[k] : ScalarType(0);

  #pragma unroll
  for (int i = 0; i < N; ++i) // N is known at compile time (4), this unrolls to: if (i < n) result += weight[i] * vector[i][k]; i<4
    if (i < n) result += wv.weights[i] * wv.ptrs[i][k];
  v[k] = result;
}

/*!
 * \brief dispatcher for the linear combination kernel on GPU
 */
template<class ScalarType>
void CSysVector<ScalarType>::LinearCombinationGPU(const unsigned long n, const std::vector<CSysVector<ScalarType>>& vs, const ScalarType* ws,
                                                  CSysVector<ScalarType>& v, bool inc)
{
  const unsigned long nElm = v.nElmDomain;
  dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE,1,1);
  int numBlocks = KernelParameters::round_up_division(KernelParameters::MVP_BLOCK_SIZE, nElm);
  dim3 gridDim(numBlocks, 1, 1);

  ScalarType* d_v = v.GetDevicePointer();

  for (unsigned long i = 0; i < n; i += 4) {
    const int rem = static_cast<int>(std::min(n - i, 4ul));
    //prepare vectors pointers and corresponding weights, passing them by value
    WeightedVecs<ScalarType, 4> vs_ws = {};
    for (int j = 0; j < rem; ++j) {
      vs_ws.ptrs[j] = vs[i + j].GetDevicePointer();  // get the pointer
      vs_ws.weights[j] = ws[i + j];                  // plain array indexing, not ws(k)
    }
    //calculate the linear combination on GPU, handle more than 4 vectors through inc || i > 0
    LinearCombinationKernel<ScalarType, 4><<<gridDim, blockDim>>>(d_v, vs_ws, rem, nElm, inc || i > 0);
    gpuErrChk(cudaPeekAtLastError());
  }

}

template class CSysVector<su2double>; //This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.
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

#if defined(USE_MIXED_PRECISION)
template void CSysVector<su2mixedfloat>::HtDTransfer(bool trigger) const;
template void CSysVector<su2mixedfloat>::DtHTransfer(bool trigger) const;
template su2mixedfloat CSysVector<su2mixedfloat>::GPUDot(const CSysVector<su2mixedfloat>& other) const;
template su2mixedfloat CSysVector<su2mixedfloat>::GPUNorm() const;
#endif

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template void CSysVector<passivedouble>::HtDTransfer(bool trigger) const;
template void CSysVector<passivedouble>::DtHTransfer(bool trigger) const;
template passivedouble CSysVector<passivedouble>::GPUDot(const CSysVector<passivedouble>& other) const;
template passivedouble CSysVector<passivedouble>::GPUNorm() const;
#endif
