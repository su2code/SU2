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

/*!
 * \brief block-level reduction of elementwise products, accumulated into a single device scalar.
 */
template <class ScalarType>
__global__ void GPUDotKernel(const ScalarType* __restrict__ a, const ScalarType* __restrict__ b,
                             unsigned long n, ScalarType* __restrict__ result) {
    //shared memory
    extern __shared__ unsigned char smem_raw[];
    ScalarType* sdata = reinterpret_cast<ScalarType*>(smem_raw);

    //local and global thread indexes for access and block reduction
    const unsigned long tid = threadIdx.x;
    unsigned long idx = blockIdx.x * blockDim.x + tid;
    const unsigned long stride = blockDim.x * gridDim.x;

    //thread reduction
    ScalarType local = ScalarType(0);
    for (; idx < n; idx += stride) local += a[idx] * b[idx];

    sdata[tid] = local;
    __syncthreads();

    //block reduction
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }
    
    //final atomic add per block
    if (tid == 0) atomicAdd(result, sdata[0]);
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::GPUDot(const CSysVector& other) const {

    dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE, 1, 1);
    int numBlocks = KernelParameters::round_up_division(KernelParameters::MVP_BLOCK_SIZE, this->nElmDomain);
    dim3 gridDim(numBlocks, 1, 1);

    // allocate and zero the result scalar
    ScalarType* d_dot_result;
    gpuErrChk(cudaMalloc(&d_dot_result, sizeof(ScalarType)));
    gpuErrChk(cudaMemset(d_dot_result, 0, sizeof(ScalarType)));
  
    const size_t sharedBytes = KernelParameters::MVP_BLOCK_SIZE * sizeof(ScalarType);
    GPUDotKernel<<<gridDim, blockDim, sharedBytes>>>(this->d_vec_val, other.d_vec_val, this->nElmDomain, d_dot_result);
    gpuErrChk(cudaPeekAtLastError());

    ScalarType result;
    gpuErrChk(cudaMemcpy(&result, d_dot_result, sizeof(ScalarType), cudaMemcpyDeviceToHost));
    gpuErrChk(cudaFree(d_dot_result));



    return result;
}

/*!
 * \brief namespace defining the scalar operators for GPU
 */
namespace {
template <class T> struct OpSetScalar { T val; __device__ __forceinline__ T operator()(T x) const { return val; } };
template <class T> struct OpAddScalar { T val; __device__ __forceinline__ T operator()(T x) const { return x + val; } };
template <class T> struct OpSubScalar { T val; __device__ __forceinline__ T operator()(T x) const { return x - val; } };
template <class T> struct OpMulScalar { T val; __device__ __forceinline__ T operator()(T x) const { return x * val; } };
template <class T> struct OpDivScalar { T val; __device__ __forceinline__ T operator()(T x) const { return x / val; } };
}  // namespace

/*!
 * \brief generic Unary Operation Kernel to apply given functors on GPU
 */
template <class ScalarType, class Operator>
__global__ void GPUUnaryOperationKernel(ScalarType* __restrict__ vec, unsigned long n, Operator op) {
    const unsigned long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) vec[idx] = op(vec[idx]);
}

template <class ScalarType>
void CSysVector<ScalarType>::GPUUnaryOperation(GPUScalarOp op, ScalarType val) {

    dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE, 1, 1);
    int numBlocks = KernelParameters::round_up_division(KernelParameters::MVP_BLOCK_SIZE, this->nElmDomain);
    dim3 gridDim(numBlocks, 1, 1);

    switch (op) {
      case GPUScalarOp::SET:
        GPUUnaryOperationKernel<<<gridDim, blockDim>>>(this->d_vec_val, this->nElmDomain, OpSetScalar<ScalarType>{val});
        break;
      case GPUScalarOp::ADD:
        GPUUnaryOperationKernel<<<gridDim, blockDim>>>(this->d_vec_val, this->nElmDomain, OpAddScalar<ScalarType>{val});
        break;
      case GPUScalarOp::SUB:
        GPUUnaryOperationKernel<<<gridDim, blockDim>>>(this->d_vec_val, this->nElmDomain, OpSubScalar<ScalarType>{val});
        break;
      case GPUScalarOp::MUL:
        GPUUnaryOperationKernel<<<gridDim, blockDim>>>(this->d_vec_val, this->nElmDomain, OpMulScalar<ScalarType>{val});
        break;
      case GPUScalarOp::DIV:
        GPUUnaryOperationKernel<<<gridDim, blockDim>>>(this->d_vec_val, this->nElmDomain, OpDivScalar<ScalarType>{val});
        break;
    }
    gpuErrChk(cudaPeekAtLastError());  
    gpuErrChk(cudaDeviceSynchronize());


}

template<class ScalarType>
void CSysVector<ScalarType>::HtDTransfer(bool trigger) const
{
   if(trigger) gpuErrChk(cudaMemcpy((void*)(d_vec_val), (void*)&vec_val[0], (sizeof(ScalarType)*nElm), cudaMemcpyHostToDevice));
}

template<class ScalarType>
void CSysVector<ScalarType>::DtHTransfer(bool trigger) const
{
   if(trigger) gpuErrChk(cudaMemcpy((void*)(&vec_val[0]), (void*)d_vec_val, (sizeof(ScalarType)*nElm), cudaMemcpyDeviceToHost));
}

template<class ScalarType>
void CSysVector<ScalarType>::GPUSetVal(ScalarType val, bool trigger) const
{
   if(trigger) gpuErrChk(cudaMemset((void*)(d_vec_val), val, (sizeof(ScalarType)*nElm)));
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

template <class ScalarType>
const su2matrix<ScalarType>& CSysVector<ScalarType>::multiDotGPU(const std::vector<CSysVector<ScalarType>>& V,
                                                                 const size_t i0, const size_t n,
                                                                 const std::vector<CSysVector<ScalarType>>& W,
                                                                 const size_t m) {
  SU2_ZONE_SCOPED

  static su2matrix<ScalarType> shared;
  if (n == 0 || m == 0) return shared;

  const size_t size = V[0].nElmDomain;

  // get all the device pointers
  std::vector<const ScalarType*> h_V_ptrs(n), h_W_ptrs(m);
  for (size_t i = 0; i < n; ++i){
    h_V_ptrs[i] = V[i0 + i].data();
  }
  for (size_t j = 0; j < m; ++j){
    //gpuErrChk(cudaMemAdvise(vec_val, nElm * sizeof(ScalarType), cudaMemAdviseSetReadMostly, device_id)); //if read only, could be good
    h_W_ptrs[j] = W[j].data();
  }

  // copy the pointers to the device arrays of pointers
  const ScalarType** d_V_ptrs;
  const ScalarType** d_W_ptrs;
  gpuErrChk(cudaMalloc(&d_V_ptrs, n * sizeof(ScalarType*)));
  gpuErrChk(cudaMalloc(&d_W_ptrs, m * sizeof(ScalarType*)));
  gpuErrChk(cudaMemcpy(d_V_ptrs, h_V_ptrs.data(), n * sizeof(ScalarType*), cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy(d_W_ptrs, h_W_ptrs.data(), m * sizeof(ScalarType*), cudaMemcpyHostToDevice));

  // allocate result buffer, zero it
  ScalarType* d_local;
  gpuErrChk(cudaMalloc(&d_local, n * m * sizeof(ScalarType)));
  gpuErrChk(cudaMemset(d_local, 0, n * m * sizeof(ScalarType)));

  dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE,1,1);
  int numBlocksPerPair = KernelParameters::round_up_division(KernelParameters::MVP_BLOCK_SIZE, size);
  dim3 gridDim(numBlocksPerPair, n * m, 1);

  GPUmultiDot<<<gridDim, blockDim, KernelParameters::MVP_BLOCK_SIZE * sizeof(ScalarType)>>>(d_V_ptrs, n, d_W_ptrs, m, size, d_local);
  gpuErrChk(cudaGetLastError());

  // copy result to host for MPI reduce
  su2matrix<ScalarType> local(n,m);
  gpuErrChk(cudaMemcpy(local.data(), d_local, n * m * sizeof(ScalarType), cudaMemcpyDeviceToHost));

  /*--- Single AllReduce of the result, only the master thread communicates. ---*/
  SU2_OMP_MASTER {
    shared.resize(n, m);

    const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
    SelectMPIWrapper<ScalarType>::W::Allreduce(local.data(), shared.data(), n * m, mpi_type, MPI_SUM,
                                               SU2_MPI::GetComm());
  }
  END_SU2_OMP_MASTER

  /*--- All threads have the same view of the result. ---*/
  SU2_OMP_BARRIER


  // clean allocations
  gpuErrChk(cudaFree(d_local));
  gpuErrChk(cudaFree(d_V_ptrs));
  gpuErrChk(cudaFree(d_W_ptrs));

  return shared;
}

template<class ScalarType, int N>
struct WeightedVecs {
  const ScalarType* ptrs[N];
  ScalarType weights[N];
};

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

template<class ScalarType>
void CSysVector<ScalarType>::LinearCombinationGPU(const unsigned long n, const std::vector<CSysVector<ScalarType>>& vs, const ScalarType* ws,
                                                  CSysVector<ScalarType>& v, bool inc)
{

  const unsigned long nElm = v.nElmDomain;
  dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE,1,1);
  int numBlocks = KernelParameters::round_up_division(KernelParameters::MVP_BLOCK_SIZE, nElm);
  dim3 gridDim(numBlocks, 1, 1);

  for (unsigned long i = 0; i < n; i += 4) {
    const int rem = static_cast<int>(std::min(n - i, 4ul));
    //prepare vectors pointers and corresponding weights, passing them by value
    WeightedVecs<ScalarType, 4> vs_ws = {};
    for (int j = 0; j < rem; ++j) {
      vs_ws.ptrs[j] = vs[i + j].data();  // already on device from multiDot
      vs_ws.weights[j] = ws[i + j];      // plain array indexing, not ws(k)
    }
    //calculate the linear combination on GPU, handle more than 4 vectors through inc || i > 0
    LinearCombinationKernel<ScalarType, 4><<<gridDim, blockDim>>>(v.data(), vs_ws, rem, nElm, inc || i > 0);
    gpuErrChk(cudaPeekAtLastError());
  }

  gpuErrChk(cudaDeviceSynchronize());

}

template class CSysVector<su2double>; //This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.
