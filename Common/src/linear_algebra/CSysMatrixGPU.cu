/*!
 * \file CSysMatrixGPU.cu
 * \brief Implementations of Kernels and Functions for Matrix Operations on the GPU
 *
 * The kernel implementations will feature a lot of if-statements.
 * The reason for such heavy usage of conditionals is to do a check
 * whether the memory locations being accessed by the threads are
 * within bounds or not. Usually the entire kernel is "wrapped" in
 * a single conditional for these checks. But, in our case, it is
 * necessary for us to use intermittent synchronization barriers like
 * __syncthreads() which will lead to thread divergence issues if used
 * inside a conditional.
 * \author A. Raj
 * \version 8.3.0 "Harrier"
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

#include "../../include/linear_algebra/CSysMatrix.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/linear_algebra/GPUComms.cuh"

using namespace cudaKernelParameters;

template<typename matrixType, typename vectorType>
__device__ void DeviceGaussElimination(matrixType* matrixCopy, vectorType* prod, unsigned long row, unsigned int threadNo, bool rowInPartition, matrixParameters matrixParam)
{
   auto matrixIndex = [=](unsigned short x, unsigned short y){
      return (x * matrixParam.blockColSize + y);
   };

   auto vectorIndex = [=](unsigned short row, unsigned short elemNo){
      return (row * matrixParam.blockRowSize + elemNo);
   };

   unsigned short x = threadNo/matrixParam.blockRowSize;
   unsigned short y = threadNo % matrixParam.blockColSize;

   matrixType pivot, weight = 0.0;

   __syncthreads();

   for(unsigned short currentRow=0; currentRow < matrixParam.blockRowSize; currentRow++)
   {
      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockSize)){
         pivot = matrixCopy[matrixIndex(currentRow, currentRow)];
      }
      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockSize)){

         if(x==currentRow)
         {
            matrixCopy[matrixIndex(currentRow, y)] /= pivot;
            if(y==0) prod[vectorIndex(row, currentRow)] /= pivot;
         }

         weight = matrixCopy[matrixIndex(x, currentRow)];

         if(x!=currentRow)
         {
            matrixCopy[matrixIndex(x, y)] -= weight * matrixCopy[matrixIndex(currentRow, y)];
            if(y==0) prod[vectorIndex(row, x)] -= weight * prod[vectorIndex(row, currentRow)];
         }
      }

   __syncthreads();
   }

}

template<typename matrixType, typename vectorType>
__global__ void FirstSymmetricIterationKernel(matrixType* matrix, vectorType* vec, vectorType* prod, const unsigned long* d_partition_offsets, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, const unsigned long* d_dia_ptr, matrixParameters matrixParam)
{
   /*--- First part of the symmetric iteration: (D+L).x* = b ---*/

   auto matrixIndex = [=](unsigned long row, unsigned short threadNo, const unsigned long* ptr_array){
      return (ptr_array[row] * matrixParam.blockSize + threadNo);
   };

   auto vectorIndex = [=](unsigned long row, unsigned short elemNo){
      return (row * matrixParam.blockRowSize + elemNo);
   };

   auto vectorMultIndex = [=](unsigned short blockNo, unsigned short blockCol){
      return (d_col_ind[blockNo] * matrixParam.blockColSize + blockCol);
   };

   unsigned long origRow = (blockIdx.x * blockDim.x + threadIdx.x)/CUDA_WARP_SIZE;
   unsigned short localRow = origRow % matrixParam.rowsPerBlock;
   unsigned short threadNo = threadIdx.x % CUDA_WARP_SIZE;

   unsigned short blockCol = threadNo % matrixParam.blockColSize;

   extern __shared__ matrixType tempBuffer[];

   /* Forward Substitution Algorithm with Gaussian Elimination. */
   for(auto iPartition = matrixParam.nChainStart; iPartition < matrixParam.nChainEnd; iPartition++)
   {

      /* We move to the set of rows present in the current partition. */
      unsigned long row = origRow + d_partition_offsets[iPartition];
      bool rowInPartition = (row < d_partition_offsets[iPartition + 1]);

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockRowSize)) {

         tempBuffer[matrixParam.shrdMemIndex(localRow, threadNo)] = 0.0;
      }

      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.activeThreads)) {

         for(unsigned long index = matrixIndex(row, threadNo, d_row_ptr); index < matrixIndex(row, 0, d_dia_ptr); index+=matrixParam.activeThreads)
         {

            unsigned short blockNo = index/matrixParam.blockSize;
            unsigned short blockRow = (index/matrixParam.blockColSize) % matrixParam.blockRowSize;

            atomicAdd(&tempBuffer[matrixParam.shrdMemIndex(localRow, blockRow)], matrix[index] * prod[vectorMultIndex(blockNo, blockCol)]);  // Compute L.x*
         }
      }

      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockSize)) {

         if(threadNo < matrixParam.blockRowSize) prod[vectorIndex(row, threadNo)] = vec[vectorIndex(row, threadNo)] - tempBuffer[matrixParam.shrdMemIndex(localRow, threadNo)];  // Compute y = b - L.x*

         /* Reinitialize the shared memory to hold the matrix diagonal elements now. */
         tempBuffer[matrixParam.shrdMemIndex(localRow, threadNo)] = matrix[matrixIndex(row, threadNo, d_dia_ptr)];
      }

      DeviceGaussElimination(&tempBuffer[matrixParam.shrdMemIndex(localRow, 0)], prod, row, threadNo, rowInPartition, matrixParam); // Solve D.x* = y

      __syncthreads();

   }
}


template<typename matrixType, typename vectorType>
__global__ void SecondSymmetricIterationKernel(matrixType* matrix, vectorType* prod, const unsigned long* d_partition_offsets, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, const unsigned long* d_dia_ptr, matrixParameters matrixParam)
{

   /*--- Second part of the symmetric iteration: (D+U).x_(1) = D.x* ---*/

   auto matrixIndex = [=](unsigned long row, unsigned short threadNo, const unsigned long* ptr_array){
      return (ptr_array[row] * matrixParam.blockSize + threadNo);
   };

   auto vectorIndex = [=](unsigned long row, unsigned short elemNo){
      return (row * matrixParam.blockRowSize + elemNo);
   };

   auto vectorMultIndex = [=](unsigned short blockNo, unsigned short blockCol){
      return (d_col_ind[blockNo] * matrixParam.blockColSize + blockCol);
   };

   unsigned long origRow = (blockIdx.x * blockDim.x + threadIdx.x)/CUDA_WARP_SIZE;
   unsigned short localRow = origRow % matrixParam.rowsPerBlock;
   unsigned short threadNo = threadIdx.x % CUDA_WARP_SIZE;

   unsigned short blockCol = threadNo % matrixParam.blockColSize;

   extern __shared__ matrixType tempBuffer[];

   /* Backward Substitution Algorithm with Gaussian Elimination. */
   for(auto iPartition = matrixParam.nChainStart; iPartition > matrixParam.nChainEnd; iPartition--)
   {

      /* We move to the set of rows present in the current partition. */
      unsigned long row = origRow + d_partition_offsets[iPartition - 1];
      bool rowInPartition = (row < d_partition_offsets[iPartition]);

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockRowSize)) {

         tempBuffer[matrixParam.shrdMemIndex(localRow, threadNo)] = 0.0;
      }

      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.activeThreads)) {

         for(unsigned long index = matrixIndex(row, threadNo, d_dia_ptr); index < matrixIndex(row, matrixParam.blockSize, d_dia_ptr); index+=matrixParam.activeThreads)
         {

            unsigned short blockNo = index/matrixParam.blockSize;
            unsigned short blockRow = (index/matrixParam.blockColSize) % matrixParam.blockRowSize;

            atomicAdd(&tempBuffer[matrixParam.shrdMemIndex(localRow, blockRow)], matrix[index] * prod[vectorMultIndex(blockNo, blockCol)]);  // Compute D.x*
         }

         for(unsigned long index = matrixIndex(row, threadNo + matrixParam.blockSize, d_dia_ptr); index < matrixIndex(row + 1, 0, d_row_ptr); index+=matrixParam.activeThreads)
         {

            unsigned short blockNo = index/matrixParam.blockSize;
            unsigned short blockRow = (index/matrixParam.blockColSize) % matrixParam.blockRowSize;

            atomicAdd(&tempBuffer[matrixParam.shrdMemIndex(localRow, blockRow)], -(matrix[index] * prod[vectorMultIndex(blockNo, blockCol)]));  // Compute y = D.x*-U.x_(n+1)
         }
      }

      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockSize)) {

         if(threadNo < matrixParam.blockRowSize) prod[vectorIndex(row, threadNo)] = tempBuffer[matrixParam.shrdMemIndex(localRow, threadNo)];

         /* Reinitialize the memory to hold the matrix diagonal elements now. */
         tempBuffer[matrixParam.shrdMemIndex(localRow, threadNo)] = matrix[matrixIndex(row, threadNo, d_dia_ptr)];
      }

      DeviceGaussElimination(&tempBuffer[matrixParam.shrdMemIndex(localRow, 0)], prod, row, threadNo, rowInPartition, matrixParam); // Solve D.x* = y

      __syncthreads();
   }
}

template<typename matrixType, typename vectorType>
__global__ void MatrixVectorProductKernel(matrixType* matrix, vectorType* vec, vectorType* prod, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, matrixParameters matrixParam)
{
   auto matrixIndex = [=](unsigned long row, unsigned short threadNo) {
      return (d_row_ptr[row] * matrixParam.blockSize + threadNo);
   };

   auto vectorIndex = [=](unsigned long row, unsigned short elemNo) {
      return (row * matrixParam.blockRowSize + elemNo);
   };

   unsigned long row = (blockIdx.x * blockDim.x + threadIdx.x)/CUDA_WARP_SIZE;
   unsigned short threadNo = threadIdx.x % CUDA_WARP_SIZE;
   unsigned short localRow = row % matrixParam.rowsPerBlock;

   unsigned short blockCol = threadNo % matrixParam.blockColSize;

   extern __shared__ vectorType res[];

   if(matrixParam.validAccess(row, threadNo, matrixParam.blockRowSize))
   {
      prod[vectorIndex(row, threadNo)] = 0.0;
      res[vectorIndex(localRow, threadNo)] = 0.0;
   }

   __syncthreads();

   if(matrixParam.validAccess(row, threadNo, matrixParam.activeThreads))
   {

      for(unsigned long index = matrixIndex(row, threadNo); index < matrixIndex(row+1, 0); index+=matrixParam.activeThreads)
      {
         /*Calculating the position of the new element the thread has shifted to.*/
         unsigned short blockNo = index/matrixParam.blockSize;
         unsigned short blockRow = (index/matrixParam.blockColSize) % matrixParam.blockRowSize;

         atomicAdd(&res[vectorIndex(localRow, blockRow)], matrix[index] * vec[(d_col_ind[blockNo]) * matrixParam.blockColSize + blockCol]);
      }

   }

   __syncthreads();

   if(matrixParam.validAccess(row, threadNo, matrixParam.blockRowSize)) prod[vectorIndex(row, threadNo)] = res[vectorIndex(localRow, threadNo)];

}

template<class ScalarType>
void CSysMatrix<ScalarType>::HtDTransfer(bool trigger) const
{
   if(trigger) gpuErrChk(cudaMemcpy((void*)(d_matrix), (void*)&matrix[0], (sizeof(ScalarType)*nnz*nVar*nEqn), cudaMemcpyHostToDevice));
}

template<class ScalarType>
void CSysMatrix<ScalarType>::GPUMatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                 CGeometry* geometry, const CConfig* config) const
                                                 {

   ScalarType* d_vec = vec.GetDevicePointer();
   ScalarType* d_prod = prod.GetDevicePointer();

   HtDTransfer();
   vec.HtDTransfer();
   prod.GPUSetVal(0.0);

   matrixParameters matrixParam(nPointDomain, nEqn, nVar, geometry->nPartition, config->GetRows_Per_Cuda_Block());

  dim3 blockDim(config->GetCuda_Block_Size(),1,1);
  unsigned int gridx = rounded_up_division(config->GetCuda_Block_Size(), matrixParam.totalRows * CUDA_WARP_SIZE);
  dim3 gridDim(gridx, 1, 1);

  MatrixVectorProductKernel<<<gridDim, blockDim, matrixParam.rowsPerBlock * matrixParam.blockRowSize * sizeof(ScalarType)>>>(d_matrix, d_vec, d_prod, d_row_ptr, d_col_ind, matrixParam);
  gpuErrChk( cudaPeekAtLastError() );

  prod.DtHTransfer();

}

template<class ScalarType>
void CSysMatrix<ScalarType>::GPUComputeLU_SGSPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
   CGeometry* geometry, const CConfig* config) const
   {
      ScalarType* d_vec = vec.GetDevicePointer();
      ScalarType* d_prod = prod.GetDevicePointer();

      HtDTransfer();
      vec.HtDTransfer();
      prod.HtDTransfer();

      matrixParameters matrixParam(nPointDomain, nEqn, nVar, geometry->nPartition, config->GetRows_Per_Cuda_Block());

      dim3 blockDim(matrixParam.rowsPerBlock * CUDA_WARP_SIZE,1,1);
      unsigned int gridx = rounded_up_division(matrixParam.rowsPerBlock, geometry->maxPartitionSize);
      dim3 gridDim(gridx, 1, 1);

      for(auto elem = geometry->chainPtr.begin(); elem != geometry->chainPtr.end() - 1; elem++)
      {
         matrixParam.nChainStart = *(elem);
         matrixParam.nChainEnd = *(elem + 1);

         FirstSymmetricIterationKernel<<<gridDim, blockDim, matrixParam.rowsPerBlock * matrixParam.blockSize * sizeof(ScalarType)>>>(d_matrix, d_vec, d_prod, d_partition_offsets, d_row_ptr, d_col_ind, d_dia_ptr, matrixParam);
         gpuErrChk( cudaPeekAtLastError() );
      }

      for(auto elem = geometry->chainPtr.rbegin(); elem != geometry->chainPtr.rend() - 1; elem++)
      {
         matrixParam.nChainStart = *(elem);
         matrixParam.nChainEnd = *(elem + 1);

         SecondSymmetricIterationKernel<<<gridDim, blockDim, matrixParam.rowsPerBlock * matrixParam.blockSize * sizeof(ScalarType)>>>(d_matrix, d_prod, d_partition_offsets, d_row_ptr, d_col_ind, d_dia_ptr, matrixParam);
         gpuErrChk( cudaPeekAtLastError() );
      }

      prod.DtHTransfer();
   }

template class CSysMatrix<su2mixedfloat>; // This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.
