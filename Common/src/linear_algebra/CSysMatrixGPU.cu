/*!
 * \file CSysMatrixGPU.cu
 * \brief Implementations of Kernels and Functions for Matrix Operations on the GPU
 * \author A. Raj
 * \version 8.2.0 "Harrier"
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

using namespace kernelParameters;


template<typename matrixType, typename vectorType>
__device__ void GaussEliminationKernel(matrixType* matrixCopy, vectorType* prod, unsigned int threadNo, matrixParameters matrixParam)
{
   #define A(I, J) matrixCopy[(I) * matrixParam.blockColSize + (J)]

   int x = threadNo/matrixParam.blockRowSize;                  
   int y = threadNo % matrixParam.blockColSize;

   matrixType pivot, weight = 0.0;
   
   __syncthreads();

   for(int currentRow=0; currentRow < matrixParam.blockRowSize; currentRow++)
   {
      pivot = A(currentRow, currentRow);

      __syncthreads();

      A(currentRow, y) /= pivot;
      if(y==0) prod[currentRow] /= pivot;

      weight = A(x, currentRow);

      if(x!=currentRow)
      {
         A(x, y) -= weight * A(currentRow, y);
         if(y==0) prod[x] -= weight * prod[currentRow];
      }

   __syncthreads();
   }

   #undef A
}

template<typename matrixType, typename vectorType>
__global__ void FirstSymmetricIterationKernel(matrixType* matrix, vectorType* vec, vectorType* prod, const unsigned long* d_partition_offsets, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, const unsigned long* d_dia_ptr, matrixParameters matrixParam)
{

   auto matrixIndex = [=](unsigned long row, unsigned short threadNo){
      return (d_row_ptr[row] * matrixParam.blockSize + threadNo);
   };

   auto matrixDiagonalIndex = [=](unsigned long row, unsigned short threadNo){
      return (d_row_ptr[row] * matrixParam.blockSize + threadNo);
   };

   auto vectorIndex = [=](unsigned long row, unsigned short elemNo){
      return (row * matrixParam.blockRowSize + elemNo);
   };

   unsigned long row = (blockIdx.x * blockDim.x + threadIdx.x)/MVP_WARP_SIZE; 
   unsigned short localRow = row % matrixParam.rowsPerBlock;
   unsigned short threadNo = threadIdx.x % MVP_WARP_SIZE;

   unsigned short blockCol = threadNo % matrixParam.blockColSize;

   extern __shared__ vectorType tempBuffer[];

   for(auto iPartition = 1ul; iPartition < matrixParam.nPartition; iPartition++)
   {  

      /*We move to the set of rows present in the current partition.*/
      row += d_partition_offsets[iPartition-1]; 
      
      bool rowInPartition = (row < d_partition_offsets[iPartition]);

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockRowSize)) {

         prod[vectorIndex(row, threadNo)] = 0.0;
         tempBuffer[vectorIndex(localRow, threadNo)] = 0.0;
      }
      
      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.activeThreads)) {

         for(unsigned long index = matrixIndex(row, threadNo); index < matrixDiagonalIndex(row, 0); index+=matrixParam.activeThreads)
         {

            unsigned short blockNo = index/matrixParam.blockSize;
            unsigned short blockRow = (index/matrixParam.blockColSize) % matrixParam.blockRowSize;    

            atomicAdd(&tempBuffer[vectorIndex(localRow, blockRow)], matrix[index] * vec[(d_col_ind[blockNo]) * matrixParam.blockColSize + blockCol]);
         }
      }

      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockRowSize)) {

         prod[vectorIndex(row, threadNo)] = vec[vectorIndex(row, threadNo)] - tempBuffer[vectorIndex(localRow, threadNo)];
      } 
     
      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockSize)) {

         tempBuffer[vectorIndex(localRow, threadNo)] = matrix[matrixDiagonalIndex(row, threadNo)];
         GaussEliminationKernel(&tempBuffer[vectorIndex(localRow, threadNo)], &prod[vectorIndex(row, threadNo/matrixParam.blockColSize)], threadNo, matrixParam);
      } 
      
      __syncthreads();

   } 
}


template<typename matrixType, typename vectorType>
__global__ void SecondSymmetricIterationKernel(matrixType* matrix, vectorType* prod, const unsigned long* d_partition_offsets, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, const unsigned long* d_dia_ptr, matrixParameters matrixParam)
{

   auto matrixIndex = [=](unsigned long row, unsigned short threadNo){
      return (d_row_ptr[row] * matrixParam.blockSize + threadNo);
   };

   auto matrixDiagonalIndex = [=](unsigned long row, unsigned short threadNo){
      return (d_dia_ptr[row] * matrixParam.blockSize + threadNo);
   };

   auto vectorIndex = [=](unsigned long row, unsigned short elemNo){
      return (row * matrixParam.blockRowSize + elemNo);
   };

   unsigned long row = (blockIdx.x * blockDim.x + threadIdx.x)/MVP_WARP_SIZE; 
   unsigned short localRow = row % matrixParam.rowsPerBlock;
   unsigned short threadNo = threadIdx.x % MVP_WARP_SIZE;

   unsigned short blockCol = threadNo % matrixParam.blockColSize;

   extern __shared__ vectorType tempBuffer[];

   for(auto iPartition = 1ul; iPartition < matrixParam.nPartition; iPartition++)
   {  

      /*We move to the set of rows present in the current partition.*/
      row += d_partition_offsets[iPartition-1]; 
      
      bool rowInPartition = (row < d_partition_offsets[iPartition]);

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockRowSize)) {

         tempBuffer[vectorIndex(localRow, threadNo)] = 0.0;
      }
      
      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.activeThreads)) {

         for(unsigned long index = matrixDiagonalIndex(row, threadNo); index < matrixDiagonalIndex(row, 0) + matrixParam.blockSize; index+=matrixParam.activeThreads)
         {

            unsigned short blockNo = index/matrixParam.blockSize;
            unsigned short blockRow = (index/matrixParam.blockColSize) % matrixParam.blockRowSize;    

            atomicAdd(&tempBuffer[vectorIndex(localRow, blockRow)], matrix[index] * prod[(d_col_ind[blockNo]) * matrixParam.blockColSize + blockCol]);
         }
      }

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.activeThreads)) {

         for(unsigned long index = matrixDiagonalIndex(row, threadNo); index < matrixIndex(row+1, 0); index+=matrixParam.activeThreads)
         {

            unsigned short blockNo = index/matrixParam.blockSize;
            unsigned short blockRow = (index/matrixParam.blockColSize) % matrixParam.blockRowSize;    

            atomicAdd(&tempBuffer[vectorIndex(localRow, blockRow)], -(matrix[index] * prod[(d_col_ind[blockNo]) * matrixParam.blockColSize + blockCol]));
         }
      }

      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockRowSize)) {

         prod[vectorIndex(row, threadNo)] = tempBuffer[vectorIndex(localRow, threadNo)];
      } 
     
      __syncthreads();

      if(matrixParam.validParallelAccess(rowInPartition, threadNo, matrixParam.blockSize)) {

         tempBuffer[vectorIndex(localRow, threadNo)] = matrix[(d_dia_ptr[row] * matrixParam.blockSize) + threadNo];
         GaussEliminationKernel(&tempBuffer[vectorIndex(localRow, threadNo)], &prod[vectorIndex(row, threadNo/matrixParam.blockColSize)], threadNo, matrixParam);
      } 
      
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

   unsigned long row = (blockIdx.x * blockDim.x + threadIdx.x)/MVP_WARP_SIZE;
   unsigned short threadNo = threadIdx.x % MVP_WARP_SIZE;
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

   matrixParameters matrixParam(nPointDomain, nEqn, nVar, geometry->nPartition);

  dim3 blockDim(MVP_BLOCK_SIZE,1,1);
  unsigned int gridx = rounded_up_division(MVP_BLOCK_SIZE, matrixParam.totalRows * MVP_WARP_SIZE);
  dim3 gridDim(gridx, 1, 1);

  MatrixVectorProductKernel<<<gridDim, blockDim>>>(d_matrix, d_vec, d_prod, d_row_ptr, d_col_ind, matrixParam);
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

      matrixParameters matrixParam(nPointDomain, nEqn, nVar, geometry->nPartition);

      dim3 blockDim(matrixParam.rowsPerBlock * MVP_WARP_SIZE,1,1);
      unsigned int gridx = rounded_up_division(matrixParam.rowsPerBlock, geometry->maxPartitionSize);
      dim3 gridDim(gridx, 1, 1);

      FirstSymmetricIterationKernel<<<gridDim, blockDim, matrixParam.rowsPerBlock * matrixParam.blockSize * sizeof(ScalarType)>>>(d_matrix, d_vec, d_prod, d_partition_offsets, d_row_ptr, d_col_ind, d_dia_ptr, matrixParam);
      gpuErrChk( cudaPeekAtLastError() );
      SecondSymmetricIterationKernel<<<gridDim, blockDim, matrixParam.rowsPerBlock * matrixParam.blockSize * sizeof(ScalarType)>>>(d_matrix, d_prod, d_partition_offsets, d_row_ptr, d_col_ind, d_dia_ptr, matrixParam);
      gpuErrChk( cudaPeekAtLastError() );
      prod.DtHTransfer();
   }

template class CSysMatrix<su2mixedfloat>; //This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.
