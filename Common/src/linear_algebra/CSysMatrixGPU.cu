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
__global__ void GaussEliminationKernel(matrixType* matrix, vectorType* prod, const unsigned long* d_dia_ptr, matrixParameters matrixParam)
{
   #define A(I, J) matrixCopy[(I) * matrixParam.blockColSize + (J)]

   int x = threadIdx.x/matrixParam.blockRowSize;                                         //  A single block is in use so no remainder step necessary like in MVP
   int y = threadIdx.x % matrixParam.blockColSize;

   matrixType pivot, weight = 0.0;

  extern __shared__ matrixType matrixCopy[];

   A(x, y) = matrix[x * matrixParam.blockColSize + y];

   __syncthreads();

   for(int currentRow=0; currentRow< matrixParam.blockRowSize; currentRow++)
   {
      pivot = A(currentRow, currentRow);

      __syncthreads();

      A(currentRow, y) /= pivot;
      if(y==0) prod[currentRow] /= pivot;

      __syncthreads();

      weight = A(x, currentRow);

      if(x!=currentRow)
      {
         A(x, y) -= weight * A(currentRow, y);
         if(y==0) prod[x] -= weight * prod[currentRow];
      }

   __syncthreads();
   }

   #undef A
   #undef validRow
}

template<typename matrixType, typename vectorType>
__global__ void FirstSymmetricIterationKernel(matrixType* matrix, vectorType* vec, vectorType* prod, unsigned long* d_partition_offsets, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, const unsigned long* d_dia_ptr, matrixParameters matrixParam, precondParameters precondParam)
{
   for(auto iPartition = 1ul; iPartition < matrixParam.nPartition - 1; iPartition++)
   {  
      #define matrixIndex(row, threadNo) d_row_ptr[row] * matrixParam.blockSize + threadNo
      #define vectorIndex(row, elemNo) row * matrixParam.blockColSize + elemNo

      int row = (blockIdx.x * blockDim.x + threadIdx.x)/MVP_WARP_SIZE + d_partition_offsets[iPartition-1];
      bool rowInLevel = (row < d_partition_offsets[iPartition]);
   
      int threadNo = threadIdx.x % MVP_WARP_SIZE;
      int activeThreads = matrixParam.blockColSize * (MVP_WARP_SIZE/matrixParam.blockRowSize);

      int blockRow = (threadNo/matrixParam.blockRowSize) % matrixParam.blockColSize;

      extern __shared__ vectorType lowProd[];

      vectorType res = 0.0;

      if(rowInLevel && threadNo<activeThreads)
      {
         for(int index = matrixIndex(row, threadNo); index < d_dia_ptr[row] * matrixParam.blockSize; index+=activeThreads)
         {
            int blockCol = index % matrixParam.blockColSize;
            int blockNo = index/(matrixParam.blockSize);
            res += matrix[index] * prod[(d_col_ind[blockNo]) * matrixParam.blockColSize + blockCol];    // Compute L.x*
         }
      }

      if(rowInLevel && threadNo<activeThreads) atomicAdd(&lowProd[vectorIndex(row - d_partition_offsets[iPartition], threadNo)], res);    // Compute y = b - L.x* (subtracts the previously obtained lower product)

      __syncthreads();

      if(rowInLevel && threadNo<matrixParam.blockColSize) prod[vectorIndex(row, threadNo)] = vec[vectorIndex(row, threadNo)] - lowProd[vectorIndex(row - d_partition_offsets[iPartition], threadNo)];   // First condition just selects the first nVar threads of the warp - this is to exclude multiple writes
      
      __syncthreads();

      if(rowInLevel && threadNo==0) GaussEliminationKernel<<<precondParam.gaussElimGridDim, precondParam.gaussElimBlockDim, matrixParam.blockSize * sizeof(matrixType)>>>(&matrix[d_dia_ptr[row]], &prod[vectorIndex(row, 0)], d_dia_ptr, matrixParam);
      
      __syncthreads();
      
      #undef matrixIndex
      #undef vectorIndex
   } 
}


template<typename matrixType, typename vectorType>
__global__ void SecondSymmetricIterationKernel(matrixType* matrix, vectorType* prod, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, const unsigned long* d_dia_ptr, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn)
{
   int row = (blockIdx.x * blockDim.x + threadIdx.x)/MVP_WARP_SIZE;
   int threadNo = threadIdx.x%MVP_WARP_SIZE;
   int activeThreads = nEqn * (MVP_WARP_SIZE/nEqn);

   int blockRow = (threadNo/nEqn)%nVar;

   vectorType res = 0.0;

   if(row<nPointDomain && threadNo<activeThreads)
   {
      // Diagonal Product
      for(int index = d_dia_ptr[row] * nVar * nEqn + threadNo; index < (d_dia_ptr[row] + 1) * nVar * nEqn; index+=activeThreads)
      {
         int blockCol = index%nEqn;
         int blockNo = index/(nVar * nEqn);
         res += matrix[index] * prod[(d_col_ind[blockNo]) * nVar + blockCol];    // Compute D.x*
      }

      // Upper Product along with Simultaneous Subtraction
      for(int index = (d_dia_ptr[row] + 1) * nVar * nEqn + threadNo; index < d_row_ptr[row+1] * nVar * nEqn; index+=activeThreads)
      {
         int blockCol = index%nVar;
         int blockNo = index/(nVar * nEqn);
         res -= matrix[index] * prod[(d_col_ind[blockNo])*nVar + blockCol];    // Compute U.x_(n+1)
      }
   }

   __syncthreads();  // Syncthreads outside conditional statements to prevent exceptions

   if(row<nPointDomain && threadNo<nVar) prod[row * nVar + threadNo] = 0.0;    // First condition just selects the first nVar threads of the warp - this is to exclude multiple writes

   __syncthreads();

   if(row<nPointDomain && threadNo<activeThreads) atomicAdd(&prod[row * nVar + blockRow], res);   // Compute y = D.x*-U.x_(n+1)
}

template<typename matrixType, typename vectorType>
__global__ void MatrixVectorProductKernel(matrixType* matrix, vectorType* vec, vectorType* prod, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, matrixParameters matrixParam)
{
   #define matrixIndex(row, threadNo) d_row_ptr[row] * matrixParam.blockSize + threadNo
   #define vectorIndex(row, elemNo) row * matrixParam.blockColSize + elemNo

   int row = (blockIdx.x * blockDim.x + threadIdx.x)/MVP_WARP_SIZE;
   int threadNo = threadIdx.x % MVP_WARP_SIZE;
   int activeThreads = matrixParam.blockRowSize * (MVP_WARP_SIZE/matrixParam.blockRowSize);

   int blockRow = (threadNo/matrixParam.blockRowSize) % matrixParam.blockColSize;      // The remainder step is necessary as the blockRow sequnce should reset when we move to the next block

   if(row<matrixParam.totalRows && threadNo<activeThreads) prod[vectorIndex(row, blockRow)] = 0.0;

   __syncthreads();

   if(row<matrixParam.totalRows && threadNo<activeThreads)
   {
      vectorType res = 0.0;

      for(int index = matrixIndex(row, threadNo); index < matrixIndex(row+1, threadNo); index+=activeThreads)
      {
         int blockCol = index % matrixParam.blockColSize;
         int blockNo = index/(matrixParam.blockSize);
         res += matrix[index] * vec[(d_col_ind[blockNo]) * matrixParam.blockColSize + blockCol];
      }

      atomicAdd(&prod[vectorIndex(row, blockRow)], res);
   }

   #undef matrixIndex
   #undef vectorIndex
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
  int gridx = rounded_up_division(MVP_BLOCK_SIZE, matrixParam.totalRows * MVP_WARP_SIZE);
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
      precondParameters precondParam(matrixParam);

      dim3 blockDim(MVP_BLOCK_SIZE,1,1);
      int gridx = rounded_up_division(MVP_BLOCK_SIZE, geometry->maxPartitionSize * MVP_WARP_SIZE);
      dim3 gridDim(gridx, 1, 1);

      FirstSymmetricIterationKernel<<<gridDim, blockDim, geometry->maxPartitionSize * MVP_WARP_SIZE * sizeof(ScalarType)>>>(d_matrix, d_vec, d_prod, d_partition_offsets, d_row_ptr, d_col_ind, d_dia_ptr, matrixParam, precondParam);
      gpuErrChk( cudaPeekAtLastError() );
      SecondSymmetricIterationKernel<<<gridDim, blockDim>>>(d_matrix, d_prod, d_row_ptr, d_col_ind, d_dia_ptr, nPointDomain, nVar, nEqn);
      gpuErrChk( cudaPeekAtLastError() );
      prod.DtHTransfer();
   }

template class CSysMatrix<su2mixedfloat>; //This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.
