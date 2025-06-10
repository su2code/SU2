/*!
 * \file CSysMatrixGPU.cu
 * \brief Implementations of Kernels and Functions for Matrix Operations on the GPU
 * \author A. Raj
 * \version 8.1.0 "Harrier"
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
#include "../../include/linear_algebra/GPUComms.cuh"


template<typename matrixType, typename vectorType>
__global__ void GaussEliminationKernel(matrixType* matrix, vectorType* prod, const unsigned long* d_dia_ptr, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn, unsigned long size)
{
   #define A(I, J) matrixCopy[localRowNo * size + (I)*nVar + (J)]
   #define validRow row<nPointDomain

   int row = (blockDim.x*blockIdx.x + threadIdx.x)/size;
   int threadNo = threadIdx.x%size;

   /*Each block caters to blockDim.x/(nVar * nEqn) rows. This variable tells us which row out 
      of the blockDim.x/(nVar * nEqn) rows  in the block does the individual thread belong to*/
   int localRowNo = threadIdx.x/size;           
  
   int x = threadNo/nEqn;                                         //  A single block is in use so no remainder step necessary like in MVP
   int y = threadNo%nVar;

   matrixType pivot, weight = 0.0;

  extern __shared__ matrixType matrixCopy[];

   if(validRow) A(x, y) = matrix[d_dia_ptr[row] * size + x * nVar + y];

   __syncthreads();

   for(int currentRow=0; currentRow<nEqn; currentRow++)
   {
      pivot = A(currentRow, currentRow);             

      __syncthreads();

      if(validRow)
      {
         A(currentRow, y) /= pivot;                              
         if(y==0) prod[row * nVar + currentRow] /= pivot;
      }

      __syncthreads();
      
      weight = A(x, currentRow);                            

      if(x!=currentRow && validRow) 
      {
         A(x, y) -= weight * A(currentRow, y);       
         if(y==0) prod[row * nVar + x] -= weight * prod[row * nVar + currentRow];
      }

   __syncthreads();
   }
}

template<typename matrixType, typename vectorType>
__global__ void FirstSymmetricIterationKernel(matrixType* matrix, vectorType* vec, vectorType* prod, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, const unsigned long* d_dia_ptr, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn)
{
   int row = (blockIdx.x * blockDim.x + threadIdx.x)/KernelParameters::MVP_WARP_SIZE;
   int threadNo = threadIdx.x%KernelParameters::MVP_WARP_SIZE;
   int activeThreads = nEqn * (KernelParameters::MVP_WARP_SIZE/nEqn);

   int blockRow = (threadNo/nEqn)%nVar;
   
   vectorType res = 0.0;

   if(row<nPointDomain && threadNo<activeThreads)
   {
      for(int index = d_row_ptr[row] * nVar * nEqn + threadNo; index < d_dia_ptr[row] * nVar * nEqn; index+=activeThreads)
      {
         int blockCol = index%nVar;
         int blockNo = index/(nVar * nEqn);
         res -= matrix[index] * prod[(d_col_ind[blockNo])*nVar + blockCol];    // Compute -L.x*
      }
   }

   __syncthreads();   // Syncthreads outside the conditional statements to prevent excecptions

   if(row<nPointDomain && threadNo<nVar) prod[row * nVar + threadNo]  = vec[row * nVar + threadNo];   // First condition just selects the first nVar threads of the warp - this is to exclude multiple writes
   
   __syncthreads(); 

   if(row<nPointDomain && threadNo<activeThreads) atomicAdd(&prod[row * nVar + blockRow], res);    // Compute y = b - L.x* (subtracts the previously obtained lower product)
}

template<typename matrixType, typename vectorType>
__global__ void SecondSymmetricIterationKernel(matrixType* matrix, vectorType* prod, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, const unsigned long* d_dia_ptr, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn)
{
   int row = (blockIdx.x * blockDim.x + threadIdx.x)/KernelParameters::MVP_WARP_SIZE;
   int threadNo = threadIdx.x%KernelParameters::MVP_WARP_SIZE;
   int activeThreads = nEqn * (KernelParameters::MVP_WARP_SIZE/nEqn);

   int blockRow = (threadNo/nEqn)%nVar;

   vectorType res = 0.0;

   if(row==nPointDomain-3) 
   {
      int hi = 1;
   }

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
__global__ void MatrixVectorProductKernel(matrixType* matrix, vectorType* vec, vectorType* prod, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn)
{
   int row = (blockIdx.x * blockDim.x + threadIdx.x)/KernelParameters::MVP_WARP_SIZE;
   int threadNo = threadIdx.x%KernelParameters::MVP_WARP_SIZE;
   int activeThreads = nEqn * (KernelParameters::MVP_WARP_SIZE/nEqn);

   int blockRow = (threadNo/nEqn)%nVar;      // The remainder step is necessary as the blockRow sequnce should reset when we move to the next block

   if(row<nPointDomain && threadNo<activeThreads) prod[row * nVar + blockRow] = 0.0;

   __syncthreads();

   if(row<nPointDomain && threadNo<activeThreads)
   {
      vectorType res = 0.0;

      for(int index = d_row_ptr[row] * nVar * nEqn + threadNo; index < d_row_ptr[row+1] * nVar * nEqn; index+=activeThreads)
      {
         int blockCol = index%nVar;
         int blockNo = index/(nVar * nEqn);
         res += matrix[index] * vec[(d_col_ind[blockNo]) * nVar + blockCol];
      }

      atomicAdd(&prod[row * nVar + blockRow], res);
   }
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

  dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE,1,1);
  int gridx = KernelParameters::rounded_up_division(KernelParameters::MVP_BLOCK_SIZE, nPointDomain * KernelParameters::MVP_WARP_SIZE);
  dim3 gridDim(gridx, 1, 1);

  MatrixVectorProductKernel<<<gridDim, blockDim>>>(d_matrix, d_vec, d_prod, d_row_ptr, d_col_ind, nPointDomain, nVar, nEqn);
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

      unsigned long size = nVar * nEqn;

      dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE,1,1);
      int gridx = KernelParameters::rounded_up_division(KernelParameters::MVP_BLOCK_SIZE, nPointDomain * KernelParameters::MVP_WARP_SIZE);
      dim3 gridDim(gridx, 1, 1);

      int geBlockx = size * (KernelParameters::MVP_BLOCK_SIZE/size);
      dim3 gaussElimBlockDim(geBlockx, 1, 1);
      int geGridx = KernelParameters::rounded_up_division(geBlockx, size * nPointDomain);
      dim3 gaussElimGridDim(geGridx,1,1);

      FirstSymmetricIterationKernel<<<gridDim, blockDim>>>(d_matrix, d_vec, d_prod, d_row_ptr, d_col_ind, d_dia_ptr, nPointDomain, nVar, nEqn);
      gpuErrChk( cudaPeekAtLastError() );
     GaussEliminationKernel<<<gaussElimGridDim, gaussElimBlockDim, sizeof(ScalarType) * geBlockx>>>(d_matrix, d_prod, d_dia_ptr, nPointDomain, nVar, nEqn, size);
      gpuErrChk( cudaPeekAtLastError() );
     SecondSymmetricIterationKernel<<<gridDim, blockDim>>>(d_matrix, d_prod, d_row_ptr, d_col_ind, d_dia_ptr, nPointDomain, nVar, nEqn);
      gpuErrChk( cudaPeekAtLastError() );
      GaussEliminationKernel<<<gaussElimGridDim, gaussElimBlockDim, sizeof(ScalarType) * geBlockx>>>(d_matrix, d_prod, d_dia_ptr, nPointDomain, nVar, nEqn, size);
      gpuErrChk( cudaPeekAtLastError() );

      prod.DtHTransfer();
   }

template class CSysMatrix<su2mixedfloat>; //This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.
