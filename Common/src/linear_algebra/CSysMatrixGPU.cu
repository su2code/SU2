/*!
 * \file CSysMatrixGPU.cu
 * \brief Implementations of Kernels and Functions for Matrix Operations on the GPU
 * \author A. Raj
 * \version 8.4.0 "Harrier"
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
__global__ void GPUMatrixVectorProductAdd(matrixType* matrix, vectorType* vec, vectorType* prod, const unsigned long* d_row_ptr, const unsigned long* d_col_ind, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn)
{
   int row = (blockIdx.x * blockDim.x + threadIdx.x)/32;
   int threadNo = threadIdx.x%32;
   int activeThreads = nVar * (32/nVar);

   int blockRow = (threadNo/nVar)%nVar;

   if(row<nPointDomain && threadNo<activeThreads) prod[row * nVar + blockRow] = 0.0;

   __syncthreads();

   if(row<nPointDomain && threadNo<activeThreads)
   {
      vectorType res = 0.0;

      for(int index = d_row_ptr[row] * nVar * nEqn + threadNo; index < d_row_ptr[row+1] * nVar * nEqn; index+=activeThreads)
      {
         int blockCol = index%nEqn;
         int blockNo = index/(nVar * nEqn);
         res += matrix[index] * vec[(d_col_ind[blockNo])*nVar + blockCol];
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
  int gridx = KernelParameters::round_up_division(KernelParameters::MVP_WARP_SIZE, nPointDomain);
  dim3 gridDim(gridx, 1, 1);

  GPUMatrixVectorProductAdd<<<gridDim, blockDim>>>(d_matrix, d_vec, d_prod, d_row_ptr, d_col_ind, nPointDomain, nVar, nEqn);
  gpuErrChk( cudaPeekAtLastError() );

  prod.DtHTransfer();

}

template class CSysMatrix<su2mixedfloat>; //This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.
