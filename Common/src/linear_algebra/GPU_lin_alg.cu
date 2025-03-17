/*!
 * \file GPU_lin_alg.cu
 * \brief Implementation of Matrix Vector Product CUDA Kernel
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
#include "../../include/linear_algebra/GPU_lin_alg.cuh"

#ifndef gpuErrChk
#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
#endif

template<typename matrixType, typename vectorType>
__global__ void GPUMatrixVectorProductAdd(matrixType* matrix, vectorType* vec, vectorType* prod, unsigned long* d_row_ptr, unsigned long* d_col_ind, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn)
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
void CSysMatrix<ScalarType>::GPUMatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                 CGeometry* geometry, const CConfig* config) const
                                                 {

  ScalarType* d_vec;
  ScalarType* d_prod;

  unsigned long vec_size = nPointDomain*nVar;

  gpuErrChk(cudaMalloc((void**)(&d_vec), (sizeof(ScalarType)*vec_size)));
  gpuErrChk(cudaMalloc((void**)(&d_prod), (sizeof(ScalarType)*vec_size)));

  gpuErrChk(cudaMemcpy((void*)(d_vec), (void*)&vec[0], (sizeof(ScalarType)*vec_size), cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemset((void*)(d_prod), 0.0, (sizeof(ScalarType)*vec_size)));

  dim3 blockDim(1024,1,1);
  double gridx = (double) nPointDomain/32.0;
  gridx = double(ceil(gridx));
  dim3 gridDim(gridx, 1.0, 1.0);

  GPUMatrixVectorProductAdd<<<gridDim, blockDim>>>(matrix, d_vec, d_prod, d_row_ptr, d_col_ind, nPointDomain, nVar, nEqn);
  gpuErrChk( cudaPeekAtLastError() );

  gpuErrChk(cudaMemcpy((void*)(&prod[0]), (void*)d_prod, (sizeof(ScalarType)*vec_size), cudaMemcpyDeviceToHost));

  gpuErrChk(cudaFree(d_vec));
  gpuErrChk(cudaFree(d_prod));

}

template class CSysMatrix<su2mixedfloat>;
