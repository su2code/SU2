/*!
 * \file GPU_lin_alg.cu
 * \brief Implementation of Matrix Vector Product CUDA Kernel
 * \author A. Raj
 * \version 8.0.1 "Harrier"
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

   int i = blockIdx.x * blockDim.x + threadIdx.x;
   int j = threadIdx.y;
   int k = threadIdx.z;

   int prod_index = i * nVar;

   if(i<nPointDomain) 
   {
      prod[prod_index + j] = 0.0;
   }

   __syncthreads();

   vectorType res = 0.0;

   if(i<nPointDomain)
   {
      for(int index = d_row_ptr[i]; index<d_row_ptr[i+1]; index++)
      {
        int matrix_index = index * nVar * nEqn;
        int vec_index = d_col_ind[index] * nEqn;
      
        res += matrix[matrix_index + (j * nEqn + k)] * vec[vec_index + k];
      }

      atomicAdd(&prod[prod_index + j],res);
   }

}


template<class ScalarType>
void CSysMatrix<ScalarType>::GPUMatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                 CGeometry* geometry, const CConfig* config) const
                                                 {

  ScalarType* d_vec;
  ScalarType* d_prod;

  unsigned long mat_size = nnz*nVar*nEqn;
  unsigned long vec_size = nPointDomain*nVar;

  gpuErrChk(cudaMalloc((void**)(&d_vec), (sizeof(ScalarType)*vec_size)));
  gpuErrChk(cudaMalloc((void**)(&d_prod), (sizeof(ScalarType)*vec_size)));
  
  gpuErrChk(cudaMemcpy((void*)(d_matrix), (void*)&matrix[0], (sizeof(ScalarType)*mat_size), cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy((void*)(d_vec), (void*)&vec[0], (sizeof(ScalarType)*vec_size), cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy((void*)(d_prod), (void*)&prod[0], (sizeof(ScalarType)*vec_size), cudaMemcpyHostToDevice));

  double xDim = (double) 1024.0/(nVar*nEqn);
  dim3 blockDim(floor(xDim), nVar, nEqn);
  double gridx = (double) nPointDomain/xDim;
  dim3 gridDim(ceil(gridx), 1, 1);

  GPUMatrixVectorProductAdd<<<gridDim, blockDim>>>(d_matrix, d_vec, d_prod, d_row_ptr, d_col_ind, nPointDomain, nVar, nEqn);
  gpuErrChk( cudaPeekAtLastError() );

  gpuErrChk(cudaMemcpy((void*)(&prod[0]), (void*)d_prod, (sizeof(ScalarType)*vec_size), cudaMemcpyDeviceToHost));
 
  gpuErrChk(cudaFree(d_vec));
  gpuErrChk(cudaFree(d_prod));

}

template class CSysMatrix<su2mixedfloat>;