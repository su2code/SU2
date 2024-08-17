#include "../../include/linear_algebra/CSysMatrix.inl"
#include "../../include/linear_algebra/CSysMatrix.hpp"
#include "../../include/linear_algebra/GPU_lin_alg.cuh"
#include "../../include/geometry/CGeometry.hpp"

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

    prod[prod_index + j] = 0.0;

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

  unsigned long nPointDomain = geometry->GetnPointDomain();
  unsigned long nDim = geometry->GetnDim();
  unsigned long nVar = nDim + 2 ;
  unsigned long nEqn = nVar;

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

  gpuErrChk(cudaMemcpy((void*)(&prod[0]), (void*)d_prod, (sizeof(ScalarType)*vec_size), cudaMemcpyDeviceToHost));
 
  gpuErrChk(cudaFree(d_vec));
  gpuErrChk(cudaFree(d_prod));

}

template class CSysMatrix<su2mixedfloat>;