#include<cuda_runtime.h>
#include "../../include/linear_algebra/CSysMatrix.inl"
#include "../../include/linear_algebra/CSysMatrix.hpp"
#include "../../include/geometry/CGeometry.hpp"

template<class ScalarType>
void CSysMatrix<ScalarType>::GPUMatrixStructMemAlloc()
{
    
}

__global__ void GPUMatrixVectorProductAdd(su2mixedfloat* matrix, double* vec, double* prod, unsigned long* d_row_ptr, unsigned long* d_col_ind, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   int j = threadIdx.y;
   int k = threadIdx.z;

   if(i<nPointDomain)
   {
      for(int index = d_row_ptr[i]; index<d_row_ptr[i+1]; index++)
      {
         int matrix_index = index * nVar * nEqn;
         int vec_index = d_col_ind[index] * nEqn;
         int prod_index = i * nVar;
      
         prod[prod_index + k] += matrix[ matrix_index + (i * nEqn + j)] * vec[vec_index + j];

      }
   }
}


template<class ScalarType>
void CSysMatrix<ScalarType>::GPUMatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                 CGeometry* geometry, const CConfig* config) const
                                                 {


  double* d_vec;
  double* d_prod;

  unsigned long nPointDomain = geometry->GetnPointDomain();
  unsigned long nDim = geometry->GetnDim();
  unsigned long nVar = nDim + 2 ;
  unsigned long nEqn = nVar;

  unsigned long mat_size = nnz*nVar*nEqn;
  unsigned long vec_size = nPointDomain*nVar;

  cudaMalloc((void**)(&d_vec), (sizeof(&vec[0])*vec_size));
  cudaMalloc((void**)(&d_prod), (sizeof(&prod[0])*vec_size));

  cudaError_t code = cudaGetLastError();
        if(code != cudaSuccess)
        {
            std::cerr << code << " Error Code " << std::endl;
        }

  cudaMemcpy((void*)(d_matrix), (void*)&matrix[0], (sizeof(&matrix[0])*mat_size), cudaMemcpyHostToDevice);
  cudaMemcpy((void*)(d_vec), (void*)&vec[0], (sizeof(&vec[0])*vec_size), cudaMemcpyHostToDevice);
  cudaMemcpy((void*)(d_prod), (void*)&prod[0], (sizeof(&prod[0])*vec_size), cudaMemcpyHostToDevice);

  code = cudaGetLastError();
        if(code != cudaSuccess)
        {
            std::cerr << code << " Error Code " << std::endl;
        }

  long xDim = floor(512.0/(nVar*nEqn));
  dim3 blockDim(xDim, nEqn, nVar);
  dim3 gridDim(ceil(nPointDomain/xDim), 1, 1);

  GPUMatrixVectorProductAdd<<<gridDim, blockDim>>>(d_matrix, d_vec, d_prod, d_row_ptr, d_col_ind, nPointDomain, nVar, nEqn);

  code = cudaGetLastError();
        if(code != cudaSuccess)
        {
            std::cerr << code << " Error Code " << std::endl;
        }

  cudaMemcpy((void*)(&prod[0]), (void*)d_prod, (sizeof(&prod[0])*vec_size), cudaMemcpyDeviceToHost);

  
  cudaFree(d_vec);
  cudaFree(d_prod);

}

template class CSysMatrix<su2mixedfloat>;