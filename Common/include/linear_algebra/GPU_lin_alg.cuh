#include<cuda_runtime.h>
#include"../../include/linear_algebra/CSysMatrix.hpp"

__global__ void GPUMatrixVectorProductAdd(su2mixedfloat* matrix, double* vec, double* prod, unsigned long* row_ptr, unsigned long* col_ind);