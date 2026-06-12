/*!
 * \file CSysMatrixGPU.cu
 * \brief Implementations of Kernels and Functions for Matrix Operations on the GPU
 * \author A. Raj
 * \version 8.5.0 "Harrier"
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

#include <cusparse.h>
#include <cstdint>
#include <type_traits>

inline void cusparseAssert(cusparseStatus_t code, const char* file, int line, bool abort = true) {
  if (code != CUSPARSE_STATUS_SUCCESS) {
    fprintf(stderr, "cuSPARSEassert: %s %s %d\n", cusparseGetErrorString(code), file, line);
    if (abort) exit(static_cast<int>(code));
  }
}

#define cusparseErrChk(ans)                    \
  {                                            \
    cusparseAssert((ans), __FILE__, __LINE__); \
  }

inline cusparseIndexType_t GetCusparseIndexType() {
  if constexpr (sizeof(unsigned long) == 4) {
    return CUSPARSE_INDEX_32I;
  } else if constexpr (sizeof(unsigned long) == 8) {
    return CUSPARSE_INDEX_64I;
  } else {
    static_assert(sizeof(unsigned long) == 4 || sizeof(unsigned long) == 8,
                  "cuSPARSE BSR SpMV only supports 32-bit or 64-bit index arrays in this path.");
  }
}

template <class ScalarType>
constexpr cudaDataType GetCudaDataType() {
  if constexpr (std::is_same<ScalarType, float>::value) {
    return CUDA_R_32F;
  } else if constexpr (std::is_same<ScalarType, double>::value) {
    return CUDA_R_64F;
  } else {
    static_assert(std::is_same<ScalarType, float>::value || std::is_same<ScalarType, double>::value,
                  "cuSPARSE BSR SpMV only supports float and double in this path.");
  }
}

template<class ScalarType>
void CSysMatrix<ScalarType>::HtDTransfer(bool trigger) const
{
   if(trigger) gpuErrChk(cudaMemcpy((void*)(d_matrix), (void*)&matrix[0], (sizeof(ScalarType)*nnz*nVar*nEqn), cudaMemcpyHostToDevice));
}

template <class ScalarType>
void CSysMatrix<ScalarType>::GPUMatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                    CGeometry* geometry, const CConfig* config) const {
  if (nVar != nEqn) {
    SU2_MPI::Error("CUDA CSysMatrix matvec with cuSPARSE BSR requires square blocks.", CURRENT_FUNCTION);
  }

  ScalarType* d_vec = vec.GetDevicePointer();
  ScalarType* d_prod = prod.GetDevicePointer();

  vec.HtDTransfer();

  const auto indexType = GetCusparseIndexType();
  const auto valueType = GetCudaDataType<ScalarType>();

  const std::int64_t blockSize = static_cast<std::int64_t>(nVar);

  const std::int64_t brows = static_cast<std::int64_t>(nPointDomain);
  const std::int64_t bcols = static_cast<std::int64_t>(nPoint);
  const std::int64_t bnnz = static_cast<std::int64_t>(nnz);

  const std::int64_t xSize = static_cast<std::int64_t>(nPoint) * blockSize;
  const std::int64_t ySize = static_cast<std::int64_t>(nPointDomain) * blockSize;

  const ScalarType alpha = 1.0;
  const ScalarType beta = 0.0;

  cusparseHandle_t handle = nullptr;
  cusparseConstSpMatDescr_t matA = nullptr;
  cusparseDnVecDescr_t vecX = nullptr;
  cusparseDnVecDescr_t vecY = nullptr;

  cusparseErrChk(cusparseCreate(&handle));

  cusparseErrChk(cusparseCreateConstBsr(&matA, brows, bcols, bnnz, blockSize, blockSize, d_row_ptr, d_col_ind, d_matrix,
                                        indexType, indexType, CUSPARSE_INDEX_BASE_ZERO, valueType, CUSPARSE_ORDER_ROW));

  cusparseErrChk(cusparseCreateDnVec(&vecX, xSize, d_vec, valueType));
  cusparseErrChk(cusparseCreateDnVec(&vecY, ySize, d_prod, valueType));

  size_t bufferSize = 0;
  void* dBuffer = nullptr;

  cusparseErrChk(cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY,
                                         valueType, CUSPARSE_SPMV_BSR_ALG1, &bufferSize));

  if (bufferSize > 0) {
    gpuErrChk(cudaMalloc(&dBuffer, bufferSize));
  }

  cusparseErrChk(cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, valueType,
                              CUSPARSE_SPMV_BSR_ALG1, dBuffer));

  if (dBuffer != nullptr) {
    gpuErrChk(cudaFree(dBuffer));
  }

  cusparseErrChk(cusparseDestroyDnVec(vecY));
  cusparseErrChk(cusparseDestroyDnVec(vecX));
  cusparseErrChk(cusparseDestroySpMat(matA));
  cusparseErrChk(cusparseDestroy(handle));

  prod.DtHTransfer();
}

template class CSysMatrix<su2mixedfloat>; //This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.
