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

struct CudaSpMVResources {
  cusparseHandle_t handle = nullptr;
  cusparseConstSpMatDescr_t mat = nullptr;
  size_t buffer_size = 0;
  void* buffer = nullptr;
};

void ReleaseCudaSpMVResources(CudaSpMVResources*& resources) {
  if (resources == nullptr) {
    return;
  }

  if (resources->buffer != nullptr) {
    gpuErrChk(cudaFree(resources->buffer));
    resources->buffer = nullptr;
    resources->buffer_size = 0;
  }

  if (resources->mat != nullptr) {
    cusparseErrChk(cusparseDestroySpMat(resources->mat));
    resources->mat = nullptr;
  }

  if (resources->handle != nullptr) {
    cusparseErrChk(cusparseDestroy(resources->handle));
    resources->handle = nullptr;
  }

  delete resources;
  resources = nullptr;
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

  if (spmv_resources == nullptr) {
    spmv_resources = new CudaSpMVResources;

    cusparseErrChk(cusparseCreate(&spmv_resources->handle));

    cusparseErrChk(cusparseCreateConstBsr(&spmv_resources->mat, brows, bcols, bnnz, blockSize, blockSize, d_row_ptr,
                                          d_col_ind, d_matrix, indexType, indexType, CUSPARSE_INDEX_BASE_ZERO,
                                          valueType, CUSPARSE_ORDER_ROW));
  }

  cusparseDnVecDescr_t vecX = nullptr;
  cusparseDnVecDescr_t vecY = nullptr;

  cusparseErrChk(cusparseCreateDnVec(&vecX, xSize, d_vec, valueType));
  cusparseErrChk(cusparseCreateDnVec(&vecY, ySize, d_prod, valueType));

  size_t required_buffer_size = 0;

  cusparseErrChk(cusparseSpMV_bufferSize(spmv_resources->handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
                                         spmv_resources->mat, vecX, &beta, vecY, valueType, CUSPARSE_SPMV_BSR_ALG1,
                                         &required_buffer_size));

  if (required_buffer_size > spmv_resources->buffer_size) {
    if (spmv_resources->buffer != nullptr) {
      gpuErrChk(cudaFree(spmv_resources->buffer));
    }

    if (required_buffer_size > 0) {
      gpuErrChk(cudaMalloc(&spmv_resources->buffer, required_buffer_size));
    } else {
      spmv_resources->buffer = nullptr;
    }

    spmv_resources->buffer_size = required_buffer_size;
  }

  cusparseErrChk(cusparseSpMV(spmv_resources->handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, spmv_resources->mat,
                              vecX, &beta, vecY, valueType, CUSPARSE_SPMV_BSR_ALG1, spmv_resources->buffer));

  cusparseErrChk(cusparseDestroyDnVec(vecY));
  cusparseErrChk(cusparseDestroyDnVec(vecX));
}
template void CSysMatrix<su2mixedfloat>::HtDTransfer(bool trigger) const;
template void CSysMatrix<su2mixedfloat>::GPUMatrixVectorProduct(const CSysVector<su2mixedfloat>& vec,
                                                                CSysVector<su2mixedfloat>& prod, CGeometry* geometry,
                                                                const CConfig* config) const;

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template void CSysMatrix<passivedouble>::HtDTransfer(bool trigger) const;
template void CSysMatrix<passivedouble>::GPUMatrixVectorProduct(const CSysVector<passivedouble>& vec,
                                                                CSysVector<passivedouble>& prod, CGeometry* geometry,
                                                                const CConfig* config) const;
#endif
