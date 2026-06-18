/*!
 * \file CSysVectorGPU.cu
 * \brief Implementations of Kernels and Functions for Vector Operations on the GPU
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

#include "../../include/linear_algebra/CSysVector.hpp"
#include "../../include/linear_algebra/GPUComms.cuh"
#include <cublas_v2.h>
#include <type_traits>

namespace {

thread_local cublasHandle_t active_solver_blas_handle = nullptr;
thread_local unsigned active_solver_blas_depth = 0;

cublasHandle_t GetBlasHandle(const char* creation_error, bool& owns_handle) {
  if (active_solver_blas_handle != nullptr) {
    owns_handle = false;
    return active_solver_blas_handle;
  }

  cublasHandle_t handle = nullptr;
  owns_handle = true;
  if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS) {
    SU2_MPI::Error(creation_error, CURRENT_FUNCTION);
  }
  return handle;
}

void ReleaseBlasHandle(cublasHandle_t handle, bool owns_handle, const char* destruction_error) {
  if (owns_handle && handle != nullptr && cublasDestroy(handle) != CUBLAS_STATUS_SUCCESS) {
    SU2_MPI::Error(destruction_error, CURRENT_FUNCTION);
  }
}

}  // namespace

void SU2_GPU_BeginSolverBLASContext() {
  if (active_solver_blas_depth == 0) {
    cublasHandle_t handle = nullptr;
    if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS) {
      SU2_MPI::Error("cuBLAS handle creation failed for the GPU linear solver context.", CURRENT_FUNCTION);
    }
    active_solver_blas_handle = handle;
  }
  ++active_solver_blas_depth;
}

void SU2_GPU_EndSolverBLASContext() {
  if (active_solver_blas_depth == 0) {
    SU2_MPI::Error("GPU linear solver BLAS context ended without a matching begin.", CURRENT_FUNCTION);
  }

  --active_solver_blas_depth;
  if (active_solver_blas_depth == 0) {
    auto status = cublasDestroy(active_solver_blas_handle);
    active_solver_blas_handle = nullptr;
    if (status != CUBLAS_STATUS_SUCCESS) {
      SU2_MPI::Error("cuBLAS handle destruction failed for the GPU linear solver context.", CURRENT_FUNCTION);
    }
  }
}

template <class ScalarType>
void CSysVector<ScalarType>::HtDTransfer(bool trigger) const {
  if (trigger)
    gpuErrChk(cudaMemcpy((void*)(d_vec_val), (void*)&vec_val[0], (sizeof(ScalarType) * nElm), cudaMemcpyHostToDevice));
}

template <class ScalarType>
void CSysVector<ScalarType>::DtHTransfer(bool trigger) const {
  if (trigger)
    gpuErrChk(cudaMemcpy((void*)(&vec_val[0]), (void*)d_vec_val, (sizeof(ScalarType) * nElm), cudaMemcpyDeviceToHost));
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::GPUDot(const CSysVector& other) const {
  bool owns_handle = false;
  cublasHandle_t handle = GetBlasHandle("cuBLAS handle creation failed in CSysVector::GPUDot.", owns_handle);
  cublasStatus_t status = CUBLAS_STATUS_SUCCESS;

  ScalarType local_dot = ScalarType(0);

  if constexpr (std::is_same_v<ScalarType, float>) {
    status = cublasSdot(handle, static_cast<int>(nElmDomain), GetDevicePointer(), 1, other.GetDevicePointer(), 1,
                        &local_dot);
  } else if constexpr (std::is_same_v<ScalarType, double>) {
    status = cublasDdot(handle, static_cast<int>(nElmDomain), GetDevicePointer(), 1, other.GetDevicePointer(), 1,
                        &local_dot);
  } else {
    ReleaseBlasHandle(handle, owns_handle, "cuBLAS handle destruction failed in CSysVector::GPUDot.");
    SU2_MPI::Error("Unsupported ScalarType in CSysVector::GPUDot.", CURRENT_FUNCTION);
    return ScalarType(0);
  }

  if (status != CUBLAS_STATUS_SUCCESS) {
    ReleaseBlasHandle(handle, owns_handle, "cuBLAS handle destruction failed in CSysVector::GPUDot.");
    SU2_MPI::Error("cuBLAS dot failed in CSysVector::GPUDot.", CURRENT_FUNCTION);
    return ScalarType(0);
  }

  ReleaseBlasHandle(handle, owns_handle, "cuBLAS handle destruction failed in CSysVector::GPUDot.");

  ScalarType global_dot = ScalarType(0);
  const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
  SelectMPIWrapper<ScalarType>::W::Allreduce(&local_dot, &global_dot, 1, mpi_type, MPI_SUM, SU2_MPI::GetComm());

  return global_dot;
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::GPUNorm() const {
  return sqrt(GPUDot(*this));
}

template void CSysVector<su2mixedfloat>::HtDTransfer(bool trigger) const;
template void CSysVector<su2mixedfloat>::DtHTransfer(bool trigger) const;
template su2mixedfloat CSysVector<su2mixedfloat>::GPUDot(const CSysVector<su2mixedfloat>& other) const;
template su2mixedfloat CSysVector<su2mixedfloat>::GPUNorm() const;

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template void CSysVector<passivedouble>::HtDTransfer(bool trigger) const;
template void CSysVector<passivedouble>::DtHTransfer(bool trigger) const;
template passivedouble CSysVector<passivedouble>::GPUDot(const CSysVector<passivedouble>& other) const;
template passivedouble CSysVector<passivedouble>::GPUNorm() const;
#endif

#ifdef CODI_REVERSE_TYPE
template void CSysVector<su2double>::HtDTransfer(bool trigger) const;
template void CSysVector<su2double>::DtHTransfer(bool trigger) const;
template su2double CSysVector<su2double>::GPUDot(const CSysVector<su2double>& other) const;
template su2double CSysVector<su2double>::GPUNorm() const;
#endif
