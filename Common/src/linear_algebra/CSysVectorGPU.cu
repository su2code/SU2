/*!
 * \file CSysVectorGPU.cu
 * \brief Implementations of Kernels and Functions for Vector Operations on the GPU
 * \author A. Raj, D. Di giusto
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
#include <algorithm>
#include <type_traits>
#include <vector>

namespace {

/*--- cuBLAS handle for the reductions. Created on first use and kept for the lifetime of
 * the program, matching the fact that CUDA is either on or off for the whole run. ---*/
cublasHandle_t solver_blas_handle = nullptr;

cublasHandle_t GetBlasHandle() {
  if (solver_blas_handle == nullptr) {
    if (cublasCreate(&solver_blas_handle) != CUBLAS_STATUS_SUCCESS) {
      SU2_MPI::Error("cuBLAS handle creation failed for the GPU linear algebra.", CURRENT_FUNCTION);
    }
  }
  return solver_blas_handle;
}

}  // namespace

namespace VecExpr {

namespace {
/*--- Deliberately not thread local: every thread of an OpenMP team has to agree on this,
 * otherwise the team splits over the worksharing constructs in CSysVector. It is only
 * written by CSysSolve, outside any parallel region over the linear system. ---*/
bool use_device_expressions = false;
}  // namespace

bool UseDeviceExpressions() { return use_device_expressions; }

void SetUseDeviceExpressions(bool use) { use_device_expressions = use; }

}  // namespace VecExpr

template <class ScalarType>
void CSysVector<ScalarType>::HtDTransfer(bool trigger) const {
  SU2_ZONE_SCOPED
  if (trigger)
    gpuErrChk(cudaMemcpy((void*)(d_vec_val), (void*)&vec_val[0], (sizeof(ScalarType) * nElm), cudaMemcpyHostToDevice));
}

template <class ScalarType>
void CSysVector<ScalarType>::DtHTransfer(bool trigger) const {
  SU2_ZONE_SCOPED
  if (trigger)
    gpuErrChk(cudaMemcpy((void*)(&vec_val[0]), (void*)d_vec_val, (sizeof(ScalarType) * nElm), cudaMemcpyDeviceToHost));
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::dotGPU(const CSysVector& other) const {
  SU2_ZONE_SCOPED
  /*--- Both operands are already on the device, the caller owns the transfers. This
   * reduces over MPI, so it must be called by a single thread (see SU2_DEVICE_REGION). ---*/
  cublasHandle_t handle = GetBlasHandle();
  cublasStatus_t status = CUBLAS_STATUS_SUCCESS;

  ScalarType local_dot = ScalarType(0);

  if constexpr (std::is_same_v<ScalarType, float>) {
    status = cublasSdot(handle, static_cast<int>(nElmDomain), GetDevicePointer(), 1, other.GetDevicePointer(), 1,
                        &local_dot);
  } else if constexpr (std::is_same_v<ScalarType, double>) {
    status = cublasDdot(handle, static_cast<int>(nElmDomain), GetDevicePointer(), 1, other.GetDevicePointer(), 1,
                        &local_dot);
  } else {
    SU2_MPI::Error("Unsupported ScalarType in CSysVector::dotGPU.", CURRENT_FUNCTION);
    return ScalarType(0);
  }

  if (status != CUBLAS_STATUS_SUCCESS) {
    SU2_MPI::Error("cuBLAS dot failed in CSysVector::dotGPU.", CURRENT_FUNCTION);
    return ScalarType(0);
  }

  ScalarType global_dot = ScalarType(0);
  const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
  SelectMPIWrapper<ScalarType>::W::Allreduce(&local_dot, &global_dot, 1, mpi_type, MPI_SUM, SU2_MPI::GetComm());

  return global_dot;
}

/*!
 * \brief multi vector product with cublas<t>gemmBatched
 */
template <class ScalarType>
su2matrix<ScalarType> CSysVector<ScalarType>::multiDotGPU(const std::vector<CSysVector<ScalarType>>& V, const size_t i0,
                                                          const size_t n, const std::vector<CSysVector<ScalarType>>& W,
                                                          const size_t m) {
  /*--- The multiDot product between n V[size] and m W[size] vectors is performed as
   * a General Matrix Multiplication between two tall-skinny matrices:
   * C = \alpha * A^T * B + \beta * C
   * being A = V[ size * n ] and B = W[ size * m ] the batched vectors ---*/
  cublasHandle_t handle = GetBlasHandle();
  cublasStatus_t status = CUBLAS_STATUS_SUCCESS;

  const size_t size = V[0].nElmDomain;
  const size_t batch = n * m;

  /*--- Persistent device workspace, cached across calls and freed automatically
   * when the program exits (static local destruction), instead of leaking. ---*/
  struct Workspace {
    ScalarType* d_local = nullptr;
    const ScalarType** d_A = nullptr;
    const ScalarType** d_B = nullptr;
    ScalarType** d_C = nullptr;
    size_t capacity = 0;

    void EnsureCapacity(size_t batch) {
      if (batch <= capacity) return;
      cudaFree(d_local);
      cudaFree(d_A);
      cudaFree(d_B);
      cudaFree(d_C);
      gpuErrChk(cudaMalloc(&d_local, batch * sizeof(ScalarType)));
      gpuErrChk(cudaMalloc(&d_A, batch * sizeof(ScalarType*)));
      gpuErrChk(cudaMalloc(&d_B, batch * sizeof(ScalarType*)));
      gpuErrChk(cudaMalloc(&d_C, batch * sizeof(ScalarType*)));
      capacity = batch;
    }

    ~Workspace() {
      cudaFree(d_local);
      cudaFree(d_A);
      cudaFree(d_B);
      cudaFree(d_C);
    }
  };
  static Workspace ws;

  // allocate persistent result buffer local on host and device, is resized if needed
  su2matrix<ScalarType> local;
  local.resize(n, m);
  ws.EnsureCapacity(batch);

  // zero out the result buffer
  gpuErrChk(cudaMemset(ws.d_local, 0, batch * sizeof(ScalarType)));

  // prepare the arrays A,B,C on host
  static std::vector<const ScalarType*> h_A, h_B;
  static std::vector<ScalarType*> h_C;
  h_A.resize(batch); h_B.resize(batch); h_C.resize(batch);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j =0; j < m; ++j) {
      const size_t idx = i * m + j;
      h_A[idx] = V[i0 + i].GetDevicePointer();
      h_B[idx] = W[j].GetDevicePointer();
      h_C[idx] = ws.d_local + idx; // C maps to d_local to store the coefficients in the 2D array
    }
  }

  // copy pointers to device
  gpuErrChk(cudaMemcpy(ws.d_A, h_A.data(), batch * sizeof(ScalarType*), cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy(ws.d_B, h_B.data(), batch * sizeof(ScalarType*), cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy(ws.d_C, h_C.data(), batch * sizeof(ScalarType*), cudaMemcpyHostToDevice));

  // define alpha = 1.0 and beta = 0.0
  const auto alpha = ScalarType(1.0);
  const auto beta = ScalarType(0.0);

  if constexpr (std::is_same_v<ScalarType, float>) {
    status = cublasSgemmBatched(handle, CUBLAS_OP_T, CUBLAS_OP_N, 1, 1, size, &alpha, ws.d_A, static_cast<int>(size),
                                ws.d_B, static_cast<int>(size), &beta, ws.d_C, 1, static_cast<int>(batch));
  } else if constexpr (std::is_same_v<ScalarType, double>) {
    status = cublasDgemmBatched(handle, CUBLAS_OP_T, CUBLAS_OP_N, 1, 1, size, &alpha, ws.d_A, static_cast<int>(size),
                                ws.d_B, static_cast<int>(size), &beta, ws.d_C, 1, static_cast<int>(batch));
  } else {
    SU2_MPI::Error("Unsupported ScalarType in CSysVector::multiDotGPU.", CURRENT_FUNCTION);
    return local;
  }

  if (status != CUBLAS_STATUS_SUCCESS) {
    SU2_MPI::Error("cuBLAS cublas<t>gemmBatched failed in CSysVector::multiDotGPU.", CURRENT_FUNCTION);
    return local;
  }

  // copy result to host for MPI reduce
  gpuErrChk(cudaMemcpy(local.data(), ws.d_local, batch * sizeof(ScalarType), cudaMemcpyDeviceToHost));

  return local;
}

/*--- Every expression the solvers assign to a CSysVector needs its assignment kernel
 * instantiated here; the host compiler cannot emit one. A shape that is missing shows up
 * as an undefined reference to VecExpr::AssignDeviceExpression at link time, and is fixed
 * by adding a line to DEVICE_EXPRESSION_SHAPES below. The aliases use CSysVector (not
 * CVectorView) because that is how the operator overloads name their operands; store_t
 * turns it into a view when the node is built. ---*/
namespace {

template <class S>
using Vec = CSysVector<S>;
template <class S>
using Sca = VecExpr::Bcast<S>;

/*--- Leaves and the shapes of the FGMRES/GMRES basis updates. ---*/
template <class S>
using DeviceBcast = Sca<S>;
template <class S>
using DeviceView = VecExpr::CVectorView<S>;
template <class S>
using DeviceNeg = VecExpr::minus_<Vec<S>, S>;

/*--- vector * scalar and scalar * vector are distinct types, both are used. ---*/
template <class S>
using DeviceScale = VecExpr::mul_<Vec<S>, Sca<S>, S>;
template <class S>
using DeviceLScale = VecExpr::mul_<Sca<S>, Vec<S>, S>;
template <class S>
using DeviceDivScale = VecExpr::div_<Vec<S>, Sca<S>, S>;

/*--- Linear combinations, CSysSolve unrolls them up to four terms. ---*/
template <class S>
using DeviceScale2 = VecExpr::add_<DeviceScale<S>, DeviceScale<S>, S>;
template <class S>
using DeviceScale3 = VecExpr::add_<DeviceScale2<S>, DeviceScale<S>, S>;
template <class S>
using DeviceScale4 = VecExpr::add_<DeviceScale3<S>, DeviceScale<S>, S>;

/*--- r = b - A_x, in CG, BCGSTAB, Smoother and FGCRODR. ---*/
template <class S>
using DeviceSub = VecExpr::sub_<Vec<S>, Vec<S>, S>;

/*--- p = beta * p + z, in CG. ---*/
template <class S>
using DeviceLScalePlus = VecExpr::add_<DeviceLScale<S>, Vec<S>, S>;

/*--- p = beta * (p - omega * v) + r, in BCGSTAB. ---*/
template <class S>
using DeviceSubLScale = VecExpr::sub_<Vec<S>, DeviceLScale<S>, S>;
template <class S>
using DeviceLScaleSub = VecExpr::mul_<Sca<S>, DeviceSubLScale<S>, S>;
template <class S>
using DeviceBcgsDir = VecExpr::add_<DeviceLScaleSub<S>, Vec<S>, S>;

}  // namespace

#define DEVICE_EXPRESSION_SHAPES(SCALAR)                  \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceBcast);    \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceView);     \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceNeg);      \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceScale);    \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceLScale);   \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceDivScale); \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceScale2);   \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceScale3);   \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceScale4);   \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceSub);      \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceLScalePlus); \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceSubLScale);  \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceLScaleSub);  \
  INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, DeviceBcgsDir)

#define INSTANTIATE_DEVICE_ASSIGN(SCALAR, OP, EXPR)                                          \
  template void VecExpr::AssignDeviceExpression<VecExpr::DeviceAssignOp::OP, SCALAR, EXPR<SCALAR>>( \
      SCALAR*, unsigned long, const VecExpr::CVecExpr<EXPR<SCALAR>, SCALAR>&)

#define INSTANTIATE_DEVICE_ASSIGN_EXPR(SCALAR, EXPR) \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Assign, EXPR);   \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Add, EXPR);      \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Subtract, EXPR); \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Multiply, EXPR); \
  INSTANTIATE_DEVICE_ASSIGN(SCALAR, Divide, EXPR)

DEVICE_EXPRESSION_SHAPES(su2mixedfloat);

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
DEVICE_EXPRESSION_SHAPES(passivedouble);
#endif

#undef DEVICE_EXPRESSION_SHAPES
#undef INSTANTIATE_DEVICE_ASSIGN_EXPR
#undef INSTANTIATE_DEVICE_ASSIGN


template void CSysVector<su2mixedfloat>::HtDTransfer(bool trigger) const;
template void CSysVector<su2mixedfloat>::DtHTransfer(bool trigger) const;
template su2mixedfloat CSysVector<su2mixedfloat>::dotGPU(const CSysVector<su2mixedfloat>& other) const;
template su2matrix<su2mixedfloat> CSysVector<su2mixedfloat>::multiDotGPU(
    const std::vector<CSysVector<su2mixedfloat>>& V, size_t i0, size_t n,
    const std::vector<CSysVector<su2mixedfloat>>& W, size_t m);

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template void CSysVector<passivedouble>::HtDTransfer(bool trigger) const;
template void CSysVector<passivedouble>::DtHTransfer(bool trigger) const;
template passivedouble CSysVector<passivedouble>::dotGPU(const CSysVector<passivedouble>& other) const;
template su2matrix<passivedouble> CSysVector<passivedouble>::multiDotGPU(
    const std::vector<CSysVector<passivedouble>>& V, size_t i0, size_t n,
    const std::vector<CSysVector<passivedouble>>& W, size_t m);
#endif
