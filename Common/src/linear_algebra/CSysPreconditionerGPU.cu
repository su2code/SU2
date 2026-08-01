/*!
 * \file CSysPreconditionerGPU.cu
 * \brief CUDA/GPU skeleton implementations for matrix-based preconditioners.
 * \author Jesse Li
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

#include "../../include/linear_algebra/CSysMatrix.inl"
#include "../../include/linear_algebra/GPUComms.cuh"

namespace {

template <class ScalarType>
__global__ void ApplyJacobiPreconditionerKernel(const ScalarType* invM, const ScalarType* vec, ScalarType* prod,
                                                unsigned long nPointDomain, unsigned long nVar) {
  const auto iPoint = static_cast<unsigned long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (iPoint >= nPointDomain) return;

  const auto block = &invM[iPoint * nVar * nVar];
  const auto rhs = &vec[iPoint * nVar];
  auto out = &prod[iPoint * nVar];

  for (auto iVar = 0ul; iVar < nVar; ++iVar) {
    ScalarType sum = ScalarType(0);
    for (auto jVar = 0ul; jVar < nVar; ++jVar) {
      sum += block[iVar * nVar + jVar] * rhs[jVar];
    }
    out[iVar] = sum;
  }
}

}  // namespace

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildJacobiPreconditionerGPU() {
  SU2_ZONE_SCOPED

  if (nVar != nEqn) {
    SU2_MPI::Error("CUDA Jacobi preconditioner requires square blocks.", CURRENT_FUNCTION);
  }

  if (invM == nullptr) {
    SU2_MPI::Error("CUDA Jacobi preconditioner was requested without host inverse block storage.", CURRENT_FUNCTION);
  }

  if (d_invM == nullptr) {
    d_invM = GPUMemoryAllocation::gpu_alloc<ScalarType, true>(nPointDomain * nVar * nEqn * sizeof(ScalarType));
  }

  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    InverseDiagonalBlock(iPoint, &(invM[iPoint * nVar * nVar]));
  }

  gpuErrChk(cudaMemcpy(d_invM, invM, nPointDomain * nVar * nVar * sizeof(ScalarType), cudaMemcpyHostToDevice));
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeJacobiPreconditionerGPU(const CSysVector<ScalarType>& vec,
                                                            CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                            const CConfig* config) const {
  (void)geometry;
  (void)config;

  SU2_ZONE_SCOPED

  if (d_invM == nullptr) {
    SU2_MPI::Error("CUDA Jacobi preconditioner used before BuildJacobiPreconditionerGPU.", CURRENT_FUNCTION);
  }

  constexpr unsigned threadsPerBlock = 128;
  const auto blocks = static_cast<unsigned>((nPointDomain + threadsPerBlock - 1) / threadsPerBlock);
  ApplyJacobiPreconditionerKernel<<<blocks, threadsPerBlock>>>(d_invM, vec.GetDevicePointer(), prod.GetDevicePointer(),
                                                               nPointDomain, nVar);
  gpuErrChk(cudaPeekAtLastError());
}

template void CSysMatrix<su2mixedfloat>::BuildJacobiPreconditionerGPU();
template void CSysMatrix<su2mixedfloat>::ComputeJacobiPreconditionerGPU(const CSysVector<su2mixedfloat>& vec,
                                                                        CSysVector<su2mixedfloat>& prod,
                                                                        CGeometry* geometry,
                                                                        const CConfig* config) const;

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template void CSysMatrix<passivedouble>::BuildJacobiPreconditionerGPU();
template void CSysMatrix<passivedouble>::ComputeJacobiPreconditionerGPU(const CSysVector<passivedouble>& vec,
                                                                        CSysVector<passivedouble>& prod,
                                                                        CGeometry* geometry,
                                                                        const CConfig* config) const;
#endif
