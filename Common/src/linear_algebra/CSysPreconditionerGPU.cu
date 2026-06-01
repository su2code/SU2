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

#include "../../include/linear_algebra/CSysMatrix.hpp"
#include "../../include/linear_algebra/GPUComms.cuh"

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildJacobiPreconditionerGPU() {
  /*--- Implementation area for building the GPU/device Jacobi preconditioner.
   * This skeleton intentionally leaves the algorithm unimplemented so the
   * surrounding dispatch and file structure can be reviewed independently of
   * the CUDA preconditioner details. ---*/
  SU2_MPI::Error("CSysMatrix::BuildJacobiPreconditionerGPU skeleton reached without an implementation.",
                 CURRENT_FUNCTION);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeJacobiPreconditionerGPU(const CSysVector<ScalarType>& vec,
                                                            CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                            const CConfig* config) const {
  (void)vec;
  (void)prod;
  (void)geometry;
  (void)config;

  /*--- Implementation area for applying the GPU/device Jacobi preconditioner.
   * This skeleton intentionally leaves the algorithm unimplemented so the
   * surrounding dispatch and file structure can be reviewed independently of
   * the CUDA preconditioner details. ---*/
  SU2_MPI::Error("CSysMatrix::ComputeJacobiPreconditionerGPU skeleton reached without an implementation.",
                 CURRENT_FUNCTION);
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
