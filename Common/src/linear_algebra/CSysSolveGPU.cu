/*!
 * \file CSysSolveGPU.cu
 * \brief CUDA/GPU skeleton implementations for Krylov solvers.
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

#include "../../include/linear_algebra/CSysSolve.hpp"
#include "../../include/linear_algebra/GPUComms.cuh"

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::FGMRES_LinSolver_GPU(const VectorType& b, VectorType& x,
                                                          const ProductType& mat_vec, const PrecondType& precond,
                                                          ScalarType tol, unsigned long m, ScalarType& residual,
                                                          bool monitoring, const CConfig* config) const {
  SU2_ZONE_SCOPED

  (void)b;
  (void)x;
  (void)mat_vec;
  (void)precond;
  (void)tol;
  (void)m;
  (void)residual;
  (void)monitoring;
  (void)config;

  /*--- Implementation area for the GPU-resident FGMRES solve path.
   * This skeleton intentionally leaves the algorithm unimplemented so the
   * surrounding dispatch, build integration, and file structure can be
   * reviewed independently of the CUDA solver details. ---*/
  SU2_MPI::Error("FGMRES GPU solver skeleton reached without an implementation.", CURRENT_FUNCTION);
  return 0;
}

/*--- Explicit instantiations for the GPU solver helper.
 * Keep these aligned with the scalar types instantiated for CSysSolve. ---*/
template unsigned long CSysSolve<su2mixedfloat>::FGMRES_LinSolver_GPU(
    const CSysVector<su2mixedfloat>& b, CSysVector<su2mixedfloat>& x,
    const CMatrixVectorProduct<su2mixedfloat>& mat_vec, const CPreconditioner<su2mixedfloat>& precond,
    su2mixedfloat tol, unsigned long m, su2mixedfloat& residual, bool monitoring, const CConfig* config) const;

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template unsigned long CSysSolve<passivedouble>::FGMRES_LinSolver_GPU(
    const CSysVector<passivedouble>& b, CSysVector<passivedouble>& x,
    const CMatrixVectorProduct<passivedouble>& mat_vec, const CPreconditioner<passivedouble>& precond,
    passivedouble tol, unsigned long m, passivedouble& residual, bool monitoring, const CConfig* config) const;
#endif
