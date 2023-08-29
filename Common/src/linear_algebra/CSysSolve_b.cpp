/*!
 * \file CSysSolve_b.cpp
 * \brief Routines for the linear solver used in the reverse sweep of AD.
 * \author T. Albring, J. Blühdorn
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/linear_algebra/CSysSolve_b.hpp"
#include "../../include/linear_algebra/CSysSolve.hpp"
#include "../../include/linear_algebra/CSysMatrix.hpp"
#include "../../include/linear_algebra/CSysVector.hpp"

#ifdef CODI_REVERSE_TYPE
template <class ScalarType>
void CSysSolve_b<ScalarType>::Solve_b(const su2double::Real* x, su2double::Real* x_b, size_t m,
                                      const su2double::Real* y, const su2double::Real* y_b, size_t n,
                                      codi::ExternalFunctionUserData* d) {
  CSysVector<su2double>* LinSysRes_b = nullptr;
  d->getDataByIndex(LinSysRes_b, 0);

  CSysVector<su2double>* LinSysSol_b = nullptr;
  d->getDataByIndex(LinSysSol_b, 1);

  CSysMatrix<ScalarType>* Jacobian = nullptr;
  d->getDataByIndex(Jacobian, 2);

  CGeometry* geometry = nullptr;
  d->getDataByIndex(geometry, 3);

  const CConfig* config = nullptr;
  d->getDataByIndex(config, 4);

  CSysSolve<ScalarType>* solver = nullptr;
  d->getDataByIndex(solver, 5);

  /*--- Initialize the right-hand side with the gradient of the solution of the primal linear system ---*/

  SU2_OMP_BARRIER
  SU2_OMP_FOR_STAT(roundUpDiv(n, omp_get_num_threads()))
  for (unsigned long i = 0; i < n; i++) {
    (*LinSysRes_b)[i] = y_b[i];
    (*LinSysSol_b)[i] = 0.0;
  }
  END_SU2_OMP_FOR

  solver->Solve_b(*Jacobian, *LinSysRes_b, *LinSysSol_b, geometry, config, false);

  SU2_OMP_FOR_STAT(roundUpDiv(n, omp_get_num_threads()))
  for (unsigned long i = 0; i < n; i++) {
    x_b[i] = SU2_TYPE::GetValue((*LinSysSol_b)[i]);
  }
  END_SU2_OMP_FOR
}

template class CSysSolve_b<su2mixedfloat>;
#ifdef USE_MIXED_PRECISION
template class CSysSolve_b<passivedouble>;
#endif
#endif
