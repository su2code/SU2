/*!
 * \file ComputeLinSysResRMS.hpp
 * \brief Helper to compute the global RMS of LinSysRes across all variables and domain points.
 * \author Nijso Beishuizen
 * \version 8.4.0 "Harrier"
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

#pragma once
#include "../../include/solvers/CSolver.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include <cmath>

/*!
 * \brief Compute the global (MPI-reduced) RMS of LinSysRes over all variables and domain points.
 *
 * \note Thread-safety: This function MUST be called by ALL threads in the current
 *       OpenMP parallel region, because squaredNorm() uses parallel for + barriers
 *       internally.  Do NOT call from inside BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS.
 *       The return value is correct only on the master thread (thread 0).
 *
 * \param[in] solver - Solver whose LinSysRes is evaluated.
 * \return Global RMS value (valid on master thread; other threads return 0).
 */
inline passivedouble ComputeLinSysResRMS(const CSolver* solver) {

  /*--- squaredNorm() -> dot() uses OMP parallel for + barriers internally,
   *    so all threads must participate. ---*/
  const su2double sqNorm = solver->LinSysRes.squaredNorm();

  /*--- The MPI reduction for nElmDomain must be single-threaded. ---*/
  passivedouble result = 0.0;
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    unsigned long nElmDomain = solver->LinSysRes.GetNElmDomain();
    unsigned long globalNElmDomain = 0;
    SU2_MPI::Allreduce(&nElmDomain, &globalNElmDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
    if (globalNElmDomain > 0)
      result = std::sqrt(SU2_TYPE::GetValue(sqNorm) / static_cast<passivedouble>(globalNElmDomain));
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  return result;
}
