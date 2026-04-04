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
#include <cmath>

/*!
 * \brief Compute the global (MPI-reduced) RMS of LinSysRes over all variables and domain points.
 * \param[in] solver - Solver whose LinSysRes is evaluated.
 * \return Global RMS value.
 */
inline su2double ComputeLinSysResRMS(const CSolver* solver) {
  unsigned long nElmDomain = solver->LinSysRes.GetNElmDomain();
  unsigned long globalNElmDomain = 0;
  SU2_MPI::Allreduce(&nElmDomain, &globalNElmDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  if (globalNElmDomain == 0) return 0.0;

  const su2double squaredNorm = solver->LinSysRes.squaredNorm();
  return std::sqrt(squaredNorm / static_cast<su2double>(globalNElmDomain));
}
