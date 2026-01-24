/*!
 * \file ComputeLinSysResRMS.hpp
 * \brief Helper function to compute global RMS of linear system residual with MPI support.
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

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"
#include "../../../Common/include/CConfig.hpp"
#include <vector>
#include <cmath>

/*!
 * \brief Compute the global RMS of the linear system residual with MPI consistency.
 * \param[in] solver - Pointer to the solver containing the residual.
 * \param[in] geometry - Pointer to the geometry for domain point count.
 * \return Global RMS of the linear system residual across all MPI ranks.
 *
 * This function computes the RMS residual in an MPI-consistent manner by:
 * 1. Computing local sum of squared residuals over domain points
 * 2. Using MPI_Allreduce to sum across all ranks
 * 3. Using MPI_Allreduce to count total domain points
 * 4. Computing sqrt(globalSum / (globalPoints * nVar))
 *
 * All ranks will compute the identical global RMS value.
 */
inline su2double ComputeLinSysResRMS(const CSolver* solver, const CGeometry* geometry) {
  const unsigned short nVar = solver->GetnVar();
  const unsigned long nPointDomain = geometry->GetnPointDomain();
  
  /*--- Compute local sum of squared residuals ---*/
  std::vector<su2double> sumRes(nVar, 0.0);
  su2double localSum = 0.0;
  
  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    const su2double* res = solver->LinSysRes.GetBlock(iPoint);
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      sumRes[iVar] += res[iVar] * res[iVar];
    }
  }
  
  for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
    localSum += sumRes[iVar];
  }
  
  /*--- Global reduction of squared residuals ---*/
  su2double globalSum = 0.0;
  SU2_MPI::Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  
  /*--- Global reduction of domain point count ---*/
  unsigned long globalNPointDomain = 0;
  SU2_MPI::Allreduce(&nPointDomain, &globalNPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  
  /*--- Compute global RMS ---*/
  if (globalNPointDomain == 0) return 0.0;
  
  return std::sqrt(globalSum / (static_cast<su2double>(globalNPointDomain) * static_cast<su2double>(nVar)));
}
