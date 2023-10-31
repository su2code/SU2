/*!
 * \file CConservativeVarsInterface.cpp
 * \brief Declaration and inlines of the class to transfer conservative variables
 *        from a generic zone into another one.
 * \author Ruben Sanchez
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

#include "../../../include/interfaces/cfd/CConservativeVarsInterface.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

CConservativeVarsInterface::CConservativeVarsInterface(unsigned short val_nVar, unsigned short val_nConst) :
  CInterface(val_nVar, val_nConst) {
}

void CConservativeVarsInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                                                   const CConfig *donor_config, unsigned long Marker_Donor,
                                                   unsigned long Vertex_Donor, unsigned long Point_Donor) {
  /*--- Retrieve solution and set it as the donor variable ---*/
  auto Solution = donor_solution->GetNodes()->GetSolution(Point_Donor);

  for (auto iVar = 0u; iVar < nVar; iVar++)
    Donor_Variable[iVar] = Solution[iVar];
}

void CConservativeVarsInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                                                    const CConfig *target_config, unsigned long Marker_Target,
                                                    unsigned long Vertex_Target, unsigned long Point_Target) {

  /*--- Set the target solution with the value of the Target Variable ---*/
  target_solution->GetNodes()->SetSolution(Point_Target,Target_Variable);
}
