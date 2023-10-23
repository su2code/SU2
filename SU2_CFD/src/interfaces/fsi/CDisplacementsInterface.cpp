/*!
 * \file CDisplacementsInterface.cpp
 * \brief Main subroutines for transferring boundary displacements.
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

#include "../../../include/interfaces/fsi/CDisplacementsInterface.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

CDisplacementsInterface::CDisplacementsInterface(unsigned short val_nVar,
                                                 unsigned short val_nConst) :
  CInterface(val_nVar, val_nConst) {
}

void CDisplacementsInterface::GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry,
                                                const CConfig *struct_config, unsigned long Marker_Struct,
                                                unsigned long Vertex_Struct, unsigned long Point_Struct) {
  /*--- The displacements come from the predicted solution, but they are no longer incremental ---*/
  auto DisplacementDonor = struct_solution->GetNodes()->GetSolution_Pred(Point_Struct);

  for (auto iVar = 0u; iVar < nVar; iVar++)
    Donor_Variable[iVar] = DisplacementDonor[iVar];

  if (struct_config->GetTime_Domain()) {
    auto VelocityDonor = struct_solution->GetNodes()->GetSolution_Vel_Pred(Point_Struct);
    for (auto iVar = nVar/2; iVar < nVar; iVar++)//Assuming dynamic interface always has nVar = 2*nDim, 2D: 4, 3D: 6
      Donor_Variable[iVar] = VelocityDonor[iVar-nVar/2];
  }
}

void CDisplacementsInterface::SetTarget_Variable(CSolver *mesh_solver, CGeometry *flow_geometry,
                                                 const CConfig *flow_config, unsigned long Marker_Flow,
                                                 unsigned long Vertex_Flow, unsigned long Point_Mesh) {

  if (!flow_config->GetTime_Domain()) {
    /*--- Impose the boundary displacements ---*/
    mesh_solver->GetNodes()->SetBound_Disp(Point_Mesh,Target_Variable);
  } else {
    /*--- Impose the boundary displacements ---*/
    for (auto iVar = 0u; iVar < nVar/2; iVar++)
      mesh_solver->GetNodes()->SetBound_Disp(Point_Mesh,iVar,Target_Variable[iVar]);

    /*--- Impose the boundary velocities ---*/
    for (auto iVar = nVar/2; iVar < nVar; iVar++)
      mesh_solver->GetNodes()->SetBound_Vel(Point_Mesh,iVar-nVar/2,Target_Variable[iVar]);
  }
}
