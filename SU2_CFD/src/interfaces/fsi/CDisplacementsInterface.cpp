/*!
 * \file CDisplacementsInterface.cpp
 * \brief Main subroutines for transferring boundary displacements.
 * \author Ruben Sanchez
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

CDisplacementsInterface::CDisplacementsInterface(unsigned short val_nVar,
                                                 unsigned short val_nConst,
                                                 CConfig *config) :
  CInterface(val_nVar, val_nConst, config) {

}

CDisplacementsInterface::~CDisplacementsInterface(void) {

}


void CDisplacementsInterface::GetPhysical_Constants(CSolver *struct_solution, CSolver *flow_solution,
                                                    CGeometry *struct_geometry, CGeometry *flow_geometry,
                                                    CConfig *struct_config, CConfig *flow_config) {
}

void CDisplacementsInterface::GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry,
                                                CConfig *struct_config, unsigned long Marker_Struct,
                                                unsigned long Vertex_Struct, unsigned long Point_Struct) {

  su2double *DisplacementDonor;
  unsigned short iVar;

  /*--- The displacements come from the predicted solution, but they are no longer incremental ---*/
  DisplacementDonor = struct_solution->GetNodes()->GetSolution_Pred(Point_Struct);

  for (iVar = 0; iVar < nVar; iVar++)
    Donor_Variable[iVar] = DisplacementDonor[iVar];

}

void CDisplacementsInterface::SetTarget_Variable(CSolver *mesh_solver, CGeometry *flow_geometry,
                                                 CConfig *flow_config, unsigned long Marker_Flow,
                                                 unsigned long Vertex_Flow, unsigned long Point_Mesh) {

  /*--- Impose the boundary displacements ---*/

  mesh_solver->GetNodes()->SetBound_Disp(Point_Mesh,Target_Variable);

}
