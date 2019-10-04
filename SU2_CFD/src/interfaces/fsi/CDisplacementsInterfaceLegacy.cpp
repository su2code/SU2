/*!
 * \file CDisplacementsInterfaceLegacy.cpp
 * \brief Declaration and inlines of the class to transfer structural displacements
 *        from a structural zone into a fluid zone.
 * \author Ruben Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../../include/interfaces/fsi/CDisplacementsInterfaceLegacy.hpp"

CDisplacementsInterfaceLegacy::CDisplacementsInterfaceLegacy(void) : CInterface() {

}

CDisplacementsInterfaceLegacy::CDisplacementsInterfaceLegacy(unsigned short val_nVar,
                                                             unsigned short val_nConst, CConfig *config) :
  CInterface(val_nVar, val_nConst, config) {

}

CDisplacementsInterfaceLegacy::~CDisplacementsInterfaceLegacy(void) {

}


void CDisplacementsInterfaceLegacy::GetPhysical_Constants(CSolver *struct_solution, CSolver *flow_solution,
                                                          CGeometry *struct_geometry, CGeometry *flow_geometry,
                                                          CConfig *struct_config, CConfig *flow_config) {
}

void CDisplacementsInterfaceLegacy::GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry,
                                                      CConfig *struct_config, unsigned long Marker_Struct,
                                                      unsigned long Vertex_Struct, unsigned long Point_Struct) {

  su2double *DisplacementDonor, *DisplacementDonor_Prev;
  unsigned short iVar;

  /*--- The displacements come from the predicted solution ---*/
  DisplacementDonor = struct_solution->GetNodes()->GetSolution_Pred(Point_Struct);

  DisplacementDonor_Prev = struct_solution->GetNodes()->GetSolution_Pred_Old(Point_Struct);

  for (iVar = 0; iVar < nVar; iVar++)
    Donor_Variable[iVar] = DisplacementDonor[iVar] - DisplacementDonor_Prev[iVar];
}

void CDisplacementsInterfaceLegacy::SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
                                                       CConfig *flow_config, unsigned long Marker_Flow,
                                                       unsigned long Vertex_Flow, unsigned long Point_Flow) {

  flow_geometry->vertex[Marker_Flow][Vertex_Flow]->SetVarCoord(Target_Variable);
}
