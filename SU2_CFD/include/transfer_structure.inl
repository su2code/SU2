/*!
 * \file transfer_structure.inl
 * \brief In-Line subroutines of the <i>transfer_structure.hpp</i> file.
 * \author R. Sanchez
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

#pragma once

inline void CTransfer::GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                                                    CGeometry *donor_geometry, CGeometry *target_geometry,
                       CConfig *donor_config, CConfig *target_config) { }

inline void CTransfer::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, 
                          CConfig *donor_config, unsigned long Marker_Donor, 
                     unsigned long Vertex_Donor, unsigned long Point_Donor) { }

inline void CTransfer::InitializeTarget_Variable(CSolver *target_solution, unsigned long Marker_Target,
                                                 unsigned long Vertex_Target, unsigned short nDonorPoints){
  for (unsigned short iVar = 0; iVar < nVar; iVar++) Target_Variable[iVar] = 0.0;
}

inline void CTransfer::RecoverTarget_Variable(long indexPoint_iVertex, su2double *Buffer_Bcast_Variables,
                                              su2double donorCoeff){
  for (unsigned short iVar = 0; iVar < nVar; iVar++) Target_Variable[iVar] += donorCoeff * Buffer_Bcast_Variables[indexPoint_iVertex*nVar+iVar];

}

inline void CTransfer::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
										  CConfig *target_config, unsigned long Marker_Target,
										  unsigned long Vertex_Target, unsigned long Point_Target) { }

inline void CTransfer::SetAverageValues(CSolver *donor_solution, CSolver *target_solution, unsigned short donorZone) { }

inline void CTransfer::SetAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry, unsigned short donorZone) { }

inline void CTransfer::SetSpanWiseLevels(CConfig *donor_config, CConfig *target_config) { }
