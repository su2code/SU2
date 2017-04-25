/*!
 * \file transfer_structure.inl
 * \brief In-Line subroutines of the <i>transfer_structure.hpp</i> file.
 * \author R. Sanchez
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

inline void CTransfer::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                      CConfig *target_config, unsigned long Marker_Target,
                      unsigned long Vertex_Target, unsigned long Point_Target) { }
