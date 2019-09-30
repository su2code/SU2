/*!
 * \file transport_model.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
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

inline su2double CViscosityModel::GetViscosity() { return Mu; }
inline su2double CViscosityModel::Getdmudrho_T () { return dmudrho_T; }
inline su2double CViscosityModel::GetdmudT_rho() { return dmudT_rho; }
inline void CViscosityModel::SetViscosity(su2double T, su2double rho) {}
inline void CViscosityModel::SetDerViscosity(su2double T, su2double rho) {}

inline su2double CConductivityModel::GetConductivity() { return Kt; }
inline su2double CConductivityModel::Getdktdrho_T () { return dktdrho_T; }
inline su2double CConductivityModel::GetdktdT_rho () { return dktdT_rho; }
inline void CConductivityModel::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {}
inline void CConductivityModel::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {}
