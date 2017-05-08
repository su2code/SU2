/*!
 * \file transport_model.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
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

inline su2double CViscosityModel::GetViscosity() { return Mu; }
inline su2double CViscosityModel::Getdmudrho_T () { return dmudrho_T; }
inline su2double CViscosityModel::GetdmudT_rho() { return dmudT_rho; }
inline void CViscosityModel::SetViscosity(su2double T, su2double rho) {}
inline void CViscosityModel::SetDerViscosity(su2double T, su2double rho) {}

inline su2double CConductivityModel::GetConductivity() { return Kt; }
inline su2double CConductivityModel::Getdktdrho_T () { return dktdrho_T; }
inline su2double CConductivityModel::GetdktdT_rho () { return dktdT_rho; }
inline void CConductivityModel::SetConductivity(su2double T, su2double rho, su2double mu, su2double cp) {}
inline void CConductivityModel::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {}
