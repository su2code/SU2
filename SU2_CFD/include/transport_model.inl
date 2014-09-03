/*!
 * \file transport_model.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna
 * \version 3.2.1 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

inline double CViscosityModel::GetViscosity() { return Mu; }
inline double CViscosityModel::Getdmudrho_T () { return dmudrho_T; }
inline double CViscosityModel::GetdmudT_rho() { return dmudT_rho; }
inline void CViscosityModel::SetViscosity(double T, double rho) {}

inline double CThermalConductivityModel::GetThermalConductivity() { return Kt; }
inline double CThermalConductivityModel::GetDerThermalConductivity_rho_T () { return dktdrho_T; }
inline double CThermalConductivityModel::GetDerThermalConductivity_T_rho () { return dktdT_rho; }
inline void CThermalConductivityModel::SetThermalConductivity(double par1, double par2) {}
