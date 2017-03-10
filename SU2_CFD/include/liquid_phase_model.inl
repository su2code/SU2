/*!
 * \file liquid_phase_model.inl
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

inline su2double CLiquidModel::GetLiquidDensity_PT() { return rho_l; }
inline su2double CLiquidModel::GetLiquidEnthalpy_PT() { return h_l; }
inline su2double CLiquidModel::GetSurfaceTension_T() { return sigma; }
inline su2double CLiquidModel::GetTsat_P() { return Tsat; }
inline su2double CLiquidModel::GetPsat_T() { return Psat; }
inline su2double CLiquidModel::GetdGibbs_PT() { return dGibbs; }

inline su2double CLiquidModel::SetLiquidDensity_PT(su2double P, su2double T, su2double Tstar) { return rho_l; }
inline su2double CLiquidModel::SetLiquidEnthalpy_Prho(su2double P, su2double T, su2double h_v, su2double Tstar) { return h_l; }
inline su2double CLiquidModel::SetSurfaceTension_T(su2double T, su2double Tstar) { return sigma; }
inline su2double CLiquidModel::SetTsat_P(su2double P) { return Tsat; }
inline su2double CLiquidModel::SetdGibbs_PT(su2double P, su2double T, su2double Gas_Constant) { return dGibbs, Psat; }