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

inline su2double CLiquidModel::GetLiquidDensity() { return rho_l; }
inline su2double CLiquidModel::GetMixtureDensity() { return rho_m; }
inline su2double CLiquidModel::GetLiquidTemperature() { return T_l; }
inline su2double CLiquidModel::GetLiquidEnthalpy() { return h_l; }
inline su2double CLiquidModel::GetSurfaceTension() { return sigma; }
inline su2double CLiquidModel::GetTsat() { return Tsat; }
inline su2double CLiquidModel::GetPsat() { return Psat; }
inline su2double CLiquidModel::GetCriticalRadius() { return Rc; }
inline su2double CLiquidModel::GetRadius() { return R; }

inline void CLiquidModel::SetLiquidProp(su2double P, su2double T, su2double rho, su2double h_v, su2double *Two_Phase_Var) {}
inline void CLiquidModel::SetTsat(su2double P) {}  
inline void CLiquidModel::SetPsat (su2double T) {}  
inline void CLiquidModel::SetLiquidDensity() {}
inline void CLiquidModel::SetTLiquid( su2double T) {}
inline void CLiquidModel::SetLiquidEnthalpy(su2double h_v) {} 
inline void CLiquidModel::SetSurfaceTension(su2double T) {}
inline void CLiquidModel::SetRadius(su2double *Two_Phase_Var) {}
inline void CLiquidModel::SetRCritical (su2double P, su2double T) {}
inline void CLiquidModel::SetDensity_Mixture (su2double rho, su2double *Two_Phase_Var) {}