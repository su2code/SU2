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

inline su2double CLiquidModel::Get_LiquidDensity() { return rho_l; }
inline su2double CLiquidModel::Get_MixtureDensity() { return rho_m; }
inline su2double CLiquidModel::Get_LiquidTemperature() { return T_l; }
inline su2double CLiquidModel::Get_LiquidEnthalpy() { return h_l; }
inline su2double CLiquidModel::Get_SurfaceTension() { return sigma; }
inline su2double CLiquidModel::Get_Tsat() { return Tsat; }
inline su2double CLiquidModel::Get_Psat() { return Psat; }
inline su2double CLiquidModel::Get_CriticalRadius() { return Rc; }
inline su2double CLiquidModel::Get_Radius() { return R; }

// inline void CLiquidModel::Set_LiquidProp(su2double P, su2double T, su2double rho, su2double h_v, su2double Rcritical, su2double Radius, su2double mom3) {}
inline void CLiquidModel::SetTsat(su2double P)  {}
inline void CLiquidModel::SetPsat (su2double T) {} 
inline void CLiquidModel::SetLiquidDensity()    {}
inline void CLiquidModel::SetTLiquid( su2double T, su2double Rcritical, su2double Radius)  {}
inline void CLiquidModel::SetLiquidEnthalpy(su2double h_v) {}
inline void CLiquidModel::SetSurfaceTension(su2double T, su2double Rdroplet)   {}
// inline void CLiquidModel::SetRadius(su2double *Two_Phase_Var) {}
// inline void CLiquidModel::SetRCritical (su2double P, su2double T)  {}
// inline void CLiquidModel::SetDensity_Mixture (su2double rho, su2double mom3)  {}