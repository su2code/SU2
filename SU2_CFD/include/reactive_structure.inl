/*!
 * \file transport_model.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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


inline void       CReactive::InitializeMixture(CConfig *config) {}
inline bool       CReactive::Get_Ionization() {}
inline bool       CReactiveMutation::Get_Ionization() { return ionization; }
inline bool       CReactiveHardCode::Get_Ionization() { return ionization; }
inline vector<su2double> CReactive::Get_MolarMass() { }
inline vector<su2double> CReactiveMutation::Get_MolarMass(){ return Ms; }
inline vector<su2double> CReactiveHardCode::Get_MolarMass(){ return Ms; }
inline vector<su2double> CReactive::Get_CvTraRotSpecies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline vector<su2double> CReactive::Get_CvVibElSpecies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double* CReactive::Get_CpTraRotSpecies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double* CReactive::Get_CpVibElSpecies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double  CReactive::Get_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve){}
inline su2double  CReactive::Get_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){}
inline su2double  CReactive::Get_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve){}
inline su2double  CReactive::Get_MixtureEnergy(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline vector<su2double> CReactive::Get_MixtureEnergies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline vector<su2double> CReactive::Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {}
inline vector<su2double> CReactive::Get_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double CReactive::Get_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline su2double CReactive::Get_CVEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double  CReactive::Get_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline vector<su2double> CReactive::Get_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double* CReactive::Get_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double* CReactive::Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double  CReactive::Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline vector<su2double> CReactive::Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline vector<su2double> CReactive::Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve) {}
inline su2double  CReactive::Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve) {}
inline su2double  CReactive::Get_Density(su2double T, su2double *Xs, su2double P){}
