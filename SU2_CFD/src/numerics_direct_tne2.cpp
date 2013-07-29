/*!
 * \file numerics_direct_tne2.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>


CUpwRoe_TNE2::CUpwRoe_TNE2(unsigned short val_nDim, unsigned short val_nVar,
                           CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
  /*--- Read configuration parameters ---*/
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  ionization = config->GetIonization();
  
  /*--- Define useful constants ---*/
  nVar     = val_nVar;
  nDim     = val_nDim;
  nSpecies = val_nVar - val_nDim - 2;
  
  /*--- Allocate arrays ---*/
	Diff_U      = new double [nVar];
  Density_i   = new double [nSpecies];
  Density_j   = new double [nSpecies];
  RoeDensity  = new double [nSpecies];
  dPdrhos     = new double [nSpecies];
	Velocity_i  = new double [nDim];
	Velocity_j  = new double [nDim];
	RoeVelocity = new double [nDim];
  l           = new double [nDim];
  m           = new double [nDim];
	Lambda      = new double [nVar];
	Epsilon     = new double [nVar];
	P_Tensor    = new double* [nVar];
	invP_Tensor = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar] = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
	}
  Proj_flux_tensor_i = new double [nVar];
	Proj_flux_tensor_j = new double [nVar];
}

CUpwRoe_TNE2::~CUpwRoe_TNE2(void) {
	unsigned short iVar;
  
	delete [] Diff_U;
  delete [] Density_i;
  delete [] Density_j;
  delete [] RoeDensity;
  delete [] dPdrhos;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] RoeVelocity;
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
  delete [] l;
  delete [] m;
	delete [] Lambda;
	delete [] Epsilon;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
  
}

void CUpwRoe_TNE2::ComputeResidual(double *val_residual, double **val_Jacobian_i,
                                   double **val_Jacobian_j, CConfig *config) {
  
  // NOTE: VIBRATIONAL TEMPERATURES NOT BEING CALCULATED CORRECTLY!!!  Cv-VE NOT INCLUDING ELECTRONIC STATES!!!
  unsigned short iDim, iSpecies, iVar, jVar, kVar;
  double *hf, *Tref, *rotmodes, *molarmass, *charvibtemp;
  double DensityMix_i, DensityMix_j, DensityEnergyForm, DensityEnergyRef;
  double DensityCvtr, DensityCvve, conc;
  double dPdrhoE, dPdrhoEve;
  bool zero_order;
  
  /*--- Face area (norm or the normal vector) ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Unit Normal ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Load variables from the config class --*/
  hf        = config->GetEnthalpy_Formation();
  Tref      = config->GetRefTemperature();
  rotmodes  = config->GetRotationModes();
  molarmass = config->GetMolar_Mass();
  charvibtemp = config->GetCharVibTemp();
  
	/*--- Conserved variables at point i,
   Need to recompute SoundSpeed / Pressure / Enthalpy in
   case of 2nd order reconstruction ---*/
  DensityMix_i      = 0.0;
  DensityEnergyRef  = 0.0;
  DensityEnergyForm = 0.0;
  DensityCvtr       = 0.0;
  conc              = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Density_i[iSpecies] = U_i[iSpecies];
    DensityMix_i       += U_i[iSpecies];
    DensityEnergyForm  += U_i[iSpecies] * (hf[iSpecies] - UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies]);
    DensityEnergyRef   +=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    DensityCvtr        +=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               += U_i[iSpecies] / molarmass[iSpecies];
    if (Density_i[iSpecies] < 0.0) zero_order = true;
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    DensityEnergyRef   -=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    DensityCvtr        -=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               -= U_i[iSpecies] / molarmass[iSpecies];
  }
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[nSpecies+iDim] / DensityMix_i;
		sq_vel          += Velocity_i[iDim]*Velocity_i[iDim];
	}
  Temperature_i    = (U_i[nSpecies+nDim] - U_i[nSpecies+nDim+1] - DensityEnergyForm
                      + DensityEnergyRef - 0.5*DensityMix_i*sq_vel) / DensityCvtr;
  Temperature_ve_i = Temperature_i;
  
	Energy_i     = U_i[nSpecies+nDim] / DensityMix_i;
  Energy_ve_i  = U_i[nSpecies+nDim+1] / DensityMix_i;
  
  Pressure_i = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Pressure_i += U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_i;
  if (ionization) {
    iSpecies = nSpecies - 1;
    Pressure_i -= U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_i;
    Pressure_i += U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_ve_i;
  }
  if (Pressure_i < 0.0) zero_order = true;
  
  dPdrhoE = UNIVERSAL_GAS_CONSTANT/DensityCvtr * conc;
  SoundSpeed_i = sqrt((1.0 + dPdrhoE) * Pressure_i/DensityMix_i);
  Enthalpy_i = (U_i[nSpecies+nDim] + Pressure_i) / DensityMix_i;
  
  
  /*--- If it is not a physical solution, then
   use the zero order reconstruction ---*/
  if (zero_order) {
    for (iVar = 0; iVar < nVar; iVar++)
      U_i[iVar] = UZeroOrder_i[iVar];
    
    DensityMix_i      = 0.0;
    DensityEnergyRef  = 0.0;
    DensityEnergyForm = 0.0;
    DensityCvtr       = 0.0;
    conc              = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Density_i[iSpecies] = U_i[iSpecies];
      DensityMix_i       += U_i[iSpecies];
      DensityEnergyForm  += U_i[iSpecies] * (hf[iSpecies] - UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies]);
      DensityEnergyRef   +=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
      * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
      DensityCvtr        +=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
      * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
      conc               += U_i[iSpecies] / molarmass[iSpecies];
    }
    if (ionization) {
      iSpecies = nSpecies-1;
      DensityEnergyRef   -=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
      * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
      DensityCvtr        -=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
      * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
      conc               -= U_i[iSpecies] / molarmass[iSpecies];
    }
    sq_vel = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = U_i[nSpecies+iDim] / DensityMix_i;
      sq_vel          += Velocity_i[iDim]*Velocity_i[iDim];
    }
    Temperature_i    = (U_i[nSpecies+nDim] - U_i[nSpecies+nDim+1] - DensityEnergyForm
                        + DensityEnergyRef - 0.5*DensityMix_i*sq_vel) / DensityCvtr;
    Temperature_ve_i = Temperature_i;
    
    Energy_i     = U_i[nSpecies+nDim] / DensityMix_i;
    Energy_ve_i  = U_i[nSpecies+nDim+1] / DensityMix_i;
    
    Pressure_i = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Pressure_i += U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_i;
    if (ionization) {
      iSpecies = nSpecies - 1;
      Pressure_i -= U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_i;
      Pressure_i += U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_ve_i;
    }
    
    dPdrhoE = UNIVERSAL_GAS_CONSTANT/DensityCvtr * conc;
    SoundSpeed_i = sqrt((1.0 + dPdrhoE) * Pressure_i/DensityMix_i);
    Enthalpy_i = (U_i[nSpecies+nDim] + Pressure_i) / DensityMix_i;
  }
  
	/*--- Conserved variables at point j,
   Need to recompute SoundSpeed / Pressure / Enthalpy in
   case of 2nd order reconstruction ---*/
  zero_order        = false;
  DensityMix_j      = 0.0;
  DensityEnergyRef  = 0.0;
  DensityEnergyForm = 0.0;
  DensityCvtr       = 0.0;
  conc              = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Density_j[iSpecies] = U_j[iSpecies];
    DensityMix_i       += U_j[iSpecies];
    DensityEnergyForm  += U_j[iSpecies] * (hf[iSpecies] - UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies]);
    DensityEnergyRef   +=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    DensityCvtr        +=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               += U_j[iSpecies] / molarmass[iSpecies];
    if (Density_j[iSpecies] < 0.0) zero_order = true;
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    DensityEnergyRef   -=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    DensityCvtr        -=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               -= U_j[iSpecies] / molarmass[iSpecies];
  }
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[nSpecies+iDim] / DensityMix_j;
		sq_vel          += Velocity_j[iDim]*Velocity_j[iDim];
	}
  Temperature_j    = (U_j[nSpecies+nDim] - U_j[nSpecies+nDim+1] - DensityEnergyForm
                      + DensityEnergyRef - 0.5*DensityMix_j*sq_vel) / DensityCvtr;
  Temperature_ve_j = Temperature_j;
  
	Energy_j     = U_j[nSpecies+nDim] / DensityMix_j;
  Energy_ve_j  = U_j[nSpecies+nDim+1] / DensityMix_j;
  
  Pressure_j = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Pressure_j += U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_j;
  if (ionization) {
    iSpecies = nSpecies - 1;
    Pressure_j -= U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_j;
    Pressure_j += U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_ve_j;
  }
  if (Pressure_j < 0.0) zero_order = true;
  
  dPdrhoE = UNIVERSAL_GAS_CONSTANT/DensityCvtr * conc;
  SoundSpeed_j = sqrt((1.0 + dPdrhoE) * Pressure_j/DensityMix_j);
  Enthalpy_j = (U_j[nSpecies+nDim] + Pressure_j) / DensityMix_j;
  
  
  /*--- If it is not a physical solution, then
   use the zero order reconstruction ---*/
  if (zero_order) {
    for (iVar = 0; iVar < nVar; iVar++)
      U_i[iVar] = UZeroOrder_i[iVar];
    
    DensityMix_j      = 0.0;
    DensityEnergyRef  = 0.0;
    DensityEnergyForm = 0.0;
    DensityCvtr       = 0.0;
    conc              = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Density_j[iSpecies] = U_j[iSpecies];
      DensityMix_j       += U_j[iSpecies];
      DensityEnergyForm  += U_j[iSpecies] * (hf[iSpecies] - UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies]);
      DensityEnergyRef   +=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
      * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
      DensityCvtr        +=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
      * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
      conc               += U_j[iSpecies] / molarmass[iSpecies];
    }
    if (ionization) {
      iSpecies = nSpecies-1;
      DensityEnergyRef   -=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
      * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
      DensityCvtr        -=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
      * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
      conc               -= U_j[iSpecies] / molarmass[iSpecies];
    }
    sq_vel = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_j[iDim] = U_j[nSpecies+iDim] / DensityMix_j;
      sq_vel          += Velocity_j[iDim]*Velocity_j[iDim];
    }
    Temperature_j    = (U_j[nSpecies+nDim] - U_j[nSpecies+nDim+1] - DensityEnergyForm
                        + DensityEnergyRef - 0.5*DensityMix_j*sq_vel) / DensityCvtr;
    Temperature_ve_j = Temperature_j;
    
    Energy_j     = U_j[nSpecies+nDim] / DensityMix_j;
    Energy_ve_j  = U_j[nSpecies+nDim+1] / DensityMix_j;
    
    Pressure_j = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Pressure_j += U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_j;
    if (ionization) {
      iSpecies = nSpecies - 1;
      Pressure_j -= U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_j;
      Pressure_j += U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_ve_j;
    }
    
    dPdrhoE = UNIVERSAL_GAS_CONSTANT/DensityCvtr * conc;
    SoundSpeed_j = sqrt((1.0 + dPdrhoE) * Pressure_j/DensityMix_j);
    Enthalpy_j = (U_j[nSpecies+nDim] + Pressure_j) / DensityMix_j;
  }
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  dPdrhoE = 0.0;
  dPdrhoEve = 0.0;
  conc        = 0.0;
  DensityCvtr = 0.0;
  DensityCvve = 0.0;
  R = sqrt(abs(DensityMix_j/DensityMix_i));
  RoeTemperature_ve = (R*Temperature_ve_j + Temperature_ve_i) / (R+1);
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    RoeDensity[iSpecies] = (R*(Density_j[iSpecies]/DensityMix_j)
                            + (Density_i[iSpecies]/DensityMix_i)) / (R+1);
    DensityCvtr        +=  RoeDensity[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    DensityCvve        +=  UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * (charvibtemp[iSpecies]/RoeTemperature_ve)
    * (charvibtemp[iSpecies]/RoeTemperature_ve) * exp(charvibtemp[iSpecies]/RoeTemperature_ve)
    / ((exp(charvibtemp[iSpecies]/RoeTemperature_ve) - 1.0)*(exp(charvibtemp[iSpecies]/RoeTemperature_ve) - 1.0));
    conc               += RoeDensity[iSpecies] / molarmass[iSpecies];
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    DensityCvtr        -=  RoeDensity[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               -= RoeDensity[iSpecies] / molarmass[iSpecies];
  }
  dPdrhoE   = UNIVERSAL_GAS_CONSTANT/DensityCvtr * conc;
  dPdrhoEve =
  sq_vel = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i) / (R+1);
  RoeEnergy_ve = (R*Energy_ve_j + Energy_ve_i) / (R+1);
  RoePressure = (R*Pressure_j + Pressure_i) / (R+1);
  RoeSoundSpeed = sqrt((1.0 + dPdrhoE) * RoePressure/(DensityMix_i*R));
  
  /*--- Compute Proj_flux_tensor_i ---*/
  GetInviscidProjFlux(Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, &Energy_ve_i, Normal, Proj_flux_tensor_i);
  
  /*--- Compute Proj_flux_tensor_j ---*/
  GetInviscidProjFlux(Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, &Energy_ve_j, Normal, Proj_flux_tensor_j);
  
  /*--- Compute projected P, invP, and Lambda ---*/
  GetPMatrix(RoeDensity, RoeVelocity, &RoeEnthalpy, &RoeEnergy_ve, &RoeSoundSpeed, dPdrhos,
             dPdrhoE, dPdrhoEve, UnitaryNormal, l, m, P_Tensor);
  GetPMatrix_inv(RoeDensity, RoeVelocity, &RoeEnergy_ve, &RoeSoundSpeed, dPdrhos,
                 dPdrhoE, dPdrhoEve, UnitaryNormal, l, m, invP_Tensor);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitaryNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
  }
  
  /*--- Flow eigenvalues and entropy correctors ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
  //	/*--- Harten and Hyman (1983) entropy correction ---*/
  //	for (iDim = 0; iDim < nDim; iDim++)
  //		Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,ProjVelocity_j-Lambda[iDim]));
  //
  //	Epsilon[nVar-2] = 4.0*max(0.0, max(Lambda[nVar-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar-2]));
  //	Epsilon[nVar-1] = 4.0*max(0.0, max(Lambda[nVar-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar-1]));
  //
  //	for (iVar = 0; iVar < nVar; iVar++)
  //		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
  //			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
  //		else
  //			Lambda[iVar] = fabs(Lambda[iVar]);
  
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  /*--- Jacobians of the inviscid flux, scaled by
   0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
  GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
  GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
  
  /*--- Diference variables iPoint and jPoint ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];
  
  /*--- Roe's Flux approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_ij = 0.0;
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
      val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
      val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
    }
  }
}

CUpwAUSM_TNE2::CUpwAUSM_TNE2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Diff_U = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	RoeVelocity = new double [nDim];
	delta_vel  = new double [nDim];
	delta_wave = new double [nVar];
	Proj_flux_tensor_i = new double [nVar];
	Proj_flux_tensor_j = new double [nVar];
	Lambda = new double [nVar];
	Epsilon = new double [nVar];
	P_Tensor = new double* [nVar];
	invP_Tensor = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar] = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
	}
}

CUpwAUSM_TNE2::~CUpwAUSM_TNE2(void) {
	unsigned short iVar;
  
	delete [] Diff_U;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] RoeVelocity;
	delete [] delta_vel;
	delete [] delta_wave;
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
	delete [] Lambda;
	delete [] Epsilon;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
  
}

void CUpwAUSM_TNE2::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	/*--- Point i, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	Density_i = U_i[0];
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
	}
	Energy_i = U_i[nDim+1] / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel));
	Pressure_i = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	Enthalpy_i = (U_i[nDim+1] + Pressure_i) / Density_i;
  
	/*--- Point j, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	Density_j = U_j[0];
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_j = U_j[nDim+1] / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel));
	Pressure_j = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	Enthalpy_j = (U_j[nDim+1] + Pressure_j) / Density_j;
  
	/*--- Projected velocities ---*/
	ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
	}
  
	double mL	= ProjVelocity_i/SoundSpeed_i;
	double mR	= ProjVelocity_j/SoundSpeed_j;
  
	double mLP;
	if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0);
	else mLP = 0.5*(mL+fabs(mL));
  
	double mRM;
	if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0);
	else mRM = 0.5*(mR-fabs(mR));
  
	double mF = mLP + mRM;
  
	double pLP;
	if (fabs(mL) <= 1.0) pLP = 0.25*Pressure_i*(mL+1.0)*(mL+1.0)*(2.0-mL);
	else pLP = 0.5*Pressure_i*(mL+fabs(mL))/mL;
  
	double pRM;
	if (fabs(mR) <= 1.0) pRM = 0.25*Pressure_j*(mR-1.0)*(mR-1.0)*(2.0+mR);
	else pRM = 0.5*Pressure_j*(mR-fabs(mR))/mR;
  
	double pF = pLP + pRM;
	double Phi = fabs(mF);
  
	val_residual[0] = 0.5*(mF*((Density_i*SoundSpeed_i)+(Density_j*SoundSpeed_j))-Phi*((Density_j*SoundSpeed_j)-(Density_i*SoundSpeed_i)));
	for (iDim = 0; iDim < nDim; iDim++)
		val_residual[iDim+1] = 0.5*(mF*((Density_i*SoundSpeed_i*Velocity_i[iDim])+(Density_j*SoundSpeed_j*Velocity_j[iDim]))
                                -Phi*((Density_j*SoundSpeed_j*Velocity_j[iDim])-(Density_i*SoundSpeed_i*Velocity_i[iDim])))+UnitaryNormal[iDim]*pF;
	val_residual[nVar-1] = 0.5*(mF*((Density_i*SoundSpeed_i*Enthalpy_i)+(Density_j*SoundSpeed_j*Enthalpy_j))-Phi*((Density_j*SoundSpeed_j*Enthalpy_j)-(Density_i*SoundSpeed_i*Enthalpy_i)));
  
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] *= Area;
  
  
  
	if (implicit) {
    
		/*--- Promediate Roe variables iPoint and jPoint ---*/
		R = sqrt(Density_j/Density_i);
		RoeDensity = R*Density_i;
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
			sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
		}
		RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
		RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));
    
		/*--- Compute P and Lambda (do it with the Normal) ---*/
		GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, P_Tensor);
    
		ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity   += RoeVelocity[iDim]*UnitaryNormal[iDim];
			ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
			ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
		}
    
		/*--- Flow eigenvalues and Entropy correctors ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			Lambda[iDim] = ProjVelocity;
		Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
		Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
		/*--- Compute inverse P ---*/
		GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, invP_Tensor);
    
		/*--- Jacobias of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
		GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
		/*--- Diference variables iPoint and jPoint ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    
		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				Proj_ModJac_Tensor_ij = 0.0;
				/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
				for (kVar = 0; kVar < nVar; kVar++)
					Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
				val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
				val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
			}
		}
	}
}


CCentLax_TNE2::CCentLax_TNE2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  /*--- Read configuration parameters ---*/
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  ionization = config->GetIonization();
  
  /*--- Define useful constants ---*/
  nVar     = val_nVar;
  nDim     = val_nDim;
  nSpecies = val_nVar - val_nDim - 2;
  
	/*--- Artifical dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_Flow();
  
	/*--- Allocate some structures ---*/
	Diff_U           = new double [nVar];
  Density_i        = new double[nSpecies];
  Density_j        = new double[nSpecies];
  MeanDensity      = new double[nSpecies];
	Velocity_i       = new double [nDim];
	Velocity_j       = new double [nDim];
	MeanVelocity     = new double [nDim];
	Proj_flux_tensor = new double [nVar];
  
}

CCentLax_TNE2::~CCentLax_TNE2(void) {
	delete [] Diff_U;
  delete [] Density_i;
  delete [] Density_j;
  delete [] MeanDensity;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] MeanVelocity;
	delete [] Proj_flux_tensor;
  
}

void CCentLax_TNE2::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i,
                                    double **val_Jacobian_j, CConfig *config) {
  
  // NOTE: Vibrational-electronic temperature not calculated correctly!!!
  // NOTE: Electronic energy not included in Cvve!!!
  // NOTE: Calculating visc contribution to Jacobian using mean values -- should use i & j!
  unsigned short iSpecies, iDim;
  double DensityMix_i, DensityMix_j;
  double *hf, *Tref, *rotmodes, *molarmass, *charvibtemp;
  double conc, DensityEnergyRef, DensityEnergyForm, DensityCvtr, DensityCvve;
  double energy, energy_ve, energy_v, energy_e, energy_formation, sqvel;
  
  /*--- Load variables from the config class --*/
  hf        = config->GetEnthalpy_Formation();
  Tref      = config->GetRefTemperature();
  rotmodes  = config->GetRotationModes();
  molarmass = config->GetMolar_Mass();
  charvibtemp = config->GetCharVibTemp();
  
	/*--- Conserved variables at point i,
   Need to recompute SoundSpeed / Pressure / Enthalpy in
   case of 2nd order reconstruction ---*/
  DensityMix_i      = 0.0;
  DensityEnergyRef  = 0.0;
  DensityEnergyForm = 0.0;
  DensityCvtr       = 0.0;
  conc              = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Density_i[iSpecies] = U_i[iSpecies];
    DensityMix_i       += U_i[iSpecies];
    DensityEnergyForm  += U_i[iSpecies] * (hf[iSpecies] - UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies]);
    DensityEnergyRef   +=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    DensityCvtr        +=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               += U_i[iSpecies] / molarmass[iSpecies];
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    DensityEnergyRef   -=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    DensityCvtr        -=  U_i[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               -= U_i[iSpecies] / molarmass[iSpecies];
  }
	sqvel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[nSpecies+iDim] / DensityMix_i;
		sqvel          += Velocity_i[iDim]*Velocity_i[iDim];
	}
  Temperature_i    = (U_i[nSpecies+nDim] - U_i[nSpecies+nDim+1] - DensityEnergyForm
                      + DensityEnergyRef - 0.5*DensityMix_i*sqvel) / DensityCvtr;
  Temperature_ve_i = Temperature_i;
	Energy_i     = U_i[nSpecies+nDim] / DensityMix_i;
  Energy_ve_i  = U_i[nSpecies+nDim+1] / DensityMix_i;
  Pressure_i = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Pressure_i += U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_i;
  if (ionization) {
    iSpecies = nSpecies - 1;
    Pressure_i -= U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_i;
    Pressure_i += U_i[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_ve_i;
  }
  dPdrhoE = UNIVERSAL_GAS_CONSTANT/DensityCvtr * conc;
  SoundSpeed_i = sqrt((1.0 + dPdrhoE) * Pressure_i/DensityMix_i);
  Enthalpy_i = (U_i[nSpecies+nDim] + Pressure_i) / DensityMix_i;
  
  /*--- Conserved variables at point j,
   Need to recompute SoundSpeed / Pressure / Enthalpy in
   case of 2nd order reconstruction ---*/
  DensityMix_j      = 0.0;
  DensityEnergyRef  = 0.0;
  DensityEnergyForm = 0.0;
  DensityCvtr       = 0.0;
  conc              = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Density_j[iSpecies] = U_j[iSpecies];
    DensityMix_j       += U_j[iSpecies];
    DensityEnergyForm  += U_j[iSpecies] * (hf[iSpecies] - UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies]);
    DensityEnergyRef   +=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    DensityCvtr        +=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               += U_j[iSpecies] / molarmass[iSpecies];
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    DensityEnergyRef   -=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    DensityCvtr        -=  U_j[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc               -= U_j[iSpecies] / molarmass[iSpecies];
  }
	sqvel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[nSpecies+iDim] / DensityMix_j;
		sqvel          += Velocity_j[iDim]*Velocity_j[iDim];
	}
  Temperature_j    = (U_j[nSpecies+nDim] - U_j[nSpecies+nDim+1] - DensityEnergyForm
                      + DensityEnergyRef - 0.5*DensityMix_j*sqvel) / DensityCvtr;
  Temperature_ve_j = Temperature_j;
	Energy_j     = U_j[nSpecies+nDim] / DensityMix_j;
  Energy_ve_j  = U_j[nSpecies+nDim+1] / DensityMix_j;
  Pressure_j = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Pressure_j += U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_j;
  if (ionization) {
    iSpecies = nSpecies - 1;
    Pressure_j -= U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_j;
    Pressure_j += U_j[iSpecies] * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Temperature_ve_j;
  }
  dPdrhoE = UNIVERSAL_GAS_CONSTANT/DensityCvtr * conc;
  SoundSpeed_j = sqrt((1.0 + dPdrhoE) * Pressure_j/DensityMix_j);
  Enthalpy_j = (U_j[nSpecies+nDim] + Pressure_j) / DensityMix_j;
  
	/*--- Compute mean values of the variables ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MeanDensity[iSpecies] = 0.5 * (Density_i[iSpecies]+Density_j[iSpecies]);
	MeanPressure = 0.5 * (Pressure_i+Pressure_j);
	MeanEnthalpy = 0.5 * (Enthalpy_i+Enthalpy_j);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanVelocity[iDim] =  0.5 * (Velocity_i[iDim]+Velocity_j[iDim]);
	MeanEnergy         = 0.5 * (Energy_i+Energy_j);
  MeanEnergy_ve      = 0.5 * (Energy_ve_i+Energy_ve_j);
  MeanTemperature_ve = 0.5 * (Temperature_ve_i + Temperature_ve_j);
  MeanTemperature    = 0.5 * (Temperature_i + Temperature_j);
  
  /*--- Compute pressure derivatives using mean variables ---*/
  DensityCvtr       = 0.0;
  DensityCvve       = 0.0;
  conc              = 0.0;
  dPdrhoEve         = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    DensityCvtr +=  MeanDensity[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    DensityCvve +=  UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * (charvibtemp[iSpecies]/MeanTemperature_ve)
    * (charvibtemp[iSpecies]/MeanTemperature_ve) * exp(charvibtemp[iSpecies]/MeanTemperature_ve)
    / ((exp(charvibtemp[iSpecies]/MeanTemperature_ve) - 1.0)*
       (exp(charvibtemp[iSpecies]/MeanTemperature_ve) - 1.0));
    conc        +=  MeanDensity[iSpecies] / molarmass[iSpecies];
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    DensityCvtr -=  MeanDensity[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc        -=  MeanDensity[iSpecies] / molarmass[iSpecies];
    DensityCvve -=  UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * (charvibtemp[iSpecies]/MeanTemperature_ve)
    * (charvibtemp[iSpecies]/MeanTemperature_ve) * exp(charvibtemp[iSpecies]/MeanTemperature_ve)
    / ((exp(charvibtemp[iSpecies]/MeanTemperature_ve) - 1.0)*
       (exp(charvibtemp[iSpecies]/MeanTemperature_ve) - 1.0));
    dPdrhoEve    = UNIVERSAL_GAS_CONSTANT/DensityCvve * MeanDensity[iSpecies]/molarmass[iSpecies];
  }
  dPdrhoE = UNIVERSAL_GAS_CONSTANT/DensityCvtr * conc;
  dPdrhoEve -= dPdrhoE;
  sqvel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += MeanVelocity[iDim]*MeanVelocity[iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    energy_v  =  UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]
    * charvibtemp[iSpecies] / ( exp(charvibtemp[iSpecies]/MeanTemperature_ve) -1.0 );
    energy_e = 0.0;
    energy_formation =  hf[iSpecies] - UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    energy =  (3.0/2.0 + rotmodes[iSpecies]/2.0) * (MeanTemperature - Tref[iSpecies])
    + energy_v + energy_e + energy_formation + 0.5*sqvel;
    energy_ve = energy_v + energy_e;
    dPdrhos[iSpecies] = UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * MeanTemperature
    + dPdrhoE*0.5*sqvel - dPdrhoE*energy - dPdrhoEve*energy_ve;
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    dPdrhos[iSpecies] -= UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * MeanTemperature;
    dPdrhos[iSpecies] += UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * MeanTemperature_ve;
  }
  
	/*--- Get projected flux tensor ---*/
	GetInviscidProjFlux(MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy,
                      &MeanEnergy_ve, Normal, Proj_flux_tensor);
  
	/*--- Compute inviscid residual ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}
  
	/*--- Jacobians of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	if (implicit) {
    GetInviscidProjJac(MeanDensity, MeanVelocity, &MeanEnthalpy, &MeanEnergy_ve,
                       dPdrhos, dPdrhoE, dPdrhoEve, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}
  
	/*--- Computes differences btw. conservative variables ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Diff_U[iVar] = U_i[iVar]-U_j[iVar];
	Diff_U[nDim+1] = DensityMix_i*Enthalpy_i-DensityMix_j*Enthalpy_j;
  
	/*--- Compute the local spectral radius and the stretching factor ---*/
	ProjVelocity_i = 0; ProjVelocity_j = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	sc0 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_0 = Param_Kappa_0*sc0*double(nDim)/3.0;
  
	/*--- Compute viscous part of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resvisc[iVar] = Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
	}
  
	if (implicit) {
		cte = Epsilon_0*StretchingFactor*MeanLambda;
    
		for (iVar = 0; iVar < (nVar-1); iVar++) {
			val_Jacobian_i[iVar][iVar] += cte;
			val_Jacobian_j[iVar][iVar] -= cte;
		}
    
		/*--- Last rows: CAREFUL!! You have differences of \rho_Enthalpy, not differences of \rho_Energy ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Jacobian_i[nSpecies+nDim][iSpecies] += cte*dPdrhos[iSpecies];
		for (iDim = 0; iDim < nDim; iDim++)
			val_Jacobian_i[nSpecies+nDim][nSpecies+iDim] -= cte*UNIVERSAL_GAS_CONSTANT/DensityCvtr*conc * Velocity_i[iDim];
		val_Jacobian_i[nSpecies+nDim][nSpecies+nDim] += cte*(1+dPdrhoE);
    val_Jacobian_i[nSpecies+nDim][nSpecies+nDim+1] += cte*dPdrhoEve;
    val_Jacobian_i[nSpecies+nDim+1][nSpecies+nDim+1] += cte;
    
		/*--- Last row of Jacobian_j ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Jacobian_j[nSpecies+nDim][iSpecies] -= cte*dPdrhos[iSpecies];
		for (iDim = 0; iDim < nDim; iDim++)
			val_Jacobian_j[nSpecies+nDim][nSpecies+iDim] += cte*UNIVERSAL_GAS_CONSTANT/DensityCvtr*conc * Velocity_j[iDim];
		val_Jacobian_j[nSpecies+nDim][nSpecies+nDim] -= cte*(1+dPdrhoE);
    val_Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] -= cte*dPdrhoEve;
    val_Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] -= cte;
	}
  cout << "CCent_Lax Test!" << endl;
  cin.get();
}
