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
  
  unsigned short iDim, iSpecies, iVar, jVar, kVar, nHeavy, nEl;
  double *Ms;
  double DensityMix_i, DensityMix_j, conc, conc_i, conc_j;
  double dPdrhoE_i, dPdrhoE_j, dPdrhoEve_i, dPdrhoEve_j, dPdrhoE, dPdrhoEve;
  double Ru, rhoCvtr_i, rhoCvtr_j, rhoCvve_i, rhoCvve_j, rho_el_i, rho_el_j, rho_el;
  
  /*--- Face area (norm or the normal vector) ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Unit Normal ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Determine the number of heavy particle species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, c]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Density_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    Density_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[VEL_INDEX+iDim];
    Velocity_j[iDim] = V_j[VEL_INDEX+iDim];
  }
  Pressure_i       = V_i[P_INDEX];
  Pressure_j       = V_j[P_INDEX];
  Enthalpy_i       = V_i[H_INDEX];
  Enthalpy_j       = V_j[H_INDEX];
  SoundSpeed_i     = V_i[A_INDEX];
  SoundSpeed_j     = V_j[A_INDEX];
  DensityMix_i     = V_i[RHO_INDEX];
  DensityMix_j     = V_j[RHO_INDEX];
  rhoCvtr_i        = V_i[RHOCVTR_INDEX];
  rhoCvtr_j        = V_j[RHOCVTR_INDEX];
  rhoCvve_i        = V_i[RHOCVVE_INDEX];
  rhoCvve_j        = V_j[RHOCVVE_INDEX];
  Temperature_i    = V_i[T_INDEX];
  Temperature_j    = V_j[T_INDEX];
  Temperature_ve_i = V_i[TVE_INDEX];
  Temperature_ve_j = V_j[TVE_INDEX];
  
  /*--- Read parameters from CConfig ---*/
  Ms = config->GetMolar_Mass();
  
  /*--- Rename for convenience ---*/
  Ru = UNIVERSAL_GAS_CONSTANT;
  
  /*--- Calculate Roe variables ---*/
  R    = sqrt(abs(DensityMix_j/DensityMix_i));
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    RoeDensity[iSpecies] = (R*Density_j[iSpecies] + Density_i[iSpecies]) / (R+1);
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim] + Velocity_i[iDim]) / (R+1);
  }
  RoeEnthalpy       = (R*Enthalpy_j       + Enthalpy_i)       / (R+1);
  RoeEnergy_ve      = (R*Energy_ve_j      + Energy_ve_i)      / (R+1);
  RoeTemperature_ve = (R*Temperature_ve_j + Temperature_ve_i) / (R+1);
  RoePressure       = (R*Pressure_j       + Pressure_i)       / (R+1);
  
  if (ionization) {
    rho_el   = RoeDensity[nSpecies-1];
    rho_el_i = Density_i[nSpecies-1];
    rho_el_j = Density_j[nSpecies-1];
  } else {
    rho_el   = 0.0;
    rho_el_i = 0.0;
    rho_el_j = 0.0;
  }

  /*--- Calculate quantities using Roe variables ---*/
  conc = 0.0;
  conc_i = 0.0;
  conc_j = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    conc_i += Density_i[iSpecies] / Ms[iSpecies];
    conc_j += Density_j[iSpecies] / Ms[iSpecies];
    conc   += RoeDensity[iSpecies] / Ms[iSpecies];
    dPdrhos[iSpecies] = (R*dPdrhos_j[iSpecies] + dPdrhos_i[iSpecies]) / (R+1);
  }
  dPdrhoE_i   = conc_i*Ru / rhoCvtr_i;
  dPdrhoE_j   = conc_j*Ru / rhoCvtr_j;
  dPdrhoEve_i = -dPdrhoE_i + rho_el_i * Ru/Ms[nSpecies-1] * 1.0/rhoCvve_i;
  dPdrhoEve_j = -dPdrhoE_j + rho_el_j * Ru/Ms[nSpecies-1] * 1.0/rhoCvve_j;
  
  dPdrhoE   = Ru*(R+1) / (R*rhoCvtr_j+rhoCvtr_i) * conc;
  dPdrhoEve = -dPdrhoE + rho_el * Ru/Ms[nSpecies-1] * (R+1)/(R*rhoCvve_j+rhoCvve_i);
  
  RoeSoundSpeed = sqrt((1.0 + dPdrhoE) * RoePressure/(DensityMix_i*R));
  
  /*--- Calculate dual grid tangent vectors for P & invP ---*/
  CreateBasis(UnitaryNormal);
  
  /*--- Compute Proj_flux_tensor_i ---*/
  GetInviscidProjFlux(Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, &Energy_ve_i, Normal, Proj_flux_tensor_i);
  
  /*--- Compute Proj_flux_tensor_j ---*/
  GetInviscidProjFlux(Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, &Energy_ve_j, Normal, Proj_flux_tensor_j);
  
  /*--- Compute projected P, invP, and Lambda ---*/
  GetPMatrix(RoeDensity, RoeVelocity, &RoeEnthalpy, &RoeEnergy_ve, &RoeSoundSpeed, dPdrhos,
             dPdrhoE, dPdrhoEve, UnitaryNormal, l, m, P_Tensor);
  GetPMatrix_inv(RoeDensity, RoeVelocity, &RoeEnergy_ve, &RoeSoundSpeed, dPdrhos,
                 dPdrhoE, dPdrhoEve, UnitaryNormal, l, m, invP_Tensor);
  
  /*--- Compute projected velocities ---*/
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitaryNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
  }
  
  /*--- Calculate eigenvalues ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Lambda[iSpecies] = ProjVelocity;
  for (iDim = 0; iDim < nDim-1; iDim++)
    Lambda[nSpecies+iDim] = ProjVelocity;
  Lambda[nSpecies+nDim-1] = ProjVelocity + RoeSoundSpeed;
  Lambda[nSpecies+nDim]   = ProjVelocity - RoeSoundSpeed;
  Lambda[nSpecies+nDim+1] = ProjVelocity;
  
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
  
  /*--- Calculate inviscid projected Jacobians ---*/
  // Note: Scaling value is 0.5 because inviscid flux is based on 0.5*(Fc_i+Fc_j)
  GetInviscidProjJac(Density_i, Velocity_i, &Enthalpy_i, &Energy_ve_i,
                     dPdrhos_i, dPdrhoE_i, dPdrhoEve_i, Normal, 0.5, val_Jacobian_i);
  GetInviscidProjJac(Density_j, Velocity_j, &Enthalpy_j, &Energy_ve_j,
                     dPdrhos_j, dPdrhoE_j, dPdrhoEve_j, Normal, 0.5, val_Jacobian_j);
  
  /*--- Difference of conserved variables at iPoint and jPoint ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];
  
  /*--- Roe's Flux approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5 * (Proj_flux_tensor_i[iVar] + Proj_flux_tensor_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {

      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
      val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
      val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
    }
  }
  
/*  cout << endl << endl << "TNE2 Roe Jacobian_i: " << endl;
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      cout << val_Jacobian_i[iVar][jVar] << "\t";
    }
    cout << endl;
  }
  cout << endl << endl;
  
  cout << endl << endl << "TNE2 Roe P: " << endl;
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      cout << P_Tensor[iVar][jVar] << "\t";
    }
    cout << endl;
  }
  cout << endl << endl;
  
  cout << endl << endl << "TNE2 Roe invP: " << endl;
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      cout << invP_Tensor[iVar][jVar] << "\t";
    }
    cout << endl;
  }
  cout << endl << endl;
  
  cout << endl << endl << "TNE2 Lambda: " << endl;
  for (iVar = 0; iVar < nVar; iVar++) {
    cout << Lambda[iVar] << endl;
  }
  cout << endl << endl;
  cin.get();*/
  
  /* Visualization
  cout << "normal: " << endl;
  for (iDim = 0; iDim < nDim; iDim++)
    cout << UnitaryNormal[iDim] << endl;
  
  cout << "l: " << endl;
  for (iDim = 0; iDim < nDim; iDim++)
    cout << l[iDim] << endl;
  
  cout << "m: " << endl;
  for (iDim = 0; iDim < nDim; iDim++)
    cout << m[iDim] << endl;
  cout << endl << endl;
  
  cout << "lambda: " << endl;
  for(iVar =0 ;iVar<nVar; iVar++)
    cout << Lambda[iVar] << endl;
  
  cout << "dPdrhos_i[0]: " << dPdrhos_i[0] << endl;
  cout << "dPdrhoE: " << dPdrhoE << endl;
  cout << "dPdrhoEve: " << dPdrhoEve << endl;
  
  cout << "Projected Jacobian: " << endl;
  for (iVar =0; iVar < nVar; iVar++){
    for (jVar =0 ; jVar < nVar; jVar++) {
      cout << val_Jacobian_i[iVar][jVar] << "\t" ;
    }
    cout << endl;
  }
  cout << endl << endl;
  cout << "Projected Upwind: " << endl;
  for (iVar =0; iVar < nVar; iVar++){
    for (jVar =0 ; jVar < nVar; jVar++) {
      cout << val_Jacobian_j[iVar][jVar] << "\t" ;
    }
    cout << endl;
  }
  cin.get();
   */
  
}

void CUpwRoe_TNE2::CreateBasis(double *val_Normal) {
  unsigned short iDim;
  double modm, modl;
  
  /*--- Define l as a vector in the plane normal to the supplied vector ---*/
  l[0] = 0.0;
  l[1] = -val_Normal[2];
  l[2] = val_Normal[1];
  
  /*--- Check for the zero vector and re-assign if needed ---*/
  if (l[0] == 0.0 && l[1] == 0.0 && l[2] == 0.0) {
    l[0] = -val_Normal[2];
    l[1] = 0.0;
    l[2] = val_Normal[0];
  }
  
  /*--- Take vector product of n * l to make m ---*/
  m[0] = val_Normal[1]*l[2] - val_Normal[2]*l[1];
	m[1] = val_Normal[2]*l[0] - val_Normal[0]*l[2];
	m[2] = val_Normal[0]*l[1] - val_Normal[1]*l[0];
  
  /*--- Normalize ---*/
  modm =0 ; modl = 0;
  for (iDim =0 ; iDim < nDim; iDim++) {
    modm += m[iDim]*m[iDim];
    modl += l[iDim]*l[iDim];
  }
  modm = sqrt(modm);
  modl = sqrt(modl);
  for (iDim =0 ; iDim < nDim; iDim++) {
    l[iDim] = l[iDim]/modl;
    m[iDim] = m[iDim]/modm;
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
  dPdrhos          = new double[nSpecies];
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
  delete [] dPdrhos;
  
}

void CCentLax_TNE2::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i,
                                    double **val_Jacobian_j, CConfig *config) {
  
  unsigned short iDim, iSpecies, iVar, jVar, kVar, nHeavy, nEl;
  double *Ms;
  double DensityMix_i, DensityMix_j, conc;
  double Ru, rhoCvtr, rhoCvtr_i, rhoCvtr_j, rhoCvve, rhoCvve_i, rhoCvve_j, rho_el, dPdrhoE, dPdrhoEve;
  
  /*--- Face area (norm or the normal vector) ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Unit Normal ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Determine the number of heavy particle species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Read from CConfig ---*/
  Ms = config->GetMolar_Mass();
  
  /*--- Rename for convenience ---*/
  Ru = UNIVERSAL_GAS_CONSTANT;
  
  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, c]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Density_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    Density_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[VEL_INDEX+iDim];
    Velocity_j[iDim] = V_j[VEL_INDEX+iDim];
  }
  Pressure_i       = V_i[P_INDEX];
  Pressure_j       = V_j[P_INDEX];
  Enthalpy_i       = V_i[H_INDEX];
  Enthalpy_j       = V_j[H_INDEX];
  DensityMix_i     = V_i[RHO_INDEX];
  DensityMix_j     = V_j[RHO_INDEX];
  rhoCvtr_i        = V_i[RHOCVTR_INDEX];
  rhoCvtr_j        = V_j[RHOCVTR_INDEX];
  rhoCvve_i        = V_i[RHOCVVE_INDEX];
  rhoCvve_j        = V_j[RHOCVVE_INDEX];
  Temperature_i    = V_i[T_INDEX];
  Temperature_j    = V_j[T_INDEX];
  Temperature_ve_i = V_i[TVE_INDEX];
  Temperature_ve_j = V_j[TVE_INDEX];
  
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
    /*--- Calculate additional mean quantities ---*/
    rhoCvtr = 0.5*(rhoCvtr_i+rhoCvtr_j);
    rhoCvve = 0.5*(rhoCvve_i+rhoCvve_j);
    conc    = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      dPdrhos[iSpecies] = 0.5 * (dPdrhos_i[iSpecies] + dPdrhos_j[iSpecies]);
      conc             += MeanDensity[iSpecies]/Ms[iSpecies];
    }
    
    if (ionization) rho_el = MeanDensity[nSpecies-1];
    else            rho_el = 0.0;
    
    dPdrhoE   = Ru/rhoCvtr * conc;
    dPdrhoEve = -dPdrhoE + rho_el * Ru/Ms[nSpecies-1] * 1.0/rhoCvve;
    
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
			val_Jacobian_i[nSpecies+nDim][nSpecies+iDim] -= cte*UNIVERSAL_GAS_CONSTANT/rhoCvtr*conc * Velocity_i[iDim];
		val_Jacobian_i[nSpecies+nDim][nSpecies+nDim] += cte*(1+dPdrhoE);
    val_Jacobian_i[nSpecies+nDim][nSpecies+nDim+1] += cte*dPdrhoEve;
    val_Jacobian_i[nSpecies+nDim+1][nSpecies+nDim+1] += cte;
    
		/*--- Last row of Jacobian_j ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Jacobian_j[nSpecies+nDim][iSpecies] -= cte*dPdrhos[iSpecies];
		for (iDim = 0; iDim < nDim; iDim++)
			val_Jacobian_j[nSpecies+nDim][nSpecies+iDim] += cte*UNIVERSAL_GAS_CONSTANT/rhoCvtr*conc * Velocity_j[iDim];
		val_Jacobian_j[nSpecies+nDim][nSpecies+nDim] -= cte*(1+dPdrhoE);
    val_Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] -= cte*dPdrhoEve;
    val_Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] -= cte;
	}
}



CSource_TNE2::CSource_TNE2(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  X = new double[nSpecies];

}

CSource_TNE2::~CSource_TNE2(void) {
  delete [] X;
  
}

void CSource_TNE2::ComputeVibRelaxation(double *val_residual, double **val_Jacobian_i, CConfig *config) {
  /*--- Translational-rotational & vibrational energy exchange via inelastic collisions ---*/
  // Note: Electronic energy not implemented
	// Note: Landau-Teller formulation
  // Note: Millikan & White relaxation time (requires P in Atm.)
	// Note: Park limiting cross section
  
  bool ionization;
  unsigned short iDim, iSpecies, jSpecies, nEv, nHeavy, nEl;
  double rhos, evib, P, T, Tve, u, v, w, rhoCvtr, rhoCvve, Ru, conc, sqvel, N;
  double Qtv, estar, tau, tauMW, tauP;
  double tau_sr, mu, A_sr, B_sr, num, denom;
  double thoTve, exptv, evibs, eels;
  double thoT, expt, Cvvs, Cvvst, Cvtrs;
  double *dTdrhos, dTdrhou, dTdrhov, dTdrhow, dTdrhoE, dTdrhoEve;
  double *dTvedrhos, dTvedrhou, dTvedrhov, dTvedrhow, dTvedrhoE, dTvedrhoEve;
  double sigma, ws;
  double *Ms, *thetav, *Tref, *hf, *xi, ef;
  
  /*--- Initialize ---*/
  dTdrhos = NULL;
  dTvedrhos = NULL;
  
  /*--- Determine the number of heavy particle species ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Rename for convenience ---*/
  Ru      = UNIVERSAL_GAS_CONSTANT;
  P       = V_i[P_INDEX];
  T       = V_i[T_INDEX];
  Tve     = V_i[TVE_INDEX];
  u       = V_i[VEL_INDEX];
  v       = V_i[VEL_INDEX+1];
  w       = V_i[VEL_INDEX+2];
  rhoCvtr = V_i[RHOCVTR_INDEX];
  rhoCvve = V_i[RHOCVVE_INDEX];
  nEv     = nSpecies+nDim+1;
  
  /*--- Read from CConfig ---*/
  Ms     = config->GetMolar_Mass();
  thetav = config->GetCharVibTemp();
  Tref   = config->GetRefTemperature();
  hf     = config->GetEnthalpy_Formation();
  xi     = config->GetRotationModes();
  
  /*--- Calculate mole fractions ---*/
  N    = 0.0;
  conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    conc += V_i[RHOS_INDEX+iSpecies] / Ms[iSpecies];
    N    += V_i[RHOS_INDEX+iSpecies] / Ms[iSpecies] * AVOGAD_CONSTANT;
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    X[iSpecies] = (V_i[RHOS_INDEX+iSpecies] / Ms[iSpecies]) / conc;
  
  sqvel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel += V_i[VEL_INDEX+iDim]*V_i[VEL_INDEX+iDim];
  }
  
  /*--- Calculate partial derivatives of T & Tve ---*/
  if (implicit) {
    dTdrhos   = new double[nSpecies];
    dTvedrhos = new double[nSpecies];
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      Cvtrs = (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
      ef    = hf[iSpecies] - Ru/Ms[iSpecies] * Tref[iSpecies];
      evibs = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve) - 1.0);
      eels  = 0.0;
      dTdrhos[iSpecies]   = (-Cvtrs*(T - Tref[iSpecies]) - ef + 0.5*sqvel) / rhoCvtr;
      dTvedrhos[iSpecies] = -(evibs + eels) / rhoCvve;
    }
    for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
      ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];
      dTdrhos[nSpecies-1] = (-ef + 0.5*sqvel) / rhoCvtr;
      dTvedrhos[nSpecies-1] = -(3.0/2.0) * Ru/Ms[nSpecies-1] * Tve / rhoCvve;
    }
    dTdrhou     = -u / rhoCvtr;
    dTdrhov     = -v / rhoCvtr;
    dTdrhow     = -w / rhoCvtr;
    dTdrhoE     = 1.0 / rhoCvtr;
    dTdrhoEve   = -1.0 / rhoCvtr;
    dTvedrhou   = 0.0;
    dTvedrhov   = 0.0;
    dTvedrhow   = 0.0;
    dTvedrhoE   = 0.0;
    dTvedrhoEve = 1.0 / rhoCvve;
  }
  

  /*--- Loop over species to calculate source term --*/
  Qtv = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (thetav[iSpecies] != 0.0) {
      /*--- Rename ---*/
      rhos = V_i[RHOS_INDEX+iSpecies];
      thoT   = thetav[iSpecies]/T;
      expt   = exp(thetav[iSpecies]);
      thoTve = thetav[iSpecies]/Tve;
      exptv = exp(thetav[iSpecies]/Tve);
      
      /*--- Millikan & White relaxation time ---*/
      num   = 0.0;
      denom = 0.0;
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        mu     = Ms[iSpecies]*Ms[jSpecies] / (Ms[iSpecies] + Ms[jSpecies]);
        A_sr   = 1.16 * 1E-3 * sqrt(mu) * pow(thetav[iSpecies], 4.0/3.0);
        B_sr   = 0.015 * pow(mu, 0.25);
        tau_sr = 101325.0/P * exp(A_sr*(pow(T,-1.0/3.0) - B_sr) - 18.42);
        num   += X[iSpecies];
        denom += X[iSpecies] / tau_sr;
      }
      tauMW = num / denom;
      
      /*--- Park limiting cross section ---*/
      sigma = 1E-20 * (5E4/T)*(5E4/T);
      ws    = sqrt(8.0*Ru*T / (PI_NUMBER*Ms[iSpecies]));
      tauP  = 1.0 / (sigma * ws * N);
      
      /*--- Species relaxation time ---*/
      tau = tauMW + tauP;
      
      /*--- Vibrational energy terms ---*/
      estar = Ru/Ms[iSpecies] * thetav[iSpecies] / (expt - 1.0);
      evib  = Ru/Ms[iSpecies] * thetav[iSpecies] / (exptv - 1.0);
      
      /*--- Add species contribution to residual ---*/
      val_residual[nEv] += rhos * (estar - evib) / tau * Volume;
      
      if (implicit) {
        /*--- Calculate species specific heats ---*/
        Cvvst = Ru/Ms[iSpecies] * thoT*thoT * expt / ((expt-1.0)*(expt-1.0));
        Cvvs  = Ru/Ms[iSpecies] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));
        
        /*--- Density ---*/
        val_Jacobian_i[nEv][iSpecies] += (estar - evib)/tau * Volume;
        for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
          val_Jacobian_i[nEv][jSpecies] += U_i[iSpecies] * (Cvvst*dTdrhos[iSpecies] - Cvvs*dTvedrhos[iSpecies]) * Volume;
        
        /*--- Momentum ---*/
        val_Jacobian_i[nEv][nSpecies]      += U_i[iSpecies]/tau * (Cvvst*dTdrhou - Cvvs*dTvedrhou) * Volume;
        val_Jacobian_i[nEv][nSpecies+1]    += U_i[iSpecies]/tau * (Cvvst*dTdrhov - Cvvs*dTvedrhov) * Volume;
        val_Jacobian_i[nEv][nSpecies+2]    += U_i[iSpecies]/tau * (Cvvst*dTdrhow - Cvvs*dTvedrhow) * Volume;
        
        /*--- Energy ---*/
        val_Jacobian_i[nEv][nSpecies+nDim] += U_i[iSpecies]/tau * (Cvvst*dTdrhoE - Cvvs*dTvedrhoE) * Volume;
        
        /*--- Vibrational energy ---*/
        val_Jacobian_i[nEv][nEv]           += U_i[iSpecies]/tau * (Cvvst*dTdrhoEve - Cvvs*dTvedrhoEve) * Volume;
      }
    }
  }
  
  if (dTdrhos != NULL)   delete[] dTdrhos;
  if (dTvedrhos != NULL) delete [] dTvedrhos;
}