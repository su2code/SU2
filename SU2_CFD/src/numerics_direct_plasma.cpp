/*!
 * \file numerics_direct_plasma.cpp
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

CUpwRoe_Plasma::CUpwRoe_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	unsigned short iVar, iSpecies;
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	nVar_Species = nDim + 2;
	Energy_vib = 0.0;
	Energy_el = 0.0;
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Enthalpy_formation = new double [nSpecies];
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
		Enthalpy_formation[iSpecies] = config->GetEnthalpy_Formation(iSpecies);
  
	Diff_U = new double [nVar_Species];
	Velocity_i	= new double [nDim];
	Velocity_j	= new double [nDim];
	RoeVelocity = new double [nDim];
	delta_vel = new double [nDim];
	delta_wave			= new double [nVar_Species];
	Proj_flux_tensor_i	= new double [nVar_Species];
	Proj_flux_tensor_j	= new double [nVar_Species];
	Lambda				= new double [nVar_Species];
	Epsilon				= new double [nVar_Species];
  
	P_Tensor			= new double* [nVar_Species];
	invP_Tensor			= new double* [nVar_Species];
	ProjJac_i			= new double* [nVar_Species];
	ProjJac_j			= new double* [nVar_Species];
  
	for (iVar = 0; iVar < nVar_Species; iVar++) {
		P_Tensor[iVar]  = new double [nVar_Species];
		invP_Tensor[iVar] = new double [nVar_Species];
		ProjJac_i[iVar]  = new double [nVar_Species];
		ProjJac_j[iVar] = new double [nVar_Species];
	}
}

void CUpwRoe_Plasma::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	unsigned short iSpecies,  loc = 0;
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);                    /*! Area of the face*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/
  
	/*--- Point i, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = iSpecies*nVar_Species;
		Density_i= U_i[loc + 0];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[loc + iDim+1] / Density_i;
			sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
		}
    
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Gamma_Minus_One = Gamma - 1.0;
		Energy_i		= U_i[loc + nDim+1] / Density_i;
    //		Pressure_i	= Gamma_Minus_One*(Energy_i-0.5*sq_vel))*Density_i;
		Pressure_i  = Gamma_Minus_One*(Energy_i - 0.5*sq_vel - Enthalpy_formation[iSpecies] - Energy_vib - Energy_el)*Density_i;
    
		SoundSpeed_i	= sqrt(Pressure_i*Gamma/Density_i);
		Enthalpy_i		= (U_i[loc + nDim + 1] + Pressure_i) / Density_i;
    
		/*--- Point j, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
    
		Density_j		= U_j[loc + 0];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_j[iDim] = U_j[loc + iDim+1] / Density_j;
			sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
		}
    
		Energy_j		= U_j[loc + nDim+1] / Density_j;
    //		Pressure_j	= (Gamma_Minus_One*(Energy_j-0.5*sq_vel))*Density_j;
		Pressure_j  = Gamma_Minus_One*(Energy_j - 0.5*sq_vel - Enthalpy_formation[iSpecies] - Energy_vib - Energy_el)*Density_j;
    
		SoundSpeed_j	= sqrt(Pressure_j*Gamma/Density_j);
		Enthalpy_j		= (U_j[loc + nDim+1] + Pressure_j) / Density_j;
    
		/*--- Promediate Roe variables iPoint and jPoint ---*/
		R = sqrt(fabs(Density_j/Density_i));
		RoeDensity = R*Density_i;
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1.0);
			sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
		}
		RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
		RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
    
		/*--- Compute Proj_flux_tensor_i ---*/
		GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, Proj_flux_tensor_i);
		/*--- Compute Proj_flux_tensor_j ---*/
		GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, Proj_flux_tensor_j);
		/*--- Compute P and Lambda (do it with the Normal) ---*/
		GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, P_Tensor);
    
		ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity += RoeVelocity[iDim]*UnitaryNormal[iDim];
			ProjVelocity_i += Velocity_i[iDim] *UnitaryNormal[iDim];
			ProjVelocity_j += Velocity_j[iDim] *UnitaryNormal[iDim];
		}
    
		/*--- Flow eigenvalues and entropy correctors ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			Lambda[iDim] = ProjVelocity;
    
		Lambda[nVar_Species-2] = ProjVelocity + RoeSoundSpeed;
		Lambda[nVar_Species-1] = ProjVelocity - RoeSoundSpeed;
    
		/*--- Harten and Hyman (1983) entropy correction ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,ProjVelocity_j-Lambda[iDim]));
    
		Epsilon[nVar_Species-2] = 4.0*max(0.0, max(Lambda[nVar_Species-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar_Species-2]));
		Epsilon[nVar_Species-1] = 4.0*max(0.0, max(Lambda[nVar_Species-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar_Species-1]));
    
		for (iVar = 0; iVar < nVar_Species; iVar++) {
			if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
				Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
			else Lambda[iVar] = fabs(Lambda[iVar]);
		}
    
		for (iVar = 0; iVar < nVar_Species; iVar++)
			Lambda[iVar] = fabs(Lambda[iVar]);
    
		if (!implicit) {
      
			/*--- Compute wave amplitudes (characteristics) ---*/
			proj_delta_vel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
				proj_delta_vel += delta_vel[iDim]*Normal[iDim];
			}
			delta_p = Pressure_j - Pressure_i;
			delta_rho = Density_j - Density_i;
			proj_delta_vel = proj_delta_vel/Area;
      
			if (nDim == 2) {
				delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
				delta_wave[1] = UnitaryNormal[1]*delta_vel[0]-UnitaryNormal[0]*delta_vel[1];
				delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
				delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
			} else {
				delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
				delta_wave[1] = UnitaryNormal[0]*delta_vel[2]-UnitaryNormal[2]*delta_vel[0];
				delta_wave[2] = UnitaryNormal[1]*delta_vel[0]-UnitaryNormal[0]*delta_vel[1];
				delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
				delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
			}
      
			/*--- Roe's Flux approximation ---*/
			for (iVar = 0; iVar < nVar_Species; iVar++) {
				val_residual[loc + iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
				for (jVar = 0; jVar < nVar_Species; jVar++)
					val_residual[loc + iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
			}
		}
    
		else {
      
			/*--- Compute inverse P ---*/
			GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, invP_Tensor);
      
			for (iVar = 0; iVar < nVar_Species; iVar++) {
				for (jVar = 0; jVar < nVar_Species; jVar++) {
					ProjJac_i[iVar][jVar] = 0.0;
					ProjJac_j[iVar][jVar] = 0.0;
				}
			}
			/*--- Jacobians of the inviscid flux, scaled by
       0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
			GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, ProjJac_i);
			GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, ProjJac_j);
      
			/*--- Diference variables iPoint and jPoint ---*/
			for (iVar = 0; iVar < nVar_Species; iVar++)
				Diff_U[iVar] = U_j[loc + iVar]-U_i[loc + iVar];
      
      
		}
		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar_Species; iVar++) {
			val_residual[loc + iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
			for (jVar = 0; jVar < nVar_Species; jVar++) {
				Proj_ModJac_Tensor_ij = 0.0;
				/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
				for (kVar = 0; kVar < nVar_Species; kVar++)
					Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
				val_residual[loc + iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
				val_Jacobian_i[loc + iVar][loc + jVar] = ProjJac_i[iVar][jVar] + 0.5*Proj_ModJac_Tensor_ij*Area;
				val_Jacobian_j[loc + iVar][loc + jVar] = ProjJac_j[iVar][jVar] - 0.5*Proj_ModJac_Tensor_ij*Area;
			}
		}
    
	}
}


CUpwRoe_Plasma::~CUpwRoe_Plasma(void) {
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
	delete [] Enthalpy_formation;
  
	for (iVar = 0; iVar < nVar_Species; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
		delete [] ProjJac_i[iVar];
		delete [] ProjJac_j[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
	delete [] ProjJac_i;
	delete [] ProjJac_j;
}



CUpwRoe_Turkel_Plasma::CUpwRoe_Turkel_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) :
CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Beta_min = config->GetminTurkelBeta();
	Beta_max = config->GetmaxTurkelBeta();
	nVar_Species = nDim + 2;
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	Energy_vib = 0.0;
	Energy_el = 0.0;
  
	Diff_U = new double [nVar_Species];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	RoeVelocity = new double [nDim];
	Proj_flux_tensor_i = new double [nVar_Species];
	Proj_flux_tensor_j = new double [nVar_Species];
	Lambda = new double [nVar_Species];
	Epsilon = new double [nVar];
	absPeJac = new double* [nVar_Species];
	invRinvPe = new double* [nVar_Species];
	R_Tensor  = new double* [nVar_Species];
	Matrix    = new double* [nVar_Species];
	Art_Visc  = new double* [nVar_Species];
	Jac_ProjFlux_i  = new double* [nVar_Species];
	Jac_ProjFlux_j  = new double* [nVar_Species];
	for (iVar = 0; iVar < nVar_Species; iVar++) {
		absPeJac[iVar] = new double [nVar_Species];
		invRinvPe[iVar] = new double [nVar_Species];
		Matrix[iVar] = new double [nVar_Species];
		Art_Visc[iVar] = new double [nVar_Species];
		R_Tensor[iVar] = new double [nVar_Species];
		Jac_ProjFlux_i[iVar] = new double [nVar_Species];
		Jac_ProjFlux_j[iVar] = new double [nVar_Species];
	}
}

CUpwRoe_Turkel_Plasma::~CUpwRoe_Turkel_Plasma(void) {
	unsigned short iVar;
  
	delete [] Diff_U;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] RoeVelocity;
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
	delete [] Lambda;
	delete [] Epsilon;
	for (iVar = 0; iVar < nVar_Species; iVar++) {
		delete [] absPeJac[iVar];
		delete [] invRinvPe[iVar];
		delete [] Matrix[iVar];
		delete [] Art_Visc[iVar];
		delete [] R_Tensor[iVar];
		delete [] Jac_ProjFlux_i[iVar];
		delete [] Jac_ProjFlux_j[iVar];
	}
	delete [] Matrix;
	delete [] Art_Visc;
	delete [] absPeJac;
	delete [] invRinvPe;
	delete [] R_Tensor;
	delete [] Jac_ProjFlux_i;
	delete [] Jac_ProjFlux_j;
  
  
}

void CUpwRoe_Turkel_Plasma::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	double r_hat, s_hat, t_hat, rhoB2a2, sqr_one_m_Betasqr_Lam1;
	double Beta2, one_m_Betasqr, one_p_Betasqr, sqr_two_Beta_c_Area;
	unsigned short iSpecies, loc = 0;
	bool precondition = config->Low_Mach_Preconditioning();
  
	/*--- Face area (norm or the normal vector) ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Unit Normal ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = nVar_Species*iSpecies;
		/*--- Conserved variables at point i,
     Need to recompute SoundSpeed / Pressure / Enthalpy in
     case of 2nd order reconstruction ---*/
		Density_i = U_i[loc + 0];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[loc + iDim+1] / Density_i;
			sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
		}
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Gamma_Minus_One = Gamma - 1.0;
		Energy_i		= U_i[loc + nDim+1] / Density_i;
    //		Pressure_i	= Gamma_Minus_One*(Energy_i-0.5*sq_vel))*Density_i;
		Pressure_i  = Gamma_Minus_One*(Energy_i - 0.5*sq_vel - Enthalpy_formation[iSpecies] - Energy_vib - Energy_el)*Density_i;
    
		SoundSpeed_i	= sqrt(Pressure_i*Gamma/Density_i);
		Enthalpy_i		= (U_i[loc + nDim + 1] + Pressure_i) / Density_i;
    
		/*--- Conserved variables at point j,
     Need to recompute SoundSpeed / Pressure / Enthalpy in
     case of 2nd order reconstruction ---*/
		Density_j		= U_j[loc + 0];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_j[iDim] = U_j[loc + iDim+1] / Density_j;
			sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
		}
    
		Energy_j		= U_j[loc + nDim+1] / Density_j;
    //		Pressure_j	= (Gamma_Minus_One*(Energy_j-0.5*sq_vel))*Density_j;
		Pressure_j  = Gamma_Minus_One*(Energy_j - 0.5*sq_vel - Enthalpy_formation[iSpecies] - Energy_vib - Energy_el)*Density_j;
    
		SoundSpeed_j	= sqrt(Pressure_j*Gamma/Density_j);
		Enthalpy_j		= (U_j[loc + nDim+1] + Pressure_j) / Density_j;
    
		/*--- Roe-averaged variables at interface between i & j ---*/
		R = sqrt(Density_j/Density_i);
		RoeDensity = R*Density_i;
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
			sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
		}
		RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
		RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
		RoePressure = RoeDensity/Gamma*RoeSoundSpeed*RoeSoundSpeed;
    
		/*--- Compute Proj_flux_tensor_i ---*/
		GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, Proj_flux_tensor_i);
    
		/*--- Compute Proj_flux_tensor_j ---*/
		GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, Proj_flux_tensor_j);
    
		ProjVelocity = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjVelocity   += RoeVelocity[iDim]*Normal[iDim];
    
		/*--- First few flow eigenvalues of A.Normal with the normal---*/
		for (iDim = 0; iDim < nDim; iDim++)
			Lambda[iDim] = ProjVelocity;
    
		/*-- Default value of preconditioning parameter Beta ---*/
		Beta = 1.0;
		if (precondition && (iSpecies == nSpecies-1)) {
			double local_Mach = sqrt(sq_vel)/RoeSoundSpeed;
			Beta 		    = max(Beta_min,min(local_Mach,Beta_max));
		}
    
		Beta2 				   = Beta*Beta;
		one_m_Betasqr 		   = 1.0 - Beta2;  // 1-Beta*Beta
		one_p_Betasqr 		   = 1.0 + Beta2;  // 1+ Beta*Beta
		sqr_one_m_Betasqr_Lam1 = pow((one_m_Betasqr*Lambda[0]),2); // [(1-Beta^2)*Lambda[0]]^2
		sqr_two_Beta_c_Area    = pow(2.0*Beta*RoeSoundSpeed*Area,2); // [2*Beta*c*Area]^2
    
		/*--- The rest of the flow eigenvalues of preconditioned matrix---*/
		Lambda[nVar_Species-2] = 0.5 * ( one_p_Betasqr*Lambda[0] + sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
		Lambda[nVar_Species-1] = 0.5 * ( one_p_Betasqr*Lambda[0] - sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
    
		s_hat = 1.0/Area * (Lambda[nVar_Species-1] - Lambda[0]*Beta2);
		r_hat = 1.0/Area * (Lambda[nVar_Species-2] - Lambda[0]*Beta2);
		t_hat = 0.5/Area * (Lambda[nVar_Species-1] - Lambda[nVar_Species-2]);
		rhoB2a2 = RoeDensity*Beta2*RoeSoundSpeed*RoeSoundSpeed;
    
		/*--- Diference variables iPoint and jPoint and absolute value of the eigen values---*/
		for (iVar = 0; iVar < nVar_Species; iVar++) {
			Lambda[iVar] = fabs(Lambda[iVar]);
			Diff_U[iVar] = U_j[loc + iVar] - U_i[loc + iVar];
		}
    
		/*--- Compute the absolute Preconditioned Jacobian in entropic Variables (do it with the Unitary Normal) ---*/
		GetPrecondJacobian(Beta2, r_hat, s_hat, t_hat, rhoB2a2, Lambda, UnitaryNormal, absPeJac);
    
		/*--- Compute the matrix from entropic variables to conserved variables ---*/
		GetinvRinvPe(Beta2, RoeEnthalpy, RoeSoundSpeed, RoeDensity, RoeVelocity, invRinvPe);
    
		/*--- Compute the matrix from entropic variables to conserved variables ---*/
		GetRMatrix(RoePressure, RoeSoundSpeed, RoeDensity, RoeVelocity, R_Tensor);
    
    
		for (iVar = 0; iVar < nVar_Species; iVar ++){
			for (jVar = 0; jVar < nVar_Species; jVar ++) {
				Matrix[iVar][jVar] = 0.0;
				for (kVar = 0; kVar < nVar_Species; kVar++)
					Matrix[iVar][jVar]  += absPeJac[iVar][kVar]*R_Tensor[kVar][jVar];
			}
		}
    
		for (iVar = 0; iVar < nVar_Species; iVar ++){
			for (jVar = 0; jVar < nVar_Species; jVar ++) {
				Art_Visc[iVar][jVar] = 0.0;
				for (kVar = 0; kVar < nVar_Species; kVar++)
					Art_Visc[iVar][jVar]  += invRinvPe[iVar][kVar]*Matrix[kVar][jVar];
			}
		}
    
		if (implicit) {
			/*--- Jacobians of the inviscid flux, scaled by
       0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
			GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Jac_ProjFlux_i);
			GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Jac_ProjFlux_j);
		}
    
		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar_Species; iVar++) {
			val_residual[loc + iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
			for (jVar = 0; jVar < nVar_Species; jVar++) {
				val_residual[loc + iVar] -= 0.5*Art_Visc[iVar][jVar]*Diff_U[jVar];
				if (implicit) {
					val_Jacobian_i[loc + iVar][loc + jVar] = Jac_ProjFlux_i[iVar][jVar] + 0.5*Art_Visc[iVar][jVar];
					val_Jacobian_j[loc + iVar][loc + jVar] = Jac_ProjFlux_j[iVar][jVar] - 0.5*Art_Visc[iVar][jVar];
				}
			}
		}
	}
  
}

CUpwRoe_PlasmaDiatomic::CUpwRoe_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	unsigned short iVar, iSpecies;
  
	nMonatomics = val_nMonatomics;
	nDiatomics  = val_nDiatomics;
  
	GammaMonatomic = config->GetGammaMonatomic();
	GammaDiatomic = config->GetGammaDiatomic();
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Diff_U = new double [nVar];
  
	Density_i		= new double[nSpecies];
	Energy_i		= new double[nSpecies];
	Energy_vib_i = new double [nSpecies];
	Energy_el_i = new double[nSpecies];
	SoundSpeed_i	= new double[nSpecies];
	Pressure_i		= new double[nSpecies];
	Enthalpy_i		= new double[nSpecies];
  
	Density_j		= new double[nSpecies];
	Energy_j		= new double[nSpecies];
	Energy_vib_j = new double[nSpecies];
	Energy_el_j = new double[nSpecies];
	SoundSpeed_j	= new double[nSpecies];
	Pressure_j		= new double[nSpecies];
	Enthalpy_j		= new double[nSpecies];
  
	RoeDensity		= new double[nSpecies];
	RoeEnthalpy		= new double[nSpecies];
	RoeSoundSpeed	= new double[nSpecies];
	RoeEnergy_vib = new double[nSpecies];
  
	ProjVelocity	= new double[nSpecies];
	ProjVelocity_i	= new double[nSpecies];
	ProjVelocity_j	= new double[nSpecies];
  
	Velocity_i		= new double* [nSpecies];
	Velocity_j		= new double* [nSpecies];
	RoeVelocity		= new double* [nSpecies];
  
	delta_vel		= new double* [nSpecies];
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Velocity_i[iSpecies]	= new double [nDim];
		Velocity_j[iSpecies]	= new double [nDim];
		RoeVelocity[iSpecies]	= new double [nDim];
		delta_vel[iSpecies]		= new double [nDim];
    
	}
  
	delta_wave			= new double [nVar];
	Proj_flux_tensor_i	= new double [nVar];
	Proj_flux_tensor_j	= new double [nVar];
	Lambda				= new double [nVar];
	Epsilon				= new double [nVar];
  
	P_Tensor			= new double* [nVar];
	invP_Tensor			= new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]  = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
	}
}

CUpwRoe_PlasmaDiatomic::~CUpwRoe_PlasmaDiatomic(void) {
	unsigned short iVar, iSpecies;
  
  
	delete [] Diff_U;
  
	delete [] Density_i;
	delete [] Energy_i;
	delete [] Energy_vib_i;
	delete [] Energy_el_i;
	delete [] SoundSpeed_i;
	delete [] Pressure_i;
	delete [] Enthalpy_i;
  
	delete [] Density_j;
	delete [] Energy_j;
	delete [] Energy_vib_j;
	delete [] Energy_el_j;
	delete [] SoundSpeed_j;
	delete [] Pressure_j;
	delete [] Enthalpy_j;
  
	delete [] RoeDensity;
	delete [] RoeEnthalpy;
	delete [] RoeSoundSpeed;
	delete [] RoeEnergy_vib;
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		delete [] Velocity_i[iSpecies];
		delete [] Velocity_j[iSpecies];
		delete [] RoeVelocity[iSpecies];
		delete [] delta_vel[iSpecies];
	}
  
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] RoeVelocity;
	delete [] delta_vel[iSpecies];
  
	delete [] ProjVelocity;
	delete [] ProjVelocity_i;
	delete [] ProjVelocity_j;
  
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

void CUpwRoe_PlasmaDiatomic::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	unsigned short iSpecies,  loc = 0;
	Area = 0;
  
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
  
	Area = sqrt(Area);                    /*! Area of the face*/
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/
  
	/*--- Point i, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		Density_i[iSpecies]	= U_i[loc + 0];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i[iSpecies];
			Vel2 += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
		}
		Energy_i[iSpecies]		= U_i[loc+nDim+1] / Density_i[iSpecies];
		Energy_el_i[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_i[iSpecies]  = U_i[loc+nDim+2] / Density_i[iSpecies];
			SoundSpeed_i[iSpecies] 	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_i[iSpecies] - 0.5*Vel2 - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_i[iSpecies] = (GammaDiatomic-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies]);
		}
		else {
			SoundSpeed_i[iSpecies] 	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_i[iSpecies] - 0.5*Vel2 - Energy_el_i[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_i[iSpecies] = (GammaMonatomic-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el_i[iSpecies]);
		}
		Enthalpy_i[iSpecies] = Energy_i[iSpecies] + Pressure_i[iSpecies] / Density_i[iSpecies];
    
		/*--- Point j, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
		Density_j[iSpecies]		= U_j[loc + 0];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_j[iSpecies][iDim] = U_j[loc+iDim+1] / Density_j[iSpecies];
			Vel2 += Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
		}
		Energy_j[iSpecies]		= U_j[loc+nDim+1] / Density_j[iSpecies];
		if (iSpecies < nDiatomics) {
			Energy_vib_j[iSpecies]  = U_j[loc+nDim+2] / Density_j[iSpecies];
			SoundSpeed_j[iSpecies] 	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_j[iSpecies] - 0.5*Vel2 - Energy_vib_j[iSpecies] - Energy_el_j[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_j[iSpecies] = (GammaDiatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_vib_j[iSpecies] - Energy_el_j[iSpecies]);
		}
		else {
			SoundSpeed_j[iSpecies] 	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_j[iSpecies] - 0.5*Vel2 - Energy_el_j[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_j[iSpecies] = (GammaMonatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el_j[iSpecies]);
		}
		Enthalpy_j[iSpecies] = Energy_j[iSpecies] + Pressure_j[iSpecies] / Density_j[iSpecies];
    
		/*--- Average Roe variables iPoint and jPoint ---*/
		R = sqrt(fabs(Density_j[iSpecies]/Density_i[iSpecies]));
		RoeDensity[iSpecies] = R*Density_i[iSpecies];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			RoeVelocity[iSpecies][iDim] = (R*Velocity_j[iSpecies][iDim]+Velocity_i[iSpecies][iDim])/(R+1.0);
			Vel2 += RoeVelocity[iSpecies][iDim]*RoeVelocity[iSpecies][iDim];
		}
		RoeEnthalpy[iSpecies] = (R*Enthalpy_j[iSpecies]+Enthalpy_i[iSpecies])/(R+1);
		RoeEnergy_vib[iSpecies] = (R*Energy_vib_j[iSpecies] + Energy_vib_i[iSpecies])/(R+1);
		if (iSpecies < nDiatomics)
			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaDiatomic-1.0)*(RoeEnthalpy[iSpecies]-0.5*Vel2 - RoeEnergy_vib[iSpecies] - config->GetEnthalpy_Formation(iSpecies))));
		else
			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaMonatomic-1.0)*(RoeEnthalpy[iSpecies] - 0.5*Vel2 - config->GetEnthalpy_Formation(iSpecies))));
	}
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		ProjVelocity[iSpecies]		= 0.0;
		ProjVelocity_i[iSpecies]		= 0.0;
		ProjVelocity_j[iSpecies]		= 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity[iSpecies]   += RoeVelocity[iSpecies][iDim]*UnitaryNormal[iDim];
			ProjVelocity_i[iSpecies] += Velocity_i[iSpecies][iDim] *UnitaryNormal[iDim];
			ProjVelocity_j[iSpecies] += Velocity_j[iSpecies][iDim] *UnitaryNormal[iDim];
		}
	}
	/*--- Flow eigenvalues and Entropy correctors ---*/
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		for (iDim = 0; iDim < nDim; iDim++) {
      
			Lambda[loc + iDim] = ProjVelocity[iSpecies];
		}
		Lambda[loc + nDim]  = ProjVelocity[iSpecies] + RoeSoundSpeed[iSpecies];
		Lambda[loc + nDim+1] = ProjVelocity[iSpecies] - RoeSoundSpeed[iSpecies];
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = ProjVelocity[iSpecies];
	}
  
	for (iVar = 0; iVar < nVar; iVar++) {
		Lambda[iVar] = fabs(Lambda[iVar]);
	}
  
	/*--- Compute Proj_flux_tensor_i ---*/
	GetInviscidProjFlux_(Density_i, Velocity_i, Pressure_i, Enthalpy_i, Energy_vib_i, Normal, Proj_flux_tensor_i);
  
	/*--- Compute Proj_flux_tensor_j ---*/
	GetInviscidProjFlux_(Density_j, Velocity_j, Pressure_j, Enthalpy_j, Energy_vib_j, Normal, Proj_flux_tensor_j);
  
	/*--- Compute P and Lambda (do it with the Normal) ---*/
	GetPMatrix_(RoeDensity, RoeVelocity, RoeEnthalpy, RoeSoundSpeed, RoeEnergy_vib, Energy_el_i, config, UnitaryNormal, P_Tensor);
  
	/*--- Compute inverse P ---*/
	GetPMatrix_inv_(RoeDensity, RoeVelocity, RoeSoundSpeed, RoeEnergy_vib, Energy_el_i, config, UnitaryNormal, invP_Tensor);
  
	/*--- Jacobians of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidProjJac_(Velocity_i, Energy_i, Energy_vib_i, Enthalpy_i, Normal, 0.5, val_Jacobian_i, config);
	GetInviscidProjJac_(Velocity_j, Energy_j, Energy_vib_j, Enthalpy_j, Normal, 0.5, val_Jacobian_j, config);
  
	/*--- Difference conserved variables between iPoint and jPoint ---*/
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

CUpwSW_PlasmaDiatomic::CUpwSW_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	unsigned short iVar, iSpecies;
  
	nMonatomics = val_nMonatomics;
	nDiatomics  = val_nDiatomics;
  
	GammaMonatomic = config->GetGammaMonatomic();
	GammaDiatomic = config->GetGammaDiatomic();
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Diff_U = new double [nVar];
  
	Density_i		= new double[nSpecies];
	Energy_i		= new double[nSpecies];
	Energy_vib_i = new double [nSpecies];
	Energy_el_i = new double[nSpecies];
	SoundSpeed_i	= new double[nSpecies];
	Pressure_i		= new double[nSpecies];
	Enthalpy_i		= new double[nSpecies];
  
	Density_j		= new double[nSpecies];
	Energy_j		= new double[nSpecies];
	Energy_vib_j = new double[nSpecies];
	Energy_el_j = new double[nSpecies];
	SoundSpeed_j	= new double[nSpecies];
	Pressure_j		= new double[nSpecies];
	Enthalpy_j		= new double[nSpecies];
  
	Density_ij		= new double[nSpecies];
	Enthalpy_ij		= new double[nSpecies];
	SoundSpeed_ij	= new double[nSpecies];
	Energy_vib_ij = new double[nSpecies];
	Energy_el_ij  = new double[nSpecies];
  
	ProjVelocity_ij	= new double[nSpecies];
	ProjVelocity_i	= new double[nSpecies];
	ProjVelocity_j	= new double[nSpecies];
  
	Velocity_i		= new double* [nSpecies];
	Velocity_j		= new double* [nSpecies];
	Velocity_ij		= new double* [nSpecies];
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Velocity_i[iSpecies]	= new double [nDim];
		Velocity_j[iSpecies]	= new double [nDim];
		Velocity_ij[iSpecies]	= new double [nDim];
	}
  
	Proj_flux_tensor_i	= new double [nVar];
	Proj_flux_tensor_j	= new double [nVar];
	Lambda				= new double [nVar];
  
	P_Tensor			= new double* [nVar];
	invP_Tensor			= new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]  = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
	}
}

CUpwSW_PlasmaDiatomic::~CUpwSW_PlasmaDiatomic(void) {
	unsigned short iVar, iSpecies;
  
  
	delete [] Diff_U;
  
	delete [] Density_i;
	delete [] Energy_i;
	delete [] Energy_vib_i;
	delete [] Energy_el_i;
	delete [] SoundSpeed_i;
	delete [] Pressure_i;
	delete [] Enthalpy_i;
  
	delete [] Density_j;
	delete [] Energy_j;
	delete [] Energy_vib_j;
	delete [] Energy_el_j;
	delete [] SoundSpeed_j;
	delete [] Pressure_j;
	delete [] Enthalpy_j;
  
	delete [] Density_ij;
	delete [] Enthalpy_ij;
	delete [] SoundSpeed_ij;
	delete [] Energy_vib_ij;
	delete [] Energy_el_ij;
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		delete [] Velocity_i[iSpecies];
		delete [] Velocity_j[iSpecies];
		delete [] Velocity_ij[iSpecies];
	}
  
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Velocity_ij;
  
	delete [] ProjVelocity_ij;
	delete [] ProjVelocity_i;
	delete [] ProjVelocity_j;
  
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
	delete [] Lambda;
  
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
}

void CUpwSW_PlasmaDiatomic::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	unsigned short iSpecies,  loc = 0;
	double epsilon;
	epsilon = 1E-4;
	Area = 0;
  
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
  
	Area = sqrt(Area);                    /*! Area of the face*/
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/
  
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_flux_tensor_i[iVar] = 0.0;
		Proj_flux_tensor_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[iVar][jVar] = 0.0;
			val_Jacobian_j[iVar][jVar] = 0.0;
		}
	}
  
	/*--- Point i, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		Density_i[iSpecies]	= U_i[loc + 0];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i[iSpecies];
			Vel2 += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
		}
		Energy_i[iSpecies]		= U_i[loc+nDim+1] / Density_i[iSpecies];
		Energy_vib_i[iSpecies] = 0.0;
		Energy_el_i[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_i[iSpecies]  = U_i[loc+nDim+2] / Density_i[iSpecies];
			SoundSpeed_i[iSpecies] 	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_i[iSpecies] - 0.5*Vel2 - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_i[iSpecies] = (GammaDiatomic-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies]);
		}
		else {
			SoundSpeed_i[iSpecies] 	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_i[iSpecies] - 0.5*Vel2 - Energy_el_i[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_i[iSpecies] = (GammaMonatomic-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el_i[iSpecies]);
		}
		//		Enthalpy_i[iSpecies] = (U_i[loc + nDim+1] + Pressure_i[iSpecies]) / Density_i[iSpecies];
		Enthalpy_i[iSpecies] = Energy_i[iSpecies] + Pressure_i[iSpecies] / Density_i[iSpecies];
    
    
		/*--- Point j, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
		Density_j[iSpecies]		= U_j[loc + 0];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_j[iSpecies][iDim] = U_j[loc+iDim+1] / Density_j[iSpecies];
			Vel2 += Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
		}
		Energy_j[iSpecies]		= U_j[loc+nDim+1] / Density_j[iSpecies];
		Energy_vib_j[iSpecies] = 0.0;
		Energy_el_j[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_j[iSpecies]  = U_j[loc+nDim+2] / Density_j[iSpecies];
			SoundSpeed_j[iSpecies] 	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_j[iSpecies] - 0.5*Vel2 - Energy_vib_j[iSpecies] - Energy_el_j[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_j[iSpecies] = (GammaDiatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_vib_j[iSpecies] - Energy_el_j[iSpecies]);
		}
		else {
			SoundSpeed_j[iSpecies] 	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_j[iSpecies] - 0.5*Vel2 - Energy_el_j[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_j[iSpecies] = (GammaMonatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el_j[iSpecies]);
		}
		Enthalpy_j[iSpecies] = Energy_j[iSpecies] + Pressure_j[iSpecies] / Density_j[iSpecies];
	}
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		ProjVelocity_i[iSpecies]  = 0.0;
		ProjVelocity_j[iSpecies]	= 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity_i[iSpecies]  += Velocity_i[iSpecies][iDim] *UnitaryNormal[iDim];
			ProjVelocity_j[iSpecies]  += Velocity_j[iSpecies][iDim] *UnitaryNormal[iDim];
		}
	}
  
	/*--- Flow eigenvalues at i (Lambda+) --- */
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda[loc + iDim] = 0.5*(ProjVelocity_i[iSpecies] + fabs(ProjVelocity_i[iSpecies]));
		}
		Lambda[loc + nDim]   = 0.5*(ProjVelocity_i[iSpecies] + SoundSpeed_i[iSpecies]
                                + fabs(ProjVelocity_i[iSpecies] + SoundSpeed_i[iSpecies]));
		Lambda[loc + nDim+1] = 0.5*(ProjVelocity_i[iSpecies] - SoundSpeed_i[iSpecies]
                                + fabs(ProjVelocity_i[iSpecies] - SoundSpeed_i[iSpecies]));
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = 0.5*(ProjVelocity_i[iSpecies] + fabs(ProjVelocity_i[iSpecies]));
	}
  
	/*--- Correct for the sonic glitch ---*/
	/*  for (iVar = 0; iVar < nDim+2; iVar++)
   Lambda[loc+iVar] = 0.5* (Lambda[loc+iVar] + sqrt(Lambda[loc+iVar]*Lambda[loc+iVar] + epsilon*epsilon));
   if (iSpecies < nDiatomics)
   Lambda[loc+nDim+2] = 0.5* (Lambda[loc+nDim+2] + sqrt(Lambda[loc+nDim+2]*Lambda[loc+nDim+2] + epsilon*epsilon));*/
  
	/*--- Compute P & invP at i ---*/
	GetPMatrix_(Density_i, Velocity_i, Enthalpy_i, SoundSpeed_i, Energy_vib_i, Energy_el_i, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_i, Velocity_i, SoundSpeed_i, Energy_vib_i, Energy_el_i, config, UnitaryNormal, invP_Tensor);
  
	/*--- Projected flux (f+) at i ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_i = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_i += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_flux_tensor_i[iVar] += Proj_ModJac_Tensor_i*U_i[jVar]*Area;
			val_Jacobian_i[iVar][jVar] += Proj_ModJac_Tensor_i*Area;
		}
	}
  
	/*--- Flow eigenvalues at j (Lambda-) --- */
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda[loc + iDim] = 0.5*(ProjVelocity_j[iSpecies] - fabs(ProjVelocity_j[iSpecies]));
		}
		Lambda[loc + nDim]   = 0.5*(ProjVelocity_j[iSpecies] + SoundSpeed_j[iSpecies]
                                - fabs(ProjVelocity_j[iSpecies] + SoundSpeed_j[iSpecies]));
		Lambda[loc + nDim+1] = 0.5*(ProjVelocity_j[iSpecies] - SoundSpeed_j[iSpecies]
                                - fabs(ProjVelocity_j[iSpecies] - SoundSpeed_j[iSpecies]));
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = 0.5*(ProjVelocity_i[iSpecies] - fabs(ProjVelocity_i[iSpecies]));
	}
  
	/*--- Correct for the sonic glitch ---*/
	/*  for (iVar = 0; iVar < nDim+2; iVar++)
   Lambda[loc+iVar] = 0.5* (Lambda[loc+iVar] - sqrt(Lambda[loc+iVar]*Lambda[loc+iVar] + epsilon*epsilon));
   if (iSpecies < nDiatomics)
   Lambda[loc+nDim+2] = 0.5* (Lambda[loc+nDim+2] - sqrt(Lambda[loc+nDim+2]*Lambda[loc+nDim+2] + epsilon*epsilon));*/
  
	/*--- Compute P & invP at j ---*/
	GetPMatrix_(Density_j, Velocity_j, Enthalpy_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_j, Velocity_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, invP_Tensor);
  
	/*--- Projected flux (f-) ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_j = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_j += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_flux_tensor_j[iVar] += Proj_ModJac_Tensor_j*U_j[jVar]*Area;
			val_Jacobian_j[iVar][jVar] += Proj_ModJac_Tensor_j*Area;
		}
	}
  
	/*--- Flux splitting ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar];
	}
}

CUpwMSW_PlasmaDiatomic::CUpwMSW_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	unsigned short iVar, iSpecies;
  
	nMonatomics = val_nMonatomics;
	nDiatomics  = val_nDiatomics;
  
	GammaMonatomic = config->GetGammaMonatomic();
	GammaDiatomic = config->GetGammaDiatomic();
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Diff_U = new double [nVar];
  
	Density_i		 = new double[nSpecies];
	Energy_i		 = new double[nSpecies];
	Energy_vib_i = new double [nSpecies];
	Energy_el_i  = new double[nSpecies];
	SoundSpeed_i = new double[nSpecies];
	Pressure_i   = new double[nSpecies];
	Enthalpy_i	 = new double[nSpecies];
  
	Density_j    = new double[nSpecies];
	Energy_j     = new double[nSpecies];
	Energy_vib_j = new double[nSpecies];
	Energy_el_j  = new double[nSpecies];
	SoundSpeed_j = new double[nSpecies];
	Pressure_j	 = new double[nSpecies];
	Enthalpy_j	 = new double[nSpecies];
  
  Densityst_i    = new double[nSpecies];
  Velocityst_i   = new double*[nSpecies];
  Soundspeedst_i = new double[nSpecies];
  Enthalpyst_i   = new double[nSpecies];
  Energy_vibst_i = new double[nSpecies];
  Energy_elst_i  = new double[nSpecies];
  
  Densityst_j    = new double[nSpecies];
  Velocityst_j   = new double*[nSpecies];
  Soundspeedst_j = new double[nSpecies];
  Enthalpyst_j   = new double[nSpecies];
  Energy_vibst_j = new double[nSpecies];
  Energy_elst_j  = new double[nSpecies];
  
	ProjVelocity_i	= new double[nSpecies];
	ProjVelocity_j	= new double[nSpecies];
  ProjVelocityst_i = new double[nSpecies];
  ProjVelocityst_j = new double[nSpecies];
  
	Velocity_i		= new double* [nSpecies];
	Velocity_j		= new double* [nSpecies];
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Velocity_i[iSpecies]	= new double [nDim];
		Velocity_j[iSpecies]	= new double [nDim];
    
    Velocityst_i[iSpecies] = new double[nDim];
    Velocityst_j[iSpecies] = new double[nDim];
	}
  
	Proj_flux_tensor_i	= new double [nVar];
	Proj_flux_tensor_j	= new double [nVar];
	Lambda_i				    = new double [nVar];
  Lambda_j				    = new double [nVar];
  
	P_Tensor			= new double* [nVar];
	invP_Tensor			= new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]  = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
	}
}

CUpwMSW_PlasmaDiatomic::~CUpwMSW_PlasmaDiatomic(void) {
	unsigned short iVar, iSpecies;
  
  
	delete [] Diff_U;
  
	delete [] Density_i;      delete [] Density_j;
	delete [] Energy_i;       delete [] Energy_j;
	delete [] Energy_vib_i;   delete [] Energy_vib_j;
	delete [] Energy_el_i;    delete [] Energy_el_j;
	delete [] SoundSpeed_i;   delete [] SoundSpeed_j;
	delete [] Pressure_i;     delete [] Pressure_j;
	delete [] Enthalpy_i;     delete [] Enthalpy_j;
  
  delete [] Densityst_i;    delete [] Densityst_j;
  delete [] Soundspeedst_i; delete [] Soundspeedst_j;
  delete [] Enthalpyst_i;   delete [] Enthalpyst_j;
  delete [] Energy_vibst_i; delete [] Energy_vibst_j;
  delete [] Energy_elst_i;  delete [] Energy_elst_j;
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		delete [] Velocity_i[iSpecies];
		delete [] Velocity_j[iSpecies];
    delete [] Velocityst_i[iSpecies];
    delete [] Velocityst_j[iSpecies];
	}
  
	delete [] Velocity_i;
	delete [] Velocity_j;
  delete [] Velocityst_i;
  delete [] Velocityst_j;
  
	delete [] ProjVelocity_i;
	delete [] ProjVelocity_j;
  delete [] ProjVelocityst_i;
  delete [] ProjVelocityst_j;
  
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
	delete [] Lambda_i;
  delete [] Lambda_j;
  
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
}

void CUpwMSW_PlasmaDiatomic::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	unsigned short iSpecies,  loc = 0;
	double epsilon, alpha, w, dp, onemw;
	epsilon = 1E-4;
  alpha = 6.0;
	Area = 0;
  
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
  
	Area = sqrt(Area);                    /*! Area of the face*/
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/
  
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_flux_tensor_i[iVar] = 0.0;
		Proj_flux_tensor_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[iVar][jVar] = 0.0;
			val_Jacobian_j[iVar][jVar] = 0.0;
		}
	}
  
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
    /*--- Point i: recompute sound speed, pressure and enthalpy in case of 2nd order reconstruction ---*/
		Density_i[iSpecies]	= U_i[loc + 0];
		Vel2 = 0;
    ProjVelocity_i[iSpecies]  = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i[iSpecies];
      ProjVelocity_i[iSpecies]  += Velocity_i[iSpecies][iDim] *UnitaryNormal[iDim];
			Vel2                      += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
		}
		Energy_i[iSpecies]		= U_i[loc+nDim+1] / Density_i[iSpecies];
		Energy_vib_i[iSpecies] = 0.0;
		Energy_el_i[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_i[iSpecies]  = U_i[loc+nDim+2] / Density_i[iSpecies];
			SoundSpeed_i[iSpecies] 	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_i[iSpecies] - 0.5*Vel2 - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_i[iSpecies] = (GammaDiatomic-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies]);
		}
		else {
			SoundSpeed_i[iSpecies] 	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_i[iSpecies] - 0.5*Vel2 - Energy_el_i[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_i[iSpecies] = (GammaMonatomic-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el_i[iSpecies]);
		}
		Enthalpy_i[iSpecies] = Energy_i[iSpecies] + Pressure_i[iSpecies] / Density_i[iSpecies];
    
    
		/*--- Point j: recompute sound speed, pressure and enthalpy in case of 2nd order reconstruction ---*/
		Density_j[iSpecies]		= U_j[loc + 0];
		Vel2 = 0;
		ProjVelocity_j[iSpecies]	= 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_j[iSpecies][iDim] = U_j[loc+iDim+1] / Density_j[iSpecies];
      ProjVelocity_j[iSpecies]  += Velocity_j[iSpecies][iDim] *UnitaryNormal[iDim];
			Vel2                      += Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
		}
		Energy_j[iSpecies]		= U_j[loc+nDim+1] / Density_j[iSpecies];
		Energy_vib_j[iSpecies] = 0.0;
		Energy_el_j[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_j[iSpecies]  = U_j[loc+nDim+2] / Density_j[iSpecies];
			SoundSpeed_j[iSpecies] 	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_j[iSpecies] - 0.5*Vel2 - Energy_vib_j[iSpecies] - Energy_el_j[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_j[iSpecies] = (GammaDiatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_vib_j[iSpecies] - Energy_el_j[iSpecies]);
		}
		else {
			SoundSpeed_j[iSpecies] 	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_j[iSpecies] - 0.5*Vel2 - Energy_el_j[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			Pressure_j[iSpecies] = (GammaMonatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el_j[iSpecies]);
		}
		Enthalpy_j[iSpecies] = Energy_j[iSpecies] + Pressure_j[iSpecies] / Density_j[iSpecies];
	}
  
  /*--- Calculate weighted state vector ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
    /*--- Calculate the weighting function ---*/
    dp = fabs(Pressure_j[iSpecies] - Pressure_i[iSpecies]) / min(Pressure_j[iSpecies],Pressure_i[iSpecies]);
    w = 0.5 * (1.0/(pow(alpha*dp,2.0) +1.0));
    onemw = 1.0 - w;
    
    /*--- Calculate weighted state vector, U_star, for i ---*/
    Densityst_i[iSpecies] = onemw*Density_i[iSpecies] + w*Density_j[iSpecies];
    Densityst_j[iSpecies] = onemw*Density_j[iSpecies] + w*Density_i[iSpecies];
    
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocityst_i[iSpecies][iDim] = onemw*Velocity_i[iSpecies][iDim] + w*Velocity_j[iSpecies][iDim];
      Velocityst_j[iSpecies][iDim] = onemw*Velocity_j[iSpecies][iDim] + w*Velocity_i[iSpecies][iDim];
      ProjVelocityst_i[iSpecies] = onemw*ProjVelocity_i[iSpecies] + w*ProjVelocity_j[iSpecies];
      ProjVelocityst_j[iSpecies] = onemw*ProjVelocity_j[iSpecies] + w*ProjVelocity_i[iSpecies];
    }
    
    Enthalpyst_i[iSpecies] = onemw*Enthalpy_i[iSpecies] + w*Enthalpy_j[iSpecies];
    Enthalpyst_j[iSpecies] = onemw*Enthalpy_j[iSpecies] + w*Enthalpy_i[iSpecies];
    
    Soundspeedst_i[iSpecies] = onemw*SoundSpeed_i[iSpecies] + w*SoundSpeed_j[iSpecies];
    Soundspeedst_j[iSpecies] = onemw*SoundSpeed_j[iSpecies] + w*SoundSpeed_i[iSpecies];
    
    Energy_vibst_i[iSpecies] = onemw*Energy_vib_i[iSpecies] + w*Energy_vib_j[iSpecies];
    Energy_vibst_j[iSpecies] = onemw*Energy_vib_j[iSpecies] + w*Energy_vib_i[iSpecies];
    
    Energy_elst_i[iSpecies] = onemw*Energy_el_i[iSpecies] + w*Energy_el_j[iSpecies];
    Energy_elst_j[iSpecies] = onemw*Energy_el_j[iSpecies] + w*Energy_el_i[iSpecies];
  }
  
	/*--- Flow eigenvalues at i (Lambda+) --- */
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda_i[loc+iDim] = 0.5*(ProjVelocityst_i[iSpecies] + fabs(ProjVelocityst_i[iSpecies]));
		}
		Lambda_i[loc+nDim]   = 0.5*(ProjVelocityst_i[iSpecies] + Soundspeedst_i[iSpecies]
                                + fabs(ProjVelocityst_i[iSpecies] + Soundspeedst_i[iSpecies]));
		Lambda_i[loc+nDim+1] = 0.5*(ProjVelocityst_i[iSpecies] - Soundspeedst_i[iSpecies]
                                + fabs(ProjVelocityst_i[iSpecies] - Soundspeedst_i[iSpecies]));
		if (iSpecies < nDiatomics) Lambda_i[loc+nDim+2] = 0.5*(ProjVelocityst_i[iSpecies] + fabs(ProjVelocityst_i[iSpecies]));
	}
  
	/*--- Compute P & invP at i ---*/
	GetPMatrix_(Densityst_i, Velocityst_i, Enthalpyst_i, Soundspeedst_i, Energy_vibst_i, Energy_elst_i, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Densityst_i, Velocityst_i, Soundspeedst_i, Energy_vibst_i, Energy_elst_i, config, UnitaryNormal, invP_Tensor);
  
	/*--- Projected flux (f+) at i ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_i = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_i += P_Tensor[iVar][kVar]*Lambda_i[kVar]*invP_Tensor[kVar][jVar];
			Proj_flux_tensor_i[iVar] += Proj_ModJac_Tensor_i*U_i[jVar]*Area;
			val_Jacobian_i[iVar][jVar] += Proj_ModJac_Tensor_i*Area;
		}
	}
  
	/*--- Flow eigenvalues at j (Lambda-) --- */
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda_j[loc+iDim] = 0.5*(ProjVelocityst_j[iSpecies] - fabs(ProjVelocityst_j[iSpecies]));
		}
		Lambda_j[loc + nDim]   = 0.5*(ProjVelocityst_j[iSpecies] + Soundspeedst_j[iSpecies]
                                  - fabs(ProjVelocityst_j[iSpecies] + Soundspeedst_j[iSpecies]));
		Lambda_j[loc + nDim+1] = 0.5*(ProjVelocityst_j[iSpecies] - Soundspeedst_j[iSpecies]
                                  - fabs(ProjVelocityst_j[iSpecies] - Soundspeedst_j[iSpecies]));
		if (iSpecies < nDiatomics) Lambda_j[loc + nDim+2] = 0.5*(ProjVelocityst_i[iSpecies] - fabs(ProjVelocityst_i[iSpecies]));
	}
  
	/*--- Compute P & invP at j ---*/
	GetPMatrix_(Densityst_j, Velocityst_j, Enthalpyst_j, Soundspeedst_j, Energy_vibst_j, Energy_elst_j, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Densityst_j, Velocityst_j, Soundspeedst_j, Energy_vibst_j, Energy_elst_j, config, UnitaryNormal, invP_Tensor);
  
	/*--- Projected flux (f-) ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_j = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_j += P_Tensor[iVar][kVar]*Lambda_j[kVar]*invP_Tensor[kVar][jVar];
			Proj_flux_tensor_j[iVar] += Proj_ModJac_Tensor_j*U_j[jVar]*Area;
			val_Jacobian_j[iVar][jVar] += Proj_ModJac_Tensor_j*Area;
		}
	}
  
	/*--- Flux splitting ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar];
	}
}



CUpwHLLC_PlasmaDiatomic::CUpwHLLC_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
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

CUpwHLLC_PlasmaDiatomic::~CUpwHLLC_PlasmaDiatomic(void) {
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

void CUpwHLLC_PlasmaDiatomic::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	/*--- Point i, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	Density_i = U_i[0];
	double sq_vel_i = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		sq_vel_i += Velocity_i[iDim]*Velocity_i[iDim];
	}
	Energy_i = U_i[nDim+1] / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel_i));
	Pressure_i = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	Enthalpy_i = (U_i[nDim+1] + Pressure_i) / Density_i;
  
	/*--- Point j, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	Density_j = U_j[0];
	double sq_vel_j = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel_j += Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_j = U_j[nDim+1] / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel_j));
	Pressure_j = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	Enthalpy_j = (U_j[nDim+1] + Pressure_j) / Density_j;
  
	/*--- Projected velocities ---*/
	ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
	}
  
	/*--- Roe's aveaging ---*/
	double Rrho = sqrt(Density_j/Density_i);
	double tmp = 1.0/(1.0+Rrho);
	double velRoe[3];
	for (iDim = 0; iDim < nDim; iDim++)
		velRoe[iDim] = tmp*(Velocity_i[iDim] + Velocity_j[iDim]*Rrho);
  
	double uRoe  = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		uRoe += velRoe[iDim]*UnitaryNormal[iDim];
  
	double gamPdivRho = tmp*((Gamma*Pressure_i/Density_i+0.5*(Gamma-1.0)*sq_vel_i) + (Gamma*Pressure_j/Density_j+0.5*(Gamma-1.0)*sq_vel_j)*Rrho);
	double sq_velRoe = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		sq_velRoe += velRoe[iDim]*velRoe[iDim];
  
	double cRoe  = sqrt(gamPdivRho - ((Gamma+Gamma)*0.5-1.0)*0.5*sq_velRoe);
  
	/*--- speed of sound at L and R ---*/
	double sL = min(uRoe-cRoe, ProjVelocity_i-SoundSpeed_i);
	double sR = max(uRoe+cRoe, ProjVelocity_j+SoundSpeed_j);
  
	/*--- speed of contact surface ---*/
	double sM = (Pressure_i-Pressure_j
               - Density_i*ProjVelocity_i*(sL-ProjVelocity_i)
               + Density_j*ProjVelocity_j*(sR-ProjVelocity_j))
  /(Density_j*(sR-ProjVelocity_j)-Density_i*(sL-ProjVelocity_i));
  
	/*--- Pressure at right and left (Pressure_j=Pressure_i) side of contact surface ---*/
	double pStar = Density_j * (ProjVelocity_j-sR)*(ProjVelocity_j-sM) + Pressure_j;
  
	if (sM >= 0.0) {
		if (sL > 0.0) {
			val_residual[0] = Density_i*ProjVelocity_i;
			for (iDim = 0; iDim < nDim; iDim++)
				val_residual[iDim+1] = Density_i*Velocity_i[iDim]*ProjVelocity_i + Pressure_i*UnitaryNormal[iDim];
			val_residual[nVar-1] = Energy_i*Density_i*ProjVelocity_i + Pressure_i*ProjVelocity_i;
		}
		else {
			double invSLmSs = 1.0/(sL-sM);
			double sLmuL = sL-ProjVelocity_i;
			double rhoSL = Density_i*sLmuL*invSLmSs;
			double rhouSL[3];
			for (iDim = 0; iDim < nDim; iDim++)
				rhouSL[iDim] = (Density_i*Velocity_i[iDim]*sLmuL+(pStar-Pressure_i)*UnitaryNormal[iDim])*invSLmSs;
			double eSL = (sLmuL*Energy_i*Density_i-Pressure_i*ProjVelocity_i+pStar*sM)*invSLmSs;
      
			val_residual[0] = rhoSL*sM;
			for (iDim = 0; iDim < nDim; iDim++)
				val_residual[iDim+1] = rhouSL[iDim]*sM + pStar*UnitaryNormal[iDim];
			val_residual[nVar-1] = (eSL+pStar)*sM;
		}
	}
	else {
		if (sR >= 0.0) {
			double invSRmSs = 1.0/(sR-sM);
			double sRmuR = sR-ProjVelocity_j;
			double rhoSR = Density_j*sRmuR*invSRmSs;
			double rhouSR[3];
			for (iDim = 0; iDim < nDim; iDim++)
				rhouSR[iDim] = (Density_j*Velocity_j[iDim]*sRmuR+(pStar-Pressure_j)*UnitaryNormal[iDim])*invSRmSs;
			double eSR = (sRmuR*Energy_j*Density_j-Pressure_j*ProjVelocity_j+pStar*sM)*invSRmSs;
      
			val_residual[0] = rhoSR*sM;
			for (iDim = 0; iDim < nDim; iDim++)
				val_residual[iDim+1] = rhouSR[iDim]*sM + pStar*UnitaryNormal[iDim];
			val_residual[nVar-1] = (eSR+pStar)*sM;
		}
		else {
			val_residual[0] = Density_j*ProjVelocity_j;
			for (iDim = 0; iDim < nDim; iDim++)
				val_residual[iDim+1] = Density_j*Velocity_j[iDim]*ProjVelocity_j + Pressure_j*UnitaryNormal[iDim];
			val_residual[nVar-1] = Energy_j*Density_j*ProjVelocity_j + Pressure_j*ProjVelocity_j;
		}
	}
  
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
    cout << "Need to change the call to InvProjJacobians, to InvProjJacobians_ and include effects of diatomic species" << endl;
    cin.get();
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

CCentJST_Plasma::CCentJST_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Artifical dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_Flow();
	Param_Kappa_4 = config->GetKappa_4th_Flow();
  
	/*--- Allocate some structures ---*/
	Diff_U = new double [nVar];
	Diff_Lapl = new double [nVar];
	Velocity_i = new double* [nSpecies];
	Velocity_j = new double* [nSpecies];
  
	Proj_flux_tensor = new double [nVar];
  
	Pressure_i = new double [nSpecies];
	Pressure_j = new double [nSpecies];
  
	MeanEnergy = new double [nSpecies];
	MeanDensity = new double [nSpecies];
	MeanPressure = new double [nSpecies];
	MeanEnthalpy = new double [nSpecies];
	MeanVelocity = new double* [nSpecies];
	MeanLambda   = new double [nSpecies];
	StretchingFactor = new double [nSpecies];
	Epsilon_2 = new double [nSpecies];
	Epsilon_4 = new double [nSpecies];
  
	for(iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		MeanVelocity[iSpecies] = new double [nDim];
		Velocity_i[iSpecies] = new double [nDim];
		Velocity_j[iSpecies] = new double [nDim];
    
	}
  
	SoundSpeed_i = new double [nSpecies];
	SoundSpeed_j = new double [nSpecies];
  
	Enthalpy_i = new double [nSpecies];
	Enthalpy_j = new double [nSpecies];
  
	Lambda_i = new double [nSpecies];
	Lambda_j = new double [nSpecies];
  
	Sensor_i = new double [nSpecies];
	Sensor_j = new double [nSpecies];
  
}

CCentJST_Plasma::~CCentJST_Plasma(void) {
  
	delete [] Diff_U; 	delete [] Diff_Lapl;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		delete [] MeanVelocity[iSpecies];
		delete [] Velocity_i[iSpecies];
		delete [] Velocity_j[iSpecies];
    
	}
	delete [] Velocity_i;		delete [] Velocity_j;
	delete [] Pressure_i;		delete [] Pressure_j;
	delete [] MeanVelocity;		delete [] Proj_flux_tensor;
	delete [] MeanEnergy;		delete [] MeanDensity;	delete [] 	MeanPressure;
	delete [] MeanEnthalpy;		delete [] MeanVelocity;
	delete [] MeanLambda;		delete [] StretchingFactor;
	delete [] Epsilon_2;		delete [] Epsilon_4;
	delete [] SoundSpeed_i;		delete [] SoundSpeed_j ;
	delete [] Enthalpy_i;		delete [] Enthalpy_j;
	delete [] Lambda_i;			delete [] Lambda_j;
	delete [] Sensor_i;			delete [] Sensor_j;
  
}

void CCentJST_Plasma::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
                                  CConfig *config) {
  
  
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*--- Conservative variables at point i and 1 ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = iSpecies*(nDim+2);
		Density_i	= U_i[loc + 0];
		Density_j	= U_j[loc + 0];
		MeanDensity[iSpecies]	 = 0.5*(Density_i+Density_j);
    
		Energy_i = U_i[loc + nDim+1] / Density_i;
		Energy_j = U_j[loc + nDim+1] / Density_j;
		MeanEnergy[iSpecies] = 0.5*(Energy_i+Energy_j);
    
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i;
			Velocity_j[iSpecies][iDim] = U_j[loc + iDim+1] / Density_j;
			MeanVelocity[iSpecies][iDim] =  0.5*(Velocity_i[iSpecies][iDim]+Velocity_j[iSpecies][iDim]);
		}
		MeanPressure[iSpecies] = 0.5*(Pressure_i[iSpecies]+Pressure_j[iSpecies]);
		MeanEnthalpy[iSpecies] = 0.5*(Enthalpy_i[iSpecies]+Enthalpy_j[iSpecies]);
    
	}
  
	/*--- Get projected flux tensor ---*/
	GetInviscidProjFlux(MeanDensity, MeanVelocity, MeanPressure, MeanEnthalpy, Normal, Proj_flux_tensor);
  
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		GetInviscidProjJac(MeanVelocity, MeanEnergy, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}
  
	/*--- Computes differences btw. Laplacians and conservative variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_U[iVar] = U_i[iVar]-U_j[iVar];
	}
  
	nVar_Species = (nDim+2);
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = iSpecies*(nDim+2);
		Density_i	= U_i[loc + 0];
		Density_j	= U_j[loc + 0];
		Diff_U[loc + nVar_Species-1] = Density_i*Enthalpy_i[iSpecies]-Density_j*Enthalpy_j[iSpecies];
	}
  
	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
  
	/*--- Compute the local espectral radius and the stretching factor ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
    
		ProjVelocity_i = 0; ProjVelocity_j = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity_i += Velocity_i[iSpecies][iDim]*Normal[iDim];
			ProjVelocity_j += Velocity_j[iSpecies][iDim]*Normal[iDim];
		}
    
		Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i[iSpecies]*Area);
		Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j[iSpecies]*Area);
		MeanLambda[iSpecies] = 0.5*(Local_Lambda_i+Local_Lambda_j);
		Phi_i = pow(Lambda_i[iSpecies]/(4.0*MeanLambda[iSpecies]+EPS), Param_p);
		Phi_j = pow(Lambda_j[iSpecies]/(4.0*MeanLambda[iSpecies]+EPS), Param_p);
		StretchingFactor[iSpecies] = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
		Epsilon_2[iSpecies] = Param_Kappa_2*0.5*(Sensor_i[iSpecies]+Sensor_j[iSpecies])*sc2;
		Epsilon_4[iSpecies] = max(0.0, Param_Kappa_4-Epsilon_2[iSpecies])*sc4;
	}
  
	/*--- Compute viscous part of the residual ---*/
	nVar_Species = (nDim+2);
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = iSpecies*(nDim+2);
		for (iVar = 0; iVar < nVar_Species; iVar++)
			val_resvisc[loc + iVar] = (Epsilon_2[iSpecies]*Diff_U[loc + iVar] - Epsilon_4[iSpecies]*Diff_Lapl[loc + iVar])*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
	}
  
	if (implicit) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			loc = iSpecies*(nDim+2);
			cte_0 = (Epsilon_2[iSpecies] + Epsilon_4[iSpecies]*double(Neighbor_i+1))*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			cte_1 = (Epsilon_2[iSpecies] + Epsilon_4[iSpecies]*double(Neighbor_j+1))*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			for (iVar = 0; iVar < (nVar_Species-1); iVar++) {
				val_Jacobian_i[loc + iVar][loc + iVar] += cte_0;
				val_Jacobian_j[loc + iVar][loc + iVar] -= cte_1;
			}
			sq_vel_i = 0.0; sq_vel_j = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				sq_vel_i += 0.5*Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
				sq_vel_j += 0.5*Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
			}
			Gamma = Vector_Gamma[iSpecies];
			Gamma_Minus_One = Gamma - 1.0;
      
			val_Jacobian_i[loc + nVar_Species-1][loc + 0] += cte_0*Gamma_Minus_One*sq_vel_i;
			for (iDim = 0; iDim < nDim; iDim++)
				val_Jacobian_i[loc + nVar_Species-1][loc + iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iSpecies][iDim];
			val_Jacobian_i[loc + nVar_Species-1][loc + nVar_Species-1] += cte_0*Gamma;
      
      
			/*--- Last row of Jacobian_j ---*/
			val_Jacobian_j[loc + nVar_Species-1][loc + 0] -= cte_1*Gamma_Minus_One*sq_vel_j;
			for (iDim = 0; iDim < nDim; iDim++)
				val_Jacobian_j[loc + nVar_Species-1][loc + iDim+1] += cte_1*Gamma_Minus_One*Velocity_j[iSpecies][iDim];
			val_Jacobian_j[loc + nVar_Species-1][loc + nVar_Species-1] -= cte_1*Gamma;
		}
	}
}


CCentJST_PlasmaDiatomic::CCentJST_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Artifical dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_Flow();
	Param_Kappa_4 = config->GetKappa_4th_Flow();
  
	/*--- Allocate some structures ---*/
	Diff_U = new double [nVar];
	Diff_Lapl = new double [nVar];
	Velocity_i = new double* [nSpecies];
	Velocity_j = new double* [nSpecies];
  
	Proj_flux_tensor = new double [nVar];
  
	Pressure_i = new double [nSpecies];
	Pressure_j = new double [nSpecies];
  
	MeanEnergy = new double [nSpecies];
	MeanEnergy_vib = new double [nSpecies];
	MeanDensity = new double [nSpecies];
	MeanPressure = new double [nSpecies];
	MeanEnthalpy = new double [nSpecies];
	MeanVelocity = new double* [nSpecies];
	MeanLambda   = new double [nSpecies];
	StretchingFactor = new double [nSpecies];
	Epsilon_2 = new double [nSpecies];
	Epsilon_4 = new double [nSpecies];
  
	for(iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		MeanVelocity[iSpecies] = new double [nDim];
		Velocity_i[iSpecies] = new double [nDim];
		Velocity_j[iSpecies] = new double [nDim];
    
	}
  
	SoundSpeed_i = new double [nSpecies];
	SoundSpeed_j = new double [nSpecies];
  
	Enthalpy_i = new double [nSpecies];
	Enthalpy_j = new double [nSpecies];
  
	Lambda_i = new double [nSpecies];
	Lambda_j = new double [nSpecies];
  
	Sensor_i = new double [nSpecies];
	Sensor_j = new double [nSpecies];
  
}

CCentJST_PlasmaDiatomic::~CCentJST_PlasmaDiatomic(void) {
  
	delete [] Diff_U; 	delete [] Diff_Lapl;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		delete [] MeanVelocity[iSpecies];
		delete [] Velocity_i[iSpecies];
		delete [] Velocity_j[iSpecies];
    
	}
	delete [] Velocity_i;		delete [] Velocity_j;
	delete [] Pressure_i;		delete [] Pressure_j;
	delete [] MeanVelocity;		delete [] Proj_flux_tensor;
	delete [] MeanEnergy;		delete [] MeanDensity;	delete [] 	MeanPressure;
	delete [] MeanEnthalpy;		delete [] MeanVelocity;
	delete [] MeanLambda;		delete [] StretchingFactor;
	delete [] Epsilon_2;		delete [] Epsilon_4;
	delete [] SoundSpeed_i;		delete [] SoundSpeed_j ;
	delete [] Enthalpy_i;		delete [] Enthalpy_j;
	delete [] Lambda_i;			delete [] Lambda_j;
	delete [] Sensor_i;			delete [] Sensor_j;
  
}

void CCentJST_PlasmaDiatomic::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
                                          CConfig *config) {
  
  
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*--- Conservative variables at point i and 1 ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		Density_i	= U_i[loc + 0];
		Density_j	= U_j[loc + 0];
		MeanDensity[iSpecies]	 = 0.5*(Density_i+Density_j);
    
		Energy_i = U_i[loc + nDim+1] / Density_i;
		Energy_j = U_j[loc + nDim+1] / Density_j;
		MeanEnergy[iSpecies] = 0.5*(Energy_i+Energy_j);
    
		Energy_vib_i = 0.0;
		Energy_vib_j = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_i = U_i[loc+nDim+2] / Density_i;
			Energy_vib_j = U_j[loc+nDim+2] / Density_j;
		}
		MeanEnergy_vib[iSpecies] = 0.5*(Energy_vib_i+Energy_vib_j);
    
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i;
			Velocity_j[iSpecies][iDim] = U_j[loc + iDim+1] / Density_j;
			MeanVelocity[iSpecies][iDim] =  0.5*(Velocity_i[iSpecies][iDim]+Velocity_j[iSpecies][iDim]);
		}
    
		MeanPressure[iSpecies] = 0.5*(Pressure_i[iSpecies]+Pressure_j[iSpecies]);
		MeanEnthalpy[iSpecies] = 0.5*(Enthalpy_i[iSpecies]+Enthalpy_j[iSpecies]);
    
	}
  
	/*--- Get projected flux tensor ---*/
	GetInviscidProjFlux_(MeanDensity, MeanVelocity, MeanPressure, MeanEnthalpy, MeanEnergy_vib, Normal, Proj_flux_tensor);
  
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		GetInviscidProjJac_(MeanVelocity, MeanEnergy, MeanEnergy_vib, MeanEnthalpy, Normal, 0.5, val_Jacobian_i, config);
		//		GetInviscidProjJac_(MeanVelocity, MeanEnergy, Normal, 0.5, val_Jacobian_i);
		//		GetInviscidProjJac(MeanVelocity, MeanEnergy, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}
  
	/*--- Computes differences btw. Laplacians and conservative variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_U[iVar] = U_i[iVar]-U_j[iVar];
	}
  
	//	nVar_Species = (nDim+2);
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		Density_i	= U_i[loc + 0];
		Density_j	= U_j[loc + 0];
		Diff_U[loc+nDim+1] = Density_i*Enthalpy_i[iSpecies]-Density_j*Enthalpy_j[iSpecies];
	}
  
	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
  
	/*--- Compute the local espectral radius and the stretching factor ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    
		ProjVelocity_i = 0; ProjVelocity_j = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity_i += Velocity_i[iSpecies][iDim]*Normal[iDim];
			ProjVelocity_j += Velocity_j[iSpecies][iDim]*Normal[iDim];
		}
		Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i[iSpecies]*Area);
		Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j[iSpecies]*Area);
		MeanLambda[iSpecies] = 0.5*(Local_Lambda_i+Local_Lambda_j);
		Phi_i = pow(Lambda_i[iSpecies]/(4.0*MeanLambda[iSpecies]+EPS), Param_p);
		Phi_j = pow(Lambda_j[iSpecies]/(4.0*MeanLambda[iSpecies]+EPS), Param_p);
		StretchingFactor[iSpecies] = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
		Epsilon_2[iSpecies] = Param_Kappa_2*0.5*(Sensor_i[iSpecies]+Sensor_j[iSpecies])*sc2;
		Epsilon_4[iSpecies] = max(0.0, Param_Kappa_4-Epsilon_2[iSpecies])*sc4;
	}
  
	/*--- Compute viscous part of the residual ---*/
	//	nVar_Species = (nDim+2);
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (iVar = 0; iVar < nVar_Species; iVar++)
			val_resvisc[loc + iVar] = (Epsilon_2[iSpecies]*Diff_U[loc + iVar] - Epsilon_4[iSpecies]*Diff_Lapl[loc + iVar])*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
	}
  
	if (implicit) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			cte_0 = (Epsilon_2[iSpecies] + Epsilon_4[iSpecies]*double(Neighbor_i+1))*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			cte_1 = (Epsilon_2[iSpecies] + Epsilon_4[iSpecies]*double(Neighbor_j+1))*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			for (iVar = 0; iVar < (nVar_Species-1); iVar++) {
				val_Jacobian_i[loc + iVar][loc + iVar] += cte_0;
				val_Jacobian_j[loc + iVar][loc + iVar] -= cte_1;
			}
			sq_vel_i = 0.0; sq_vel_j = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				sq_vel_i += 0.5*Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
				sq_vel_j += 0.5*Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
			}
			Gamma = Vector_Gamma[iSpecies];
			Gamma_Minus_One = Gamma - 1.0;
      
			val_Jacobian_i[loc+nDim+1][loc+0] += cte_0*Gamma_Minus_One*sq_vel_i;
			for (iDim = 0; iDim < nDim; iDim++)
				val_Jacobian_i[loc+nDim+1][loc+iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iSpecies][iDim];
			val_Jacobian_i[loc+nDim+1][loc+nDim+1] += cte_0*Gamma;
      
      
			/*--- Last row of Jacobian_j ---*/
			val_Jacobian_j[loc + nDim+1][loc+0] -= cte_1*Gamma_Minus_One*sq_vel_j;
			for (iDim = 0; iDim < nDim; iDim++)
				val_Jacobian_j[loc+nDim+1][loc+iDim+1] += cte_1*Gamma_Minus_One*Velocity_j[iSpecies][iDim];
			val_Jacobian_j[loc+nDim+1][loc+nDim+1] -= cte_1*Gamma;
		}
	}
}

CCentLax_PlasmaDiatomic::CCentLax_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	/*--- Artifical dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_Plasma();
  
	/*--- Allocate some structures ---*/
	Diff_U = new double [nVar];
	Velocity_i = new double* [nSpecies];
	Velocity_j = new double* [nSpecies];
	MeanVelocity = new double* [nSpecies];
	for(iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		MeanVelocity[iSpecies] = new double [nDim];
		Velocity_i[iSpecies] = new double [nDim];
		Velocity_j[iSpecies] = new double [nDim];
	}
  
	Proj_flux_tensor = new double [nVar];
  
	Pressure_i = new double [nSpecies];
	Pressure_j = new double [nSpecies];
  
	Energy_i = new double[nSpecies];
	Energy_j = new double[nSpecies];
	Energy_vib_i = new double[nSpecies];
	Energy_vib_j = new double[nSpecies];
	MeanEnergy = new double [nSpecies];
	MeanEnergy_vib = new double [nSpecies];
	MeanDensity = new double [nSpecies];
	MeanPressure = new double [nSpecies];
	MeanEnthalpy = new double [nSpecies];
	MeanLambda   = new double [nSpecies];
	StretchingFactor = new double [nSpecies];
  
	Epsilon_0 = new double [nSpecies];
  
	SoundSpeed_i = new double [nSpecies];
	SoundSpeed_j = new double [nSpecies];
  
	Enthalpy_i = new double [nSpecies];
	Enthalpy_j = new double [nSpecies];
  
	Lambda_i = new double [nSpecies];
	Lambda_j = new double [nSpecies];
  
	Sensor_i = new double [nSpecies];
	Sensor_j = new double [nSpecies];
}

CCentLax_PlasmaDiatomic::~CCentLax_PlasmaDiatomic(void) {
  
	delete [] Diff_U; 	delete [] Diff_Lapl;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		delete [] MeanVelocity[iSpecies];
		delete [] Velocity_i[iSpecies];
		delete [] Velocity_j[iSpecies];
    
	}
	delete [] Velocity_i;		delete [] Velocity_j;
	delete [] Pressure_i;		delete [] Pressure_j;
	delete [] Energy_i;			delete [] Energy_j;
	delete [] Energy_vib_i; delete [] Energy_vib_j;
	delete [] MeanVelocity;		delete [] Proj_flux_tensor;
	delete [] MeanEnergy;		delete [] MeanDensity;	delete [] 	MeanPressure;
	delete [] MeanEnthalpy;		delete [] MeanVelocity;
	delete [] MeanLambda;		delete [] StretchingFactor;
	delete [] Epsilon_0;
	delete [] SoundSpeed_i;		delete [] SoundSpeed_j ;
	delete [] Enthalpy_i;		delete [] Enthalpy_j;
	delete [] Lambda_i;			delete [] Lambda_j;
	delete [] Sensor_i;			delete [] Sensor_j;
  
}

void CCentLax_PlasmaDiatomic::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
                                          CConfig *config) {
  
	/*--- Conservative variables at point i and 1 ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		Density_i	= U_i[loc + 0];
		Density_j	= U_j[loc + 0];
		MeanDensity[iSpecies]	 = 0.5*(Density_i+Density_j);
    
		Energy_i[iSpecies] = U_i[loc + nDim+1] / Density_i;
		Energy_j[iSpecies] = U_j[loc + nDim+1] / Density_j;
		MeanEnergy[iSpecies] = 0.5*(Energy_i[iSpecies]+Energy_j[iSpecies]);
    
		Energy_vib_i[iSpecies] = 0.0;
		Energy_vib_j[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_i[iSpecies] = U_i[loc+nDim+2] / Density_i;
			Energy_vib_j[iSpecies] = U_j[loc+nDim+2] / Density_j;
		}
		MeanEnergy_vib[iSpecies] = 0.5*(Energy_vib_i[iSpecies]+Energy_vib_j[iSpecies]);
    
		sq_vel_i = 0.0; sq_vel_j = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i;
			Velocity_j[iSpecies][iDim] = U_j[loc + iDim+1] / Density_j;
			sq_vel_i += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
			sq_vel_j += Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
			MeanVelocity[iSpecies][iDim] =  0.5*(Velocity_i[iSpecies][iDim]+Velocity_j[iSpecies][iDim]);
		}
		double Gamma_minus_one;
		Gamma_minus_one = config->GetSpecies_Gamma(iSpecies) - 1.0;
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Energy_el_i = 0.0;
		Energy_el_j = 0.0;
    
		Pressure_i[iSpecies] = (Gamma_minus_one)*Density_i*(Energy_i[iSpecies] - 0.5*sq_vel_i - Enthalpy_formation - Energy_vib_i[iSpecies] - Energy_el_i);
		Pressure_j[iSpecies] = (Gamma_minus_one)*Density_j*(Energy_j[iSpecies] - 0.5*sq_vel_j - Enthalpy_formation - Energy_vib_j[iSpecies] - Energy_el_j);
		MeanPressure[iSpecies] = 0.5*(Pressure_i[iSpecies]+Pressure_j[iSpecies]);
    
		Enthalpy_i[iSpecies] = Energy_i[iSpecies] + Pressure_i[iSpecies]/Density_i;
		Enthalpy_j[iSpecies] = Energy_j[iSpecies] + Pressure_j[iSpecies]/Density_j;
		MeanEnthalpy[iSpecies] = 0.5*(Enthalpy_i[iSpecies]+Enthalpy_j[iSpecies]);
	}
  
	/*--- Get projected flux tensor ---*/
	GetInviscidProjFlux_(MeanDensity, MeanVelocity, MeanPressure, MeanEnthalpy, MeanEnergy_vib, Normal, Proj_flux_tensor);
  
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		GetInviscidProjJac_(MeanVelocity, MeanEnergy, MeanEnergy_vib, MeanEnthalpy, Normal, 0.5, val_Jacobian_i, config);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}
  
	/*--- Computes differences btw. conservative variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_U[iVar] = U_i[iVar]-U_j[iVar];
	}
	for (iSpecies = 0; iSpecies<nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		Diff_U[loc+nDim+1] = U_i[loc+0]*Enthalpy_i[iSpecies] - U_j[loc+0]*Enthalpy_j[iSpecies];
	}
  
	sc0 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
  
	/*--- Compute the local spectral radius and the stretching factor ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		ProjVelocity_i = 0; ProjVelocity_j = 0; Area = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity_i += Velocity_i[iSpecies][iDim]*Normal[iDim];
			ProjVelocity_j += Velocity_j[iSpecies][iDim]*Normal[iDim];
			Area += Normal[iDim]*Normal[iDim];
		}
		Area = sqrt(Area);
		Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i[iSpecies]*Area);
		Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j[iSpecies]*Area);
    
		MeanLambda[iSpecies] = 0.5*(Local_Lambda_i+Local_Lambda_j);
		Phi_i = pow(Lambda_i[iSpecies]/(4.0*MeanLambda[iSpecies]+EPS), Param_p);
		Phi_j = pow(Lambda_j[iSpecies]/(4.0*MeanLambda[iSpecies]+EPS), Param_p);
    
		StretchingFactor[iSpecies] = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
		Epsilon_0[iSpecies] = Param_Kappa_0*sc0*double(nDim)/3.0;
	}
  
	/*--- Compute viscous part of the residual ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
    
		if ( iSpecies < nDiatomics ) {
			loc = (nDim+3)*iSpecies;
			nVar_Species = nDim+3;
		} else {
			loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			nVar_Species = nDim+2;
		}
    
		for (iVar = 0; iVar < nVar_Species; iVar++)
			val_resvisc[loc + iVar] = (Epsilon_0[iSpecies]*Diff_U[loc + iVar])*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
	}
  
	if (implicit) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			if ( iSpecies < nDiatomics ) {
				loc = (nDim+3)*iSpecies;
				nVar_Species = nDim+3;
			} else {
				loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				nVar_Species = nDim+2;
			}
			cte_0 = Epsilon_0[iSpecies]*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			for (iVar = 0; iVar < nDim+1; iVar++) {
				val_Jacobian_i[loc+iVar][loc+iVar] += cte_0;
				val_Jacobian_j[loc+iVar][loc+iVar] -= cte_0;
			}
      
			sq_vel_i = 0.0;
			sq_vel_j = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				sq_vel_i += 0.5*Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
				sq_vel_j += 0.5*Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
			}
      
			/*--- Energy rows: CAREFUL!! You have differences of \rho_Enthalpy, not differences of \rho_Energy ---*/
			val_Jacobian_i[loc+nDim+1][loc+0] += cte_0*(Vector_Gamma[iSpecies]-1.0)*(sq_vel_i-config->GetEnthalpy_Formation(iSpecies));
			for (iDim = 0; iDim < nDim; iDim++)
				val_Jacobian_i[loc+nDim+1][loc+iDim+1] -= cte_0*(Vector_Gamma[iSpecies]-1.0)*Velocity_i[iSpecies][iDim];
			val_Jacobian_i[loc+nDim+1][loc+nDim+1] += cte_0*Vector_Gamma[iSpecies];
			if (iSpecies < nDiatomics)
				val_Jacobian_i[loc+nDim+1][loc+nDim+2] -= cte_0*(Vector_Gamma[iSpecies]-1.0);
      
			/*--- Energy row of Jacobian_j ---*/
			val_Jacobian_j[loc+nDim+1][loc+0] -= cte_0*(Vector_Gamma[iSpecies]-1.0)*(sq_vel_j-config->GetEnthalpy_Formation(iSpecies));
			for (iDim = 0; iDim < nDim; iDim++)
				val_Jacobian_j[loc+nDim+1][loc+iDim+1] += cte_0*(Vector_Gamma[iSpecies]-1.0)*Velocity_j[iSpecies][iDim];
			val_Jacobian_j[loc+nDim+1][loc+nDim+1] -= cte_0*Vector_Gamma[iSpecies];
			if (iSpecies < nDiatomics)
				val_Jacobian_j[loc+nDim+1][loc+nDim+2] += cte_0*(Vector_Gamma[iSpecies]-1.0);
      
			/*--- Vibrational Energy row of Jacobian_i & Jacobian_j ---*/
			if (iSpecies < nDiatomics) {
				val_Jacobian_i[loc+nDim+2][loc+nDim+2] += cte_0;
				val_Jacobian_j[loc+nDim+2][loc+nDim+2] -= cte_0;
			}
		}
	}
}

CAvgGrad_Plasma::CAvgGrad_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies,unsigned short val_nDiatomics, unsigned short val_nMonatomics,  CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Allocate some useful vectors ---*/
  
	PrimVar_i        = new double [nDim+6];
	PrimVar_j        = new double [nDim+6];
  Mean_PrimVar     = new double [nDim+6];
  Mean_GradPrimVar = new double*[nDim+3];
	for (iVar = 0; iVar < nDim+3; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
  
	Mean_Laminar_Viscosity         = new double[nSpecies];
	Mean_Eddy_Viscosity            = new double[nSpecies];
  Mean_Thermal_Conductivity      = new double[nSpecies];
  Mean_Thermal_Conductivity_vib  = new double[nSpecies];
}

CAvgGrad_Plasma::~CAvgGrad_Plasma(void) {
  
	delete [] PrimVar_i;
	delete [] PrimVar_j;
	delete [] Mean_PrimVar;
	for (iVar = 0; iVar < nDim+3; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
  
  delete [] Mean_Laminar_Viscosity;
  delete [] Mean_Eddy_Viscosity;
  delete [] Mean_Thermal_Conductivity;
  delete [] Mean_Thermal_Conductivity_vib;
}

void CAvgGrad_Plasma::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	unsigned short loc, iVar, jVar, nVar_species, iSpecies;
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar =0 ; jVar < nVar; jVar++) {
      val_Jacobian_i[iVar][jVar] = 0.0;
      val_Jacobian_j[iVar][jVar] = 0.0;
    }
  
	/*--- Normalized normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Calculate distance between nodes (for use in implicit calculations only) ---*/
	dist_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
	dist_ij = sqrt(dist_ij);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if ( iSpecies < nDiatomics ) {
      loc = (nDim+3)*iSpecies;
      nVar_species = nDim+3;
    }	else {
      loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
      nVar_species = nDim+2;
    }
    
    /*--- Set the primitive variables (T_tr,    vel,  T_vib,      p,    rho,      h,      c) and compute the mean ---*/
    /*--- Indices:                    (   0, iDim+1, nDim+1, nDim+2, nDim+3, nDim+4, nDim+5)                      ---*/
    for (iVar = 0; iVar < nDim+6; iVar++) {
      PrimVar_i[iVar] = Varray_i[iSpecies][iVar];
      PrimVar_j[iVar] = Varray_j[iSpecies][iVar];
      Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar] + PrimVar_j[iVar]);
    }
    
    /*--- Calculate mean values of the transport coefficients ---*/
    Mean_Laminar_Viscosity[iSpecies] = 0.5*(Laminar_Viscosity_MultipleSpecies_i[iSpecies] + Laminar_Viscosity_MultipleSpecies_j[iSpecies]);
    Mean_Eddy_Viscosity[iSpecies] = 0.5*(Eddy_Viscosity_MultipleSpecies_i[iSpecies] + Eddy_Viscosity_MultipleSpecies_j[iSpecies]);
    Mean_Thermal_Conductivity[iSpecies] = 0.5*(Thermal_Conductivity_MultipleSpecies_i[iSpecies] + Thermal_Conductivity_MultipleSpecies_j[iSpecies]);
    Mean_Thermal_Conductivity_vib[iSpecies] = 0.5*(Thermal_Conductivity_vib_MultipleSpecies_i[iSpecies] + Thermal_Conductivity_vib_MultipleSpecies_j[iSpecies]);
    
    /*--- Mean gradient approximation ---*/
    for (iVar = 0; iVar < nDim+3; iVar++)
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i_array[iSpecies][iVar][iDim] + PrimVar_Grad_j_array[iSpecies][iVar][iDim]);
      }
    
    /*--- Get projected viscous flux tensor ---*/
    GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_Thermal_Conductivity, Mean_Thermal_Conductivity_vib, iSpecies);
    
    /*--- Update viscous residual with the species projected flux ---*/
    for (iVar = 0; iVar < nVar_species; iVar++)
      val_residual[loc+iVar] = Proj_Flux_Tensor[iVar];
    
    /*--- Implicit part ---*/
    if (implicit) {
      
      if (config->GetKind_GasModel() == ARGON)
        GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitaryNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j, iSpecies);
      
      else {
        
        GetViscousProjJacsDiatomics(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_Thermal_Conductivity, Mean_Thermal_Conductivity_vib, dist_ij, UnitaryNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j, iSpecies);
      }
    }
  }
}

CAvgGradCorrected_Plasma::CAvgGradCorrected_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies,unsigned short val_nDiatomics, unsigned short val_nMonatomics,  CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Allocate some useful vectors ---*/
	Edge_Vector = new double [nDim];
	Proj_Mean_GradPrimVar_Edge = new double [nDim+3];
  
	PrimVar_i = new double [nDim+6];
	PrimVar_j = new double [nDim+6];
	Mean_Laminar_Viscosity = new double [nSpecies];
	Mean_Eddy_Viscosity = new double [nSpecies];
	Mean_PrimVar = new double [nDim+6];
	Mean_Density = new double [nSpecies];
	Mean_GradPrimVar = new double* [nDim+3];
	for (iVar = 0; iVar < nDim+3; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGradCorrected_Plasma::~CAvgGradCorrected_Plasma(void) {
  
	delete [] PrimVar_i;
	delete [] PrimVar_j;
	delete [] Mean_PrimVar;
	for (iVar = 0; iVar < nDim+6; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
	delete [] Proj_Mean_GradPrimVar_Edge;
	delete [] Edge_Vector;
}

void CAvgGradCorrected_Plasma::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	unsigned short loc, iVar, nVar_species, iSpecies;
  
	/*--- Normalized normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
	}
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) {
			loc = (nDim+3)*iSpecies;
			nVar_species = nDim+3;
		}	else {
			loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			nVar_species = nDim+2;
		}
		/*--- Set the primitive variables (T_tr,    vel,  T_vib,      p,    rho,      h,      c) and compute the mean ---*/
		/*--- Indices:                    (   0, iDim+1, nDim+1, nDim+2, nDim+3, nDim+4, nDim+5)                      ---*/
		for (iVar = 0; iVar < nDim+6; iVar++) {
			PrimVar_i[iVar] = Varray_i[iSpecies][iVar];
			PrimVar_j[iVar] = Varray_j[iSpecies][iVar];
			Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar] + PrimVar_j[iVar]);
		}
		Mean_Laminar_Viscosity[iSpecies] = 0.5*(Laminar_Viscosity_MultipleSpecies_i[iSpecies] + Laminar_Viscosity_MultipleSpecies_j[iSpecies]);
		Mean_Eddy_Viscosity[iSpecies] = 0.5*(Eddy_Viscosity_MultipleSpecies_i[iSpecies] + Eddy_Viscosity_MultipleSpecies_j[iSpecies]);
    
		/*--- Projection of the mean gradient in the direction of the edge ---*/
		for (iVar = 0; iVar < nDim+3; iVar++) {
			Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i_array[iSpecies][iVar][iDim] + PrimVar_Grad_j_array[iSpecies][iVar][iDim]);
				Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
			}
			for (iDim = 0; iDim < nDim; iDim++) {
				Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                         (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
			}
		}
    
		/*--- Get projected viscous flux tensor ---*/
		GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, iSpecies);
    
		/*--- Update viscous residual with the species projected flux ---*/
		for (iVar = 0; iVar < nVar_species; iVar++)
			val_residual[loc+iVar] = Proj_Flux_Tensor[iVar];
    
		/*--- Implicit part ---*/
		if (implicit) {
			GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, sqrt(dist_ij_2), UnitaryNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j, iSpecies);
		}
	}
}

CSourcePieceWise_Plasma::CSourcePieceWise_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics,
                                                 unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar, val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	unsigned short iReaction;
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	U_id = new double [nVar];
	Temp_tr_id = new double[nSpecies];
  
	nSpecies 		= config->GetnSpecies();
	nReactions 		= config->GetnReactions();
	//	GammaMonatomic 	= config->GetGammaMonatomic();
	//	GammaDiatomic 	= config->GetGammaDiatomic();
	nMonatomics 	= config->GetnMonatomics();
	nDiatomics 		= config->GetnDiatomics();
  
	Molar_Mass 		= new double[nSpecies];
	ChargeNumber 	= new double[nSpecies];
	w_dot 			= new double[nSpecies];					w_dotd = new double[nSpecies];
	Q_tv 			= new double[nSpecies];
	Q_elastic 		= new double [nSpecies];			Q_elasticd	  = new double [nSpecies];
	Enthalpy_formation = new double [nSpecies];
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Molar_Mass[iSpecies] = config->GetMolar_Mass(iSpecies);
		ChargeNumber[iSpecies] = config->GetParticle_ChargeNumber(iSpecies);
		Enthalpy_formation[iSpecies] = config->GetEnthalpy_Formation(iSpecies);
	}
  
	ExtentOfReaction = new double[nReactions];
	ReactionRateFwd  = new double[nReactions];	ReactionRateFwdd = new double[nReactions];
	ReactionRateBkw  = new double[nReactions];	ReactionRateBkwd = new double[nReactions];
	fwdRxn			 = new double[nReactions];			fwdRxnd			= new double[nReactions];
	bkwRxn			 = new double[nReactions];			bkwRxnd		= new double[nReactions];
	Keq				 = new double[nReactions];				Keqd = new double[nReactions];
	Cf				 = new double[nReactions];
	eta				 = new double[nReactions];
	theta			 = new double[nReactions];
	RxnReactants	 = new int*[nReactions];
	RxnProducts 	 = new int*[nReactions];
	EqRxnConstants	 = new double*[nReactions];
	T_rxnf = new double[nReactions];						T_rxnfd = new double[nReactions];
	T_rxnb = new double[nReactions];						T_rxnbd = new double[nReactions];
  
	SpeciesPressure_id = new double[nSpecies];
  
	RxnConstantTable = new double*[6];
	for (iVar = 0; iVar < 6; iVar++)
		RxnConstantTable[iVar] = new double[5];
  
	for (iReaction = 0; iReaction < nReactions; iReaction++) {
		RxnReactants[iReaction] = new int [nSpecies];
		RxnProducts[iReaction]  = new int [nSpecies];
		EqRxnConstants[iReaction] = new double[5];
	}
  
	if (config->GetKind_GasModel() == N2   || config->GetKind_GasModel() == O2 ||
      config->GetKind_GasModel() == AIR5 || config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == ARGON_SID) {
		for (iReaction = 0; iReaction < nReactions; iReaction++) {
			Cf[iReaction] = config->GetArrheniusCoeff(iReaction);
			eta[iReaction] = config->GetArrheniusEta(iReaction);
			theta[iReaction] = config->GetArrheniusTheta(iReaction);
		}
	}
  
	Omega00 = config->GetCollisionIntegral00();
	Omega11 = config->GetCollisionIntegral11();
  
	Molecular_Mass = new double [nSpecies];
	Molecular_Diameter = new double [nSpecies];
	Energy_vib = 0.0; Energy_el = 0.0;
	P = new double*[nSpecies];								Pd = new double*[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		P[iSpecies] = new double [nDim];				Pd[iSpecies] = new double [nDim];
	}
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Molecular_Mass[iSpecies] = config->GetMolar_Mass(iSpecies)/AVOGAD_CONSTANT;
		Molecular_Diameter[iSpecies] = config->GetMolecular_Diameter(iSpecies);
	}
  
	Reactions = config->GetReaction_Map();
	nReactions = config->GetnReactions();
  
	CharVibTemp = new double [nSpecies];
	CharVibTempd = new double [nSpecies];
  
	Residual_New = new double[nVar];
	Residual_Baseline = new double[nVar];
	U_Baseline = new double[nVar];
  
  
	unsigned short iDim, iVar, iSpecies;
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Kb = 1.38E-23;
  
	SourceVector = new double[nVar];
	SourceJacobian = new double*[nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar ++) {
		SourceVector[iVar] = 10.0;
		SourceJacobian[iVar] = new double [nVar];
		for (unsigned short jVar = 0; jVar < nVar; jVar ++)
			SourceJacobian[iVar][jVar] = 0.0;
	}
  
  
	MagneticField = new double [3];
	Mag_Force 	  = new double* [nSpecies];
	VcrossB = new double* [nSpecies];
	velocity = new double *[nSpecies];
  
	/* For the more accurate magnetic field model */
	VioncrossB = new double [nDim];
	Current_Density = new double [3];
	JcrossB = new double [3];
	Electric_Conductivity = config->GetElec_Conductivity();
	vector_r = new double [nDim];
	dpcenter = new double [nDim];
  
	for (iVar = 0; iVar < 3; iVar ++) {
		Current_Density[iVar] = 0.0;
		JcrossB[iVar] = 0.0;
	}
  
	/* Till here For the more accurate magnetic field model */
  
	MagneticDipole = new double [3];
	ElectricField = new double [nDim];
  
	for (iVar = 0; iVar < 3; iVar ++) {
		MagneticField[iVar] = 0.0;
		MagneticDipole[iVar] = 0.0;
	}
	Stagnation_B = config->GetStagnation_B();
	dpcenter[0] = config->GetDipoleDist();
	dpcenter[1] = 0.0;
	dpcenter[2] = 0.0;
  
	MagneticDipole[0] = Stagnation_B*2*PI_NUMBER*pow((0.01825+dpcenter[0]),3)/MAGNETIC_CONSTANT;
	MagneticDipole[1] = 0.0;
	MagneticDipole[2] = 0.0;
  
	for (iDim = 0; iDim < nDim; iDim++)
		ElectricField[iDim] = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		velocity[iSpecies] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim ++) velocity[iSpecies][iDim]  = 0.0;
		VcrossB[iSpecies] = new double [3];
		for (iVar = 0; iVar < 3; iVar ++) VcrossB[iSpecies][iVar] = 0.0;
		Mag_Force[iSpecies] = new double [3];
		for (iVar = 0; iVar < 3; iVar ++) Mag_Force[iSpecies][iVar] = 0.0;
    
	}
  
	AvgNum = AVOGAD_CONSTANT;
	M1Avg = config->GetMolar_Mass(0);
	M2Avg = config->GetMolar_Mass(1);
	M3Avg = config->GetMolar_Mass(2);
	M1M2M3Avg3 = M1Avg * M2Avg* M3Avg;
  
	M1 = M1Avg/AVOGAD_CONSTANT;
	M3 = M3Avg/AVOGAD_CONSTANT;
	M2 = M1-M3;
  
	ec = ELECTRON_CHARGE;    eps0 = FREE_PERMITTIVITY;
	Te = 10000;
	Rc = 8314.462175;
	r12 = 4E-10;     r13 = 2E-10;     r23 = ec*ec/(32.0*eps0*Kb*Te);
	sigma12 = PI_NUMBER*r12*r12;     sigma13 = PI_NUMBER*r13*r13;     sigma23 = PI_NUMBER*r23*r23;
  
	/*--- Specific heats at constant volume for each species ---*/
	Cv1 = 1.5*Rc/(M1Avg);                                                       // Specific heat of Argon gas
	Cv2 = 1.5*Rc/(M2Avg);                                                       // Specific heat of Argon positive ion
	Cv3 = 1.5*Rc/(M3Avg);                                                       // Specific heat of electron
  
	gam1 = 1 + Rc/(M1Avg*Cv1);
	gam2 = 1 + Rc/(M2Avg*Cv2);
	gam3 = 1 + Rc/(M3Avg*Cv3);
	tol = 1E-15;
	Tstart = 1000.0;
  
  
}

CSourcePieceWise_Plasma::~CSourcePieceWise_Plasma(void) {
	delete [] Molar_Mass;
	delete [] Molecular_Mass;
	delete [] Molecular_Diameter;
	delete [] w_dot;						delete [] w_dotd;
	delete [] Q_tv;
	delete [] Q_elastic;				delete [] Q_elasticd;
	delete [] ChargeNumber;
	delete [] Enthalpy_formation;
  
	delete [] ExtentOfReaction;
	delete [] ReactionRateFwd;	delete [] ReactionRateFwdd;
	delete [] ReactionRateBkw;	delete [] ReactionRateBkwd;
	delete [] fwdRxn;						delete [] fwdRxnd;
	delete [] bkwRxn;						delete [] bkwRxnd;
	delete [] Keq;							delete [] Keqd;
	delete [] T_rxnf;						delete [] T_rxnfd;
	delete [] T_rxnb;						delete [] T_rxnbd;
	delete [] Cf;
	delete [] eta;
	delete [] theta;
	delete [] CharVibTemp;
	delete [] CharVibTempd;
  
	for (iVar = 0; iVar < 6; iVar++)
		delete [] RxnConstantTable[iVar];
	delete [] RxnConstantTable;
  
	unsigned short iReaction;
	for (iReaction = 0; iReaction < nReactions; iReaction++) {
		delete [] RxnReactants[iReaction];
		delete [] RxnProducts[iReaction];
		delete [] EqRxnConstants[iReaction];
	}
	delete [] RxnReactants;	delete [] RxnProducts;	delete [] EqRxnConstants;
  
	unsigned short iSpecies;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		delete [] P[iSpecies];
		delete [] Pd[iSpecies];
	}
	delete [] P;
	delete [] Pd;
  
	delete [] Temp_tr_id;
	delete [] SpeciesPressure_id;
  
	delete [] Residual_New;	delete [] Residual_Baseline;	delete [] U_Baseline;
  
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] SourceJacobian[iVar];
	delete [] SourceVector;	delete[] SourceJacobian;
  
	//delete [] ElectricField;
	delete[] MagneticField;
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		delete [] velocity[iSpecies];
		delete [] VcrossB[iSpecies];
	}
	delete [] velocity;		delete [] VcrossB;
	delete [] VioncrossB;	delete [] Current_Density; delete [] JcrossB;
	delete [] MagneticDipole;
	delete [] dpcenter;
	delete [] vector_r;
  
}

void CSourcePieceWise_Plasma::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
  
	tol = 1E-60;
	double zero = 1E-15;
	double dens_ratio;
  
  //	cout << " Electric Field = " << ElectricField[0] << " , " << ElectricField[1] << endl;
  //    cin.get();
  
	/* brief: defining conservative variables, residual and jacobian for point i first: */
	if (nDim ==2 ) {
		r1 = U_i[0];	r2 = U_i[4];	r3 = U_i[8];	// density of each species
		m1 = U_i[1];	m2 = U_i[5];	m3 = U_i[9];	// density*horizontal velocity of each species
		n1 = U_i[2];	n2 = U_i[6];	n3 = U_i[10];	// density*vertical velocity of each species
		e1 = U_i[3];	e2 = U_i[7];	e3 = U_i[11];	// energy per unit volume of each species
    
	}
	l1 = 0;	l2 = 0;	l3 = 0;	// density*vertical velocity of each species
  
	if (nDim ==3 ) {
		r1 = U_i[0];	r2 = U_i[5];	r3 = U_i[10];	// density of each species
		m1 = U_i[1];	m2 = U_i[6];	m3 = U_i[11];	// density*horizontal velocity of each species
		n1 = U_i[2];	n2 = U_i[7];	n3 = U_i[12];	// density*vertical velocity of each species
		l1 = U_i[3];	l2 = U_i[8];	l3 = U_i[13];	// density*vertical velocity of each species
		e1 = U_i[4];	e2 = U_i[9];	e3 = U_i[14];	// energy per unit volume of each species
	}
  
	P1 = (gam1-1)*(e1 - 0.5*(m1*m1 + n1*n1 + l1*l1)/r1 - (Enthalpy_formation[0] + Energy_vib + Energy_el)*r1 );	// Partial Pressure of species 1
	P2 = (gam2-1)*(e2 - 0.5*(m2*m2 + n2*n2 + l2*l2)/r2 - (Enthalpy_formation[1] + Energy_vib + Energy_el)*r2 );	// Partial Pressure of species 1
	P3 = (gam3-1)*(e3 - 0.5*(m3*m3 + n3*n3 + l3*l3)/r3 - (Enthalpy_formation[2] + Energy_vib + Energy_el)*r3 );	// Partial Pressure of species 1
  
	//	P1 = (gam1-1)*(e1 - 0.5*(m1*m1 + n1*n1 + l1*l1)/r1);	// Partial Pressure of species 2
	//	P2 = (gam2-1)*(e2 - 0.5*(m2*m2 + n2*n2 + l2*l2)/r2);	// Partial Pressure of species 2
	//	P3 = (gam3-1)*(e3 - 0.5*(m3*m3 + n3*n3 + l3*l3)/r3);	// Partial Pressure of species 3
  
	T1 = P1/(Rc/(M1*AvgNum)*r1);	// Temperature of species 1
	T2 = P2/(Rc/(M2*AvgNum)*r2);	// Temperature of species 2
	T3 = P3/(Rc/(M3*AvgNum)*r3);	// Temperature of species 3
  
	r1 = max(r1, zero);	r2 = max(r2, zero);	r3 = max(r3, zero);
	P1 = max(P1, zero);	P2 = max(P2, zero);	P3 = max(P3, zero);
	T1 = max(T1, zero);	T2 = max(T2, zero);	T3 = max(T3, zero);
  
	/*Brief: Partial derivative of species temperature with the conservative variables */
  
	dT1_dr1 = (-Cv1*T1 + .5*(m1*m1 + n1*n1 + l1*l1)/(r1*r1))/(r1*Cv1);
	dT2_dr2 = (-Cv2*T2 + .5*(m2*m2 + n2*n2 + l2*l2)/(r2*r2))/(r2*Cv2);
	dT3_dr3 = (-Cv3*T3 + .5*(m3*m3 + n3*n3 + l3*l3)/(r3*r3))/(r3*Cv3);
  
	dT1_dm1 = - m1/r1 / (r1*Cv1) ;
	dT2_dm2 = - m2/r2 / (r2*Cv2);
	dT3_dm3 = - m3/r3 / (r3*Cv3);
  
	dT1_dn1 = - n1/r1 / (r1*Cv1) ;
	dT2_dn2 = - n2/r2 / (r2*Cv2);
	dT3_dn3 = - n3/r3 / (r3*Cv3);
  
	dT1_dl1 = - l1/r1 / (r1*Cv1) ;
	dT2_dl2 = - l2/r2 / (r2*Cv2);
	dT3_dl3 = - l3/r3 / (r3*Cv3);
  
	dT1_de1 = 	1.0/(Cv1*r1);
	dT2_de2 = 	1.0/(Cv2*r2);
	dT3_de3 = 	1.0/(Cv3*r3);
  
	Tstart = 1000.0;
	T3 = max(T3,Tstart);
	T1 = max(T1,Tstart);
	T2 = max(T2,Tstart);
  
	/* Relaxation factors f12, f13, f23 */
  
	double max_limit = 0.0;
  
	//	if ( T1 > 3000.0)
	max_limit = 3E-15 ;//* sqrt(1E-8 * 1000/(T3* r3));
  
	f23 = r3/M3*sqrt(8*Kb*T3/(PI_NUMBER*M3))*max(1.95E-10/(T3*T3)*log(1.53E14*T3*T3*T3*M3/r3), max_limit);
  
	if (T3 < 10000) f13 = r1/M1*sqrt(8*Kb*T3/(PI_NUMBER*M3)) * (0.39 - 0.551E-4*T3 + 0.595E-8*T3*T3)*1E-20;
	else f13 = r1/M1*sqrt(8*Kb*T3/(PI_NUMBER*M3))* (-0.35 + 0.775E-4*T3)*1E-20;
  
	f12 = pow((sigma12/sigma13),1.5)*f13;	/* f12 was estimated to be equal to f13 */
  
	/* Heat Transfer Terms, following Hoffert and Lien here */
	QT1 = 2*r2*Cv2*f12*M2/M1*(T2-T1) + 2*r3*Cv3*f13*M3/M1*(T3-T1);	// Heat transfered to argon atoms
	QT2 = 2*r2*Cv2*f12*M2/M1*(T1-T2) + 2*r3*Cv3*f23*M3/M2*(T3-T2);	// Heat transfered to argon positive ions
	QT3 = 2*r3*Cv3*f23*M3/M2*(T2-T3) + 2*r3*Cv3*f13*M3/M1*(T1-T3);	// Heat transfered to electrons
  
	/* Derivative of the heat transfer terms with respect to species temperature */
	dQT1_dT1 = -2*r2*Cv2*f12*M2/M1 - 2*r3*Cv3*f13*M3/M1;	//derivative of QT1 wrt T1
  
	dQT2_dT2 =-2*r2*Cv2*f12*M2/M1 - 2*r3*Cv3*f23*M3/M2;	//derivative of QT2 wrt T2
  
	dQT3_dT3 = - 2*r3*Cv3*f23*M3/M2 - 2*r3*Cv3*f13*M3/M1;	//derivative of QT3 wrt T3
  
	T1 = max(T1,Tstart);
	T2 = T1;
  
	/*Brief: Partial derivative of collision frequencies wrt temperature variables */
  
	/*
   nu12 = frequency of momentum transfer from species 2 to species 1
   nu21 = frequency of momentum transfer from species 1 to species 2
   nu13 = frequency of momentum transfer from species 3 to species 1
   nu31 = frequency of momentum transfer from species 1 to species 3
   nu23 = frequency of momentum transfer from species 3 to species 2
   nu32 = frequency of momentum transfer from species 2 to species 3
   
	 */
  
	C12   = sqrt(8.0*Kb/PI_NUMBER* (T1/M1 + T2/M2));	// collision velocity of species 1,2
	C13   = sqrt(8.0*Kb/PI_NUMBER* (T1/M1 + T3/M3));	// collision velocity of species 1,3
	C23   = sqrt(8.0*Kb/PI_NUMBER* (T2/M2 + T3/M3));	// collision velocity of species 2,3
  
	r23 = ec*ec/(32.0*eps0*Kb*T3);sigma23 = PI_NUMBER*r23*r23;
	nu12 = r2/(M2+M1) * sigma12 * C12; 	 nu21 = r1/(M2+M1) * sigma12 * C12;
	nu13 = r3/(M3+M1) * sigma13 * C13; 	 nu31 = r1/(M3+M1) * sigma13 * C13;
	nu23 = r3/(M2+M3) * sigma23 * C23; 	 nu32 = r2/(M2+M3) * sigma23 * C23;
  
	dv12_dT1 = 0; 	dv12_dT2 = 0;	dv21_dT1 = 0;	dv21_dT2 = 0;
	dv13_dT1 = 0;	dv13_dT3 = 0;	dv31_dT1 = 0;	dv31_dT3 = 0;
	dv23_dT2 = 0;	dv23_dT3 = 0;	dv32_dT2 = 0;	dv32_dT3 = 0;
  
	dC12_dT1 = 1/(2*M1) * sqrt(8*Kb/PI_NUMBER) * 1/ (sqrt(T1/M1 + T2/M2)); 	dC12_dT2 = 1/(2*M2) * sqrt(8*Kb/PI_NUMBER) * 1/ (sqrt(T1/M1 + T2/M2));
  
	dC13_dT1 = 1/(2*M1) * sqrt(8*Kb/PI_NUMBER) * 1/ (sqrt(T1/M1 + T3/M3)); 	dC13_dT3 = 1/(2*M3) * sqrt(8*Kb/PI_NUMBER) * 1/ (sqrt(T1/M1 + T3/M3));
  
	dC23_dT2 = 1/(2*M2) * sqrt(8*Kb/PI_NUMBER) * 1/ (sqrt(T2/M2 + T3/M3)); 	dC23_dT3 = 1/(2*M3) * sqrt(8*Kb/PI_NUMBER) * 1/ (sqrt(T2/M2 + T3/M3));
  
	dv12_dT1 = r2/(M2+M1) * sigma12 * dC12_dT1; 	dv12_dT2 = r2/(M2+M1) * sigma12 * dC12_dT2;
	dv21_dT1 = r1/(M2+M1) * sigma12 * dC12_dT1; 	dv21_dT2 = r1/(M2+M1) * sigma12 * dC12_dT2;
  
	dv13_dT1 = r3/(M3+M1) * sigma13 * dC13_dT1; 	dv13_dT3 = r3/(M3+M1) * sigma13 * dC13_dT3;
	dv31_dT1 = r1/(M3+M1) * sigma13 * dC13_dT1; 	dv31_dT3 = r1/(M3+M1) * sigma13 * dC13_dT3;
  
	dv23_dT2 = r3/(M2+M3) * sigma23 * dC23_dT2; 	dv23_dT3 = r3/(M2+M3) * sigma23 * dC23_dT3;
	dv32_dT2 = r2/(M2+M3) * sigma23 * dC23_dT2; 	dv32_dT3 = r2/(M2+M3) * sigma23 * dC23_dT3;
  
	/* Mass conservation terms, species formation and destruction formulae  */
	/*
   kf1 = rate of formation of species 1
   kf2 = rate of formation of species 2
   kf3 = rate of formation of species 3
   
   kb1 = rate of backward reaction of species 1
   kb2 = rate of backward reaction of species 2
   kb3 = rate of backward reaction of species 3
   
	 */
  
	/* Table of reaction constants  */
  
	/* ****************************************************************/
	/***/ C1     = 10.12;	 C2     = 10.12;    C3	   = 22.60E4;   /***/
	/***/ Ck1    = 2.9E22;   Ck2    = 2.9E22;   Ck3    = 2.9E22;   /***/
	/***/ eta1   = 1.5;	 eta2   = 1.5;	    eta3   = 1.5;	   /***/
	/***/ zeta1  = 1.5;	 zeta2  = 1.5;	    zeta3  = 1.5;	   /***/
	/***/ theta1 = 135300.0; theta2 = 135300.0; theta3 = 135300.0; /***/
	/***/ phi1   = 183100.0; phi2   = 183100.0; phi3   = 183100.0; /***/
	/******************************************************************/
  
	//T3 = sqrt(T3*T1);
	T1 = max(T1,Tstart);
	T2 = T1;
  
	kf1 = C1*pow(T1,eta1)*(theta1/T1 + 2.0)*exp(-theta1/T1);
	kf2 = C2*pow(T2,eta2)*(theta2/T2 + 2.0)*exp(-theta2/T2);
	kf3 = C3*pow(T3,eta3)*(theta3/T3 + 2.0)*exp(-theta3/T3);
  
	ke1 = Ck1*pow(T1,zeta1)*exp(-phi1/T1);
	ke2 = Ck2*pow(T2,zeta2)*exp(-phi2/T2);
	ke3 = Ck3*pow(T3,zeta3)*exp(-phi3/T3);
  
	kb1 = kf1/ke1; 	kb2 = kf2/ke2; 	kb3 = kf3/ke3;
  
	if (T1 > Tstart) kb1 = kf1/ke1;
	else kb1 = 0.0;
	if ( T2 > Tstart)  kb2 = kf2/ke2;
	else kb2 = 0.0;
	if ( T3 > Tstart)  kb3 = kf3/ke3;
	else kb3 = 0.0;
  
	R =  -kf1*r1*r1/(M1*M1*AvgNum*AvgNum) + kb1*r2*r3*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum) +
  -kf2*r1*r2/(M1*M2*AvgNum*AvgNum) + kb2*r2*r3*r2/(M2*M3*M2*AvgNum*AvgNum*AvgNum) +
  -kf3*r1*r3/(M1*M3*AvgNum*AvgNum) + kb3*r2*r3*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum);
  
	dkf1_dT1 = C1*eta1*pow(T1,eta1-1.0)*(theta1/T1 + 2.0)*exp(-theta1/T1) +
  C1*pow(T1,eta1)*(-theta1/(T1*T1))*exp(-theta1/T1)          +
  C1*pow(T1,eta1)*(theta1/T1 + 2.0)*exp(-theta1/T1)*(theta1/(T1*T1));
  
	dkf2_dT2 = C2*eta2*pow(T2,eta2-1.0)*(theta2/T2 + 2.0)*exp(-theta2/T2) +
  C2*pow(T2,eta2)*(-theta2/(T2*T2))*exp(-theta2/T2)          +
  C2*pow(T2,eta2)*(theta2/T2 + 2.0)*exp(-theta2/T2)*(theta2/(T2*T2));
  
	dkf3_dT3 = C3*eta3*pow(T3,eta3-1.0)*(theta3/T3 + 2.0)*exp(-theta3/T3) +
  C3*pow(T3,eta3)*(-theta3/(T3*T3))*exp(-theta3/T3)          +
  C3*pow(T3,eta3)*(theta3/T3 + 2.0)*exp(-theta3/T3)*(theta3/(T3*T3));
  
	dke1_dT1 = Ck1*zeta1*pow(T1,zeta1-1.0)*exp(-phi1/T1) + Ck1*pow(T1,zeta1)*exp(-phi1/T1)*(phi1/(T1*T1));
	dke2_dT2 = Ck2*zeta2*pow(T2,zeta2-1.0)*exp(-phi2/T2) + Ck2*pow(T2,zeta2)*exp(-phi2/T2)*(phi2/(T2*T2));
	dke3_dT3 = Ck3*zeta3*pow(T3,zeta3-1.0)*exp(-phi3/T3) + Ck3*pow(T3,zeta3)*exp(-phi3/T3)*(phi3/(T3*T3));
  
	dkb1_dT1 = (ke1*dkf1_dT1 - kf1*dke1_dT1)/(ke1*ke1);
	dkb2_dT2 = (ke2*dkf2_dT2 - kf2*dke2_dT2)/(ke2*ke2);
	dkb3_dT3 = (ke3*dkf3_dT3 - kf3*dke3_dT3)/(ke3*ke3);
  
	/* derivative of expression for rate of reaction R wrt to first four conservative variables */
	dR_dr1 = ( - 2*kf1*r1/(M1*M1*AvgNum*AvgNum)	 + kb1*r2*r3/(M2*M3*M1*AvgNum*AvgNum*AvgNum)
            - kf2*r2/(M1*M2*AvgNum*AvgNum)	 + 0.0
            - kf3*r3/(M1*M3*AvgNum*AvgNum)	 + 0.0);
	dR_dm1 = 0.0;
	dR_dn1 = 0.0;
	dR_dl1 = 0.0;
	dR_de1 =  ( -r1*r1/(M1*M1*AvgNum*AvgNum)*dkf1_dT1*dT1_de1   + r2*r3*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)*dkb1_dT1*dT1_de1);
  
	/* derivative of expression for rate of reaction R wrt to 5th,6th,7th and 8th conservative variables */
	dR_dr2 = ( + 0.0	   + kb1*r3*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)
            - kf2*r1/(M1*M2*AvgNum*AvgNum)	   + 2*kb2*r2*r3/(M2*M3*M2*AvgNum*AvgNum*AvgNum)
            + 0.0	   + kb3*r3*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum));
	dR_dm2 = 0.0;
	dR_dn2 = 0.0;
	dR_dl2 = 0.0;
	dR_de2 = ( - r1*r2/(M1*M2*AvgNum*AvgNum)*dkf2_dT2*dT2_de2  + r2*r3*r2/(M2*M3*M2*AvgNum*AvgNum*AvgNum)*dkb2_dT2*dT2_de2);
  
  
	/* derivative of expression for rate of reaction R wrt to last four conservative variables */
	dR_dr3 = ( + 0.0	 + kb1*r2*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)
            + 0.0	 	 + kb2*r2*r2/(M2*M3*M2*AvgNum*AvgNum*AvgNum)
            - kf3*r1/(M1*M3*AvgNum*AvgNum)	 + 2*kb3*r2*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum));
	dR_dr3 = dR_dr2 * M3/M2;
	dR_dm3 = 0.0;
	dR_dn3 = 0.0;
	dR_dl3 = 0.0;
	dR_de3 = ( - r1*r3/(M1*M3*AvgNum*AvgNum)*dkf3_dT3*dT3_de3  + r2*r3*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum)*dkb3_dT3*dT3_de3);
	dR_de3 = dR_de2*M3/M2;
  
	velocity[0][0] = m1/r1; 	velocity[0][1] = n1/r1;
	velocity[1][0] = m2/r2; 	velocity[1][1] = n2/r2;
	velocity[2][0] = m3/r3; 	velocity[2][1] = n3/r3;
	if (nDim ==3)
		velocity[0][2] = l1/r1; 	velocity[1][2] = l2/r2;  velocity[2][2] = l3/r3;
  
	mdotr = 0.0;
	distance = 0.0;
  
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		vector_r[iDim] = Coord_i[iDim] - dpcenter[iDim];
		mdotr +=MagneticDipole[iDim]*vector_r[iDim];
		distance += vector_r[iDim]*vector_r[iDim];
	}
	distance = sqrt(distance);
	rpower5 = pow(distance,5);
	rcubed  = pow(distance,3);
	for (iDim = 0; iDim < nDim; iDim ++ )
		MagneticField[iDim] = MAGNETIC_CONSTANT/(4.0*PI_NUMBER) * (3.0 * vector_r[iDim]* mdotr/rpower5 - MagneticDipole[iDim]/rcubed);
  
  
	//
	//
	//	mdotr = 0.0;
	//	distance = 0.0;
	//	for (iDim = 0; iDim < nDim; iDim ++ ) {
	//		mdotr +=MagneticDipole[iDim]*Coord_i[iDim];
	//		distance += Coord_i[iDim]*Coord_i[iDim];
	//	}
	//
	//	distance = sqrt(distance);
	//	rpower5 = pow(distance,5);
	//	rcubed  = pow(distance,3);
	//
	//	for (iDim = 0; iDim < nDim; iDim ++ )
	//		MagneticField[iDim] = MAGNETIC_CONSTANT/(4.0*PI_NUMBER) * (3.0 * Coord_i[iDim]* mdotr/rpower5 - MagneticDipole[iDim]/rcubed);
  
  
  
  
  
	double Bx = MagneticField[0];
	double By = MagneticField[1];
	double Bz = MagneticField[2];
	double s0 = Electric_Conductivity;
  
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		if (nDim ==2) {
			VcrossB[iSpecies][0] = 0.0;	VcrossB[iSpecies][1] = 0.0;
			Mag_Force[iSpecies][0] = 0.0; Mag_Force[iSpecies][1] = 0.0;
		}
		else {
			VcrossB[iSpecies][0] = velocity[iSpecies][1]*MagneticField[2] - velocity[iSpecies][2]*MagneticField[1];
			VcrossB[iSpecies][1] = velocity[iSpecies][2]*MagneticField[0] - velocity[iSpecies][0]*MagneticField[2];
			VcrossB[iSpecies][2] = velocity[iSpecies][0]*MagneticField[1] - velocity[iSpecies][1]*MagneticField[0];
			Mag_Force[iSpecies][0] = 0.0; Mag_Force[iSpecies][1] = 0.0;	Mag_Force[iSpecies][2] = 0.0;
      
		}
	}
  
	for (iDim = 0; iDim < nDim; iDim ++ )
		Current_Density[iDim] = Electric_Conductivity*(ElectricField[iDim] + VcrossB[0][iDim]);
  
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		JcrossB[0] = Current_Density[1]*MagneticField[2] - Current_Density[2]*MagneticField[1];
		JcrossB[1] = Current_Density[2]*MagneticField[0] - Current_Density[0]*MagneticField[2];
		JcrossB[2] = Current_Density[0]*MagneticField[1] - Current_Density[1]*MagneticField[0];
	}
  
	double JdotJoversgima = 0.0;
  
	for (iDim = 0; iDim < nDim; iDim ++ )
		JdotJoversgima += (Current_Density[iDim]*Current_Density[iDim]);
  
	JdotJoversgima = JdotJoversgima/s0 ;
	/*
	 * Magnetic Forces using the Jcross B formulation result in the following expression:
	 * JcrossB[0] = Fx = ihat  [ Bz * ( Uz*Bx - Ux*Bz ) - By * ( Ux*By - Uy*Bx ) ] * ElectricalConductivity;
	 * JcrossB[1] = Fy = jhat  [ Bx * ( Ux*By - Uy*Bx ) - Bz * ( Uy*Bz - Uz*By ) ] * ElectricalConductivity;
	 * JcrossB[2] = Fz = khat  [ By * ( Uy*Bz - Uz*By ) - Bx * ( Uz*Bx - Ux*Bz ) ] * ElectricalConductivity;
	 *
	 * U dot JcrossB[0] =  [ Bz * ( Ux*Uz*Bx - Ux*Ux*Bz ) - By * ( Ux*Ux*By - Ux*Uy*Bx ) ] +
	 *					   [ Bx * ( Uy*Ux*By - Uy*Uy*Bx ) - Bz * ( Uy*Uy*Bz - Uy*Uz*By ) ] +
	 *	                   [ By * ( Uz*Uy*Bz - Uz*Uz*By ) - Bx * ( Uz*Uz*Bx - Uz*Ux*Bz ) ] * ElectricalConductivity
	 */
  
	double delJcrossBdelm = s0/r1/r1 * ( Bx * (l1*Bz +   n1*By) + By * (n1*Bx - 2*m1*By) + Bz * (l1*Bx - 2*m1*Bz) );
	double delJcrossBdeln = s0/r1/r1 * ( Bx * (m1*By - 2*n1*Bx) + By * (m1*Bx +   l1*Bz) + Bz * (l1*By - 2*n1*Bz) );
	double delJcrossBdell = s0/r1/r1 * ( Bx * (m1*Bz - 2*l1*Bx) + By * (n1*Bz - 2*l1*By) + Bz * (m1*Bx +   n1*By) );
  
	double delJdotJdelm   =  1.0/s0*( Current_Density[0] * 0  	  - Bz*s0/r1*Current_Density[1] + By*s0/r1*Current_Density[2] );
	double delJdotJdeln   =  1.0/s0*( Current_Density[0] * Bz*s0/r1 - 0*Current_Density[0] 	    - Bx*s0/r1*Current_Density[2] );
	double delJdotJdell   =  1.0/s0*(-Current_Density[0] * By*s0/r1 + Bx*s0/r1*Current_Density[1] + 0*Current_Density[2] );
  
	/* Derivative of the first term in the source terms wrt all the conservative variables */
	unsigned short loc;
  
	//	JcrossB[0] = 0.0;	JcrossB[1] = 0.0;	JcrossB[2] = 0.0;	JdotJoversgima = 0.0;
	//	delJcrossBdelm = 0.0; delJcrossBdeln = 0.0; delJcrossBdell = 0.0;
	//	delJdotJdelm = 0.0;   delJdotJdeln = 0.0;   delJdotJdell = 0.0;
  
  
	loc = 0;
	dens_ratio = r1/(r1+r2+r3);
	JcrossB[0] = dens_ratio*JcrossB[0]; JcrossB[1] = dens_ratio*JcrossB[1];	JcrossB[2] = dens_ratio*JcrossB[2];
	JdotJoversgima = dens_ratio*JdotJoversgima;
	delJcrossBdelm = dens_ratio*delJcrossBdelm; delJcrossBdeln = dens_ratio*delJcrossBdeln;	delJcrossBdell = dens_ratio*delJcrossBdell;
	delJdotJdelm   = dens_ratio*delJdotJdelm;	delJdotJdeln   = dens_ratio*delJdotJdeln;	delJdotJdell   = dens_ratio*delJdotJdell;
  
  
  
	//	for (iDim = 0; iDim < nDim; iDim ++)
	//		Mag_Force[0][iDim] = dens_ratio*JcrossB[iDim];
  
	SourceVector[loc+0] = M1*AvgNum*R;
  
	SourceJacobian[loc+0][loc+0] =  M1Avg*dR_dr1;
	SourceJacobian[loc+0][loc+1] = M1Avg*dR_dm1;
	SourceJacobian[loc+0][loc+2] = M1Avg*dR_dn1;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] = M1Avg*dR_dl1;
  
	SourceJacobian[loc+0][loc+nDim+1] = M1Avg*dR_de1;
  
	/* Derivative of the second term in the source terms wrt all the conservative variables */
  
	SourceVector[loc+1] = JcrossB[0] + r1*nu12*(m2/r2 - m1/r1) + r1*nu13*(m3/r3 - m1/r1);
  
	SourceJacobian[loc+1][loc+0] = nu12*m2/r2 + nu13*m3/r3 - JcrossB[0]/r1;
	SourceJacobian[loc+1][loc+1] = -nu12 - nu13  - s0 * ( Bz*Bz + By*By)/r1;
	SourceJacobian[loc+1][loc+2] = s0*By*Bx/r1;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] =  s0*Bx*Bz/r1;
	SourceJacobian[loc+1][loc+nDim+1] = JcrossB[0] + r1*dv12_dT1*dT1_de1*(m2/r2 - m1/r1) + r1*dv13_dT1*dT1_de1*(m3/r3 - m1/r1);
  
	/* Derivative of the third term in the source terms wrt all the conservative variables */
  
	SourceVector[loc+2] = JcrossB[1] + r1*nu12*(n2/r2 - n1/r1) + r1*nu13*(n3/r3 - n1/r1);
  
	SourceJacobian[loc+2][loc+0]  = nu12*n2/r2 + nu13*n3/r3 - JcrossB[1]/r1;
	SourceJacobian[loc+2][loc+1]  = s0*Bx*By/r1;
	SourceJacobian[loc+2][loc+2]  = -nu12 - nu13 -s0*(Bx*Bx + Bz*Bz)/r1;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] =  s0*Bz*By/r1;
	SourceJacobian[loc+2][loc+nDim+1] = JcrossB[1] + r1*(n2/r2 -n1/r1)*dv12_dT1*dT1_de1 + r1*(n3/r3 - n1/r1)*dv13_dT1*dT1_de1;
  
  
	if (nDim == 3) {
		SourceVector[loc+3] = JcrossB[2] + r1*nu12*(l2/r2 - l1/r1) + r1*nu13*(l3/r3 - l1/r1);
    
		SourceJacobian[loc+3][loc+0]  = nu12*l2/r2 + nu13*l3/r3 - JcrossB[2]/r1;
		SourceJacobian[loc+3][loc+1]  = s0*Bx*Bz/r1;
		SourceJacobian[loc+3][loc+2]  = s0*By*Bz/r1;
		SourceJacobian[loc+3][loc+3] = -nu12 - nu13 - s0*(By*By + Bx*Bx)/r1;
		SourceJacobian[loc+3][loc+nDim+1] = JcrossB[2] + r1*(l2/r2 -l1/r1)*dv12_dT1*dT1_de1 + r1*(l3/r3 - l1/r1)*dv13_dT1*dT1_de1;
	}
	/* Derivative of the Fourth term in the source terms wrt all the conservative variables */
  
	SourceVector[loc+nDim+1] = JcrossB[0]*m1/r1 + JcrossB[1]*n1/r1 + JcrossB[2]*l1/r1 + JdotJoversgima +
  m1*nu12*(m2/r2 - m1/r1) + m1*nu13*(m3/r3 - m1/r1) +
  n1*nu12*(n2/r2 - n1/r1) + n1*nu13*(n3/r3 - n1/r1) +
  l1*nu12*(l2/r2 - l1/r1) + l1*nu13*(l3/r3 - l1/r1) + QT1;
  
	SourceJacobian[loc+nDim+1][loc+0]  =  -2/r1*(JcrossB[0]*m1/r1 + JcrossB[1]*n1/r1 + JcrossB[2]*l1/r1) - 2*JdotJoversgima/r1 + (m1*m1 + n1*n1 + l1*l1)/(r1*r1) * (nu12 + nu13);
	SourceJacobian[loc+nDim+1][loc+1]  = nu12*(m2/r2 -2*m1/r1) + nu13*(m3/r3 - 2*m1/r1) + delJcrossBdelm + delJdotJdelm + dQT1_dT1*dT1_dm1;
	SourceJacobian[loc+nDim+1][loc+2]  = nu12*(n2/r2 -2*n1/r1) + nu13*(n3/r3 - 2*n1/r1) + delJcrossBdeln + delJdotJdeln + dQT1_dT1*dT1_dn1;
	if (nDim == 3) SourceJacobian[loc+nDim+1][loc+3]  = nu12*(l2/r2 -2*l1/r1) + nu13*(l3/r3 - 2*l1/r1)   + delJcrossBdell + delJdotJdell + dQT1_dT1*dT1_dl1;
  
	SourceJacobian[loc+nDim+1][loc+nDim+1]  = JdotJoversgima/e1 +
  (m1*(m2/r2 - m1/r1) + n1*(n2/r2 - n1/r1) + l1*(l2/r2 - l1/r1))*dv12_dT1*dT1_de1 +
  (m1*(m3/r3 - m1/r1) + n1*(n3/r3 - n1/r1) + l1*(l3/r3 - l1/r1))*dv13_dT1*dT1_de1 + dQT1_dT1*dT1_de1;
  
  
  
  
  
  
  
	/* Derivative of the fifth term in the source terms wrt all the conservative variables */
	loc = (nDim+2);
  
	for (iDim = 0; iDim < nDim; iDim ++ )
		Current_Density[iDim] = Electric_Conductivity*(ElectricField[iDim] + VcrossB[1][iDim]);
  
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		JcrossB[0] = Current_Density[1]*MagneticField[2] - Current_Density[2]*MagneticField[1];
		JcrossB[1] = Current_Density[2]*MagneticField[0] - Current_Density[0]*MagneticField[2];
		JcrossB[2] = Current_Density[0]*MagneticField[1] - Current_Density[1]*MagneticField[0];
	}
	JdotJoversgima = 0.0;
	for (iDim = 0; iDim < nDim; iDim ++ )
		JdotJoversgima += (Current_Density[iDim]*Current_Density[iDim]);
  
	JdotJoversgima = JdotJoversgima/s0;
	delJcrossBdelm = s0/r2/r2 * ( Bx * (l2*Bz +   n2*By) + By * (n2*Bx - 2*m2*By) + Bz * (l2*Bx - 2*m2*Bz) );
	delJcrossBdeln = s0/r2/r2 * ( Bx * (m2*By - 2*n2*Bx) + By * (m2*Bx +   l2*Bz) + Bz * (l2*By - 2*n2*Bz) );
	delJcrossBdell = s0/r2/r2 * ( Bx * (m2*Bz - 2*l2*Bx) + By * (n2*Bz - 2*l2*By) + Bz * (m2*Bx +   n2*By) );
  
	delJdotJdelm   =  1/s0*( Current_Density[0] * 0  	  - Bz*s0/r2*Current_Density[1] + By*s0/r2*Current_Density[2] );
	delJdotJdeln   =  1/s0*( Current_Density[0] * Bz*s0/r2 - 0*Current_Density[0] 	    - Bx*s0/r2*Current_Density[2] );
	delJdotJdell   =  1/s0*(-Current_Density[0] * By*s0/r2 + Bx*s0/r2*Current_Density[1] + 0*Current_Density[2] );
  
	//	JcrossB[0] = 0.0;	JcrossB[1] = 0.0;	JcrossB[2] = 0.0;	JdotJoversgima = 0.0;
	//	delJcrossBdelm = 0.0; delJcrossBdeln = 0.0; delJcrossBdell = 0.0;
	//	delJdotJdelm = 0.0;   delJdotJdeln = 0.0;   delJdotJdell = 0.0;
  
	dens_ratio = r2/(r1+r2+r3);
	JcrossB[0] = dens_ratio*JcrossB[0]; JcrossB[1] = dens_ratio*JcrossB[1];	JcrossB[2] = dens_ratio*JcrossB[2];
	JdotJoversgima = dens_ratio*JdotJoversgima;
	delJcrossBdelm = dens_ratio*delJcrossBdelm; delJcrossBdeln = dens_ratio*delJcrossBdeln;	delJcrossBdell = dens_ratio*delJcrossBdell;
	delJdotJdelm   = dens_ratio*delJdotJdelm;	delJdotJdeln   = dens_ratio*delJdotJdeln;	delJdotJdell   = dens_ratio*delJdotJdell;
  
  
  
	//	for (iDim = 0; iDim < nDim; iDim ++)
	//		Mag_Force[1][iDim] = dens_ratio*JcrossB[iDim];
  
	SourceVector[loc+0]= -M2*AvgNum*R;
  
	SourceJacobian[loc+0][loc+0] = -M2*AvgNum*dR_dr2;
	SourceJacobian[loc+0][loc+1] = -M2*AvgNum*dR_dm2;
	SourceJacobian[loc+0][loc+2] = -M2*AvgNum*dR_dn2;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] =  -M2*AvgNum*dR_dl2;
	SourceJacobian[loc+0][loc+nDim+1] = -M2*AvgNum*dR_de2;
  
	/* Derivative of the Sixth term in the source terms with all the conservative variables */
  
	SourceVector[loc+1] = JcrossB[0] + r2*nu21*(m1/r1 - m2/r2) + r2*nu23*(m3/r3 - m2/r2);
  
	SourceJacobian[loc+1][loc+0] = nu21*m1/r1 + nu23*m3/r3 - JcrossB[0]/r2;
	SourceJacobian[loc+1][loc+1] = -nu21 - nu23 - s0 * ( Bz*Bz + By*By)/r2;
	SourceJacobian[loc+1][loc+2] = s0*By*Bx/r2;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] =  s0*Bx*Bz/r2;
	SourceJacobian[loc+1][loc+nDim+1]  = JcrossB[0] + r2*(m1/r1 - m2/r2)*dv21_dT2*dT2_de2 + r2*(m3/r3 - m2/r2)*dv23_dT2*dT2_de2;
  
  
	/* Derivative of the Seventh term in the source terms wrt all the conservative variables */
  
	SourceVector[loc+2]  =  JcrossB[1]  + r2*nu21*(n1/r1 - n2/r2) + r2*nu23*(n3/r3 - n2/r2);
  
	SourceJacobian[loc+2][loc+0] = nu21*n1/r1 + nu23*n3/r3 - JcrossB[1]/r2;
	SourceJacobian[loc+2][loc+1] = s0*Bx*By/r2;
	SourceJacobian[loc+2][loc+2] = -nu21 - nu23 -s0*(Bx*Bx + Bz*Bz)/r2;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] = s0*Bz*By/r2;
	SourceJacobian[loc+2][loc+nDim+1] = JcrossB[1] + r2*(n1/r1 - n2/r2)*dv21_dT2*dT2_de2 + r2*(n3/r3 - n2/r2)*dv23_dT2*dT2_de2;
  
	if (nDim == 3) {
		SourceVector[loc+3] = JcrossB[2]  + r2*nu21*(l1/r1 - l2/r2) + r2*nu23*(l3/r3 - l2/r2);
    
		SourceJacobian[loc+3][loc+0]  = nu21*l1/r1 + nu23*l3/r3 - JcrossB[2]/r2;
		SourceJacobian[loc+3][loc+1]  = s0*Bx*Bz/r2;
		SourceJacobian[loc+3][loc+2]  = s0*By*Bz/r2;
		SourceJacobian[loc+3][loc+3] = -nu21 - nu23 - s0*(By*By + Bx*Bx)/r2;
		SourceJacobian[loc+3][loc+nDim+1] = JcrossB[2]+ r2*(l1/r1 -l2/r2)*dv21_dT2*dT2_de2 + r2*(l3/r3 - l2/r2)*dv23_dT2*dT2_de2;
	}
  
	/* Derivative of the Eight term in the source terms wrt all the conservative variables */
  
	SourceVector[loc+nDim+1] = JcrossB[0]*m2/r2 + JcrossB[1]*n2/r2 + JcrossB[2]*l2/r2 + JdotJoversgima +
  m2*nu21*(m1/r1 - m2/r2) + m2*nu23*(m3/r3 - m2/r2)+
  n2*nu21*(n1/r1 - n2/r2) + n2*nu23*(n3/r3 - n2/r2) +
  l2*nu21*(l1/r1 - l2/r2) + l2*nu23*(l3/r3 - l2/r2) + QT2;
  
	SourceJacobian[loc+nDim+1][loc+0] = -2/r2*(JcrossB[0]*m2/r2 + JcrossB[1]*n2/r2 + JcrossB[2]*l2/r2) - 2*JdotJoversgima/r2 + (nu21+nu23)*(m2*m2 + n2*n2 + l2*l2)/(r2*r2) ;// +  dQT2_dT2*dT2_dr2;
	SourceJacobian[loc+nDim+1][loc+1] = nu21*(m1/r1 - 2*m2/r2) + nu23*(m3/r3 - 2*m2/r2) + delJcrossBdelm + delJdotJdelm + dQT2_dT2*dT2_dm2;
	SourceJacobian[loc+nDim+1][loc+2] = nu21*(n1/r1 - 2*n2/r2) + nu23*(n3/r3 - 2*n2/r2) + delJcrossBdeln + delJdotJdeln + dQT2_dT2*dT2_dn2;
	if (nDim ==3) SourceJacobian[loc+nDim+1][loc+3] = nu21*(l1/r1 - 2*l2/r2) + nu23*(l3/r3 - 2*l2/r2) + delJcrossBdell + delJdotJdell + dQT2_dT2*dT2_dl2;
  
	SourceJacobian[loc+nDim+1][loc+nDim+1] = (JcrossB[0]*m2/r2 + JcrossB[1]*n2/r2 + JcrossB[2]*l2/r2 + JdotJoversgima)/e2 +
  (n2*(n1/r1 - n2/r2) + m2*(m1/r1 - m2/r2) + l2*(l1/r1 - l2/r2))*dv21_dT2*dT2_de2 +
  (n2*(n3/r3 - n2/r2) + m2*(m3/r3 - m2/r2) + l2*(l3/r3 - l2/r2))*dv23_dT2*dT2_de2 + dQT2_dT2*dT2_de2;
  
  
  
  
  
  
  
	/* Derivative of the ninth term in the source terms wrt all the conservative variables */
	loc = 2*(nDim+2);
  
	for (iDim = 0; iDim < nDim; iDim ++ )
		Current_Density[iDim] = Electric_Conductivity*(ElectricField[iDim] + VcrossB[2][iDim]);
  
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		JcrossB[0] = Current_Density[1]*MagneticField[2] - Current_Density[2]*MagneticField[1];
		JcrossB[1] = Current_Density[2]*MagneticField[0] - Current_Density[0]*MagneticField[2];
		JcrossB[2] = Current_Density[0]*MagneticField[1] - Current_Density[1]*MagneticField[0];
	}
	JdotJoversgima = 0.0;
	for (iDim = 0; iDim < nDim; iDim ++ )
		JdotJoversgima += (Current_Density[iDim]*Current_Density[iDim]);
  
	JdotJoversgima = JdotJoversgima/s0 ;
  
	delJcrossBdelm = s0/r3/r3 * ( Bx * (l3*Bz +   n3*By) + By * (n3*Bx - 2*m3*By) + Bz * (l3*Bx - 2*m3*Bz) );
	delJcrossBdeln = s0/r3/r3 * ( Bx * (m3*By - 2*n3*Bx) + By * (m3*Bx +   l3*Bz) + Bz * (l3*By - 2*n3*Bz) );
	delJcrossBdell = s0/r3/r3 * ( Bx * (m3*Bz - 2*l3*Bx) + By * (n3*Bz - 2*l3*By) + Bz * (m3*Bx +   n3*By) );
  
	delJdotJdelm   = 1.0/s0*( Current_Density[0] * 0  	    - Bz*s0/r3*Current_Density[1]  + By*s0/r3*Current_Density[2] );
	delJdotJdeln   = 1.0/s0*( Current_Density[0] * Bz*s0/r3 - 0*Current_Density[0] 	   	   - Bx*s0/r3*Current_Density[2] );
	delJdotJdell   = 1.0/s0*(-Current_Density[0] * By*s0/r3 + Bx*s0/r3*Current_Density[1]  + 0*Current_Density[2] );
  
	for (iDim = 0; iDim < nDim; iDim ++)
		Mag_Force[2][iDim] = JcrossB[iDim];
  
	//	JcrossB[0] = 0.0;	JcrossB[1] = 0.0;	JcrossB[2] = 0.0;	JdotJoversgima = 0.0;
	//	delJcrossBdelm = 0.0; delJcrossBdeln = 0.0; delJcrossBdell = 0.0;
	//	delJdotJdelm = 0.0;   delJdotJdeln = 0.0;   delJdotJdell = 0.0;
  
	dens_ratio = r3/(r1+r2+r3);
	JcrossB[0] = dens_ratio*JcrossB[0]; JcrossB[1] = dens_ratio*JcrossB[1];	JcrossB[2] = dens_ratio*JcrossB[2];
	JdotJoversgima = dens_ratio*JdotJoversgima;
	delJcrossBdelm = dens_ratio*delJcrossBdelm; delJcrossBdeln = dens_ratio*delJcrossBdeln;	delJcrossBdell = dens_ratio*delJcrossBdell;
	delJdotJdelm   = dens_ratio*delJdotJdelm;	delJdotJdeln   = dens_ratio*delJdotJdeln;	delJdotJdell   = dens_ratio*delJdotJdell;
  
	SourceVector[loc+0] = -M3*AvgNum*R;
  
	SourceJacobian[loc+0][loc+0] = -M3*AvgNum*dR_dr3;
	SourceJacobian[loc+0][loc+1] = -M3*AvgNum*dR_dm3;
	SourceJacobian[loc+0][loc+2] = -M3*AvgNum*dR_dn3;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] = -M3*AvgNum*dR_dl3;
	SourceJacobian[loc+0][loc+nDim+1] = -M3*AvgNum*dR_de3;
  
	/* Derivative of the Tenth term in the source terms wrt all the conservative variables */
  
	SourceVector[loc+1] =  JcrossB[0]  + r3*nu31*(m1/r1 - m3/r3) + r3*nu32*(m2/r2 - m3/r3);
  
	SourceJacobian[loc+1][loc+0] = nu31*m1/r1 + nu32*m2/r2 - JcrossB[0]/r3;
	SourceJacobian[loc+1][loc+1] =  -nu31 - nu32- s0 * ( Bz*Bz + By*By)/r3;
	SourceJacobian[loc+1][loc+2] =   s0*Bx*Bz/r3;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] =  s0*Bx*Bz/r3;
	SourceJacobian[loc+1][loc+nDim+1] =  JcrossB[0] + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_de3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_de3;
  
	/* Derivative of the Eleventh term in the source terms wrt all the conservative variables */
	SourceVector[loc+2]  = JcrossB[1]  + r3*nu31*(n1/r1 - n3/r3) + r3*nu32*(n2/r2 - n3/r3);
  
	SourceJacobian[loc+2][loc+0] = nu31*n1/r1 + nu32*n2/r2 - JcrossB[1]/r3;
	SourceJacobian[loc+2][loc+1] = s0*Bx*By/r3;
	SourceJacobian[loc+2][loc+2] =  -nu31 - nu32 -s0*(Bx*Bx + Bz*Bz)/r3;
	if (nDim ==3)  SourceJacobian[loc+2][loc+3] =  s0*Bz*By/r3 ;
	SourceJacobian[loc+2][loc+nDim+1] = JcrossB[1] + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_de3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_de3;
  
	/* Derivative of the Twelfth term in the source terms wrt all the conservative variables */
	if (nDim ==3) {
		SourceVector[loc+3]  = JcrossB[2] + r3*nu31*(l1/r1 - l3/r3) + r3*nu32*(l2/r2 - l3/r3);
		SourceJacobian[loc+3][loc+0] = nu31*l1/r1 + nu32*l2/r2 - JcrossB[2]/r3;
		SourceJacobian[loc+3][loc+1] =  s0*Bx*Bz/r3;
		SourceJacobian[loc+3][loc+2] =  s0*By*Bz/r3;
		SourceJacobian[loc+3][loc+3] =  -nu31 - nu32 - s0*(By*By + Bx*Bx)/r3 ;
		SourceJacobian[loc+3][loc+nDim+1] =  JcrossB[2] + r3*(l1/r1 - l3/r3)*dv31_dT3*dT3_de3 + r3*(l2/r2 - l3/r3)*dv32_dT3*dT3_de3;
	}
  
	SourceVector[loc+nDim+1]  = JcrossB[0]*m3/r3 + JcrossB[1]*n3/r3 + JcrossB[2]*l3/r3 + JdotJoversgima +
  m3*nu31*(m1/r1 - m3/r3) + m3*nu32*(m2/r2 - m3/r3) +
  n3*nu31*(n1/r1 - n3/r3) + n3*nu32*(n2/r2 - n3/r3) +
  l3*nu31*(l1/r1 - l3/r3) + l3*nu32*(l2/r2 - l3/r3) + QT3;
  
	SourceJacobian[loc+nDim+1][loc+0] =  -2/r3*(JcrossB[0]*m3/r3 + JcrossB[1]*n3/r3 + JcrossB[2]*l3/r3) - 2*JdotJoversgima/r3 + (nu31+nu32)*(m3*m3 + n3*n3 +l3*l3)/(r3*r3);// + dQT3_dT3*dT3_dr3
	SourceJacobian[loc+nDim+1][loc+1]  = nu31*(m1/r1 - 2*m3/r3) + nu32*(m2/r2 - 2*m3/r3) + delJcrossBdelm  + delJdotJdelm + dQT3_dT3*dT3_dm3;
	SourceJacobian[loc+nDim+1][loc+2]  = nu31*(n1/r1 - 2*n3/r3) + nu32*(n2/r2 - 2*n3/r3) + delJcrossBdeln  + delJdotJdeln + dQT3_dT3*dT3_dn3;
	if (nDim ==3)SourceJacobian[loc+nDim+1][loc+3]  = nu31*(l1/r1 - 2*l3/r3) + nu32*(l2/r2 - 2*l3/r3) + delJcrossBdell + delJdotJdell  + dQT3_dT3*dT3_dl3;
	SourceJacobian[loc+nDim+1][loc+nDim+1]  = (JcrossB[0]*m3/r3 + JcrossB[1]*n3/r3 + JcrossB[2]*l3/r3 + JdotJoversgima)/e3 +
  (m3*(m1/r1 - m3/r3) +  n3*(n1/r1 - n3/r3) +  l3*(l1/r1 - l3/r3))*dv31_dT3*dT3_de3 +
  (m3*(m2/r2 - m3/r3)+  n3*(n2/r2 - n3/r3) + l3*(l2/r2 - l3/r3))*dv32_dT3*dT3_de3 + dQT3_dT3*dT3_de3;
  
  
  
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		JcrossB[0] = Current_Density[1]*MagneticField[2] - Current_Density[2]*MagneticField[1];
		JcrossB[1] = Current_Density[2]*MagneticField[0] - Current_Density[0]*MagneticField[2];
		JcrossB[2] = Current_Density[0]*MagneticField[1] - Current_Density[1]*MagneticField[0];
	}
  
	for (unsigned short iVar = 0; iVar < nVar; iVar ++) {
		val_residual[iVar] = SourceVector[iVar]*Volume;
		if (implicit) {
			for (unsigned short jVar = 0; jVar < nVar; jVar ++) {
				val_Jacobian_i[iVar][jVar] = SourceJacobian[iVar][jVar]*Volume;
			}
		}
	}
}

void CSourcePieceWise_Plasma::ComputeResidual_Axisymmetric(double *val_residual, CConfig *config) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CSourcePieceWise_Plasma::ComputeResidual_Axisymmetric
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *Coord_i Volume
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim nVar nSpecies nDiatomics
  
	//SU2_CPP2C DECL_LIST START
  
	//SU2_CPP2C DECL_LIST END
  
  
	double yinv, Enthalpy_i, Velocity_i, sq_vel;
	unsigned short iDim, iSpecies, loc;
  
	if (Coord_i[1] > 0.0) yinv = 1.0/Coord_i[1];
	else yinv = 0.0;
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		sq_vel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i = U_i[loc+iDim+1] / U_i[loc+0];
			sq_vel += Velocity_i *Velocity_i;
		}
    
		/*--- Set the primitive variables (T_tr,    vel,  T_vib,      p,    rho,      h,      c) and compute the mean ---*/
		/*--- Indices:                    (   0, iDim+1, nDim+1, nDim+2, nDim+3, nDim+4, nDim+5)                      ---*/
		Enthalpy_i = Varray_i[iSpecies][nDim+4];
    
		val_residual[loc+0] = yinv*Volume*U_i[loc+2];
		val_residual[loc+1] = yinv*Volume*U_i[loc+1]*U_i[loc+2]/U_i[loc+0];
		val_residual[loc+2] = yinv*Volume*U_i[loc+2]*U_i[loc+2]/U_i[loc+0];
		val_residual[loc+3] = yinv*Volume*Enthalpy_i*U_i[loc+2];
		if (iSpecies < nDiatomics)
			val_residual[loc+4] = yinv*Volume*U_i[loc+4]*U_i[loc+2]/U_i[loc+0];
	}
	//SU2_CPP2C END CSourcePieceWise_Plasma::ComputeResidual_Axisymmetric
}


void CSourcePieceWise_Plasma::SetJacobian_Axisymmetric(double **val_Jacobian_i, CConfig *config) {
  
	unsigned short iVar, jVar;
  
	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
		Residual_New[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++)
			val_Jacobian_i[iVar][jVar] = 0.0;
	}
  
	/*--- Calculate jacobian using automatic differentiation ---*/
	if (config->GetKind_SourJac_Plasma() == AUTO_DIFF) {
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) U_id[jVar] = 0.0;
			U_id[iVar] = 1.0;
			ComputeResidual_Axisymmetric_ad(Residual_Baseline, Residual_New, config);
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_i[jVar][iVar] = Residual_New[jVar];
		}
	}
  
	/*--- Calculate jacobian using finite differencing ---*/
	else if (config->GetKind_SourJac_Plasma() == FINITE_DIFF) {
		double FDEpsilon;
    
		/*--- Compute the baseline line residual ---*/
		ComputeResidual_Axisymmetric(Residual_Baseline, config);
    
		/*--- Compute forward finite differences ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			FDEpsilon = 1E-6*U_Baseline[iVar];
			if (FDEpsilon < 1E-14)
				FDEpsilon = 1E-14;
      
			/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
			U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
			for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
			ComputeResidual_Axisymmetric(Residual_New, config);
      
			/*--- Undo the change in U_i to keep constant the solution ---*/
			U_i[iVar] = U_i[iVar] - FDEpsilon;
      
			/*--- Save the new Jacobian (per column) ---*/
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
			}
		}
	}
}


void CSourcePieceWise_Plasma::SetJacobian_Chemistry(double **val_Jacobian_i, CConfig *config) {
  
	unsigned short iVar, jVar;
  
	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
		Residual_New[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++)
			val_Jacobian_i[iVar][jVar] = 0.0;
	}
  
	/*--- Calculate jacobian using automatic differentiation ---*/
	if (config->GetKind_SourJac_Plasma() == AUTO_DIFF) {
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) U_id[jVar] = 0.0;
			U_id[iVar] = 1.0;
			ComputeResidual_Chemistry_ad(Residual_Baseline, Residual_New, config);
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_i[jVar][iVar] = Residual_New[jVar];
		}
	}
  
	/*--- Calculate jacobian using finite differencing ---*/
	else if (config->GetKind_SourJac_Plasma() == FINITE_DIFF) {
		double FDEpsilon;
    
		/*--- Compute the baseline line residual ---*/
		ComputeResidual_Chemistry(Residual_Baseline, config);
    
		/*--- Compute forward finite differences ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			FDEpsilon = 1E-6*U_Baseline[iVar];
			if (FDEpsilon < 1E-14)
				FDEpsilon = 1E-14;
      
			/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
			U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
			for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
			ComputeResidual_Chemistry(Residual_New, config);
      
			/*--- Undo the change in U_i to keep constant the solution ---*/
			U_i[iVar] = U_i[iVar] - FDEpsilon;
      
			/*--- Save the new Jacobian (per column) ---*/
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
			}
		}
	}
}


void CSourcePieceWise_Plasma::GetEq_Rxn_Coefficients(double **EqRxnConstants, CConfig *config) {
	unsigned short iSpecies, iReaction, loc, ii;
	double mixNumDensity, spNumDensity, interpFact;
  
	// Park 90 Chemistry model implemented in Scalabrin.
  
	/*--- Calculate the mixture number density at the current node ---*/
	mixNumDensity = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		spNumDensity = U_i[loc+0]/Molar_Mass[iSpecies]*AVOGAD_CONSTANT;
		mixNumDensity += spNumDensity;
	}
  
	/*--- Retrieve constants for a specified reaction ---*/
	for (iReaction = 0; iReaction < nReactions; iReaction++) {
		config->GetChemistryEquilConstants(RxnConstantTable, iReaction);
    
		/*--- Interpolate table values for current mixture number density ---*/
		if (mixNumDensity < 1E14) {
			for (ii = 0; ii < 5; ii++)
				EqRxnConstants[iReaction][ii] = RxnConstantTable[0][ii];
		} else if (mixNumDensity >= 1E14 && mixNumDensity < 1E15) {
			interpFact = (mixNumDensity - 1E14) / (1E15 - 1E14);
			for (ii = 0; ii < 5; ii++)
				EqRxnConstants[iReaction][ii] = interpFact*(RxnConstantTable[1][ii]-RxnConstantTable[0][ii]) + RxnConstantTable[0][ii];
		} else if (mixNumDensity >= 1E15 && mixNumDensity < 1E16) {
			interpFact = (mixNumDensity - 1E15) / (1E16 - 1E15);
			for (ii = 0; ii < 5; ii++)
				EqRxnConstants[iReaction][ii] = interpFact*(RxnConstantTable[2][ii]-RxnConstantTable[1][ii]) + RxnConstantTable[1][ii];
		} else if (mixNumDensity >= 1E16 && mixNumDensity < 1E17) {
			interpFact = (mixNumDensity - 1E16) / (1E17 - 1E16);
			for (ii = 0; ii < 5; ii++)
				EqRxnConstants[iReaction][ii] = interpFact*(RxnConstantTable[3][ii]-RxnConstantTable[2][ii]) + RxnConstantTable[2][ii];
		} else if (mixNumDensity >= 1E17 && mixNumDensity < 1E18) {
			interpFact = (mixNumDensity - 1E17) / (1E18 - 1E17);
			for (ii = 0; ii < 5; ii++)
				EqRxnConstants[iReaction][ii] = interpFact*(RxnConstantTable[4][ii]-RxnConstantTable[3][ii]) + RxnConstantTable[3][ii];
		} else if (mixNumDensity >= 1E18 && mixNumDensity < 1E19) {
			interpFact = (mixNumDensity - 1E18) / (1E19 - 1E18);
			for (ii = 0; ii < 5; ii++)
				EqRxnConstants[iReaction][ii] = interpFact*(RxnConstantTable[5][ii]-RxnConstantTable[4][ii]) + RxnConstantTable[4][ii];
		} else if (mixNumDensity >= 1E19) {
			for (ii = 0; ii < 5; ii++)
				EqRxnConstants[iReaction][ii] = RxnConstantTable[5][ii];
		}
    
	}
}

void CSourcePieceWise_Plasma::ComputeResidual_Chemistry(double *val_residual, CConfig *config) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CSourcePieceWise_Plasma::ComputeResidual_Chemistry
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *Molecular_Diameter *Molar_Mass Volume *Cf *eta *theta **EqRxnConstants
	//SU2_CPP2C VARS INT ***Reactions
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim nVar nSpecies nReactions nDiatomics
  
	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nSpecies Temp_tr_i
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nReactions SIZE=5 EqRxnConstants
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nReactions T_rxnf T_rxnb ReactionRateFwd ReactionRateBkw Keq fwdRxn bkwRxn
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nSpecies w_dot
  
	//SU2_CPP2C DECL_LIST END
  
	double T_min;						// Minimum temperature for the modified temperature calculations.
	double epsilon;					// Parameter for the modified temperature calculations.
	unsigned short iSpecies, jSpecies, iReaction, iVar, iDim, ii;
	unsigned short iLoc, jLoc;
	unsigned short counterFwd, counterBkw;
	double T_tr;
  
	T_min = 800;
	epsilon = 80;
  
	/*--- Initialize all components of the residual and Jacobian to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.0;
	}
  
	for (iReaction = 0; iReaction < nReactions; iReaction++) {
    
		/*--- Calculate the rate-controlling temperatures ---*/
		//NOTE:  This implementation takes the geometric mean of the TR temps of all particpants.
		//       This is NOT consistent with the Park chemical models, and is only a first-cut approach.
		T_rxnf[iReaction] = 1.0;
		T_rxnb[iReaction] = 1.0;
		counterFwd = 0;
		counterBkw = 0;
		for (ii = 0; ii < 3; ii++) {
			iSpecies = Reactions[iReaction][0][ii];
			jSpecies = Reactions[iReaction][1][ii];
			/*--- Reactants ---*/
			if (iSpecies != nSpecies) {
				T_tr = Varray_i[iSpecies][0];
				T_rxnf[iReaction] *= T_tr;
				counterFwd++;
			}
			/*--- Products ---*/
			if (jSpecies != nSpecies) {
				T_tr = Varray_i[jSpecies][0];
				T_rxnb[iReaction] *= T_tr;
				counterBkw++;
			}
		}
		T_rxnf[iReaction] = exp(1.0/counterFwd*log(T_rxnf[iReaction]));
		T_rxnb[iReaction] = exp(1.0/counterBkw*log(T_rxnb[iReaction]));
    
		/*--- Apply a modified temperature to ease the stiffness at low temperatures ---*/
		T_rxnf[iReaction] = 0.5 * (T_rxnf[iReaction]+T_min + sqrt((T_rxnf[iReaction]-T_min)*(T_rxnf[iReaction]-T_min) + epsilon*epsilon));
		T_rxnb[iReaction] = 0.5 * (T_rxnb[iReaction]+T_min + sqrt((T_rxnb[iReaction]-T_min)*(T_rxnb[iReaction]-T_min) + epsilon*epsilon));
    
		GetEq_Rxn_Coefficients(EqRxnConstants, config);
    
		/*--- Calculate equilibrium extent of reaction ---*/
		//NOTE: Scalabrin implementation
		Keq[iReaction] = exp(EqRxnConstants[iReaction][0]*(T_rxnb[iReaction]/10000.0)
                         + EqRxnConstants[iReaction][1]
                         + EqRxnConstants[iReaction][2]*log(10000.0/T_rxnb[iReaction])
                         + EqRxnConstants[iReaction][3]*(10000.0/T_rxnb[iReaction])
                         + EqRxnConstants[iReaction][4]*(10000.0/T_rxnb[iReaction])*(10000.0/T_rxnb[iReaction]));
    
		/*--- Calculate reaction rate coefficients ---*/
		ReactionRateFwd[iReaction] = Cf[iReaction] * exp(eta[iReaction]*log(T_rxnf[iReaction])) * exp(-theta[iReaction]/T_rxnf[iReaction]);
		ReactionRateBkw[iReaction] = Cf[iReaction] * exp(eta[iReaction]*log(T_rxnb[iReaction])) * exp(-theta[iReaction]/T_rxnb[iReaction]) / Keq[iReaction];
    
		fwdRxn[iReaction] = 1.0;
		bkwRxn[iReaction] = 1.0;
		for (ii = 0; ii < 3; ii++) {
			/*--- Reactants ---*/
			iSpecies = Reactions[iReaction][0][ii];
			if ( iSpecies != nSpecies) {
				if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
				else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				fwdRxn[iReaction] *= 0.001*U_i[iLoc+0]/Molar_Mass[iSpecies];
			}
			/*--- Products ---*/
			jSpecies = Reactions[iReaction][1][ii];
			if (jSpecies != nSpecies) {
				if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
				else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
				bkwRxn[iReaction] *= 0.001*U_i[jLoc+0]/Molar_Mass[jSpecies];
			}
		}
		fwdRxn[iReaction] = 1000.0 * ReactionRateFwd[iReaction] * fwdRxn[iReaction];
		bkwRxn[iReaction] = 1000.0 * ReactionRateBkw[iReaction] * bkwRxn[iReaction];
	}
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) w_dot[iSpecies] = 0.0;
  
	for (iReaction = 0; iReaction < nReactions; iReaction++) {
		for (ii = 0; ii < 3; ii++) {
			/*--- Products ---*/
			iSpecies = Reactions[iReaction][1][ii];
			if (iSpecies != nSpecies)
				w_dot[iSpecies] += Molar_Mass[iSpecies] * (fwdRxn[iReaction] - bkwRxn[iReaction]);
			/*--- Reactants ---*/
			iSpecies = Reactions[iReaction][0][ii];
			if (iSpecies != nSpecies)
				w_dot[iSpecies] -= Molar_Mass[iSpecies] * (fwdRxn[iReaction] - bkwRxn[iReaction]);
		}
	}
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		val_residual[iLoc] = w_dot[iSpecies]*Volume;
		for (iDim = 0; iDim < nDim; iDim++)
			val_residual[iLoc+1+iDim] = w_dot[iSpecies] * U_i[iLoc+1+iDim]/U_i[iLoc] * Volume;
		val_residual[iLoc+1+nDim] = w_dot[iSpecies] * U_i[iLoc+1+nDim]/U_i[iLoc] * Volume;
		if (iSpecies < nDiatomics) {
			val_residual[iLoc+2+nDim] = w_dot[iSpecies] * U_i[iLoc+2+nDim]/U_i[iLoc] * Volume;
		}
	}
	//SU2_CPP2C END CSourcePieceWise_Plasma::ComputeResidual_Chemistry
}


void CSourcePieceWise_Plasma::SetJacobian_ElecForce(double **val_Jacobian_i, CConfig *config) {
  
	unsigned short iVar, jVar;
  
	/*--- In a near future AD should be used instead of FD ---*/
	double FDEpsilon;
  
	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
		Residual_New[iVar] = 0.0;
	}
  
	/*--- Compute the baseline line residual ---*/
  
	ComputeResidual_ElecForce(Residual_Baseline, config);
  
	/*--- Compute forward finite differences ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		FDEpsilon = 0.01*U_Baseline[iVar];
		if (FDEpsilon == 0)
			FDEpsilon = 1E-9;
    
		/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
		U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
		for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
		ComputeResidual_ElecForce(Residual_New, config);
    
		/*--- Undo the change in U_i to keep constant the solution ---*/
		U_i[iVar] = U_i[iVar] - FDEpsilon;
    
		/*--- Save the new Jacobian (per column) ---*/
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
		}
	}
}

void CSourcePieceWise_Plasma::ComputeResidual_ElecForce(double *val_residual, CConfig *config) {
  
}

void CSourcePieceWise_Plasma::SetJacobian_MomentumExch(double **val_Jacobian_i, CConfig *config) {
  
	unsigned short iVar, jVar;
  
	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
		Residual_New[iVar] = 0.0;
	}
  
	/*--- Calculate jacobian using automatic differentiation ---*/
	if (config->GetKind_SourJac_Plasma() == AUTO_DIFF) {
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) U_id[jVar] = 0.0;
			U_id[iVar] = 1.0;
			ComputeResidual_MomentumExch_ad(Residual_Baseline, Residual_New, config);
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_i[jVar][iVar] = Residual_New[jVar];
		}
	}
  
	/*--- Calculate jacobian using finite differencing ---*/
	else if (config->GetKind_SourJac_Plasma() == FINITE_DIFF) {
    
		double FDEpsilon;
    
		/*--- Compute the baseline line residual ---*/
		ComputeResidual_MomentumExch(Residual_Baseline, config);
    
		/*--- Compute forward finite differences ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			FDEpsilon = 1E-6*U_Baseline[iVar];
			if (FDEpsilon < 1E-14)
				FDEpsilon = 1E-14;
      
			/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
			U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
			for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
			ComputeResidual_MomentumExch(Residual_New, config);
      
			/*--- Undo the change in U_i to keep constant the solution ---*/
			U_i[iVar] = U_i[iVar] - FDEpsilon;
      
			/*--- Save the new Jacobian (per column) ---*/
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
			}
		}
	}
}

void CSourcePieceWise_Plasma::ComputeResidual_MomentumExch(double *val_residual, CConfig *config) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CSourcePieceWise_Plasma::ComputeResidual_MomentumExch
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *Molecular_Diameter *Molecular_Mass ***Omega11 Volume
	//SU2_CPP2C VARS INT *ChargeNumber
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim nVar nSpecies nDiatomics PI_NUMBER ELECTRON_CHARGE FREE_PERMITTIVITY BOLTZMANN_CONSTANT
  
	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nSpecies Temp_tr_i
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nSpecies SIZE=nDim P
	//SU2_CPP2C DECL_LIST END
  
	unsigned short iDim, iSpecies, jSpecies, iVar, iLoc, jLoc;
	double collisionFreq, collisionArea, velocity_Species_i, velocity_Species_j, collisionVelocity;
	double T_control;
	double radius_electronIonCollision, electron_temperature;
	double T_tr_i, T_tr_j;
  
	/*--- Initialization ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		for (iDim = 0; iDim < nDim; iDim++)
			P[iSpecies][iDim] = 0.0;
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = 0.0;
  
	electron_temperature = Varray_i[nSpecies-1][0];
  
	/*--- Solve for momentum exchange between species due to collisions ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
			if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
			else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
      
			if (iSpecies != jSpecies) {
				T_tr_i = Varray_i[iSpecies][0];
				T_tr_j = Varray_i[jSpecies][0];
				T_control = sqrt(T_tr_i*T_tr_j);
				collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] * pow(T_control, Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2]);
        
				/*----- An if-condition for the case when a positive charge collides with an electron -----*/
				if( (ChargeNumber[iSpecies] == 1 && jSpecies == nSpecies -1) || (iSpecies == nSpecies-1 && ChargeNumber[jSpecies] == 1)) {
					radius_electronIonCollision = ELECTRON_CHARGE*ELECTRON_CHARGE/(32.0*FREE_PERMITTIVITY*BOLTZMANN_CONSTANT*electron_temperature);
					collisionArea = PI_NUMBER * radius_electronIonCollision * radius_electronIonCollision;
				}
        
				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*T_tr_i/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*T_tr_j/(PI_NUMBER*Molecular_Mass[jSpecies]));
				collisionVelocity  = sqrt(velocity_Species_i*velocity_Species_i + velocity_Species_j*velocity_Species_j);
				collisionFreq = (U_i[jLoc+0] * collisionArea * collisionVelocity)/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);
        
				for (iDim = 0; iDim < nDim; iDim++)
					P[iSpecies][iDim]  +=  U_i[iLoc+0]* collisionFreq * (U_i[jLoc+1+iDim]/U_i[jLoc+0] - U_i[iLoc+1+iDim]/U_i[iLoc+0]);
        
			}
		}
		for (iDim = 0; iDim < nDim; iDim++) {
			val_residual[iLoc+1+iDim] = P[iSpecies][iDim] * Volume;
			val_residual[iLoc+nDim + 1] += P[iSpecies][iDim]* U_i[iLoc+1+iDim]/U_i[iLoc+0] * Volume;
		}
	}
	//SU2_CPP2C END CSourcePieceWise_Plasma::ComputeResidual_MomentumExch
}


void CSourcePieceWise_Plasma::SetJacobian_EnergyExch(double **val_Jacobian_i, CConfig *config) {
  
	unsigned short iVar, jVar;
  
	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
		Residual_New[iVar] = 0.0;
	}
  
	/*--- Calculate jacobian using automatic differentiation ---*/
	if (config->GetKind_SourJac_Plasma() == AUTO_DIFF) {
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) U_id[jVar] = 0.0;
			U_id[iVar] = 1.0;
			ComputeResidual_EnergyExch_ad(Residual_Baseline, Residual_New, config);
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_i[jVar][iVar] = Residual_New[jVar];
		}
	}
  
	/*--- Calculate jacobian using finite differencing ---*/
	else if (config->GetKind_SourJac_Plasma() == FINITE_DIFF) {
		double FDEpsilon;
    
		/*--- Compute the baseline line residual ---*/
		ComputeResidual_EnergyExch(Residual_Baseline, config);
    
		/*--- Compute forward finite differences ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			FDEpsilon = 1E-6*U_Baseline[iVar];
			if (FDEpsilon < 1E-14)
				FDEpsilon = 1E-14;
      
			/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
			U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
			for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
			ComputeResidual_EnergyExch(Residual_New, config);
      
			/*--- Undo the change in U_i to keep constant the solution ---*/
			U_i[iVar] = U_i[iVar] - FDEpsilon;
      
			/*--- Save the new Jacobian (per column) ---*/
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
			}
		}
	}
}


void CSourcePieceWise_Plasma::ComputeResidual_EnergyExch(double *val_residual, CConfig *config) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CSourcePieceWise_Plasma::ComputeResidual_EnergyExch
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *Molecular_Diameter *Molecular_Mass *Molar_Mass Volume ***Omega11
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim nVar nSpecies nDiatomics PI_NUMBER BOLTZMANN_CONSTANT UNIVERSAL_GAS_CONSTANT
  
	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nSpecies SpeciesPressure_i Temp_tr_i CharVibTemp
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nSpecies Temp_tr_i
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nSpecies Q_elastic Q_tv
	//SU2_CPP2C DECL_LIST END
  
	unsigned short iDim, iSpecies, jSpecies, iLoc, jLoc, iVar;
	double collisionFreq, collisionArea, velocity_Species_i, velocity_Species_j, collisionVelocity, vel_dot_prod;
	double coefficient;
	double T_control;
	double T_tr_i, T_tr_j;
  
	/*--- Energy transport via elastic collisions ---*/
	//Comment: From Lee 1985, and originally from Engineering Magnetohydrodynamics by Sutton (1965)
	/*Comment: Two terms, both from kinetic theory.  The first term accounts for changes in species temperature from
	 elastic encounters between particles, the second term accounts for frictional heating between species*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Q_elastic[iSpecies] = 0.0;
		Q_tv[iSpecies] = 0.0;
	}
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.0;
	}
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){
			if (jSpecies != iSpecies) {
				if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
				else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
        
				T_tr_i = Varray_i[iSpecies][0];
				T_tr_j = Varray_i[jSpecies][0];
				T_control = sqrt(T_tr_i*T_tr_j);
				collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] * pow(T_control, Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2]);
        
				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*T_tr_i/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*T_tr_j/(PI_NUMBER*Molecular_Mass[jSpecies]));
				collisionVelocity = sqrt(velocity_Species_i*velocity_Species_i + velocity_Species_j*velocity_Species_j);
				collisionFreq = (U_i[jLoc+0] * collisionArea * collisionVelocity)/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);
				vel_dot_prod = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					vel_dot_prod += (U_i[iLoc+1+iDim] - U_i[jLoc+1+iDim])*U_i[iLoc+1+iDim];
				}
        
				/* Exchange from Lee and Sutton, heavy particles */
				coefficient = 2.0*U_i[iLoc+0] / (Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]) * collisionFreq * 3.0/2.0*BOLTZMANN_CONSTANT;
				Q_elastic[iSpecies] += coefficient*(T_tr_j - T_tr_i);
			}
		}
		val_residual[iLoc+1+nDim] =  Q_elastic[iSpecies] * Volume;
	}
  
	/*--- Translational-rotational & vibrational energy exchange via inelastic collisions ---*/
	//Comment: Landau-Teller formulation
	//Comment: Based on Scalabrin and Candler.  May need to re-visit with much more detail and derive for my formulation
  
	double tau_sr, tau_ps, LimitingXSection, AvgMolecularSpeed, ReducedMass, A_sr, estar_vs, e_vs, q_tr_vs;
	double MixturePressure, MixtureNumDensity;
	double P_i, P_j;
  
	/*--- Calculate mixture quantities ---*/
	MixturePressure   = 0.0;
	MixtureNumDensity = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) {
			iLoc = (nDim+3)*iSpecies;
			CharVibTemp[iSpecies] = config->GetCharVibTemp(iSpecies);
		} else
			iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		MixturePressure += Varray_i[iSpecies][nDim+2];//SpeciesPressure_i[iSpecies];
		MixtureNumDensity += U_i[iLoc+0]/Molecular_Mass[iSpecies];
	}
  
	for (iSpecies = 0; iSpecies < nDiatomics; iSpecies++) {
		iLoc = (nDim+3)*iSpecies;
		Q_tv[iSpecies] = 0.0;
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){
			if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
			else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
      
			T_tr_i = Varray_i[iSpecies][0];
			T_tr_j = Varray_i[jSpecies][0];
			P_i = Varray_i[iSpecies][nDim+2];
			P_j = Varray_i[jSpecies][nDim+2];
      
			/*--- Calculate Landau-Teller relaxation time ---*/
			ReducedMass = Molar_Mass[iSpecies]*Molar_Mass[jSpecies]/(Molar_Mass[iSpecies]+Molar_Mass[jSpecies]);
			LimitingXSection = 1E-20*pow(50000/T_tr_i, 2.0);
      
			AvgMolecularSpeed = sqrt(8.0*UNIVERSAL_GAS_CONSTANT*T_tr_i/(PI_NUMBER*Molar_Mass[iSpecies]));
      
			A_sr = 1.16 * 1E-3 * sqrt(ReducedMass) * pow(CharVibTemp[iSpecies], 4.0/3.0);
			tau_sr =   101325/(P_i+P_j)* exp(A_sr*(pow(sqrt(T_tr_i*T_tr_j),-1.0/3.0) - 0.015*pow(ReducedMass,0.25)) - 18.42);
			tau_ps = 1.0/(LimitingXSection*AvgMolecularSpeed*MixtureNumDensity);
      
			estar_vs = UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies] * CharVibTemp[iSpecies]/(exp(CharVibTemp[iSpecies]/T_tr_j)-1.0);
			e_vs = U_i[iLoc+nDim+2]/U_i[iLoc+0];
      
			/*--- Energy transferred per-molecule from species r to species s ---*/
			q_tr_vs = Molecular_Mass[iSpecies] * (estar_vs - e_vs)/(tau_sr + tau_ps);
      
			/*--- Convert to energy per volume for r and s species and multiply by volume for residual ---*/
			val_residual[iLoc+nDim+2] += q_tr_vs*U_i[iLoc+0]/Molecular_Mass[iSpecies]*Volume;
			val_residual[iLoc+nDim+1] += q_tr_vs*U_i[iLoc+0]/Molecular_Mass[iSpecies]*Volume;
			val_residual[jLoc+nDim+1] -= q_tr_vs*U_i[jLoc+0]/Molecular_Mass[jSpecies]*Volume;
		}
	}
	//SU2_CPP2C END CSourcePieceWise_Plasma::ComputeResidual_EnergyExch
}


CSource_Magnet::CSource_Magnet(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
	MagneticField = new double [3];
	MagneticDipole = new double [3];
	Current_Density = new double [3];
	JcrossB = new double [3];
	VcrossB = new double [3];
	Electric_Conductivity = config->GetElec_Conductivity();
	vector_r = new double [nDim];
	dpcenter = new double [nDim];
	velocity = new double[nDim];
  
	for (iVar = 0; iVar < 3; iVar ++) {
		MagneticField[iVar] = 0.0;
		MagneticDipole[iVar] = 0.0;
		Current_Density[iVar] = 0.0;
		JcrossB[iVar] = 0.0;
		VcrossB[iVar] = 0.0;
	}
  
	Stagnation_B = config->GetStagnation_B();
	dpcenter[0] = config->GetDipoleDist();
	dpcenter[1] = 0.0;
	dpcenter[2] = 0.0;
  
	MagneticDipole[0] = Stagnation_B*2*PI_NUMBER*pow((0.01825+dpcenter[0]),3)/MAGNETIC_CONSTANT;
	MagneticDipole[1] = 0.0;
	MagneticDipole[2] = 0.0;
  
	if (nDim ==2) VcrossB[iDim] = 0.0;
  
  
}

CSource_Magnet::~CSource_Magnet(void) {
	delete []	MagneticField;
	delete []	MagneticDipole;
	delete []	Current_Density;
	delete []	JcrossB;
	delete []   VcrossB;
	delete []   velocity;
}

void CSource_Magnet::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
  
	double mdotr, distance, rpower5, rcubed;
	double r1,m1,n1,l1,e1;
	r1 = U_i[0];	// density of each species
	m1 = U_i[1];	// density*horizontal velocity
	n1 = U_i[2];	// density*vertical velocity
	e1 = U_i[3];	// energy per unit volume
	l1 = 0.0;			// density*vertical velocity of each species
  
	if (nDim ==3 ) {
		l1 = U_i[3];	// density*vertical velocity
		e1 = U_i[4];
	}
  
	velocity[0] = m1/r1; 	velocity[1] = n1/r1;
	if (nDim ==3) velocity[2] = l1/r1;
  
  
	mdotr = 0.0;
	distance = 0.0;
  
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		vector_r[iDim] = Coord_i[iDim] - dpcenter[iDim];
		mdotr +=MagneticDipole[iDim]*vector_r[iDim];
		distance += vector_r[iDim]*vector_r[iDim];
	}
	distance = sqrt(distance);
	rpower5 = pow(distance,5);
	rcubed  = pow(distance,3);
	for (iDim = 0; iDim < nDim; iDim ++ )
		MagneticField[iDim] = MAGNETIC_CONSTANT/(4.0*PI_NUMBER) * (3.0 * vector_r[iDim]* mdotr/rpower5 - MagneticDipole[iDim]/rcubed);
  
	double Bx = MagneticField[0];
	double By = MagneticField[1];
	double Bz = MagneticField[2];
	double s0 = Electric_Conductivity;
  
	if (nDim ==3) {
		VcrossB[0] = velocity[1]*MagneticField[2] - velocity[2]*MagneticField[1];
		VcrossB[1] = velocity[2]*MagneticField[0] - velocity[0]*MagneticField[2];
		VcrossB[2] = velocity[0]*MagneticField[1] - velocity[1]*MagneticField[0];
	}
  
	for (iDim = 0; iDim < nDim; iDim ++ )
		Current_Density[iDim] = Electric_Conductivity*VcrossB[iDim];
  
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		JcrossB[0] = Current_Density[1]*MagneticField[2] - Current_Density[2]*MagneticField[1];
		JcrossB[1] = Current_Density[2]*MagneticField[0] - Current_Density[0]*MagneticField[2];
		JcrossB[2] = Current_Density[0]*MagneticField[1] - Current_Density[1]*MagneticField[0];
	}
  
	double JdotJoversgima = 0.0;
  
	for (iDim = 0; iDim < nDim; iDim ++ )
		JdotJoversgima += (Current_Density[iDim]*Current_Density[iDim]);
  
	JdotJoversgima = JdotJoversgima/s0 ;
  
	double delJcrossBdelm = s0/r1/r1 * ( Bx * (l1*Bz +   n1*By) + By * (n1*Bx - 2*m1*By) + Bz * (l1*Bx - 2*m1*Bz) );
	double delJcrossBdeln = s0/r1/r1 * ( Bx * (m1*By - 2*n1*Bx) + By * (m1*Bx +   l1*Bz) + Bz * (l1*By - 2*n1*Bz) );
	double delJcrossBdell = s0/r1/r1 * ( Bx * (m1*Bz - 2*l1*Bx) + By * (n1*Bz - 2*l1*By) + Bz * (m1*Bx +   n1*By) );
  
	double delJdotJdelm   =  1.0/s0*( Current_Density[0] * 0  	  - Bz*s0/r1*Current_Density[1] + By*s0/r1*Current_Density[2] );
	double delJdotJdeln   =  1.0/s0*( Current_Density[0] * Bz*s0/r1 - 0*Current_Density[0] 	    - Bx*s0/r1*Current_Density[2] );
	double delJdotJdell   =  1.0/s0*(-Current_Density[0] * By*s0/r1 + Bx*s0/r1*Current_Density[1] + 0*Current_Density[2] );
  
	/* Derivative of the second term in the source terms wrt all the conservative variables */
  
	val_residual[1] = JcrossB[0]*Volume;
	val_residual[2] = JcrossB[1]*Volume;
	if (nDim == 3)	val_residual[3] = JcrossB[2]*Volume;
	val_residual[nDim+1] = (JcrossB[0]*m1/r1 + JcrossB[1]*n1/r1 + JcrossB[2]*l1/r1 + JdotJoversgima)*Volume;
  
	if (implicit) {
		val_Jacobian_i[1][0] = -JcrossB[0]/r1*Volume;
		val_Jacobian_i[1][1] =- s0 * ( Bz*Bz + By*By)/r1*Volume;
		val_Jacobian_i[1][2] = s0*By*Bx/r1*Volume;
		if (nDim ==3) val_Jacobian_i[1][3] =  s0*Bx*Bz/r1*Volume;
		val_Jacobian_i[1][nDim+1] = JcrossB[0]*Volume;
    
		val_Jacobian_i[2][0]  = - JcrossB[1]/r1*Volume;
		val_Jacobian_i[2][1]  = s0*Bx*By/r1*Volume;
		val_Jacobian_i[2][2]  = -s0*(Bx*Bx + Bz*Bz)/r1*Volume;
		if (nDim ==3) val_Jacobian_i[2][3] =  s0*Bz*By/r1*Volume;
		val_Jacobian_i[2][nDim+1] = JcrossB[1]*Volume;
    
		if (nDim == 3) {
			val_Jacobian_i[3][0]  = - JcrossB[2]/r1*Volume;
			val_Jacobian_i[3][1]  = s0*Bx*Bz/r1*Volume;
			val_Jacobian_i[3][2]  = s0*By*Bz/r1*Volume;
			val_Jacobian_i[3][3] = - s0*(By*By + Bx*Bx)/r1*Volume;
			val_Jacobian_i[3][nDim+1] = JcrossB[2]*Volume;
		}
    
		val_Jacobian_i[nDim+1][0]  =  (-2/r1*(JcrossB[0]*m1/r1 + JcrossB[1]*n1/r1 + JcrossB[2]*l1/r1) - 2*JdotJoversgima/r1)*Volume;
		val_Jacobian_i[nDim+1][1]  = (delJcrossBdelm + delJdotJdelm)*Volume;
		val_Jacobian_i[nDim+1][2]  = (delJcrossBdeln + delJdotJdeln)*Volume;
		if (nDim == 3) val_Jacobian_i[nDim+1][3]  = (delJcrossBdell + delJdotJdell)*Volume;
		val_Jacobian_i[nDim+1][nDim+1]  = JdotJoversgima/e1*Volume;
	}
}
