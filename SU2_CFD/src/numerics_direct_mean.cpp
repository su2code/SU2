/*!
 * \file numerics_direct_mean.cpp
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


CUpwRoe_Flow::CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();
  
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

CUpwRoe_Flow::~CUpwRoe_Flow(void) {
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

void CUpwRoe_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	/*--- Face area (norm or the normal vector) ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Unit Normal ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	/*--- Conserved variables at point i,
   Need to recompute SoundSpeed / Pressure / Enthalpy in
   case of 2nd order reconstruction ---*/
	Density_i = U_i[0];
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
	}
	Energy_i     = U_i[nDim+1] / Density_i;
	Pressure_i   = (Gamma_Minus_One*(Energy_i-0.5*sq_vel))*Density_i;
  
  /*--- If it is not a physical solution, then
   use the zero order reconstruction ---*/
  if (((Density_i < 0.0) || (Pressure_i < 0.0)) &&  UZeroOrder_i != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      U_i[iVar] = UZeroOrder_i[iVar];
    
    Density_i = U_i[0];
    sq_vel = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = U_i[iDim+1] / Density_i;
      sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
    }
    Energy_i     = U_i[nDim+1] / Density_i;
    Pressure_i   = (Gamma_Minus_One*(Energy_i-0.5*sq_vel))*Density_i;
  }
  
	SoundSpeed_i = sqrt(Pressure_i*Gamma/Density_i);
	Enthalpy_i   = (U_i[nDim+1] + Pressure_i) / Density_i;
  
	/*--- Conserved variables at point j,
   Need to recompute SoundSpeed / Pressure / Enthalpy in
   case of 2nd order reconstruction ---*/
	Density_j = U_j[0];
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_j     = U_j[nDim+1] / Density_j;
	Pressure_j   = (Gamma_Minus_One*(Energy_j-0.5*sq_vel))*Density_j;
  
  /*--- If it is not a physical solution, then
   use the zero order reconstruction ---*/
  if (((Density_j < 0.0) || (Pressure_j < 0.0)) &&  UZeroOrder_j != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      U_j[iVar] = UZeroOrder_j[iVar];
    
    Density_j = U_j[0];
    sq_vel = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_j[iDim] = U_j[iDim+1] / Density_j;
      sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
    }
    Energy_j     = U_j[nDim+1] / Density_j;
    Pressure_j   = (Gamma_Minus_One*(Energy_j-0.5*sq_vel))*Density_j;
  }
  
	SoundSpeed_j = sqrt(Pressure_j*Gamma/Density_j);
	Enthalpy_j   = (U_j[nDim+1] + Pressure_j) / Density_j;
  
	/*--- Roe-averaged variables at interface between i & j ---*/
	R = sqrt(fabs(Density_j/Density_i));
	RoeDensity = R*Density_i;
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
		sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
	}
	RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
	RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));
  
	/*--- Compute Proj_flux_tensor_i ---*/
	GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, Proj_flux_tensor_i);
  
	/*--- Compute Proj_flux_tensor_j ---*/
	GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, Proj_flux_tensor_j);
  
	/*--- Compute P and Lambda (do it with the Normal) ---*/
	GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, P_Tensor);
  
	ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity   += RoeVelocity[iDim]*UnitaryNormal[iDim];
		ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
	}
  
	/*--- Adjustment for a rotating frame ---*/
	if (rotating_frame) {
		ProjVelocity   -= Rot_Flux/Area;
		ProjVelocity_i -= Rot_Flux/Area;
		ProjVelocity_j -= Rot_Flux/Area;
	}
  
	/*--- Adjustment due to mesh motion (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0; double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitaryNormal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*UnitaryNormal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*UnitaryNormal[iDim];
		}
		ProjVelocity   -= ProjGridVel;
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j -= ProjGridVel;
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
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
			for (jVar = 0; jVar < nVar; jVar++)
				val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
		}
    
		/*--- Flux contribution for a rotating frame ---*/
		if (rotating_frame) {
			ProjVelocity = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			}
		}
    
		/*--- Flux contribution due to grid motion ---*/
		if (grid_movement) {
			ProjVelocity = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			}
		}
	}
	else {
    
		/*--- Compute inverse P ---*/
		GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, invP_Tensor);
    
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
    
		/*--- Jacobian contributions due to a rotating frame ---*/
		if (rotating_frame) {
			ProjVelocity = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
				/*--- Implicit terms ---*/
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
    
		/*--- Jacobian contributions due to grid motion ---*/
		if (grid_movement) {
			ProjVelocity = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
				/*--- Implicit terms ---*/
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
    
	}
  
}

CUpwRoePrim_Flow::CUpwRoePrim_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();
  
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

CUpwRoePrim_Flow::~CUpwRoePrim_Flow(void) {
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

void CUpwRoePrim_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	/*--- Face area (norm or the normal vector) ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Unit Normal ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	/*--- Primitive variables at point i ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity_i[iDim] = V_i[iDim+1];
	Pressure_i = V_i[nDim+1];
	Density_i = V_i[nDim+2];
	Enthalpy_i = V_i[nDim+3];
  
	/*--- Conserved variables at point j,
	 Need to recompute SoundSpeed / Pressure / Enthalpy in
	 case of 2nd order reconstruction ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity_j[iDim] = V_j[iDim+1];
	Pressure_j = V_j[nDim+1];
	Density_j = V_j[nDim+2];
	Enthalpy_j = V_j[nDim+3];
  
	/*--- Roe-averaged variables at interface between i & j ---*/
	R = sqrt(fabs(Density_j/Density_i));
	RoeDensity = R*Density_i;
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
		sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
	}
	RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
	RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));
  
	/*--- Compute Proj_flux_tensor_i ---*/
	GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, Proj_flux_tensor_i);
  
	/*--- Compute Proj_flux_tensor_j ---*/
	GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, Proj_flux_tensor_j);
  
	/*--- Compute P and Lambda (do it with the Normal) ---*/
	GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, P_Tensor);
  
	ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity   += RoeVelocity[iDim]*UnitaryNormal[iDim];
		ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
	}
  
	/*--- Adjustment for a rotating frame ---*/
	if (rotating_frame) {
		ProjVelocity   -= Rot_Flux/Area;
		ProjVelocity_i -= Rot_Flux/Area;
		ProjVelocity_j -= Rot_Flux/Area;
	}
  
	/*--- Adjustment due to mesh motion (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0; double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitaryNormal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*UnitaryNormal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*UnitaryNormal[iDim];
		}
		ProjVelocity   -= ProjGridVel;
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j -= ProjGridVel;
	}
  
	/*--- Flow eigenvalues and entropy correctors ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		Lambda[iDim] = ProjVelocity;
	Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
	Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
	/*--- Harten and Hyman entropy correction ---*/
	SoundSpeed_i = sqrt(Pressure_i*Gamma/Density_i); SoundSpeed_j = sqrt(Pressure_j*Gamma/Density_j);
	for (iDim = 0; iDim < nDim; iDim++)
		Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,ProjVelocity_j-Lambda[iDim]));
	Epsilon[nVar-2] = 4.0*max(0.0, max(Lambda[nVar-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar-2]));
	Epsilon[nVar-1] = 4.0*max(0.0, max(Lambda[nVar-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar-1]));
  
	for (iVar = 0; iVar < nVar; iVar++)
		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
		else
			Lambda[iVar] = fabs(Lambda[iVar]);
  
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
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
		for (jVar = 0; jVar < nVar; jVar++)
			val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
	}
  
	/*--- Flux contribution for a rotating frame ---*/
	if (rotating_frame) {
		ProjVelocity = Rot_Flux;
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
		}
	}
  
	/*--- Flux contribution due to grid motion ---*/
	if (grid_movement) {
		ProjVelocity = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
		}
	}
  
	if (implicit) {
    
		/*--- Compute inverse P ---*/
		GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, invP_Tensor);
    
		/*--- Jacobians of the inviscid flux ---*/
		Energy_i = Enthalpy_i - Pressure_i/Density_i; Energy_j = Enthalpy_j - Pressure_j/Density_j;
		GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				Proj_ModJac_Tensor_ij = 0.0;
				for (kVar = 0; kVar < nVar; kVar++)
					Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
				val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
				val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
			}
		}
    
		/*--- Jacobian contributions due to a rotating frame ---*/
		if (rotating_frame) {
			ProjVelocity = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
    
		/*--- Jacobian contributions due to grid motion ---*/
		if (grid_movement) {
			ProjVelocity = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
    
	}
  
}

CUpwRoe_Turkel_Flow::CUpwRoe_Turkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Beta_min = config->GetminTurkelBeta();
	Beta_max = config->GetmaxTurkelBeta();
  
	Diff_U = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	RoeVelocity = new double [nDim];
	Proj_flux_tensor_i = new double [nVar];
	Proj_flux_tensor_j = new double [nVar];
	Lambda = new double [nVar];
	Epsilon = new double [nVar];
	absPeJac = new double* [nVar];
	invRinvPe = new double* [nVar];
	R_Tensor  = new double* [nVar];
	Matrix    = new double* [nVar];
	Art_Visc  = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		absPeJac[iVar] = new double [nVar];
		invRinvPe[iVar] = new double [nVar];
		Matrix[iVar] = new double [nVar];
		Art_Visc[iVar] = new double [nVar];
		R_Tensor[iVar] = new double [nVar];
	}
}

CUpwRoe_Turkel_Flow::~CUpwRoe_Turkel_Flow(void) {
	unsigned short iVar;
  
	delete [] Diff_U;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] RoeVelocity;
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
	delete [] Lambda;
	delete [] Epsilon;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] absPeJac[iVar];
		delete [] invRinvPe[iVar];
		delete [] Matrix[iVar];
		delete [] Art_Visc[iVar];
		delete [] R_Tensor[iVar];
	}
	delete [] Matrix;
	delete [] Art_Visc;
	delete [] absPeJac;
	delete [] invRinvPe;
	delete [] R_Tensor;
  
  
}

void CUpwRoe_Turkel_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	double r_hat, s_hat, t_hat, rhoB2a2, sqr_one_m_Betasqr_Lam1;
	double Beta2, one_m_Betasqr, one_p_Betasqr, sqr_two_Beta_c_Area;
	double local_Mach;
  
	/*--- Face area (norm or the normal vector) ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Unit Normal ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	/*--- Conserved variables at point i,
   Need to recompute SoundSpeed / Pressure / Enthalpy in
   case of 2nd order reconstruction ---*/
	Density_i = U_i[0];
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
	}
	Energy_i     = U_i[nDim+1] / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel));
	Pressure_i   = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	Enthalpy_i   = (U_i[nDim+1] + Pressure_i) / Density_i;
  
	/*--- Conserved variables at point j,
   Need to recompute SoundSpeed / Pressure / Enthalpy in
   case of 2nd order reconstruction ---*/
	Density_j = U_j[0];
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_j     = U_j[nDim+1] / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel));
	Pressure_j   = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	Enthalpy_j   = (U_j[nDim+1] + Pressure_j) / Density_j;
  
	/*--- Roe-averaged variables at interface between i & j ---*/
	R = sqrt(Density_j/Density_i);
	RoeDensity = R*Density_i;
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
		sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
	}
	RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
	RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));
	RoePressure = RoeDensity/Gamma*RoeSoundSpeed*RoeSoundSpeed;
  
	/*--- Compute Proj_flux_tensor_i ---*/
	GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, Proj_flux_tensor_i);
  
	/*--- Compute Proj_flux_tensor_j ---*/
	GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, Proj_flux_tensor_j);
  
	ProjVelocity = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVelocity   += RoeVelocity[iDim]*Normal[iDim];
  
  /*--- Adjustment for a rotating frame ---*/
  if (rotating_frame) {
		ProjVelocity -= Rot_Flux/Area;
	}
  
	/*--- First few flow eigenvalues of A.Normal with the normal---*/
	for (iDim = 0; iDim < nDim; iDim++)
		Lambda[iDim] = ProjVelocity;
  
	local_Mach = sqrt(sq_vel)/RoeSoundSpeed;
	Beta 	   = max(Beta_min,min(local_Mach,Beta_max));
	Beta2 	   = Beta*Beta;
  
	one_m_Betasqr 		   = 1.0 - Beta2;  // 1-Beta*Beta
	one_p_Betasqr 		   = 1.0 + Beta2;  // 1+Beta*Beta
	sqr_one_m_Betasqr_Lam1 = pow((one_m_Betasqr*Lambda[0]),2); // [(1-Beta^2)*Lambda[0]]^2
	sqr_two_Beta_c_Area    = pow(2.0*Beta*RoeSoundSpeed*Area,2); // [2*Beta*c*Area]^2
  
	/*--- The rest of the flow eigenvalues of preconditioned matrix---*/
	Lambda[nVar-2] = 0.5 * ( one_p_Betasqr*Lambda[0] + sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
	Lambda[nVar-1] = 0.5 * ( one_p_Betasqr*Lambda[0] - sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
  
	s_hat = 1.0/Area * (Lambda[nVar-1] - Lambda[0]*Beta2);
	r_hat = 1.0/Area * (Lambda[nVar-2] - Lambda[0]*Beta2);
	t_hat = 0.5/Area * (Lambda[nVar-1] - Lambda[nVar-2]);
	rhoB2a2 = RoeDensity*Beta2*RoeSoundSpeed*RoeSoundSpeed;
  
	/*--- Diference variables iPoint and jPoint and absolute value of the eigen values---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_U[iVar] = U_j[iVar]-U_i[iVar];
		Lambda[iVar] = fabs(Lambda[iVar]);
	}
  
	/*--- Compute the absolute Preconditioned Jacobian in entropic Variables (do it with the Unitary Normal) ---*/
	GetPrecondJacobian(Beta2, r_hat, s_hat, t_hat, rhoB2a2, Lambda, UnitaryNormal, absPeJac);
  
	/*--- Compute the matrix from entropic variables to conserved variables ---*/
	GetinvRinvPe(Beta2, RoeEnthalpy, RoeSoundSpeed, RoeDensity, RoeVelocity, invRinvPe);
  
	/*--- Compute the matrix from entropic variables to conserved variables ---*/
	GetRMatrix(RoePressure, RoeSoundSpeed, RoeDensity, RoeVelocity, R_Tensor);
  
	if (implicit) {
		/*--- Jacobians of the inviscid flux, scaled by
     0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
		GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
	}
  
	for (iVar = 0; iVar < nVar; iVar ++){
		for (jVar = 0; jVar < nVar; jVar ++) {
			Matrix[iVar][jVar] = 0.0;
			for (kVar = 0; kVar < nVar; kVar++)
				Matrix[iVar][jVar]  += absPeJac[iVar][kVar]*R_Tensor[kVar][jVar];
		}
	}
  
	for (iVar = 0; iVar < nVar; iVar ++){
		for (jVar = 0; jVar < nVar; jVar ++) {
			Art_Visc[iVar][jVar] = 0.0;
			for (kVar = 0; kVar < nVar; kVar++)
				Art_Visc[iVar][jVar]  += invRinvPe[iVar][kVar]*Matrix[kVar][jVar];
		}
	}
  
	/*--- Roe's Flux approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
		for (jVar = 0; jVar < nVar; jVar++) {
			val_residual[iVar] -= 0.5*Art_Visc[iVar][jVar]*Diff_U[jVar];
			if (implicit) {
				val_Jacobian_i[iVar][jVar] += 0.5*Art_Visc[iVar][jVar];
				val_Jacobian_j[iVar][jVar] -= 0.5*Art_Visc[iVar][jVar];
			}
		}
	}
  
  /*--- Contributions due to a rotating frame ---*/
  if (rotating_frame) {
    ProjVelocity = Rot_Flux;
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      /*--- Implicit terms ---*/
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }
  
}


CUpwRoeArtComp_Flow::CUpwRoeArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	gravity = config->GetGravityForce();
	Froude = config->GetFroude();
  
	Diff_U = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	MeanVelocity = new double [nDim];
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

CUpwRoeArtComp_Flow::~CUpwRoeArtComp_Flow(void) {
	unsigned short iVar;
  
	delete [] Diff_U;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] MeanVelocity;
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

void CUpwRoeArtComp_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	/*--- Compute face area ---*/
	Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
  /*--- Compute and unitary normal vector ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		UnitaryNormal[iDim] = Normal[iDim]/Area;
    if (fabs(UnitaryNormal[iDim]) < EPS) UnitaryNormal[iDim] = EPS;
  }
  
	/*--- Set velocity and pressure variables at points iPoint and jPoint ---*/
	Pressure_i = U_i[0]; Pressure_j = U_j[0];
  
  ProjVelocity = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
    Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
		ProjVelocity += MeanVelocity[iDim]*Normal[iDim];
  }
  
	/*--- Mean variables at points iPoint and jPoint ---*/
	MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
	MeanPressure = 0.5*(Pressure_i + Pressure_j);
	MeanBetaInc2 = 0.5*(BetaInc2_i + BetaInc2_j);
	MeanSoundSpeed = sqrt(ProjVelocity*ProjVelocity + (MeanBetaInc2/MeanDensity) * Area * Area);
  
	/*--- Compute Proj_flux_tensor_i ---*/
	GetInviscidArtCompProjFlux(&DensityInc_i, Velocity_i, &Pressure_i, &BetaInc2_i, Normal, Proj_flux_tensor_i);
  
	/*--- Compute Proj_flux_tensor_j ---*/
	GetInviscidArtCompProjFlux(&DensityInc_j, Velocity_j, &Pressure_j, &BetaInc2_j, Normal, Proj_flux_tensor_j);
  
	/*--- Compute P and Lambda (matrix of eigenvalues) ---*/
	GetPArtCompMatrix(&MeanDensity, MeanVelocity, &MeanBetaInc2, UnitaryNormal, P_Tensor);
  
	/*--- Flow eigenvalues ---*/
	if (nDim == 2) {
		Lambda[0] = ProjVelocity;
		Lambda[1] = ProjVelocity + MeanSoundSpeed;
		Lambda[2] = ProjVelocity - MeanSoundSpeed;
	}
	if (nDim == 3) {
		Lambda[0] = ProjVelocity;
		Lambda[1] = ProjVelocity;
		Lambda[2] = ProjVelocity + MeanSoundSpeed;
		Lambda[3] = ProjVelocity - MeanSoundSpeed;
	}
  
  /*--- Absolute value of the eigenvalues ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Lambda[iVar] = fabs(Lambda[iVar]);
  
	/*--- Compute inverse P ---*/
	GetPArtCompMatrix_inv(&MeanDensity, MeanVelocity, &MeanBetaInc2, UnitaryNormal, invP_Tensor);
  
	if (implicit) {
		/*--- Jacobian of the inviscid flux ---*/
		GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, Normal, 0.5, val_Jacobian_j);
	}
  
	/*--- Diference variables iPoint and jPoint ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Diff_U[iVar] = U_j[iVar] - U_i[iVar];
  
	/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar];
			if (implicit) {
				val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij;
				val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij;
			}
		}
	}
  
}

CUpwAUSM_Flow::CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
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

CUpwAUSM_Flow::~CUpwAUSM_Flow(void) {
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

void CUpwAUSM_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
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

CUpwHLLC_Flow::CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
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

CUpwHLLC_Flow::~CUpwHLLC_Flow(void) {
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

void CUpwHLLC_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
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

CCentJST_Flow::CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	grid_movement = config->GetGrid_Movement();
	rotating_frame = config->GetRotating_Frame();
  
	/*--- Artifical dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_Flow();
	Param_Kappa_4 = config->GetKappa_4th_Flow();
  
	/*--- Allocate some structures ---*/
	Diff_U = new double [nVar];
	Diff_Lapl = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	MeanVelocity = new double [nDim];
	Proj_flux_tensor = new double [nVar];
  
}

CCentJST_Flow::~CCentJST_Flow(void) {
	delete [] Diff_U;
	delete [] Diff_Lapl;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] MeanVelocity;
	delete [] Proj_flux_tensor;
}

void CCentJST_Flow::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
                                CConfig *config) {
  
	/*--- Density, velocity and energy at points i and j ---*/
	Density_i = U_i[0]; Density_j = U_j[0];
  
	sq_vel_i = 0.0; sq_vel_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_i = U_i[nDim+1] / Density_i;
	Energy_j = U_j[nDim+1] / Density_j;
  
	/*--- Mean Values ---*/
	MeanDensity = 0.5*(Density_i+Density_j);
	MeanPressure = 0.5*(Pressure_i+Pressure_j);
	MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
	MeanEnergy = 0.5*(Energy_i+Energy_j);
  
	/*--- Get projected flux tensor ---*/
	GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, Proj_flux_tensor);
  
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}
  
	/*--- Adjustment for a rotating frame ---*/
	if (rotating_frame) {
		ProjVelocity = Rot_Flux;
		for (iVar = 0; iVar < nVar; iVar++) {
			val_resconv[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			if (implicit) {
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
	}
  
	/*--- Adjustment due to mesh motion ---*/
	if (grid_movement) {
		ProjVelocity = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		for (iVar = 0; iVar < nVar; iVar++) {
			val_resconv[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			if (implicit) {
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
	}
  
	/*--- Computes differences btw. Laplacians and conservative variables,
   with a correction for the enthalpy ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_U[iVar] = U_i[iVar]-U_j[iVar];
	}
	Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
	/*--- Compute the local espectral radius and the stretching factor ---*/
	ProjVelocity_i = 0; ProjVelocity_j = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
  
  /*--- Adjustment for a rotating frame ---*/
	if (rotating_frame) {
		ProjVelocity_i -= Rot_Flux;
		ProjVelocity_j -= Rot_Flux;
	}
  
	/*--- Adjustment due to mesh motion ---*/
	if (grid_movement) {
		ProjGridVel_i = 0.0; ProjGridVel_j = 0.0;  ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*Normal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*Normal[iDim];
		}
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j -= ProjGridVel;
	}
  
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
	
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
  
	Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
	Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
	/*--- Compute viscous part of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_resvisc[iVar] = (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/
	if (implicit) {
		cte_0 = (Epsilon_2 + Epsilon_4*double(Neighbor_i+1))*StretchingFactor*MeanLambda;
		cte_1 = (Epsilon_2 + Epsilon_4*double(Neighbor_j+1))*StretchingFactor*MeanLambda;
    
		for (iVar = 0; iVar < (nVar-1); iVar++) {
			val_Jacobian_i[iVar][iVar] += cte_0;
			val_Jacobian_j[iVar][iVar] -= cte_1;
		}
    
		/*--- Last row of Jacobian_i ---*/
		val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*sq_vel_i;
		for (iDim = 0; iDim < nDim; iDim++)
			val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
		val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;
    
		/*--- Last row of Jacobian_j ---*/
		val_Jacobian_j[nVar-1][0] -= cte_1*Gamma_Minus_One*sq_vel_j;
		for (iDim = 0; iDim < nDim; iDim++)
			val_Jacobian_j[nVar-1][iDim+1] += cte_1*Gamma_Minus_One*Velocity_j[iDim];
		val_Jacobian_j[nVar-1][nVar-1] -= cte_1*Gamma;
	}
}

CCentJSTArtComp_Flow::CCentJSTArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	grid_movement = config->GetGrid_Movement();
	rotating_frame = config->GetRotating_Frame();
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	gravity = config->GetGravityForce();
	Froude = config->GetFroude();
  
	/*--- Artifical dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_Flow();
	Param_Kappa_4 = config->GetKappa_4th_Flow();
  
	/*--- Allocate some structures ---*/
	Diff_U = new double [nVar];
	Diff_Lapl = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	MeanVelocity = new double [nDim];
	Proj_flux_tensor = new double [nVar];
  
}

CCentJSTArtComp_Flow::~CCentJSTArtComp_Flow(void) {
	delete [] Diff_U;
	delete [] Diff_Lapl;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] MeanVelocity;
	delete [] Proj_flux_tensor;
}

void CCentJSTArtComp_Flow::ComputeResidual(double *val_resconv, double *val_resvisc,
                                       double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	/*--- Conservative variables at point i and j ---*/
	Pressure_i = U_i[0]; Pressure_j = U_j[0];
	sq_vel_i = 0.0; sq_vel_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
		sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
	}
  
	/*--- Compute mean values of the variables ---*/
	MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
	MeanPressure = 0.5*(Pressure_i + Pressure_j);
	MeanBetaInc2 = 0.5*(BetaInc2_i + BetaInc2_j);
  
	/*--- Get projected flux tensor ---*/
	GetInviscidArtCompProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, Normal, Proj_flux_tensor);
  
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		GetInviscidArtCompProjJac(&MeanDensity, MeanVelocity, &MeanBetaInc2, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}
  
	/*--- Computes differences between Laplacians and conservative variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_U[iVar] = U_i[iVar]-U_j[iVar];
	}
  
	/*--- Compute the local espectral radius and the stretching factor ---*/
	ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
  
	SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
	SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);
  
	Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
	Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
  
	Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
	Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
	/*--- Compute viscous part of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_resvisc[iVar] = (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
  
	if (implicit) {
		cte_0 = (Epsilon_2 + Epsilon_4*double(Neighbor_i+1))*StretchingFactor*MeanLambda;
		cte_1 = (Epsilon_2 + Epsilon_4*double(Neighbor_j+1))*StretchingFactor*MeanLambda;
    
		for (iVar = 0; iVar < nVar; iVar++) {
			val_Jacobian_i[iVar][iVar] += cte_0;
			val_Jacobian_j[iVar][iVar] -= cte_1;
		}
	}
}

CCentLax_Flow::CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	grid_movement = config->GetGrid_Movement();
	rotating_frame = config->GetRotating_Frame();
  
	/*--- Artifical dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_Flow();
  
	/*--- Allocate some structures ---*/
	Diff_U = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	MeanVelocity = new double [nDim];
	Proj_flux_tensor = new double [nVar];
  
}

CCentLax_Flow::~CCentLax_Flow(void) {
	delete [] Diff_U;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] MeanVelocity;
	delete [] Proj_flux_tensor;
  
}

void CCentLax_Flow::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
                                CConfig *config) {
  
	/*--- Evaluate points 0 and 1 ---*/
	Density_i = U_i[0]; Density_j = U_j[0];
	sq_vel_i = 0.0; sq_vel_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_i = U_i[nVar-1]/Density_i;
	Energy_j = U_j[nVar-1]/Density_j;
  
	/*--- Compute mean values of the variables ---*/
	MeanDensity = 0.5*(Density_i+Density_j);
	MeanPressure = 0.5*(Pressure_i+Pressure_j);
	MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
	MeanEnergy = 0.5*(Energy_i+Energy_j);
  
	/*--- Get projected flux tensor ---*/
	GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, Proj_flux_tensor);
  
	/*--- Compute inviscid residual ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}
  
	/*--- Jacobians of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	if (implicit) {
		GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}
  
	/*--- Adjustment for a rotating frame ---*/
	if (rotating_frame) {
		ProjVelocity = Rot_Flux;
		for (iVar = 0; iVar < nVar; iVar++) {
			val_resconv[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			if (implicit) {
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
	}
  
	/*--- Adjustment due to grid motion ---*/
	if (grid_movement) {
		ProjVelocity = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		for (iVar = 0; iVar < nVar; iVar++) {
			val_resconv[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			if (implicit) {
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
	}
  
	/*--- Computes differences btw. conservative variables ---*/
	for (iVar = 0; iVar < nDim+1; iVar++)
		Diff_U[iVar] = U_i[iVar]-U_j[iVar];
	Diff_U[nDim+1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
	/*--- Compute the local espectral radius and the stretching factor ---*/
	ProjVelocity_i = 0; ProjVelocity_j = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
  
	if (rotating_frame) {
		ProjVelocity_i -= Rot_Flux;
		ProjVelocity_j -= Rot_Flux;
	}
  
	/*--- TDE ---*/
	if (grid_movement) {
		double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0; double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*Normal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*Normal[iDim];
		}
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j -= ProjGridVel;
	}
  
	Area = sqrt(Area);
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	/*	cout << "EULER LocalLambda_i " << Local_Lambda_i << endl;
   cout << "EULER LocalLambda_j " << Local_Lambda_j << endl;
   cout << "SoundSpeed_i " << SoundSpeed_i << endl;
   cout << "Area " << Area << endl;*/
  
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	sc0 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_0 = Param_Kappa_0*sc0*double(nDim)/3.0;
  
	/*	cout << "Euler sc0: " << sc0 << endl;
   cout << "Epsilon_0: " << Epsilon_0 << endl;*/
  
	/*--- Compute viscous part of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resvisc[iVar] = Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
		/*		cout << "artificial dissipation: " << (Epsilon_0*Diff_U[iVar])*StretchingFactor*MeanLambda << endl;
     cout << "epsilon0: " << Epsilon_0 << endl;
     cout << "DiffU: " << Diff_U[iVar] << endl;
     cout << "stretching: " << StretchingFactor << endl;
     cout << "MeanLambda: " << MeanLambda << endl;
     cin.get();*/
		//		cout << val_resconv[iVar] << endl;
	}
	/*	cout << "EULER val_resvisc:" << endl;
   for (iVar = 0; iVar < nVar; iVar++)
   cout << val_resvisc[iVar] << endl;
   cin.get();*/
  
	if (implicit) {
		cte = Epsilon_0*StretchingFactor*MeanLambda;
		//		cout << "cte " << cte << endl;
    
		for (iVar = 0; iVar < (nVar-1); iVar++) {
			val_Jacobian_i[iVar][iVar] += cte;
			val_Jacobian_j[iVar][iVar] -= cte;
		}
    
		/*--- Last rows: CAREFUL!! You have differences of \rho_Enthalpy, not differences of \rho_Energy ---*/
		val_Jacobian_i[nVar-1][0] += cte*Gamma_Minus_One*sq_vel_i;
		for (iDim = 0; iDim < nDim; iDim++)
			val_Jacobian_i[nVar-1][iDim+1] -= cte*Gamma_Minus_One*Velocity_i[iDim];
		val_Jacobian_i[nVar-1][nVar-1] += cte*Gamma;
    
		/*--- Last row of Jacobian_j ---*/
		val_Jacobian_j[nVar-1][0] -= cte*Gamma_Minus_One*sq_vel_j;
		for (iDim = 0; iDim < nDim; iDim++)
			val_Jacobian_j[nVar-1][iDim+1] += cte*Gamma_Minus_One*Velocity_j[iDim];
		val_Jacobian_j[nVar-1][nVar-1] -= cte*Gamma;
	}
	/*	cout << "EULER val_Jacobian_i: " << endl;
   for (iVar = 0; iVar < nVar; iVar++) {
   for (jVar =0; jVar < nVar; jVar++) {
   cout << "\t" << val_Jacobian_i[iVar][jVar];
   }
   cout << endl;
   }
   cin.get();*/
}

CCentLaxArtComp_Flow::CCentLaxArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	grid_movement = config->GetGrid_Movement();
	rotating_frame = config->GetRotating_Frame();
	gravity = config->GetGravityForce();
	Froude = config->GetFroude();
  
	/*--- Artificial dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_Flow();
  
	/*--- Allocate some structures ---*/
	Diff_U = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	MeanVelocity = new double [nDim];
	Proj_flux_tensor = new double [nVar];
  
}

CCentLaxArtComp_Flow::~CCentLaxArtComp_Flow(void) {
	delete [] Diff_U;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] MeanVelocity;
	delete [] Proj_flux_tensor;
  
}

void CCentLaxArtComp_Flow::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
                                       CConfig *config) {
  
	/*--- Conservative variables at point i and j ---*/
	Pressure_i = U_i[0]; Pressure_j = U_j[0];
	sq_vel_i = 0.0; sq_vel_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
		sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
  
	/*--- Correction ---*/
	//	double epsilon = - STANDART_GRAVITY * ((0.5*DensityInc_i*(Coord_i[nDim-1]-Coord_j[nDim-1])
	//																					+ 0.5*DensityInc_j*(Coord_j[nDim-1]-Coord_i[nDim-1]))
	//																				 /(2.0*config->GetFroude()*config->GetFroude()));
  
	/*--- Compute mean values of the variables ---*/
	MeanDensity = 0.5*(DensityInc_i+DensityInc_j);
	MeanPressure = 0.5*(Pressure_i+Pressure_j);
	MeanBetaInc2 = 0.5*(BetaInc2_i+BetaInc2_j);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  
	/*--- Get projected flux tensor ---*/
	GetInviscidArtCompProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, Normal, Proj_flux_tensor);
  
	/*--- Compute inviscid residual ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		GetInviscidArtCompProjJac(&MeanDensity, MeanVelocity, &MeanBetaInc2, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}
  
	/*--- Computes differences btw. conservative variables ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  
	/*--- Compute the local espectral radius and the stretching factor ---*/
	ProjVelocity_i = 0; ProjVelocity_j = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
  
	Area = sqrt(Area);
  
	SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
	SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);
  
	Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
	Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;
	MeanLambda = 0.5*(Local_Lambda_i + Local_Lambda_j);
  
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	sc0 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_0 = Param_Kappa_0*sc0*double(nDim)/3.0;
  
	/*--- Compute viscous part of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_resvisc[iVar] = Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			val_Jacobian_i[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_j[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
		}
	}
  
}

CAvgGrad_Flow::CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Compressible flow, primitive variables nDim+3, (T,vx,vy,vz,P,rho) ---*/
	PrimVar_i = new double [nDim+3];
	PrimVar_j = new double [nDim+3];
	Mean_PrimVar = new double [nDim+3];
  
  /*--- Compressible flow, primitive gradient variables nDim+3, (T,vx,vy,vz) ---*/
	Mean_GradPrimVar = new double* [nDim+1];
	for (iVar = 0; iVar < nDim+1; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGrad_Flow::~CAvgGrad_Flow(void) {
  
	delete [] PrimVar_i;
	delete [] PrimVar_j;
	delete [] Mean_PrimVar;
  
	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGrad_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CAvgGrad_Flow::ComputeResidual
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i *U_j **PrimVar_Grad_i **PrimVar_Grad_j
	//SU2_CPP2C INVARS Laminar_Viscosity_i Laminar_Viscosity_j Eddy_Viscosity_i Eddy_Viscosity_j
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *Normal Gamma Gamma_Minus_One Gas_Constant
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim nVar
  
	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS INT SCALAR iDim jDim iVar
	//SU2_CPP2C VARS DOUBLE SCALAR Area Mean_Laminar_Viscosity Mean_Eddy_Viscosity dist_ij
	//SU2_CPP2C VARS DOUBLE SCALAR total_viscosity heat_flux_factor div_vel cp Density sq_vel
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nVar PrimVar_i PrimVar_j Mean_PrimVar Proj_Flux_Tensor
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nVar SIZE=nDim Mean_GradPrimVar
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nDim SIZE=nDim tau
	//SU2_CPP2C DECL_LIST END
  
  
	/*--- Normalized normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	for (iVar = 0; iVar < nDim+3; iVar++) {
		PrimVar_i[iVar] = V_i[iVar];
		PrimVar_j[iVar] = V_j[iVar];
		Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
	}
  
	/*--- Mean Viscosities and turbulent kinetic energy---*/
	Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
	Mean_Eddy_Viscosity    = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
	Mean_turb_ke           = 0.5*(turb_ke_i + turb_ke_j);
  
	/*--- Mean gradient approximation ---*/
	for (iVar = 0; iVar < nDim+1; iVar++) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
		}
	}
  
	/*--- Get projected flux tensor ---*/
	//SU2_CPP2C SUBROUTINE START GetViscousProjFlux
	GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
	//SU2_CPP2C SUBROUTINE LOCATION numerics_structure.cpp
	//SU2_CPP2C SUBROUTINE VARS Mean_PrimVar Mean_GradPrimVar Normal Mean_Laminar_Viscosity Mean_Eddy_Viscosity
	//SU2_CPP2C SUBROUTINE END
  
	/*--- Update viscous residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
	//SU2_CPP2C COMMENT START
  
	/*--- Compute the implicit part ---*/
	if (implicit) {
		dist_ij = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
		dist_ij = sqrt(dist_ij);
    
		GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                       dist_ij, UnitaryNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
	}
  
	//SU2_CPP2C COMMENT END
  
	//SU2_CPP2C END CAvgGrad_Flow::ComputeResidual
}

CAvgGradArtComp_Flow::CAvgGradArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Incompressible flow, primitive variables nDim+1, (P,vx,vy,vz) ---*/
	PrimVar_i = new double [nVar];
	PrimVar_j = new double [nVar];
	Mean_PrimVar = new double [nVar];
	Mean_GradPrimVar = new double* [nVar];
  
  /*--- Incompressible flow, gradient primitive variables nDim+1, (P,vx,vy,vz) ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGradArtComp_Flow::~CAvgGradArtComp_Flow(void) {
  
	delete [] PrimVar_i;
	delete [] PrimVar_j;
	delete [] Mean_PrimVar;
  
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGradArtComp_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i,
                                       double **val_Jacobian_j, CConfig *config) {
  
	/*--- Normalized normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	/*--- Conversion to Primitive Variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		PrimVar_i[iVar] = V_i[iVar];
		PrimVar_j[iVar] = V_j[iVar];
		Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
	}
  
	/*--- Mean Viscosities ---*/
	Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
	Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  
	/*--- Mean gradient approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
  
	/*--- Get projected flux tensor ---*/
	GetViscousArtCompProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
	/*--- Update viscous residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
	/*--- Implicit part ---*/
	if (implicit) {
		dist_ij = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
		dist_ij = sqrt(dist_ij);
    
		GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitaryNormal,
                              Area, val_Jacobian_i, val_Jacobian_j);
	}
}

CAvgGradCorrected_Flow::CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Compressible flow, primitive variables nDim+3, (T,vx,vy,vz,P,rho) ---*/
	PrimVar_i = new double [nDim+3];
	PrimVar_j = new double [nDim+3];
	Mean_PrimVar = new double [nDim+3];
  
  /*--- Compressible flow, primitive gradient variables nDim+1, (T,vx,vy,vz) ---*/
	Proj_Mean_GradPrimVar_Edge = new double [nDim+1];
	Mean_GradPrimVar = new double* [nDim+1];
	for (iVar = 0; iVar < nDim+1; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
  
  Edge_Vector = new double [nDim];
  
}

CAvgGradCorrected_Flow::~CAvgGradCorrected_Flow(void) {
  
	delete [] PrimVar_i;
	delete [] PrimVar_j;
	delete [] Mean_PrimVar;
	delete [] Proj_Mean_GradPrimVar_Edge;
	delete [] Edge_Vector;
  
	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGradCorrected_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
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
  
	for (iVar = 0; iVar < nDim+3; iVar++) {
		PrimVar_i[iVar] = V_i[iVar];
		PrimVar_j[iVar] = V_j[iVar];
		Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
	}
  
	/*--- Mean Viscosities and turbulent kinetic energy ---*/
	Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
	Mean_Eddy_Viscosity    = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
	Mean_turb_ke           = 0.5*(turb_ke_i + turb_ke_j);
  
	/*--- Projection of the mean gradient in the direction of the edge ---*/
	for (iVar = 0; iVar < nDim+1; iVar++) {
		Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
		}
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                       (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
		}
	}
  
	/*--- Get projected flux tensor ---*/
	GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
	/*--- Save residual value ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
	/*--- Compute the implicit part ---*/
	if (implicit) {
		GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                       sqrt(dist_ij_2), UnitaryNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
	}
}

CAvgGradCorrectedArtComp_Flow::CAvgGradCorrectedArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
	PrimVar_i = new double [nVar];
	PrimVar_j = new double [nVar];
	Mean_PrimVar = new double [nVar];
	Proj_Mean_GradPrimVar_Edge = new double [nVar];
	Edge_Vector = new double [nDim];
  
	Mean_GradPrimVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGradCorrectedArtComp_Flow::~CAvgGradCorrectedArtComp_Flow(void) {
  
	delete [] PrimVar_i;
	delete [] PrimVar_j;
	delete [] Mean_PrimVar;
	delete [] Proj_Mean_GradPrimVar_Edge;
	delete [] Edge_Vector;
  
	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGradCorrectedArtComp_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	/*--- Normalized normal vector ---*/
	Area = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;
  
	/*--- Conversion to Primitive Variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		PrimVar_i[iVar] = V_i[iVar];
		PrimVar_j[iVar] = V_j[iVar];
		Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
	}
  
	/*--- Mean Viscosities ---*/
	Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
	Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  
	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
	}
  
	/*--- Projection of the mean gradient in the direction of the edge ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
		}
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                       (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
		}
	}
  
	/*--- Get projected flux tensor ---*/
	GetViscousArtCompProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
	/*--- Update viscous residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
	/*--- Implicit part ---*/
	if (implicit) {
		GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, sqrt(dist_ij_2), UnitaryNormal,
                              Area, val_Jacobian_i, val_Jacobian_j);
	}
}

CSourcePieceWise_Gravity::CSourcePieceWise_Gravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	incompressible = (config->GetIncompressible() == YES);
  
}

CSourcePieceWise_Gravity::~CSourcePieceWise_Gravity(void) { }

void CSourcePieceWise_Gravity::ComputeResidual(double *val_residual, CConfig *config) {
	unsigned short iVar;
  
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = 0.0;
  
	if (incompressible) {
    
		/*--- Compute the Froude number  ---*/
		Froude = config->GetFroude();
    
		/*--- Evaluate the source term  ---*/
		val_residual[nDim] = Volume * DensityInc_i / (Froude * Froude);
    
	}
  else {
    
		/*--- Evaluate the source term  ---*/
		val_residual[nDim] = Volume * U_i[0] * STANDART_GRAVITY;
    
	}
  
}

CSourceRotatingFrame_Flow::CSourceRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
}

CSourceRotatingFrame_Flow::~CSourceRotatingFrame_Flow(void) { }

void CSourceRotatingFrame_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
  
	double CrossProduct[3], vel[3] = {0,0,0};
  
	/*--- Retrieve the angular velocity vector ---*/
	double *Omega = config->GetOmega_FreeStreamND();
  
	/*--- Calculate momentum source terms as: rho * ( Omega X V ) ---*/
	for(unsigned short iDim = 0; iDim < nDim; iDim++)
		vel[iDim] = U_i[iDim+1];
  
	CrossProduct[0] = Omega[1]*vel[2] - Omega[2]*vel[1];
	CrossProduct[1] = Omega[2]*vel[0] - Omega[0]*vel[2];
	CrossProduct[2] = Omega[0]*vel[1] - Omega[1]*vel[0];
  
	if (nDim == 2) {
		val_residual[0] = 0.0;
		val_residual[1] = CrossProduct[0]*Volume;
		val_residual[2] = CrossProduct[1]*Volume;
		val_residual[3] = 0.0;
	}
  
	if (nDim == 3) {
		val_residual[0] = 0.0;
		val_residual[1] = CrossProduct[0]*Volume;
		val_residual[2] = CrossProduct[1]*Volume;
		val_residual[3] = CrossProduct[2]*Volume;
		val_residual[4] = 0.0;
	}
  
}

CSourceAxisymmetric_Flow::CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
}

CSourceAxisymmetric_Flow::~CSourceAxisymmetric_Flow(void) { }

void CSourceAxisymmetric_Flow::ComputeResidual(double *val_residual, double **Jacobian_i, CConfig *config) {
  
	double yinv, Pressure_i, Enthalpy_i, Velocity_i, sq_vel;
	unsigned short iDim;
  
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool incompressible  = config->GetIncompressible();
  
	if (Coord_i[1] > 0.0) yinv = 1.0/Coord_i[1];
	else yinv = 0.0;
  
	if (incompressible) {
		val_residual[0] = yinv*Volume*U_i[2]*BetaInc2_i;
		val_residual[1] = yinv*Volume*U_i[1]*U_i[2]/DensityInc_i;
		val_residual[2] = yinv*Volume*U_i[2]*U_i[2]/DensityInc_i;
	}
	else {
		sq_vel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i = U_i[iDim+1] / U_i[0];
			sq_vel += Velocity_i *Velocity_i;
		}
    
		Pressure_i = (Gamma-1.0)*U_i[0]*(U_i[nDim+1]/U_i[0]-0.5*sq_vel);
		Enthalpy_i = (U_i[nDim+1] + Pressure_i) / U_i[0];
    
		val_residual[0] = yinv*Volume*U_i[2];
		val_residual[1] = yinv*Volume*U_i[1]*U_i[2]/U_i[0];
		val_residual[2] = yinv*Volume*(U_i[2]*U_i[2]/U_i[0]);
		val_residual[3] = yinv*Volume*Enthalpy_i*U_i[2];
	}
  
  if (implicit) {
    Jacobian_i[0][0] = 0;
    Jacobian_i[0][1] = 0;
    Jacobian_i[0][2] = 1.;
    Jacobian_i[0][3] = 0;
    
    Jacobian_i[1][0] = -U_i[1]*U_i[2]/(U_i[0]*U_i[0]);
    Jacobian_i[1][1] = U_i[2]/U_i[0];
    Jacobian_i[1][2] = U_i[1]/U_i[0];
    Jacobian_i[1][3] = 0;
    
    Jacobian_i[2][0] = -U_i[2]*U_i[2]/(U_i[0]*U_i[0]);
    Jacobian_i[2][1] = 0;
    Jacobian_i[2][2] = 2*U_i[2]/U_i[0];
    Jacobian_i[2][3] = 0;
    
    Jacobian_i[3][0] = -Gamma*U_i[2]*U_i[3]/(U_i[0]*U_i[0]) + (Gamma-1)*U_i[2]*(U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]*U_i[0]);
    Jacobian_i[3][1] = -(Gamma-1)*U_i[2]*U_i[1]/(U_i[0]*U_i[0]);
    Jacobian_i[3][2] = Gamma*U_i[3]/U_i[0] - 1/2*(Gamma-1)*( (U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]) + 2*U_i[2]*U_i[2]/(U_i[0]*U_i[0]) );
    Jacobian_i[3][3] = Gamma*U_i[2]/U_i[0];
    
    for (int iVar=0; iVar<4; iVar++)
      for (int jVar=0; jVar<4; jVar++)
        Jacobian_i[iVar][jVar] *= yinv*Volume;
  }
}

CSource_JouleHeating::CSource_JouleHeating(unsigned short val_nDim, unsigned short val_nVar,
                                           CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	Velocity = new double [nDim];
	Gamma = config->GetGamma();
	Gas_Constant = config->GetGas_ConstantND();
  
  
}

CSource_JouleHeating::~CSource_JouleHeating(void) {
  
}

void CSource_JouleHeating::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
  
	double Current = 100.0;//config->GetCurrent();
	for (unsigned short iVar = 0; iVar < nVar; iVar ++) {
		val_residual[iVar] = 0.0;
		for (unsigned short jVar = 0; jVar < nVar; jVar ++)
			val_Jacobian_i[iVar][jVar] = 0.0;
    
	}
	//		if (fabs(Integralsqr < 1E-16) || (Integralsqr != Integralsqr)) {
	//			cout << " Integral = "<< Integralsqr << endl;
	//			cin.get();
	//		}
	val_residual[nVar-1] = Current*Current*Elec_Conduct*Volume/(4.0*PI_NUMBER*PI_NUMBER*Integralsqr);
}


void CSource_JouleHeating::SetElec_Cond() {
  
	Density = U_i[0];
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity[iDim] = U_i[iDim+1]/Density;
		sq_vel += Velocity[iDim]*Velocity[iDim];
	}
	Energy = U_i[nDim+1]/Density;
	SoundSpeed = sqrt(Gamma*Gamma_Minus_One*(Energy-0.5*sq_vel));
	Pressure = (Gamma-1.0) * Density * (Energy - 0.5*sq_vel);
	Temperature = Pressure/(Gas_Constant*Density);
	double Patm = Pressure / 101325.0;
  
	double *coeff_a, *coeff_c, *coeff_d;
	coeff_a = new double[8];
	coeff_c = new double[8];
	coeff_d = new double[8];
  
	double w, sigma, x0,x1,x2, x3,x4;
  
	x0 = 1.0; x1 = Patm; x2 = Patm*Patm ; x3 = pow(Patm,3); x4 = pow(Patm,4);
	x1 = log(x1); x2 = log(x2); x3 = log(x3); x4 = log(x4);
  
	coeff_a[0] = exp(1.635045e0*x0 + 4.450390e-2*x1 - 5.928863e-4*x2  + 0.0*x3 + 0.0*x4);
	coeff_c[0] = exp(5.748398e0*x0 + 6.411299e-2*x1 + 0.0*x2    	  + 0.0*x3 + 0.0*x4);
	coeff_d[0] = exp(1.786355e0*x0 - 1.212690e-2*x1 - 2.375673e-4*x2  + 0.0*x3 + 0.0*x4);
	w		   = exp(1.419925e0*x0 - 3.875497e-2*x1 + 0.0*x2	  	  + 0.0*x3 + 0.0*x4);
  
	sigma = coeff_a[0] - coeff_c[0]*exp(-pow(Temperature/coeff_d[0],w));
  
	coeff_c[1] = exp(8.930803e0*x0 + 5.718843e-2*x1 + 1.093759e-3*x2  + 0.0*x3 			+ 0.0*x4);
	coeff_c[2] = exp(8.576847e0*x0 + 1.004174e-1*x1 + 7.406762e-3*x2  - 1.095186e-3*x3  + 0.0*x4);
	coeff_c[3] = exp(1.023493e1*x0 + 6.651575e-2*x1 + 1.090308e-3*x2  - 6.576415e-5*x3 	+ 4.715318e-7*x4);
	coeff_c[4] = exp(1.072380e1*x0 - 5.671452e-2*x1 + 1.468210e-4*x2  + 2.608196e-5*x3 	+ 6.511719e-6*x4);
	coeff_c[5] = exp(1.106431e1*x0 + 5.578774e-2*x1 + 6.639655e-4*x2  - 0.0*x3 			+ 0.0*x4);
	coeff_c[6] = exp(1.023203e1*x0 + 8.703300e-2*x1 + 5.007602e-3*x2  + 0.0*x3 			+ 0.0*x4);
	coeff_c[7] = exp(2.946755e1*x0 - 4.289010e0*x1  - 3.224136e-1*x2  + 9.371814e-2*x3 	+ 0.0*x4);
  
  
	coeff_d[1] = exp(7.014976e0*x0 +  7.625175e-2*x1 + 3.011941e-4*x2  + 0.0*x3 		 + 0.0*x4);
	coeff_d[2] = exp(9.113182e0*x0 -  8.202725e-2*x1 + 6.299430e-3*x2  + 9.099828e-4*x3  + 0.0*x4);
	coeff_d[3] = exp(8.039563e0*x0 +  1.435966e-1*x1 + 8.862611e-3*x2  - 3.478227e-4*x3  - 3.614711e-5*x4);
	coeff_d[4] = exp(8.556977e0*x0 +  2.227207e-1*x1 - 2.773160e-3*x2  - 1.684219e-3*x3  + 1.878188e-4*x4);
	coeff_d[5] = exp(9.309043e0*x0 +  1.366510e-1*x1 - 2.599317e-3*x2  + 0.0*x3 		 + 0.0*x4);
	coeff_d[6] = exp(1.130562e1*x0 -  2.184155e-2*x1 - 1.865543e-4*x2  + 0.0*x3 		 + 0.0*x4);
	coeff_d[7] = exp(2.430324e1*x0 -  2.653523e0*x1  - 3.309222e-1*x2  + 4.769061e-2*x3  + 0.0*x4);
  
  
	coeff_a[1] =  exp(4.493934e-2*x0 -  9.063708e-3*x1 - 2.489500e-3*x2  + 0.0*x3 		   + 0.0*x4);
	x0 = 1.0; 	  x1 = log(Patm); x2 = x1*x1 ; x3 = pow(x1,3); x4 = pow(x1,4);
  
	coeff_a[2] =  (1.593153e0*x0  +  4.137850e-2*x1 + 1.430491e-2*x2  - 4.403957e-7*x3  + 0.0*x4);
	coeff_a[3] = -(2.627897e-1*x0 +  2.917558e-3*x1 + 3.036205e-3*x2  - 1.926651e-4*x3  - 2.917018e-5*x4);
	coeff_a[4] = -(1.707216e-1*x0 +  2.035164e-2*x1 + 1.809127e-3*x2  - 9.630175e-5*x3  + 1.781249e-5*x4);
	coeff_a[5] = -(2.480007e-1*x0 +  2.217818e-2*x1 + 9.094614e-4*x2  + 0.0*x3 		   + 0.0*x4);
	coeff_a[6] =  (3.636707e0*x0  -  1.436268e-1*x1 - 2.934926e-3*x2  + 0.0*x3 		   + 0.0*x4);
	coeff_a[7] =  coeff_a[3] + coeff_a[4] + coeff_a[5] - coeff_a[1] - coeff_a[2] - coeff_a[6];
  
  
	double q = 0;
	for (int i = 1; i < 8; ++i) {
		q = (Temperature - coeff_c[i]) / coeff_d[i];
		sigma = sigma + coeff_a[i] * q/(exp(q) + exp(-q));
	}
  
	Elec_Conduct = exp(sigma);
  
	//	if ( Elec_Conduct != Elec_Conduct) {
	//	cout << " Elec_Cond in calculation = " << Elec_Conduct << endl;
	//	cout << "Density = " << Density << endl;
	//	cout << "Energy = " << Energy << endl;
	//	cout << "Patm = " << Patm << endl;
	//	cout << "Temperature = " << Temperature << endl;
	//	cout << "SoundSpeed = " << SoundSpeed << endl;
	//	cout << "sigma = " << sigma << endl;
	//	for (int i = 0; i < 8; ++i) {
	//		q = (Temperature - coeff_c[i]) / coeff_d[i];
	//		cout << " i = " << i << ",  q = " << q  << " coeff_a[i] = " << coeff_a[i];
	//		cout << " coeff_c[i] = " << coeff_c[i] << " coeff_d[i] = " << coeff_d[i] << endl;
	//	}
	//	cin.get();
	//	}
}



