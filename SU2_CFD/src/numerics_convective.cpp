/*!
 * \file numerics_convective.cpp
 * \brief This file contains all the convective term discretization.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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

void CUpwRoe_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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

	/*--- Adjustment for a rotating frame: the velocity projection
        must include a contribution from rotation when computing
        the eigenvalues ---*/
	if (rotating_frame) {
		//		double ProjRotVel = 0.0; double ProjRotVel_i = 0.0; double ProjRotVel_j = 0.0;
		//		for (iDim = 0; iDim < nDim; iDim++) {
		//			ProjRotVel   += 0.5*(RotVel_i[iDim]+RotVel_j[iDim])*UnitaryNormal[iDim];
		//			ProjRotVel_i += RotVel_i[iDim]*UnitaryNormal[iDim];
		//			ProjRotVel_j += RotVel_j[iDim]*UnitaryNormal[iDim];
		//		}
		//		ProjVelocity   -= ProjRotVel;
		//    ProjVelocity_i -= ProjRotVel_i;
		//		ProjVelocity_j -= ProjRotVel_j;
		ProjVelocity   -= Rot_Flux/Area;
		ProjVelocity_i -= Rot_Flux/Area;
		ProjVelocity_j -= Rot_Flux/Area;
	}

	/*--- Flow eigenvalues and Entropy correctors ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Lambda[iDim] = ProjVelocity;
		Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,ProjVelocity_j-Lambda[iDim]));
	}
	Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
	Epsilon[nVar-2] = 4.0*max(0.0, max(Lambda[nVar-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar-2]));
	Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
	Epsilon[nVar-1] = 4.0*max(0.0, max(Lambda[nVar-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar-1]));

	/*--- Entropy correction ---*/
	/*	for (iVar = 0; iVar < nVar; iVar++)
		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
		else
			Lambda[iVar] = fabs(Lambda[iVar]);*/

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

		if (nDim == 3) {
			delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
			delta_wave[1] = UnitaryNormal[0]*delta_vel[2]-UnitaryNormal[2]*delta_vel[0];
			delta_wave[2] = UnitaryNormal[1]*delta_vel[0]-UnitaryNormal[0]*delta_vel[1];
			delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
			delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
		}
		else {
			delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
			delta_wave[1] = UnitaryNormal[1]*delta_vel[0]-UnitaryNormal[0]*delta_vel[1];
			delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
			delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
		}

		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
			for (jVar = 0; jVar < nVar; jVar++)
				val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
		}

		/*--- Flux contribution for a rotating frame ---*/
		if (rotating_frame) {
			//			ProjVelocity = 0.0;
			//			for (iDim = 0; iDim < nDim; iDim++)
			//				ProjVelocity += 0.5*(RotVel_i[iDim]+RotVel_j[iDim])*Normal[iDim];
			ProjVelocity = Rot_Flux;
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
		GetInviscidProjJac(Velocity_i, Energy_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidProjJac(Velocity_j, Energy_j, Normal, 0.5, val_Jacobian_j); 

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

		/*--- Contributions for a rotating frame ---*/
		if (rotating_frame) {
			//			ProjVelocity = 0.0;
			//			for (iDim = 0; iDim < nDim; iDim++)
			//				ProjVelocity += 0.5*(RotVel_i[iDim]+RotVel_j[iDim])*Normal[iDim];
			ProjVelocity = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
				/*--- Implicit terms ---*/
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

void CUpwRoeArtComp_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	/*--- Compute area nad unitary normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);

	double eps = numeric_limits<double>::epsilon();

	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = eps + Normal[iDim]/Area;

	/*--- Set the variables at point i ---*/
	Pressure_i = U_i[0];
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;

	/*--- Variables at point j ---*/
	Pressure_j = U_j[0];
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;

	/*--- Mean variables at points iPoint and jPoint ---*/
	MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
	MeanPressure = 0.5*(Pressure_i + Pressure_j);
	MeanBetaInc2 = 0.5*(BetaInc2_i + BetaInc2_j);

	ProjVelocity = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);	
		ProjVelocity += MeanVelocity[iDim]*Normal[iDim];
	}

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

	for (iVar = 0; iVar < nVar; iVar++)
		Lambda[iVar] = fabs(Lambda[iVar]);

	/*--- Compute inverse P ---*/
	GetPArtCompMatrix_inv(&MeanDensity, MeanVelocity, &MeanBetaInc2, UnitaryNormal, invP_Tensor);

	if (implicit) {
		/*--- Jacobian of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
		GetInviscidArtCompProjJac(DensityInc_i, Velocity_i, BetaInc2_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidArtCompProjJac(DensityInc_j, Velocity_j, BetaInc2_j, Normal, 0.5, val_Jacobian_j); 
	}

	/*--- TEST JACOBIAN DECOMPOSITION ----*/
	/*	ProjVelocity_i = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
	SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i));
	GetInviscidArtCompProjJac(DensityInc_i, Velocity_i, BetaInc2_i, UnitaryNormal, 1.0, val_Jacobian_i);
	cout <<" Jacobian "<< val_Jacobian_i[0][0]<<" "<< val_Jacobian_i[0][1]<<" "<< val_Jacobian_i[0][2]<<" "<< val_Jacobian_i[0][3]<<endl;
	cout <<" Jacobian "<< val_Jacobian_i[1][0]<<" "<< val_Jacobian_i[1][1]<<" "<< val_Jacobian_i[1][2]<<" "<< val_Jacobian_i[1][3]<<endl;
	cout <<" Jacobian "<< val_Jacobian_i[2][0]<<" "<< val_Jacobian_i[2][1]<<" "<< val_Jacobian_i[2][2]<<" "<< val_Jacobian_i[2][3]<<endl;
	cout <<" Jacobian "<< val_Jacobian_i[3][0]<<" "<< val_Jacobian_i[3][1]<<" "<< val_Jacobian_i[3][2]<<" "<< val_Jacobian_i[3][3]<<endl;
	GetPArtCompMatrix(&DensityInc_i, Velocity_i, &BetaInc2_i, UnitaryNormal, P_Tensor);
	GetPArtCompMatrix_inv(&DensityInc_i, Velocity_i, &BetaInc2_i, UnitaryNormal, invP_Tensor);
	Lambda[0] = ProjVelocity_i;
	Lambda[1] = ProjVelocity_i;
	Lambda[2] = ProjVelocity_i + SoundSpeed_i;
	Lambda[3] = ProjVelocity_i - SoundSpeed_i;
	double JacobianNew[4][4];
	for (iVar = 0; iVar < nVar; iVar++)
		for (jVar = 0; jVar < nVar; jVar++) 
			JacobianNew[iVar][jVar] = 0.0;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) { 
			Proj_ModJac_Tensor_ij = 0.0;
			for (kVar = 0; kVar < nVar; kVar++)
				JacobianNew[iVar][jVar] += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
		}
	}
	cout <<" JacobianNew "<< JacobianNew[0][0]<<" "<< JacobianNew[0][1]<<" "<< JacobianNew[0][2] <<" "<< JacobianNew[0][3]<<endl;
	cout <<" JacobianNew "<< JacobianNew[1][0]<<" "<< JacobianNew[1][1]<<" "<< JacobianNew[1][2] <<" "<< JacobianNew[1][3]<<endl;
	cout <<" JacobianNew "<< JacobianNew[2][0]<<" "<< JacobianNew[2][1]<<" "<< JacobianNew[2][2] <<" "<< JacobianNew[2][3]<<endl;
	cout <<" JacobianNew "<< JacobianNew[3][0]<<" "<< JacobianNew[3][1]<<" "<< JacobianNew[3][2] <<" "<< JacobianNew[3][3]<<endl;
	cout <<" Normal Vector "<< UnitaryNormal[0] <<" "<< UnitaryNormal[1] <<" "<<UnitaryNormal[2] << endl;

	cout << endl;

	cin.get();*/

	/*--- Diference variables iPoint and jPoint ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Diff_U[iVar] = U_j[iVar] - U_i[iVar];

	/*--- Flux approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
		for (jVar = 0; jVar < nVar; jVar++) { 
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
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

CUpwRoe_AdjFlow::CUpwRoe_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Residual_Roe = new double [nVar];
	RoeVelocity = new double [nDim];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	Lambda = new double [nVar];
	P_Tensor = new double* [nVar];
	invP_Tensor = new double* [nVar];
	Proj_flux_tensor_i = new double*[nVar];
	Proj_flux_tensor_j = new double*[nVar];
	Proj_ModJac_Tensor = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar] = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
		Proj_flux_tensor_i[iVar] = new double[nVar];
		Proj_flux_tensor_j[iVar] = new double[nVar];
		Proj_ModJac_Tensor[iVar] = new double[nVar];
	}

}

CUpwRoe_AdjFlow::~CUpwRoe_AdjFlow(void) {

	delete [] Residual_Roe;
	delete [] RoeVelocity;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Lambda;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
		delete [] Proj_flux_tensor_i[iVar];
		delete [] Proj_flux_tensor_j[iVar];
		delete [] Proj_ModJac_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
	delete [] Proj_ModJac_Tensor;

}

void CUpwRoe_AdjFlow::SetResidual (double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, 
		double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,CConfig *config) {

	/*--- Compute the area ---*/
	area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		area += Normal[iDim]*Normal[iDim];
	area = sqrt(area);
	rarea = 1.0 / area;

	/*--- Components of the normal & unit normal vector of the current face ---*/
	Sx = Normal[0]; 
	Sy = Normal[1];
	Sz = 0.0; if (nDim == 3) Sz = Normal[2];
	nx = Sx * rarea; 
	ny = Sy * rarea;
	nz = Sz * rarea;

	/*--- Flow variable states at point i (left, _l) and j (right, _r)---*/
	rho_l  = U_i[0]; rho_r  = U_j[0];
	u_l = U_i[1]/U_i[0]; v_l = U_i[2]/U_i[0]; w_l = 0.0;
	u_r = U_j[1]/U_j[0]; v_r = U_j[2]/U_j[0]; w_r = 0.0;
	if (nDim == 3) w_l = U_i[3]/U_i[0];
	if (nDim == 3) w_r = U_j[3]/U_j[0];
	h_l = Enthalpy_i; h_r = Enthalpy_j;

	/*--- One-half speed squared ---*/
	q_l = ONE2 * ((u_l*u_l) + (v_l*v_l) + (w_l*w_l));
	q_r = ONE2 * ((u_r*u_r) + (v_r*v_r) + (w_r*w_r));

	/*--- Projected velocity ---*/
	Q_l = (u_l * Sx) + (v_l * Sy) + (w_l * Sz);
	Q_r = (u_r * Sx) + (v_r * Sy) + (w_r * Sz);

	/*--- Mean adjoint variables ---*/
	psi1 = ONE2 * (Psi_i[0] + Psi_j[0]);
	psi2 = ONE2 * (Psi_i[1] + Psi_j[1]);
	psi3 = ONE2 * (Psi_i[2] + Psi_j[2]);
	psi4 = 0.0; if (nDim == 3) psi4 = ONE2 * (Psi_i[3] + Psi_j[3]);
	psi5 = ONE2 * (Psi_i[nVar-1] + Psi_j[nVar-1]);

	/*--- Left state ---*/
	l1psi = (Sx * psi2) + (Sy * psi3) + (Sz * psi4) + (Q_l * psi5);
	l2psi = psi1 + (u_l * psi2) + (v_l * psi3) + (w_l * psi4) + (h_l * psi5);

	val_residual_i[0] = Q_l * psi1 - l2psi * Q_l + l1psi * Gamma_Minus_One * q_l;
	val_residual_i[1] = Q_l * psi2 + l2psi * Sx  - l1psi * Gamma_Minus_One * u_l;
	val_residual_i[2] = Q_l * psi3 + l2psi * Sy  - l1psi * Gamma_Minus_One * v_l;
	if (nDim == 3) val_residual_i[3] = Q_l * psi4 + l2psi * Sz  - l1psi * Gamma_Minus_One * w_l;
	val_residual_i[nVar-1] = Q_l * psi5 + l1psi * Gamma_Minus_One;

	/*--- Right state ---*/
	l1psi = (Sx * psi2) + (Sy * psi3) + (Sz * psi4) + (Q_r * psi5);
	l2psi = psi1 + (u_r * psi2) + (v_r * psi3) + (w_r * psi4) + (h_r * psi5);

	val_residual_j[0] = -(Q_r * psi1 - l2psi * Q_r + l1psi * Gamma_Minus_One * q_r);
	val_residual_j[1] = -(Q_r * psi2 + l2psi * Sx  - l1psi * Gamma_Minus_One * u_r);
	val_residual_j[2] = -(Q_r * psi3 + l2psi * Sy  - l1psi * Gamma_Minus_One * v_r);
	if (nDim == 3) val_residual_j[3] = -(Q_r * psi4 + l2psi * Sz  - l1psi * Gamma_Minus_One * w_r);
	val_residual_j[nVar-1] = -(Q_r * psi5 + l1psi * Gamma_Minus_One);


	/*--- f_{roe} = P^{-T} |lambda| P^T \delta \psi ---*/ 
	psi1_l = Psi_i[0];
	psi2_l = Psi_i[1];
	psi3_l = Psi_i[2];
	psi4_l = 0.0; if (nDim == 3) psi4_l = Psi_i[3];
	psi5_l = Psi_i[nVar-1];

	psi1_r = Psi_j[0];
	psi2_r = Psi_j[1];
	psi3_r = Psi_j[2];
	psi4_r = 0.0; if (nDim == 3) psi4_r = Psi_j[3];
	psi5_r = Psi_j[nVar-1];

	/*--- Roe averaging ---*/
	rrho_l   = 1.0 / rho_l;
	weight   = sqrt(rho_r * rrho_l);
	rweight1 = 1.0 / (1.0 + weight);
	weight  *= rweight1;

	h = h_l * rweight1 + weight * h_r;
	u = u_l * rweight1 + weight * u_r;
	v = v_l * rweight1 + weight * v_r;
	w = w_l * rweight1 + weight * w_r;

	psi1 = ONE2 * (psi1_r - psi1_l);
	psi2 = ONE2 * (psi2_r - psi2_l);
	psi3 = ONE2 * (psi3_r - psi3_l);
	psi4 = ONE2 * (psi4_r - psi4_l);
	psi5 = ONE2 * (psi5_r - psi5_l);

	q2 = (u*u) + (v*v) + (w*w);
	Q  = (u * Sx) + (v * Sy) + (w * Sz);
	vn = nx * u   + ny * v   + nz * w;
	cc = Gamma_Minus_One * h - 0.5 * Gamma_Minus_One * q2;
	c  = sqrt(cc);

	/*--- Adjustment for a rotating frame: the velocity projection
   must include a contribution from rotation when computing
   the eigenvalues ---*/
	if (rotating_frame) {
		double ProjRotVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjRotVel += 0.5*(RotVel_i[iDim]+RotVel_j[iDim])*Normal[iDim];
		ProjRotVel = Rot_Flux;
		Q -= ProjRotVel;
	}

	/*--- Eigenvalues from the primal solution ---*/
	absQ  = fabs(Q);
	absQp = fabs(Q + c * area);
	absQm = fabs(Q - c * area);

	alpha  = ONE2 * Gamma_Minus_One * q2 / cc;
	beta_u = psi2 + u * psi5;
	beta_v = psi3 + v * psi5;
	beta_w = psi4 + w * psi5;
	eta    = Gamma_Minus_One / cc;
	l1psi  = (nx * psi2) + (ny * psi3) + (nz * psi4) + (vn * psi5);
	l2psi  = psi1 + (u * psi2) + (v * psi3) + (w * psi4) + (h * psi5);
	l1l2p  = (l2psi + c * l1psi) * absQp;
	l1l2m  = (l2psi - c * l1psi) * absQm;

	/*--- adjoint flux computation in the x,y and z coordinate system ---*/
	Residual_Roe[0] = ((1.0-alpha)*l2psi - (1.0-alpha)*cc/Gamma_Minus_One*psi5
			- u*beta_u*(1.0-(nx*nx)) - v*beta_v*(1.0-(ny*ny))
			- w*beta_w*(1.0-(nz*nz)) + ny*nz*(w*beta_v + v*beta_w)
			+ nx*nz*(w*beta_u + u*beta_w) + ny*nx*(v*beta_u + u*beta_v) ) * absQ
			- ONE2 / c * vn * (l1l2p - l1l2m) + ONE2 * alpha *  (l1l2p + l1l2m);

	Residual_Roe[1] = (l2psi*u*eta - u*psi5 + beta_u*(1.0-(nx*nx))
			- nx*(beta_v*ny + beta_w*nz) ) * absQ + ONE2*nx/c  * (l1l2p - l1l2m )
			- ONE2*eta*u * (l1l2p + l1l2m );

	Residual_Roe[2] = (l2psi*v*eta - v*psi5 + beta_v*(1.0-(ny*ny))
			- ny*(beta_w*nz + beta_u*nx) ) * absQ + ONE2*ny/c  * (l1l2p - l1l2m )
			- ONE2*eta*v * (l1l2p + l1l2m );

	if (nDim == 3) Residual_Roe[3] = (l2psi*w*eta - w*psi5 + beta_w*(1.0-(nz*nz)) - nz*(beta_u*nx + beta_v*ny) ) * absQ 
			+ ONE2*nz/c  * (l1l2p - l1l2m ) - ONE2*eta*w * (l1l2p + l1l2m );

	Residual_Roe[nVar-1] = (psi5 - l2psi*eta) * absQ + ONE2*eta*(l1l2p + l1l2m);

	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual_i[iVar]   += Residual_Roe[iVar];
		val_residual_j[iVar]   -= Residual_Roe[iVar];
	}

	/*--- Flux contribution for a rotating frame ---*/
	if (rotating_frame) {
		double ProjVelocity = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjVelocity += 0.5*(RotVel_i[iDim]+RotVel_j[iDim])*Normal[iDim];
		ProjVelocity = Rot_Flux;
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual_i[iVar] -= ProjVelocity * 0.5*(Psi_i[iVar]+Psi_j[iVar]);
			val_residual_j[iVar] += ProjVelocity * 0.5*(Psi_i[iVar]+Psi_j[iVar]);
		}
	}


	/*--- Implicit Contributions ---*/
	if (implicit) {

		/*--- Prepare variables for use in matrix routines ---*/
		RoeDensity = U_i[0]*sqrt(U_j[0]/U_i[0]);
		RoeSoundSpeed = c;
		UnitaryNormal[0] = nx;  UnitaryNormal[1] = ny;  if (nDim == 3 ) UnitaryNormal[2] = nz;
		RoeVelocity[0]   = u;   RoeVelocity[1]   = v;   if (nDim == 3 ) RoeVelocity[2]   = w;
		Velocity_i[0]    = u_l; Velocity_i[1]    = v_l; if (nDim == 3 ) Velocity_i[2]    = w_l;
		Velocity_j[0]    = u_r; Velocity_j[1]    = v_r; if (nDim == 3 ) Velocity_j[2]    = w_r;
		Energy_i = U_i[nDim+1] / U_i[0]; Energy_j = U_j[nDim+1] / U_j[0];

		/*--- Jacobians of the inviscid flux, scaled by
		 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
		GetInviscidProjJac(Velocity_i, Energy_i, Normal, 0.5, Proj_flux_tensor_i);
		GetInviscidProjJac(Velocity_j, Energy_j, Normal, 0.5, Proj_flux_tensor_j);

		/*--- Compute P, inverse P, and store eigenvalues ---*/
		GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, invP_Tensor);
		GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, P_Tensor);

		/*--- Flow eigenvalues ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			Lambda[iDim] = absQ;
		Lambda[nVar-2] = absQp;
		Lambda[nVar-1] = absQm;

		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) { 
				Proj_ModJac_Tensor_ij = 0.0;
				/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
				for (kVar = 0; kVar < nVar; kVar++) 
					Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
				Proj_ModJac_Tensor[iVar][jVar] = 0.5*Proj_ModJac_Tensor_ij*area;
			}
		}

		/*--- Transpose the matrices and store the Jacobians. Note the negative
		 sign for the ji and jj Jacobians bc the normal direction is flipped. ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) { 
				val_Jacobian_ii[jVar][iVar] = Proj_flux_tensor_i[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ij[jVar][iVar] = Proj_flux_tensor_i[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ji[jVar][iVar] = -(Proj_flux_tensor_j[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar]);
				val_Jacobian_jj[jVar][iVar] = -(Proj_flux_tensor_j[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar]);
			}
		}

		/*--- Contributions for a rotating frame ---*/
		if (rotating_frame) {
			double ProjVelocity = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjVelocity += 0.5*(RotVel_i[iDim]+RotVel_j[iDim])*Normal[iDim];
			ProjVelocity = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				/*--- Adjust Jacobian main diagonal ---*/
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjVelocity;
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjVelocity;
			}
		}

	}
}


CUpwRoeArtComp_AdjFlow::CUpwRoeArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

	MeanVelocity = new double [nDim];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	Lambda = new double [nVar];
	P_Tensor = new double* [nVar];
	invP_Tensor = new double* [nVar];
	Proj_Jac_Tensor_i = new double*[nVar];
	Proj_Jac_Tensor_j = new double*[nVar];
	Proj_ModJac_Tensor = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar] = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
		Proj_Jac_Tensor_i[iVar] = new double[nVar];
		Proj_Jac_Tensor_j[iVar] = new double[nVar];
		Proj_ModJac_Tensor[iVar] = new double[nVar];
	}

}

CUpwRoeArtComp_AdjFlow::~CUpwRoeArtComp_AdjFlow(void) {

	delete [] MeanVelocity;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Lambda;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
		delete [] Proj_Jac_Tensor_i[iVar];
		delete [] Proj_Jac_Tensor_j[iVar];
		delete [] Proj_ModJac_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
	delete [] Proj_Jac_Tensor_i;
	delete [] Proj_Jac_Tensor_j;
	delete [] Proj_ModJac_Tensor;

}

void CUpwRoeArtComp_AdjFlow::SetResidual (double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, 
		double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,CConfig *config) {

	/*--- Compute area nad unitary normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);

	/*--- Components of the normal & unit normal vector of the current face ---*/	
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = numeric_limits<double>::epsilon() + Normal[iDim]/Area;

	/*--- Set the variables at point i, and j ---*/
	Pressure_i = U_i[0]; Pressure_j = U_j[0];
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
	}

	/*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidArtCompProjJac(DensityInc_i, Velocity_i, BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
	GetInviscidArtCompProjJac(DensityInc_j, Velocity_j, BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);

	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual_i[iVar] = 0.0; val_residual_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_residual_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
			val_residual_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
		}
	}

	/*--- Mean variables at points iPoint and jPoint ---*/
	MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
	MeanPressure = 0.5*(Pressure_i + Pressure_j);
	MeanBetaInc2 = 0.5*(BetaInc2_i + BetaInc2_j);

	ProjVelocity = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);	
		ProjVelocity += MeanVelocity[iDim]*Normal[iDim];
	}

	MeanSoundSpeed = sqrt(ProjVelocity*ProjVelocity + (MeanBetaInc2/MeanDensity) * Area * Area);

	/*--- Compute P, inverse P, and store eigenvalues ---*/
	GetPArtCompMatrix_inv(&MeanDensity, MeanVelocity, &MeanBetaInc2, UnitaryNormal, invP_Tensor);
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

	for (iVar = 0; iVar < nVar; iVar++)
		Lambda[iVar] = fabs(Lambda[iVar]);


	/*--- Flux approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) { 
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++) 
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_ModJac_Tensor[iVar][jVar] = 0.5*Proj_ModJac_Tensor_ij;
		}
	}

	for (iVar = 0; iVar < nVar; iVar++)
		for (jVar = 0; jVar < nVar; jVar++) {
			val_residual_i[iVar] -= Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar]);
			val_residual_j[iVar] += Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar]);
		}

	/*--- Implicit contributions, Transpose the matrices and store the Jacobians. Note the negative
	 sign for the ji and jj Jacobians bc the normal direction is flipped. ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) { 
				val_Jacobian_ii[jVar][iVar] = Proj_Jac_Tensor_i[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ij[jVar][iVar] = Proj_Jac_Tensor_i[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ji[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar]);
				val_Jacobian_jj[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar]);
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

void CUpwAUSM_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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
		GetInviscidProjJac(Velocity_i, Energy_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidProjJac(Velocity_j, Energy_j, Normal, 0.5, val_Jacobian_j); 

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

void CUpwHLLC_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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
		GetInviscidProjJac(Velocity_i, Energy_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidProjJac(Velocity_j, Energy_j, Normal, 0.5, val_Jacobian_j); 

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

CUpwLin_LevelSet::CUpwLin_LevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];

}

CUpwLin_LevelSet::~CUpwLin_LevelSet(void) {

	delete [] Velocity_i;
	delete [] Velocity_j;

}

void CUpwLin_LevelSet::SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	unsigned short iDim;
	double a0, a1, q_ij;

	q_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i; 
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
		q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
	}

	a0 = 0.5*(q_ij+fabs(q_ij)); a1 = 0.5*(q_ij-fabs(q_ij));

	val_residual[0] = a0*LevelSetVar_i[0]+a1*LevelSetVar_j[0];

	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_j[0][0] = a1;
	}

}

CUpwLin_AdjLevelSet::CUpwLin_AdjLevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];

}

CUpwLin_AdjLevelSet::~CUpwLin_AdjLevelSet(void) {

	delete [] Velocity_i;
	delete [] Velocity_j;

}

void CUpwLin_AdjLevelSet::SetResidual (double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, 
		double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config)  {

	unsigned short iDim;
	double a0, a1, q_ij;

	/*--- i point ---*/
	q_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i; Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
		q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
	}

	a0 = 0.5*(q_ij+fabs(q_ij)); a1 = 0.5*(q_ij-fabs(q_ij));

	val_residual_i[0] = a0*LevelSetVar_i[0]+a1*LevelSetVar_j[0];

	if (implicit) {
		val_Jacobian_ii[0][0] = a0;
		val_Jacobian_ij[0][0] = a1;
	}

	/*--- j point ---*/
	q_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i; Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
		q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
	}

	a0 = 0.5*(q_ij+fabs(q_ij)); a1 = 0.5*(q_ij-fabs(q_ij));

	val_residual_j[0] = -(a0*LevelSetVar_i[0]+a1*LevelSetVar_j[0]);

	if (implicit) {
		val_Jacobian_ji[0][0] = -a0;
		val_Jacobian_jj[0][0] = -a1;
	}

}

CUpwLin_TurbSA::CUpwLin_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];

}

CUpwLin_TurbSA::~CUpwLin_TurbSA(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwLin_TurbSA::SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	unsigned short iDim;

	switch (config->GetKind_Turb_Model()) {
	case SA : case SA_COMP :
		double Density_i = U_i[0];
		double Density_j = U_j[0];

		double q_ij = 0;
		if (!rotating_frame) {
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity_i[iDim] = U_i[iDim+1]/Density_i;
				Velocity_j[iDim] = U_j[iDim+1]/Density_j;
				q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
			}
		}
		else {
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity_i[iDim] = U_i[iDim+1]/Density_i - RotVel_i[iDim];
				Velocity_j[iDim] = U_j[iDim+1]/Density_j - RotVel_j[iDim];
				q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
			}
		}

		double a0 = 0.5*(q_ij+fabs(q_ij));
		double a1 = 0.5*(q_ij-fabs(q_ij));
		val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];

		if (implicit)
			val_Jacobian_i[0][0] = a0;
		break;
	}
}

CUpwLin_TurbSST::CUpwLin_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
}

CUpwLin_TurbSST::~CUpwLin_TurbSST(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwLin_TurbSST::SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);

	double Density_i = U_i[0];
	double Density_j = U_j[0];

	double q_ij = 0;
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/Density_i;
		Velocity_j[iDim] = U_j[iDim+1]/Density_j;
		q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
	}

	double a0 = 0.5*(q_ij+fabs(q_ij));
	double a1 = 0.5*(q_ij-fabs(q_ij));
	val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];

	if (implicit)
		val_Jacobian_i[0][0] = a0;

}


CUpwLin_AdjTurb::CUpwLin_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
}

CUpwLin_AdjTurb::~CUpwLin_AdjTurb(void) {
	delete [] Velocity_i;
}

void CUpwLin_AdjTurb::SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);

	/*--- Non-conservative term  -->  -\nabla \psi_\mu  B^{cv}
	 B^{cv} = -v ---*/

	unsigned short iDim;
	double proj_conv_flux = 0;

	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/U_i[0];
		proj_conv_flux += Velocity_i[iDim]*Normal[iDim]; // projection of convective flux at iPoint
	}
	double psinu0 = TurbPsi_i[0];
	double psinu1;
	if (proj_conv_flux > 0)
		psinu1 = psinu0 + proj_conv_flux;
	else
		psinu1 = psinu0;

	val_residual[0] = 0.5*( proj_conv_flux*(psinu0+psinu1)-fabs(proj_conv_flux)*(psinu1-psinu0));
	if (implicit) {
		val_Jacobian_i[0][0] = 0.5*( proj_conv_flux + fabs(proj_conv_flux));
	}
}

CUpwSca_TurbSA::CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
}

CUpwSca_TurbSA::~CUpwSca_TurbSA(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwSca_TurbSA::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	bool rotating_frame = config->GetRotating_Frame();

	Density_i = U_i[0];
	Density_j = U_j[0];

	q_ij = 0;
	if (!rotating_frame) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i;
			Velocity_j[iDim] = U_j[iDim+1]/Density_j;
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	}
	else {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - RotVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - RotVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	}

	a0 = 0.5*(q_ij+fabs(q_ij));
	a1 = 0.5*(q_ij-fabs(q_ij));
	val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];

	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_j[0][0] = a1;
	}

}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
}

CUpwSca_TurbSST::~CUpwSca_TurbSST(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwSca_TurbSST::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	Density_i = U_i[0];
	Density_j = U_j[0];
	q_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/Density_i;
		Velocity_j[iDim] = U_j[iDim+1]/Density_j;
		q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
	}

	a0 = 0.5*(q_ij+fabs(q_ij));
	a1 = 0.5*(q_ij-fabs(q_ij));

	val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
	val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];

	if (implicit) {
		val_Jacobian_i[0][0] = a0;		val_Jacobian_i[0][1] = 0.0;
		val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = a0;

		val_Jacobian_j[0][0] = a1;		val_Jacobian_j[0][1] = 0.0;
		val_Jacobian_j[1][0] = 0.0;		val_Jacobian_j[1][1] = a1;
	}

}

CUpwSca_AdjTurb::CUpwSca_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
}

CUpwSca_AdjTurb::~CUpwSca_AdjTurb(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwSca_AdjTurb::SetResidual (double *val_residual_i, double *val_residual_j, 
		double **val_Jacobian_ii, double **val_Jacobian_ij,
		double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) {

	bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);

	/*--- Non-conservative term  -->  -\nabla \psi_\mu  B^{cv}
	 B^{cv} = -\nabla \hat{nu}/\sigma + v ---*/

	unsigned short iDim;
	double proj_conv_flux_i = 0, proj_conv_flux_j = 0, proj_conv_flux_ij = 0;
	double sigma = 2./3.;

	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/U_i[0];
		Velocity_j[iDim] = U_j[iDim+1]/U_j[0];
		proj_conv_flux_i += (TurbVar_Grad_i[0][iDim]/sigma - Velocity_i[iDim])*Normal[iDim]; // projection of convective flux at iPoint
		proj_conv_flux_j += (TurbVar_Grad_j[0][iDim]/sigma - Velocity_j[iDim])*Normal[iDim]; // projection of convective flux at jPoint
	}
	proj_conv_flux_ij = 0.5*fabs(proj_conv_flux_i+proj_conv_flux_j); // projection of average convective flux

	val_residual_i[0] = 0.5*( proj_conv_flux_i*(TurbPsi_i[0]+TurbPsi_j[0])-proj_conv_flux_ij*(TurbPsi_j[0]-TurbPsi_i[0]));
	val_residual_j[0] = 0.5*(-proj_conv_flux_j*(TurbPsi_j[0]+TurbPsi_i[0])-proj_conv_flux_ij*(TurbPsi_i[0]-TurbPsi_j[0]));
	if (implicit) {
		val_Jacobian_ii[0][0] = 0.5*( proj_conv_flux_i+proj_conv_flux_ij);
		val_Jacobian_ij[0][0] = 0.5*( proj_conv_flux_i-proj_conv_flux_ij);
		val_Jacobian_ji[0][0] = 0.5*(-proj_conv_flux_j-proj_conv_flux_ij);
		val_Jacobian_jj[0][0] = 0.5*(-proj_conv_flux_j+proj_conv_flux_ij);
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

void CCentJST_Flow::SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
		bool art_diss, CConfig *config) {

	/*--- Conservative variables at point i and 1 ---*/
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
		GetInviscidProjJac(MeanVelocity, MeanEnergy, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}

	/*--- If there is any grid movement ---*/
	if (grid_movement) {
		ProjVelocity = 0;
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

	/*--- If we are in a rotating frame ---*/
	if (rotating_frame) {
		//		ProjVelocity = 0.0;
		//		for (iDim = 0; iDim < nDim; iDim++)
		//			ProjVelocity += 0.5*(RotVel_i[iDim]+RotVel_j[iDim])*Normal[iDim];
		ProjVelocity = Rot_Flux;
		for (iVar = 0; iVar < nVar; iVar++) {
			val_resconv[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			if (implicit) {
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
	}

	if (!art_diss) return;

	/*--- Computes differences btw. Laplacians and conservative variables ---*/
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

	if (rotating_frame) {
		//		double ProjRotVel_i = 0.0; double ProjRotVel_j = 0.0;
		//		for (iDim = 0; iDim < nDim; iDim++) {
		//			ProjRotVel_i += RotVel_i[iDim]*Normal[iDim];
		//			ProjRotVel_j += RotVel_j[iDim]*Normal[iDim];
		//		}
		//		ProjVelocity_i -= ProjRotVel_i;
		//		ProjVelocity_j -= ProjRotVel_j;
		ProjVelocity_i -= Rot_Flux;
		ProjVelocity_j -= Rot_Flux;
	}

	Area = sqrt(Area);
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS), Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS), Param_p);	
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

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

		for (iVar = 0; iVar < (nVar-1); iVar++) {
			val_Jacobian_i[iVar][iVar] += cte_0;
			val_Jacobian_j[iVar][iVar] -= cte_1;
		}

		/*--- Last rows: CAREFUL!! You have differences of \rho_Enthalpy, not differences of \rho_Energy ---*/
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

void CCentJSTArtComp_Flow::SetResidual(double *val_resconv, double *val_resvisc, 
		double **val_Jacobian_i, double **val_Jacobian_j,
		bool art_diss, CConfig *config) {

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
	//																				+ 0.5*DensityInc_j*(Coord_j[nDim-1]-Coord_i[nDim-1]))
	//																			 /(2.0*config->GetFroude()*config->GetFroude()));

	/*--- Compute mean values of the variables ---*/
	MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
	MeanPressure = 0.5*(Pressure_i + Pressure_j);
	MeanBetaInc2 = 0.5*(BetaInc2_i + BetaInc2_j);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);

	/*--- Get projected flux tensor ---*/
	GetInviscidArtCompProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, Normal, Proj_flux_tensor);

	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = Proj_flux_tensor[iVar];
		val_resvisc[iVar] = 0.0;
	}

	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		GetInviscidArtCompProjJac(MeanDensity, MeanVelocity, MeanBetaInc2, Normal, 0.5, val_Jacobian_i);
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
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS), Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS), Param_p);	
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

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

CCentJST_AdjFlow::CCentJST_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	grid_movement = config->GetGrid_Movement();
	rotating_frame = config->GetRotating_Frame();

	Diff_Psi = new double [nVar]; Diff_Lapl = new double [nVar];
	Und_Lapl_i = new double [nVar]; Und_Lapl_j = new double [nVar];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];
	MeanPhi = new double [nDim];

	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_AdjFlow();
	Param_Kappa_4 = config->GetKappa_4th_AdjFlow();
	implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
}

CCentJST_AdjFlow::~CCentJST_AdjFlow(void) {

	delete [] Diff_Psi; delete [] Diff_Lapl;
	delete [] Und_Lapl_i; delete [] Und_Lapl_j;
	delete [] Velocity_i; delete [] Velocity_j;
	delete [] MeanPhi;
}

void CCentJST_AdjFlow::SetResidual (double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
		double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
		unsigned short local_art_diss, CConfig *config) {

	/*--- Mean Values ---*/
	MeanPsiRho =  0.5*(Psi_i[0]+Psi_j[0]);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanPhi[iDim] =  0.5*(Psi_i[iDim+1]+Psi_j[iDim+1]);
	MeanPsiE =  0.5*(Psi_i[nVar-1]+Psi_j[nVar-1]);

	/*--- Point i convective residual evaluation ---*/
	ProjVelocity_i = 0; ProjPhi = 0; ProjPhi_Vel = 0; sq_vel = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_i[iDim] = U_i[iDim+1] / U_i[0];
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjPhi += MeanPhi[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_i[iDim];
		sq_vel += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	phis1 = ProjPhi + ProjVelocity_i*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_i*MeanPsiE;

	val_resconv_i[0] = ProjVelocity_i*MeanPsiRho - phis2*ProjVelocity_i + Gamma_Minus_One*phis1*sq_vel;
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_i[iDim+1] = ProjVelocity_i*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_i[iDim];
	val_resconv_i[nVar-1] = ProjVelocity_i*MeanPsiE + Gamma_Minus_One*phis1;

	/*--- Rotating Frame ---*/
	if (rotating_frame) {
		double ProjRotVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjRotVel += RotVel_i[iDim]*Normal[iDim];
		ProjRotVel = Rot_Flux;
		val_resconv_i[0] -= ProjRotVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjRotVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjRotVel*MeanPsiE;
	}

	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		val_Jacobian_ii[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[0][jDim+1] = -0.5*ProjVelocity_i*Velocity_i[jDim] + Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_ii[0][nVar-1] = 0.5*ProjVelocity_i*(Gamma_Minus_One*sq_vel - Enthalpy_i);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_ii[iDim+1][0] = 0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_ii[iDim+1][jDim+1] = 0.5*Normal[iDim]*Velocity_i[jDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*Normal[jDim];
			val_Jacobian_ii[iDim+1][iDim+1] += 0.5*ProjVelocity_i;
			val_Jacobian_ii[iDim+1][nVar-1] = 0.5*Enthalpy_i*Normal[iDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*ProjVelocity_i;
		}
		val_Jacobian_ii[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[nVar-1][jDim+1] = 0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_ii[nVar-1][nVar-1] = 0.5*Gamma*ProjVelocity_i;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ij[iVar][jVar] = val_Jacobian_ii[iVar][jVar];

		/*--- Contribution from a rotating frame ---*/
		if (rotating_frame) {
			double ProjRotVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjRotVel += 0.5*RotVel_i[iDim]*Normal[iDim];
			ProjRotVel = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjRotVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjRotVel;
			}
		}
	}


	/*--- Point j convective residual evaluation ---*/
	ProjVelocity_j = 0; ProjPhi_Vel = 0; sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_j[iDim] = U_j[iDim+1] / U_j[0];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_j[iDim];
		sq_vel += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}

	phis1 = ProjPhi + ProjVelocity_j*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_j*MeanPsiE;

	val_resconv_j[0] = -(ProjVelocity_j*MeanPsiRho - phis2*ProjVelocity_j + Gamma_Minus_One*phis1*sq_vel);
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_j[iDim+1] = -(ProjVelocity_j*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_j[iDim]);
	val_resconv_j[nVar-1] = -(ProjVelocity_j*MeanPsiE + Gamma_Minus_One*phis1);

	/*--- Rotating Frame ---*/
	if (rotating_frame) {
		double ProjRotVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjRotVel += RotVel_j[iDim]*Normal[iDim];
		ProjRotVel = Rot_Flux;
		val_resconv_j[0] += ProjRotVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjRotVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjRotVel*MeanPsiE;
	}

	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		val_Jacobian_jj[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[0][jDim+1] = 0.5*ProjVelocity_j*Velocity_j[jDim] - Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_jj[0][nVar-1] = -0.5*ProjVelocity_j*(Gamma_Minus_One*sq_vel - Enthalpy_j);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_jj[iDim+1][0] = -0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_jj[iDim+1][jDim+1] = -0.5*Normal[iDim]*Velocity_j[jDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*Normal[jDim];
			val_Jacobian_jj[iDim+1][iDim+1] -= 0.5*ProjVelocity_j;
			val_Jacobian_jj[iDim+1][nVar-1] = -0.5*Enthalpy_j*Normal[iDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*ProjVelocity_j;
		}
		val_Jacobian_jj[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[nVar-1][jDim+1] = -0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_jj[nVar-1][nVar-1] = -0.5*Gamma*ProjVelocity_j;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ji[iVar][jVar] = val_Jacobian_jj[iVar][jVar];

		/*--- Contribution from a rotating frame ---*/
		if (rotating_frame) {
			double ProjRotVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjRotVel += 0.5*RotVel_j[iDim]*Normal[iDim];
			ProjRotVel = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjRotVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjRotVel;
			}
		}
	}

	/*--- Computes differences btw. variables and Laplacians ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_Psi[iVar]  = Psi_i[iVar]-Psi_j[iVar];
	}

	/*--- Account for effects of rotation. ---*/
	if (rotating_frame) {
		//		double ProjRotVel_i = 0.0; double ProjRotVel_j = 0.0;
		//		for (iDim = 0; iDim < nDim; iDim++) {
		//			ProjRotVel_i += RotVel_i[iDim]*Normal[iDim];
		//			ProjRotVel_j += RotVel_j[iDim]*Normal[iDim];
		//		}
		//		ProjVelocity_i -= ProjRotVel_i;
		//		ProjVelocity_j += ProjRotVel_j;
		ProjVelocity_i -= Rot_Flux;
		ProjVelocity_j += Rot_Flux;
	}

	/*--- Compute the spectral radius and stretching factor ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

	double sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
	Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
	Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;

	/*--- Compute viscous residual 1st- & 3rd-order dissipation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = (Epsilon_2*Diff_Psi[iVar]-Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
		if (implicit) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
		}
	}
}

CCentJSTArtComp_AdjFlow::CCentJSTArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Diff_Psi = new double [nVar]; Diff_Lapl = new double [nVar];
	Und_Lapl_i = new double [nVar]; Und_Lapl_j = new double [nVar];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];
	Proj_Jac_Tensor_i = new double*[nVar];
	Proj_Jac_Tensor_j = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Jac_Tensor_i[iVar] = new double[nVar];
		Proj_Jac_Tensor_j[iVar] = new double[nVar];
	}

	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_AdjFlow();
	Param_Kappa_4 = config->GetKappa_4th_AdjFlow();
	implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

}

CCentJSTArtComp_AdjFlow::~CCentJSTArtComp_AdjFlow(void) {

	delete [] Diff_Psi; delete [] Diff_Lapl;
	delete [] Und_Lapl_i; delete [] Und_Lapl_j;
	delete [] Velocity_i; delete [] Velocity_j;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Proj_Jac_Tensor_i[iVar];
		delete [] Proj_Jac_Tensor_j[iVar];
	}
	delete [] Proj_Jac_Tensor_i;
	delete [] Proj_Jac_Tensor_j;

}

void CCentJSTArtComp_AdjFlow::SetResidual (double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
		double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
		unsigned short local_art_diss, CConfig *config) {

	/*--- Set the variables at point i ---*/
	Pressure_i = U_i[0]; Pressure_j = U_j[0];
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
	}

	/*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidArtCompProjJac(DensityInc_i, Velocity_i, BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
	GetInviscidArtCompProjJac(DensityInc_j, Velocity_j, BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);

	for (iVar = 0; iVar < nDim; iVar++) {
		val_resconv_i[iVar]  = 0.0; val_resconv_j[iVar]  = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_resconv_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
			val_resconv_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
		}
	}

	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_ii[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
				val_Jacobian_ij[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
				val_Jacobian_ji[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
				val_Jacobian_jj[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
			}
	}

	/*--- Computes differences btw. variables and Laplacians ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_Psi[iVar]  = Psi_i[iVar]-Psi_j[iVar];
	}

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

	/*--- Compute the spectral radius and stretching factor ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

	/*--- Compute streching factor ---*/
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
	Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
	Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;

	/*--- Compute viscous residual 1st- & 3rd-order dissipation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = (Epsilon_2*Diff_Psi[iVar]-Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
		if (implicit) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
		}
	}
}

CCentJST_LinFlow::CCentJST_LinFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Param_p = 0.3;
	Param_Kappa_4 = config->GetKappa_4th_LinFlow();

	Diff_DeltaU = new double [nVar];
	Diff_Lapl = new double [nVar];
	Und_Lapl_i = new double [nVar];
	Und_Lapl_j = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	MeanVelocity = new double [nDim];
	MeanDeltaVel = new double [nDim];
	MeanJacobian = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		MeanJacobian[iVar] = new double [nVar];
	Jacobian_i = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Jacobian_i[iVar] = new double [nVar];
	Jacobian_j = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Jacobian_j[iVar] = new double [nVar];
}

CCentJST_LinFlow::~CCentJST_LinFlow(void) {
	unsigned short iVar;

	delete [] Diff_DeltaU;
	delete [] Diff_Lapl;
	delete [] Und_Lapl_i;
	delete [] Und_Lapl_j;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] MeanDeltaVel;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] MeanJacobian[iVar];
	delete [] MeanJacobian;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Jacobian_i[iVar];
	delete [] Jacobian_i;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Jacobian_j[iVar];
	delete [] Jacobian_j;
}

void CCentJST_LinFlow::SetResidual (double *val_resconv, double *val_resvisc, double **val_Jacobian_i, 
		double **val_Jacobian_j, bool art_diss, CConfig *config) {

	/*--- Mean Values of the linealized variables ---*/
	MeanDeltaRho =  0.5*(DeltaU_i[0]+DeltaU_j[0]);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanDeltaVel[iDim] =  0.5*(DeltaU_i[iDim+1]+DeltaU_j[iDim+1]);
	MeanDeltaE =  0.5*(DeltaU_i[nVar-1]+DeltaU_j[nVar-1]);

	/*--- Values of the flow variables at point i ---*/
	Density_i = U_i[0];
	ProjVelocity_i = 0; Area = 0;
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	DensityEnergy_i = U_i[nDim+1];
	Energy_i = DensityEnergy_i / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel));
	Pressure_i = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	Enthalpy_i = (DensityEnergy_i + Pressure_i) / Density_i;

	/*--- Values of the flow variables at point j ---*/
	Density_j = U_j[0];
	ProjVelocity_j = 0; 
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		sq_vel += Velocity_j[iDim]*Velocity_j[iDim]; 
	}
	DensityEnergy_j = U_j[nDim+1];
	Energy_j = DensityEnergy_j / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel));
	Pressure_j = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	Enthalpy_j = (DensityEnergy_j + Pressure_j) / Density_j;

	/*--- Mean values the flow variables ---*/
	MeanDensity = 0.5*(Density_i + Density_j);
	for (iDim = 0; iDim < nDim; iDim++) MeanVelocity[iDim] = 0.5*(Velocity_j[iDim] + Velocity_i[iDim]);
	MeanPressure = 0.5*(Pressure_i + Pressure_j);
	MeanEnthalpy = 0.5*(Enthalpy_j+Enthalpy_i);
	MeanEnergy = 0.5*(Energy_j+Energy_i);

	/*--- Compute projected inviscid jacobian (scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal) ---*/
	GetInviscidProjJac(Velocity_i, Energy_i, Normal, 0.5, Jacobian_i);
	GetInviscidProjJac(Velocity_j, Energy_j, Normal, 0.5, Jacobian_j);

	/*--- Compute inviscid flux $Jacobian x DeltaU$ ---*/

	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++)
			val_resconv[iVar] += Jacobian_i[iVar][jVar] * DeltaU_i[jVar] + Jacobian_j[iVar][jVar] * DeltaU_j[jVar];
	}

	if (!art_diss) return;

	/*--- Computes differences btw. variables and Laplacians ---*/
	for (iVar = 0; iVar < nVar; iVar++) 
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];

	/*--- Calcula el radio espectral local, Factor de stretching factor ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

	Phi_i = pow(0.5*max(0.0,(Lambda_i - MeanLambda)/(MeanLambda+EPS)), Param_p);
	Phi_j = pow(0.5*max(0.0,(Lambda_j - MeanLambda)/(MeanLambda+EPS)), Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

	sc4 = 9.0/(double(Neighbor_i*(1+Neighbor_i))) + 9.0/(double(Neighbor_j*(1+Neighbor_j)));

	Epsilon_4 = Param_Kappa_4*sc4;

	for (iVar = 0; iVar < nVar; iVar++)
		val_resvisc[iVar] = -Epsilon_4*Diff_Lapl[iVar]*StretchingFactor*MeanLambda;	
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

void CCentLax_Flow::SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
		bool art_diss, CConfig *config) {

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
		GetInviscidProjJac(MeanVelocity, MeanEnergy, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}

	/*--- In case there is a grid movement ---*/
	if (grid_movement) {
		ProjVelocity = 0;
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

	/*--- If we are in a rotating frame ---*/
	if (rotating_frame) {
		//		ProjVelocity = 0;
		//		for (iDim = 0; iDim < nDim; iDim++)
		//			ProjVelocity += 0.5*(RotVel_i[iDim]+RotVel_j[iDim])*Normal[iDim];
		ProjVelocity = Rot_Flux;
		for (iVar = 0; iVar < nVar; iVar++) {
			val_resconv[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			if (implicit) {
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}
	}

	if (!art_diss) return;

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
		//		double ProjRotVel_i = 0.0; double ProjRotVel_j = 0.0;
		//		for (iDim = 0; iDim < nDim; iDim++) {
		//			ProjRotVel_i += RotVel_i[iDim]*Normal[iDim];
		//			ProjRotVel_j += RotVel_j[iDim]*Normal[iDim];
		//		}
		//		ProjVelocity_i -= ProjRotVel_i;
		//		ProjVelocity_j -= ProjRotVel_j;
		ProjVelocity_i -= Rot_Flux;
		ProjVelocity_j -= Rot_Flux;
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
	for (iVar = 0; iVar < nVar; iVar++)
		val_resvisc[iVar] = Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;

	if (implicit) {	
		cte = Epsilon_0*StretchingFactor*MeanLambda;

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

void CCentLaxArtComp_Flow::SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
		bool art_diss, CConfig *config) {

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
		GetInviscidArtCompProjJac(MeanDensity, MeanVelocity, MeanBetaInc2, Normal, 0.5, val_Jacobian_i);
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
	}

	if (!art_diss) return;

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

CCentLax_AdjFlow::CCentLax_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Diff_Psi = new double [nVar]; 	MeanPhi = new double [nDim];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];

	implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_AdjFlow();
	rotating_frame = config->GetRotating_Frame();

}

CCentLax_AdjFlow::~CCentLax_AdjFlow(void) {

	delete [] Diff_Psi; delete [] MeanPhi;
	delete [] Velocity_i; delete [] Velocity_j;

}

void CCentLax_AdjFlow::SetResidual (double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
		double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
		unsigned short local_art_diss, CConfig *config) {

	/*--- Mean value of the adjoint variables ---*/
	MeanPsiRho =  0.5*(Psi_i[0]+Psi_j[0]);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanPhi[iDim] =  0.5*(Psi_i[iDim+1]+Psi_j[iDim+1]);
	MeanPsiE =  0.5*(Psi_i[nVar-1]+Psi_j[nVar-1]);

	/*--- Evaluation at point i ---*/	
	ProjVelocity_i = 0; ProjPhi = 0; ProjPhi_Vel = 0; sq_vel = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_i[iDim] = U_i[iDim+1] / U_i[0];
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjPhi += MeanPhi[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_i[iDim];
		sq_vel += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	phis1 = ProjPhi + ProjVelocity_i*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_i*MeanPsiE;

	/*--- Compute inviscid residual at point i ---*/	
	val_resconv_i[0] = ProjVelocity_i*MeanPsiRho - phis2*ProjVelocity_i + Gamma_Minus_One*phis1*sq_vel;
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_i[iDim+1] = ProjVelocity_i*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_i[iDim];
	val_resconv_i[nVar-1] = ProjVelocity_i*MeanPsiE + Gamma_Minus_One*phis1;

	/*--- Rotating Frame ---*/
	if (rotating_frame) {
		double ProjRotVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjRotVel += RotVel_i[iDim]*Normal[iDim];
		ProjRotVel = Rot_Flux;
		val_resconv_i[0] -= ProjRotVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjRotVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjRotVel*MeanPsiE;
	}

	/*--- Inviscid contribution to the implicit part ---*/	
	if (implicit) {
		val_Jacobian_ii[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[0][jDim+1] = -0.5*ProjVelocity_i*Velocity_i[jDim] + Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_ii[0][nVar-1] = 0.5*ProjVelocity_i*(Gamma_Minus_One*sq_vel - Enthalpy_i);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_ii[iDim+1][0] = 0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_ii[iDim+1][jDim+1] = 0.5*Normal[iDim]*Velocity_i[jDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*Normal[jDim];
			val_Jacobian_ii[iDim+1][iDim+1] += 0.5*ProjVelocity_i;
			val_Jacobian_ii[iDim+1][nVar-1] = 0.5*Enthalpy_i*Normal[iDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*ProjVelocity_i;
		}
		val_Jacobian_ii[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[nVar-1][jDim+1] = 0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_ii[nVar-1][nVar-1] = 0.5*Gamma*ProjVelocity_i;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ij[iVar][jVar] = val_Jacobian_ii[iVar][jVar];

		/*--- Contribution from a rotating frame ---*/
		if (rotating_frame) {
			double ProjRotVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjRotVel += 0.5*RotVel_i[iDim]*Normal[iDim];
			ProjRotVel = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjRotVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjRotVel;
			}
		}
	}

	/*--- Evaluation at point j ---*/	
	ProjVelocity_j = 0; ProjPhi_Vel = 0; sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_j[iDim] = U_j[iDim+1] / U_j[0];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_j[iDim];
		sq_vel += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}

	phis1 = ProjPhi + ProjVelocity_j*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_j*MeanPsiE;

	/*--- Compute inviscid residual at point j ---*/	
	val_resconv_j[0] = -(ProjVelocity_j*MeanPsiRho - phis2*ProjVelocity_j + Gamma_Minus_One*phis1*sq_vel);
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_j[iDim+1] = -(ProjVelocity_j*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_j[iDim]);
	val_resconv_j[nVar-1] = -(ProjVelocity_j*MeanPsiE + Gamma_Minus_One*phis1);

	/*--- Rotating Frame ---*/
	if (rotating_frame) {
		double ProjRotVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjRotVel += RotVel_j[iDim]*Normal[iDim];
		ProjRotVel = Rot_Flux;
		val_resconv_j[0] += ProjRotVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjRotVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjRotVel*MeanPsiE;
	}

	/*--- Inviscid contribution to the implicit part ---*/	
	if (implicit) {
		val_Jacobian_jj[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[0][jDim+1] = 0.5*ProjVelocity_j*Velocity_j[jDim] - Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_jj[0][nVar-1] = -0.5*ProjVelocity_j*(Gamma_Minus_One*sq_vel - Enthalpy_j);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_jj[iDim+1][0] = -0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_jj[iDim+1][jDim+1] = -0.5*Normal[iDim]*Velocity_j[jDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*Normal[jDim];
			val_Jacobian_jj[iDim+1][iDim+1] -= 0.5*ProjVelocity_j;
			val_Jacobian_jj[iDim+1][nVar-1] = -0.5*Enthalpy_j*Normal[iDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*ProjVelocity_j;
		}
		val_Jacobian_jj[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[nVar-1][jDim+1] = -0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_jj[nVar-1][nVar-1] = -0.5*Gamma*ProjVelocity_j;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ji[iVar][jVar] = val_Jacobian_jj[iVar][jVar];

		/*--- Contribution from a rotating frame ---*/
		if (rotating_frame) {
			double ProjRotVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjRotVel += 0.5*RotVel_j[iDim]*Normal[iDim];
			ProjRotVel = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjRotVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjRotVel;
			}
		}
	}

	if (local_art_diss == 0) return;

	/*--- Computes differences btw. variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) 
		Diff_Psi[iVar] = Psi_i[iVar]-Psi_j[iVar];

	/*--- Account for effects of rotation. ---*/
	if (rotating_frame) {
		//		double ProjRotVel_i = 0.0; double ProjRotVel_j = 0.0;
		//		for (iDim = 0; iDim < nDim; iDim++) {
		//			ProjRotVel_i += RotVel_i[iDim]*Normal[iDim];
		//			ProjRotVel_j += RotVel_j[iDim]*Normal[iDim];
		//		}
		//		ProjVelocity_i -= ProjRotVel_i;
		//		ProjVelocity_j += ProjRotVel_j;
		ProjVelocity_i -= Rot_Flux;
		ProjVelocity_j += Rot_Flux;
	}

	/*--- Compute spectral radius ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

	/*--- Compute streching factor ---*/
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_0 = Param_Kappa_0*sc2*double(nDim)/3.0;

	/*--- Artifical dissipation evaluation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = Epsilon_0*StretchingFactor*MeanLambda*Diff_Psi[iVar];
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
	}

	/*--- Contribution to implicit part ---*/
	if (implicit) {	
		for (iVar = 0; iVar < nVar; iVar++) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
		}
	}

}

CCentLaxArtComp_AdjFlow::CCentLaxArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Diff_Psi = new double [nVar]; 	MeanPhi = new double [nDim];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];
	Proj_Jac_Tensor_i = new double*[nVar];
	Proj_Jac_Tensor_j = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Jac_Tensor_i[iVar] = new double[nVar];
		Proj_Jac_Tensor_j[iVar] = new double[nVar];
	}

	implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_AdjFlow();

}

CCentLaxArtComp_AdjFlow::~CCentLaxArtComp_AdjFlow(void) {

	delete [] Diff_Psi; delete [] MeanPhi;
	delete [] Velocity_i; delete [] Velocity_j;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Proj_Jac_Tensor_i[iVar];
		delete [] Proj_Jac_Tensor_j[iVar];
	}
	delete [] Proj_Jac_Tensor_i;
	delete [] Proj_Jac_Tensor_j;

}

void CCentLaxArtComp_AdjFlow::SetResidual (double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
		double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
		unsigned short local_art_diss, CConfig *config) {

	/*--- Set the variables at point i ---*/
	Pressure_i = U_i[0]; Pressure_j = U_j[0];
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
	}

	/*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidArtCompProjJac(DensityInc_i, Velocity_i, BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
	GetInviscidArtCompProjJac(DensityInc_j, Velocity_j, BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);

	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv_i[iVar]  = 0.0; val_resconv_j[iVar]  = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_resconv_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*0.5*(Psi_i[jVar]+Psi_j[jVar]);
			val_resconv_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*0.5*(Psi_i[jVar]+Psi_j[jVar]);
		}
	}

	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_ii[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
				val_Jacobian_ij[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
				val_Jacobian_ji[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
				val_Jacobian_jj[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
			}
	}

	/*--- Computes differences btw. variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) 
		Diff_Psi[iVar] = Psi_i[iVar]-Psi_j[iVar];

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

	/*--- Compute spectral radius ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

	/*--- Compute streching factor ---*/
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_0 = Param_Kappa_0*sc2*double(nDim)/3.0;

	/*--- Artifical dissipation evaluation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = Epsilon_0*StretchingFactor*MeanLambda*Diff_Psi[iVar];
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
	}

	/*--- Contribution to implicit part ---*/
	if (implicit) {	
		for (iVar = 0; iVar < nVar; iVar++) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
		}
	}

}

CCentLax_LinFlow::CCentLax_LinFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Artificial Dissipation coefficients ---*/
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_LinFlow();

	Diff_DeltaU = new double [nVar];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	MeanVelocity = new double [nDim];
	MeanDeltaVel = new double [nDim];
	MeanJacobian = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		MeanJacobian[iVar] = new double [nVar];
	Jacobian_i = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Jacobian_i[iVar] = new double [nVar];
	Jacobian_j = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Jacobian_j[iVar] = new double [nVar];
}

CCentLax_LinFlow::~CCentLax_LinFlow(void) {
	unsigned short iVar;

	delete [] Diff_DeltaU;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] MeanVelocity;
	delete [] MeanDeltaVel;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] MeanJacobian[iVar];
	delete [] MeanJacobian;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Jacobian_i[iVar];
	delete [] Jacobian_i;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Jacobian_j[iVar];
	delete [] Jacobian_j;
}

void CCentLax_LinFlow::SetResidual (double *val_resconv, double *val_resvisc, double **val_Jacobian_i,
		double **val_Jacobian_j, bool art_diss, CConfig *config) {

	/*--- Mean Values of the linealized variables ---*/
	MeanDeltaRho =  0.5*(DeltaU_i[0]+DeltaU_j[0]);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanDeltaVel[iDim] =  0.5*(DeltaU_i[iDim+1]+DeltaU_j[iDim+1]);
	MeanDeltaE =  0.5*(DeltaU_i[nVar-1]+DeltaU_j[nVar-1]);

	/*--- Values of the flow variables at point i ---*/
	Density_i = U_i[0];
	ProjVelocity_i = 0; Area = 0;
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	DensityEnergy_i = U_i[nDim+1];
	Energy_i = DensityEnergy_i / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel));
	Pressure_i = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	Enthalpy_i = (DensityEnergy_i + Pressure_i) / Density_i;

	/*--- Values of the flow variables at point j ---*/
	Density_j = U_j[0];
	ProjVelocity_j = 0; sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		sq_vel += Velocity_j[iDim]*Velocity_j[iDim]; }
	DensityEnergy_j = U_j[nDim+1];
	Energy_j = DensityEnergy_j / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel));
	Pressure_j = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	Enthalpy_j = (DensityEnergy_j + Pressure_j) / Density_j;

	/*--- Mean values the flow variables ---*/
	MeanDensity = 0.5*(Density_i + Density_j);
	for (iDim = 0; iDim < nDim; iDim++) MeanVelocity[iDim] = 0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
	MeanPressure = 0.5*(Pressure_i + Pressure_j);
	MeanEnthalpy = 0.5*(Enthalpy_i + Enthalpy_j);
	MeanEnergy = 0.5*(Energy_i + Energy_j);

	/*--- Compute projected inviscid jacobian (scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal) ---*/
	GetInviscidProjJac(Velocity_i, Energy_i, Normal, 0.5, Jacobian_i);
	GetInviscidProjJac(Velocity_j, Energy_j, Normal, 0.5, Jacobian_j);

	/*--- Compute inviscid flux $Jacobian x DeltaU$ ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++)
			val_resconv[iVar] += Jacobian_i[iVar][jVar] * DeltaU_i[jVar] + Jacobian_j[iVar][jVar] * DeltaU_j[jVar];
	}

	if (!art_diss) return;

	/*--- Computes differences btw. variables and dS ---*/
	for (iVar = 0; iVar < nVar; iVar++) 
		Diff_DeltaU[iVar] = DeltaU_i[iVar]-DeltaU_j[iVar];

	/*--- Compute spectral radius Factor de stretching factor ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

	/*--- Compute stretching factor ---*/
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_i = Param_Kappa_0*sc2*double(nDim)/3.0;

	/*--- Evaluate artificial dissipation ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_resvisc[iVar] = Epsilon_i*Diff_DeltaU[iVar]*StretchingFactor*MeanLambda;
}

CUpwRoe_Plasma::~CUpwRoe_Plasma(void) {
	unsigned short iVar, iFluids;


	delete [] Diff_U;

	delete [] Density_i;
	delete [] Energy_i;
	delete [] SoundSpeed_i;
	delete [] Pressure_i;
	delete [] Enthalpy_i;

	delete [] Density_j;
	delete [] Energy_j;
	delete [] SoundSpeed_j;
	delete [] Pressure_j;
	delete [] Enthalpy_j;

	delete [] RoeDensity;
	delete [] RoeEnthalpy;
	delete [] RoeSoundSpeed;

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		delete [] Velocity_i[iFluids];
		delete [] Velocity_j[iFluids];
		delete [] RoeVelocity[iFluids];
		delete [] delta_vel[iFluids];
	}

	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] RoeVelocity;
	delete [] delta_vel[iFluids];

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

CUpwRoe_Plasma::CUpwRoe_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies,unsigned short val_nFluids ,CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies,val_nFluids, config) {

	unsigned short iVar, iFluids;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Diff_U = new double [nVar];

	Density_i		= new double[nFluids];
	Energy_i		= new double[nFluids];
	SoundSpeed_i	= new double[nFluids];
	Pressure_i		= new double[nFluids];
	Enthalpy_i		= new double[nFluids];

	Density_j		= new double[nFluids];
	Energy_j		= new double[nFluids];
	SoundSpeed_j	= new double[nFluids];
	Pressure_j		= new double[nFluids];
	Enthalpy_j		= new double[nFluids];

	RoeDensity		= new double[nFluids];
	RoeEnthalpy		= new double[nFluids];
	RoeSoundSpeed	= new double[nFluids];

	ProjVelocity	= new double[nFluids];
	ProjVelocity_i	= new double[nFluids];
	ProjVelocity_j	= new double[nFluids];

	Velocity_i		= new double* [nFluids];
	Velocity_j		= new double* [nFluids];
	RoeVelocity		= new double* [nFluids];

	delta_vel		= new double* [nFluids];

	for (iFluids = 0; iFluids < nFluids; iFluids++) {
		Velocity_i[iFluids]	= new double [nDim];
		Velocity_j[iFluids]	= new double [nDim];
		RoeVelocity[iFluids]	= new double [nDim];
		delta_vel[iFluids]		= new double [nDim];

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

void CUpwRoe_Plasma::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	unsigned short iFluids,  loc = 0; 
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);                    /*! Area of the face*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/

	/*--- Point i, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	for (iFluids = 0; iFluids < nFluids; iFluids ++) {

		loc = iFluids*(nDim+2);
		Density_i[iFluids]	= U_i[loc + 0];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) { 
			Velocity_i[iFluids][iDim] = U_i[loc + iDim+1] / Density_i[iFluids];
			sq_vel += Velocity_i[iFluids][iDim]*Velocity_i[iFluids][iDim];
		}

		Energy_i[iFluids]		= U_i[loc + nDim+1] / Density_i[iFluids];
		SoundSpeed_i[iFluids]	= sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i[iFluids]-0.5*sq_vel)));
		Pressure_i[iFluids]		= (SoundSpeed_i[iFluids] * SoundSpeed_i[iFluids] * Density_i[iFluids]) / Gamma;
		Enthalpy_i[iFluids]		= (U_i[loc + nDim+1] + Pressure_i[iFluids]) / Density_i[iFluids];

		/*--- Point j, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/

		Density_j[iFluids]		= U_j[loc + 0];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) { 
			Velocity_j[iFluids][iDim] = U_j[loc + iDim+1] / Density_j[iFluids];
			sq_vel += Velocity_j[iFluids][iDim]*Velocity_j[iFluids][iDim];
		}

		Energy_j[iFluids]		= U_j[loc + nDim+1] / Density_j[iFluids];
		SoundSpeed_j[iFluids]	= sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j[iFluids]-0.5*sq_vel)));
		Pressure_j[iFluids]		= (SoundSpeed_j[iFluids] * SoundSpeed_j[iFluids] * Density_j[iFluids]) / Gamma;
		Enthalpy_j[iFluids]		= (U_j[loc + nDim+1] + Pressure_j[iFluids]) / Density_j[iFluids];

		/*--- Promediate Roe variables iPoint and jPoint ---*/
		R = sqrt(fabs(Density_j[iFluids]/Density_i[iFluids]));
		RoeDensity[iFluids] = R*Density_i[iFluids];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) { 
			RoeVelocity[iFluids][iDim] = (R*Velocity_j[iFluids][iDim]+Velocity_i[iFluids][iDim])/(R+1.0);
			sq_vel += RoeVelocity[iFluids][iDim]*RoeVelocity[iFluids][iDim];
		}
		RoeEnthalpy[iFluids] = (R*Enthalpy_j[iFluids]+Enthalpy_i[iFluids])/(R+1);
		RoeSoundSpeed[iFluids] = sqrt(fabs((Gamma-1)*(RoeEnthalpy[iFluids]-0.5*sq_vel)));

	}
	/*--- Compute Proj_flux_tensor_i ---*/
	GetInviscidProjFlux(Density_i, Velocity_i, Pressure_i, Enthalpy_i, Normal, Proj_flux_tensor_i);
	/*--- Compute Proj_flux_tensor_j ---*/
	GetInviscidProjFlux(Density_j, Velocity_j, Pressure_j, Enthalpy_j, Normal, Proj_flux_tensor_j);
	/*--- Compute P and Lambda (do it with the Normal) ---*/
	GetPMatrix(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitaryNormal, P_Tensor);

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		ProjVelocity[iFluids]		= 0.0; 
		ProjVelocity_i[iFluids]		= 0.0; 
		ProjVelocity_j[iFluids]		= 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity[iFluids]   += RoeVelocity[iFluids][iDim]*UnitaryNormal[iDim];
			ProjVelocity_i[iFluids] += Velocity_i[iFluids][iDim] *UnitaryNormal[iDim];
			ProjVelocity_j[iFluids] += Velocity_j[iFluids][iDim] *UnitaryNormal[iDim];	
		}
	}
	/*--- Flow eigenvalues and Entropy correctors ---*/

	for (iFluids = 0; iFluids < nFluids; iFluids++) {

		loc = iFluids*(nDim+2);

		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda[loc + iDim] = ProjVelocity[iFluids];
			Epsilon[loc + iDim] = 4.0*max(EPS, max(Lambda[loc + iDim]-ProjVelocity_i[iFluids],ProjVelocity_j[iFluids]-Lambda[loc + iDim]));
		}

		Lambda[loc + nDim]  = ProjVelocity[iFluids] + RoeSoundSpeed[iFluids];
		Epsilon[loc + nDim] = 4.0*max(EPS, max(Lambda[loc + nDim]-(ProjVelocity_i[iFluids]+SoundSpeed_i[iFluids]),(ProjVelocity_j[iFluids]+SoundSpeed_j[iFluids])-Lambda[loc + nDim]));
		Lambda[loc + nDim+1] = ProjVelocity[iFluids] - RoeSoundSpeed[iFluids];
		Epsilon[loc + nDim+1] = 4.0*max(EPS, max(Lambda[loc + nDim+1]-(ProjVelocity_i[iFluids]-SoundSpeed_i[iFluids]),(ProjVelocity_j[iFluids]-SoundSpeed_j[iFluids])-Lambda[loc + nDim+1]));
	}

	/*--- Entropy correction ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
		else
			Lambda[iVar] = fabs(Lambda[iVar]);

	/*	for (iVar = 0; iVar < nVar; iVar++)
			Lambda[iVar] = fabs(Lambda[iVar]);
	 */

	if (!implicit) {

		/*--- Compute wave amplitudes (characteristics) ---*/
		for (iFluids = 0; iFluids < nFluids; iFluids++) {
			proj_delta_vel = 0.0;
			loc = iFluids * (nDim+2);
			for (iDim = 0; iDim < nDim; iDim++) {
				delta_vel[iFluids][iDim] = Velocity_j[iFluids][iDim] - Velocity_i[iFluids][iDim];
				proj_delta_vel += delta_vel[iFluids][iDim]*Normal[iDim];
			}
			delta_p = Pressure_j[iFluids] - Pressure_i[iFluids];
			delta_rho = Density_j[iFluids] - Density_i[iFluids];
			proj_delta_vel = proj_delta_vel/Area;

			if (nDim == 3) {
				delta_wave[loc+0] = delta_rho - delta_p/(RoeSoundSpeed[iFluids]*RoeSoundSpeed[iFluids]);
				delta_wave[loc+1] = UnitaryNormal[0]*delta_vel[iFluids][2]-UnitaryNormal[2]*delta_vel[iFluids][0];
				delta_wave[loc+2] = UnitaryNormal[1]*delta_vel[iFluids][0]-UnitaryNormal[0]*delta_vel[iFluids][1];
				delta_wave[loc+3] = proj_delta_vel + delta_p/(RoeDensity[iFluids]*RoeSoundSpeed[iFluids]);
				delta_wave[loc+4] = -proj_delta_vel + delta_p/(RoeDensity[iFluids]*RoeSoundSpeed[iFluids]);
			}
			else {
				delta_wave[loc+0] = delta_rho - delta_p/(RoeSoundSpeed[iFluids]*RoeSoundSpeed[iFluids]);
				delta_wave[loc+1] = UnitaryNormal[1]*delta_vel[iFluids][0]-UnitaryNormal[0]*delta_vel[iFluids][1];
				delta_wave[loc+2] = proj_delta_vel + delta_p/(RoeDensity[iFluids]*RoeSoundSpeed[iFluids]);
				delta_wave[loc+3] = -proj_delta_vel + delta_p/(RoeDensity[iFluids]*RoeSoundSpeed[iFluids]);
			}
		}
		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
			for (jVar = 0; jVar < nVar; jVar++) {
				val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
			}
		}
	}

	else {
		/*--- Compute inverse P ---*/
		GetPMatrix_inv(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitaryNormal, invP_Tensor);

		/*--- Jacobias of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
		GetInviscidProjJac(Velocity_i, Energy_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidProjJac(Velocity_j, Energy_j, Normal, 0.5, val_Jacobian_j); 

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

void CUpwRoe_PlasmaDiatomic::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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
			//			SoundSpeed_i[iSpecies]	= sqrt(fabs(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_i[iSpecies]-0.5*Vel2)));
			SoundSpeed_i[iSpecies] 	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_i[iSpecies] - 0.5*Vel2 - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			//			Pressure_i[iSpecies]		= (SoundSpeed_i[iSpecies] * SoundSpeed_i[iSpecies] * Density_i[iSpecies]) / GammaDiatomic;
			Pressure_i[iSpecies] = (GammaDiatomic-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies]);
		}
		else {
			//			SoundSpeed_i[iSpecies]	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_i[iSpecies]-0.5*Vel2));
			SoundSpeed_i[iSpecies] 	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_i[iSpecies] - 0.5*Vel2 - Energy_el_i[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			//			Pressure_i[iSpecies]		= (SoundSpeed_i[iSpecies] * SoundSpeed_i[iSpecies] * Density_i[iSpecies]) / GammaMonatomic;
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
		if (iSpecies < nDiatomics) {
			Energy_vib_j[iSpecies]  = U_j[loc+nDim+2] / Density_j[iSpecies];
			//			SoundSpeed_j[iSpecies]	= sqrt(fabs(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_j[iSpecies]-0.5*Vel2)));
			SoundSpeed_j[iSpecies] 	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_j[iSpecies] - 0.5*Vel2 - Energy_vib_j[iSpecies] - Energy_el_j[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			//			Pressure_j[iSpecies]		= (SoundSpeed_j[iSpecies] * SoundSpeed_j[iSpecies] * Density_j[iSpecies]) / GammaDiatomic;
			Pressure_j[iSpecies] = (GammaDiatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_vib_j[iSpecies] - Energy_el_j[iSpecies]);
		}
		else {
			//			SoundSpeed_j[iSpecies]	= sqrt(fabs(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_j[iSpecies]-0.5*Vel2)));
			SoundSpeed_j[iSpecies] 	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_j[iSpecies] - 0.5*Vel2 - Energy_el_j[iSpecies] - config->GetEnthalpy_Formation(iSpecies)));
			//			Pressure_j[iSpecies]		= (SoundSpeed_j[iSpecies] * SoundSpeed_j[iSpecies] * Density_j[iSpecies]) / GammaMonatomic;
			Pressure_j[iSpecies] = (GammaMonatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 1.0/2.0*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el_j[iSpecies]);
		}
		//		Enthalpy_j[iSpecies]		= (U_j[loc + nDim+1] + Pressure_j[iSpecies]) / Density_j[iSpecies];
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
		if (iSpecies < nDiatomics)
			//			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaDiatomic-1.0)*(RoeEnthalpy[iSpecies]-0.5*Vel2)));
			//			RoeSoundSpeed[iSpecies] = (R*SoundSpeed_j[iSpecies]+SoundSpeed_i[iSpecies])/(R+1.0);
			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaDiatomic-1.0)*(RoeEnthalpy[iSpecies]-0.5*Vel2 - config->GetEnthalpy_Formation(iSpecies))));

		else
			//			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaMonatomic-1.0)*(RoeEnthalpy[iSpecies]-0.5*Vel2)));
			//			RoeSoundSpeed[iSpecies] = (R*SoundSpeed_j[iSpecies]+SoundSpeed_i[iSpecies])/(R+1.0);
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
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = 1.0;
	}

	for (iVar = 0; iVar < nVar; iVar++) {
		Lambda[iVar] = fabs(Lambda[iVar]);
	}

	/*--- Compute Proj_flux_tensor_i ---*/
	GetInviscidProjFlux_(Density_i, Velocity_i, Pressure_i, Enthalpy_i, Energy_vib_i, Normal, Proj_flux_tensor_i);

	/*--- Compute Proj_flux_tensor_j ---*/
	GetInviscidProjFlux_(Density_j, Velocity_j, Pressure_j, Enthalpy_j, Energy_vib_j, Normal, Proj_flux_tensor_j);

	/*--- Compute P and Lambda (do it with the Normal) ---*/
	//	GetPMatrix_(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitaryNormal, P_Tensor);
	GetPMatrix_(RoeDensity, RoeVelocity, RoeEnthalpy, RoeSoundSpeed, Energy_vib_i, Energy_el_i, config, UnitaryNormal, P_Tensor);

	/*	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << "\t" << P_Tensor[iVar][jVar];
		}
		cout << endl;
	}
	cin.get(); */
	/*	cout <<"P_Tensor"<< P_Tensor[0][0] << " " <<P_Tensor[0][1] << " " <<P_Tensor[0][2] << " " <<P_Tensor[0][3] << endl;
	cout <<"P_Tensor"<< P_Tensor[1][0] << " " <<P_Tensor[1][1] << " " <<P_Tensor[1][2] << " " <<P_Tensor[1][3] << endl;
	cout <<"P_Tensor"<< P_Tensor[2][0] << " " <<P_Tensor[2][1] << " " <<P_Tensor[2][2] << " " <<P_Tensor[2][3] << endl;
	cout <<"P_Tensor"<< P_Tensor[3][0] << " " <<P_Tensor[3][1] << " " <<P_Tensor[3][2] << " " <<P_Tensor[3][3] << endl;*/



	/*--- Compute inverse P ---*/
	GetPMatrix_inv_(RoeDensity, RoeVelocity, RoeSoundSpeed, Energy_vib_i, Energy_el_i, config, UnitaryNormal, invP_Tensor);
	/*	cout <<"Pinv_Tensor"<< invP_Tensor[0][0] << " " <<invP_Tensor[0][1] << " " <<invP_Tensor[0][2] << " " <<invP_Tensor[0][3] << endl;
	cout <<"Pinv_Tensor"<< invP_Tensor[1][0] << " " <<invP_Tensor[1][1] << " " <<invP_Tensor[1][2] << " " <<invP_Tensor[1][3] << endl;
	cout <<"Pinv_Tensor"<< invP_Tensor[2][0] << " " <<invP_Tensor[2][1] << " " <<invP_Tensor[2][2] << " " <<invP_Tensor[2][3] << endl;
	cout <<"Pinv_Tensor"<< invP_Tensor[3][0] << " " <<invP_Tensor[3][1] << " " <<invP_Tensor[3][2] << " " <<invP_Tensor[3][3] << endl;*/

	/*--- Jacobias of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidProjJac_(Velocity_i, Energy_i, Normal, 0.5, val_Jacobian_i);
	/*	cout <<"val_Jacobian_i"<< val_Jacobian_i[0][0] << " " <<val_Jacobian_i[0][1] << " " <<val_Jacobian_i[0][2] << " " <<val_Jacobian_i[0][3] << endl;
	cout <<"val_Jacobian_i"<< val_Jacobian_i[1][0] << " " <<val_Jacobian_i[1][1] << " " <<val_Jacobian_i[1][2] << " " <<val_Jacobian_i[1][3] << endl;
	cout <<"val_Jacobian_i"<< val_Jacobian_i[2][0] << " " <<val_Jacobian_i[2][1] << " " <<val_Jacobian_i[2][2] << " " <<val_Jacobian_i[2][3] << endl;
	cout <<"val_Jacobian_i"<< val_Jacobian_i[3][0] << " " <<val_Jacobian_i[3][1] << " " <<val_Jacobian_i[3][2] << " " <<val_Jacobian_i[3][3] << endl;*/


	/*	unsigned short kVar;
	double **tmp1, **tmp2;
	tmp1 = new double *[nVar];
	tmp2 = new double *[nVar]; 
	for (iVar = 0; iVar < nVar; iVar++){
		tmp1[iVar] = new double [nVar];
		tmp2[iVar] = new double [nVar];
	}

	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			tmp1[iVar][jVar] = 0.0;
			for (kVar = 0; kVar < nVar; kVar++) {
				tmp1[iVar][jVar] += invP_Tensor[iVar][kVar] * val_Jacobian_i[kVar][jVar];
			}
		}
	}
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			tmp2[iVar][jVar] = 0.0;
			for (kVar = 0; kVar < nVar; kVar++) {
				tmp2[iVar][jVar] += tmp1[iVar][kVar] * P_Tensor[kVar][jVar];
			}
		}
	}
	cout << endl << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << "\t" << tmp2[iVar][jVar];
		}
		cout << endl;
	}
	cin.get(); */


	GetInviscidProjJac_(Velocity_j, Energy_j, Normal, 0.5, val_Jacobian_j);
	/*	cout <<"val_Jacobian_j"<< val_Jacobian_j[0][0] << " " <<val_Jacobian_j[0][1] << " " <<val_Jacobian_j[0][2] << " " <<val_Jacobian_j[0][3] << endl;
	cout <<"val_Jacobian_j"<< val_Jacobian_j[1][0] << " " <<val_Jacobian_j[1][1] << " " <<val_Jacobian_j[1][2] << " " <<val_Jacobian_j[1][3] << endl;
	cout <<"val_Jacobian_j"<< val_Jacobian_j[2][0] << " " <<val_Jacobian_j[2][1] << " " <<val_Jacobian_j[2][2] << " " <<val_Jacobian_j[2][3] << endl;
	cout <<"val_Jacobian_j"<< val_Jacobian_j[3][0] << " " <<val_Jacobian_j[3][1] << " " <<val_Jacobian_j[3][2] << " " <<val_Jacobian_j[3][3] << endl;


	cout <<"Lambda"<< Lambda[0] << " " <<Lambda[1] << " " <<Lambda[2] << " " <<Lambda[3] << endl; cin.get(); */

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

CUpwRoe_AdjPlasmaDiatomic::CUpwRoe_AdjPlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {

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
	SoundSpeed_i	= new double[nSpecies];
	Pressure_i		= new double[nSpecies];
	Enthalpy_i		= new double[nSpecies];

	Density_j		= new double[nSpecies];
	Energy_j		= new double[nSpecies];
	Energy_vib_j = new double [nSpecies];
	SoundSpeed_j	= new double[nSpecies];
	Pressure_j		= new double[nSpecies];
	Enthalpy_j		= new double[nSpecies];

	RoeDensity		= new double[nSpecies];
	RoeEnthalpy		= new double[nSpecies];
	RoeSoundSpeed	= new double[nSpecies];

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
	Lambda				= new double [nVar];
	Epsilon				= new double [nVar];

	P_Tensor			= new double* [nVar];
	invP_Tensor			= new double* [nVar];
	Proj_Jac_Tensor_i	= new double* [nVar];
	Proj_Jac_Tensor_j	= new double* [nVar];
	Proj_ModJac_Tensor = new double*[nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]  = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
		Proj_Jac_Tensor_i[iVar]	= new double [nVar];
		Proj_Jac_Tensor_j[iVar]	= new double [nVar];
		Proj_ModJac_Tensor[iVar] = new double[nVar];
	}
}

CUpwRoe_AdjPlasmaDiatomic::~CUpwRoe_AdjPlasmaDiatomic(void) {
	unsigned short iVar, iSpecies;


	delete [] Diff_U;

	delete [] Density_i;
	delete [] Energy_i;
	delete [] SoundSpeed_i;
	delete [] Pressure_i;
	delete [] Enthalpy_i;

	delete [] Density_j;
	delete [] Energy_j;
	delete [] SoundSpeed_j;
	delete [] Pressure_j;
	delete [] Enthalpy_j;

	delete [] RoeDensity;
	delete [] RoeEnthalpy;
	delete [] RoeSoundSpeed;

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
	delete [] Lambda;
	delete [] Epsilon;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
		delete [] Proj_Jac_Tensor_i[iVar];
		delete [] Proj_Jac_Tensor_j[iVar];
		delete [] Proj_ModJac_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
	delete [] Proj_Jac_Tensor_i;
	delete [] Proj_Jac_Tensor_j;
	delete [] Proj_ModJac_Tensor;
}

void CUpwRoe_AdjPlasmaDiatomic::SetResidual(double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, 
		double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) {

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
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i[iSpecies];
			sq_vel += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
		}
		Energy_i[iSpecies]		= U_i[loc + nDim+1] / Density_i[iSpecies];
		if (iSpecies < nDiatomics) {
			SoundSpeed_i[iSpecies]	= sqrt(fabs(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_i[iSpecies]-0.5*sq_vel)));
			Pressure_i[iSpecies]		= (SoundSpeed_i[iSpecies] * SoundSpeed_i[iSpecies] * Density_i[iSpecies]) / GammaDiatomic;
			Energy_vib_i[iSpecies]  = U_i[loc + nDim+2];
		}
		else {
			SoundSpeed_i[iSpecies]	= sqrt(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_i[iSpecies]-0.5*sq_vel));
			Pressure_i[iSpecies]		= (SoundSpeed_i[iSpecies] * SoundSpeed_i[iSpecies] * Density_i[iSpecies]) / GammaMonatomic;
		}
		Enthalpy_i[iSpecies]		= (U_i[loc + nDim+1] + Pressure_i[iSpecies]) / Density_i[iSpecies];

		/*--- Point j, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
		Density_j[iSpecies]		= U_j[loc + 0];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_j[iSpecies][iDim] = U_j[loc + iDim+1] / Density_j[iSpecies];
			sq_vel += Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
		}
		Energy_j[iSpecies]		= U_j[loc + nDim+1] / Density_j[iSpecies];
		if (iSpecies < nDiatomics) {
			SoundSpeed_j[iSpecies]	= sqrt(fabs(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_j[iSpecies]-0.5*sq_vel)));
			Pressure_j[iSpecies]		= (SoundSpeed_j[iSpecies] * SoundSpeed_j[iSpecies] * Density_j[iSpecies]) / GammaDiatomic;
			Energy_vib_j[iSpecies]	= U_j[loc + nDim+2];
		}
		else {
			SoundSpeed_j[iSpecies]	= sqrt(fabs(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_j[iSpecies]-0.5*sq_vel)));
			Pressure_j[iSpecies]		= (SoundSpeed_j[iSpecies] * SoundSpeed_j[iSpecies] * Density_j[iSpecies]) / GammaMonatomic;
		}
		Enthalpy_j[iSpecies]		= (U_j[loc + nDim+1] + Pressure_j[iSpecies]) / Density_j[iSpecies];

		/*--- Promediate Roe variables iPoint and jPoint ---*/
		R = sqrt(fabs(Density_j[iSpecies]/Density_i[iSpecies]));
		RoeDensity[iSpecies] = R*Density_i[iSpecies];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			RoeVelocity[iSpecies][iDim] = (R*Velocity_j[iSpecies][iDim]+Velocity_i[iSpecies][iDim])/(R+1.0);
			sq_vel += RoeVelocity[iSpecies][iDim]*RoeVelocity[iSpecies][iDim];
		}
		RoeEnthalpy[iSpecies] = (R*Enthalpy_j[iSpecies]+Enthalpy_i[iSpecies])/(R+1);
		if (iSpecies < nDiatomics)
			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaDiatomic-1.0)*(RoeEnthalpy[iSpecies]-0.5*sq_vel)));
		else
			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaMonatomic-1.0)*(RoeEnthalpy[iSpecies]-0.5*sq_vel)));
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

	/*--- Jacobias of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidProjJac_(Velocity_i, Energy_i, Normal, 0.5, Proj_Jac_Tensor_i);		
	GetInviscidProjJac_(Velocity_j, Energy_j, Normal, 0.5, Proj_Jac_Tensor_j);

	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual_i[iVar] = 0.0; val_residual_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_residual_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
			val_residual_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
		}
	}

	/*--- Compute P, inverse P, and store eigenvalues ---*/
	//	GetPMatrix_(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitaryNormal, P_Tensor);
	double *Energy_el_i;
	Energy_el_i = new double[nSpecies];
	for (iSpecies=0; iSpecies<nSpecies; iSpecies++)
		Energy_el_i[iSpecies] = 0.0;

	GetPMatrix_(RoeDensity, RoeVelocity, RoeEnthalpy, RoeSoundSpeed, Energy_vib_i, Energy_el_i, config, UnitaryNormal, P_Tensor);		

	GetPMatrix_inv_(RoeDensity, RoeVelocity, RoeSoundSpeed, Energy_vib_i, Energy_el_i, config, UnitaryNormal, invP_Tensor);

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (iDim = 0; iDim < nDim; iDim++)
			Lambda[loc + iDim] = ProjVelocity[iSpecies];
		Lambda[loc + nDim]  = ProjVelocity[iSpecies] + RoeSoundSpeed[iSpecies];
		Lambda[loc + nDim+1] = ProjVelocity[iSpecies] - RoeSoundSpeed[iSpecies];
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = 1.0;
	}

	for (iVar = 0; iVar < nVar; iVar++)
		Lambda[iVar] = fabs(Lambda[iVar]);

	/*--- Flux approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) { 
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++) 
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_ModJac_Tensor[iVar][jVar] = 0.5*Proj_ModJac_Tensor_ij;
		}
	}

	for (iVar = 0; iVar < nVar; iVar++)
		for (jVar = 0; jVar < nVar; jVar++) {
			val_residual_i[iVar] -= Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar]);
			val_residual_j[iVar] += Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar]);
		}

	/*--- Implicit contributions, Transpose the matrices and store the Jacobians. Note the negative
	 sign for the ji and jj Jacobians bc the normal direction is flipped. ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) { 
				val_Jacobian_ii[jVar][iVar] = Proj_Jac_Tensor_i[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ij[jVar][iVar] = Proj_Jac_Tensor_i[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ji[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar]);
				val_Jacobian_jj[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar]);
			}
		}
	}
}

CUpwRoe_PlasmaDiatomic_AM::CUpwRoe_PlasmaDiatomic_AM(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {

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
	SoundSpeed_i	= new double[nSpecies];
	Pressure_i		= new double[nSpecies];
	Enthalpy_i		= new double[nSpecies];

	Density_j		= new double[nSpecies];
	Energy_j		= new double[nSpecies];
	Energy_vib_j = new double [nSpecies];
	SoundSpeed_j	= new double[nSpecies];
	Pressure_j		= new double[nSpecies];
	Enthalpy_j		= new double[nSpecies];

	RoeDensity		= new double[nSpecies];
	RoeEnthalpy		= new double[nSpecies];
	RoeVibEnergy	= new double[nDiatomics];
	RoeSoundSpeed	= new double[nSpecies];

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

CUpwRoe_PlasmaDiatomic_AM::~CUpwRoe_PlasmaDiatomic_AM(void) {
	unsigned short iVar, iSpecies;


	delete [] Diff_U;

	delete [] Density_i;
	delete [] Energy_i;
	delete [] SoundSpeed_i;
	delete [] Pressure_i;
	delete [] Enthalpy_i;

	delete [] Density_j;
	delete [] Energy_j;
	delete [] SoundSpeed_j;
	delete [] Pressure_j;
	delete [] Enthalpy_j;

	delete [] RoeDensity;
	delete [] RoeEnthalpy;
	delete [] RoeVibEnergy;

	delete [] RoeSoundSpeed;

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

void CUpwRoe_PlasmaDiatomic_AM::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) { 
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i[iSpecies];
			sq_vel += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
		}
		Energy_i[iSpecies]		= U_i[loc + nDim+1] / Density_i[iSpecies];
		if (iSpecies < nDiatomics) {
			SoundSpeed_i[iSpecies]	= sqrt(fabs(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_i[iSpecies]-0.5*sq_vel)));
			Pressure_i[iSpecies]	= (SoundSpeed_i[iSpecies] * SoundSpeed_i[iSpecies] * Density_i[iSpecies]) / GammaDiatomic;
			Energy_vib_i[iSpecies]  = U_i[loc + nDim+2];
		}	
		else { 
			SoundSpeed_i[iSpecies]	= sqrt(fabs(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_i[iSpecies]-0.5*sq_vel)));
			Pressure_i[iSpecies]	= (SoundSpeed_i[iSpecies] * SoundSpeed_i[iSpecies] * Density_i[iSpecies]) / GammaMonatomic;
		}
		Enthalpy_i[iSpecies]		= (U_i[loc + nDim+1] + Pressure_i[iSpecies]) / Density_i[iSpecies];

		/*--- Point j, Needs to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/

		Density_j[iSpecies]		= U_j[loc + 0];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) { 
			Velocity_j[iSpecies][iDim] = U_j[loc + iDim+1] / Density_j[iSpecies];
			sq_vel += Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
		}
		Energy_j[iSpecies]		= U_j[loc + nDim+1] / Density_j[iSpecies];
		if (iSpecies < nDiatomics) {
			SoundSpeed_j[iSpecies]	= sqrt(fabs(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_j[iSpecies]-0.5*sq_vel)));
			Pressure_j[iSpecies]	= (SoundSpeed_j[iSpecies] * SoundSpeed_j[iSpecies] * Density_j[iSpecies]) / GammaDiatomic;
			Energy_vib_j[iSpecies]	= U_j[loc + nDim+2];
		}
		else {
			SoundSpeed_j[iSpecies]	= sqrt(fabs(GammaMonatomic*(GammaMonatomic-1.0)*(Energy_j[iSpecies]-0.5*sq_vel)));
			Pressure_j[iSpecies]	= (SoundSpeed_j[iSpecies] * SoundSpeed_j[iSpecies] * Density_j[iSpecies]) / GammaMonatomic;
		}
		Enthalpy_j[iSpecies]		= (U_j[loc + nDim+1] + Pressure_j[iSpecies]) / Density_j[iSpecies];

		/*--- Calculate Roe variables iPoint and jPoint ---*/
		R = sqrt(fabs(Density_j[iSpecies]/Density_i[iSpecies]));
		RoeDensity[iSpecies] 		= R*Density_i[iSpecies];
		sq_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) { 
			RoeVelocity[iSpecies][iDim] = (R*Velocity_j[iSpecies][iDim]+Velocity_i[iSpecies][iDim])/(R+1.0);
			sq_vel += RoeVelocity[iSpecies][iDim]*RoeVelocity[iSpecies][iDim];
		}
		RoeEnthalpy[iSpecies] = (R*Enthalpy_j[iSpecies]+Enthalpy_i[iSpecies])/(R+1);

		if (iSpecies < nDiatomics) {
			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaDiatomic-1.0)*(RoeEnthalpy[iSpecies]-0.5*sq_vel)));
			RoeVibEnergy[iSpecies] = (R*Energy_vib_j[iSpecies]+Energy_vib_i[iSpecies])/(R+1);
		}
		else 
			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaMonatomic-1.0)*(RoeEnthalpy[iSpecies]-0.5*sq_vel)));

		/*--- Flow eigenvalues and Entropy correctors ---*/
		ProjVelocity		= 0.0;
		ProjVelocity_i		= 0.0;
		ProjVelocity_j		= 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjVelocity   += RoeVelocity[iSpecies][iDim] *UnitaryNormal[iDim];
			ProjVelocity_i += Velocity_i[iSpecies][iDim] *UnitaryNormal[iDim];
			ProjVelocity_j += Velocity_j[iSpecies][iDim] *UnitaryNormal[iDim];
		}

		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda[loc + iDim] = ProjVelocity;
			Epsilon[loc + iDim] = 4.0*max(0.0, max(Lambda[loc + iDim]-ProjVelocity_i,ProjVelocity_j-Lambda[loc + iDim]));
		}

		Lambda[loc + nDim]  = ProjVelocity + RoeSoundSpeed[iSpecies];
		Epsilon[loc + nDim] = 4.0*max(0.0, max(Lambda[loc + nDim]-(ProjVelocity_i+SoundSpeed_i[iSpecies]),(ProjVelocity_j+SoundSpeed_j[iSpecies])-Lambda[loc + nDim]));
		Lambda[loc + nDim+1] = ProjVelocity - RoeSoundSpeed[iSpecies];
		Epsilon[loc + nDim+1] = 4.0*max(0.0, max(Lambda[loc + nDim+1]-(ProjVelocity_i-SoundSpeed_i[iSpecies]),(ProjVelocity_j-SoundSpeed_j[iSpecies])-Lambda[loc + nDim+1]));
		if (iSpecies < nDiatomics) {
			Lambda[loc + nDim+2] = ProjVelocity;
			Epsilon[loc + nDim+2] = 4.0*max(0.0, max(Lambda[loc + nDim+2]-ProjVelocity_i,ProjVelocity_j-Lambda[loc + nDim+2]));
		}
	}

	for (iVar = 0; iVar < nVar; iVar++) {
		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
		else
			Lambda[iVar] = fabs(Lambda[iVar]);
	}

	if (!implicit) {

		/*--- Compute wave amplitudes (characteristics) ---*/
		for (iFluids = 0; iFluids < nFluids; iFluids++) {
			proj_delta_vel = 0.0;
			loc = iFluids * (nDim+2);
			for (iDim = 0; iDim < nDim; iDim++) {
				delta_vel[iFluids][iDim] = Velocity_j[iFluids][iDim] - Velocity_i[iFluids][iDim];
				proj_delta_vel += delta_vel[iFluids][iDim]*Normal[iDim];
			}
			delta_p = Pressure_j[iFluids] - Pressure_i[iFluids];
			delta_rho = Density_j[iFluids] - Density_i[iFluids];
			proj_delta_vel = proj_delta_vel/Area;

			if (nDim == 3) {
				delta_wave[loc+0] = delta_rho - delta_p/(RoeSoundSpeed[iFluids]*RoeSoundSpeed[iFluids]);
				delta_wave[loc+1] = UnitaryNormal[0]*delta_vel[iFluids][2]-UnitaryNormal[2]*delta_vel[iFluids][0];
				delta_wave[loc+2] = UnitaryNormal[1]*delta_vel[iFluids][0]-UnitaryNormal[0]*delta_vel[iFluids][1];
				delta_wave[loc+3] = proj_delta_vel + delta_p/(RoeDensity[iFluids]*RoeSoundSpeed[iFluids]);
				delta_wave[loc+4] = -proj_delta_vel + delta_p/(RoeDensity[iFluids]*RoeSoundSpeed[iFluids]);
			}
			else {
				delta_wave[loc+0] = delta_rho - delta_p/(RoeSoundSpeed[iFluids]*RoeSoundSpeed[iFluids]);
				delta_wave[loc+1] = UnitaryNormal[1]*delta_vel[iFluids][0]-UnitaryNormal[0]*delta_vel[iFluids][1];
				delta_wave[loc+2] = proj_delta_vel + delta_p/(RoeDensity[iFluids]*RoeSoundSpeed[iFluids]);
				delta_wave[loc+3] = -proj_delta_vel + delta_p/(RoeDensity[iFluids]*RoeSoundSpeed[iFluids]);
			}
		}
		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
			for (jVar = 0; jVar < nVar; jVar++) {
				val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
			}
		}
	}
	else {
		/*--- Compute Proj_flux_tensor_i ---*/
		GetInviscidProjFlux_(Density_i, Velocity_i, Pressure_i, Enthalpy_i, Energy_vib_i, Normal, Proj_flux_tensor_i);
		/*--- Compute Proj_flux_tensor_j ---*/
		GetInviscidProjFlux_(Density_j, Velocity_j, Pressure_j, Enthalpy_j, Energy_vib_j, Normal, Proj_flux_tensor_j);
		/*--- Compute P and Lambda (do it with the Normal) ---*/
		GetPMatrix_AM(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitaryNormal, P_Tensor, RoeVibEnergy);


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
}


CConvective_Template::CConvective_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {


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

CConvective_Template::~CConvective_Template(void) {
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
void CConvective_Template::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		/*!< \brief Normal: Normal vector, it norm is the area of the face. */
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);                    /*! Area of the face*/

	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/

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

	/*--- Flow eigenvalues and Entropy correctors ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Lambda[iDim] = ProjVelocity;
		Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,ProjVelocity_j-Lambda[iDim]));
	}
	Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
	Epsilon[nVar-2] = 4.0*max(0.0, max(Lambda[nVar-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar-2]));
	Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
	Epsilon[nVar-1] = 4.0*max(0.0, max(Lambda[nVar-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar-1]));

	/*--- Entropy correction ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
		else
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

		if (nDim == 3) {
			delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
			delta_wave[1] = UnitaryNormal[0]*delta_vel[2]-UnitaryNormal[2]*delta_vel[0];
			delta_wave[2] = UnitaryNormal[1]*delta_vel[0]-UnitaryNormal[0]*delta_vel[1];
			delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
			delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
		}
		else {
			delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
			delta_wave[1] = UnitaryNormal[1]*delta_vel[0]-UnitaryNormal[0]*delta_vel[1];
			delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
			delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
		}

		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
			for (jVar = 0; jVar < nVar; jVar++)
				val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
		}
	}
	else {

		/*--- Compute inverse P ---*/
		GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, invP_Tensor);

		/*--- Jacobias of the inviscid flux, scale = 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
		GetInviscidProjJac(Velocity_i, Energy_i, Normal, 0.5, val_Jacobian_i);
		GetInviscidProjJac(Velocity_j, Energy_j, Normal, 0.5, val_Jacobian_j);

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
}
