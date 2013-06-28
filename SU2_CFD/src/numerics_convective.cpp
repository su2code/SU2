/*!
 * \file numerics_convective.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.3
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

void CUpwRoe_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CUpwRoe_Flow::SetResidual
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i *U_j
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *Normal Gamma Gamma_Minus_One
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim nVar

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS INT SCALAR iDim iVar jVar
	//SU2_CPP2C VARS DOUBLE SCALAR Area sq_vel
	//SU2_CPP2C VARS DOUBLE SCALAR Density_i Energy_i SoundSpeed_i Pressure_i Enthalpy_i
	//SU2_CPP2C VARS DOUBLE SCALAR Density_j Energy_j SoundSpeed_j Pressure_j Enthalpy_j
	//SU2_CPP2C VARS DOUBLE SCALAR R RoeDensity RoeEnthalpy RoeSoundSpeed
	//SU2_CPP2C VARS DOUBLE SCALAR ProjVelocity ProjVelocity_i ProjVelocity_j
	//SU2_CPP2C VARS DOUBLE SCALAR proj_delta_vel delta_p delta_rho
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nDim UnitaryNormal Velocity_i Velocity_j RoeVelocity delta_vel
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nVar Lambda delta_wave Proj_flux_tensor_i Proj_flux_tensor_j
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nVar SIZE=nVar P_Tensor
	//SU2_CPP2C DECL_LIST END

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
	if (Pressure_i < 0.0) Pressure_i = Pressure_Old_i;

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
	if (Pressure_j < 0.0) Pressure_j = Pressure_Old_j;

	SoundSpeed_j = sqrt(Pressure_j*Gamma/Density_j);
	Enthalpy_j   = (U_j[nDim+1] + Pressure_j) / Density_j;

	/*--- Roe-averaged variables at interface between i & j ---*/
	R = sqrt(abs(Density_j/Density_i));
	RoeDensity = R*Density_i;
	sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
		sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
	}
	RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
	RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));

	/*--- Compute Proj_flux_tensor_i ---*/
	//SU2_CPP2C SUBROUTINE START GetInviscidProjFlux
	GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, Proj_flux_tensor_i);
	//SU2_CPP2C SUBROUTINE LOCATION numerics_structure.cpp
	//SU2_CPP2C SUBROUTINE VARS Density_i Velocity_i Pressure_i Enthalpy_i Proj_flux_tensor_i Normal
	//SU2_CPP2C SUBROUTINE END

	/*--- Compute Proj_flux_tensor_j ---*/
	//SU2_CPP2C SUBROUTINE START GetInviscidProjFlux
	GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, Proj_flux_tensor_j);
	//SU2_CPP2C SUBROUTINE LOCATION numerics_structure.cpp
	//SU2_CPP2C SUBROUTINE VARS Density_j Velocity_j Pressure_j Enthalpy_j Proj_flux_tensor_j Normal
	//SU2_CPP2C SUBROUTINE END

	/*--- Compute P and Lambda (do it with the Normal) ---*/
	//SU2_CPP2C SUBROUTINE START GetPMatrix
	GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitaryNormal, P_Tensor);
	//SU2_CPP2C SUBROUTINE LOCATION numerics_structure.cpp
	//SU2_CPP2C SUBROUTINE VARS RoeDensity RoeVelocity RoeSoundSpeed P_Tensor UnitaryNormal
	//SU2_CPP2C SUBROUTINE END

	ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity   += RoeVelocity[iDim]*UnitaryNormal[iDim];
		ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
	}

	//SU2_CPP2C COMMENT START
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
	//SU2_CPP2C COMMENT END


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

	//SU2_CPP2C COMMENT START
	if (!implicit) {
		//SU2_CPP2C COMMENT END

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

		//SU2_CPP2C COMMENT START

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
		//SU2_CPP2C COMMENT END
		//SU2_CPP2C COMMENT START
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
	//SU2_CPP2C COMMENT END

	//SU2_CPP2C END CUpwRoe_Flow::SetResidual
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

void CUpwRoePrim_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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
	R = sqrt(abs(Density_j/Density_i));
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

void CUpwRoe_Turkel_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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

	/*--- First few flow eigenvalues of A.Normal with the normal---*/
	for (iDim = 0; iDim < nDim; iDim++)
		Lambda[iDim] = ProjVelocity;

	local_Mach = sqrt(sq_vel)/RoeSoundSpeed;
	Beta 	   = max(Beta_min,min(local_Mach,Beta_max));
	Beta2 	   = Beta*Beta;

	one_m_Betasqr 		   = 1.0 - Beta2;  // 1-Beta*Beta
	one_p_Betasqr 		   = 1.0 + Beta2;  // 1+ Beta*Beta
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

	/*--- Compute face area ---*/
	Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);

  /*--- Compute and unitary normal vector ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		UnitaryNormal[iDim] = Normal[iDim]/Area;
    if (UnitaryNormal[iDim] == 0.0) UnitaryNormal[iDim] = EPS;
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

CUpwRoe_AdjFlow::CUpwRoe_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();

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

	/*--- Contribution to velocity projection due to a rotating frame ---*/
	if (rotating_frame) {
		double ProjRotVel = Rot_Flux;
		Q -= ProjRotVel;
	}

	/*--- Contribution to velocity projection due to grid movement ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		Q -= ProjGridVel;
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

	/*--- Flux contribution due to a rotating frame ---*/
	if (rotating_frame) {
		double ProjVelocity = Rot_Flux;
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual_i[iVar] -= ProjVelocity * 0.5*(Psi_i[iVar]+Psi_j[iVar]);
			val_residual_j[iVar] += ProjVelocity * 0.5*(Psi_i[iVar]+Psi_j[iVar]);
		}
	}

	/*--- Flux contribution due to grid movement ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual_i[iVar] -= ProjGridVel * 0.5*(Psi_i[iVar]+Psi_j[iVar]);
			val_residual_j[iVar] += ProjGridVel * 0.5*(Psi_i[iVar]+Psi_j[iVar]);
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
		GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Proj_flux_tensor_i);
		GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Proj_flux_tensor_j);

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

		/*--- Jacobian contributions for a rotating frame ---*/
		if (rotating_frame) {
			double ProjVelocity = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				/*--- Adjust Jacobian main diagonal ---*/
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjVelocity;
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjVelocity;
			}
		}

		/*--- Jacobian contribution due to grid movement ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				/*--- Adjust Jacobian main diagonal ---*/
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjGridVel;
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjGridVel;
			}
		}

	}
}

CUpwRoe_AdjDiscFlow::CUpwRoe_AdjDiscFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

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

	Velocity_id = new double [nDim];
	Velocity_jd = new double [nDim];
	RoeVelocityd = new double [nDim];
	delta_veld  = new double [nDim];
	delta_waved = new double [nVar];
	Proj_flux_tensor_id = new double [nVar];
	Proj_flux_tensor_jd = new double [nVar];
	Lambdad = new double [nVar];
	Epsilond = new double [nVar];
	P_Tensord = new double* [nVar];
	invP_Tensord = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensord[iVar] = new double [nVar];
		invP_Tensord[iVar] = new double [nVar];
	}

	val_residual = new double [nVar];
	val_residuald = new double [nVar];
    
    U_id = new double [nVar];
	U_jd = new double [nVar];



}

CUpwRoe_AdjDiscFlow::~CUpwRoe_AdjDiscFlow(void) {
	unsigned short iVar;

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

	delete [] Velocity_id;
	delete [] Velocity_jd;
	delete [] RoeVelocityd;
	delete [] delta_veld;
	delete [] delta_waved;
	delete [] Proj_flux_tensor_id;
	delete [] Proj_flux_tensor_jd;
	delete [] Lambdad;
	delete [] Epsilond;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensord[iVar];
		delete [] invP_Tensord[iVar];
	}
	delete [] P_Tensord;
	delete [] invP_Tensord;

	delete [] val_residual;
	delete [] val_residuald;
    
    delete [] U_id;
    delete [] U_jd;

}

void CUpwRoe_AdjDiscFlow::SetResidual(double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	//unsigned short iVar, jVar;
	unsigned short iPos, jPos; // iVar, jVar are global and used in SetResidual_ad, so cannot be used here

	// U_i sensitivity:

	for (iPos = 0; iPos < nVar; iPos++){

		for (jPos = 0; jPos < nVar; jPos++){
			U_id[jPos] = 0.0;
			U_jd[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		U_id[iPos] = 1.0;

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			// Transpose each block: [jPos][iPos] -> [iPos][jPos]
			val_Jacobian_i[iPos][jPos] = val_residuald[jPos];

		}

	}

	//	 U_j sensitivity:

	for (iPos = 0; iPos < nVar; iPos++){

		for (jPos = 0; jPos < nVar; jPos++){
			U_id[jPos] = 0.0;
			U_jd[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		U_jd[iPos] = 1.0;

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			// Transpose each block: [jPos][iPos] -> [iPos][jPos]
			val_Jacobian_j[iPos][jPos] = val_residuald[jPos];

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
	GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
	GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);

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

CUpwRoeArtComp_AdjDiscFlow::CUpwRoeArtComp_AdjDiscFlow() {

}

CUpwRoeArtComp_AdjDiscFlow::~CUpwRoeArtComp_AdjDiscFlow(void) {

}

void CUpwRoeArtComp_AdjDiscFlow::SetResidual () {

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

CUpwLin_LevelSet::CUpwLin_LevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];

}

CUpwLin_LevelSet::~CUpwLin_LevelSet(void) {

	delete [] Velocity_i;
	delete [] Velocity_j;

}

void CUpwLin_LevelSet::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                   double **val_JacobianMeanFlow_i, double **val_JacobianMeanFlow_j, CConfig *config) {
	unsigned short iDim;
	double a0, a1, q_ij, Velocity_i[3], Velocity_j[3]; //, dqij_dvi[3], dqij_dvj[3], dabsqij_dvi[3], dabsqij_dvj[3], da0_dvi[3], da0_dvj[3], da1_dvi[3], da1_dvj[3];
  
	q_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = V_i[iDim+1]; //U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] =  V_i[iDim+1]; //U_j[iDim+1]/DensityInc_j;
		q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
	}
  
	a0 = 0.5*(q_ij+fabs(q_ij)); a1 = 0.5*(q_ij-fabs(q_ij));
  
	val_residual[0] = a0*LevelSetVar_i[0]+a1*LevelSetVar_j[0];
  
	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_j[0][0] = a1;
    
//    for (iDim = 0; iDim < nDim; iDim++) {
//      dqij_dvi[iDim] = 0.5 * Normal[iDim]/DensityInc_i;
//      dqij_dvj[iDim] = 0.5 * Normal[iDim]/DensityInc_j;
//      if ( q_ij >= 0.0 ) {
//        dabsqij_dvi[iDim] = dqij_dvi[iDim];
//        dabsqij_dvj[iDim] = dqij_dvj[iDim];
//      }
//      else {
//        dabsqij_dvi[iDim] = -dqij_dvi[iDim];
//        dabsqij_dvj[iDim] = -dqij_dvj[iDim];
//      }
//      da0_dvi[iDim] = 0.5 * (dqij_dvi[iDim] + dabsqij_dvi[iDim]);
//      da1_dvi[iDim] = 0.5 * (dqij_dvi[iDim] - dabsqij_dvi[iDim]);
//      
//      da0_dvj[iDim] = 0.5 * (dqij_dvj[iDim] + dabsqij_dvj[iDim]);
//      da1_dvj[iDim] = 0.5 * (dqij_dvj[iDim] - dabsqij_dvj[iDim]);
//    }
//    
//    val_JacobianMeanFlow_i[0][0] = 0.0; // No pressure contribution;
//    for (iDim = 0; iDim < nDim; iDim++) {
//      val_JacobianMeanFlow_i[0][iDim+1] = da0_dvi[iDim]*LevelSetVar_i[0]+da1_dvi[iDim]*LevelSetVar_j[0];
//      val_JacobianMeanFlow_j[0][iDim+1] = da0_dvj[iDim]*LevelSetVar_i[0]+da1_dvj[iDim]*LevelSetVar_j[0];
//    }
    
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

CUpwLin_AdjDiscLevelSet::CUpwLin_AdjDiscLevelSet() {

}

CUpwLin_AdjDiscLevelSet::~CUpwLin_AdjDiscLevelSet(void) {


}

void CUpwLin_AdjDiscLevelSet::SetResidual ()  {

}

CUpwLin_TransLM::CUpwLin_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement  = config->GetGrid_Movement();
	incompressible = config->GetIncompressible();

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];

}

CUpwLin_TransLM::~CUpwLin_TransLM(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwLin_TransLM::SetResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {


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
	val_residual[0] = a0*TransVar_i[0]+a1*TransVar_j[0];
	val_residual[1] = a0*TransVar_i[1]+a1*TransVar_j[1];
//	cout << "Velicity x: " << Velocity_i[0] << ", " << Velocity_j[0] << endl;
//	cout << "Velicity y: " << Velocity_i[1] << ", " << Velocity_j[1] << endl;
//	cout << "val_resid: " << val_residual[0] << ", " << val_residual[1] << endl;


	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_i[1][1] = a0;
	}
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

CUpwLin_AdjDiscTurb::CUpwLin_AdjDiscTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

}

CUpwLin_AdjDiscTurb::~CUpwLin_AdjDiscTurb(void) {

}

void CUpwLin_AdjDiscTurb::SetResidual () {


}

CUpwLin_AdjDiscTurbSA::CUpwLin_AdjDiscTurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

}

CUpwLin_AdjDiscTurbSA::~CUpwLin_AdjDiscTurbSA(void) {

}

void CUpwLin_AdjDiscTurbSA::SetResidual (double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {


}

CUpwSca_TurbSA::CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement  = config->GetGrid_Movement();
	incompressible = config->GetIncompressible();

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

	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CUpwSca_TurbSA::SetResidual
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i *U_j *TurbVar_i *TurbVar_j
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *Normal
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS INT SCALAR iDim
	//SU2_CPP2C VARS DOUBLE SCALAR Density_i Density_j q_ij a0 a1
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nDim Velocity_i Velocity_j
	//SU2_CPP2C DECL_LIST END

	//SU2_CPP2C COMMENT START
	if (incompressible) {
		Density_i = DensityInc_i;
		Density_j = DensityInc_j;
	}
	else {

		//SU2_CPP2C COMMENT END
		Density_i = U_i[0];
		Density_j = U_j[0];
		//SU2_CPP2C COMMENT START
	}
	//SU2_CPP2C COMMENT END

	q_ij = 0;
	//SU2_CPP2C COMMENT START
	if (rotating_frame) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - RotVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - RotVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	} else if (grid_movement) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - GridVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - GridVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	} else {
	//SU2_CPP2C COMMENT END
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i;
			Velocity_j[iDim] = U_j[iDim+1]/Density_j;
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	//SU2_CPP2C COMMENT START
	}
	//SU2_CPP2C COMMENT END

	a0 = 0.5*(q_ij+fabs(q_ij));
	a1 = 0.5*(q_ij-fabs(q_ij));
	val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];

	//SU2_CPP2C COMMENT START
	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_j[0][0] = a1;
	}
	//SU2_CPP2C COMMENT END

	//SU2_CPP2C END CUpwSca_TurbSA::SetResidual

}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement  = config->GetGrid_Movement();
	incompressible = config->GetIncompressible();

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

	if (incompressible) {
		Density_i = DensityInc_i;
		Density_j = DensityInc_j;
	}
	else {
		Density_i = U_i[0];
		Density_j = U_j[0];
	}

	q_ij = 0;
	if (rotating_frame) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - RotVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - RotVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	} else if (grid_movement) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - GridVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - GridVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	} else {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i;
			Velocity_j[iDim] = U_j[iDim+1]/Density_j;
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
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

CUpwSca_TransLM::CUpwSca_TransLM(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
}

CUpwSca_TransLM::~CUpwSca_TransLM(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwSca_TransLM::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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
	val_residual[0] = a0*TransVar_i[0]+a1*TransVar_j[0];
	val_residual[1] = a0*TransVar_i[1]+a1*TransVar_j[1];

	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_j[0][0] = a1;
		val_Jacobian_i[1][1] = a0;
		val_Jacobian_j[1][1] = a1;
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

CUpwSca_AdjDiscTurb::CUpwSca_AdjDiscTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config)  {

}

CUpwSca_AdjDiscTurb::~CUpwSca_AdjDiscTurb(void) {

}

void CUpwSca_AdjDiscTurb::SetResidual () {


}

CUpwSca_AdjDiscTurbSA::CUpwSca_AdjDiscTurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config)  {
    
    unsigned short nFlowVar = nDim + 2;
    
    TurbVar_id = new double[nVar];
    TurbVar_jd = new double[nVar];
    val_residual = new double [nVar];
    val_residuald = new double [nVar];
    Velocity_i = new double [nDim];
    Velocity_j = new double [nDim];
    Velocity_id = new double [nDim];
    Velocity_jd = new double [nDim];
    
    U_id = new double [nFlowVar];
	U_jd = new double [nFlowVar];

}

CUpwSca_AdjDiscTurbSA::~CUpwSca_AdjDiscTurbSA(void) {
    
    delete [] TurbVar_i;
    delete [] TurbVar_j;
    delete [] val_residual;
    delete [] val_residuald;
    delete [] Velocity_i;
    delete [] Velocity_j;
    delete [] Velocity_id;
    delete [] Velocity_jd;
    
    delete [] U_id;
    delete [] U_jd;

}

void CUpwSca_AdjDiscTurbSA::SetResidual (double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	//unsigned short iVar, jVar;
	unsigned short iPos, jPos; // iVar, jVar are global and used in SetResidual_ad, so cannot be used here

	unsigned short nTotalVar, nFlowVar;
	nFlowVar = nDim + 2;
	nTotalVar = nVar + nFlowVar;

	// U_i sensitivity:

	for (iPos = 0; iPos < nFlowVar; iPos++){

		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
			U_jd[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
			TurbVar_jd[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		U_id[iPos] = 1.0;
        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//			cout << "U_id: " << U_id[jPos] << endl;
//			cout << "U_jd: " << U_jd[jPos] << endl;
//		}

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			// Transpose each block: [jPos][iPos] -> [iPos][jPos]
			val_Jacobian_i[iPos][jPos] = val_residuald[jPos];

		}

//        // TEST AD *****************
//        for (jPos = 0; jPos < nVar; jPos++) {
//            cout << "AD: " << val_residuald[jPos] << endl;
//        }
//        cout << "--" << endl;
//        // AD *****************
	}
    
//    // TEST FD *****************
//    double *temp_U_i, *temp_U_j;
//    double *temp_TurbVar_i, *temp_TurbVar_j;
//    temp_U_i = new double[nFlowVar];
//    temp_U_j = new double[nFlowVar];
//    temp_TurbVar_i  = new double[nVar];
//    temp_TurbVar_j =  new double[nVar];
//    double *temp1_val_residual, *temp2_val_residual;
//    temp1_val_residual = new double[nVar];
//    temp2_val_residual = new double[nVar];
//    
//    double delta;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_U_i[jPos] = U_i[jPos];
//        temp_U_j[jPos] = U_j[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//        temp_TurbVar_j[jPos] = TurbVar_j[jPos];
//    }
//    
//    for (iPos = 0; iPos < nFlowVar; iPos++){
//        
//		for (jPos = 0; jPos < nFlowVar; jPos++){
//			U_i[jPos] = temp_U_i[jPos];
//			U_j[jPos] = temp_U_j[jPos];
//		}
//        
//		for (jPos = 0; jPos < nVar; jPos++){
//			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//			TurbVar_j[jPos] = temp_TurbVar_j[jPos];
//		}
//        
//        if (fabs(temp_U_i[iPos] + temp_U_j[iPos]) > 1e-15)
//            delta = 0.01*(temp_U_i[iPos] + temp_U_j[iPos])/2;
//        else
//            delta = 1e-15;
//        
//        U_i[iPos] = temp_U_i[iPos] - delta;
//                
//		this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//			U_i[jPos] = temp_U_i[jPos];
//			U_j[jPos] = temp_U_j[jPos];
//		}
//        
//		for (jPos = 0; jPos < nVar; jPos++){
//			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//			TurbVar_j[jPos] = temp_TurbVar_j[jPos];
//		}
//        
//        U_i[iPos] = temp_U_i[iPos] + delta;
//        
//        
//        
//		this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//        
//        cout << "U_i: " << temp_U_i[iPos] << endl;
//        
//        for (jPos = 0; jPos < nVar; jPos++) {
//            cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//        }
//	}
//    delete [] temp_U_i;
//    delete [] temp_U_j;
//    delete [] temp_TurbVar_i;
//    delete [] temp_TurbVar_j;
//    delete [] temp1_val_residual;
//    delete [] temp2_val_residual;
//    cin.get();
//    // FD *****************

	// TurbVar_i sensitivity

	for (iPos = 0; iPos < nVar; iPos++){

		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
			U_jd[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
			TurbVar_jd[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		TurbVar_id[iPos] = 1.0;

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			// Transpose each block: [jPos][iPos] -> [iPos][jPos]
			val_Jacobian_i[iPos+nFlowVar][jPos] = val_residuald[jPos];

		}
        
//      // TEST AD *****************
//        for (jPos = 0; jPos < nVar; jPos++) {
//                    cout << "AD: " << val_residuald[jPos] << endl;
//                }
//                cout << "--" << endl;
//      // AD *********

	}
    
//    //    // TEST FD *****************
//        double *temp_U_i, *temp_U_j;
//        double *temp_TurbVar_i, *temp_TurbVar_j;
//        temp_U_i = new double[nFlowVar];
//        temp_U_j = new double[nFlowVar];
//        temp_TurbVar_i  = new double[nVar];
//        temp_TurbVar_j =  new double[nVar];
//        double *temp1_val_residual, *temp2_val_residual;
//        temp1_val_residual = new double[nVar];
//        temp2_val_residual = new double[nVar];
//    
//        double delta;
//    
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//            temp_U_i[jPos] = U_i[jPos];
//            temp_U_j[jPos] = U_j[jPos];
//        }
//    
//        for (jPos = 0; jPos < nVar; jPos++){
//            temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//            temp_TurbVar_j[jPos] = TurbVar_j[jPos];
//        }
//    
//        for (iPos = 0; iPos < nVar; iPos++){
//    
//    		for (jPos = 0; jPos < nFlowVar; jPos++){
//    			U_i[jPos] = temp_U_i[jPos];
//    			U_j[jPos] = temp_U_j[jPos];
//    		}
//    
//    		for (jPos = 0; jPos < nVar; jPos++){
//    			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//    			TurbVar_j[jPos] = temp_TurbVar_j[jPos];
//    		}
//    
//            if (fabs(temp_TurbVar_i[iPos] + temp_TurbVar_j[iPos]) > 1e-15)
//                delta = 0.01*(temp_TurbVar_i[iPos] + temp_TurbVar_j[iPos])/2;
//            else
//                delta = 1e-15;
//    
//            TurbVar_i[iPos] = temp_TurbVar_i[iPos] - delta;
//    
//    		this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//    
//            for (jPos = 0; jPos < nFlowVar; jPos++){
//    			U_i[jPos] = temp_U_i[jPos];
//    			U_j[jPos] = temp_U_j[jPos];
//    		}
//    
//    		for (jPos = 0; jPos < nVar; jPos++){
//    			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//    			TurbVar_j[jPos] = temp_TurbVar_j[jPos];
//    		}
//    
//            TurbVar_i[iPos] = temp_TurbVar_i[iPos] + delta;
//    
//    
//    
//    		this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//    
//           // cout << "U_i: " << temp_U_i[iPos] << endl;
//    
//            for (jPos = 0; jPos < nVar; jPos++) {
//                cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//            }
//    	}
//        delete [] temp_U_i;
//        delete [] temp_U_j;
//        delete [] temp_TurbVar_i;
//        delete [] temp_TurbVar_j;
//        delete [] temp1_val_residual;
//        delete [] temp2_val_residual;
//        cin.get();
//        // FD ****

	// U_j sensitivity:
    
//    double *temp_U_i, *temp_U_j;
//    temp_U_i = new double[nFlowVar];
//    temp_U_j = new double[nFlowVar];
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_U_i[jPos] = U_j[jPos];
//        temp_U_j[jPos] = U_i[jPos];
//    }
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//            U_i[jPos] = temp_U_i[jPos];
//            U_j[jPos] = temp_U_j[jPos];
//        }

	for (iPos = 0; iPos < nFlowVar; iPos++){

		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
			U_jd[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
			TurbVar_jd[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		U_jd[iPos] = 1.0;

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			// Transpose each block: [jPos][iPos] -> [iPos][jPos]
			val_Jacobian_j[iPos][jPos] = val_residuald[jPos];

		}
        
        //        // TEST AD *****************
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//			cout << U_jd[jPos] << endl;
//		}
//                for (jPos = 0; jPos < nVar; jPos++) {
//                    cout << "AD: " << val_residuald[jPos] << endl;
//                }
               //cout << "--" << endl;
        //        // AD *****************

	}
//    cout << "--" << endl;
    
//    //    // TEST FD *****************
//    //    double *temp_U_i, *temp_U_j;
//        double *temp_TurbVar_i, *temp_TurbVar_j;
//    //    temp_U_i = new double[nFlowVar];
//     //   temp_U_j = new double[nFlowVar];
//        temp_TurbVar_i  = new double[nVar];
//        temp_TurbVar_j =  new double[nVar];
//        double *temp1_val_residual, *temp2_val_residual;
//        temp1_val_residual = new double[nVar];
//        temp2_val_residual = new double[nVar];
//    
//        double delta;
//    
////        for (jPos = 0; jPos < nFlowVar; jPos++){
////            temp_U_i[jPos] = U_j[jPos];
////            temp_U_j[jPos] = U_i[jPos];
////        }
//    
//        for (jPos = 0; jPos < nVar; jPos++){
//            temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//            temp_TurbVar_j[jPos] = TurbVar_j[jPos];
//        }
//    
//        for (iPos = 0; iPos < nFlowVar; iPos++){
//    
//    		for (jPos = 0; jPos < nFlowVar; jPos++){
//    			U_i[jPos] = temp_U_i[jPos];
//    			U_j[jPos] = temp_U_j[jPos];
//    		}
//    
//    		for (jPos = 0; jPos < nVar; jPos++){
//    			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//    			TurbVar_j[jPos] = temp_TurbVar_j[jPos];
//    		}
//    
//            if (fabs(temp_U_j[iPos] + temp_U_j[iPos]) > 1e-15)
//                delta = 0.01*(temp_U_j[iPos] + temp_U_j[iPos])/2;
//            else
//                delta = 1e-15;
//    
//            U_j[iPos] = temp_U_j[iPos] - delta;
//    
//    		this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//    
//            for (jPos = 0; jPos < nFlowVar; jPos++){
//    			U_i[jPos] = temp_U_i[jPos];
//    			U_j[jPos] = temp_U_j[jPos];
//    		}
//    
//    		for (jPos = 0; jPos < nVar; jPos++){
//    			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//    			TurbVar_j[jPos] = temp_TurbVar_j[jPos];
//    		}
//    
//            U_j[iPos] = temp_U_j[iPos] + delta;
//    
//    
//    
//    		this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//    
//            cout << "U_j: " << temp_U_j[iPos] << endl;
//    
//            for (jPos = 0; jPos < nVar; jPos++) {
//                cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//            }
//    	}
//        delete [] temp_U_i;
//        delete [] temp_U_j;
//        delete [] temp_TurbVar_i;
//        delete [] temp_TurbVar_j;
//        delete [] temp1_val_residual;
//        delete [] temp2_val_residual;
//        cin.get();
//    //    // FD *****************

	// TurbVar_j sensitivity

	for (iPos = 0; iPos < nVar; iPos++){

		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
			U_jd[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
			TurbVar_jd[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		TurbVar_jd[iPos] = 1.0;

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			// Transpose each block: [jPos][iPos] -> [iPos][jPos]
			val_Jacobian_j[iPos+nFlowVar][jPos] = val_residuald[jPos];

		}
        
//        // TEST AD *****************
//        for (jPos = 0; jPos < nVar; jPos++) {
//            cout << "AD: " << val_residuald[jPos] << endl;
//        }
//        cout << "--" << endl;
//        // AD *********

	}
    
//    //    // TEST FD *****************
//    double *temp_U_i, *temp_U_j;
//    double *temp_TurbVar_i, *temp_TurbVar_j;
//    temp_U_i = new double[nFlowVar];
//    temp_U_j = new double[nFlowVar];
//    temp_TurbVar_i  = new double[nVar];
//    temp_TurbVar_j =  new double[nVar];
//    double *temp1_val_residual, *temp2_val_residual;
//    temp1_val_residual = new double[nVar];
//    temp2_val_residual = new double[nVar];
//    
//    double delta;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_U_i[jPos] = U_i[jPos];
//        temp_U_j[jPos] = U_j[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//        temp_TurbVar_j[jPos] = TurbVar_j[jPos];
//    }
//    
//    for (iPos = 0; iPos < nVar; iPos++){
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//            U_i[jPos] = temp_U_i[jPos];
//            U_j[jPos] = temp_U_j[jPos];
//        }
//        
//        for (jPos = 0; jPos < nVar; jPos++){
//            TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//            TurbVar_j[jPos] = temp_TurbVar_j[jPos];
//        }
//        
//        if (fabs(temp_TurbVar_i[iPos] + temp_TurbVar_j[iPos]) > 1e-15)
//            delta = 0.01*(temp_TurbVar_i[iPos] + temp_TurbVar_j[iPos])/2;
//        else
//            delta = 1e-15;
//        
//        TurbVar_j[iPos] = temp_TurbVar_j[iPos] - delta;
//        
//        this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//            U_i[jPos] = temp_U_i[jPos];
//            U_j[jPos] = temp_U_j[jPos];
//        }
//        
//        for (jPos = 0; jPos < nVar; jPos++){
//            TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//            TurbVar_j[jPos] = temp_TurbVar_j[jPos];
//        }
//        
//        TurbVar_j[iPos] = temp_TurbVar_j[iPos] + delta;
//        
//        
//        
//        this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//        
//        // cout << "U_i: " << temp_U_i[iPos] << endl;
//        
//        for (jPos = 0; jPos < nVar; jPos++) {
//            cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//        }
//    }
//    delete [] temp_U_i;
//    delete [] temp_U_j;
//    delete [] temp_TurbVar_i;
//    delete [] temp_TurbVar_j;
//    delete [] temp1_val_residual;
//    delete [] temp2_val_residual;
//    cin.get();
//    // FD ****
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

void CCentJSTArtComp_Flow::SetResidual(double *val_resconv, double *val_resvisc, 
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
		CConfig *config) {

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

	/*--- Flux contributions due to a rotating frame at point i ---*/
	if (rotating_frame) {
		double ProjRotVel = Rot_Flux;
		val_resconv_i[0] -= ProjRotVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjRotVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjRotVel*MeanPsiE;
	}

	/*--- Flux contributions due to grid movement at point i (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_i[0] -= ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjGridVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjGridVel*MeanPsiE;
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

		/*--- Jacobian contributions due to a rotating frame at point i ---*/
		if (rotating_frame) {
			double ProjRotVel = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjRotVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjRotVel;
			}
		}

		/*--- Jacobian contributions due to grid movement at point i (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjGridVel;
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

	/*--- Flux contributions due to a rotating frame at point j ---*/
	if (rotating_frame) {
		double ProjRotVel = Rot_Flux;
		val_resconv_j[0] += ProjRotVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjRotVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjRotVel*MeanPsiE;
	}

	/*--- Flux contributions due to grid movement at point j (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_j[0] += ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjGridVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjGridVel*MeanPsiE;
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

		/*--- Jacobian contributions due to a rotating frame at point j ---*/
		if (rotating_frame) {
			double ProjRotVel = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjRotVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjRotVel;
			}
		}

		/*--- Jacobian contributions due to grid movement at point j (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjGridVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjGridVel;
			}
		}
	}

	/*--- Computes differences btw. variables and Laplacians ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_Psi[iVar]  = Psi_i[iVar]-Psi_j[iVar];
	}

	/*--- Adjustment to projected velocity due to a rotating frame ---*/
	if (rotating_frame) {
		ProjVelocity_i -= Rot_Flux;
		ProjVelocity_j += Rot_Flux;
	}

	/*--- Adjustment to projected velocity due to mesh motion (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0; double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*Normal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*Normal[iDim];
		}
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j += ProjGridVel;
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

CCentJST_AdjDiscFlow::CCentJST_AdjDiscFlow() {

}

CCentJST_AdjDiscFlow::~CCentJST_AdjDiscFlow(void) {

}

void CCentJST_AdjDiscFlow::SetResidual () {

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
		CConfig *config) {

	/*--- Set the variables at point i ---*/
	Pressure_i = U_i[0]; Pressure_j = U_j[0];
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
	}

	/*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
	GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);

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
		double **val_Jacobian_j, CConfig *config) {

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
	GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Jacobian_i);
	GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Jacobian_j);

	/*--- Compute inviscid flux $Jacobian x DeltaU$ ---*/

	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++)
			val_resconv[iVar] += Jacobian_i[iVar][jVar] * DeltaU_i[jVar] + Jacobian_j[iVar][jVar] * DeltaU_j[jVar];
	}

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

void CCentLaxArtComp_Flow::SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j, 
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

CCentLax_AdjFlow::CCentLax_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Diff_Psi = new double [nVar]; 	MeanPhi = new double [nDim];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];

	implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();

	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_AdjFlow();

}

CCentLax_AdjFlow::~CCentLax_AdjFlow(void) {

	delete [] Diff_Psi; delete [] MeanPhi;
	delete [] Velocity_i; delete [] Velocity_j;

}

void CCentLax_AdjFlow::SetResidual (double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
		double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
		CConfig *config) {

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

	/*--- Flux contributions due to a rotating frame at point i ---*/
	if (rotating_frame) {
		double ProjRotVel = Rot_Flux;
		val_resconv_i[0] -= ProjRotVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjRotVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjRotVel*MeanPsiE;
	}

	/*--- Flux contributions due to grid movement at point i (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_i[0] -= ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjGridVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjGridVel*MeanPsiE;
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

		/*--- Jacobian contributions due to a rotating frame at point i ---*/
		if (rotating_frame) {
			double ProjRotVel = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjRotVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjRotVel;
			}
		}

		/*--- Jacobian contributions due to grid movement at point i (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjGridVel;
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

	/*--- Flux contributions due to a rotating frame at point j ---*/
	if (rotating_frame) {
		double ProjRotVel = Rot_Flux;
		val_resconv_j[0] += ProjRotVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjRotVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjRotVel*MeanPsiE;
	}

	/*--- Flux contributions due to grid movement at point j (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_j[0] += ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjGridVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjGridVel*MeanPsiE;
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

		/*--- Jacobian contributions due to a rotating frame at point j ---*/
		if (rotating_frame) {
			double ProjRotVel = Rot_Flux;
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjRotVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjRotVel;
			}
		}

		/*--- Jacobian contributions due to grid movement at point j (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjGridVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjGridVel;
			}
		}
	}

	/*--- Computes differences btw. variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) 
		Diff_Psi[iVar] = Psi_i[iVar]-Psi_j[iVar];

	/*--- Adjustment to projected velocity due to a rotating frame ---*/
	if (rotating_frame) {
		ProjVelocity_i -= Rot_Flux;
		ProjVelocity_j += Rot_Flux;
	}

	/*--- Adjustment to projected velocity due to mesh motion (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0; double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*Normal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*Normal[iDim];
		}
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j += ProjGridVel;
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

CCentLax_AdjDiscFlow::CCentLax_AdjDiscFlow() {

}

CCentLax_AdjDiscFlow::~CCentLax_AdjDiscFlow(void) {

}

void CCentLax_AdjDiscFlow::SetResidual () {

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
		CConfig *config) {

	/*--- Set the variables at point i ---*/
	Pressure_i = U_i[0]; Pressure_j = U_j[0];
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
	}

	/*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
	GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);

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

CCentLaxArtComp_AdjDiscFlow::CCentLaxArtComp_AdjDiscFlow() {

}

CCentLaxArtComp_AdjDiscFlow::~CCentLaxArtComp_AdjDiscFlow(void) {

}

void CCentLaxArtComp_AdjDiscFlow::SetResidual () {

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
		double **val_Jacobian_j, CConfig *config) {

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
	GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Jacobian_i);
	GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Jacobian_j);

	/*--- Compute inviscid flux $Jacobian x DeltaU$ ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_resconv[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++)
			val_resconv[iVar] += Jacobian_i[iVar][jVar] * DeltaU_i[jVar] + Jacobian_j[iVar][jVar] * DeltaU_j[jVar];
	}

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

void CUpwRoe_Plasma::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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
		if (Pressure_i < 0.0) Pressure_i = Pressure_Old_i_MS[iSpecies];

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
		if (Pressure_j < 0.0) Pressure_j= Pressure_Old_j_MS[iSpecies];

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

void CUpwRoe_Turkel_Plasma::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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
		if (Pressure_i < 0.0) Pressure_i = Pressure_Old_i_MS[iSpecies];

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
		if (Pressure_j < 0.0) Pressure_j= Pressure_Old_j_MS[iSpecies];

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

void CUpwSW_PlasmaDiatomic::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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

void CUpwMSW_PlasmaDiatomic::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
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

void CUpwHLLC_PlasmaDiatomic::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

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



CUpwRoe_AdjPlasmaDiatomic::CUpwRoe_AdjPlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {

	unsigned short iVar, iSpecies;  

	nMonatomics = val_nMonatomics;
	nDiatomics  = val_nDiatomics;

	GammaMonatomic = config->GetGammaMonatomic();
	GammaDiatomic = config->GetGammaDiatomic();

	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Diff_U = new double [nVar];

	Density_i	  	= new double[nSpecies];
	Energy_i	  	= new double[nSpecies];
	Energy_el_i   = new double[nSpecies];
	Energy_vib_i  = new double[nSpecies];
	SoundSpeed_i	= new double[nSpecies];
	Pressure_i		= new double[nSpecies];
	Enthalpy_i		= new double[nSpecies];

	Density_j	  	= new double[nSpecies];
	Energy_j	  	= new double[nSpecies];
	Energy_el_j	  = new double[nSpecies];
	Energy_vib_j  = new double[nSpecies];
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
	Lambda				= new double [nVar];
	Epsilon				= new double [nVar];

	P_Tensor			= new double* [nVar];
	invP_Tensor			= new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]  = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
	}

	Proj_Jac_Tensor_i = new double*[nVar];
	Proj_Jac_Tensor_j = new double*[nVar];
	Proj_ModJac_Tensor = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Jac_Tensor_i[iVar] = new double[nVar];
		Proj_Jac_Tensor_j[iVar] = new double[nVar];
		Proj_ModJac_Tensor[iVar] = new double[nVar];
	}

}

CUpwRoe_AdjPlasmaDiatomic::~CUpwRoe_AdjPlasmaDiatomic(void) {

	unsigned short iVar, iSpecies;


	delete [] Diff_U;

	delete [] Density_i;
	delete [] Energy_i;
	delete [] Energy_el_i;
	delete [] Energy_vib_i;
	delete [] SoundSpeed_i;
	delete [] Pressure_i;
	delete [] Enthalpy_i;

	delete [] Density_j;
	delete [] Energy_j;
	delete [] Energy_el_j;
	delete [] Energy_vib_j;
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
	delete [] Lambda;
	delete [] Epsilon;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Proj_Jac_Tensor_i[iVar];
		delete [] Proj_Jac_Tensor_j[iVar];
		delete [] Proj_ModJac_Tensor[iVar];
	}
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

	/*--- Point i, Need to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		Density_i[iSpecies]	= U_i[loc + 0];
		double Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i[iSpecies];
			Vel2 += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
		}
		Energy_i[iSpecies]		= U_i[loc+nDim+1] / Density_i[iSpecies];
		Energy_vib_i[iSpecies] = 0.0;
		Energy_el_i[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_i[iSpecies]  = U_i[loc + nDim+2] / Density_i[iSpecies];
			SoundSpeed_i[iSpecies]	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_i[iSpecies]-0.5*Vel2-Energy_vib_i[iSpecies]-Energy_el_i[iSpecies]-config->GetEnthalpy_Formation(iSpecies)));
			Pressure_i[iSpecies]		= (GammaDiatomic-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies]-0.5*Vel2-Energy_vib_i[iSpecies]-Energy_el_i[iSpecies]-config->GetEnthalpy_Formation(iSpecies));
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
		Energy_vib_j[iSpecies] = 0.0;
		Energy_el_j[iSpecies] = 0.0;
		if (iSpecies < nDiatomics) {
			Energy_vib_j[iSpecies]	= U_j[loc + nDim+2]/Density_j[iSpecies];
			SoundSpeed_j[iSpecies]	= sqrt(GammaDiatomic*(GammaDiatomic-1.0)*(Energy_j[iSpecies]-0.5*Vel2-Energy_vib_j[iSpecies]-Energy_el_j[iSpecies]-config->GetEnthalpy_Formation(iSpecies)));
			Pressure_j[iSpecies]		= (GammaDiatomic-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies]-0.5*Vel2-Energy_vib_j[iSpecies]-Energy_el_j[iSpecies]-config->GetEnthalpy_Formation(iSpecies));
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
			RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaDiatomic-1.0)*(RoeEnthalpy[iSpecies] - 0.5*Vel2 - RoeEnergy_vib[iSpecies] - config->GetEnthalpy_Formation(iSpecies))));
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
		for (iDim = 0; iDim < nDim; iDim++)
			Lambda[loc + iDim] = ProjVelocity[iSpecies];
		Lambda[loc + nDim]  = ProjVelocity[iSpecies] + RoeSoundSpeed[iSpecies];
		Lambda[loc + nDim+1] = ProjVelocity[iSpecies] - RoeSoundSpeed[iSpecies];
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = ProjVelocity[iSpecies];
	}

	for (iVar = 0; iVar < nVar; iVar++)
		Lambda[iVar] = fabs(Lambda[iVar]);

	/*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidProjJac_(Velocity_i, Energy_i, Energy_vib_i, Enthalpy_i,  Normal, 0.5, Proj_Jac_Tensor_i, config);
	GetInviscidProjJac_(Velocity_j, Energy_j, Energy_vib_j, Enthalpy_j,  Normal, 0.5, Proj_Jac_Tensor_j, config);

	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual_i[iVar] = 0.0; val_residual_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_residual_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
			val_residual_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
		}
	}

	/*--- Compute P and Lambda (do it with the Normal) ---*/
	GetPMatrix_(RoeDensity, RoeVelocity, RoeEnthalpy, RoeSoundSpeed, RoeEnergy_vib, Energy_el_i, config, UnitaryNormal, P_Tensor);

	/*--- Compute inverse P ---*/
	GetPMatrix_inv_(RoeDensity, RoeVelocity, RoeSoundSpeed, RoeEnergy_vib, Energy_el_i, config, UnitaryNormal, invP_Tensor);

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

	/*	cout << "Residual i: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << val_residual_i[iVar] << endl;
	cout << endl << "************" << endl;
	cout << "Residual j: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << val_residual_j[iVar] << endl;
	cin.get();

	cout << endl << "************" << endl;	
	cout << "Jacobian_ii: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ii[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_ij: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ij[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_ji: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ji[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_jj: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_jj[iVar][jVar] << "\t";
		}
		cout << endl;
	}
	cin.get();*/

}

CUpwRoe_AdjDiscPlasmaDiatomic::CUpwRoe_AdjDiscPlasmaDiatomic() {

}

CUpwRoe_AdjDiscPlasmaDiatomic::~CUpwRoe_AdjDiscPlasmaDiatomic(void) {

}

void CUpwRoe_AdjDiscPlasmaDiatomic::SetResidual() {

}

CUpwSW_AdjPlasmaDiatomic::CUpwSW_AdjPlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {

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

	Proj_ModJac_Tensor = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Proj_ModJac_Tensor[iVar] = new double [nVar];

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

CUpwSW_AdjPlasmaDiatomic::~CUpwSW_AdjPlasmaDiatomic(void) {
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
		delete [] Proj_ModJac_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
	delete [] Proj_ModJac_Tensor;
}

void CUpwSW_AdjPlasmaDiatomic::SetResidual(double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, 
		double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) {

	unsigned short iSpecies,  loc = 0;
	Area = 0;

	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual_i[iVar] = 0.0;
		val_residual_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_ii[iVar][jVar] = 0.0;
			val_Jacobian_ij[iVar][jVar] = 0.0;
			val_Jacobian_ji[iVar][jVar] = 0.0;
			val_Jacobian_jj[iVar][jVar] = 0.0;
		}
	}

	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];

	Area = sqrt(Area);                    /*! Area of the face*/

	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/

	/*--- Point i, Need to recompute SoundSpeed / Pressure / Enthalpy in case of 2nd order reconstruction ---*/
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

	///////////////////////////////////////////////////////////////////////////////////////////////////////
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

	/*--- Compute P & invP at i ---*/
	GetPMatrix_(Density_i, Velocity_i, Enthalpy_i, SoundSpeed_i, Energy_vib_i, Energy_el_i, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_i, Velocity_i, SoundSpeed_i, Energy_vib_i, Energy_el_i, config, UnitaryNormal, invP_Tensor);

	/*--- Projected Jacobian (A+) at i ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_ModJac_Tensor[iVar][jVar] = Proj_ModJac_Tensor_ij;
		}
	}

	/*--- Projected flux (i->j) = transpose(A+_i) * Psi_i ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_flux_tensor_i[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_flux_tensor_i[iVar] += Proj_ModJac_Tensor[jVar][iVar]*Psi_i[jVar]*Area;
			val_Jacobian_ii[iVar][jVar] = Proj_ModJac_Tensor[jVar][iVar]*Area;
		}
		val_residual_i[iVar] += Proj_flux_tensor_i[iVar];
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
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = 0.5*(ProjVelocity_j[iSpecies] - fabs(ProjVelocity_j[iSpecies]));
	}

	/*--- Compute P & invP at j ---*/
	GetPMatrix_(Density_j, Velocity_j, Enthalpy_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_j, Velocity_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, invP_Tensor);

	/*--- Projected Jacobian (A-) at j ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_ModJac_Tensor[iVar][jVar] = Proj_ModJac_Tensor_ij;
		}
	}

	/*--- Projected flux = transpose(A-_j) * Psi_j ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_flux_tensor_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {			
			Proj_flux_tensor_j[iVar] += Proj_ModJac_Tensor[jVar][iVar]*Psi_j[jVar]*Area;
			val_Jacobian_ij[iVar][jVar] = Proj_ModJac_Tensor[jVar][iVar]*Area;
		}
		val_residual_i[iVar] += Proj_flux_tensor_j[iVar];
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	/*--- j->i Calculation ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = -UnitaryNormal[iDim];

	/*--- Flow eigenvalues at (Lambda-) i --- */
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);		
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda[loc + iDim] = 0.5*(-ProjVelocity_i[iSpecies] - fabs(-ProjVelocity_i[iSpecies]));
		}
		Lambda[loc + nDim]   = 0.5*(-ProjVelocity_i[iSpecies] + SoundSpeed_i[iSpecies]
		                                                                     - fabs(-ProjVelocity_i[iSpecies] + SoundSpeed_i[iSpecies]));
		Lambda[loc + nDim+1] = 0.5*(-ProjVelocity_i[iSpecies] - SoundSpeed_i[iSpecies]
		                                                                     - fabs(-ProjVelocity_i[iSpecies] - SoundSpeed_i[iSpecies]));
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = 0.5*(-ProjVelocity_i[iSpecies] - fabs(-ProjVelocity_i[iSpecies]));
	}

	/*--- Compute P & invP at i ---*/
	GetPMatrix_(Density_i, Velocity_i, Enthalpy_i, SoundSpeed_i, Energy_vib_i, Energy_el_i, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_i, Velocity_i, SoundSpeed_i, Energy_vib_i, Energy_el_i, config, UnitaryNormal, invP_Tensor);

	/*--- Projected Jacobian (A) at i ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_ModJac_Tensor[iVar][jVar] = Proj_ModJac_Tensor_ij;
		}
	}

	/*--- Projected flux (j->i) = transpose(A-_i) * Psi_i ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_flux_tensor_i[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_flux_tensor_i[iVar] += Proj_ModJac_Tensor[jVar][iVar]*Psi_i[jVar]*Area;
			val_Jacobian_ji[iVar][jVar] = Proj_ModJac_Tensor[jVar][iVar]*Area;
		}
		val_residual_j[iVar] += Proj_flux_tensor_i[iVar];
	}

	/*--- Flow eigenvalues at j (Lambda+) --- */
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);		
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda[loc + iDim] = 0.5*(-ProjVelocity_j[iSpecies] + fabs(-ProjVelocity_j[iSpecies]));
		}
		Lambda[loc + nDim]   = 0.5*(-ProjVelocity_j[iSpecies] + SoundSpeed_j[iSpecies]
		                                                                     + fabs(-ProjVelocity_j[iSpecies] + SoundSpeed_j[iSpecies]));
		Lambda[loc + nDim+1] = 0.5*(-ProjVelocity_j[iSpecies] - SoundSpeed_j[iSpecies]
		                                                                     + fabs(-ProjVelocity_j[iSpecies] - SoundSpeed_j[iSpecies]));
		if (iSpecies < nDiatomics) Lambda[loc + nDim+2] = 0.5*(-ProjVelocity_j[iSpecies] + fabs(-ProjVelocity_j[iSpecies]));
	}

	/*--- Compute P & invP at j ---*/
	GetPMatrix_(Density_j, Velocity_j, Enthalpy_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_j, Velocity_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, invP_Tensor);

	/*--- Projected Jacobian (A-) at j ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
			Proj_ModJac_Tensor[iVar][jVar] = Proj_ModJac_Tensor_ij;
		}
	}

	/*--- Projected flux = transpose(A+_j) * Psi_j ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_flux_tensor_j[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {			
			Proj_flux_tensor_j[iVar] += Proj_ModJac_Tensor[jVar][iVar]*Psi_j[jVar]*Area;
			val_Jacobian_jj[iVar][jVar] = Proj_ModJac_Tensor[jVar][iVar]*Area;
		}
		val_residual_j[iVar] += Proj_flux_tensor_j[iVar];
	}
	/*	cout << "Proj Flux j-i: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << Proj_flux_tensor_i[iVar] << endl;
	cout << endl << "************" << endl;
	cout << "Proj Flux j-j: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << Proj_flux_tensor_j[iVar] << endl;
	cin.get();*/

	///////////////////////////////////////////////////////////////////////////////////////////////////////

	/*	cout << "Residual i: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << val_residual_i[iVar] << endl;
	cout << endl << "************" << endl;
	cout << "Residual j: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << val_residual_j[iVar] << endl;
	cin.get();

	cout << endl << "************" << endl;	
	cout << "Jacobian_ii: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ii[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_ij: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ij[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_ji: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ji[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_jj: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_jj[iVar][jVar] << "\t";
		}
		cout << endl;
	}
	cin.get();*/

	/*--- Reset to i->j ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = -UnitaryNormal[iDim];
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

void CCentJST_Plasma::SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
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

void CCentJST_PlasmaDiatomic::SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
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

void CCentLax_PlasmaDiatomic::SetResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, double **val_Jacobian_j,
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

CCentLax_AdjPlasmaDiatomic::CCentLax_AdjPlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {

	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Artifical dissipation part ---*/
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_AdjPlasma();

	/*--- Allocate some structures ---*/
	Diff_Psi = new double [nVar];
	Velocity_i = new double* [nSpecies];
	Velocity_j = new double* [nSpecies];
	MeanVelocity = new double* [nSpecies];
	for(iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		MeanVelocity[iSpecies] = new double [nDim];
		Velocity_i[iSpecies] = new double [nDim];
		Velocity_j[iSpecies] = new double [nDim];		
	}

	Proj_Jac_Tensor_i = new double*[nVar];
	Proj_Jac_Tensor_j = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Jac_Tensor_i[iVar] = new double[nVar];
		Proj_Jac_Tensor_j[iVar] = new double[nVar];
	}

	Pressure_i = new double [nSpecies];
	Pressure_j = new double [nSpecies];

	MeanEnergy = new double [nSpecies];
	Energy_i = new double [nSpecies];
	Energy_j = new double [nSpecies];
	Energy_vib_i = new double [nSpecies];
	Energy_vib_j = new double [nSpecies];
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

CCentLax_AdjPlasmaDiatomic::~CCentLax_AdjPlasmaDiatomic(void) {

	delete [] Diff_Psi;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		delete [] MeanVelocity[iSpecies];
		delete [] Velocity_i[iSpecies];
		delete [] Velocity_j[iSpecies];
	}

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Proj_Jac_Tensor_i[iVar];
		delete [] Proj_Jac_Tensor_j[iVar];
	}
	delete [] Proj_Jac_Tensor_i;
	delete [] Proj_Jac_Tensor_j;

	delete [] Velocity_i;		delete [] Velocity_j;
	delete [] Pressure_i;		delete [] Pressure_j;
	delete [] MeanVelocity;		delete [] Proj_flux_tensor;
	delete [] MeanEnergy;		delete [] MeanDensity;	delete [] 	MeanPressure;
	delete [] Energy_i;	delete [] Energy_j;	
	delete [] Energy_vib_i;	delete [] Energy_vib_j;	
	delete [] MeanEnthalpy;		delete [] MeanVelocity;
	delete [] MeanLambda;		delete [] StretchingFactor;
	delete [] Epsilon_0;
	delete [] SoundSpeed_i;		delete [] SoundSpeed_j ;
	delete [] Enthalpy_i;		delete [] Enthalpy_j;
	delete [] Lambda_i;			delete [] Lambda_j;
	delete [] Sensor_i;			delete [] Sensor_j;

}

void CCentLax_AdjPlasmaDiatomic::SetResidual(double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j, 
		double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
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

		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i;
			Velocity_j[iSpecies][iDim] = U_j[loc + iDim+1] / Density_j;
			MeanVelocity[iSpecies][iDim] =  0.5*(Velocity_i[iSpecies][iDim]+Velocity_j[iSpecies][iDim]);
		}

		MeanPressure[iSpecies] = 0.5*(Pressure_i[iSpecies]+Pressure_j[iSpecies]);
		MeanEnthalpy[iSpecies] = 0.5*(Enthalpy_i[iSpecies]+Enthalpy_j[iSpecies]);		
	}

	/*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
	GetInviscidProjJac_(Velocity_i, Energy_i, Energy_vib_i, Enthalpy_i,  Normal, 0.5, Proj_Jac_Tensor_i, config);
	GetInviscidProjJac_(Velocity_j, Energy_j, Energy_vib_j, Enthalpy_j,  Normal, 0.5, Proj_Jac_Tensor_j, config);


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


	/*--- Computes differences btw. conservative variables ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Diff_Psi[iVar] = Psi_i[iVar]-Psi_j[iVar];

	sc0 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));

	/*--- Compute the local espectral radius and the stretching factor ---*/
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

		for (iVar = 0; iVar < nVar_Species; iVar++) {
			val_resvisc_i[loc + iVar] = -Epsilon_0[iSpecies]*Diff_Psi[loc + iVar]*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			val_resvisc_j[loc + iVar] = Epsilon_0[iSpecies]*Diff_Psi[loc + iVar]*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
		}		
	}

	/*--- Contribution to implicit part ---*/
	if (implicit) {	
		for (iVar = 0; iVar < nVar; iVar++) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_0[iSpecies]*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			val_Jacobian_ij[iVar][iVar] += Epsilon_0[iSpecies]*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			val_Jacobian_ji[iVar][iVar] += Epsilon_0[iSpecies]*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
			val_Jacobian_jj[iVar][iVar] -= Epsilon_0[iSpecies]*StretchingFactor[iSpecies]*MeanLambda[iSpecies];
		}
	}
	/*	cout << "Conv Residual i: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << val_resconv_i[iVar] << endl;
	cout << endl << "************" << endl;
	cout << "Conv Residual j: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << val_resconv_j[iVar] << endl;
	cin.get();

	cout << "Visc Residual i: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << val_resvisc_i[iVar] << endl;
	cout << endl << "************" << endl;
	cout << "Visc Residual j: " << endl;
	for (iVar = 0; iVar < nVar; iVar++)
		cout << val_resvisc_j[iVar] << endl;
	cin.get();*/

	/*	cout << endl << "************" << endl;
	cout << "Jacobian_ii: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ii[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_ij: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ij[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_ji: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_ji[iVar][jVar] << "\t";
		}
		cout << endl;
	}

	cout << endl << "************" << endl;	
	cout << "Jacobian_jj: " << endl;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			cout << val_Jacobian_jj[iVar][jVar] << "\t";
		}
		cout << endl;
	}
	cin.get();*/

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
}
