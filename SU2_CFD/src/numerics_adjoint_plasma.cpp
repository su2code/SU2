/*!
 * \file numerics_adjoint_plasma.cpp
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

void CUpwRoe_AdjPlasmaDiatomic::ComputeResidual(double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii,
                                            double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) {
  
  double Gamma, Vel2, hf;
	unsigned short iSpecies,  loc = 0;
	Area = 0;
  
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);                    /*! Area of the face*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/
  
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    Gamma = config->GetSpecies_Gamma(iSpecies);
    hf    = config->GetEnthalpy_Formation(iSpecies);
    
    /*--- Point i: Recompute quantities in case of flux reconstruction ---*/
		Density_i[iSpecies]	   = U_i[loc + 0];
    Energy_i[iSpecies]		 = U_i[loc+nDim+1] / Density_i[iSpecies];
    Energy_el_i[iSpecies]  = 0.0;
    Energy_vib_i[iSpecies] = 0.0;
    if (iSpecies < nDiatomics)
      Energy_vib_i[iSpecies] = U_i[loc+nDim+2] / Density_i[iSpecies];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i[iSpecies];
			Vel2 += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
		}
    Pressure_i[iSpecies] = (Gamma-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 0.5*Vel2 - hf -
                                                                Energy_vib_i[iSpecies] - Energy_el_i[iSpecies]);
    SoundSpeed_i[iSpecies] = sqrt(Gamma*(Gamma - 1.0)*(Energy_i[iSpecies] - 0.5*Vel2  - hf
                                                       - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies]));
		Enthalpy_i[iSpecies] = Energy_i[iSpecies] + Pressure_i[iSpecies] / Density_i[iSpecies];
    
		/*--- Point j: Recompute quantities in case of flux reconstruction ---*/
		Density_j[iSpecies]	   = U_j[loc + 0];
    Energy_j[iSpecies]		 = U_j[loc+nDim+1] / Density_j[iSpecies];
		Energy_el_j[iSpecies]  = 0.0;
		Energy_vib_j[iSpecies] = 0.0;
		if (iSpecies < nDiatomics)
			Energy_vib_j[iSpecies]	= U_j[loc + nDim+2]/Density_j[iSpecies];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_j[iSpecies][iDim] = U_j[loc+iDim+1] / Density_j[iSpecies];
			Vel2 += Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
		}
    Pressure_j[iSpecies] = (Gamma-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 0.5*Vel2 - hf -
                                                                Energy_vib_j[iSpecies] - Energy_el_j[iSpecies]);
    SoundSpeed_j[iSpecies] = sqrt(Gamma*(Gamma - 1.0)*(Energy_j[iSpecies] - 0.5*Vel2  - hf -
                                                       Energy_vib_j[iSpecies] - Energy_el_j[iSpecies]));
    Enthalpy_j[iSpecies] = Energy_j[iSpecies] + Pressure_j[iSpecies] / Density_j[iSpecies];
    
		/*--- Calculate Roe-averaged variables ---*/
		R = sqrt(fabs(Density_j[iSpecies]/Density_i[iSpecies]));
		RoeDensity[iSpecies] = R*Density_i[iSpecies];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			RoeVelocity[iSpecies][iDim] = (R*Velocity_j[iSpecies][iDim]+Velocity_i[iSpecies][iDim])/(R+1.0);
			Vel2 += RoeVelocity[iSpecies][iDim]*RoeVelocity[iSpecies][iDim];
		}
		RoeEnthalpy[iSpecies] = (R*Enthalpy_j[iSpecies]+Enthalpy_i[iSpecies])/(R+1);
		RoeEnergy_vib[iSpecies] = (R*Energy_vib_j[iSpecies] + Energy_vib_i[iSpecies])/(R+1);
    RoeSoundSpeed[iSpecies] = (R*SoundSpeed_j[iSpecies] + SoundSpeed_i[iSpecies]/(R+1));
    
    /*		if (iSpecies < nDiatomics)
     RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaDiatomic-1.0)*(RoeEnthalpy[iSpecies] - 0.5*Vel2 - RoeEnergy_vib[iSpecies] - config->GetEnthalpy_Formation(iSpecies))));
     else
     RoeSoundSpeed[iSpecies] = sqrt(fabs((GammaMonatomic-1.0)*(RoeEnthalpy[iSpecies] - 0.5*Vel2 - config->GetEnthalpy_Formation(iSpecies))));*/
    
	}
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		ProjVelocity[iSpecies]   = 0.0;
		ProjVelocity_i[iSpecies] = 0.0;
		ProjVelocity_j[iSpecies] = 0.0;
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
			val_residual_i[iVar] -= Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar])*Area;
			val_residual_j[iVar] += Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar])*Area;
		}
  
	/*--- Implicit contributions, Transpose the matrices and store the Jacobians. Note the negative
	 sign for the ji and jj Jacobians bc the normal direction is flipped. ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_ii[jVar][iVar] = (Proj_Jac_Tensor_i[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar])*Area;
				val_Jacobian_ij[jVar][iVar] = (Proj_Jac_Tensor_i[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar])*Area;
				val_Jacobian_ji[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar])*Area;
				val_Jacobian_jj[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar])*Area;
			}
		}
	}
}

CUpwSW_AdjPlasmaDiatomic::CUpwSW_AdjPlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics, unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
  
	unsigned short iVar, iSpecies;
  
	nMonatomics = val_nMonatomics;
	nDiatomics  = val_nDiatomics;
  
	GammaMonatomic = config->GetGammaMonatomic();
	GammaDiatomic = config->GetGammaDiatomic();
  
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	Diff_U = new double [nVar];
  
	Density_i    = new double[nSpecies];
	Energy_i     = new double[nSpecies];
	Energy_vib_i = new double [nSpecies];
	Energy_el_i  = new double[nSpecies];
	SoundSpeed_i = new double[nSpecies];
	Pressure_i   = new double[nSpecies];
	Enthalpy_i   = new double[nSpecies];
  
	Density_j    = new double[nSpecies];
	Energy_j     = new double[nSpecies];
	Energy_vib_j = new double[nSpecies];
	Energy_el_j  = new double[nSpecies];
	SoundSpeed_j = new double[nSpecies];
	Pressure_j   = new double[nSpecies];
	Enthalpy_j   = new double[nSpecies];
  
	ProjVelocity_i = new double[nSpecies];
	ProjVelocity_j = new double[nSpecies];
  
	Velocity_i = new double* [nSpecies];
	Velocity_j = new double* [nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Velocity_i[iSpecies] = new double [nDim];
		Velocity_j[iSpecies] = new double [nDim];
	}
  
	Proj_ModJac_Tensor = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Proj_ModJac_Tensor[iVar] = new double [nVar];
  
	Proj_flux_tensor_i = new double [nVar];
	Proj_flux_tensor_j = new double [nVar];
	Lambda_p = new double [nVar];
  Lambda_m = new double [nVar];
  
	P_Tensor			= new double* [nVar];
	invP_Tensor			= new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]  = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
	}
}

CUpwSW_AdjPlasmaDiatomic::~CUpwSW_AdjPlasmaDiatomic(void) {
	unsigned short iVar, iSpecies;
  
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
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		delete [] Velocity_i[iSpecies];
		delete [] Velocity_j[iSpecies];
	}
  
	delete [] Velocity_i;
	delete [] Velocity_j;
  
	delete [] ProjVelocity_i;
	delete [] ProjVelocity_j;
  
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
	delete [] Lambda_p;
  delete [] Lambda_m;
  
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
		delete [] Proj_ModJac_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
	delete [] Proj_ModJac_Tensor;
}

void CUpwSW_AdjPlasmaDiatomic::ComputeResidual(double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii,
                                           double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) {
  
	unsigned short iSpecies,  loc = 0;
  double hf, Gamma;
  
  /*--- Initialize residuals and Jacobians ---*/
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
  
  /*--- Geometry parameters ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);                    /*! Area of the face*/
	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;   /* ! Unit Normal*/
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    Gamma = config->GetSpecies_Gamma(iSpecies);
    hf    = config->GetEnthalpy_Formation(iSpecies);
    
    /*--- Point i: Recompute quantities in case of flux reconstruction ---*/
		Density_i[iSpecies]	   = U_i[loc + 0];
    Energy_i[iSpecies]		 = U_i[loc+nDim+1] / Density_i[iSpecies];
    Energy_el_i[iSpecies]  = 0.0;
    Energy_vib_i[iSpecies] = 0.0;
    if (iSpecies < nDiatomics)
      Energy_vib_i[iSpecies] = U_i[loc+nDim+2] / Density_i[iSpecies];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iSpecies][iDim] = U_i[loc + iDim+1] / Density_i[iSpecies];
			Vel2 += Velocity_i[iSpecies][iDim]*Velocity_i[iSpecies][iDim];
		}
    Pressure_i[iSpecies] = (Gamma-1.0) * Density_i[iSpecies] * (Energy_i[iSpecies] - 0.5*Vel2 - hf -
                                                                Energy_vib_i[iSpecies] - Energy_el_i[iSpecies]);
    SoundSpeed_i[iSpecies] = sqrt(Gamma*(Gamma - 1.0)*(Energy_i[iSpecies] - 0.5*Vel2  - hf
                                                       - Energy_vib_i[iSpecies] - Energy_el_i[iSpecies]));
		Enthalpy_i[iSpecies] = Energy_i[iSpecies] + Pressure_i[iSpecies] / Density_i[iSpecies];
    
		/*--- Point j: Recompute quantities in case of flux reconstruction ---*/
		Density_j[iSpecies]	   = U_j[loc + 0];
    Energy_j[iSpecies]		 = U_j[loc+nDim+1] / Density_j[iSpecies];
		Energy_el_j[iSpecies]  = 0.0;
		Energy_vib_j[iSpecies] = 0.0;
		if (iSpecies < nDiatomics)
			Energy_vib_j[iSpecies]	= U_j[loc + nDim+2]/Density_j[iSpecies];
		Vel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_j[iSpecies][iDim] = U_j[loc+iDim+1] / Density_j[iSpecies];
			Vel2 += Velocity_j[iSpecies][iDim]*Velocity_j[iSpecies][iDim];
		}
    Pressure_j[iSpecies] = (Gamma-1.0) * Density_j[iSpecies] * (Energy_j[iSpecies] - 0.5*Vel2 - hf -
                                                                Energy_vib_j[iSpecies] - Energy_el_j[iSpecies]);
    SoundSpeed_j[iSpecies] = sqrt(Gamma*(Gamma - 1.0)*(Energy_j[iSpecies] - 0.5*Vel2  - hf -
                                                       Energy_vib_j[iSpecies] - Energy_el_j[iSpecies]));
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
	/*--- Flow eigenvalues for Fij--- */
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
    /*--- Lambda+ at i ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda_p[loc+iDim] = 0.5*(ProjVelocity_i[iSpecies] + fabs(ProjVelocity_i[iSpecies]));
		}
		Lambda_p[loc+nDim]   = 0.5*(ProjVelocity_i[iSpecies] + SoundSpeed_i[iSpecies]
                                + fabs(ProjVelocity_i[iSpecies] + SoundSpeed_i[iSpecies]));
		Lambda_p[loc+nDim+1] = 0.5*(ProjVelocity_i[iSpecies] - SoundSpeed_i[iSpecies]
                                + fabs(ProjVelocity_i[iSpecies] - SoundSpeed_i[iSpecies]));
		if (iSpecies < nDiatomics) Lambda_p[loc+nDim+2] = 0.5*(ProjVelocity_i[iSpecies] + fabs(ProjVelocity_i[iSpecies]));
    
    /*--- Lambda- at j ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
			Lambda_m[loc+iDim] = 0.5*(ProjVelocity_j[iSpecies] - fabs(ProjVelocity_j[iSpecies]));
		}
		Lambda_m[loc+nDim]   = 0.5*(ProjVelocity_j[iSpecies] + SoundSpeed_j[iSpecies]
                                - fabs(ProjVelocity_j[iSpecies] + SoundSpeed_j[iSpecies]));
		Lambda_m[loc+nDim+1] = 0.5*(ProjVelocity_j[iSpecies] - SoundSpeed_i[iSpecies]
                                - fabs(ProjVelocity_j[iSpecies] - SoundSpeed_j[iSpecies]));
		if (iSpecies < nDiatomics) Lambda_m[loc+nDim+2] = 0.5*(ProjVelocity_j[iSpecies] - fabs(ProjVelocity_j[iSpecies]));
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
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda_p[kVar]*invP_Tensor[kVar][jVar];
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
  
	/*--- Compute P & invP at j ---*/
	GetPMatrix_(Density_j, Velocity_j, Enthalpy_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_j, Velocity_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, invP_Tensor);
  
	/*--- Projected Jacobian (A-) at j ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda_m[kVar]*invP_Tensor[kVar][jVar];
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
  
	/*--- Flow eigenvalues for Fji --- */
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
    /*--- Lambda- at i (with -nij) ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda_m[loc+iDim] = 0.5*(-ProjVelocity_i[iSpecies] - fabs(-ProjVelocity_i[iSpecies]));
		}
		Lambda_m[loc+nDim]   = 0.5*(-ProjVelocity_i[iSpecies] + SoundSpeed_i[iSpecies]
                                - fabs(-ProjVelocity_i[iSpecies] + SoundSpeed_i[iSpecies]));
		Lambda_m[loc+nDim+1] = 0.5*(-ProjVelocity_i[iSpecies] - SoundSpeed_i[iSpecies]
                                - fabs(-ProjVelocity_i[iSpecies] - SoundSpeed_i[iSpecies]));
		if (iSpecies < nDiatomics) Lambda_m[loc+nDim+2] = 0.5*(-ProjVelocity_i[iSpecies] - fabs(-ProjVelocity_i[iSpecies]));
    
    /*--- Lambda+ at j (with -nij) ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Lambda_p[loc+iDim] = 0.5*(-ProjVelocity_j[iSpecies] + fabs(-ProjVelocity_j[iSpecies]));
		}
		Lambda_p[loc+nDim]   = 0.5*(-ProjVelocity_j[iSpecies] + SoundSpeed_j[iSpecies]
                                - fabs(-ProjVelocity_j[iSpecies] + SoundSpeed_j[iSpecies]));
		Lambda_p[loc+nDim+1] = 0.5*(-ProjVelocity_j[iSpecies] - SoundSpeed_j[iSpecies]
                                - fabs(-ProjVelocity_j[iSpecies] - SoundSpeed_j[iSpecies]));
		if (iSpecies < nDiatomics) Lambda_p[loc+nDim+2] = 0.5*(-ProjVelocity_j[iSpecies] - fabs(-ProjVelocity_j[iSpecies]));
	}
  
	/*--- Compute P & invP at i ---*/
	GetPMatrix_(Density_i, Velocity_i, Enthalpy_i, SoundSpeed_i, Energy_vib_i, Energy_el_i, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_i, Velocity_i, SoundSpeed_i, Energy_vib_i, Energy_el_i, config, UnitaryNormal, invP_Tensor);
  
	/*--- Projected Jacobian (A-) at i ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda_m[kVar]*invP_Tensor[kVar][jVar];
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
  
	/*--- Compute P & invP at j ---*/
	GetPMatrix_(Density_j, Velocity_j, Enthalpy_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, P_Tensor);
	GetPMatrix_inv_(Density_j, Velocity_j, SoundSpeed_j, Energy_vib_j, Energy_el_j, config, UnitaryNormal, invP_Tensor);
  
	/*--- Projected Jacobian (A-) at j ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			Proj_ModJac_Tensor_ij = 0.0;
			/*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
			for (kVar = 0; kVar < nVar; kVar++)
				Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda_p[kVar]*invP_Tensor[kVar][jVar];
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

void CCentLax_AdjPlasmaDiatomic::ComputeResidual(double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j,
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