/*!
 * \file numerics_viscous.cpp
 * \brief This file contains all the viscous term discretization.
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

CAvgGrad_Flow::CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Allocate some useful vectors ---*/
	Velocity_i = new double [nVar];
	Velocity_j = new double [nVar];
	Prim_Var_i = new double [nVar];
	Prim_Var_j = new double [nVar];
	Mean_PrimVar = new double [nVar];
	Mean_GradPrimVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGrad_Flow::~CAvgGrad_Flow(void) {

	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Prim_Var_i;
	delete [] Prim_Var_j;
	delete [] Mean_PrimVar;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGrad_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	/*--- Normalized normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);

	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;

	/*--- Conversion to Primitive Variables ---*/
	ConsVar2PrimVar(U_i, Prim_Var_i);
	ConsVar2PrimVar(U_j, Prim_Var_j);

	/*--- Mean Primitive variable value ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_PrimVar[iVar] = 0.5*(Prim_Var_i[iVar]+Prim_Var_j[iVar]);

	/*--- Mean Viscosities ---*/
	Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
	Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);

	/*--- Mean gradient approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);

	/*--- Get projected flux tensor ---*/
	GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

	/*--- Update viscous residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Proj_Flux_Tensor[iVar];

	/*--- Implicit part ---*/
	if (implicit) {
		dist_ij = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
		dist_ij = sqrt(dist_ij);
		Mean_Density = 0.5*(U_i[0]+U_j[0]);
		GetViscousProjJacs(Mean_Density, Mean_PrimVar[nVar-1], Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, 
				dist_ij, UnitaryNormal, Area, Mean_PrimVar, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
	}
}

CAvgGradArtComp_Flow::CAvgGradArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

	/*--- Allocate some useful vectors ---*/
	Velocity_i = new double [nVar];
	Velocity_j = new double [nVar];
	Prim_Var_i = new double [nVar];
	Prim_Var_j = new double [nVar];
	Mean_PrimVar = new double [nVar];
	Mean_GradPrimVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
	
}

CAvgGradArtComp_Flow::~CAvgGradArtComp_Flow(void) {

	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Prim_Var_i;
	delete [] Prim_Var_j;
	delete [] Mean_PrimVar;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGradArtComp_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, 
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
		Prim_Var_i[iVar] = U_i[iVar];
		Prim_Var_j[iVar] = U_j[iVar];		
	}

	/*--- Mean Primitive variable value ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_PrimVar[iVar] = 0.5*(Prim_Var_i[iVar]+Prim_Var_j[iVar]);

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
		
		GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, 
															dist_ij, UnitaryNormal, Area, val_Jacobian_i, val_Jacobian_j);
	}
}

CAvgGrad_TurbSA::CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	sigma = 2./3.;

	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Kappa = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGrad_TurbSA::~CAvgGrad_TurbSA(void) {
	unsigned short iVar;

	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Kappa;
	delete [] Proj_Mean_GradTurbVar_Edge;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGrad_TurbSA::SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {

	switch (config->GetKind_Turb_Model()) {
	case SA :

		/*--- Compute mean effective viscosity ---*/
		nu_i = Laminar_Viscosity_i/U_i[0];
		nu_j = Laminar_Viscosity_j/U_j[0];
		nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);

		/*--- Compute vector going from iPoint to jPoint ---*/
		dist_ij_2 = 0; proj_vector_ij = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
			dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
			proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
		}
		proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors

		/*--- Mean gradient approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
			Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
				Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
			}
		}

		val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;

		/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
		if (implicit) {
			Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
			Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
		}
		break;

	case SA_COMP :

		/*--- Compute mean effective viscosity ---*/
		nu_e = 0.5*(Laminar_Viscosity_i+Laminar_Viscosity_j+TurbVar_i[0]+TurbVar_j[0]);
		nu_hat_i = TurbVar_i[0]/U_i[0];
		nu_hat_j = TurbVar_j[0]/U_j[0];

		/*--- Compute vector going from iPoint to jPoint ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
			dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
			proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
		}
		proj_vector_ij = proj_vector_ij/dist_ij_2;

		/*--- Mean gradient approximation. Projection of the mean
			 gradient in the direction of the edge ---*/
		Proj_Mean_GradTurbVar_Kappa[0] = 0.0;
		Proj_Mean_GradTurbVar_Edge[0] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbVar[0][iDim] =
					0.5*( (TurbVar_Grad_i[0][iDim] - ConsVar_Grad_i[0][iDim]*nu_hat_i)/U_i[0]
					                                                                       +(TurbVar_Grad_j[0][iDim] - ConsVar_Grad_j[0][iDim]*nu_hat_j)/U_j[0] );
			Proj_Mean_GradTurbVar_Kappa[0] += Mean_GradTurbVar[0][iDim]*Normal[iDim];
		}

		val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;

		/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
		if (implicit) {
			Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e/U_i[0]*proj_vector_ij)/sigma;
			Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e/U_j[0]*proj_vector_ij)/sigma;
		}
		break;
	}
}

CAvgGrad_TurbSST::CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	sigma = 2./3.;
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Kappa = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGrad_TurbSST::~CAvgGrad_TurbSST(void) {
	unsigned short iVar;

	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Kappa;
	delete [] Proj_Mean_GradTurbVar_Edge;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGrad_TurbSST::SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {

	/*--- Compute mean effective viscosity ---*/
	diff_i = Laminar_Viscosity_i + Eddy_Viscosity_i;
	diff_j = Laminar_Viscosity_j + Eddy_Viscosity_j;
	diff_e = 0.5*(diff_i + diff_j);

	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0; proj_vector_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
		proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
	}
	proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors

	/*--- Mean gradient approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
		Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
		}
	}

	val_residual[0] = diff_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;

	/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
	if (implicit) {
		Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-diff_e*proj_vector_ij)/sigma;
		Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+diff_e*proj_vector_ij)/sigma;
	}
}

CAvgGrad_AdjFlow::CAvgGrad_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iDim;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	Mean_Velocity = new double [nDim];
	Mean_GradPhi = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Mean_GradPhi[iDim] = new double [nDim];
	Mean_GradPsiE = new double [nDim];

}

CAvgGrad_AdjFlow::~CAvgGrad_AdjFlow(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Mean_Velocity;
	delete [] Mean_GradPsiE;
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		delete [] Mean_GradPhi[iDim];
}

void CAvgGrad_AdjFlow::SetResidual (double *val_residual_i, double *val_residual_j) {
	unsigned short iDim, jDim;
	double Density_i, sq_vel_i, Energy_i, SoundSpeed_i, Pressure_i, ViscDens_i, XiDens_i, 
	Density_j, sq_vel_j, Energy_j, SoundSpeed_j, Pressure_j, ViscDens_j, XiDens_j, dPhiE_dn;
	double Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
	Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5, 
	Sigma_yz5, Sigma_5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;	

	/*--- States in the point i ---*/
	Density_i = U_i[0];		
	sq_vel_i = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
	}
	Energy_i = U_i[nDim+1] / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-sq_vel_i));
	Pressure_i = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	ViscDens_i = (Laminar_Viscosity_i + Eddy_Viscosity_i) / Density_i;
	XiDens_i = Gamma *(Laminar_Viscosity_i/PRANDTL + Eddy_Viscosity_i/PRANDTL_TURB) / Density_i;

	/*--- States in the point j ---*/
	Density_j = U_j[0];
	sq_vel_j = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_j = U_j[nDim+1] / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-sq_vel_j));
	Pressure_j = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	ViscDens_j = (Laminar_Viscosity_j + Eddy_Viscosity_j) / Density_j;
	XiDens_j = Gamma *(Laminar_Viscosity_j/PRANDTL + Eddy_Viscosity_j/PRANDTL_TURB) / Density_j;

	/*--- Average of the derivatives of the adjoint variables ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Mean_GradPsiE[iDim] =  0.5*(PsiVar_Grad_i[nVar-1][iDim]+PsiVar_Grad_j[nVar-1][iDim]);
		for (jDim = 0; jDim < nDim; jDim++)
			Mean_GradPhi[iDim][jDim] =  0.5*(PsiVar_Grad_i[iDim+1][jDim]+PsiVar_Grad_j[iDim+1][jDim]);
	}

	dPhiE_dn = 0; 
	for (iDim = 0; iDim < nDim; iDim++) 
		dPhiE_dn += Mean_GradPsiE[iDim]*Normal[iDim];

	/*--- Compute the viscous residual ---*/
	if (nDim == 3) {
		Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_zz = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
		Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xz = ViscDens_i * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
		Sigma_yz = ViscDens_i * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
		Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_zz5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
		Sigma_xz5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[0]);
		Sigma_yz5 = ViscDens_i * (Velocity_i[1] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[1]);
		Sigma_5   = XiDens_i * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
		eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

		val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy + Velocity_i[2] * Normal[2] * eta_zz
				+ (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
				+ (Velocity_i[0] * Normal[2] + Velocity_i[2] * Normal[0]) * eta_xz
				+ (Velocity_i[2] * Normal[1] + Velocity_i[1] * Normal[2]) * eta_yz
				- (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
		val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_i[0] * Sigma_5);
		val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_i[1] * Sigma_5);
		val_residual_i[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_i[2] * Sigma_5);
		val_residual_i[4] = (Sigma_5);

		Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_zz = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
		Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xz = ViscDens_j * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
		Sigma_yz = ViscDens_j * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
		Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_zz5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
		Sigma_xz5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[0]);
		Sigma_yz5 = ViscDens_j * (Velocity_j[1] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[1]);
		Sigma_5   = XiDens_j * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
		eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

		val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy + Velocity_j[2] * Normal[2] * eta_zz
				+ (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
				+ (Velocity_j[0] * Normal[2] + Velocity_j[2] * Normal[0]) * eta_xz
				+ (Velocity_j[2] * Normal[1] + Velocity_j[1] * Normal[2]) * eta_yz
				- (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
		val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_j[0] * Sigma_5);
		val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_j[1] * Sigma_5);
		val_residual_j[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_j[2] * Sigma_5);
		val_residual_j[4] = (Sigma_5);
	}

	if (nDim == 2) {
		Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
		Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
		Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1]);
		Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1]);
		Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
		Sigma_5   = XiDens_i * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;

		val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy
				+ (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
				- (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
		val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] - Velocity_i[0] * Sigma_5);
		val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] - Velocity_i[1] * Sigma_5);
		val_residual_i[3] = (Sigma_5);

		Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
		Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
		Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1]);
		Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1]);
		Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
		Sigma_5   = XiDens_j * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;

		val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy
				+ (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
				- (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
		val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1]  - Velocity_j[0] * Sigma_5);
		val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1]  - Velocity_j[1] * Sigma_5);
		val_residual_j[3] = (Sigma_5);
	}
}

CAvgGradArtComp_AdjFlow::CAvgGradArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iDim;
	
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	Mean_Velocity = new double [nDim];
	Mean_GradPhi = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Mean_GradPhi[iDim] = new double [nDim];
	Mean_GradPsiE = new double [nDim];
	
}

CAvgGradArtComp_AdjFlow::~CAvgGradArtComp_AdjFlow(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Mean_Velocity;
	delete [] Mean_GradPsiE;
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		delete [] Mean_GradPhi[iDim];
}

void CAvgGradArtComp_AdjFlow::SetResidual (double *val_residual_i, double *val_residual_j) {
	unsigned short iDim, jDim;
	double Density_i, sq_vel_i, Energy_i, SoundSpeed_i, Pressure_i, ViscDens_i, XiDens_i, 
	Density_j, sq_vel_j, Energy_j, SoundSpeed_j, Pressure_j, ViscDens_j, XiDens_j, dPhiE_dn;
	double Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
	Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5, 
	Sigma_yz5, Sigma_5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;	
	
	/*--- States in the point i ---*/
	Density_i = U_i[0];		
	sq_vel_i = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
	}
	Energy_i = U_i[nDim+1] / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-sq_vel_i));
	Pressure_i = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	ViscDens_i = (Laminar_Viscosity_i + Eddy_Viscosity_i) / Density_i;
	XiDens_i = Gamma *(Laminar_Viscosity_i/PRANDTL + Eddy_Viscosity_i/PRANDTL_TURB) / Density_i;
	
	/*--- States in the point j ---*/
	Density_j = U_j[0];
	sq_vel_j = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_j = U_j[nDim+1] / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-sq_vel_j));
	Pressure_j = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	ViscDens_j = (Laminar_Viscosity_j + Eddy_Viscosity_j) / Density_j;
	XiDens_j = Gamma *(Laminar_Viscosity_j/PRANDTL + Eddy_Viscosity_j/PRANDTL_TURB) / Density_j;
	
	/*--- Average of the derivatives of the adjoint variables ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Mean_GradPsiE[iDim] =  0.5*(PsiVar_Grad_i[nVar-1][iDim]+PsiVar_Grad_j[nVar-1][iDim]);
		for (jDim = 0; jDim < nDim; jDim++)
			Mean_GradPhi[iDim][jDim] =  0.5*(PsiVar_Grad_i[iDim+1][jDim]+PsiVar_Grad_j[iDim+1][jDim]);
	}
	
	dPhiE_dn = 0; 
	for (iDim = 0; iDim < nDim; iDim++) 
		dPhiE_dn += Mean_GradPsiE[iDim]*Normal[iDim];
	
	/*--- Compute the viscous residual ---*/
	if (nDim == 3) {
		Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_zz = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
		Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xz = ViscDens_i * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
		Sigma_yz = ViscDens_i * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
		Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_zz5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
		Sigma_xz5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[0]);
		Sigma_yz5 = ViscDens_i * (Velocity_i[1] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[1]);
		Sigma_5   = XiDens_i * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
		eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;
		
		val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy + Velocity_i[2] * Normal[2] * eta_zz
													 + (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
													 + (Velocity_i[0] * Normal[2] + Velocity_i[2] * Normal[0]) * eta_xz
													 + (Velocity_i[2] * Normal[1] + Velocity_i[1] * Normal[2]) * eta_yz
													 - (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
		val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_i[0] * Sigma_5);
		val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_i[1] * Sigma_5);
		val_residual_i[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_i[2] * Sigma_5);
		val_residual_i[4] = (Sigma_5);
		
		Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_zz = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
		Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xz = ViscDens_j * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
		Sigma_yz = ViscDens_j * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
		Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_zz5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
		Sigma_xz5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[0]);
		Sigma_yz5 = ViscDens_j * (Velocity_j[1] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[1]);
		Sigma_5   = XiDens_j * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
		eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;
		
		val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy + Velocity_j[2] * Normal[2] * eta_zz
													 + (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
													 + (Velocity_j[0] * Normal[2] + Velocity_j[2] * Normal[0]) * eta_xz
													 + (Velocity_j[2] * Normal[1] + Velocity_j[1] * Normal[2]) * eta_yz
													 - (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
		val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_j[0] * Sigma_5);
		val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_j[1] * Sigma_5);
		val_residual_j[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_j[2] * Sigma_5);
		val_residual_j[4] = (Sigma_5);
	}
	
	if (nDim == 2) {
		Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
		Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
		Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1]);
		Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1]);
		Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
		Sigma_5   = XiDens_i * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;
		
		val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy
													 + (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
													 - (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
		val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] - Velocity_i[0] * Sigma_5);
		val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] - Velocity_i[1] * Sigma_5);
		val_residual_i[3] = (Sigma_5);
		
		Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
		Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
		Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1]);
		Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1]);
		Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
		Sigma_5   = XiDens_j * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;
		
		val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy
													 + (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
													 - (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
		val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1]  - Velocity_j[0] * Sigma_5);
		val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1]  - Velocity_j[1] * Sigma_5);
		val_residual_j[3] = (Sigma_5);
	}
}

CAvgGradCorrected_Flow::CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

	Velocity_i = new double [nVar];
	Velocity_j = new double [nVar];
	Prim_Var_i = new double [nVar];
	Prim_Var_j = new double [nVar];
	Mean_PrimVar = new double [nVar];
	Proj_Mean_GradPrimVar_Edge = new double [nVar];
	Edge_Vector = new double [nDim];

	Mean_GradPrimVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGradCorrected_Flow::~CAvgGradCorrected_Flow(void) {

	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Prim_Var_i;
	delete [] Prim_Var_j;
	delete [] Mean_PrimVar;
	delete [] Proj_Mean_GradPrimVar_Edge;
	delete [] Edge_Vector;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGradCorrected_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	/*--- Normalized normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);

	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;

	/*--- Conversion to Primitive Variables ---*/
	ConsVar2PrimVar(U_i, Prim_Var_i);
	ConsVar2PrimVar(U_j, Prim_Var_j);

	/*--- Mean Primitive variable value ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_PrimVar[iVar] = 0.5*(Prim_Var_i[iVar]+Prim_Var_j[iVar]);

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
		for (iDim = 0; iDim < nDim; iDim++)
			Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] - 
					(Prim_Var_j[iVar]-Prim_Var_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
	}

	/*--- Get projected flux tensor ---*/
	GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

	/*--- Save residual value ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Proj_Flux_Tensor[iVar];

	/*--- Compute the implicit part ---*/
	if (implicit) {
		GetViscousProjJacs(Mean_PrimVar[0], Mean_PrimVar[nVar-1], Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, 
				sqrt(dist_ij_2), UnitaryNormal, Area, Mean_PrimVar, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
	}
}

CAvgGradCorrectedArtComp_Flow::CAvgGradCorrectedArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;

	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

	Velocity_i = new double [nVar];
	Velocity_j = new double [nVar];
	Prim_Var_i = new double [nVar];
	Prim_Var_j = new double [nVar];
	Mean_PrimVar = new double [nVar];
	Proj_Mean_GradPrimVar_Edge = new double [nVar];
	Edge_Vector = new double [nDim];

	Mean_GradPrimVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGradCorrectedArtComp_Flow::~CAvgGradCorrectedArtComp_Flow(void) {

	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Prim_Var_i;
	delete [] Prim_Var_j;
	delete [] Mean_PrimVar;
	delete [] Proj_Mean_GradPrimVar_Edge;
	delete [] Edge_Vector;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGradCorrectedArtComp_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	/*--- Normalized normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);

	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;

	/*--- Conversion to Primitive Variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Prim_Var_i[iVar] = U_i[iVar];
		Prim_Var_j[iVar] = U_j[iVar];		
	}

	/*--- Mean Primitive variable value ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_PrimVar[iVar] = 0.5*(Prim_Var_i[iVar]+Prim_Var_j[iVar]);

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
		for (iDim = 0; iDim < nDim; iDim++)
			Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] - 
					(Prim_Var_j[iVar]-Prim_Var_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
	}

	/*--- Get projected flux tensor ---*/
	GetViscousArtCompProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

	/*--- Update viscous residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Proj_Flux_Tensor[iVar];

	/*--- Implicit part ---*/
	if (implicit) {
		double dist_ij = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
		dist_ij = sqrt(dist_ij);
		
		GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, 
															dist_ij, UnitaryNormal, Area, val_Jacobian_i, val_Jacobian_j);
	}
}

CAvgGradCorrected_TurbSA::CAvgGradCorrected_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	sigma = 2./3.;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Kappa = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Proj_Mean_GradTurbVar_Corrected = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGradCorrected_TurbSA::~CAvgGradCorrected_TurbSA(void) {
	unsigned short iVar;

	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Kappa;
	delete [] Proj_Mean_GradTurbVar_Edge;
	delete [] Proj_Mean_GradTurbVar_Corrected;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGradCorrected_TurbSA::SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {

	switch (config->GetKind_Turb_Model()) {
	case SA :
		/*--- Compute mean effective viscosity ---*/
		nu_i = Laminar_Viscosity_i/U_i[0];
		nu_j = Laminar_Viscosity_j/U_j[0];
		nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);

		/*--- Compute vector going from iPoint to jPoint ---*/
		dist_ij_2 = 0; proj_vector_ij = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
			dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
			proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
		}
		proj_vector_ij = proj_vector_ij/dist_ij_2;

		/*--- Mean gradient approximation. Projection of the mean gradient
			 in the direction of the edge ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
			Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
				Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
				Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
			}
			Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
			Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
					(TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
		}

		val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;

		/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
		if (implicit) {
			Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
			Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
		}
		break;

	case SA_COMP :
		/*--- Compute mean effective viscosity ---*/
		nu_e = 0.5*(Laminar_Viscosity_i+Laminar_Viscosity_j+TurbVar_i[0]+TurbVar_j[0]);
		nu_hat_i = TurbVar_i[0]/U_i[0];
		nu_hat_j = TurbVar_j[0]/U_j[0];

		/*--- Compute vector going from iPoint to jPoint ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
			dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
			proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
		}
		proj_vector_ij = proj_vector_ij/dist_ij_2;

		/*--- Mean gradient approximation. Projection of the mean gradient
			 in the direction of the edge ---*/
		Proj_Mean_GradTurbVar_Kappa[0] = 0.0;
		Proj_Mean_GradTurbVar_Edge[0] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbVar[0][iDim] =
					0.5*( (TurbVar_Grad_i[0][iDim] - ConsVar_Grad_i[0][iDim]*nu_hat_i)/U_i[0]
					                                                                       +(TurbVar_Grad_j[0][iDim] - ConsVar_Grad_j[0][iDim]*nu_hat_j)/U_j[0] );
			Proj_Mean_GradTurbVar_Kappa[0] += Mean_GradTurbVar[0][iDim]*Normal[iDim];
			Proj_Mean_GradTurbVar_Edge[0] += Mean_GradTurbVar[0][iDim]*Edge_Vector[iDim];
		}
		Proj_Mean_GradTurbVar_Corrected[0] = Proj_Mean_GradTurbVar_Kappa[0];
		Proj_Mean_GradTurbVar_Corrected[0] -= Proj_Mean_GradTurbVar_Edge[0]*proj_vector_ij -
				(nu_hat_j-nu_hat_i)*proj_vector_ij;

		val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;

		/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
		if (implicit) {
			Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e/U_i[0]*proj_vector_ij)/sigma;
			Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e/U_j[0]*proj_vector_ij)/sigma;
		}
		break;
	}
}

CAvgGradCorrected_TurbSST::CAvgGradCorrected_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	sigma_k1  = 0.85;
	sigma_om1 = 0.5;
	sigma_k2  = 1.0;
	sigma_om2 = 0.856;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Normal = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Proj_Mean_GradTurbVar_Corrected = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGradCorrected_TurbSST::~CAvgGradCorrected_TurbSST(void) {
	unsigned short iVar;

	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Normal;
	delete [] Proj_Mean_GradTurbVar_Edge;
	delete [] Proj_Mean_GradTurbVar_Corrected;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGradCorrected_TurbSST::SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {

	double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
	double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;

	/*--- Compute the blended constant for the viscous terms ---*/
	sigma_kine_i  = F1_i*sigma_k1+(1-F1_i)*sigma_k2;
	sigma_kine_j  = F1_j*sigma_k1+(1-F1_j)*sigma_k2;
	sigma_omega_i = F1_i*sigma_om1+(1-F1_i)*sigma_om2;
	sigma_omega_j = F1_j*sigma_om1+(1-F1_j)*sigma_om2;

	/*--- Compute mean effective viscosity ---*/
	diff_i_kine  = Laminar_Viscosity_i+sigma_kine_i*Eddy_Viscosity_i;
	diff_j_kine  = Laminar_Viscosity_j+sigma_kine_j*Eddy_Viscosity_j;
	diff_i_omega = Laminar_Viscosity_i+sigma_omega_i*Eddy_Viscosity_i;
	diff_j_omega = Laminar_Viscosity_j+sigma_omega_j*Eddy_Viscosity_j;

	diff_kine  = 0.5*(diff_i_kine+diff_j_kine);    // Could instead use weighted average!
	diff_omega = 0.5*(diff_i_omega+diff_j_omega);

	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0; proj_vector_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
		proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
	}
	proj_vector_ij = proj_vector_ij/dist_ij_2;

	/*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
		Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
			Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
		}
		Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
		Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
				(TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
	}

	val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
	val_residual[1] = diff_omega*Proj_Mean_GradTurbVar_Corrected[1];

	/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
	if (implicit) {
		Jacobian_i[0][0] = -diff_kine*proj_vector_ij/U_i[0];		Jacobian_i[0][1] = 0.0;
		Jacobian_i[1][0] = 0.0;									    Jacobian_i[1][1] = -diff_omega*proj_vector_ij/U_i[0];

		Jacobian_j[0][0] = diff_kine*proj_vector_ij/U_j[0]; 		Jacobian_j[0][1] = 0.0;
		Jacobian_j[1][0] = 0.0;									    Jacobian_j[1][1] = diff_omega*proj_vector_ij/U_j[0];
	}

}

CAvgGradCorrected_AdjFlow::CAvgGradCorrected_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	Mean_Velocity = new double [nDim];

	Mean_GradPsiVar = new double* [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Mean_GradPsiVar[iVar] = new double [nDim];

	Edge_Vector = new double [nDim];
	Proj_Mean_GradPsiVar_Edge = new double [nVar];

	Mean_GradPhi = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Mean_GradPhi[iDim] = new double [nDim];
	Mean_GradPsiE = new double [nDim];

}

CAvgGradCorrected_AdjFlow::~CAvgGradCorrected_AdjFlow(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Mean_Velocity;
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradPsiVar_Edge;

	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradPsiVar[iVar];
	delete [] Mean_GradPsiVar;

	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		delete [] Mean_GradPhi[iDim];
	delete [] Mean_GradPhi;
	delete [] Mean_GradPsiE;
}

void CAvgGradCorrected_AdjFlow::SetResidual(double *val_residual_i, double *val_residual_j) {

	unsigned short iVar, iDim, jDim;
	double Density_i, sq_vel_i, Energy_i, SoundSpeed_i, Pressure_i, ViscDens_i, XiDens_i, 
	Density_j, sq_vel_j, Energy_j, SoundSpeed_j, Pressure_j, ViscDens_j, XiDens_j, dist_ij_2, dPhiE_dn;
	double Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
	Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5, 
	Sigma_yz5, Sigma_5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;	

	/*--- States in the point i ---*/
	Density_i = U_i[0];		
	sq_vel_i = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
	}
	Energy_i = U_i[nDim+1] / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-sq_vel_i));
	Pressure_i = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	ViscDens_i = (Laminar_Viscosity_i + Eddy_Viscosity_i) / Density_i;
	XiDens_i = Gamma *(Laminar_Viscosity_i/PRANDTL + Eddy_Viscosity_i/PRANDTL_TURB) / Density_i;

	/*--- States in the point j ---*/
	Density_j = U_j[0];
	sq_vel_j = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_j = U_j[nDim+1] / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-sq_vel_j));
	Pressure_j = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	ViscDens_j = (Laminar_Viscosity_j + Eddy_Viscosity_j) / Density_j;
	XiDens_j = Gamma *(Laminar_Viscosity_j/PRANDTL + Eddy_Viscosity_j/PRANDTL_TURB) / Density_j;

	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Mean_Velocity[iDim] = 0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
	}

	/*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge, weiss correction ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GradPsiVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradPsiVar[iVar][iDim] = 0.5*(PsiVar_Grad_i[iVar][iDim] + PsiVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradPsiVar_Edge[iVar] += Mean_GradPsiVar[iVar][iDim]*Edge_Vector[iDim];
		}
		for (iDim = 0; iDim < nDim; iDim++)
			Mean_GradPsiVar[iVar][iDim] -= (Proj_Mean_GradPsiVar_Edge[iVar] - 
					(Psi_j[iVar]-Psi_i[iVar]))*Edge_Vector[iDim]/dist_ij_2;
	}

	/*--- Average of the derivatives of the adjoint variables ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Mean_GradPsiE[iDim] = Mean_GradPsiVar[nVar-1][iDim];
		for (jDim = 0; jDim < nDim; jDim++)
			Mean_GradPhi[iDim][jDim] = Mean_GradPsiVar[iDim+1][jDim];
	}

	dPhiE_dn = 0; 
	for (iDim = 0; iDim < nDim; iDim++) 
		dPhiE_dn += Mean_GradPsiE[iDim]*Normal[iDim];

	/*--- Compute the viscous residual ---*/
	if (nDim == 3) {
		Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_zz = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
		Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xz = ViscDens_i * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
		Sigma_yz = ViscDens_i * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
		Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_zz5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
		Sigma_xz5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[0]);
		Sigma_yz5 = ViscDens_i * (Velocity_i[1] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[1]);
		Sigma_5   = XiDens_i * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
		eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

		val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy + Velocity_i[2] * Normal[2] * eta_zz
				+ (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
				+ (Velocity_i[0] * Normal[2] + Velocity_i[2] * Normal[0]) * eta_xz
				+ (Velocity_i[2] * Normal[1] + Velocity_i[1] * Normal[2]) * eta_yz
				- (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
		val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_i[0] * Sigma_5);
		val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_i[1] * Sigma_5);
		val_residual_i[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_i[2] * Sigma_5);
		val_residual_i[4] = (Sigma_5);

		Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_zz = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
		Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xz = ViscDens_j * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
		Sigma_yz = ViscDens_j * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
		Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_zz5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
		Sigma_xz5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[0]);
		Sigma_yz5 = ViscDens_j * (Velocity_j[1] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[1]);
		Sigma_5   = XiDens_j * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
		eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

		val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy + Velocity_j[2] * Normal[2] * eta_zz
				+ (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
				+ (Velocity_j[0] * Normal[2] + Velocity_j[2] * Normal[0]) * eta_xz
				+ (Velocity_j[2] * Normal[1] + Velocity_j[1] * Normal[2]) * eta_yz
				- (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
		val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_j[0] * Sigma_5);
		val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_j[1] * Sigma_5);
		val_residual_j[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_j[2] * Sigma_5);
		val_residual_j[4] = (Sigma_5);
	}

	if (nDim == 2) {
		Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
		Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
		Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1]);
		Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1]);
		Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
		Sigma_5   = XiDens_i * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;

		val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy
				+ (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
				- (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
		val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] - Velocity_i[0] * Sigma_5);
		val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] - Velocity_i[1] * Sigma_5);
		val_residual_i[3] = (Sigma_5);

		Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
		Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
		Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1]);
		Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1]);
		Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
		Sigma_5   = XiDens_j * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;

		val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy
				+ (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
				- (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
		val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1]  - Velocity_j[0] * Sigma_5);
		val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1]  - Velocity_j[1] * Sigma_5);
		val_residual_j[3] = (Sigma_5);
	}
}

void CAvgGradCorrected_AdjFlow::SetResidual(double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, double **val_Jacobian_ij,
		double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) {

	unsigned short iVar, jVar, iDim, jDim;
	double Density_i, sq_vel_i, Energy_i, SoundSpeed_i, Pressure_i, ViscDens_i, XiDens_i, 
	Density_j, sq_vel_j, Energy_j, SoundSpeed_j, Pressure_j, ViscDens_j, XiDens_j, dist_ij_2, dPhiE_dn;
	double Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
	Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5, 
	Sigma_yz5, Sigma_5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;	

	/*--- States in the point i ---*/
	Density_i = U_i[0];		
	sq_vel_i = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_i[iDim] = U_i[iDim+1] / Density_i;
		sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
	}
	Energy_i = U_i[nDim+1] / Density_i;
	SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-sq_vel_i));
	Pressure_i = (SoundSpeed_i * SoundSpeed_i * Density_i) / Gamma;
	ViscDens_i = (Laminar_Viscosity_i + Eddy_Viscosity_i) / Density_i;
	XiDens_i = Gamma *(Laminar_Viscosity_i/PRANDTL + Eddy_Viscosity_i/PRANDTL_TURB) / Density_i;

	/*--- States in the point j ---*/
	Density_j = U_j[0];
	sq_vel_j = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity_j[iDim] = U_j[iDim+1] / Density_j;
		sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
	Energy_j = U_j[nDim+1] / Density_j;
	SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-sq_vel_j));
	Pressure_j = (SoundSpeed_j * SoundSpeed_j * Density_j) / Gamma;
	ViscDens_j = (Laminar_Viscosity_j + Eddy_Viscosity_j) / Density_j;
	XiDens_j = Gamma *(Laminar_Viscosity_j/PRANDTL + Eddy_Viscosity_j/PRANDTL_TURB) / Density_j;

	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Mean_Velocity[iDim] = 0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
	}

	/*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge, weiss correction ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GradPsiVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradPsiVar[iVar][iDim] = 0.5*(PsiVar_Grad_i[iVar][iDim] + PsiVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradPsiVar_Edge[iVar] += Mean_GradPsiVar[iVar][iDim]*Edge_Vector[iDim];
		}
		for (iDim = 0; iDim < nDim; iDim++)
			Mean_GradPsiVar[iVar][iDim] -= (Proj_Mean_GradPsiVar_Edge[iVar] - 
					(Psi_j[iVar]-Psi_i[iVar]))*Edge_Vector[iDim]/dist_ij_2;
	}

	/*--- Average of the derivatives of the adjoint variables ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Mean_GradPsiE[iDim] = Mean_GradPsiVar[nVar-1][iDim];
		for (jDim = 0; jDim < nDim; jDim++)
			Mean_GradPhi[iDim][jDim] = Mean_GradPsiVar[iDim+1][jDim];
	}

	dPhiE_dn = 0; 
	for (iDim = 0; iDim < nDim; iDim++) 
		dPhiE_dn += Mean_GradPsiE[iDim]*Normal[iDim];

	/*--- Compute the viscous residual ---*/
	if (nDim == 3) {

		/*--- Residual at iPoint ---*/
		Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_zz = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
		Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xz = ViscDens_i * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
		Sigma_yz = ViscDens_i * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
		Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_zz5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_i[2] * Mean_GradPsiE[2]);
		Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
		Sigma_xz5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[0]);
		Sigma_yz5 = ViscDens_i * (Velocity_i[1] * Mean_GradPsiE[2] + Velocity_i[2] * Mean_GradPsiE[1]);
		Sigma_5   = XiDens_i * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
		eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

		val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy + Velocity_i[2] * Normal[2] * eta_zz
				+ (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
				+ (Velocity_i[0] * Normal[2] + Velocity_i[2] * Normal[0]) * eta_xz
				+ (Velocity_i[2] * Normal[1] + Velocity_i[1] * Normal[2]) * eta_yz
				- (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
		val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_i[0] * Sigma_5);
		val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_i[1] * Sigma_5);
		val_residual_i[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_i[2] * Sigma_5);
		val_residual_i[4] = (Sigma_5);

		/*--- Computation of the Jacobians at Point i---*/
		double dSigmaxx_phi1, dSigmayy_phi1, dSigmazz_phi1, dSigmaxy_phi1, dSigmaxz_phi1, dSigmayz_phi1;
		double dSigmaxx_phi2, dSigmayy_phi2, dSigmazz_phi2, dSigmaxy_phi2, dSigmaxz_phi2, dSigmayz_phi2;
		double dSigmaxx_phi3, dSigmayy_phi3, dSigmazz_phi3, dSigmaxy_phi3, dSigmaxz_phi3, dSigmayz_phi3;
		double dSigmaxx5_psi5, dSigmayy5_psi5, dSigmazz5_psi5, dSigmaxy5_psi5, dSigmaxz5_psi5, dSigmayz5_psi5, dSigma5_psi5;

		dSigmaxx_phi1 = -FOUR3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
		dSigmaxx_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
		dSigmaxx_phi3 =   TWO3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
		dSigmayy_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
		dSigmayy_phi2 = -FOUR3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
		dSigmayy_phi3 =   TWO3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
		dSigmazz_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
		dSigmazz_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
		dSigmazz_phi3 = -FOUR3 * ViscDens_i * Edge_Vector[2]/dist_ij_2;
		dSigmaxy_phi1 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;
		dSigmaxy_phi2 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;
		dSigmaxy_phi3 = 0;
		dSigmaxz_phi1 = -ViscDens_i * Edge_Vector[2]/dist_ij_2;
		dSigmaxz_phi2 = 0;
		dSigmaxz_phi3 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;
		dSigmayz_phi1 = 0;
		dSigmayz_phi2 = -ViscDens_i * Edge_Vector[2]/dist_ij_2;
		dSigmayz_phi3 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;

		dSigmaxx5_psi5 = -ViscDens_i * ( FOUR3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] -  TWO3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
		dSigmayy5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] + FOUR3*Velocity_i[1]*Edge_Vector[1] -  TWO3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
		dSigmazz5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] + FOUR3*Velocity_i[2]*Edge_Vector[2])/dist_ij_2;
		dSigmaxy5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[1] + Velocity_i[1]*Edge_Vector[0] )/dist_ij_2;
		dSigmaxz5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[2] + Velocity_i[2]*Edge_Vector[0] )/dist_ij_2;
		dSigmayz5_psi5 = -ViscDens_i * ( Velocity_i[1]*Edge_Vector[2] + Velocity_i[2]*Edge_Vector[1] )/dist_ij_2;
		dSigma5_psi5   = -XiDens_i * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] + Edge_Vector[2]*Normal[2] )/dist_ij_2;

		val_Jacobian_ii[0][0] = 0;
		val_Jacobian_ii[0][1] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi1 + Velocity_i[1]*Normal[1]*dSigmayy_phi1 + Velocity_i[2]*Normal[2]*dSigmazz_phi1 
				+ (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi1
				+ (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi1
				+ (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi1 );
		val_Jacobian_ii[0][2] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi2 + Velocity_i[1]*Normal[1]*dSigmayy_phi2 + Velocity_i[2]*Normal[2]*dSigmazz_phi2 
				+ (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi2
				+ (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi2
				+ (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi2 );
		val_Jacobian_ii[0][3] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi3 + Velocity_i[1]*Normal[1]*dSigmayy_phi3 + Velocity_i[2]*Normal[2]*dSigmazz_phi3 
				+ (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi3
				+ (Velocity_i[0]*Normal[2] + Velocity_i[2]*Normal[0])*dSigmaxz_phi3
				+ (Velocity_i[2]*Normal[1] + Velocity_i[1]*Normal[2])*dSigmayz_phi3 );
		val_Jacobian_ii[0][4] = (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * dSigma5_psi5;

		val_Jacobian_ii[1][0] = 0;
		val_Jacobian_ii[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1 + Normal[2]*dSigmaxz_phi1;
		val_Jacobian_ii[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2 + Normal[2]*dSigmaxz_phi2;
		val_Jacobian_ii[1][3] = Normal[0]*dSigmaxx_phi3 + Normal[1]*dSigmaxy_phi3 + Normal[2]*dSigmaxz_phi3;
		val_Jacobian_ii[1][4] = -Velocity_i[0]*dSigma5_psi5;

		val_Jacobian_ii[2][0] = 0;
		val_Jacobian_ii[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1 + Normal[2]*dSigmayz_phi1;
		val_Jacobian_ii[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2 + Normal[2]*dSigmayz_phi2;
		val_Jacobian_ii[2][3] = Normal[0]*dSigmaxy_phi3 + Normal[1]*dSigmayy_phi3 + Normal[2]*dSigmayz_phi3;
		val_Jacobian_ii[2][4] = -Velocity_i[1]*dSigma5_psi5;

		val_Jacobian_ii[3][0] = 0;
		val_Jacobian_ii[3][1] = Normal[0]*dSigmaxz_phi1 + Normal[1]*dSigmayz_phi1 + Normal[2]*dSigmazz_phi1;
		val_Jacobian_ii[3][2] = Normal[0]*dSigmaxz_phi2 + Normal[1]*dSigmayz_phi2 + Normal[2]*dSigmazz_phi2;
		val_Jacobian_ii[3][3] = Normal[0]*dSigmaxz_phi3 + Normal[1]*dSigmayz_phi3 + Normal[2]*dSigmazz_phi3;
		val_Jacobian_ii[3][4] = -Velocity_i[2]*dSigma5_psi5;

		val_Jacobian_ii[4][0] = 0;
		val_Jacobian_ii[4][1] = 0;
		val_Jacobian_ii[4][2] = 0;
		val_Jacobian_ii[4][3] = 0;
		val_Jacobian_ii[4][4] = dSigma5_psi5;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ij[iVar][jVar] = -val_Jacobian_ii[iVar][jVar];

		/*--- Residual at jPoint ---*/
		Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1] - TWO3  * Mean_GradPhi[2][2]);
		Sigma_zz = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1] + FOUR3 * Mean_GradPhi[2][2]);
		Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xz = ViscDens_j * (Mean_GradPhi[2][0] + Mean_GradPhi[0][2]);
		Sigma_yz = ViscDens_j * (Mean_GradPhi[2][1] + Mean_GradPhi[1][2]);
		Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1] -  TWO3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_zz5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1] + FOUR3 * Velocity_j[2] * Mean_GradPsiE[2]);
		Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
		Sigma_xz5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[0]);
		Sigma_yz5 = ViscDens_j * (Velocity_j[1] * Mean_GradPsiE[2] + Velocity_j[2] * Mean_GradPsiE[1]);
		Sigma_5   = XiDens_j * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
		eta_xy = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;

		val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy + Velocity_j[2] * Normal[2] * eta_zz
				+ (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
				+ (Velocity_j[0] * Normal[2] + Velocity_j[2] * Normal[0]) * eta_xz
				+ (Velocity_j[2] * Normal[1] + Velocity_j[1] * Normal[2]) * eta_yz
				- (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
		val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] + eta_xz * Normal[2] - Velocity_j[0] * Sigma_5);
		val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] + eta_yz * Normal[2] - Velocity_j[1] * Sigma_5);
		val_residual_j[3] = (eta_xz * Normal[0] + eta_yz * Normal[1] + eta_zz * Normal[2] - Velocity_j[2] * Sigma_5);
		val_residual_j[4] = (Sigma_5);

		/*--- Computation of the Jacobians at Point j---*/
		dSigmaxx_phi1 = FOUR3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
		dSigmaxx_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
		dSigmaxx_phi3 = -TWO3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
		dSigmayy_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
		dSigmayy_phi2 = FOUR3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
		dSigmayy_phi3 = -TWO3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
		dSigmazz_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
		dSigmazz_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
		dSigmazz_phi3 = FOUR3 * ViscDens_j * Edge_Vector[2]/dist_ij_2;
		dSigmaxy_phi1 = ViscDens_j * Edge_Vector[1]/dist_ij_2;
		dSigmaxy_phi2 = ViscDens_j * Edge_Vector[0]/dist_ij_2;
		dSigmaxy_phi3 = 0;
		dSigmaxz_phi1 = ViscDens_j * Edge_Vector[2]/dist_ij_2;
		dSigmaxz_phi2 = 0;
		dSigmaxz_phi3 = ViscDens_j * Edge_Vector[0]/dist_ij_2;
		dSigmayz_phi1 = 0;
		dSigmayz_phi2 = ViscDens_j * Edge_Vector[2]/dist_ij_2;
		dSigmayz_phi3 = ViscDens_j * Edge_Vector[1]/dist_ij_2;

		dSigmaxx5_psi5 = ViscDens_j * ( FOUR3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] -  TWO3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
		dSigmayy5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] + FOUR3*Velocity_j[1]*Edge_Vector[1] -  TWO3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
		dSigmazz5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] + FOUR3*Velocity_j[2]*Edge_Vector[2])/dist_ij_2;
		dSigmaxy5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[1] + Velocity_j[1]*Edge_Vector[0] )/dist_ij_2;
		dSigmaxz5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[2] + Velocity_j[2]*Edge_Vector[0] )/dist_ij_2;
		dSigmayz5_psi5 = ViscDens_j * ( Velocity_j[1]*Edge_Vector[2] + Velocity_j[2]*Edge_Vector[1] )/dist_ij_2;
		dSigma5_psi5   = XiDens_j * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] + Edge_Vector[2]*Normal[2] )/dist_ij_2;

		val_Jacobian_jj[0][0] = 0;
		val_Jacobian_jj[0][1] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi1 + Velocity_j[1]*Normal[1]*dSigmayy_phi1 + Velocity_j[2]*Normal[2]*dSigmazz_phi1 
				+ (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi1
				+ (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi1
				+ (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi1 );
		val_Jacobian_jj[0][2] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi2 + Velocity_j[1]*Normal[1]*dSigmayy_phi2 + Velocity_j[2]*Normal[2]*dSigmazz_phi2 
				+ (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi2
				+ (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi2
				+ (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi2 );
		val_Jacobian_jj[0][3] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi3 + Velocity_j[1]*Normal[1]*dSigmayy_phi3 + Velocity_j[2]*Normal[2]*dSigmazz_phi3 
				+ (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi3
				+ (Velocity_j[0]*Normal[2] + Velocity_j[2]*Normal[0])*dSigmaxz_phi3
				+ (Velocity_j[2]*Normal[1] + Velocity_j[1]*Normal[2])*dSigmayz_phi3 );
		val_Jacobian_jj[0][4] = (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * dSigma5_psi5;

		val_Jacobian_jj[1][0] = 0;
		val_Jacobian_jj[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1 + Normal[2]*dSigmaxz_phi1;
		val_Jacobian_jj[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2 + Normal[2]*dSigmaxz_phi2;
		val_Jacobian_jj[1][3] = Normal[0]*dSigmaxx_phi3 + Normal[1]*dSigmaxy_phi3 + Normal[2]*dSigmaxz_phi3;
		val_Jacobian_jj[1][4] = -Velocity_j[0]*dSigma5_psi5;

		val_Jacobian_jj[2][0] = 0;
		val_Jacobian_jj[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1 + Normal[2]*dSigmayz_phi1;
		val_Jacobian_jj[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2 + Normal[2]*dSigmayz_phi2;
		val_Jacobian_jj[2][3] = Normal[0]*dSigmaxy_phi3 + Normal[1]*dSigmayy_phi3 + Normal[2]*dSigmayz_phi3;
		val_Jacobian_jj[2][4] = -Velocity_j[1]*dSigma5_psi5;

		val_Jacobian_jj[3][0] = 0;
		val_Jacobian_jj[3][1] = Normal[0]*dSigmaxz_phi1 + Normal[1]*dSigmayz_phi1 + Normal[2]*dSigmazz_phi1;
		val_Jacobian_jj[3][2] = Normal[0]*dSigmaxz_phi2 + Normal[1]*dSigmayz_phi2 + Normal[2]*dSigmazz_phi2;
		val_Jacobian_jj[3][3] = Normal[0]*dSigmaxz_phi3 + Normal[1]*dSigmayz_phi3 + Normal[2]*dSigmazz_phi3;
		val_Jacobian_jj[3][4] = -Velocity_j[2]*dSigma5_psi5;

		val_Jacobian_jj[4][0] = 0;
		val_Jacobian_jj[4][1] = 0;
		val_Jacobian_jj[4][2] = 0;
		val_Jacobian_jj[4][3] = 0;
		val_Jacobian_jj[4][4] = dSigma5_psi5;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ji[iVar][jVar] = -val_Jacobian_jj[iVar][jVar];
	}

	if (nDim == 2) {
		/*--- Residual at iPoint ---*/
		Sigma_xx = ViscDens_i * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
		Sigma_yy = ViscDens_i * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
		Sigma_xy = ViscDens_i * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xx5 = ViscDens_i * ( FOUR3 * Velocity_i[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_i[1] * Mean_GradPsiE[1]);
		Sigma_yy5 = ViscDens_i * (- TWO3 * Velocity_i[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_i[1] * Mean_GradPsiE[1]);
		Sigma_xy5 = ViscDens_i * (Velocity_i[0] * Mean_GradPsiE[1] + Velocity_i[1] * Mean_GradPsiE[0]);
		Sigma_5   = XiDens_i * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;

		val_residual_i[0] = - (Velocity_i[0] * Normal[0] * eta_xx  + Velocity_i[1] * Normal[1] * eta_yy
				+ (Velocity_i[0] * Normal[1] + Velocity_i[1] * Normal[0]) * eta_xy
				- (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * Sigma_5);
		val_residual_i[1] = (eta_xx * Normal[0] + eta_xy * Normal[1] - Velocity_i[0] * Sigma_5);
		val_residual_i[2] = (eta_xy * Normal[0] + eta_yy * Normal[1] - Velocity_i[1] * Sigma_5);
		val_residual_i[3] = (Sigma_5);

		/*--- Computation of the Jacobians at Point i---*/
		double dSigmaxx_phi1, dSigmayy_phi1, dSigmaxy_phi1;
		double dSigmaxx_phi2, dSigmayy_phi2, dSigmaxy_phi2;
		double dSigmaxx5_psi5, dSigmayy5_psi5, dSigmaxy5_psi5, dSigma5_psi5;

		dSigmaxx_phi1 = -FOUR3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
		dSigmaxx_phi2 =   TWO3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
		dSigmayy_phi1 =   TWO3 * ViscDens_i * Edge_Vector[0]/dist_ij_2;
		dSigmayy_phi2 = -FOUR3 * ViscDens_i * Edge_Vector[1]/dist_ij_2;
		dSigmaxy_phi1 = -ViscDens_i * Edge_Vector[1]/dist_ij_2;
		dSigmaxy_phi2 = -ViscDens_i * Edge_Vector[0]/dist_ij_2;

		dSigmaxx5_psi5 = -ViscDens_i * ( FOUR3*Velocity_i[0]*Edge_Vector[0] -  TWO3*Velocity_i[1]*Edge_Vector[1] )/dist_ij_2;
		dSigmayy5_psi5 = -ViscDens_i * (- TWO3*Velocity_i[0]*Edge_Vector[0] + FOUR3*Velocity_i[1]*Edge_Vector[1] )/dist_ij_2;
		dSigmaxy5_psi5 = -ViscDens_i * ( Velocity_i[0]*Edge_Vector[1] + Velocity_i[1]*Edge_Vector[0] )/dist_ij_2;
		dSigma5_psi5   = -XiDens_i * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] )/dist_ij_2;

		val_Jacobian_ii[0][0] = 0;
		val_Jacobian_ii[0][1] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi1 + Velocity_i[1]*Normal[1]*dSigmayy_phi1  
				+ (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi1 );
		val_Jacobian_ii[0][2] = -( Velocity_i[0]*Normal[0]*dSigmaxx_phi2 + Velocity_i[1]*Normal[1]*dSigmayy_phi2  
				+ (Velocity_i[0]*Normal[1] + Velocity_i[1]*Normal[0])*dSigmaxy_phi2 );
		val_Jacobian_ii[0][3] = (sq_vel_i - Pressure_i/(Density_i*Gamma_Minus_One)) * dSigma5_psi5;

		val_Jacobian_ii[1][0] = 0;
		val_Jacobian_ii[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1;
		val_Jacobian_ii[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2;
		val_Jacobian_ii[1][3] = -Velocity_i[0]*dSigma5_psi5;

		val_Jacobian_ii[2][0] = 0;
		val_Jacobian_ii[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1;
		val_Jacobian_ii[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2;
		val_Jacobian_ii[2][3] = -Velocity_i[1]*dSigma5_psi5;

		val_Jacobian_ii[3][0] = 0;
		val_Jacobian_ii[3][1] = 0;
		val_Jacobian_ii[3][2] = 0;
		val_Jacobian_ii[3][3] = dSigma5_psi5;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ij[iVar][jVar] = -val_Jacobian_ii[iVar][jVar];

		/*--- Residual at jPoint ---*/
		Sigma_xx = ViscDens_j * (FOUR3 * Mean_GradPhi[0][0] -  TWO3 * Mean_GradPhi[1][1]);
		Sigma_yy = ViscDens_j * (-TWO3 * Mean_GradPhi[0][0] + FOUR3 * Mean_GradPhi[1][1]);
		Sigma_xy = ViscDens_j * (Mean_GradPhi[1][0] + Mean_GradPhi[0][1]);
		Sigma_xx5 = ViscDens_j * ( FOUR3 * Velocity_j[0] * Mean_GradPsiE[0] -  TWO3 * Velocity_j[1] * Mean_GradPsiE[1]);
		Sigma_yy5 = ViscDens_j * (- TWO3 * Velocity_j[0] * Mean_GradPsiE[0] + FOUR3 * Velocity_j[1] * Mean_GradPsiE[1]);
		Sigma_xy5 = ViscDens_j * (Velocity_j[0] * Mean_GradPsiE[1] + Velocity_j[1] * Mean_GradPsiE[0]);
		Sigma_5   = XiDens_j * dPhiE_dn;
		eta_xx = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_xy = Sigma_xy + Sigma_xy5;

		val_residual_j[0] = - (Velocity_j[0] * Normal[0] * eta_xx  + Velocity_j[1] * Normal[1] * eta_yy
				+ (Velocity_j[0] * Normal[1] + Velocity_j[1] * Normal[0]) * eta_xy
				- (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * Sigma_5);
		val_residual_j[1] = (eta_xx * Normal[0] + eta_xy * Normal[1]  - Velocity_j[0] * Sigma_5);
		val_residual_j[2] = (eta_xy * Normal[0] + eta_yy * Normal[1]  - Velocity_j[1] * Sigma_5);
		val_residual_j[3] = (Sigma_5);

		/*--- Computation of the Jacobians at Point j---*/
		dSigmaxx_phi1 = FOUR3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
		dSigmaxx_phi2 = -TWO3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
		dSigmayy_phi1 = -TWO3 * ViscDens_j * Edge_Vector[0]/dist_ij_2;
		dSigmayy_phi2 = FOUR3 * ViscDens_j * Edge_Vector[1]/dist_ij_2;
		dSigmaxy_phi1 = ViscDens_j * Edge_Vector[1]/dist_ij_2;
		dSigmaxy_phi2 = ViscDens_j * Edge_Vector[0]/dist_ij_2;

		dSigmaxx5_psi5 = ViscDens_j * ( FOUR3*Velocity_j[0]*Edge_Vector[0] -  TWO3*Velocity_j[1]*Edge_Vector[1] )/dist_ij_2;
		dSigmayy5_psi5 = ViscDens_j * (- TWO3*Velocity_j[0]*Edge_Vector[0] + FOUR3*Velocity_j[1]*Edge_Vector[1] )/dist_ij_2;
		dSigmaxy5_psi5 = ViscDens_j * ( Velocity_j[0]*Edge_Vector[1] + Velocity_j[1]*Edge_Vector[0] )/dist_ij_2;
		dSigma5_psi5   = XiDens_j * ( Edge_Vector[0]*Normal[0] + Edge_Vector[1]*Normal[1] )/dist_ij_2;

		val_Jacobian_jj[0][0] = 0;
		val_Jacobian_jj[0][1] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi1 + Velocity_j[1]*Normal[1]*dSigmayy_phi1  
				+ (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi1 );
		val_Jacobian_jj[0][2] = -( Velocity_j[0]*Normal[0]*dSigmaxx_phi2 + Velocity_j[1]*Normal[1]*dSigmayy_phi2  
				+ (Velocity_j[0]*Normal[1] + Velocity_j[1]*Normal[0])*dSigmaxy_phi2 );
		val_Jacobian_jj[0][3] = (sq_vel_j - Pressure_j/(Density_j*Gamma_Minus_One)) * dSigma5_psi5;

		val_Jacobian_jj[1][0] = 0;
		val_Jacobian_jj[1][1] = Normal[0]*dSigmaxx_phi1 + Normal[1]*dSigmaxy_phi1;
		val_Jacobian_jj[1][2] = Normal[0]*dSigmaxx_phi2 + Normal[1]*dSigmaxy_phi2;
		val_Jacobian_jj[1][3] = -Velocity_j[0]*dSigma5_psi5;

		val_Jacobian_jj[2][0] = 0;
		val_Jacobian_jj[2][1] = Normal[0]*dSigmaxy_phi1 + Normal[1]*dSigmayy_phi1;
		val_Jacobian_jj[2][2] = Normal[0]*dSigmaxy_phi2 + Normal[1]*dSigmayy_phi2;
		val_Jacobian_jj[2][3] = -Velocity_j[1]*dSigma5_psi5;

		val_Jacobian_jj[3][0] = 0;
		val_Jacobian_jj[3][1] = 0;
		val_Jacobian_jj[3][2] = 0;
		val_Jacobian_jj[3][3] = dSigma5_psi5;

		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ji[iVar][jVar] = -val_Jacobian_jj[iVar][jVar];
	}
}

CAvgGradCorrectedArtComp_AdjFlow::CAvgGradCorrectedArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	Mean_Velocity = new double [nDim];
	
	Mean_GradPsiVar = new double* [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Mean_GradPsiVar[iVar] = new double [nDim];
	
	Edge_Vector = new double [nDim];
	Proj_Mean_GradPsiVar_Edge = new double [nVar];
	
	Mean_GradPhi = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Mean_GradPhi[iDim] = new double [nDim];
	Mean_GradPsiE = new double [nDim];
	
}

CAvgGradCorrectedArtComp_AdjFlow::~CAvgGradCorrectedArtComp_AdjFlow(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Mean_Velocity;
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradPsiVar_Edge;
	
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradPsiVar[iVar];
	delete [] Mean_GradPsiVar;
	
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		delete [] Mean_GradPhi[iDim];
	delete [] Mean_GradPhi;
	delete [] Mean_GradPsiE;
}

void CAvgGradCorrectedArtComp_AdjFlow::SetResidual(double *val_residual_i, double *val_residual_j) {

}

void CAvgGradCorrectedArtComp_AdjFlow::SetResidual(double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii, double **val_Jacobian_ij,
																						double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) {

}

CAvgGradCorrected_AdjTurb::CAvgGradCorrected_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	Edge_Vector = new double [nDim];
	Mean_GradTurbPsi = new double* [nVar];
	Proj_Mean_GradTurbPsi_Kappa = new double [nVar];
	Proj_Mean_GradTurbPsi_Edge = new double [nVar];
	Proj_Mean_GradTurbPsi_Corrected = new double [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbPsi[iVar] = new double [nDim];
}

CAvgGradCorrected_AdjTurb::~CAvgGradCorrected_AdjTurb(void) {
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbPsi_Kappa;
	delete [] Proj_Mean_GradTurbPsi_Edge;
	delete [] Proj_Mean_GradTurbPsi_Corrected;
	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		delete [] Mean_GradTurbPsi[iVar];
	}
	delete [] Mean_GradTurbPsi;
}

void CAvgGradCorrected_AdjTurb::SetResidual (double *val_residual, double **val_Jacobian_i,
		double **val_Jacobian_j, CConfig *config) {

	bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);
	bool WEISS_CORR = true; 

	double sigma = 2./3.;
	double nu_i, nu_j, nu_e;
	double dist_ij_2 = 0; // squared distance btw points
	double proj_vector_ij = 0; // projection of edge vector on Kappa
	unsigned short iVar, iDim;

	// Compute mean effective viscosity
	nu_i = Laminar_Viscosity_i/U_i[0];
	nu_j = Laminar_Viscosity_j/U_j[0];
	nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0])/sigma;

	// Compute vector going from iPoint to jPoint
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
		proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
	}
	proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors

	// Mean gradient approximation
	for (iVar = 0; iVar < nVar; iVar++) {
		// Correction of the mean gradient in the direction of the edge 
		Proj_Mean_GradTurbPsi_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbPsi[iVar][iDim] = 0.5*(TurbPsi_Grad_i[iVar][iDim] + TurbPsi_Grad_j[iVar][iDim]);
			if (WEISS_CORR) Proj_Mean_GradTurbPsi_Edge[iVar] += Mean_GradTurbPsi[iVar][iDim]*Edge_Vector[iDim];
		}
		if (WEISS_CORR)
			for (iDim = 0; iDim < nDim; iDim++)
				Mean_GradTurbPsi[iVar][iDim] -= (Proj_Mean_GradTurbPsi_Edge[iVar] - 
						(TurbPsi_j[iVar]-TurbPsi_i[iVar]))*Edge_Vector[iDim]/dist_ij_2;
		// Projection of the corrected gradient
		Proj_Mean_GradTurbPsi_Corrected[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Proj_Mean_GradTurbPsi_Corrected[iVar] += Mean_GradTurbPsi[iVar][iDim]*Normal[iDim];
	}

	val_residual[0] = -nu_e*Proj_Mean_GradTurbPsi_Corrected[0];

	if (implicit) {
		// For Jacobians -> Use of TSL approx. to compute derivatives of the gradients
		val_Jacobian_i[0][0] =  nu_e*proj_vector_ij;
		val_Jacobian_j[0][0] = -nu_e*proj_vector_ij;
	}
}

void CAvgGradCorrected_AdjTurb::SetResidual (double *val_residual_i, double *val_residual_j,
		double **val_Jacobian_ii, double **val_Jacobian_ij,
		double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) {

	bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);
	bool WEISS_CORR = true; 

	double sigma = 2./3.;
	double nu_i, nu_j, nu_e_i, nu_e_j;
	double dist_ij_2 = 0; // squared distance btw points
	double proj_vector_ij = 0; // projection of edge vector on Kappa
	unsigned short iVar, iDim;

	// Compute mean effective viscosity
	nu_i = Laminar_Viscosity_i/U_i[0];
	nu_j = Laminar_Viscosity_j/U_j[0];
	nu_e_i = (nu_i+TurbVar_i[0])/sigma;
	nu_e_j = (nu_j+TurbVar_j[0])/sigma;

	// Compute vector going from iPoint to jPoint
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
		proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
	}
	proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors

	// Mean gradient approximation
	for (iVar = 0; iVar < nVar; iVar++) {
		// Correction of the mean gradient in the direction of the edge 
		Proj_Mean_GradTurbPsi_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbPsi[iVar][iDim] = 0.5*(TurbPsi_Grad_i[iVar][iDim] + TurbPsi_Grad_j[iVar][iDim]);
			if (WEISS_CORR) Proj_Mean_GradTurbPsi_Edge[iVar] += Mean_GradTurbPsi[iVar][iDim]*Edge_Vector[iDim];
		}
		if (WEISS_CORR)
			for (iDim = 0; iDim < nDim; iDim++)
				Mean_GradTurbPsi[iVar][iDim] -= (Proj_Mean_GradTurbPsi_Edge[iVar] - 
						(TurbPsi_j[iVar]-TurbPsi_i[iVar]))*Edge_Vector[iDim]/dist_ij_2;
		// Projection of the corrected gradient
		Proj_Mean_GradTurbPsi_Corrected[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Proj_Mean_GradTurbPsi_Corrected[iVar] += Mean_GradTurbPsi[iVar][iDim]*Normal[iDim];
	}

	val_residual_i[0] = -nu_e_i*Proj_Mean_GradTurbPsi_Corrected[0];
	val_residual_j[0] =  nu_e_j*Proj_Mean_GradTurbPsi_Corrected[0];

	if (implicit) {
		// For Jacobians -> Use of TSL approx. to compute derivatives of the gradients
		val_Jacobian_ii[0][0] =  nu_e_i*proj_vector_ij;
		val_Jacobian_ij[0][0] = -nu_e_i*proj_vector_ij;
		val_Jacobian_ji[0][0] = -nu_e_j*proj_vector_ij;
		val_Jacobian_jj[0][0] =  nu_e_j*proj_vector_ij;
	}
}

CGalerkin_Flow::CGalerkin_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

CGalerkin_Flow::~CGalerkin_Flow(void) { }

void CGalerkin_Flow::SetResidual(double **val_stiffmatrix_elem, CConfig *config) {

	double a[4], b[4], c[4], Area;

	if (nDim == 2) {

		for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}

		Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);	/* Norm of the normal component of area, area = 1/2*cross(a,b) */

		a[0] = 0.5 * (Coord_1[0]*Coord_2[1]-Coord_2[0]*Coord_1[1]) / Area;
		a[1] = 0.5 * (Coord_2[0]*Coord_0[1]-Coord_0[0]*Coord_2[1]) / Area;
		a[2] = 0.5 * (Coord_0[0]*Coord_1[1]-Coord_1[0]*Coord_0[1]) / Area;

		b[0] = 0.5 * (Coord_1[1]-Coord_2[1]) / Area;
		b[1] = 0.5 * (Coord_2[1]-Coord_0[1]) / Area;
		b[2] = 0.5 * (Coord_0[1]-Coord_1[1]) / Area;

		c[0] = 0.5 * (Coord_2[0]-Coord_1[0]) / Area;
		c[1] = 0.5 * (Coord_0[0]-Coord_2[0]) / Area;
		c[2] = 0.5 * (Coord_1[0]-Coord_0[0]) / Area;

		/* Compute the stiffness matrix, K & multiply it by the Area */

		val_stiffmatrix_elem[0][0] = Area * (b[0]*b[0]+c[0]*c[0]);	
		val_stiffmatrix_elem[0][1] = Area * (b[0]*b[1]+c[0]*c[1]);	
		val_stiffmatrix_elem[0][2] = Area * (b[0]*b[2]+c[0]*c[2]);	
		val_stiffmatrix_elem[1][0] = Area * (b[0]*b[1]+c[0]*c[1]);	
		val_stiffmatrix_elem[1][1] = Area * (b[1]*b[1]+c[1]*c[1]);	
		val_stiffmatrix_elem[1][2] = Area * (b[1]*b[2]+c[1]*c[2]);	
		val_stiffmatrix_elem[2][0] = Area * (b[0]*b[2]+c[0]*c[2]);	
		val_stiffmatrix_elem[2][1] = Area * (b[1]*b[2]+c[1]*c[2]);	
		val_stiffmatrix_elem[2][2] = Area * (b[2]*b[2]+c[2]*c[2]);
	}

	if (nDim == 3) {
	}
}

CAvgGrad_Plasma::CAvgGrad_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nFluids, CConfig *config) : CNumerics(val_nDim, val_nVar,val_nSpecies,val_nFluids, config) {

	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Allocate some useful vectors ---*/
	Prim_Var_i = new double [nVar];
	Prim_Var_j = new double [nVar];
	Mean_Laminar_Viscosity = new double [nSpecies];
	Mean_Eddy_Viscosity = new double [nSpecies];
	Mean_PrimVar = new double [nVar];
	Mean_Density = new double [nSpecies];
	Mean_GradPrimVar = new double* [nVar+nSpecies];
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGrad_Plasma::~CAvgGrad_Plasma(void) {

	delete [] Prim_Var_i;
	delete [] Prim_Var_j;
	delete [] Mean_PrimVar;
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
}

void CAvgGrad_Plasma::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	unsigned short loc, iSpecies;


	/*--- Normalized normal vector ---*/
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);

	for (iDim = 0; iDim < nDim; iDim++)
		UnitaryNormal[iDim] = Normal[iDim]/Area;

	/*--- Conversion to Primitive Variables ---*/
	ConsVar2PrimVar_MultiSpecies(U_i, Prim_Var_i);
	ConsVar2PrimVar_MultiSpecies(U_j, Prim_Var_j);

	/*--- Mean Primitive variable value ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_PrimVar[iVar] = 0.5*(Prim_Var_i[iVar]+Prim_Var_j[iVar]);

	/*--- Mean Viscosities ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		Mean_Laminar_Viscosity[iSpecies] = 0.5*(Laminar_Viscosity_MultipleSpecies_i[iSpecies] + Laminar_Viscosity_MultipleSpecies_j[iSpecies]);
		Mean_Eddy_Viscosity[iSpecies] = 0.5*(Eddy_Viscosity_MultipleSpecies_i[iSpecies] + Eddy_Viscosity_MultipleSpecies_j[iSpecies]);
	}

	/*--- Mean gradient approximation ---*/
	for (iVar = 0; iVar < nVar + nSpecies; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);

	/*--- Get projected flux tensor ---*/
	GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

	/*--- Update viscous residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Proj_Flux_Tensor[iVar];


	/*--- Implicit part ---*/
	if (implicit) {
		dist_ij = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
		dist_ij = sqrt(dist_ij);
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			loc = (nDim+2)* iSpecies;
			Mean_Density[iSpecies] = 0.5*(U_i[loc + 0]+U_j[loc + 0]);
		}
		GetViscousProjJacs(Mean_Density, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
				dist_ij, UnitaryNormal, Area, Mean_PrimVar, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
	}

}

CViscous_Template::CViscous_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CViscous_Template::~CViscous_Template(void) { }

void CViscous_Template::SetResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) { }
