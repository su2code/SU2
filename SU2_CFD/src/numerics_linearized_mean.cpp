/*!
 * \file numerics_linearized_mean.cpp
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

void CCentJST_LinFlow::ComputeResidual (double *val_resconv, double *val_resvisc, double **val_Jacobian_i,
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

void CCentLax_LinFlow::ComputeResidual (double *val_resconv, double *val_resvisc, double **val_Jacobian_i,
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
