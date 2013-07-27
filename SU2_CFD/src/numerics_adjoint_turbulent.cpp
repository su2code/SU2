/*!
 * \file numerics_adjoint_turbulent.cpp
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

CUpwLin_AdjTurb::CUpwLin_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Velocity_i = new double [nDim];
}

CUpwLin_AdjTurb::~CUpwLin_AdjTurb(void) {
	delete [] Velocity_i;
}

void CUpwLin_AdjTurb::ComputeResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
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

void CUpwSca_AdjTurb::ComputeResidual (double *val_residual_i, double *val_residual_j,
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

void CAvgGradCorrected_AdjTurb::ComputeResidual (double *val_residual, double **val_Jacobian_i,
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

void CAvgGradCorrected_AdjTurb::ComputeResidual (double *val_residual_i, double *val_residual_j,
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

CSourcePieceWise_AdjTurb::CSourcePieceWise_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Velocity = new double [nDim];
	tau = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		tau[iDim] = new double [nDim];
}

CSourcePieceWise_AdjTurb::~CSourcePieceWise_AdjTurb(void) {
	delete [] Velocity;
  
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		delete [] tau[iDim];
	delete [] tau;
}

void CSourcePieceWise_AdjTurb::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	unsigned short iDim, jDim;
	bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);
  
	val_residual[0] = 0.0;
	if (implicit)
		val_Jacobian_i[0][0] = 0.0;
  
	if (dist_i > 0.0) {
    
		/*--- Computation of Vorticity and Divergence of velocity ---*/
		double div_vel = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity[iDim] = U_i[iDim+1]/U_i[0];
			div_vel += PrimVar_Grad_i[iDim+1][iDim];
		}
    
		double Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
		if (nDim == 3)
			Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) +
                    (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );
		Vorticity = sqrt(Vorticity);
    
		/*--- FIRST PART: -Bs*TurbPsi_i ---*/
		/*--- CLOUSURE CONSTANTS ---*/
		double cv1 = 7.1;
		double cv1_3 = cv1*cv1*cv1;
		double k = 0.41;
		double k2 = k*k;
		double cb1 = 0.1355;
		double cw2 = 0.3;
		double cw3_6 = pow(2.0,6.0);
		double sigma = 2./3.;
		double cb2 = 0.622;
		double cw1 = cb1/k2+(1+cb2)/sigma;
    
		double nu, Ji, fv1, fv2, Shat, dist_0_2, Ji_2, Ji_3, one_o_oneplusJifv1;
		double r, g, g_6, glim, fw;
		double dTs_nuhat, dTs_Shat, dShat_nuhat, dTs_fw, dfw_g, dg_r, dr_nuhat, dr_Shat;
		double dShat_fv2, dfv2_fv1, dfv1_Ji, dJi_nuhat, dfv2_Ji;
		double Bs;
    
		dist_0_2 = dist_i*dist_i;
		nu = Laminar_Viscosity_i/U_i[0];
		Ji = TurbVar_i[0]/nu;
		Ji_2 = Ji*Ji;
		Ji_3 = Ji_2*Ji;
		fv1 = Ji_3/(Ji_3+cv1_3);
		one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
		fv2 = 1.0 - Ji*one_o_oneplusJifv1;
		Shat = max(Vorticity + TurbVar_i[0]*fv2/(k2*dist_0_2),TURB_EPS);
    
		//		r = TurbVar_i[0]/(Shat*k2*dist_0_2);
		r = min(TurbVar_i[0]/(Shat*k2*dist_0_2),10.);
		g = r + cw2*(pow(r,6.)-r);
		g_6 = pow(g,6.);
		glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
		fw = g*glim;
    
		dTs_nuhat = cb1*Shat-2.0*cw1*fw*TurbVar_i[0]/dist_0_2;
		dTs_Shat = cb1*TurbVar_i[0];
		dTs_fw = -cw1*TurbVar_i[0]*TurbVar_i[0]/dist_0_2;
		dfw_g  = glim*cw3_6/(g_6+cw3_6);
		dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
		dr_nuhat = 1.0/(Shat*k2*dist_0_2);
		dr_Shat = -dr_nuhat*TurbVar_i[0]/Shat;
    
		dShat_nuhat = fv2/(k2*dist_0_2);
		dShat_fv2 = TurbVar_i[0]/(k2*dist_0_2);
		dfv2_fv1 = Ji_2*one_o_oneplusJifv1*one_o_oneplusJifv1;
		dfv1_Ji = 3.0*cv1_3*Ji_2/((Ji_3+cv1_3)*(Ji_3+cv1_3));
		dJi_nuhat = 1.0/nu;
		dfv2_Ji = -one_o_oneplusJifv1*one_o_oneplusJifv1;
		dShat_nuhat += dShat_fv2*(dfv2_fv1*dfv1_Ji+dfv2_Ji)*dJi_nuhat;
    
		Bs = dTs_nuhat;											 // nu_hat term
		Bs += dTs_Shat*dShat_nuhat;								 // S_hat term
		Bs += dTs_fw*dfw_g*dg_r*(dr_nuhat+dr_Shat*dShat_nuhat);	 // fw terms
    
		val_residual[0] = -Bs*TurbPsi_i[0]*Volume;
    
    //        cout << val_residual[0] << endl;
    //        cin.get();
    
		if (implicit)
			val_Jacobian_i[0][0] = -Bs*Volume;
    
		/*---SECOND PART: \partial_nu_hat mu^k F^{vk} cdot \grad Psi ---*/
		double dEddyVisc_nuhat;
    if (!config->GetFrozen_Visc())
      dEddyVisc_nuhat = U_i[0]*fv1*(1.0 + 3.0*cv1_3/(Ji_3+cv1_3));
    else
      dEddyVisc_nuhat = 0;
    
		for (iDim = 0; iDim < nDim; iDim++) {
			for (jDim = 0; jDim < nDim; jDim++)
				tau[iDim][jDim] = PrimVar_Grad_i[iDim+1][jDim] + PrimVar_Grad_i[jDim+1][iDim];
			tau[iDim][iDim] -= TWO3*div_vel;
		}
    
		double Gas_Constant = config->GetGas_ConstantND();
		double Cp = (Gamma/Gamma_Minus_One)*Gas_Constant;
		double tau_gradphi = 0.0, vel_tau_gradpsi5 = 0.0, gradT_gradpsi5 = 0.0;
    
		for (iDim = 0; iDim < nDim; iDim++) {
			gradT_gradpsi5 += PrimVar_Grad_i[0][iDim]*PsiVar_Grad_i[nVar-1][iDim];
			for (jDim = 0; jDim < nDim; jDim++) {
				tau_gradphi += tau[iDim][jDim]*PsiVar_Grad_i[iDim+1][jDim];
				vel_tau_gradpsi5 += Velocity[iDim]*tau[iDim][jDim]*PsiVar_Grad_i[nVar-1][jDim];
			}
		}
		val_residual[0] += (tau_gradphi + vel_tau_gradpsi5 + Cp/PRANDTL_TURB*gradT_gradpsi5)*dEddyVisc_nuhat*Volume;
		// no contributions to the Jacobians since this term does not depend on TurbPsi
    
    //        cout << (tau_gradphi + vel_tau_gradpsi5 + Cp/PRANDTL_TURB*gradT_gradpsi5)*dEddyVisc_nuhat*Volume << endl;
    //        cin.get();
    //
    //        cout << tau_gradphi << endl;
    //        cin.get();
    //
    //        cout << Volume << endl;
    //        cin.get();
    //
    //        cout << val_residual[0] << endl;
    //        cin.get();
    
    //        for (unsigned short iVar = 1; iVar < nDim+1; iVar++)
    //            for (iDim = 0; iDim < nDim; iDim++)
    //                cout << "iVar: " << iVar << ", iDim: " << iDim << ", PsiVar_Grad_i: " << PsiVar_Grad_i[iVar][iDim] << endl;
    //        cin.get();
    
	}
}

CSourceConservative_AdjTurb::CSourceConservative_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
}

CSourceConservative_AdjTurb::~CSourceConservative_AdjTurb(void) {
}

void CSourceConservative_AdjTurb::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	/*--- SOURCE term  -->  \nabla ( \psi_\mu \B7 E^{s} )
	 E^{s} = 2 c_{b2}/\sigma \nabla \hat{nu} ---*/
  
	unsigned short iDim;
	bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);
  
	double cb2 = 0.622;
	double sigma = 2./3.;
	double coeff = 2.0*cb2/sigma;
	double E_ij, proj_TurbVar_Grad_i, proj_TurbVar_Grad_j;
  
	E_ij = 0.0;	proj_TurbVar_Grad_i = 0.0; proj_TurbVar_Grad_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		proj_TurbVar_Grad_i += coeff*TurbVar_Grad_i[0][iDim]*Normal[iDim];
		proj_TurbVar_Grad_j += coeff*TurbVar_Grad_j[0][iDim]*Normal[iDim];
		E_ij += 0.5*(TurbPsi_i[0]*proj_TurbVar_Grad_i + TurbPsi_j[0]*proj_TurbVar_Grad_j);
	}
  
	val_residual[0] = E_ij;
  
	if (implicit) {
		val_Jacobian_i[0][0] = 0.5*proj_TurbVar_Grad_i;
		val_Jacobian_j[0][0] = 0.5*proj_TurbVar_Grad_j;
	}
}
