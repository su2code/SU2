/*!
 * \file numerics_source.cpp
 * \brief This file contains all the source term discretization.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

CSourceNothing::CSourceNothing(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

CSourceNothing::~CSourceNothing(void) { }

CSourcePieceWise_TurbSA::CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	//	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Closure constants ---*/
	cv1_3 = pow(7.1,3.0);
	k2 = pow(0.41,2.0);
	cb1 = 0.1355;
	cw2 = 0.3;
	cw3_6 = pow(2.0,6.0);
	sigma = 2./3.;
	cb2 = 0.622;
	cw1 = cb1/k2+(1+cb2)/sigma;
}

CSourcePieceWise_TurbSA::~CSourcePieceWise_TurbSA(void) { }

void CSourcePieceWise_TurbSA::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	val_residual[0] = 0.0;
	if (implicit)
		val_Jacobian_i[0][0] = 0.0;

	/*--- Computation of divergence of velocity and vorticity ---*/
	DivVelocity = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		DivVelocity += PrimVar_Grad_i[iDim+1][iDim];

	Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
	if (nDim == 3)
		Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) +
				(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );
	Vorticity = sqrt(Vorticity);

	switch (config->GetKind_Turb_Model()) {
	case SA :
		if (dist_i > 0.0) {

			/*--- Production term ---*/
			dist_0_2 = dist_i*dist_i;
			nu = Laminar_Viscosity_i/U_i[0];
			Ji = TurbVar_i[0]/nu;
			Ji_2 = Ji*Ji;
			Ji_3 = Ji_2*Ji;
			fv1 = Ji_3/(Ji_3+cv1_3);
			fv2 = 1.0 - Ji/(1.0+Ji*fv1);
			Omega = Vorticity;
			Shat = max(Omega + TurbVar_i[0]*fv2/(k2*dist_0_2),TURB_EPS);
			val_residual[0] += cb1*Shat*TurbVar_i[0]*Volume;

			/*--- Destruction term ---*/
			r = min(TurbVar_i[0]/(Shat*k2*dist_0_2),10.);
			g = r + cw2*(pow(r,6.)-r);
			g_6 = pow(g,6.);
			glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
			fw = g*glim;
			val_residual[0] -= cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_0_2*Volume;

			/*--- Diffusion term ---*/
			norm2_Grad = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
			val_residual[0] += cb2/sigma*norm2_Grad*Volume;

			/*--- Implicit part ---*/
			if (implicit) {

				/*--- Production term ---*/
				dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
				dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
				if ( Shat <= TURB_EPS )
					dShat = 0.0;
				else
					dShat = (fv2+TurbVar_i[0]*dfv2)/(k2*dist_0_2);
				val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;

				/*--- Destruction term ---*/
				dr = (Shat-TurbVar_i[0]*dShat)/(Shat*Shat*k2*dist_0_2);
				if (r == 10.) dr = 0.0;
				dg = dr*(1.+cw2*(6.*pow(r,5.)-1.));
				dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
				val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] +	2.*fw)*TurbVar_i[0]/dist_0_2*Volume;
			}
		}
		break;

	case SA_COMP :

		if (dist_i > 0.0) {

			/*--- Production term ---*/
			dist_0_2 = dist_i*dist_i;
			Ji = TurbVar_i[0]/Laminar_Viscosity_i;
			Ji_2 = Ji*Ji;
			Ji_3 = Ji_2*Ji;
			fv1 = Ji_3/(Ji_3+cv1_3);
			fv2 = 1.0 - Ji/(1+Ji*fv1);
			Omega = Vorticity;
			Shat = max(Omega + TurbVar_i[0]*fv2/(U_i[0]*k2*dist_0_2),TURB_EPS);
			val_residual[0] += cb1*Shat*TurbVar_i[0]*Volume;

			/*--- Destruction term ---*/
			r = min(TurbVar_i[0]/(U_i[0]*Shat*k2*dist_0_2),10.);
			g = r + cw2*(pow(r,6.)-r);
			g_6 = pow(g,6.);
			glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
			fw = g*glim;
			val_residual[0] -= cw1*fw*TurbVar_i[0]*TurbVar_i[0]/(U_i[0]*dist_0_2)*Volume;

			/*--- Diffusion term ---*/
			nu_hat_i = TurbVar_i[0]/U_i[0];
			norm2_Grad = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				norm2_Grad += pow(TurbVar_Grad_i[0][iDim]-ConsVar_Grad_i[0][iDim]*nu_hat_i,2.0);
			val_residual[0] += cb2/(sigma*U_i[0])*norm2_Grad*Volume;

			/*--- Implicit part ---*/
			if (implicit) {

				/*--- Production term ---*/
				dfv1 = 3.0*Ji_2*cv1_3/(Laminar_Viscosity_i*pow(Ji_3+cv1_3,2.));
				dfv2 = -(1/Laminar_Viscosity_i-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
				if ( Shat <= TURB_EPS )
					dShat = 0.0;
				else
					dShat = (fv2+TurbVar_i[0]*dfv2)/(k2*dist_0_2);
				val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;

				/*--- Destruction term ---*/
				double dr, dg, dfw;;
				dr = (Shat-TurbVar_i[0]*dShat)/(Shat*Shat*U_i[0]*k2*dist_0_2);
				if (r == 10.) dr = 0.0;
				dg = dr*(1.+cw2*(6.*pow(r,5.)-1.));
				dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
				val_Jacobian_i[0][0] -= cw1/U_i[0]*(dfw*TurbVar_i[0] + 2.*fw)*TurbVar_i[0]/dist_0_2*Volume;

				/*--- Diffusion term ---*/
				prod_grads = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					grad_nu_hat = TurbVar_Grad_i[0][iDim]-ConsVar_Grad_i[0][iDim]*nu_hat_i;
					prod_grads += grad_nu_hat*ConsVar_Grad_i[0][iDim]/U_i[0];
				}
				val_Jacobian_i[0][0] -= 2.0*cb2/(sigma*U_i[0])*prod_grads*Volume;
			}
		}
		break;
	}
}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	//	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Closure constants ---*/
	beta_star  = 0.09;
	von_Karman = 0.41;
	sigma_omega_1 = 0.5;
	sigma_omega_2 = 0.856;
	beta_1 = 0.0750;
	beta_2 = 0.0828;
	gamma_1 = beta_1/beta_star-sigma_omega_1*von_Karman*von_Karman/sqrt(beta_star);
	gamma_2 = beta_2/beta_star-sigma_omega_2*von_Karman*von_Karman/sqrt(beta_star);

	/*--- --- Allocate memory for Strain rate and Reynolds stress tensors ---*/
	strain_rate = new double *[nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		strain_rate[iDim] = new double [nDim];

	reynolds_stress= new double *[nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		reynolds_stress[iDim] = new double [nDim];
}

CSourcePieceWise_TurbSST::~CSourcePieceWise_TurbSST(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		delete [] strain_rate[iDim];
	delete [] strain_rate;

	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		delete [] reynolds_stress[iDim];
	delete [] reynolds_stress;
}

void CSourcePieceWise_TurbSST::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	unsigned short iDim, jDim;
	double gamma_blended, beta_blended;
	double production = 0.0, cross_difusion = 0.0;

	val_residual[0] = 0.0;
	val_residual[1] = 0.0;
	if (implicit){
		val_Jacobian_i[0][0] = 0.0;		val_Jacobian_i[0][1] = 0.0;
		val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = 0.0;
	}

	/*--- Computation of blended constants for the source terms---*/
	gamma_blended = F1*gamma_1+(1-F1)*gamma_2;
	beta_blended  = F1*beta_1+(1-F1)*beta_2;

	if (dist_i > 0.0) {
		/*--- Production term ---*/
		DivVelocity = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			DivVelocity += PrimVar_Grad_i[iDim+1][iDim];

		for (iDim = 0; iDim < nDim; iDim++)
			for(jDim = 0; jDim < nDim; jDim++)
				strain_rate[iDim][jDim] = 0.5*(PrimVar_Grad_i[iDim+1][jDim]+PrimVar_Grad_i[jDim+1][iDim]);

		for (iDim = 0; iDim < nDim; iDim++)
			for (jDim = 0; jDim < nDim; jDim++)
				reynolds_stress[iDim][jDim] = 2*Eddy_Viscosity_i*strain_rate[iDim][jDim];
		reynolds_stress[0][0] -= 2.0/3.0*Eddy_Viscosity_i*DivVelocity + 2.0/3.0*U_i[0]*TurbVar_i[0];
		reynolds_stress[1][1] -= 2.0/3.0*Eddy_Viscosity_i*DivVelocity + 2.0/3.0*U_i[0]*TurbVar_i[0];
		reynolds_stress[2][2] -= 2.0/3.0*Eddy_Viscosity_i*DivVelocity + 2.0/3.0*U_i[0]*TurbVar_i[0];

		for (iDim = 0; iDim < nDim; iDim++)
			for (jDim = 0; jDim < nDim; jDim++)
				production += reynolds_stress[iDim][jDim]*PrimVar_Grad_i[iDim+1][jDim];
		production = min(production,20*beta_star*U_i[0]*TurbVar_i[1]*TurbVar_i[0]);

		val_residual[0] += production*Volume;
		val_residual[1] += gamma_blended/(Eddy_Viscosity_i/U_i[0])*production*Volume;

		/*--- Destruction term ---*/
		val_residual[0] -= beta_star*U_i[0]*TurbVar_i[1]*TurbVar_i[0]*Volume;
		val_residual[1] -= beta_blended*U_i[0]*TurbVar_i[1]*TurbVar_i[1]*Volume;

		/*--- Diffusion term ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			cross_difusion += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[1][iDim];
		cross_difusion *= 2*U_i[0]*(1-F1)*sigma_omega_2/TurbVar_i[1];
		val_residual[1] += cross_difusion*Volume;

		/*--- Implicit part ---*/
		if (implicit) {
			val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]/U_i[0];		val_Jacobian_i[0][1] = 0.0;
			val_Jacobian_i[1][0] = 0.0;									val_Jacobian_i[1][1] = -2*beta_blended*TurbVar_i[1]/U_i[0]-abs(cross_difusion)/TurbVar_i[1];
		}
	}

}

CSourcePieceWise_Gravity::CSourcePieceWise_Gravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	incompressible = (config->GetIncompressible() == YES);

	U_ref = config->GetVelocity_Ref();
	L_ref = config->GetLength_Ref();

}

CSourcePieceWise_Gravity::~CSourcePieceWise_Gravity(void) { }

void CSourcePieceWise_Gravity::SetResidual(double *val_residual, CConfig *config) {
	unsigned short iVar;

	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = 0.0;

	if (incompressible) {
		/*--- Compute the Froude number  ---*/
		Froude = U_ref / sqrt( L_ref * STANDART_GRAVITY);

		/*--- Evaluate the dource term  ---*/
		if (nDim == 2) val_residual[2] = 0.0; //Volume / (Froude * Froude);
		if (nDim == 3) val_residual[3] = 0.0; //Volume / (Froude * Froude);

	}

}

CSourcePieceWise_Elec::CSourcePieceWise_Elec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

CSourcePieceWise_Elec::~CSourcePieceWise_Elec(void) { }

void CSourcePieceWise_Elec::SetResidual(double *val_residual, CConfig *config) {

	double Kai_n, Kai_np1 = 0.0, ne, np;
	double diff_ru_Elec_x, diff_ru_Pos_x, diff_ru_Elec_y, diff_ru_Pos_y, diff_ru_Pos_z, diff_ru_Elec_z;
	double alpha;
	double dt = TimeStep;
	double rho_Pos = 0.0, rho_Elec = 0.0, mass_Elec, mass_Pos;
	double a[4], b[4], Area = 0.0;

	if (nDim == 2) {
		for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);			/*--- Norm of the normal component of area, area = 1/2*cross(a,b) ---*/
	}

	if (nDim == 2) 	{	
		rho_Pos  = 1.0/3.0*(U_0[4] + U_1[4] + U_2[4]) ; 
		rho_Elec = 1.0/3.0*(U_0[8] + U_1[8] + U_2[8]) ; 
	}

	if (nDim == 3) 	{	
		rho_Pos  = 1.0/3.0*(U_0[5] + U_1[5] + U_2[5]) ;  
		rho_Elec = 1.0/3.0*(U_0[10] + U_1[10] + U_2[10]) ; 
	}

	/*--- Source q ---*/
	mass_Elec = config->GetParticle_Mass(2);
	mass_Pos  = config->GetParticle_Mass(1);

	ne = rho_Elec / mass_Elec;
	np = rho_Pos  / mass_Pos;

	Kai_n = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne - np );
	alpha = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Elec / (mass_Elec*mass_Elec) + rho_Pos / (mass_Pos*mass_Pos)); 

	if (nDim == 2) 	{
		diff_ru_Pos_x  = 1.0/3.0*(ConsVar_Grad_0[5][0] + ConsVar_Grad_1[5][0] + ConsVar_Grad_2[5][0]);
		diff_ru_Elec_x = 1.0/3.0*(ConsVar_Grad_0[9][0] + ConsVar_Grad_1[9][0] + ConsVar_Grad_2[9][0]);
		diff_ru_Pos_y  = 1.0/3.0*(ConsVar_Grad_0[5][1] + ConsVar_Grad_1[5][1] + ConsVar_Grad_2[5][1]);
		diff_ru_Elec_y = 1.0/3.0*(ConsVar_Grad_0[9][1] + ConsVar_Grad_1[9][1] + ConsVar_Grad_2[9][1]);
	}

	if (nDim == 3) 	{	
		diff_ru_Pos_x  = 1.0/3.0*(ConsVar_Grad_0[6][0]	+ ConsVar_Grad_1[6][0] + ConsVar_Grad_2[6][0]);
		diff_ru_Elec_x = 1.0/3.0*(ConsVar_Grad_0[11][0]	+ ConsVar_Grad_1[11][0] + ConsVar_Grad_2[11][0]);
		diff_ru_Pos_y  = 1.0/3.0*(ConsVar_Grad_0[6][1]	+ ConsVar_Grad_1[6][1] + ConsVar_Grad_2[6][1]);
		diff_ru_Elec_y = 1.0/3.0*(ConsVar_Grad_0[11][1]	+ ConsVar_Grad_1[11][1] + ConsVar_Grad_2[11][1]);
		diff_ru_Pos_z  = 1.0/3.0*(ConsVar_Grad_0[6][2]  + ConsVar_Grad_1[6][2] + ConsVar_Grad_2[6][2]);
		diff_ru_Elec_z = 1.0/3.0*(ConsVar_Grad_0[11][2] + ConsVar_Grad_1[11][2] + ConsVar_Grad_2[11][2]);
	}

	if (nDim == 2) 			
		Kai_np1 = 1.0/(1.0+alpha) * (Kai_n - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Elec_x+diff_ru_Elec_y)/mass_Elec - (diff_ru_Pos_x+diff_ru_Pos_y)/mass_Pos) );
	if (nDim == 3) 			
		Kai_np1 = 1.0/(1.0+alpha) * (Kai_n - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Elec_x+diff_ru_Elec_y + diff_ru_Elec_z)/mass_Elec - (diff_ru_Pos_x+diff_ru_Pos_y+diff_ru_Pos_z)/mass_Pos) );


	/*--- Residual = transpose(N) * source (q) * Area  ---*/
	val_residual[0] =-1.0/3.0 * Kai_np1 * Area;
	val_residual[1] =-1.0/3.0 * Kai_np1 * Area;
	val_residual[2] =-1.0/3.0 * Kai_np1 * Area;	

}

CSourcePieceWise_AdjFlow::CSourcePieceWise_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iDim;

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity = new double [nVar];
	GradDensity = new double [nDim];
	GradInvDensity = new double [nDim];
	dPoDensity2 = new double [nDim];
	alpha = new double [nDim];
	beta = new double [nDim];
	Sigma_5_vec = new double [nDim];

	GradVel_o_Rho = new double* [nDim];
	sigma = new double* [nDim];
	Sigma_phi = new double* [nDim];
	Sigma_5_Tensor = new double* [nDim];
	Sigma = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++) {
		GradVel_o_Rho[iDim] = new double [nDim];
		sigma[iDim] = new double [nDim];
		Sigma_phi[iDim] = new double [nDim];
		Sigma_5_Tensor[iDim] = new double [nDim];
		Sigma[iDim] = new double [nDim];
	}
}

CSourcePieceWise_AdjFlow::~CSourcePieceWise_AdjFlow(void) {
	unsigned short iDim;

	for (iDim = 0; iDim < nDim; iDim++) {
		delete [] GradVel_o_Rho[iDim];
		delete [] sigma[iDim];
		delete [] Sigma_phi[iDim];
		delete [] Sigma_5_Tensor[iDim];
		delete [] Sigma[iDim];
	}
	delete [] GradVel_o_Rho;
	delete [] sigma;
	delete [] Sigma_phi;
	delete [] Sigma_5_Tensor;
	delete [] Sigma;

	delete [] Velocity;
	delete [] GradDensity;
	delete [] GradInvDensity;
	delete [] dPoDensity2;
	delete [] alpha;
	delete [] beta;
	delete [] Sigma_5_vec;
}

void CSourcePieceWise_AdjFlow::SetResidual (double *val_residual, CConfig *config) {

#ifdef DOMINO_VERSION

	unsigned short iDim;
	double rtau_xx, rtau_xy, rtau_xz, rtau_yy, rtau_yz, rtau_zz;
	double Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz, Sigma_xx5, Sigma_yy5, 
	Sigma_zz5, Sigma_xy5, Sigma_xz5, Sigma_yz5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz, 
	Sigma5_x, Sigma5_y, Sigma5_z, alpha_x, alpha_y, alpha_z, beta_x, beta_y, beta_z;
	double gm1 = Gamma_Minus_One;

	double Density = U_i[0];		
	double sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++){ 
		Velocity[iDim] = U_i[iDim+1]/Density;
		sq_vel += 0.5*Velocity[iDim]*Velocity[iDim]; 
	}
	double Energy = U_i[nDim+1]/Density;
	double SoundSpeed = sqrt(Gamma*gm1*(Energy-sq_vel));
	double P = (SoundSpeed*SoundSpeed*Density)/Gamma;
	double invDensity     = 1.0/Density;
	double invDensitysq   = invDensity*invDensity;
	double invDensitycube = invDensitysq*invDensity;
	double mue_eff = Laminar_Viscosity_i + Eddy_Viscosity_i;
	double lambda = -TWO3*mue_eff;
	double xi = Gamma*(Laminar_Viscosity_i/PRANDTL + Eddy_Viscosity_i/PRANDTL_TURB);
	double vv2 = sq_vel;

	/*--- Density gradient ---*/
	double dDensity_dx = PrimVar_Grad_i[nDim+2][0];
	double dDensity_dy = PrimVar_Grad_i[nDim+2][1];
	double dDensity_dz = 0.0;
	if (nDim == 3) dDensity_dz = PrimVar_Grad_i[nDim+2][2];

	/*--- Velocity field gradient ---*/
	double  GradVelocity[3][3];
	GradVelocity[0][0] = PrimVar_Grad_i[1][0];
	GradVelocity[0][1] = PrimVar_Grad_i[1][1];
	if (nDim == 3) GradVelocity[0][2] = PrimVar_Grad_i[1][2];
	GradVelocity[1][0] = PrimVar_Grad_i[2][0];
	GradVelocity[1][1] = PrimVar_Grad_i[2][1];
	if (nDim == 3) GradVelocity[1][2] = PrimVar_Grad_i[2][2];
	if (nDim == 3) {
		GradVelocity[2][0] = PrimVar_Grad_i[3][0];
		GradVelocity[2][1] = PrimVar_Grad_i[3][1];
		GradVelocity[2][2] = PrimVar_Grad_i[3][2];
	}

	/*--- Pressure gradient ---*/
	double dP_dx = PrimVar_Grad_i[nDim+1][0];
	double dP_dy = PrimVar_Grad_i[nDim+1][1];
	double dP_dz = 0.0;
	if (nDim == 3) dP_dz = PrimVar_Grad_i[nDim+1][2];

	/*--- Adjoint energy gradient ---*/
	double dpsi5_dx = PsiVar_Grad_i[nVar-1][0];
	double dpsi5_dy = PsiVar_Grad_i[nVar-1][1];
	double dpsi5_dz = 0.0;
	if (nDim == 3) dpsi5_dz = PsiVar_Grad_i[nVar-1][2];

	/*--- Adjoint velocity gradient ---*/
	double dpsi2_dx = PsiVar_Grad_i[1][0];
	double dpsi2_dy = PsiVar_Grad_i[1][1];
	double dpsi2_dz = 0.0;
	if (nDim == 3) dpsi2_dz = PsiVar_Grad_i[1][2];

	double dpsi3_dx = PsiVar_Grad_i[2][0];
	double dpsi3_dy = PsiVar_Grad_i[2][1];
	double dpsi3_dz = 0.0;
	if (nDim == 3) dpsi3_dz = PsiVar_Grad_i[2][2];

	double dpsi4_dx = 0.0;
	double dpsi4_dy = 0.0;
	double dpsi4_dz = 0.0;
	if (nDim == 3) {
		dpsi4_dx = PsiVar_Grad_i[3][0];
		dpsi4_dy = PsiVar_Grad_i[3][1];
		dpsi4_dz = PsiVar_Grad_i[3][2];
	}

	/*--- 1.- Computation of the gradient of (1/Density) = - gradient(Density)/(Density^2) ---*/
	double dinvDensity_dx = - dDensity_dx * invDensitysq;
	double dinvDensity_dy = - dDensity_dy * invDensitysq;
	double dinvDensity_dz = 0.0;
	if (nDim == 3) dinvDensity_dz = - dDensity_dz * invDensitysq;

	/*--- Computation of the  derivatives of P/(Density^2) It is computed as
	 dp_dk(P/(Density^2)) = (dP_dk - 2 * P * dDensity_dk) / Density^3 ---*/		 
	double dPDensity2_dx = (dP_dx * Density - 2.0 * dDensity_dx * P) * invDensitycube;
	double dPDensity2_dy = (dP_dy * Density - 2.0 * dDensity_dy * P) * invDensitycube;
	double dPDensity2_dz = 0.0;
	if (nDim == 3) dPDensity2_dz = (dP_dz * Density - 2.0 * dDensity_dz * P) * invDensitycube;

	/*--- Computation of the  derivatives of v_i/rho ---*/
	double dVelocityxDensity_dx = GradVelocity[0][0] / Density + Velocity[0] * dinvDensity_dx;
	double dVelocityxDensity_dy = GradVelocity[0][1] / Density + Velocity[0] * dinvDensity_dy;
	double dVelocityxDensity_dz = 0.0;
	if (nDim == 3) dVelocityxDensity_dz = GradVelocity[0][2] / Density + Velocity[0] * dinvDensity_dz;

	double dVelocityyDensity_dx = GradVelocity[1][0] / Density + Velocity[1] * dinvDensity_dx;
	double dVelocityyDensity_dy = GradVelocity[1][1] / Density + Velocity[1] * dinvDensity_dy;
	double dVelocityyDensity_dz = 0.0;
	if (nDim == 3) dVelocityyDensity_dz = GradVelocity[1][2] / Density + Velocity[1] * dinvDensity_dz;

	double dVelocityzDensity_dx = 0.0;
	double dVelocityzDensity_dy = 0.0;
	double dVelocityzDensity_dz = 0.0;
	if (nDim == 3) {
		dVelocityzDensity_dx = GradVelocity[2][0] / Density + Velocity[2] * dinvDensity_dx;
		dVelocityzDensity_dy = GradVelocity[2][1] / Density + Velocity[2] * dinvDensity_dy;
		dVelocityzDensity_dz = GradVelocity[2][2] / Density + Velocity[2] * dinvDensity_dz;
	}

	/*--- Abbreviations ---*/
	alpha_x = xi * dinvDensity_dx;
	alpha_y = xi * dinvDensity_dy;
	alpha_z = 0.0;
	if (nDim == 3) alpha_z = xi * dinvDensity_dz;

	beta_x = ( xi / gm1 ) * dPDensity2_dx;
	beta_y = ( xi / gm1 ) * dPDensity2_dy;
	beta_z = 0.0;
	if (nDim == 3) beta_z = ( xi / gm1 ) * dPDensity2_dz;

	/*--- Components of the effective stress tensor divided by density, point i ---*/
	if (nDim == 2) {

		rtau_xx = invDensity * lambda * (GradVelocity[1][1]  - 2.0 * GradVelocity[0][0]);
		rtau_yy = invDensity * lambda * (GradVelocity[0][0]  - 2.0 * GradVelocity[1][1]);

		rtau_xy = invDensity * mue_eff * (GradVelocity[0][1] + GradVelocity[1][0]);

		/*--- Components of the adjoint tensors ---*/

		Sigma_xx = mue_eff * (FOUR3 * dpsi2_dx -  TWO3 * dpsi3_dy);
		Sigma_yy = mue_eff * (-TWO3 * dpsi2_dx + FOUR3 * dpsi3_dy);

		Sigma_xy = mue_eff * (dpsi3_dx + dpsi2_dy);

		Sigma_xx5 = mue_eff * ( FOUR3 * Velocity[0] * dpsi5_dx - TWO3 * Velocity[1] * dpsi5_dy -  TWO3 * Velocity[2] * dpsi5_dz);
		Sigma_yy5 = mue_eff * (- TWO3 * Velocity[0] * dpsi5_dx + FOUR3 * Velocity[1] * dpsi5_dy -  TWO3 * Velocity[2] * dpsi5_dz);

		Sigma_xy5 = mue_eff * (Velocity[0] * dpsi5_dy + Velocity[1] * dpsi5_dx);
		Sigma_xz5 = mue_eff * (Velocity[0] * dpsi5_dz + Velocity[2] * dpsi5_dx);

		eta_xx = Sigma_xx + Sigma_xx5;
		eta_yy = Sigma_yy + Sigma_yy5;

		eta_xy = Sigma_xy + Sigma_xy5;

		Sigma5_x = ( xi /Density ) *  dpsi5_dx;
		Sigma5_y = ( xi /Density ) *  dpsi5_dy;

		/*--- Viscous residual ---*/

		val_residual[0] = (  (- rtau_xx * Velocity[0] - rtau_xy  * Velocity[1] 
		                                                                    + alpha_x * vv2  - beta_x) *  dpsi5_dx +
		                                                                    (- rtau_xy * Velocity[0] - rtau_yy  *  Velocity[1]
		                                                                                                                    + alpha_y * vv2  - beta_y) *  dpsi5_dy)*dV;

		val_residual[1] = (    (rtau_xx - Velocity[0] * alpha_x) * dpsi5_dx
				+ (rtau_xy - Velocity[0] * alpha_y) * dpsi5_dy)*dV;

		val_residual[2] = (    (rtau_xy - Velocity[1] * alpha_x) * dpsi5_dx
				+ (rtau_yy - Velocity[1] * alpha_y) * dpsi5_dy)*dV;

		val_residual[3]  = (alpha_x * dpsi5_dx + alpha_y * dpsi5_dy)*dV;

		/*--- Extra terms coming from grad(psi)* D *grad(invM) ---*/

		/*		val_residual[0] += (- dVelocityxDensity_dx * eta_xx - dVelocityxDensity_dy * eta_xy
							- dVelocityyDensity_dx * eta_xy - dVelocityyDensity_dy * eta_yy 
							+ Sigma5_x * Velocity[0] * GradVelocity[0][0] + Sigma5_x * Velocity[1] * GradVelocity[1][0]
							+ Sigma5_y * Velocity[0] * GradVelocity[0][1] + Sigma5_y * Velocity[1] * GradVelocity[1][1])*dV;

		val_residual[1] += (dinvDensity_dx * eta_xx + dinvDensity_dy * eta_xy 
							- GradVelocity[0][0] * Sigma5_x - GradVelocity[0][1] * Sigma5_y)*dV;

		val_residual[2] += (dinvDensity_dx * eta_xy + dinvDensity_dy * eta_yy 
							- GradVelocity[1][0] * Sigma5_x - GradVelocity[1][1] * Sigma5_y)*dV; ---*/

	}

	if (nDim == 3) {
		rtau_xx = invDensity * lambda * (GradVelocity[1][1] + GradVelocity[2][2] - 2.0 * GradVelocity[0][0]);
		rtau_yy = invDensity * lambda * (GradVelocity[0][0] + GradVelocity[2][2] - 2.0 * GradVelocity[1][1]);
		rtau_zz = invDensity * lambda * (GradVelocity[0][0] + GradVelocity[1][1] - 2.0 * GradVelocity[2][2]);

		rtau_xy = invDensity * mue_eff * (GradVelocity[0][1] + GradVelocity[1][0]);
		rtau_xz = invDensity * mue_eff * (GradVelocity[0][2] + GradVelocity[2][0]);
		rtau_yz = invDensity * mue_eff * (GradVelocity[1][2] + GradVelocity[2][1]);

		/*--- Components of the adjoint tensors ---*/

		Sigma_xx = mue_eff * (FOUR3 * dpsi2_dx -  TWO3 * dpsi3_dy 
				- TWO3  * dpsi4_dz);
		Sigma_yy = mue_eff * (-TWO3 * dpsi2_dx + FOUR3 * dpsi3_dy 
				- TWO3  * dpsi4_dz);
		Sigma_zz = mue_eff * (-TWO3 * dpsi2_dx -  TWO3 * dpsi3_dy 
				+ FOUR3 * dpsi4_dz);

		Sigma_xy = mue_eff * (dpsi3_dx + dpsi2_dy);
		Sigma_xz = mue_eff * (dpsi4_dx + dpsi2_dz);
		Sigma_yz = mue_eff * (dpsi4_dy + dpsi3_dz);

		Sigma_xx5 = mue_eff * ( FOUR3 * Velocity[0] * dpsi5_dx 
				- TWO3 * Velocity[1] * dpsi5_dy -  TWO3 * Velocity[2] * dpsi5_dz);
		Sigma_yy5 = mue_eff * (- TWO3 * Velocity[0] * dpsi5_dx 
				+ FOUR3 * Velocity[1] * dpsi5_dy -  TWO3 * Velocity[2] * dpsi5_dz);
		Sigma_zz5 = mue_eff * (- TWO3 * Velocity[0] * dpsi5_dx 
				- TWO3 * Velocity[1] * dpsi5_dy + FOUR3 * Velocity[2] * dpsi5_dz);

		Sigma_xy5 = mue_eff * (Velocity[0] * dpsi5_dy + Velocity[1] * dpsi5_dx);
		Sigma_xz5 = mue_eff * (Velocity[0] * dpsi5_dz + Velocity[2] * dpsi5_dx);
		Sigma_yz5 = mue_eff * (Velocity[1] * dpsi5_dz + Velocity[2] * dpsi5_dy);

		eta_xx = Sigma_xx + Sigma_xx5;
		eta_yy = Sigma_yy + Sigma_yy5;
		eta_zz = Sigma_zz + Sigma_zz5;

		eta_xy = Sigma_xy + Sigma_xy5;
		eta_xz = Sigma_xz + Sigma_xz5;
		eta_yz = Sigma_yz + Sigma_yz5;

		Sigma5_x = ( xi /Density ) *  dpsi5_dx;
		Sigma5_y = ( xi /Density ) *  dpsi5_dy;
		Sigma5_z = ( xi /Density ) *  dpsi5_dz;

		/*--- Viscous residual ---*/

		val_residual[0] = (  (- rtau_xx * Velocity[0] - rtau_xy  * Velocity[1] - rtau_xz * Velocity[2]
		                                                                                            + alpha_x * vv2  - beta_x) *  dpsi5_dx +
		                                                                                            (- rtau_xy * Velocity[0] - rtau_yy  *  Velocity[1] - rtau_yz * Velocity[2]
		                                                                                                                                                                    + alpha_y * vv2  - beta_y) *  dpsi5_dy +
		                                                                                                                                                                    (- rtau_xz * Velocity[0] - rtau_yz  *  Velocity[1] - rtau_zz * Velocity[2]
		                                                                                                                                                                                                                                            + alpha_z * vv2  - beta_z) *  dpsi5_dz)*dV;

		val_residual[1] = (    (rtau_xx - Velocity[0] * alpha_x) * dpsi5_dx
				+ (rtau_xy - Velocity[0] * alpha_y) * dpsi5_dy
				+ (rtau_xz - Velocity[0] * alpha_z) * dpsi5_dz)*dV;

		val_residual[2] = (    (rtau_xy - Velocity[1] * alpha_x) * dpsi5_dx
				+ (rtau_yy - Velocity[1] * alpha_y) * dpsi5_dy
				+ (rtau_yz - Velocity[1] * alpha_z) * dpsi5_dz)*dV;

		val_residual[3] = (	  (rtau_xz - Velocity[2] * alpha_x) * dpsi5_dx
				+ (rtau_yz - Velocity[2] * alpha_y) * dpsi5_dy
				+ (rtau_zz - Velocity[2] * alpha_z) * dpsi5_dz )*dV;

		val_residual[4]  = (alpha_x * dpsi5_dx + alpha_y * dpsi5_dy
				+ alpha_z * dpsi5_dz)*dV;

		/*--- Extra terms coming from grad(psi)* D *grad(invM) ---*/

		/*		val_residual[0] += (- dVelocityxDensity_dx * eta_xx - dVelocityxDensity_dy * eta_xy - dVelocityxDensity_dz * eta_xz
							- dVelocityyDensity_dx * eta_xy - dVelocityyDensity_dy * eta_yy - dVelocityyDensity_dz * eta_yz
							- dVelocityzDensity_dx * eta_xz - dVelocityzDensity_dy * eta_yz - dVelocityzDensity_dz * eta_zz
							+ Sigma5_x * Velocity[0] * GradVelocity[0][0] + Sigma5_x * Velocity[1] * GradVelocity[1][0] + Sigma5_x * Velocity[2] * GradVelocity[2][0]
							+ Sigma5_y * Velocity[0] * GradVelocity[0][1] + Sigma5_y * Velocity[1] * GradVelocity[1][1] + Sigma5_y * Velocity[2] * GradVelocity[2][1]
							+ Sigma5_z * Velocity[0] * GradVelocity[0][2] + Sigma5_z * Velocity[1] * GradVelocity[1][2] + Sigma5_z * Velocity[2] * GradVelocity[2][2])*dV;

		val_residual[1] += (dinvDensity_dx * eta_xx + dinvDensity_dy * eta_xy + dinvDensity_dz * eta_xz 
							- GradVelocity[0][0] * Sigma5_x - GradVelocity[0][1] * Sigma5_y - GradVelocity[0][2] * Sigma5_z)*dV;

		val_residual[2] += (dinvDensity_dx * eta_xy + dinvDensity_dy * eta_yy + dinvDensity_dz * eta_yz 
							- GradVelocity[1][0] * Sigma5_x - GradVelocity[1][1] * Sigma5_y - GradVelocity[1][2] * Sigma5_z)*dV;

		val_residual[3] += (dinvDensity_dx * eta_xz + dinvDensity_dy * eta_yz + dinvDensity_dz * eta_zz 
							- GradVelocity[2][0] * Sigma5_x - GradVelocity[2][1] * Sigma5_y - GradVelocity[2][2] * Sigma5_z)*dV; ---*/
	}

#else

	unsigned short iDim, jDim;
	double Density = U_i[0];		
	double sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) { 
		Velocity[iDim] = U_i[iDim+1]/Density;
		sq_vel += Velocity[iDim]*Velocity[iDim];
	}
	double Energy = U_i[nDim+1]/Density;
	double SoundSpeed = sqrt(Gamma*Gamma_Minus_One*(Energy-0.5*sq_vel));
	double Pressure = (SoundSpeed*SoundSpeed*Density)/Gamma;
	double invDensity     = 1.0/Density;
	double invDensitysq   = invDensity*invDensity;
	double invDensitycube = invDensitysq*invDensity;
	double mu_tot_1 = Laminar_Viscosity_i + Eddy_Viscosity_i;
	double mu_tot_2 = Laminar_Viscosity_i/PRANDTL + Eddy_Viscosity_i/PRANDTL_TURB;
	double invR = Gamma * config->GetMach_FreeStreamND() * config->GetMach_FreeStreamND();

	/*--- Required gradients of the flow variables, point j ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		/*--- grad density ---*/
		GradDensity[iDim] = PrimVar_Grad_i[nDim+2][iDim];
		/*--- grad (1/rho) ---*/
		GradInvDensity[iDim] = -GradDensity[iDim]*invDensitysq;
		/*--- Computation of the derivatives of P/(Density^2) ---*/						 
		dPoDensity2[iDim] = (PrimVar_Grad_i[nVar-1][iDim]*Density - 2.0*GradDensity[iDim]*Pressure)*invDensitycube;
		/*--- Abbreviations: alpha, beta, sigma_5_vec ---*/
		alpha[iDim] = Gamma*mu_tot_2*GradInvDensity[iDim];
		beta[iDim] = Gamma/Gamma_Minus_One*mu_tot_2*dPoDensity2[iDim];
		Sigma_5_vec[iDim] = Gamma*mu_tot_2*PsiVar_Grad_i[nVar-1][iDim];
	}

	/*--- Definition of tensors and derivatives of velocity over density ---*/
	double div_vel = 0, div_phi = 0, vel_gradpsi5 = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		div_vel += PrimVar_Grad_i[iDim+1][iDim];
		div_phi += PsiVar_Grad_i[iDim+1][iDim];
		vel_gradpsi5 += Velocity[iDim]*PsiVar_Grad_i[nVar-1][iDim];
		for (jDim = 0; jDim < nDim; jDim++) {
			sigma[iDim][jDim] = mu_tot_1*(PrimVar_Grad_i[iDim+1][jDim]+PrimVar_Grad_i[jDim+1][iDim]);
			Sigma_phi[iDim][jDim] = mu_tot_1*(PsiVar_Grad_i[iDim+1][jDim]+PsiVar_Grad_i[jDim+1][iDim]);
			Sigma_5_Tensor[iDim][jDim] = mu_tot_1*(Velocity[jDim]*PsiVar_Grad_i[nVar-1][iDim]+Velocity[iDim]*PsiVar_Grad_i[nVar-1][jDim]);
			GradVel_o_Rho[iDim][jDim] = (PrimVar_Grad_i[iDim+1][jDim]*Density - Velocity[iDim]*GradDensity[jDim])*invDensitysq;
		}
	}
	for (iDim = 0; iDim < nDim; iDim++) {
		sigma[iDim][iDim] -= TWO3*mu_tot_1*div_vel;
		Sigma_phi[iDim][iDim] -= TWO3*mu_tot_1*div_phi;
		Sigma_5_Tensor[iDim][iDim] -= TWO3*mu_tot_1*vel_gradpsi5;
	}
	for (iDim = 0; iDim < nDim; iDim++)
		for (jDim = 0; jDim < nDim; jDim++)
			Sigma[iDim][jDim] = Sigma_phi[iDim][jDim] + Sigma_5_Tensor[iDim][jDim];

	/*--- Vector-Tensors products ---*/
	double gradT_gradpsi5 = 0, sigma_gradpsi = 0, vel_sigma_gradpsi5 = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		gradT_gradpsi5 += PrimVar_Grad_i[0][iDim]*PsiVar_Grad_i[nVar-1][iDim];
		for (jDim = 0; jDim < nDim; jDim++) {
			sigma_gradpsi += sigma[iDim][jDim]*PsiVar_Grad_i[jDim+1][iDim];
			vel_sigma_gradpsi5 += Velocity[iDim]*sigma[iDim][jDim]*PsiVar_Grad_i[nVar-1][jDim];
		}
	}

	/*--- Residuals ---*/
	double alpha_gradpsi5 = 0, beta_gradpsi5 = 0, Sigma_gradvel_o_rho = 0, Sigma5_vel_gradvel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		alpha_gradpsi5 += alpha[iDim]*PsiVar_Grad_i[nVar-1][iDim];
		beta_gradpsi5 += beta[iDim]*PsiVar_Grad_i[nVar-1][iDim];
		for (jDim = 0; jDim < nDim; jDim++) {
			Sigma_gradvel_o_rho += Sigma[iDim][jDim]*GradVel_o_Rho[iDim][jDim];
			Sigma5_vel_gradvel += Sigma_5_vec[iDim]*(Velocity[jDim]*PrimVar_Grad_i[jDim+1][iDim]);
		}
	}
	val_residual[0] = (-vel_sigma_gradpsi5/Density - Sigma_gradvel_o_rho + 0.5*sq_vel*alpha_gradpsi5 - 
			beta_gradpsi5 + Sigma5_vel_gradvel/Density) * Volume;
	for (iDim = 0; iDim < nDim; iDim++)
		for (jDim = 0; jDim < nDim; jDim++)
			val_residual[iDim+1] = (sigma[iDim][jDim]*PsiVar_Grad_i[nVar-1][jDim]/Density + 
					Sigma[iDim][jDim]*GradInvDensity[jDim] - Velocity[iDim]*alpha_gradpsi5 -
					Sigma_5_vec[jDim]*PrimVar_Grad_i[iDim+1][jDim]/Density) * Volume;
	val_residual[nVar-1] = alpha_gradpsi5 * Volume;


	/*--- Coupling terms coming from the turbulent equations ---*/
	if (config->GetKind_Solver() == ADJ_RANS) {

		/*--- Closure constants ---*/
		double cv1_3 = 7.1*7.1*7.1;
		double k2 = 0.41*0.41;
		double cb1 = 0.1355;
		double cw2 = 0.3;
		double cw3_6 = pow(2.0,6.0);
		double sigma = 2./3.;
		double cb2 = 0.622;
		double cw1 = cb1/k2+(1+cb2)/sigma;

		double nu, Ji, Ji_2, Ji_3, fv1;
		nu = Laminar_Viscosity_i/U_i[0];
		Ji = TurbVar_i[0]/nu;
		Ji_2 = Ji*Ji;
		Ji_3 = Ji_2*Ji;
		fv1 = Ji_3/(Ji_3+cv1_3);

		/*--- Contributions due to variation of viscosities ---*/
		//		double mu1 = 1.404/config->GetReynolds();
		//		double mu2 = 0.404;
		//		double dVisc_T = Laminar_Viscosity_i*(Temp_i+3.0*mu2)/(2.0*Temp_i*(Temp_i+mu2));
		double dVisc_T = 0.0;

		double Cp = (Gamma/Gamma_Minus_One)/invR;
		double kappa_psi = (sigma_gradpsi + vel_sigma_gradpsi5)/mu_tot_1 + Cp/PRANDTL_TURB*gradT_gradpsi5;
		double cv1_const = 3.0*cv1_3/(Ji_3+cv1_3);
		double theta = (kappa_psi*(1.0-Eddy_Viscosity_i/Laminar_Viscosity_i*cv1_const) - 
				Cp/PRANDTL_TURB*gradT_gradpsi5*(1.0-PRANDTL_TURB/PRANDTL))*dVisc_T*Gamma_Minus_One*invR/Density;
		double xi = kappa_psi*(1+cv1_const)*Eddy_Viscosity_i/Density;

		val_residual[0] += (theta*(sq_vel-Energy) + xi)*Volume;
		for (iDim = 0; iDim < nDim; iDim++)
			val_residual[iDim+1] -= theta*Velocity[iDim]*Volume;
		val_residual[nVar-1] += theta*Volume;

		/*--- Coupling residuals ---*/
		if (dist_i > 0.0) {
			double fv2, Omega, Shat, dist_0_2, one_o_oneplusJifv1;
			double r, g, g_6, glim, fw;
			double dfw_g, dg_r, dr_nuhat, dr_Shat;
			double dShat_fv2, dfv2_fv1, dfv1_Ji, dJi_nu, dJi_nuhat, dfv2_Ji;

			/*--- Vorticity ---*/
			Omega = (PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0])*(PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0]);
			if (nDim == 3) Omega += (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) + 
					(PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1])*(PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1]);
			Omega = sqrt(Omega);

			dist_0_2 = dist_i*dist_i;
			one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
			fv2 = 1.0 - Ji*one_o_oneplusJifv1;
			Shat = max(Omega + TurbVar_i[0]*fv2/(k2*dist_0_2), TURB_EPS);

			// r = TurbVar_i[0]/(Shat*k2*dist_0_2);
			r = min(TurbVar_i[0]/(Shat*k2*dist_0_2), 10.);
			g = r + cw2*(pow(r,6.)-r);
			g_6 = pow(g,6.);
			glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
			fw = g*glim;

			dfw_g  = glim*cw3_6/(g_6+cw3_6); 
			dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0); 
			dr_nuhat = 1.0/(Shat*k2*dist_0_2); 
			dr_Shat = -dr_nuhat*TurbVar_i[0]/Shat; 

			dShat_fv2 = TurbVar_i[0]/(k2*dist_0_2); 
			dfv2_fv1 = Ji_2*one_o_oneplusJifv1*one_o_oneplusJifv1; 
			dfv1_Ji = 3.0*cv1_3*Ji_2/((Ji_3+cv1_3)*(Ji_3+cv1_3));
			dJi_nuhat = 1.0/nu; 
			dJi_nu = -Ji/nu; 
			dfv2_Ji = -one_o_oneplusJifv1*one_o_oneplusJifv1; 

			/*--- Terms 1 & 2: -Fcv·nabla(TurbPsi_i) - Fs·TurbPsi_i ---*/
			double gradTurbVar_gradTurbPsi = 0, vel_gradTurbPsi = 0;
			for (iDim = 0; iDim < nDim; iDim++) {
				gradTurbVar_gradTurbPsi += TurbVar_Grad_i[0][iDim]*TurbPsi_Grad_i[0][iDim];
				vel_gradTurbPsi += Velocity[iDim]*TurbPsi_Grad_i[0][iDim];
			}

			double alpha_coeff = Gamma_Minus_One*invR/Density*dVisc_T;
			double beta_coeff = alpha_coeff*(sq_vel-Energy)-Laminar_Viscosity_i/Density;
			double Fs_coeff = TurbPsi_i[0]*(cb1*TurbVar_i[0]-cw1*TurbVar_i[0]*TurbVar_i[0]/dist_0_2*dfw_g*dg_r*dr_Shat)*
					dShat_fv2*(dfv2_Ji+dfv2_fv1*dfv1_Ji)*dJi_nu;
			double Gamma = Fs_coeff - gradTurbVar_gradTurbPsi/sigma;

			val_residual[0] -= (Gamma*beta_coeff - TurbVar_i[0]*vel_gradTurbPsi)/Density*Volume;
			for (iDim = 0; iDim < nDim; iDim++)
				val_residual[iDim+1] += (Gamma*alpha_coeff*Velocity[iDim] - TurbVar_i[0]*TurbPsi_Grad_i[0][iDim])/Density*Volume;
			val_residual[nVar-1] -= (Gamma*alpha_coeff)/Density*Volume;

			/*--- Terms 3: -partial{T^s}_GradVel x GradN ---*/
			//			double Ms_coeff = (cb1*TurbVar_i[0]-cw1*TurbVar_i[0]*TurbVar_i[0]/dist_0_2*dfw_g*dg_r*dr_Shat);
			//			Ms_coeff *= TurbPsi_i[0]/(Omega + TURB_EPS);
			//
			//			for (iDim = 0; iDim < nDim; iDim++) {
			//				for (jDim = 0; jDim < nDim; jDim++) {
			//					val_residual[0] += Ms_coeff*(PrimVar_Grad_i[iDim+1][jDim]-PrimVar_Grad_i[jDim+1][iDim])*
			//					GradVel_o_Rho[iDim][jDim]*dV;
			//					val_residual[iDim+1] -= Ms_coeff*(PrimVar_Grad_i[iDim+1][jDim]-PrimVar_Grad_i[jDim+1][iDim])*
			//					GradInvDensity[jDim]*dV;
			//				}
			//			}
		}
	}

#endif

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

void CSourcePieceWise_AdjTurb::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
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

		if (implicit)
			val_Jacobian_i[0][0] = -Bs*Volume;

		/*---SECOND PART: \partial_nu_hat mu^k F^{vk} cdot \grad Psi ---*/
		double dEddyVisc_nuhat = U_i[0]*fv1*(1.0 + 3.0*cv1_3/(Ji_3+cv1_3));

		for (iDim = 0; iDim < nDim; iDim++) {
			for (jDim = 0; jDim < nDim; jDim++) 
				tau[iDim][jDim] = PrimVar_Grad_i[iDim+1][jDim] + PrimVar_Grad_i[jDim+1][iDim];
			tau[iDim][iDim] -= TWO3*div_vel;
		}

		double invR = Gamma*config->GetMach_FreeStreamND()*config->GetMach_FreeStreamND();
		double Cp = (Gamma/Gamma_Minus_One)/invR;
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

	}
}

CSourcePieceWise_AdjElec::CSourcePieceWise_AdjElec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { 

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

CSourcePieceWise_AdjElec::~CSourcePieceWise_AdjElec(void) { }

void CSourcePieceWise_AdjElec::SetResidual(double *val_residual, CConfig *config) {
	val_residual[0] = Volume*sin(PI_NUMBER*Coord_i[0])*sin(PI_NUMBER*Coord_i[1]);
}

CSourcePieceWise_LinElec::CSourcePieceWise_LinElec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

}

CSourcePieceWise_LinElec::~CSourcePieceWise_LinElec(void) { }

void CSourcePieceWise_LinElec::SetResidual(double *val_residual, CConfig *config) {
	val_residual[0] = Volume*Coord_i[0]*Coord_i[1]*
			pow((1.0-Coord_i[0]),config->GetChargeCoeff())*(1.0-Coord_i[1])*
			log(1.0-Coord_i[0]);
}

CSourceConservative_AdjFlow::CSourceConservative_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Velocity = new double [nDim];
	Residual_i = new double [nVar];
	Residual_j = new double [nVar];
	Mean_Residual = new double [nVar];

	Mean_PrimVar_Grad = new double* [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Mean_PrimVar_Grad[iVar] = new double [nDim];
}

CSourceConservative_AdjFlow::~CSourceConservative_AdjFlow(void) {
	delete [] Mean_Residual;
	delete [] Residual_j;
	delete [] Residual_i;
	delete [] Velocity;

	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_PrimVar_Grad[iVar];
	delete [] Mean_PrimVar_Grad;
}

void CSourceConservative_AdjFlow::SetResidual (double *val_residual, CConfig *config) {
	unsigned short iDim, jDim, iVar;
	double rho, nu, Ji, fv1, fv2, Omega, Shat, dist_sq, Ji_2, Ji_3, one_o_oneplusJifv1;
	double r, g, g_6, glim, dfw_g, dg_r, dr_nuhat, dr_Shat, Ms_coeff, invOmega;

	/*--- CLOUSURE CONSTANTS ---*/
	double cv1_3 = 7.1*7.1*7.1;
	double k2 = 0.41*0.41;
	double cb1 = 0.1355;
	double cw2 = 0.3;
	double cw3_6 = pow(2.0,6.0);
	double sigma = 2./3.;
	double cb2 = 0.622;
	double cw1 = cb1/k2+(1+cb2)/sigma;

	for (iVar = 0; iVar < nVar; iVar++) {
		Residual_i[iVar] = 0.0;
		Residual_j[iVar] = 0.0;
	}

	/*--- iPoint ---*/

	/*--- Density and velocities ---*/
	rho = U_i[0];
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity[iDim] = U_i[iDim+1]/rho;

	/*--- Vorticity ---*/
	Omega = (PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0])*(PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0]);
	if (nDim == 3) Omega += (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) + 
			(PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1])*(PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1]);
	Omega = sqrt(Omega);
	invOmega = 1.0/(Omega + TURB_EPS);
	//	invOmega = min(1.0/Omega, max_invOmega);

	/*--- Compute Ms_coeff -> coming from partial derivatives ---*/
	Ms_coeff = 0.0;
	if (dist_i > 0) {
		dist_sq = dist_i*dist_i;
		nu = Laminar_Viscosity_i/rho;
		Ji = TurbVar_i[0]/nu;
		Ji_2 = Ji*Ji;
		Ji_3 = Ji_2*Ji;
		fv1 = Ji_3/(Ji_3+cv1_3);
		one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
		fv2 = 1.0 - Ji*one_o_oneplusJifv1;
		Shat = max(Omega + TurbVar_i[0]*fv2/(k2*dist_sq),TURB_EPS);

		r = min(TurbVar_i[0]/(Shat*k2*dist_sq),10.);
		g = r + cw2*(pow(r,6.)-r);
		g_6 = pow(g,6.);
		glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);

		dfw_g  = glim*cw3_6/(g_6+cw3_6); 
		dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);  
		dr_nuhat = 1.0/(Shat*k2*dist_sq); 
		dr_Shat = -dr_nuhat*TurbVar_i[0]/Shat; 

		Ms_coeff = (cb1*TurbVar_i[0]-cw1*TurbVar_i[0]*TurbVar_i[0]/dist_sq*dfw_g*dg_r*dr_Shat);
	}
	Ms_coeff *= TurbPsi_i[0]*invOmega/rho;

	/*--- Compute residual of iPoint ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		for (jDim = 0; jDim < nDim; jDim++) {
			Residual_i[0] -= Ms_coeff*(Velocity[jDim]*PrimVar_Grad_i[jDim+1][iDim]*Normal[iDim] -
					Velocity[jDim]*PrimVar_Grad_i[iDim+1][jDim]*Normal[iDim]);
			Residual_i[iDim+1] += Ms_coeff*(PrimVar_Grad_i[iDim+1][jDim]*Normal[jDim] - 
					PrimVar_Grad_i[jDim+1][iDim]*Normal[jDim]);
		}
	}

	/*--- jPoint ---*/

	/*--- Density and velocities ---*/
	rho = U_j[0];
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity[iDim] = U_j[iDim+1]/rho;

	/*--- Vorticity ---*/
	Omega = (PrimVar_Grad_j[1][1]-PrimVar_Grad_j[2][0])*(PrimVar_Grad_j[1][1]-PrimVar_Grad_j[2][0]);
	if (nDim == 3) Omega += (PrimVar_Grad_j[1][2]-PrimVar_Grad_j[3][0])*(PrimVar_Grad_j[1][2]-PrimVar_Grad_j[3][0]) + 
			(PrimVar_Grad_j[2][2]-PrimVar_Grad_j[3][1])*(PrimVar_Grad_j[2][2]-PrimVar_Grad_j[3][1]);
	Omega = sqrt(Omega);
	invOmega = 1.0/(Omega + TURB_EPS);
	//	invOmega = min(1.0/Omega, max_invOmega);

	/*--- Compute Ms_coeff -> coming from partial derivatives ---*/
	Ms_coeff = 0.0;
	if (dist_j > 0) {
		dist_sq = dist_j*dist_j;
		nu = Laminar_Viscosity_j/rho;
		Ji = TurbVar_j[0]/nu;
		Ji_2 = Ji*Ji;
		Ji_3 = Ji_2*Ji;
		fv1 = Ji_3/(Ji_3+cv1_3);
		one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
		fv2 = 1.0 - Ji*one_o_oneplusJifv1;
		Shat = max(Omega + TurbVar_j[0]*fv2/(k2*dist_sq),TURB_EPS);

		r = min(TurbVar_j[0]/(Shat*k2*dist_sq),10.);
		g = r + cw2*(pow(r,6.)-r);
		g_6 = pow(g,6.);
		glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);

		dfw_g  = glim*cw3_6/(g_6+cw3_6); 
		dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);  
		dr_nuhat = 1.0/(Shat*k2*dist_sq); 
		dr_Shat = -dr_nuhat*TurbVar_j[0]/Shat; 

		Ms_coeff = (cb1*TurbVar_j[0]-cw1*TurbVar_j[0]*TurbVar_j[0]/dist_sq*dfw_g*dg_r*dr_Shat);
	}
	Ms_coeff *= TurbPsi_j[0]*invOmega/rho;

	/*--- Compute residual of jPoint ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		for (jDim = 0; jDim < nDim; jDim++) {
			Residual_j[0] -= Ms_coeff*(Velocity[jDim]*PrimVar_Grad_j[jDim+1][iDim]*Normal[iDim] -
					Velocity[jDim]*PrimVar_Grad_j[iDim+1][jDim]*Normal[iDim]);
			Residual_j[iDim+1] += Ms_coeff*(PrimVar_Grad_j[iDim+1][jDim]*Normal[jDim] - 
					PrimVar_Grad_j[jDim+1][iDim]*Normal[jDim]);
		}
	}

	/*--- MEAN RESIDUAL ---*/
	for (iVar = 0; iVar < nVar; iVar++) 
		val_residual[iVar] = 0.5*(Residual_i[iVar] + Residual_j[iVar]);

}

CSourceConservative_AdjTurb::CSourceConservative_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

}

CSourceConservative_AdjTurb::~CSourceConservative_AdjTurb(void) {
}

void CSourceConservative_AdjTurb::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	/*--- SOURCE term  -->  \nabla ( \psi_\mu · E^{s} )
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

CSourceRotationalFrame_Flow::CSourceRotationalFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

}

CSourceRotationalFrame_Flow::~CSourceRotationalFrame_Flow(void) { }

void CSourceRotationalFrame_Flow::SetResidual(double *val_residual, CConfig *config) {

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

CSourcePieceWise_Plasma::CSourcePieceWise_Plasma(unsigned short val_nDim, unsigned short val_nVar, 
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Kb = 1.38E-23; 
	M1 = config->GetParticle_Mass(0);
	M3 = config->GetParticle_Mass(2);
	M2 = M1-M3;

	AvgNum = AVOGAD_CONSTANT;
	M1Avg = M1*AVOGAD_CONSTANT;
	M2Avg = M2*AVOGAD_CONSTANT;
	M3Avg = M3*AVOGAD_CONSTANT;
	M1M2M3Avg3 = M1Avg * M2Avg* M3Avg;

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

}

void CSourcePieceWise_Plasma::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Ex = - ConsVar_Grad[0][0];
	Ey = - ConsVar_Grad[0][1];

	Ey = 0.0;
	if ( Ex < -10000)
		Ex = -10000.0;

	if (nDim == 3) 	Ex = 0.0;

	tol = 1E-60;
	double zero = 1E-15;


	/* brief: defining conservative variables, residual and jacobian for point i first: */
	if (nDim ==2 ) {
		r1 = U_i[0];	r2 = U_i[4];	r3 = U_i[8];	// density of each species
		m1 = U_i[1];	m2 = U_i[5];	m3 = U_i[9];	// density*horizontal velocity of each species
		n1 = U_i[2];	n2 = U_i[6];	n3 = U_i[10];	// density*vertical velocity of each species
		e1 = U_i[3];	e2 = U_i[7];	e3 = U_i[11];	// energy per unit volume of each species
	}

	if (nDim ==3 ) {
		r1 = U_i[0];	r2 = U_i[5];	r3 = U_i[10];	// density of each species
		m1 = U_i[1];	m2 = U_i[6];	m3 = U_i[11];	// density*horizontal velocity of each species
		n1 = U_i[2];	n2 = U_i[7];	n3 = U_i[12];	// density*vertical velocity of each species
		e1 = U_i[4];	e2 = U_i[9];	e3 = U_i[14];	// energy per unit volume of each species
	}

	P1 = (gam1-1)*(e1 - 0.5*(m1*m1 + n1*n1)/r1);	// Partial Pressure of species 1
	P2 = (gam2-1)*(e2 - 0.5*(m2*m2 + n2*n2)/r2);	// Partial Pressure of species 2
	P3 = (gam3-1)*(e3 - 0.5*(m3*m3 + n3*n3)/r3);	// Partial Pressure of species 3

	T1 = P1/(Rc/(M1*AvgNum)*r1);	// Temperature of species 1
	T2 = P2/(Rc/(M2*AvgNum)*r2);	// Temperature of species 2
	T3 = P3/(Rc/(M3*AvgNum)*r3);	// Temperature of species 3

	r1 = max(r1, zero);	r2 = max(r2, zero);	r3 = max(r3, zero);
	P1 = max(P1, zero);	P2 = max(P2, zero);	P3 = max(P3, zero);
	T1 = max(T1, zero);	T2 = max(T2, zero);	T3 = max(T3, zero);
	//m3 = max(m3,8000*1.2E-8);
	/*Brief: Partial derivative of species temperature with the conservative variables */

	//e1 = r1*Cv1*T1 + .5*(m1*m1 + n1*n1)/r1;

	dT1_dr1 = (-Cv1*T1 + .5*(m1*m1 + n1*n1)/(r1*r1))/(r1*Cv1);
	dT2_dr2 = (-Cv2*T2 + .5*(m2*m2 + n2*n2)/(r2*r2))/(r2*Cv2);
	dT3_dr3 = (-Cv3*T3 + .5*(m3*m3 + n3*n3)/(r3*r3))/(r3*Cv3);

	dT1_dm1 = - m1/r1 / (r1*Cv1) ;
	dT2_dm2 = - m2/r2 / (r2*Cv2);
	dT3_dm3 = - m3/r3 / (r3*Cv3);

	dT1_dn1 = - n1/r1 / (r1*Cv1) ;
	dT2_dn2 = - n2/r2 / (r2*Cv2);
	dT3_dn3 = - n3/r3 / (r3*Cv3);

	dT1_de1 = 1.0/(Cv1*r1);
	dT2_de2 = 1.0/(Cv2*r2);
	dT3_de3 = 1.0/(Cv3*r3);

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

	f13= 1*f13;
	f12= 1*f12;
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
	double FF =1.0;

	kf1 = FF*C1*pow(T1,eta1)*(theta1/T1 + 2.0)*exp(-theta1/T1);
	kf2 = FF*C2*pow(T2,eta2)*(theta2/T2 + 2.0)*exp(-theta2/T2);
	kf3 = FF*C3*pow(T3,eta3)*(theta3/T3 + 2.0)*exp(-theta3/T3);

	FF = 1.0;
	ke1 = FF*Ck1*pow(T1,zeta1)*exp(-phi1/T1);
	ke2 = FF*Ck2*pow(T2,zeta2)*exp(-phi2/T2);
	ke3 = FF*Ck3*pow(T3,zeta3)*exp(-phi3/T3);

	kb1 = kf1/ke1;
	kb2 = kf2/ke2;
	kb3 = kf3/ke3;

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
	dR_de1 =  ( -r1*r1/(M1*M1*AvgNum*AvgNum)*dkf1_dT1*dT1_de1   + r2*r3*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)*dkb1_dT1*dT1_de1);

	/* derivative of expression for rate of reaction R wrt to 5th,6th,7th and 8th conservative variables */
	dR_dr2 = ( + 0.0	   + kb1*r3*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)
			- kf2*r1/(M1*M2*AvgNum*AvgNum)	   + 2*kb2*r2*r3/(M2*M3*M2*AvgNum*AvgNum*AvgNum)
			+ 0.0	   + kb3*r3*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum));
	dR_dm2 = 0.0;
	dR_dn2 = 0.0;
	dR_de2 = ( - r1*r2/(M1*M2*AvgNum*AvgNum)*dkf2_dT2*dT2_de2  + r2*r3*r2/(M2*M3*M2*AvgNum*AvgNum*AvgNum)*dkb2_dT2*dT2_de2);


	/* derivative of expression for rate of reaction R wrt to last four conservative variables */
	dR_dr3 = ( + 0.0	 + kb1*r2*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)
			+ 0.0	 	 + kb2*r2*r2/(M2*M3*M2*AvgNum*AvgNum*AvgNum)
			- kf3*r1/(M1*M3*AvgNum*AvgNum)	 + 2*kb3*r2*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum));
	dR_dr3 = dR_dr2 * M3/M2;
	dR_dm3 = 0.0;
	dR_dn3 = 0.0;
	dR_de3 = ( - r1*r3/(M1*M3*AvgNum*AvgNum)*dkf3_dT3*dT3_de3  + r2*r3*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum)*dkb3_dT3*dT3_de3);
	dR_de3 = dR_de2*M3/M2;

	/* Derivative of the first term in the source terms wrt all the conservative variables */
	S1 = M1*AvgNum*R;

	dS1_dr1 = M1*AvgNum*dR_dr1;
	dS1_dm1 = M1*AvgNum*dR_dm1;
	dS1_dn1 = M1*AvgNum*dR_dn1;
	dS1_de1 = M1*AvgNum*dR_de1;

	/* Derivative of the second term in the source terms wrt all the conservative variables */

	S2 = r1*nu12*(m2/r2 - m1/r1)+ r1*nu13*(m3/r3 - m1/r1);

	dS2_dr1 = nu12*m2/r2 + nu13*m3/r3;
	dS2_dm1 = -nu12 - nu13;
	dS2_dn1 = 0.0;
	dS2_de1 = r1*dv12_dT1*dT1_de1*(m2/r2 - m1/r1) + r1*dv13_dT1*dT1_de1*(m3/r3 - m1/r1);

	/* Derivative of the third term in the source terms wrt all the conservative variables */

	S3 = r1*nu12*(n2/r2 - n1/r1) + r1*nu13*(n3/r3 - n1/r1);

	dS3_dr1 = nu12*n2/r2 + nu13*n3/r3;
	dS3_dm1 = 0.0;
	dS3_dn1 = -nu12 - nu13;
	dS3_de1 =  r1*(n2/r2 -n1/r1)*dv12_dT1*dT1_de1 + r1*(n3/r3 - n1/r1)*dv13_dT1*dT1_de1;

	/* Derivative of the Fourth term in the source terms wrt all the conservative variables */

	S4 = m1*nu12*(m2/r2 - m1/r1) + m1*nu13*(m3/r3 - m1/r1) +  n1*nu12*(n2/r2 - n1/r1) + n1*nu13*(n3/r3 - n1/r1) + QT1;

	dS4_dr1 = 0.0;//(nu12 + nu13)*(m1*m1 + n1*n1)/(r1*r1);// + dQT1_dT1*dT1_dr1;
	dS4_dm1 = nu12*(m2/r2 -2*m1/r1) + nu13*(m3/r3 - 2*m1/r1) + dQT1_dT1*dT1_dm1;
	dS4_dn1 = nu12*(n2/r2 -2*n1/r1) + nu13*(n3/r3 - 2*n1/r1) + dQT1_dT1*dT1_dn1;
	dS4_de1 = (m1*(m2/r2 - m1/r1) + n1*(n2/r2 - n1/r1))*dv12_dT1*dT1_de1 +
			(m1*(m3/r3 - m1/r1) + n1*(n3/r3 - n1/r1))*dv13_dT1*dT1_de1 + dQT1_dT1*dT1_de1;

	/* Derivative of the fifth term in the source terms wrt all the conservative variables */

	S5 = -M2*AvgNum*R;

	dS5_dr2 = -M2*AvgNum*dR_dr2;
	dS5_dm2 = -M2*AvgNum*dR_dm2;
	dS5_dn2 = -M2*AvgNum*dR_dn2;
	dS5_de2 = -M2*AvgNum*dR_de2;

	/* Derivative of the Sixth term in the source terms with all the conservative variables */

	S6 = r2*ec*Ex/M2 + r2*nu21*(m1/r1 - m2/r2) + r2*nu23*(m3/r3 - m2/r2);

	dS6_dr2 = nu21*m1/r1 + nu23*m3/r3 + ec*Ex/M2;
	dS6_dm2 = -nu21 - nu23;
	dS6_dn2 = 0;
	dS6_de2 = 0 + r2*(m1/r1 - m2/r2)*dv21_dT2*dT2_de2 + r2*(m3/r3 - m2/r2)*dv23_dT2*dT2_de2;

	/* Derivative of the Seventh term in the source terms wrt all the conservative variables */

	S7 = r2*ec*Ey/M2 + r2*nu21*(n1/r1 - n2/r2) + r2*nu23*(n3/r3 - n2/r2);

	dS7_dr2 = nu21*n1/r1 + nu23*n3/r3 + ec*Ey/M2;
	dS7_dm2 = 0;
	dS7_dn2 = -nu21 - nu23 ;
	dS7_de2 = 0 + r2*(n1/r1 - n2/r2)*dv21_dT2*dT2_de2 + r2*(n3/r3 - n2/r2)*dv23_dT2*dT2_de2;


	/* Derivative of the Eight term in the source terms wrt all the conservative variables */

	S8 = m2*ec*Ex/M2 + n2*ec*Ey/M2 + m2*nu21*(m1/r1 - m2/r2) + m2*nu23*(m3/r3 - m2/r2)+ n2*nu21*(n1/r1 - n2/r2) +
			n2*nu23*(n3/r3 - n2/r2) + QT2;

	dS8_dr2 = 0.0;//(nu21+nu23)*(m2*m2 + n2*n2)/(r2*r2) ;// +  dQT2_dT2*dT2_dr2;
	dS8_dm2 = nu21*(m1/r1 - 2*m2/r2) + nu23*(m3/r3 - 2*m2/r2) + ec*Ex/M2  + dQT2_dT2*dT2_dm2;
	dS8_dn2 = nu21*(n1/r1 - 2*n2/r2) + nu23*(n3/r3 - 2*n2/r2) + ec*Ey/M2  + dQT2_dT2*dT2_dn2;
	dS8_de2 = (n2*(n1/r1 - n2/r2)+ m2*(m1/r1 - m2/r2))*dv21_dT2*dT2_de2+
			(n2*(n3/r3 - n2/r2)+ m2*(m3/r3 - m2/r2))*dv23_dT2*dT2_de2 + dQT2_dT2*dT2_de2;


	/* Derivative of the ninth term in the source terms wrt all the conservative variables */

	S9 = -M3*AvgNum*R;

	dS9_dr3 = -M3*AvgNum*dR_dr3;
	dS9_dm3 = -M3*AvgNum*dR_dm3;
	dS9_dn3 = -M3*AvgNum*dR_dn3;
	dS9_de3 = -M3*AvgNum*dR_de3;

	/* Derivative of the Tenth term in the source terms wrt all the conservative variables */

	S10 =  -r3*ec*Ex/M3 + r3*nu31*(m1/r1 - m3/r3) + r3*nu32*(m2/r2 - m3/r3);

	dS10_dr3 = nu31*m1/r1 + nu32*m2/r2 - ec*Ex/M3 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_dr3 +  r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_dr3;
	dS10_dm3 =  -nu31 - nu32 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_dm3  + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_dm3;
	dS10_dn3 =  0 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_dn3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_dn3;
	dS10_de3 =  0 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_de3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_de3;
	//dS10_dm3 = 10*dS10_dm3;

	/* Derivative of the Eleventh term in the source terms wrt all the conservative variables */

	S11 = -r3*ec*Ey/M3 + r3*nu31*(n1/r1 - n3/r3) + r3*nu32*(n2/r2 - n3/r3);

	dS11_dr3 = nu31*n1/r1 + nu32*n2/r2 - ec*Ey/M3 ;//+ r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_dr3 +  r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_dr3;
	dS11_dm3 =  0 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_dm3  + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_dm3;
	dS11_dn3 =  -nu31 - nu32 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_dn3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_dn3;
	dS11_de3 =  0 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_de3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_de3;

	/* Derivative of the Twelfth term in the source terms wrt all the conservative variables */

	S12 = -m3*ec*Ex/M3 - n3*ec*Ey/M3 + m3*nu31*(m1/r1 - m3/r3) + m3*nu32*(m2/r2 -m3/r3) +
			n3*nu31*(n1/r1 - n3/r3) + n3*nu32*(n2/r2 - n3/r3) + QT3;

	dS12_dr3 = (nu31+nu32)*(m3*m3 + n3*n3)/(r3*r3);// + dQT3_dT3*dT3_dr3;
	dS12_dm3 = nu31*(m1/r1 - 2*m3/r3) + nu32*(m2/r2 - 2*m3/r3) -ec*Ex/M3 + dQT3_dT3*dT3_dm3;
	dS12_dn3 = nu31*(n1/r1 - 2*n3/r3) + nu32*(n2/r2 - 2*n3/r3) -ec*Ey/M3 + dQT3_dT3*dT3_dn3;
	dS12_de3 = (m3*(m1/r1 - m3/r3)+  n3*(n1/r1 - n3/r3))*dv31_dT3*dT3_de3 +
			(m3*(m2/r2 - m3/r3)+  n3*(n2/r2 - n3/r3))*dv32_dT3*dT3_de3 + dQT3_dT3*dT3_de3;


	// The factor is sqmueo
	/*	double sqmueo = 1.12E-3;
		cout << "T1 = " << T1 << " T2 = " << T2 << " T3 = " << T3 << endl;
		cout << " UUS(1) = " << m1/r1 << " VVS(1) =" << n1/r1 << endl;
		cout << " UUS(2) = " << m2/r2 << " VVS(2) =" << n2/r2 << endl;
		cout << " UUS(3) = " << m3/r3 << " VVS(3) =" << n3/r3 << endl;
		cout << " CHSCON(1) = " << 0.0 << " CHSCON(2) = " <<r2*ec/M2*1.12E-3 <<  " CHSCON(3) = " <<-r3*ec/M3*1.12E-3 << endl;
		cout << " EFX = " << Ex/sqmueo << " EFY = " << Ey/sqmueo << endl;
		cout << " Momentum Transfer VIA collision source terms" << endl;
		cout << " FMTX(1) = " << S2 << " FMTY(1) = " << S3 << endl;
		cout << " FMTX(2) = " << S6 << " FMTY(2) = " << S7 << endl;
		cout << " FMTX(3) = " << S10 << " FMTY(3) = " << S11 << endl;
		cout << " FMTIA(1) = " << r1*nu12 << " FMTEA(1) = " << r1*nu13  << " FMTEI(1) =" << r2*nu23<< endl;
		cout << " RNUEA = " << f13 << " RNUEI = " << f23 <<  endl;
		cout << " RELAX(1) = " << QT1 << " RELAX(2) = " << QT2 << " RELAX(3) = " << QT3 << endl;
	 */
	if (nDim ==2) {
		val_residual[0] = S1*Volume;
		val_residual[1] = S2*Volume;
		val_residual[2] = S3*Volume;
		val_residual[3] = S4*Volume;
		val_residual[4] = S5*Volume;
		val_residual[5] = S6*Volume;
		val_residual[6] = S7*Volume;
		val_residual[7] = S8*Volume;
		val_residual[8] = S9*Volume;
		val_residual[9] = S10*Volume;
		val_residual[10] = S11*Volume;
		val_residual[11] = S12*Volume;
	}

	if (nDim ==3) {
		val_residual[0] = S1*Volume;
		val_residual[1] = S2*Volume;
		val_residual[2] = S3*Volume;
		val_residual[3] = 0;
		val_residual[4] = S4*Volume;

		val_residual[5] = S5*Volume;
		val_residual[6] = S6*Volume;
		val_residual[7] = S7*Volume;
		val_residual[8] = 0.0;
		val_residual[9] = S8*Volume;

		val_residual[10] = S9*Volume;
		val_residual[11] = S10*Volume;
		val_residual[12] = S11*Volume;
		val_residual[13] = 0.0;
		val_residual[14] = S12*Volume;

	}

	for (unsigned short iVar = 0; iVar < nVar; iVar ++)
		if(fabs(val_residual[iVar]) < tol)
			val_residual[iVar] = 0.0;


	for (unsigned short iVar = 0; iVar < nVar; iVar ++)
		for (unsigned short jVar = 0; jVar < nVar; jVar ++)
			val_Jacobian_i[iVar][jVar] = 0.0;


	if (implicit) {

		if (nDim ==2) {
			val_Jacobian_i[0][0] = dS1_dr1*Volume;
			val_Jacobian_i[0][1] = dS1_dm1*Volume;
			val_Jacobian_i[0][2] = dS1_dn1*Volume;
			val_Jacobian_i[0][3] = dS1_de1*Volume;

			val_Jacobian_i[1][0] = dS2_dr1*Volume;
			val_Jacobian_i[1][1] = dS2_dm1*Volume;
			val_Jacobian_i[1][2] = dS2_dn1*Volume;
			val_Jacobian_i[1][3] = dS2_de1*Volume;

			val_Jacobian_i[2][0] = dS3_dr1*Volume;
			val_Jacobian_i[2][1] = dS3_dm1*Volume;
			val_Jacobian_i[2][2] = dS3_dn1*Volume;
			val_Jacobian_i[2][3] = dS3_de1*Volume;

			val_Jacobian_i[3][0] = dS4_dr1*Volume;
			val_Jacobian_i[3][1] = dS4_dm1*Volume;
			val_Jacobian_i[3][2] = dS4_dn1*Volume;
			val_Jacobian_i[3][3] = dS4_de1*Volume;

			val_Jacobian_i[4][4] = dS5_dr2*Volume;
			val_Jacobian_i[4][5] = dS5_dm2*Volume;
			val_Jacobian_i[4][6] = dS5_dn2*Volume;
			val_Jacobian_i[4][7] = dS5_de2*Volume;

			val_Jacobian_i[5][4] = dS6_dr2*Volume;
			val_Jacobian_i[5][5] = dS6_dm2*Volume;
			val_Jacobian_i[5][6] = dS6_dn2*Volume;
			val_Jacobian_i[5][7] = dS6_de2*Volume;

			val_Jacobian_i[6][4] = dS7_dr2*Volume;
			val_Jacobian_i[6][5] = dS7_dm2*Volume;
			val_Jacobian_i[6][6] = dS7_dn2*Volume;
			val_Jacobian_i[6][7] = dS7_de2*Volume;

			val_Jacobian_i[7][4] = dS8_dr2*Volume;
			val_Jacobian_i[7][5] = dS8_dm2*Volume;
			val_Jacobian_i[7][6] = dS8_dn2*Volume;
			val_Jacobian_i[7][7] = dS8_de2*Volume;

			val_Jacobian_i[8][8]   = dS9_dr3*Volume;
			val_Jacobian_i[8][9]   = dS9_dm3*Volume;
			val_Jacobian_i[8][10]  = dS9_dn3*Volume;
			val_Jacobian_i[8][11]  = dS9_de3*Volume;

			val_Jacobian_i[9][8]   = dS10_dr3*Volume;
			val_Jacobian_i[9][9]   = dS10_dm3*Volume;
			val_Jacobian_i[9][10]  = dS10_dn3*Volume;
			val_Jacobian_i[9][11]  = dS10_de3*Volume;

			val_Jacobian_i[10][8]  = dS11_dr3*Volume;
			val_Jacobian_i[10][9]  = dS11_dm3*Volume;
			val_Jacobian_i[10][10] = dS11_dn3*Volume;
			val_Jacobian_i[10][11] = dS11_de3*Volume;

			val_Jacobian_i[11][8]  = dS12_dr3*Volume;
			val_Jacobian_i[11][9]  = dS12_dm3*Volume;
			val_Jacobian_i[11][10] = dS12_dn3*Volume;
			val_Jacobian_i[11][11] = dS12_de3*Volume;

		}

		if (nDim ==3) {
			val_Jacobian_i[0][0] = dS1_dr1*Volume;
			val_Jacobian_i[0][1] = dS1_dm1*Volume;
			val_Jacobian_i[0][2] = dS1_dn1*Volume;
			val_Jacobian_i[0][3] = 0;
			val_Jacobian_i[0][4] = dS1_de1*Volume;

			val_Jacobian_i[1][0] = dS2_dr1*Volume;
			val_Jacobian_i[1][1] = dS2_dm1*Volume;
			val_Jacobian_i[1][2] = dS2_dn1*Volume;
			val_Jacobian_i[1][3] = 0;
			val_Jacobian_i[1][4] = dS2_de1*Volume;

			val_Jacobian_i[2][0] = dS3_dr1*Volume;
			val_Jacobian_i[2][1] = dS3_dm1*Volume;
			val_Jacobian_i[2][2] = dS3_dn1*Volume;
			val_Jacobian_i[2][3] = 0;
			val_Jacobian_i[2][4] = dS3_de1*Volume;

			val_Jacobian_i[3][0] = 0*Volume;
			val_Jacobian_i[3][1] = 0*Volume;
			val_Jacobian_i[3][2] = 0*Volume;
			val_Jacobian_i[3][3] = 0;
			val_Jacobian_i[3][4] = 0*Volume;

			val_Jacobian_i[4][0] = dS4_dr1*Volume;
			val_Jacobian_i[4][1] = dS4_dm1*Volume;
			val_Jacobian_i[4][2] = dS4_dn1*Volume;
			val_Jacobian_i[4][3] = 0;
			val_Jacobian_i[4][4] = dS4_de1*Volume;




			val_Jacobian_i[5][5] = dS5_dr2*Volume;
			val_Jacobian_i[5][6] = dS5_dm2*Volume;
			val_Jacobian_i[5][7] = dS5_dn2*Volume;
			val_Jacobian_i[5][8] = 0;
			val_Jacobian_i[5][9] = dS5_de2*Volume;

			val_Jacobian_i[6][5] = dS6_dr2*Volume;
			val_Jacobian_i[6][6] = dS6_dm2*Volume;
			val_Jacobian_i[6][7] = dS6_dn2*Volume;
			val_Jacobian_i[6][8] = 0*Volume;
			val_Jacobian_i[6][9] = dS6_de2*Volume;

			val_Jacobian_i[7][5] = dS7_dr2*Volume;
			val_Jacobian_i[7][6] = dS7_dm2*Volume;
			val_Jacobian_i[7][7] = dS7_dn2*Volume;
			val_Jacobian_i[7][8] = 0*Volume;
			val_Jacobian_i[7][9] = dS7_de2*Volume;

			val_Jacobian_i[8][5] = 0*Volume;
			val_Jacobian_i[8][6] = 0*Volume;
			val_Jacobian_i[8][7] = 0*Volume;
			val_Jacobian_i[8][8] = 0*Volume;
			val_Jacobian_i[8][9] = 0*Volume;

			val_Jacobian_i[9][5] = dS8_dr2*Volume;
			val_Jacobian_i[9][6] = dS8_dm2*Volume;
			val_Jacobian_i[9][7] = dS8_dn2*Volume;
			val_Jacobian_i[9][8] = 0*Volume;
			val_Jacobian_i[9][9] = dS8_de2*Volume;



			val_Jacobian_i[10][10] = dS9_dr3*Volume;
			val_Jacobian_i[10][11] = dS9_dm3*Volume;
			val_Jacobian_i[10][12] = dS9_dn3*Volume;
			val_Jacobian_i[10][13] = 0*Volume;
			val_Jacobian_i[10][14] = dS9_de3*Volume;

			val_Jacobian_i[11][10] = dS10_dr3*Volume;
			val_Jacobian_i[11][11] = dS10_dm3*Volume;
			val_Jacobian_i[11][12] = dS10_dn3*Volume;
			val_Jacobian_i[11][13] = 0*Volume;
			val_Jacobian_i[11][14] = dS10_de3*Volume;

			val_Jacobian_i[12][10] = dS11_dr3*Volume;
			val_Jacobian_i[12][11] = dS11_dm3*Volume;
			val_Jacobian_i[12][12] = dS11_dn3*Volume;
			val_Jacobian_i[12][13] = 0*Volume;
			val_Jacobian_i[12][14] = dS11_de3*Volume;

			val_Jacobian_i[13][10] = 0*Volume;
			val_Jacobian_i[13][11] = 0*Volume;
			val_Jacobian_i[13][12] = 0*Volume;
			val_Jacobian_i[13][13] = 0*Volume;
			val_Jacobian_i[13][14] = 0*Volume;

			val_Jacobian_i[14][10] = dS12_dr3*Volume;
			val_Jacobian_i[14][11] = dS12_dm3*Volume;
			val_Jacobian_i[14][12] = dS12_dn3*Volume;
			val_Jacobian_i[14][13] = 0*Volume;
			val_Jacobian_i[14][14] = dS12_de3*Volume;

		}




		for (unsigned short iVar = 0; iVar < nVar; iVar ++)
			for (unsigned short jVar = 0; jVar < nVar; jVar ++) {
				if(fabs(val_Jacobian_i[iVar][jVar]) < tol)
					val_Jacobian_i[iVar][jVar] = 0.0;
			}
	}
}

CSourcePieceWise_PlasmaDiatomic::CSourcePieceWise_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar, 
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Kb = 1.38E-23; 
	M1 = config->GetParticle_Mass(0);
	M3 = config->GetParticle_Mass(2);
	M2 = M1-M3;

	AvgNum = AVOGAD_CONSTANT;
	M1Avg = M1*AVOGAD_CONSTANT;
	M2Avg = M2*AVOGAD_CONSTANT;
	M3Avg = M3*AVOGAD_CONSTANT;
	M1M2M3Avg3 = M1Avg * M2Avg* M3Avg;

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

CSourcePieceWise_PlasmaDiatomic::~CSourcePieceWise_PlasmaDiatomic(void) {

}

void CSourcePieceWise_PlasmaDiatomic::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Ex = - ConsVar_Grad[0][0];
	Ey = - ConsVar_Grad[0][1];

	Ey = 0.0;
	if ( Ex < -10000)
		Ex = -10000.0;

	tol = 1E-60;
	double zero = 1E-15;

	/* brief: defining conservative variables, residual and jacobian for point i first: */
	r1 = U_i[0];	r2 = U_i[4];	r3 = U_i[8];	// density of each species
	m1 = U_i[1];	m2 = U_i[5];	m3 = U_i[9];	// density*horizontal velocity of each species
	n1 = U_i[2];	n2 = U_i[6];	n3 = U_i[10];	// density*vertical velocity of each species
	e1 = U_i[3];	e2 = U_i[7];	e3 = U_i[11];	// energy per unit volume of each species

	P1 = (gam1-1)*(e1 - 0.5*(m1*m1 + n1*n1)/r1);	// Partial Pressure of species 1
	P2 = (gam2-1)*(e2 - 0.5*(m2*m2 + n2*n2)/r2);	// Partial Pressure of species 2
	P3 = (gam3-1)*(e3 - 0.5*(m3*m3 + n3*n3)/r3);	// Partial Pressure of species 3

	T1 = P1/(Rc/(M1*AvgNum)*r1);	// Temperature of species 1
	T2 = P2/(Rc/(M2*AvgNum)*r2);	// Temperature of species 2
	T3 = P3/(Rc/(M3*AvgNum)*r3);	// Temperature of species 3

	r1 = max(r1, zero);	r2 = max(r2, zero);	r3 = max(r3, zero);
	P1 = max(P1, zero);	P2 = max(P2, zero);	P3 = max(P3, zero);
	T1 = max(T1, zero);	T2 = max(T2, zero);	T3 = max(T3, zero);
	//m3 = max(m3,8000*1.2E-8);
	/*Brief: Partial derivative of species temperature with the conservative variables */

	//e1 = r1*Cv1*T1 + .5*(m1*m1 + n1*n1)/r1;

	dT1_dr1 = (-Cv1*T1 + .5*(m1*m1 + n1*n1)/(r1*r1))/(r1*Cv1);
	dT2_dr2 = (-Cv2*T2 + .5*(m2*m2 + n2*n2)/(r2*r2))/(r2*Cv2);
	dT3_dr3 = (-Cv3*T3 + .5*(m3*m3 + n3*n3)/(r3*r3))/(r3*Cv3);

	dT1_dm1 = - m1/r1 / (r1*Cv1) ;
	dT2_dm2 = - m2/r2 / (r2*Cv2);
	dT3_dm3 = - m3/r3 / (r3*Cv3);

	dT1_dn1 = - n1/r1 / (r1*Cv1) ;
	dT2_dn2 = - n2/r2 / (r2*Cv2);
	dT3_dn3 = - n3/r3 / (r3*Cv3);

	dT1_de1 = 1.0/(Cv1*r1);
	dT2_de2 = 1.0/(Cv2*r2);
	dT3_de3 = 1.0/(Cv3*r3);

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

	f13= 1*f13;
	f12= 1*f12;
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
	double FF =1.0;

	kf1 = FF*C1*pow(T1,eta1)*(theta1/T1 + 2.0)*exp(-theta1/T1);
	kf2 = FF*C2*pow(T2,eta2)*(theta2/T2 + 2.0)*exp(-theta2/T2);
	kf3 = FF*C3*pow(T3,eta3)*(theta3/T3 + 2.0)*exp(-theta3/T3);

	FF = 1.0;
	ke1 = FF*Ck1*pow(T1,zeta1)*exp(-phi1/T1);
	ke2 = FF*Ck2*pow(T2,zeta2)*exp(-phi2/T2);
	ke3 = FF*Ck3*pow(T3,zeta3)*exp(-phi3/T3);

	kb1 = kf1/ke1;
	kb2 = kf2/ke2;
	kb3 = kf3/ke3;

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
	dR_de1 =  ( -r1*r1/(M1*M1*AvgNum*AvgNum)*dkf1_dT1*dT1_de1   + r2*r3*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)*dkb1_dT1*dT1_de1);

	/* derivative of expression for rate of reaction R wrt to 5th,6th,7th and 8th conservative variables */
	dR_dr2 = ( + 0.0	   + kb1*r3*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)
			- kf2*r1/(M1*M2*AvgNum*AvgNum)	   + 2*kb2*r2*r3/(M2*M3*M2*AvgNum*AvgNum*AvgNum)
			+ 0.0	   + kb3*r3*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum));
	dR_dm2 = 0.0;
	dR_dn2 = 0.0;
	dR_de2 = ( - r1*r2/(M1*M2*AvgNum*AvgNum)*dkf2_dT2*dT2_de2  + r2*r3*r2/(M2*M3*M2*AvgNum*AvgNum*AvgNum)*dkb2_dT2*dT2_de2);


	/* derivative of expression for rate of reaction R wrt to last four conservative variables */
	dR_dr3 = ( + 0.0	 + kb1*r2*r1/(M2*M3*M1*AvgNum*AvgNum*AvgNum)
			+ 0.0	 	 + kb2*r2*r2/(M2*M3*M2*AvgNum*AvgNum*AvgNum)
			- kf3*r1/(M1*M3*AvgNum*AvgNum)	 + 2*kb3*r2*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum));
	dR_dr3 = dR_dr2 * M3/M2;
	dR_dm3 = 0.0;
	dR_dn3 = 0.0;
	dR_de3 = ( - r1*r3/(M1*M3*AvgNum*AvgNum)*dkf3_dT3*dT3_de3  + r2*r3*r3/(M2*M3*M3*AvgNum*AvgNum*AvgNum)*dkb3_dT3*dT3_de3);
	dR_de3 = dR_de2*M3/M2;

	/* Derivative of the first term in the source terms wrt all the conservative variables */
	S1 = M1*AvgNum*R;

	dS1_dr1 = M1*AvgNum*dR_dr1;
	dS1_dm1 = M1*AvgNum*dR_dm1;
	dS1_dn1 = M1*AvgNum*dR_dn1;
	dS1_de1 = M1*AvgNum*dR_de1;

	/* Derivative of the second term in the source terms wrt all the conservative variables */

	S2 = r1*nu12*(m2/r2 - m1/r1)+ r1*nu13*(m3/r3 - m1/r1);

	dS2_dr1 = nu12*m2/r2 + nu13*m3/r3;
	dS2_dm1 = -nu12 - nu13;
	dS2_dn1 = 0.0;
	dS2_de1 = r1*dv12_dT1*dT1_de1*(m2/r2 - m1/r1) + r1*dv13_dT1*dT1_de1*(m3/r3 - m1/r1);

	/* Derivative of the third term in the source terms wrt all the conservative variables */

	S3 = r1*nu12*(n2/r2 - n1/r1) + r1*nu13*(n3/r3 - n1/r1);

	dS3_dr1 = nu12*n2/r2 + nu13*n3/r3;
	dS3_dm1 = 0.0;
	dS3_dn1 = -nu12 - nu13;
	dS3_de1 =  r1*(n2/r2 -n1/r1)*dv12_dT1*dT1_de1 + r1*(n3/r3 - n1/r1)*dv13_dT1*dT1_de1;

	/* Derivative of the Fourth term in the source terms wrt all the conservative variables */

	S4 = m1*nu12*(m2/r2 - m1/r1) + m1*nu13*(m3/r3 - m1/r1) +  n1*nu12*(n2/r2 - n1/r1) + n1*nu13*(n3/r3 - n1/r1) + QT1;

	dS4_dr1 = 0.0;//(nu12 + nu13)*(m1*m1 + n1*n1)/(r1*r1);// + dQT1_dT1*dT1_dr1;
	dS4_dm1 = nu12*(m2/r2 -2*m1/r1) + nu13*(m3/r3 - 2*m1/r1) + dQT1_dT1*dT1_dm1;
	dS4_dn1 = nu12*(n2/r2 -2*n1/r1) + nu13*(n3/r3 - 2*n1/r1) + dQT1_dT1*dT1_dn1;
	dS4_de1 = (m1*(m2/r2 - m1/r1) + n1*(n2/r2 - n1/r1))*dv12_dT1*dT1_de1 +
			(m1*(m3/r3 - m1/r1) + n1*(n3/r3 - n1/r1))*dv13_dT1*dT1_de1 + dQT1_dT1*dT1_de1;

	/* Derivative of the fifth term in the source terms wrt all the conservative variables */

	S5 = -M2*AvgNum*R;

	dS5_dr2 = -M2*AvgNum*dR_dr2;
	dS5_dm2 = -M2*AvgNum*dR_dm2;
	dS5_dn2 = -M2*AvgNum*dR_dn2;
	dS5_de2 = -M2*AvgNum*dR_de2;

	/* Derivative of the Sixth term in the source terms with all the conservative variables */

	S6 = r2*ec*Ex/M2 + r2*nu21*(m1/r1 - m2/r2) + r2*nu23*(m3/r3 - m2/r2);

	dS6_dr2 = nu21*m1/r1 + nu23*m3/r3 + ec*Ex/M2;
	dS6_dm2 = -nu21 - nu23;
	dS6_dn2 = 0;
	dS6_de2 = 0 + r2*(m1/r1 - m2/r2)*dv21_dT2*dT2_de2 + r2*(m3/r3 - m2/r2)*dv23_dT2*dT2_de2;

	/* Derivative of the Seventh term in the source terms wrt all the conservative variables */

	S7 = r2*ec*Ey/M2 + r2*nu21*(n1/r1 - n2/r2) + r2*nu23*(n3/r3 - n2/r2);

	dS7_dr2 = nu21*n1/r1 + nu23*n3/r3 + ec*Ey/M2;
	dS7_dm2 = 0;
	dS7_dn2 = -nu21 - nu23 ;
	dS7_de2 = 0 + r2*(n1/r1 - n2/r2)*dv21_dT2*dT2_de2 + r2*(n3/r3 - n2/r2)*dv23_dT2*dT2_de2;


	/* Derivative of the Eight term in the source terms wrt all the conservative variables */

	S8 = m2*ec*Ex/M2 + n2*ec*Ey/M2 + m2*nu21*(m1/r1 - m2/r2) + m2*nu23*(m3/r3 - m2/r2)+ n2*nu21*(n1/r1 - n2/r2) +
			n2*nu23*(n3/r3 - n2/r2) + QT2;

	dS8_dr2 = 0.0;//(nu21+nu23)*(m2*m2 + n2*n2)/(r2*r2) ;// +  dQT2_dT2*dT2_dr2;
	dS8_dm2 = nu21*(m1/r1 - 2*m2/r2) + nu23*(m3/r3 - 2*m2/r2) + ec*Ex/M2  + dQT2_dT2*dT2_dm2;
	dS8_dn2 = nu21*(n1/r1 - 2*n2/r2) + nu23*(n3/r3 - 2*n2/r2) + ec*Ey/M2  + dQT2_dT2*dT2_dn2;
	dS8_de2 = (n2*(n1/r1 - n2/r2)+ m2*(m1/r1 - m2/r2))*dv21_dT2*dT2_de2+
			(n2*(n3/r3 - n2/r2)+ m2*(m3/r3 - m2/r2))*dv23_dT2*dT2_de2 + dQT2_dT2*dT2_de2;


	/* Derivative of the ninth term in the source terms wrt all the conservative variables */

	S9 = -M3*AvgNum*R;

	dS9_dr3 = -M3*AvgNum*dR_dr3;
	dS9_dm3 = -M3*AvgNum*dR_dm3;
	dS9_dn3 = -M3*AvgNum*dR_dn3;
	dS9_de3 = -M3*AvgNum*dR_de3;

	/* Derivative of the Tenth term in the source terms wrt all the conservative variables */

	S10 =  -r3*ec*Ex/M3 + r3*nu31*(m1/r1 - m3/r3) + r3*nu32*(m2/r2 - m3/r3);

	dS10_dr3 = nu31*m1/r1 + nu32*m2/r2 - ec*Ex/M3 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_dr3 +  r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_dr3;
	dS10_dm3 =  -nu31 - nu32 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_dm3  + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_dm3;
	dS10_dn3 =  0 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_dn3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_dn3;
	dS10_de3 =  0 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_de3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_de3;
	//dS10_dm3 = 10*dS10_dm3;

	/* Derivative of the Eleventh term in the source terms wrt all the conservative variables */

	S11 = -r3*ec*Ey/M3 + r3*nu31*(n1/r1 - n3/r3) + r3*nu32*(n2/r2 - n3/r3);

	dS11_dr3 = nu31*n1/r1 + nu32*n2/r2 - ec*Ey/M3 ;//+ r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_dr3 +  r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_dr3;
	dS11_dm3 =  0 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_dm3  + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_dm3;
	dS11_dn3 =  -nu31 - nu32 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_dn3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_dn3;
	dS11_de3 =  0 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_de3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_de3;

	/* Derivative of the Twelfth term in the source terms wrt all the conservative variables */

	S12 = -m3*ec*Ex/M3 - n3*ec*Ey/M3 + m3*nu31*(m1/r1 - m3/r3) + m3*nu32*(m2/r2 -m3/r3) +
			n3*nu31*(n1/r1 - n3/r3) + n3*nu32*(n2/r2 - n3/r3) + QT3;

	dS12_dr3 = (nu31+nu32)*(m3*m3 + n3*n3)/(r3*r3);// + dQT3_dT3*dT3_dr3;
	dS12_dm3 = nu31*(m1/r1 - 2*m3/r3) + nu32*(m2/r2 - 2*m3/r3) -ec*Ex/M3 + dQT3_dT3*dT3_dm3;
	dS12_dn3 = nu31*(n1/r1 - 2*n3/r3) + nu32*(n2/r2 - 2*n3/r3) -ec*Ey/M3 + dQT3_dT3*dT3_dn3;
	dS12_de3 = (m3*(m1/r1 - m3/r3)+  n3*(n1/r1 - n3/r3))*dv31_dT3*dT3_de3 +
			(m3*(m2/r2 - m3/r3)+  n3*(n2/r2 - n3/r3))*dv32_dT3*dT3_de3 + dQT3_dT3*dT3_de3;


	// The factor is sqmueo
	/*	double sqmueo = 1.12E-3;
	 cout << "T1 = " << T1 << " T2 = " << T2 << " T3 = " << T3 << endl;
	 cout << " UUS(1) = " << m1/r1 << " VVS(1) =" << n1/r1 << endl;
	 cout << " UUS(2) = " << m2/r2 << " VVS(2) =" << n2/r2 << endl;
	 cout << " UUS(3) = " << m3/r3 << " VVS(3) =" << n3/r3 << endl;
	 cout << " CHSCON(1) = " << 0.0 << " CHSCON(2) = " <<r2*ec/M2*1.12E-3 <<  " CHSCON(3) = " <<-r3*ec/M3*1.12E-3 << endl;
	 cout << " EFX = " << Ex/sqmueo << " EFY = " << Ey/sqmueo << endl;
	 cout << " Momentum Transfer VIA collision source terms" << endl;
	 cout << " FMTX(1) = " << S2 << " FMTY(1) = " << S3 << endl;
	 cout << " FMTX(2) = " << S6 << " FMTY(2) = " << S7 << endl;
	 cout << " FMTX(3) = " << S10 << " FMTY(3) = " << S11 << endl;
	 cout << " FMTIA(1) = " << r1*nu12 << " FMTEA(1) = " << r1*nu13  << " FMTEI(1) =" << r2*nu23<< endl;
	 cout << " RNUEA = " << f13 << " RNUEI = " << f23 <<  endl;
	 cout << " RELAX(1) = " << QT1 << " RELAX(2) = " << QT2 << " RELAX(3) = " << QT3 << endl;
	 */
	val_residual[0] = S1*Volume;
	val_residual[1] = S2*Volume;
	val_residual[2] = S3*Volume;
	val_residual[3] = S4*Volume;
	val_residual[4] = S5*Volume;
	val_residual[5] = S6*Volume;
	val_residual[6] = S7*Volume;
	val_residual[7] = S8*Volume;
	val_residual[8] = S9*Volume;
	val_residual[9] = S10*Volume;
	val_residual[10] = S11*Volume;
	val_residual[11] = S12*Volume;

	for (unsigned short iVar = 0; iVar < 12; iVar ++)
		if(fabs(val_residual[iVar]) < tol)
			val_residual[iVar] = 0.0;

	if (implicit) {

		val_Jacobian_i[0][0] = dS1_dr1*Volume;
		val_Jacobian_i[0][1] = dS1_dm1*Volume;
		val_Jacobian_i[0][2] = dS1_dn1*Volume;
		val_Jacobian_i[0][3] = dS1_de1*Volume;
		val_Jacobian_i[0][4] = 0;
		val_Jacobian_i[0][5] = 0;
		val_Jacobian_i[0][6] = 0;
		val_Jacobian_i[0][7] = 0;
		val_Jacobian_i[0][8] = 0;
		val_Jacobian_i[0][9] = 0;
		val_Jacobian_i[0][10] = 0;
		val_Jacobian_i[0][11] = 0;

		val_Jacobian_i[1][0] = dS2_dr1*Volume;
		val_Jacobian_i[1][1] = dS2_dm1*Volume;
		val_Jacobian_i[1][2] = dS2_dn1*Volume;
		val_Jacobian_i[1][3] = dS2_de1*Volume;
		val_Jacobian_i[1][4] = 0.0;
		val_Jacobian_i[1][5] = 0.0;
		val_Jacobian_i[1][6] = 0.0;
		val_Jacobian_i[1][7] = 0.0;
		val_Jacobian_i[1][8] = 0.0;
		val_Jacobian_i[1][9] = 0.0;
		val_Jacobian_i[1][10] = 0.0;
		val_Jacobian_i[1][11] = 0.0;

		val_Jacobian_i[2][0] = dS3_dr1*Volume;
		val_Jacobian_i[2][1] = dS3_dm1*Volume;
		val_Jacobian_i[2][2] = dS3_dn1*Volume;
		val_Jacobian_i[2][3] = dS3_de1*Volume;
		val_Jacobian_i[2][4] = 0.0;
		val_Jacobian_i[2][5] = 0.0;
		val_Jacobian_i[2][6] = 0.0;
		val_Jacobian_i[2][7] = 0.0;
		val_Jacobian_i[2][8] = 0.0;
		val_Jacobian_i[2][9] = 0.0;
		val_Jacobian_i[2][10] = 0.0;
		val_Jacobian_i[2][11] = 0.0;


		val_Jacobian_i[3][0] = dS4_dr1*Volume;
		val_Jacobian_i[3][1] = dS4_dm1*Volume;
		val_Jacobian_i[3][2] = dS4_dn1*Volume;
		val_Jacobian_i[3][3] = dS4_de1*Volume;
		val_Jacobian_i[3][4] = 0.0;
		val_Jacobian_i[3][5] = 0.0;
		val_Jacobian_i[3][6] = 0.0;
		val_Jacobian_i[3][7] = 0.0;
		val_Jacobian_i[3][8] = 0.0;
		val_Jacobian_i[3][9] = 0.0;
		val_Jacobian_i[3][10] = 0.0;
		val_Jacobian_i[3][11] = 0.0;

		val_Jacobian_i[4][0] = 0.0;
		val_Jacobian_i[4][1] = 0.0;
		val_Jacobian_i[4][2] = 0.0;
		val_Jacobian_i[4][3] = 0.0;
		val_Jacobian_i[4][4] = dS5_dr2*Volume;
		val_Jacobian_i[4][5] = dS5_dm2*Volume;
		val_Jacobian_i[4][6] = dS5_dn2*Volume;
		val_Jacobian_i[4][7] = dS5_de2*Volume;
		val_Jacobian_i[4][8] = 0.0;
		val_Jacobian_i[4][9] = 0.0;
		val_Jacobian_i[4][10] = 0.0;
		val_Jacobian_i[4][11] = 0.0;

		val_Jacobian_i[5][0] = 0.0;
		val_Jacobian_i[5][1] = 0.0;
		val_Jacobian_i[5][2] = 0.0;
		val_Jacobian_i[5][3] = 0.0;
		val_Jacobian_i[5][4] = dS6_dr2*Volume;
		val_Jacobian_i[5][5] = dS6_dm2*Volume;
		val_Jacobian_i[5][6] = dS6_dn2*Volume;
		val_Jacobian_i[5][7] = dS6_de2*Volume;
		val_Jacobian_i[5][8] = 0.0;
		val_Jacobian_i[5][9] = 0.0;
		val_Jacobian_i[5][10] = 0.0;
		val_Jacobian_i[5][11] = 0.0;

		val_Jacobian_i[6][0] = 0.0;
		val_Jacobian_i[6][1] = 0.0;
		val_Jacobian_i[6][2] = 0.0;
		val_Jacobian_i[6][3] = 0.0;
		val_Jacobian_i[6][4] = dS7_dr2*Volume;
		val_Jacobian_i[6][5] = dS7_dm2*Volume;
		val_Jacobian_i[6][6] = dS7_dn2*Volume;
		val_Jacobian_i[6][7] = dS7_de2*Volume;
		val_Jacobian_i[6][8] = 0.0;
		val_Jacobian_i[6][9] = 0.0;
		val_Jacobian_i[6][10] = 0.0;
		val_Jacobian_i[6][11] = 0.0;

		val_Jacobian_i[7][0] = 0.0;
		val_Jacobian_i[7][1] = 0.0;
		val_Jacobian_i[7][2] = 0.0;
		val_Jacobian_i[7][3] = 0.0;
		val_Jacobian_i[7][4] = dS8_dr2*Volume;
		val_Jacobian_i[7][5] = dS8_dm2*Volume;
		val_Jacobian_i[7][6] = dS8_dn2*Volume;
		val_Jacobian_i[7][7] = dS8_de2*Volume;
		val_Jacobian_i[7][8] = 0.0;
		val_Jacobian_i[7][9] = 0.0;
		val_Jacobian_i[7][10] = 0.0;
		val_Jacobian_i[7][11] = 0.0;

		val_Jacobian_i[8][0] = 0.0;
		val_Jacobian_i[8][1] = 0.0;
		val_Jacobian_i[8][2] = 0.0;
		val_Jacobian_i[8][3] = 0.0;
		val_Jacobian_i[8][4] = 0.0;
		val_Jacobian_i[8][5] = 0.0;
		val_Jacobian_i[8][6] = 0.0;
		val_Jacobian_i[8][7] = 0.0;
		val_Jacobian_i[8][8] = dS9_dr3*Volume;
		val_Jacobian_i[8][9] = dS9_dm3*Volume;
		val_Jacobian_i[8][10] = dS9_dn3*Volume;
		val_Jacobian_i[8][11] = dS9_de3*Volume;

		val_Jacobian_i[9][0] = 0.0;
		val_Jacobian_i[9][1] = 0.0;
		val_Jacobian_i[9][2] = 0.0;
		val_Jacobian_i[9][3] = 0.0;
		val_Jacobian_i[9][4] = 0.0;
		val_Jacobian_i[9][5] = 0.0;
		val_Jacobian_i[9][6] = 0.0;
		val_Jacobian_i[9][7] = 0.0;
		val_Jacobian_i[9][8] = dS10_dr3*Volume;
		val_Jacobian_i[9][9] = dS10_dm3*Volume;
		val_Jacobian_i[9][10] = dS10_dn3*Volume;
		val_Jacobian_i[9][11] = dS10_de3*Volume;

		val_Jacobian_i[10][0] = 0.0;
		val_Jacobian_i[10][1] = 0.0;
		val_Jacobian_i[10][2] = 0.0;
		val_Jacobian_i[10][3] = 0.0;
		val_Jacobian_i[10][4] = 0.0;
		val_Jacobian_i[10][5] = 0.0;
		val_Jacobian_i[10][6] = 0.0;
		val_Jacobian_i[10][7] = 0.0;
		val_Jacobian_i[10][8] = dS11_dr3*Volume;
		val_Jacobian_i[10][9] = dS11_dm3*Volume;
		val_Jacobian_i[10][10] = dS11_dn3*Volume;
		val_Jacobian_i[10][11] = dS11_de3*Volume;

		val_Jacobian_i[11][0] = 0.0;
		val_Jacobian_i[11][1] = 0.0;
		val_Jacobian_i[11][2] = 0.0;
		val_Jacobian_i[11][3] = 0.0;
		val_Jacobian_i[11][4] = 0.0;
		val_Jacobian_i[11][5] = 0.0;
		val_Jacobian_i[11][6] = 0.0;
		val_Jacobian_i[11][7] = 0.0;
		val_Jacobian_i[11][8] = dS12_dr3*Volume;
		val_Jacobian_i[11][9] = dS12_dm3*Volume;
		val_Jacobian_i[11][10] = dS12_dn3*Volume;
		val_Jacobian_i[11][11] = dS12_de3*Volume;

		/*for (unsigned short iVar = 0; iVar < 12; iVar ++)
		 for (unsigned short jVar = 0; jVar < 12; jVar ++)
		 if(fabs(val_Jacobian_i[iVar][jVar]) < tol)
		 val_Jacobian_i[iVar][jVar] = 0.0;
		 */
	}
}

CSourcePieceWise_Combustion::CSourcePieceWise_Combustion(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	Gas_Constant = config->GetGas_Constant();

	Alpha = 1000;
	T_i = 288.1;

}

CSourcePieceWise_Combustion::~CSourcePieceWise_Combustion(void) { }

void CSourcePieceWise_Combustion::SetResidual(double *val_residual, CConfig *config) {

	val_residual[0] = 0.0;
	double rho = U_i[0];
	double rhoE = U_i[nDim+1];

	double Velocity2 = 0.0; 
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Velocity2 += U_i[iDim+1]*U_i[iDim+1]/(rho*rho);
	double Pressure = (Gamma-1.0)*rho*(rhoE/rho-0.5*Velocity2);
	double temp = Pressure / ( Gas_Constant * rho );

	double lambda = LambdaComb_i;

	/*--- Evaluate the dource term  ---*/
	if (temp > T_i) {
		val_residual[0] = Alpha*rho*(1-lambda)*Volume;
	}

}

CSource_Template::CSource_Template(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSource_Template::~CSource_Template(void) { }

void CSource_Template::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) { }
