/*!
 * \file numerics_source.cpp
 * \brief This file contains all the source term discretization.
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

}

CSourcePieceWise_TurbSST::~CSourcePieceWise_TurbSST(void) { }

void CSourcePieceWise_TurbSST::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	unsigned short iDim;
	double gamma_blended, beta_blended;
	double diverg, production;

	val_residual[0] = 0.0;
	val_residual[1] = 0.0;
	if (implicit){
		val_Jacobian_i[0][0] = 0.0;		val_Jacobian_i[0][1] = 0.0;
		val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = 0.0;
	}

	/*--- Computation of blended constants for the source terms---*/
	gamma_blended = F1_i*gamma_1+(1-F1_i)*gamma_2;
	beta_blended  = F1_i*beta_1+(1-F1_i)*beta_2;

	if (dist_i > 0.0) {
		/*--- Production ---*/
		diverg = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			diverg += PrimVar_Grad_i[iDim+1][iDim];

		production = Eddy_Viscosity_i*StrainMag*StrainMag - 2.0/3.0*U_i[0]*TurbVar_i[0]*diverg;
		production = min(production,20*beta_star*U_i[0]*TurbVar_i[1]*TurbVar_i[0]);
		production = max(production,0.0);

		val_residual[0] += production*Volume;
		val_residual[1] += gamma_blended/(Eddy_Viscosity_i/U_i[0])*production*Volume;

		/*--- Dissipation ---*/
		val_residual[0] -= beta_star*U_i[0]*TurbVar_i[1]*TurbVar_i[0]*Volume;
		val_residual[1] -= beta_blended*U_i[0]*TurbVar_i[1]*TurbVar_i[1]*Volume;

		/*--- Cross diffusion ---*/
		val_residual[1] += (1.0-F1_i)*CDkw*Volume;

		/*--- Implicit part ---*/
		if (implicit) {
			val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;		val_Jacobian_i[0][1] = 0.0;
			val_Jacobian_i[1][0] = 0.0;									val_Jacobian_i[1][1] = -2*beta_blended*TurbVar_i[1]*Volume;
		}
	}

}

CSourcePieceWise_FreeSurface::CSourcePieceWise_FreeSurface(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

}

CSourcePieceWise_FreeSurface::~CSourcePieceWise_FreeSurface(void) { }

void CSourcePieceWise_FreeSurface::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	unsigned short iVar;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = 0.0;

	if (nDim == 2) {
		val_residual[1] = -Volume * U_i[0] * AuxVar_Grad_i[0];
		val_residual[2] = -Volume * U_i[0] * AuxVar_Grad_i[1];

		if (implicit) {
			val_Jacobian_i[0][0] = 0.0;													val_Jacobian_i[0][1] = 0.0;		val_Jacobian_i[0][2] = 0.0;
			val_Jacobian_i[1][0] = -Volume*AuxVar_Grad_i[0];		val_Jacobian_i[1][1] = 0.0;		val_Jacobian_i[1][2] = 0.0;
			val_Jacobian_i[2][0] = -Volume*AuxVar_Grad_i[1];		val_Jacobian_i[2][1] = 0.0;		val_Jacobian_i[2][2] = 0.0;
		}
	} 


	/*	if (nDim == 2) {
		val_residual[0] = 0.0; //(Volume/DensityInc_i) * ( ( U_i[1] ) * AuxVar_Grad_i[0] + ( U_i[2] ) * AuxVar_Grad_i[1] );
		val_residual[1] = (Volume/DensityInc_i) * ( ( U_i[1]*U_i[1] + U_i[0]/DensityInc_i ) * AuxVar_Grad_i[0] + (U_i[1]*U_i[2]) * AuxVar_Grad_i[1] );
		val_residual[2] = (Volume/DensityInc_i) * ( (U_i[1]*U_i[2]) * AuxVar_Grad_i[0] + ( U_i[2]*U_i[2] + U_i[0]/DensityInc_i ) * AuxVar_Grad_i[1] );

		if (implicit) {
			val_Jacobian_i[0][0] = 0.0;		
			val_Jacobian_i[0][1] = 0.0; //(Volume/DensityInc_i) * AuxVar_Grad_i[0];		
			val_Jacobian_i[0][2] = 0.0; //(Volume/DensityInc_i) * AuxVar_Grad_i[1];
			val_Jacobian_i[1][0] = (Volume/DensityInc_i) * (1.0/DensityInc_i) * AuxVar_Grad_i[0];		
			val_Jacobian_i[1][1] = (Volume/DensityInc_i) * (2.0*U_i[1]*AuxVar_Grad_i[0]+U_i[2]*AuxVar_Grad_i[1]);		
			val_Jacobian_i[1][2] = (Volume/DensityInc_i) * (U_i[1]*AuxVar_Grad_i[1]);	
			val_Jacobian_i[2][0] = (Volume/DensityInc_i) * (1.0/DensityInc_i) * AuxVar_Grad_i[1];		
			val_Jacobian_i[2][1] = (Volume/DensityInc_i) * (U_i[2]*AuxVar_Grad_i[0]);	
			val_Jacobian_i[2][2] = (Volume/DensityInc_i) * (U_i[1]*AuxVar_Grad_i[0]+2.0*U_i[2]*AuxVar_Grad_i[1]);	
		}

	} */

}

CSourcePieceWise_Gravity::CSourcePieceWise_Gravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	incompressible = (config->GetIncompressible() == YES);

}

CSourcePieceWise_Gravity::~CSourcePieceWise_Gravity(void) { }

void CSourcePieceWise_Gravity::SetResidual(double *val_residual, CConfig *config) {
	unsigned short iVar;

	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = 0.0;

	if (incompressible) {

		/*--- Compute the Froude number  ---*/
		Froude = config->GetFroude();

		/*--- Evaluate the source term  ---*/
		val_residual[nDim] = Volume * DensityInc_i / (Froude * Froude);

	}

}

CSourcePieceWise_Elec::CSourcePieceWise_Elec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

CSourcePieceWise_Elec::~CSourcePieceWise_Elec(void) { }

void CSourcePieceWise_Elec::SetResidual(double *val_residual, CConfig *config) {

	if (config->GetKind_GasModel() == ARGON) {
		double Kai_n, Kai_np1 = 0.0, ne, np;
		double diff_ru_Elec_x, diff_ru_Pos_x, diff_ru_Elec_y, diff_ru_Pos_y;
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

			rho_Pos  = 1.0/3.0*(U_0[4] + U_1[4] + U_2[4]) ; 
			rho_Elec = 1.0/3.0*(U_0[8] + U_1[8] + U_2[8]) ; 

			/*--- Source q ---*/
			mass_Elec = config->GetParticle_Mass(2);
			mass_Pos  = config->GetParticle_Mass(1);

			ne = rho_Elec / mass_Elec;
			np = rho_Pos  / mass_Pos;

			Kai_n = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne - np );
			alpha = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Elec / (mass_Elec*mass_Elec) + rho_Pos / (mass_Pos*mass_Pos));

			diff_ru_Pos_x  = 1.0/3.0*(ConsVar_Grad_0[5][0] + ConsVar_Grad_1[5][0] + ConsVar_Grad_2[5][0]);
			diff_ru_Elec_x = 1.0/3.0*(ConsVar_Grad_0[9][0] + ConsVar_Grad_1[9][0] + ConsVar_Grad_2[9][0]);
			diff_ru_Pos_y  = 1.0/3.0*(ConsVar_Grad_0[5][1] + ConsVar_Grad_1[5][1] + ConsVar_Grad_2[5][1]);
			diff_ru_Elec_y = 1.0/3.0*(ConsVar_Grad_0[9][1] + ConsVar_Grad_1[9][1] + ConsVar_Grad_2[9][1]);


			Kai_np1 = 1.0/(1.0+alpha) * (Kai_n - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Elec_x+diff_ru_Elec_y)/mass_Elec - (diff_ru_Pos_x+diff_ru_Pos_y)/mass_Pos) );

			/*--- Residual = transpose(N) * source (q) * Area  ---*/
			val_residual[0] =-1.0/3.0 * Kai_np1 * Area;
			val_residual[1] =-1.0/3.0 * Kai_np1 * Area;
			val_residual[2] =-1.0/3.0 * Kai_np1 * Area;

			/*		val_residual[0] = 0.0;
			val_residual[1] = 0.0;
			val_residual[2] = 0.0;
			 */

		}


	}
	if (config->GetKind_GasModel() == AIR7) {
		double Chi_n_0, Chi_n_1, Chi_n_2; 
		double Chi_np1_0 = 0.0, Chi_np1_1 = 0.0, Chi_np1_2 = 0.0;
		double ne_0, ne_1, ne_2;
		double np_0, np_1, np_2;
		double diff_ru_Neg_x_0 = 0.0, diff_ru_Neg_x_1 = 0.0, diff_ru_Neg_x_2 = 0.0;
		double diff_ru_Pos_x_0 = 0.0, diff_ru_Pos_x_1 = 0.0, diff_ru_Pos_x_2 = 0.0; 
		double diff_ru_Neg_y_0 = 0.0, diff_ru_Neg_y_1 = 0.0, diff_ru_Neg_y_2 = 0.0; 
		double diff_ru_Pos_y_0 = 0.0, diff_ru_Pos_y_1 = 0.0, diff_ru_Pos_y_2 = 0.0;
		double alpha_0, alpha_1, alpha_2;
		double dt = TimeStep;
		double rho_Pos_0 = 0.0, rho_Pos_1 = 0.0, rho_Pos_2 = 0.0;
		double rho_Neg_0 = 0.0, rho_Neg_1 = 0.0, rho_Neg_2 = 0.0;		
		double mass_Pos_0 = 0.0, mass_Pos_1 = 0.0, mass_Pos_2 = 0.0;
		double massTotal_Pos_0 = 0.0, massTotal_Pos_1 = 0.0, massTotal_Pos_2 = 0.0;
		double mass_Neg_0 = 0.0, mass_Neg_1 = 0.0, mass_Neg_2 = 0.0;
		double massTotal_Neg_0 = 0.0, massTotal_Neg_1 = 0.0, massTotal_Neg_2 = 0.0;
		double a[4], b[4], Area = 0.0;
		unsigned short counterPos = 0, counterNeg = 0;
		unsigned short loc;
		unsigned short iSpecies;

		if (nDim == 2) {
			for (unsigned short iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);			/*--- Norm of the normal component of area, area = 1/2*cross(a,b) ---*/


			for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);			
				if (config->GetParticle_ChargeNumber(iSpecies) > 0) {
					rho_Pos_0 += U_0[loc+0];
					rho_Pos_1 += U_1[loc+0];
					rho_Pos_2 += U_2[loc+0];
					mass_Pos_0 += U_0[loc+0] * config->GetMolar_Mass(iSpecies);
					mass_Pos_1 += U_1[loc+0] * config->GetMolar_Mass(iSpecies);
					mass_Pos_2 += U_2[loc+0] * config->GetMolar_Mass(iSpecies);
					massTotal_Pos_0 += U_0[loc+0];
					massTotal_Pos_1 += U_1[loc+0];
					massTotal_Pos_2 += U_2[loc+0];
					diff_ru_Pos_x_0 += ConsVar_Grad_0[loc+1][0]; 
					diff_ru_Pos_x_1 += ConsVar_Grad_1[loc+1][0];
					diff_ru_Pos_x_2 += ConsVar_Grad_2[loc+1][0];
					diff_ru_Pos_y_0 += ConsVar_Grad_0[loc+2][1];
					diff_ru_Pos_y_1 += ConsVar_Grad_1[loc+2][1];
					diff_ru_Pos_y_2 += ConsVar_Grad_2[loc+2][1];
					counterPos++;
				} else if (config->GetParticle_ChargeNumber(iSpecies) < 0) {
					rho_Neg_0 += U_0[loc+0];
					rho_Neg_1 += U_1[loc+0];
					rho_Neg_2 += U_2[loc+0];
					mass_Neg_0 += U_0[loc+0] * config->GetMolar_Mass(iSpecies);
					mass_Neg_1 += U_1[loc+0] * config->GetMolar_Mass(iSpecies);
					mass_Neg_2 += U_2[loc+0] * config->GetMolar_Mass(iSpecies);
					massTotal_Neg_0 += U_0[loc+0];
					massTotal_Neg_1 += U_1[loc+0];
					massTotal_Neg_2 += U_2[loc+0];					
					diff_ru_Neg_x_0 += ConsVar_Grad_0[loc+1][0]; 
					diff_ru_Neg_x_1 += ConsVar_Grad_1[loc+1][0];
					diff_ru_Neg_x_2 += ConsVar_Grad_2[loc+1][0];
					diff_ru_Neg_y_0 += ConsVar_Grad_0[loc+2][1];
					diff_ru_Neg_y_1 += ConsVar_Grad_1[loc+2][1];
					diff_ru_Neg_y_2 += ConsVar_Grad_2[loc+2][1];
					counterNeg++;
				}
			}
			rho_Pos_0 = rho_Pos_0/counterPos;
			rho_Pos_1 = rho_Pos_1/counterPos;
			rho_Pos_2 = rho_Pos_2/counterPos;
			rho_Neg_0 = rho_Neg_0/counterNeg;
			rho_Neg_1 = rho_Neg_1/counterNeg;
			rho_Neg_2 = rho_Neg_2/counterNeg;

			mass_Pos_0 = mass_Pos_0 / massTotal_Pos_0;
			mass_Pos_1 = mass_Pos_1 / massTotal_Pos_1;
			mass_Pos_2 = mass_Pos_2 / massTotal_Pos_2;
			mass_Neg_0 = mass_Neg_0 / massTotal_Neg_0;
			mass_Neg_1 = mass_Neg_1 / massTotal_Neg_1;
			mass_Neg_2 = mass_Neg_2 / massTotal_Neg_2;

			diff_ru_Pos_x_0 = diff_ru_Pos_x_0 / counterPos;			
			diff_ru_Pos_x_1 = diff_ru_Pos_x_1 / counterPos;
			diff_ru_Pos_x_2 = diff_ru_Pos_x_2 / counterPos;
			diff_ru_Pos_y_0 = diff_ru_Pos_y_0 / counterPos;			
			diff_ru_Pos_y_1 = diff_ru_Pos_y_1 / counterPos;
			diff_ru_Pos_y_2 = diff_ru_Pos_y_2 / counterPos;

			diff_ru_Neg_x_0 = diff_ru_Neg_x_0 / counterNeg;			
			diff_ru_Neg_x_1 = diff_ru_Neg_x_1 / counterNeg;
			diff_ru_Neg_x_2 = diff_ru_Neg_x_2 / counterNeg;
			diff_ru_Neg_y_0 = diff_ru_Neg_y_0 / counterNeg;			
			diff_ru_Neg_y_1 = diff_ru_Neg_y_1 / counterNeg;
			diff_ru_Neg_y_2 = diff_ru_Neg_y_2 / counterNeg;


			/*--- Source q ---*/
			ne_0 = rho_Neg_0 / mass_Neg_0;
			ne_1 = rho_Neg_1 / mass_Neg_1;
			ne_2 = rho_Neg_2 / mass_Neg_2;
			np_0 = rho_Pos_0 / mass_Pos_0;
			np_1 = rho_Pos_1 / mass_Pos_1;
			np_2 = rho_Pos_2 / mass_Pos_2;

			Chi_n_0 = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne_0 - np_0);
			Chi_n_1 = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne_1 - np_1);
			Chi_n_2 = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne_2 - np_2);

			alpha_0 = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Neg_0 / (mass_Neg_0*mass_Neg_0) + rho_Pos_0 / (mass_Pos_0*mass_Pos_0));
			alpha_1 = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Neg_1 / (mass_Neg_1*mass_Neg_1) + rho_Pos_1 / (mass_Pos_1*mass_Pos_1));
			alpha_2 = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Neg_2 / (mass_Neg_2*mass_Neg_2) + rho_Pos_2 / (mass_Pos_2*mass_Pos_2));			

			Chi_np1_0 = 1.0/(1.0+alpha_0) * (Chi_n_0 - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Neg_x_0+diff_ru_Neg_y_0)/mass_Neg_0 - (diff_ru_Pos_x_0+diff_ru_Pos_y_0)/mass_Pos_0));
			Chi_np1_1 = 1.0/(1.0+alpha_1) * (Chi_n_1 - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Neg_x_1+diff_ru_Neg_y_1)/mass_Neg_1 - (diff_ru_Pos_x_1+diff_ru_Pos_y_1)/mass_Pos_1));
			Chi_np1_2 = 1.0/(1.0+alpha_2) * (Chi_n_2 - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Neg_x_2+diff_ru_Neg_y_2)/mass_Neg_2 - (diff_ru_Pos_x_2+diff_ru_Pos_y_2)/mass_Pos_2));

			/*--- Residual = transpose(N) * source (q) * Area  ---*/
			val_residual[0] =-1.0/3.0 * Chi_np1_0 * Area;
			val_residual[1] =-1.0/3.0 * Chi_np1_1 * Area;
			val_residual[2] =-1.0/3.0 * Chi_np1_2 * Area;

		}

		if (nDim == 3) {
			val_residual[0] = 0.0;
			val_residual[1] = 0.0;
			val_residual[2] = 0.0;
		}
	}
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

			/*--- Terms 1 & 2: -Fcv\B7nabla(TurbPsi_i) - Fs\B7TurbPsi_i ---*/
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

CSourcePieceWise_AdjLevelSet::CSourcePieceWise_AdjLevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { 
}

CSourcePieceWise_AdjLevelSet::~CSourcePieceWise_AdjLevelSet(void) { }

void CSourcePieceWise_AdjLevelSet::SetResidual(double *val_residual, CConfig *config) {}

CSourcePieceWise_LevelSet::CSourcePieceWise_LevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { 
}

CSourcePieceWise_LevelSet::~CSourcePieceWise_LevelSet(void) { }

void CSourcePieceWise_LevelSet::SetResidual(double *val_residual, CConfig *config) {}

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

CSourceRotationalFrame_AdjFlow::CSourceRotationalFrame_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceRotationalFrame_AdjFlow::~CSourceRotationalFrame_AdjFlow(void) { }

void CSourceRotationalFrame_AdjFlow::SetResidual(double *val_residual, CConfig *config) {

	unsigned short iDim;
	double phi[3] = {0,0,0};

	/*--- Retrieve the angular velocity vector ---*/
	double *Omega = config->GetOmega_FreeStreamND();

	/*--- Get adjoint velocity ---*/
	for(iDim = 0; iDim < nDim; iDim++)
		phi[iDim] = U_i[iDim+1];

	/*--- Compute the source term ---*/
	if (nDim == 2) {
		val_residual[0] = 0.0;
		val_residual[1] =  Omega[2]*phi[1]*Volume;
		val_residual[2] = -Omega[2]*phi[0]*Volume;
		val_residual[3] = 0.0;
	} else {
		val_residual[0] = 0.0;
		val_residual[1] = (Omega[2]*phi[1] - Omega[1]*phi[2])*Volume;
		val_residual[2] = (Omega[0]*phi[2] - Omega[2]*phi[0])*Volume;
		val_residual[3] = (Omega[1]*phi[0] - Omega[0]*phi[1])*Volume;
		val_residual[4] = 0.0;
	}

}

CSourceAxisymmetric_Flow::CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

}

CSourceAxisymmetric_Flow::~CSourceAxisymmetric_Flow(void) { }

void CSourceAxisymmetric_Flow::SetResidual(double *val_residual, CConfig *config) {

	double yinv, Pressure_i, Enthalpy_i, Velocity_i, sq_vel;
	unsigned short iDim;

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
		val_residual[2] = yinv*Volume*U_i[2]*U_i[2]/U_i[0];
		val_residual[3] = yinv*Volume*Enthalpy_i*U_i[2];
	}

}

CSourcePieceWise_Plasma::CSourcePieceWise_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics,
		unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar, val_nSpecies, val_nDiatomics, val_nMonatomics, config) {

	unsigned short iDim, iVar, iSpecies;
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Kb = 1.38E-23; 
	M1 = config->GetParticle_Mass(0);
	M3 = config->GetParticle_Mass(2);
	M2 = config->GetParticle_Mass(1);

	SourceVector = new double[nVar];
	SourceJacobian = new double*[nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar ++) {
		SourceVector[iVar] = 10.0;
		SourceJacobian[iVar] = new double [nVar];
		for (unsigned short jVar = 0; jVar < nVar; jVar ++)
			SourceJacobian[iVar][jVar] = 0.0;
	}


	ElectricField = new double [nDim];
	MagneticField = new double [3];
	VcrossB = new double* [nSpecies];
	velocity = new double *[nSpecies];

	EMF  = new double* [nSpecies];
	MagneticDipole = new double [3];
	for (iDim = 0; iDim < nDim; iDim ++)
		ElectricField[iDim] = 0.0;

	for (iVar = 0; iVar < 3; iVar ++) {
		MagneticField[iVar] = 0.0;
		MagneticDipole[iVar] = 0.0;

	}
	MagneticDipole[0] = config->GetMagneticDipole(0);;
	MagneticDipole[1] = config->GetMagneticDipole(1);
	MagneticDipole[2] = config->GetMagneticDipole(2);

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		velocity[iSpecies] = new double [nDim];
		EMF[iSpecies] = new double [nDim];

		for (iDim = 0; iDim < nDim; iDim ++) {
			velocity[iSpecies][iDim]  = 0.0;
			EMF[iSpecies][iDim]  = 0.0;
		}

		VcrossB[iSpecies] = new double [3];
		for (iVar = 0; iVar < 3; iVar ++)
			VcrossB[iSpecies][iVar] = 0.0;

	}

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
	//# ifdef NewSource

	/*	Ex = - ConsVar_Grad[0][0];
	Ey = - ConsVar_Grad[0][1];
	//Ey = 0.0;
	//Ez = Ey;
	if ( Ex < -10000)
		Ex = -10000.0;

	if (nDim == 3) 	Ex = 0.0;
	 */

	tol = 1E-60;
	double zero = 1E-15;

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

	P1 = (gam1-1)*(e1 - 0.5*(m1*m1 + n1*n1 + l1*l1)/r1);	// Partial Pressure of species 1
	P2 = (gam2-1)*(e2 - 0.5*(m2*m2 + n2*n2 + l2*l2)/r2);	// Partial Pressure of species 2
	P3 = (gam3-1)*(e3 - 0.5*(m3*m3 + n3*n3 + l3*l3)/r3);	// Partial Pressure of species 3

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

	for (iVar = 0; iVar < 3; iVar ++)
		MagneticField[iVar] = 0.0;

	mdotr = 0.0;
	distance = 0.0;

	for (iDim = 0; iDim < nDim; iDim ++ ) {
		mdotr +=MagneticDipole[iDim]*Coord_i[iDim];
		distance += Coord_i[iDim]*Coord_i[iDim];
	}

	distance = sqrt(distance);
	rpower5 = pow(distance,5);
	rcubed  = pow(distance,3);

	for (iDim = 0; iDim < nDim; iDim ++ )
		MagneticField[iDim] = MAGNETIC_CONSTANT/(4.0*PI_NUMBER) * (3.0 * Coord_i[iDim]* mdotr/rpower5 - MagneticDipole[iDim]/rcubed);

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		if (nDim ==2)
			VcrossB[iSpecies][iDim] = 0.0;
		else {
			VcrossB[iSpecies][0] = velocity[iSpecies][1]*MagneticField[2] - velocity[iSpecies][2]*MagneticField[1];
			VcrossB[iSpecies][1] = velocity[iSpecies][2]*MagneticField[0] - velocity[iSpecies][0]*MagneticField[2];
			VcrossB[iSpecies][2] = velocity[iSpecies][0]*MagneticField[1] - velocity[iSpecies][1]*MagneticField[0];
		}
	}

	for (iDim = 0; iDim < nDim; iDim ++ ) {
		if(config->GetElectricSolver()) ElectricField[iDim] = -ConsVar_Grad[0][iDim];
		else ElectricField[iDim] = 0.0;
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ )
			EMF[iSpecies][iDim] = ElectricField[iDim] + VcrossB[iSpecies][iDim];
	}

	/* Derivative of the first term in the source terms wrt all the conservative variables */
	unsigned short loc;

	loc = 0;
	SourceVector[loc+0] = M1*AvgNum*R;

	SourceJacobian[loc+0][loc+0] =  M1*AvgNum*dR_dr1;
	SourceJacobian[loc+0][loc+1] = M1*AvgNum*dR_dm1;
	SourceJacobian[loc+0][loc+2] = M1*AvgNum*dR_dn1;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] = M1*AvgNum*dR_dl1;

	SourceJacobian[loc+0][loc+nDim+1] = M1*AvgNum*dR_de1;

	/* Derivative of the second term in the source terms wrt all the conservative variables */

	SourceVector[loc+1] = r1*nu12*(m2/r2 - m1/r1)+ r1*nu13*(m3/r3 - m1/r1);

	SourceJacobian[loc+1][loc+0] = nu12*m2/r2 + nu13*m3/r3;
	SourceJacobian[loc+1][loc+1] = -nu12 - nu13;
	SourceJacobian[loc+1][loc+2] = 0.0;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] = 0;
	SourceJacobian[loc+1][loc+nDim+1] = r1*dv12_dT1*dT1_de1*(m2/r2 - m1/r1) + r1*dv13_dT1*dT1_de1*(m3/r3 - m1/r1);

	/* Derivative of the third term in the source terms wrt all the conservative variables */

	SourceVector[loc+2] = r1*nu12*(n2/r2 - n1/r1) + r1*nu13*(n3/r3 - n1/r1);

	SourceJacobian[loc+2][loc+0]  = nu12*n2/r2 + nu13*n3/r3;
	SourceJacobian[loc+2][loc+1]  = 0.0;
	SourceJacobian[loc+2][loc+2]  = -nu12 - nu13;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] = 0;
	SourceJacobian[loc+2][loc+nDim+1] =  r1*(n2/r2 -n1/r1)*dv12_dT1*dT1_de1 + r1*(n3/r3 - n1/r1)*dv13_dT1*dT1_de1;


	if (nDim == 3) {
		SourceVector[loc+3] = r1*nu12*(l2/r2 - l1/r1) + r1*nu13*(l3/r3 - l1/r1);

		SourceJacobian[loc+3][loc+0]  = nu12*l2/r2 + nu13*l3/r3;
		SourceJacobian[loc+3][loc+1]  = 0.0;
		SourceJacobian[loc+3][loc+2]  = 0.0;
		SourceJacobian[loc+3][loc+3] = -nu12 - nu13;;
		SourceJacobian[loc+3][loc+nDim+1] =  r1*(l2/r2 -l1/r1)*dv12_dT1*dT1_de1 + r1*(l3/r3 - l1/r1)*dv13_dT1*dT1_de1;
	}
	/* Derivative of the Fourth term in the source terms wrt all the conservative variables */

	SourceVector[loc+nDim+1] = m1*nu12*(m2/r2 - m1/r1) + m1*nu13*(m3/r3 - m1/r1) +
			n1*nu12*(n2/r2 - n1/r1) + n1*nu13*(n3/r3 - n1/r1) +
			l1*nu12*(l2/r2 - l1/r1) + l1*nu13*(l3/r3 - l1/r1) + QT1;

	SourceJacobian[loc+nDim+1][loc+0]  = 0.0;
	SourceJacobian[loc+nDim+1][loc+1]  = nu12*(m2/r2 -2*m1/r1) + nu13*(m3/r3 - 2*m1/r1) + dQT1_dT1*dT1_dm1;
	SourceJacobian[loc+nDim+1][loc+2]  = nu12*(n2/r2 -2*n1/r1) + nu13*(n3/r3 - 2*n1/r1) + dQT1_dT1*dT1_dn1;
	if (nDim == 3) SourceJacobian[loc+nDim+1][loc+3]  = nu12*(l2/r2 -2*l1/r1) + nu13*(l3/r3 - 2*l1/r1) + dQT1_dT1*dT1_dl1;

	SourceJacobian[loc+nDim+1][loc+nDim+1]  = (m1*(m2/r2 - m1/r1) + n1*(n2/r2 - n1/r1) + l1*(l2/r2 - l1/r1))*dv12_dT1*dT1_de1 +
			(m1*(m3/r3 - m1/r1) + n1*(n3/r3 - n1/r1) + l1*(l3/r3 - l1/r1))*dv13_dT1*dT1_de1 + dQT1_dT1*dT1_de1;




	/* Derivative of the fifth term in the source terms wrt all the conservative variables */
	loc = (nDim+2);

	SourceVector[loc+0]= -M2*AvgNum*R;

	SourceJacobian[loc+0][loc+0] = -M2*AvgNum*dR_dr2;
	SourceJacobian[loc+0][loc+1] = -M2*AvgNum*dR_dm2;
	SourceJacobian[loc+0][loc+2] = -M2*AvgNum*dR_dn2;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] =  -M2*AvgNum*dR_dl2;
	SourceJacobian[loc+0][loc+nDim+1] = -M2*AvgNum*dR_de2;

	/* Derivative of the Sixth term in the source terms with all the conservative variables */

	SourceVector[loc+1] = r2*ec*EMF[1][0]/M2 + r2*nu21*(m1/r1 - m2/r2) + r2*nu23*(m3/r3 - m2/r2);

	SourceJacobian[loc+1][loc+0] = nu21*m1/r1 + nu23*m3/r3 + ec*EMF[1][0]/M2;
	SourceJacobian[loc+1][loc+1] = -nu21 - nu23;
	SourceJacobian[loc+1][loc+2] = ec*MagneticField[2]/M2;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] =  -ec*MagneticField[1]/M2;
	SourceJacobian[loc+1][loc+nDim+1]  = 0 + r2*(m1/r1 - m2/r2)*dv21_dT2*dT2_de2 + r2*(m3/r3 - m2/r2)*dv23_dT2*dT2_de2;

	/* Derivative of the Seventh term in the source terms wrt all the conservative variables */

	SourceVector[loc+2]  = r2*ec*EMF[1][1]/M2 + r2*nu21*(n1/r1 - n2/r2) + r2*nu23*(n3/r3 - n2/r2);

	SourceJacobian[loc+2][loc+0] = nu21*n1/r1 + nu23*n3/r3 + ec*EMF[1][1]/M2;
	SourceJacobian[loc+2][loc+1] = -ec*MagneticField[2]/M2;
	SourceJacobian[loc+2][loc+2] = -nu21 - nu23 ;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] = ec*MagneticField[0]/M2;;
	SourceJacobian[loc+2][loc+nDim+1] = 0 + r2*(n1/r1 - n2/r2)*dv21_dT2*dT2_de2 + r2*(n3/r3 - n2/r2)*dv23_dT2*dT2_de2;

	if (nDim == 3) {
		SourceVector[loc+3] = r2*ec*EMF[1][2]/M2 + r2*nu21*(l1/r1 - l2/r2) + r2*nu23*(l3/r3 - l2/r2);

		SourceJacobian[loc+3][loc+0]  = nu12*l2/r2 + nu13*l3/r3;
		SourceJacobian[loc+3][loc+1]  = ec*MagneticField[1]/M2;
		SourceJacobian[loc+3][loc+2]  = -ec*MagneticField[0]/M2;
		SourceJacobian[loc+3][loc+3] = -nu12 - nu13;;
		SourceJacobian[loc+3][loc+nDim+1] =  r1*(l2/r2 -l1/r1)*dv12_dT1*dT1_de1 + r1*(l3/r3 - l1/r1)*dv13_dT1*dT1_de1;
	}

	/* Derivative of the Eight term in the source terms wrt all the conservative variables */

	SourceVector[loc+nDim+1] = m2*ec*EMF[1][0]/M2 + n2*ec*EMF[1][1]/M2  + l2*ec*EMF[1][2]/M2 +
			m2*nu21*(m1/r1 - m2/r2) + m2*nu23*(m3/r3 - m2/r2)+
			n2*nu21*(n1/r1 - n2/r2) + n2*nu23*(n3/r3 - n2/r2) +
			l2*nu21*(l1/r1 - l2/r2) + l2*nu23*(l3/r3 - l2/r2) + QT2;

	SourceJacobian[loc+nDim+1][loc+0] = -1/r2*(m2*ec*VcrossB[1][0]/M2 + n2*ec*VcrossB[1][1]/M2  + l2*ec*VcrossB[1][2]/M2) ;//(nu21+nu23)*(m2*m2 + n2*n2)/(r2*r2) ;// +  dQT2_dT2*dT2_dr2;
	SourceJacobian[loc+nDim+1][loc+1] = nu21*(m1/r1 - 2*m2/r2) + nu23*(m3/r3 - 2*m2/r2) + ec*EMF[1][0]/M2  + dQT2_dT2*dT2_dm2;
	SourceJacobian[loc+nDim+1][loc+2] = nu21*(n1/r1 - 2*n2/r2) + nu23*(n3/r3 - 2*n2/r2) + ec*EMF[1][1]/M2  + dQT2_dT2*dT2_dn2;
	if (nDim ==3) SourceJacobian[loc+nDim+1][loc+3] = nu21*(l1/r1 - 2*l2/r2) + nu23*(l3/r3 - 2*l2/r2) + ec*EMF[1][2]/M2  + dQT2_dT2*dT2_dl2;
	SourceJacobian[loc+nDim+1][loc+nDim+1] = (n2*(n1/r1 - n2/r2)+ m2*(m1/r1 - m2/r2)+ l2*(l1/r1 - l2/r2))*dv21_dT2*dT2_de2+
			(n2*(n3/r3 - n2/r2)+ m2*(m3/r3 - m2/r2)+ l2*(l3/r3 - l2/r2))*dv23_dT2*dT2_de2 + dQT2_dT2*dT2_de2;



	//cout  << Coord_i[0]/0.01905 << ", " << Coord_i[1]/0.01905 << ", " <<  Coord_i[2]/0.01905 << ",  ";

	/* Derivative of the ninth term in the source terms wrt all the conservative variables */
	loc = 2*(nDim+2);
	SourceVector[loc+0] = -M3*AvgNum*R;

	SourceJacobian[loc+0][loc+0] = -M3*AvgNum*dR_dr3;
	SourceJacobian[loc+0][loc+1] = -M3*AvgNum*dR_dm3;
	SourceJacobian[loc+0][loc+2] = -M3*AvgNum*dR_dn3;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] = -M3*AvgNum*dR_dl3;
	SourceJacobian[loc+0][loc+nDim+1] = -M3*AvgNum*dR_de3;

	/* Derivative of the Tenth term in the source terms wrt all the conservative variables */

	SourceVector[loc+1] =  -r3*ec*EMF[2][0]/M3 + r3*nu31*(m1/r1 - m3/r3) + r3*nu32*(m2/r2 - m3/r3);
	//cout  << MagneticField[0] << ",  ";

	SourceJacobian[loc+1][loc+0] = nu31*m1/r1 + nu32*m2/r2 - ec*EMF[2][0]/M3;
	SourceJacobian[loc+1][loc+1] =  -nu31 - nu32;
	SourceJacobian[loc+1][loc+2] =  -ec*MagneticField[2]/M3;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] = ec*MagneticField[1]/M3;
	SourceJacobian[loc+1][loc+nDim+1] =  0 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_de3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_de3;

	/* Derivative of the Eleventh term in the source terms wrt all the conservative variables */
	SourceVector[loc+2]  = -r3*ec*EMF[2][1]/M3 + r3*nu31*(n1/r1 - n3/r3) + r3*nu32*(n2/r2 - n3/r3);
	//cout <<  MagneticField[1] << ",  " ;

	SourceJacobian[loc+2][loc+0] = nu31*n1/r1 + nu32*n2/r2 - ec*EMF[2][1]/M3 ;
	SourceJacobian[loc+2][loc+1] =  ec*MagneticField[2]/M3;
	SourceJacobian[loc+2][loc+2] =  -nu31 - nu32 ;
	if (nDim ==3)  SourceJacobian[loc+2][loc+3] = -ec*MagneticField[0]/M3;
	SourceJacobian[loc+2][loc+nDim+1] =  0 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_de3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_de3;

	/* Derivative of the Twelfth term in the source terms wrt all the conservative variables */
	if (nDim ==3) {
		SourceVector[loc+3]  = -r3*ec*EMF[2][2]/M3 + r3*nu31*(l1/r1 - l3/r3) + r3*nu32*(l2/r2 - l3/r3);
		//cout <<  MagneticField[2] << endl;

		SourceJacobian[loc+3][loc+0] = nu31*l1/r1 + nu32*l2/r2 - ec*EMF[2][2]/M3 ;
		SourceJacobian[loc+3][loc+1] =  -ec*MagneticField[1]/M3;
		SourceJacobian[loc+3][loc+2] =   ec*MagneticField[0]/M3;
		SourceJacobian[loc+3][loc+3] =  -nu31 - nu32 ;
		SourceJacobian[loc+3][loc+nDim+1] =  0 + r3*(l1/r1 - l3/r3)*dv31_dT3*dT3_de3 + r3*(l2/r2 - l3/r3)*dv32_dT3*dT3_de3;
	}

	SourceVector[loc+nDim+1]  = -m3*ec*EMF[2][0]/M3 - n3*ec*EMF[2][1]/M3 - l3*ec*EMF[2][2]/M3 +
			m3*nu31*(m1/r1 - m3/r3) + m3*nu32*(m2/r2 - m3/r3) +
			n3*nu31*(n1/r1 - n3/r3) + n3*nu32*(n2/r2 - n3/r3) +
			l3*nu31*(l1/r1 - l3/r3) + l3*nu32*(l2/r2 - l3/r3) + QT3;

	SourceJacobian[loc+nDim+1][loc+0] = (nu31+nu32)*(m3*m3 + n3*n3 +l3*l3)/(r3*r3) -1/r3*(-m3*ec*VcrossB[2][0]/M3 - n3*ec*VcrossB[2][1]/M3 - l3*ec*VcrossB[2][2]/M3);// + dQT3_dT3*dT3_dr3;
	SourceJacobian[loc+nDim+1][loc+1]  = nu31*(m1/r1 - 2*m3/r3) + nu32*(m2/r2 - 2*m3/r3) -ec*EMF[2][0]/M3 + dQT3_dT3*dT3_dm3;
	SourceJacobian[loc+nDim+1][loc+2]  = nu31*(n1/r1 - 2*n3/r3) + nu32*(n2/r2 - 2*n3/r3) -ec*EMF[2][1]/M3 + dQT3_dT3*dT3_dn3;
	if (nDim ==3)SourceJacobian[loc+nDim+1][loc+3]  = nu31*(l1/r1 - 2*l3/r3) + nu32*(l2/r2 - 2*l3/r3) -ec*EMF[2][2]/M3 + dQT3_dT3*dT3_dl3;
	SourceJacobian[loc+nDim+1][loc+nDim+1]  = (m3*(m1/r1 - m3/r3)+  n3*(n1/r1 - n3/r3)+  l3*(l1/r1 - l3/r3))*dv31_dT3*dT3_de3 +
			(m3*(m2/r2 - m3/r3)+  n3*(n2/r2 - n3/r3) + l3*(l2/r2 - l3/r3))*dv32_dT3*dT3_de3 + dQT3_dT3*dT3_de3;



	for (unsigned short iVar = 0; iVar < nVar; iVar ++) {
		val_residual[iVar] = SourceVector[iVar]*Volume;
		if (implicit) {
			for (unsigned short jVar = 0; jVar < nVar; jVar ++) {
				val_Jacobian_i[iVar][jVar] = SourceJacobian[iVar][jVar]*Volume;
			}
		}
	}
}

CSourcePieceWise_Plasma_Air::CSourcePieceWise_Plasma_Air(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics,
		unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar, val_nSpecies, val_nDiatomics, val_nMonatomics, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Kb = BOLTZMANN_CONSTANT;

	zero = 1E-15;
	nS = nSpecies;
	Mass = new double [nSpecies+1];
	Molwt = new double [nSpecies];
	Species_Charge = new double [nSpecies];
	Cv_heatcap = new double [nSpecies];
	Gamma_species = new double [nSpecies];
	n = new double[nSpecies + 1];
	MassSource = new double [nSpecies + 1];
	MomentumSource = new double* [nSpecies];
	EnergySource = new double [nSpecies];

	Density = new double [nSpecies];
	Pressure = new double [nSpecies];
	Temperature = new double [nSpecies];
	Energy = new double [nSpecies];
	rhoU = new double* [nSpecies];
	velocity = new double *[nSpecies];

	ElectricField = new double [nDim];
	MagneticField = new double [3];
	VcrossB = new double* [nSpecies];
	EMF  = new double* [nSpecies];
	MagneticDipole = new double [3];
	for (iDim = 0; iDim < nDim; iDim ++)
		ElectricField[iDim] = 0.0;

	for (iVar = 0; iVar < 3; iVar ++) {
		MagneticField[iVar] = 0.0;
		MagneticDipole[iVar] = 0.0;

	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		VcrossB[iSpecies] = new double [3];
		for (iVar = 0; iVar < 3; iVar ++)
			VcrossB[iSpecies][iVar] = 0.0;
	}

	double AtomicOxygen, AtomicNitrogen;

	AtomicOxygen = 2.6567625500E-26; AtomicNitrogen = 2.3258669700E-26;

	Mass[0] = 	2*AtomicOxygen;
	Mass[1] = 	3*AtomicOxygen;
	Mass[2] = 	2*AtomicNitrogen;
	Mass[3] = 	AtomicNitrogen+AtomicOxygen;
	Mass[4] = 	2*AtomicNitrogen;
	Mass[5] = 	2*AtomicNitrogen;
	Mass[6] = 	2*AtomicNitrogen;
	Mass[7] = 	2*AtomicNitrogen;
	Mass[8] = 	2*AtomicNitrogen;
	Mass[9] = 	2*AtomicOxygen-ELECTRON_MASS;
	Mass[10] = 	2*AtomicNitrogen-ELECTRON_MASS;
	Mass[11] = 	2*AtomicOxygen+ELECTRON_MASS;
	Mass[12] = 	AtomicOxygen;
	Mass[13] = 	AtomicNitrogen;
	Mass[14] = 	AtomicOxygen;
	Mass[15] = 	AtomicNitrogen;
	Mass[16] = 	AtomicOxygen-ELECTRON_MASS;
	Mass[17] = 	AtomicOxygen+ELECTRON_MASS;
	Mass[18] = 	ELECTRON_MASS;

	dTemperature_drhoU = new double*[nSpecies];
	dTemperature_Energy = new double [nSpecies];
	QT =  new double [nSpecies];
	Collision_Freq_Heat_Transfer = new double* [nSpecies];
	Collision_Freq_Momentum_Tranfer = new double* [nSpecies];
	CollisionArea = new double* [nSpecies];
	CollisionVelo = new double* [nSpecies];
	dQT_dTemperature = new double *[nSpecies];


	nReactions = 64;
	Reactants = 0;
	Products = 1;
	Reactions = new unsigned short**[nReactions];
	RateofReaction = new double[nReactions];

	for (iReactions = 0; iReactions < nReactions; iReactions ++){
		Reactions[iReactions] = new unsigned short*[2];
		Reactions[iReactions][Reactants] = new unsigned short [3];
		Reactions[iReactions][Products] = new unsigned short [3];
	}

	for (iReactions = 0; iReactions < nReactions; iReactions ++){
		for (iVar = 0; iVar < 3; iVar ++ ) {
			Reactions[iReactions][Reactants][iVar] = 0;
			Reactions[iReactions][Products][iVar] = 0;
		}
		RateofReaction[iReactions] = 0.0;
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		Molwt[iSpecies] = Mass[iSpecies]*AVOGAD_CONSTANT;
		Species_Charge[iSpecies] = 0.0;
		n[iSpecies] = 1.0;
		if (iSpecies < nDiatomics) 	Gamma_species[iSpecies]  =  config->GetGammaDiatomic();
		else Gamma_species[iSpecies]  =  config->GetGammaMonatomic();
		Cv_heatcap[iSpecies] = 1/(Gamma_species[iSpecies]-1.0)*UNIVERSAL_GAS_CONSTANT/Molwt[iSpecies];

		MassSource[iSpecies] = 0.0;
		EnergySource[iSpecies] = 0.0;
		Density[iSpecies] = 0.0;
		Pressure[iSpecies] = 0.0;
		Temperature[iSpecies] = 0.0;
		Energy[iSpecies] = 0.0;
		dTemperature_Energy[iSpecies] = 0.0;
		QT[iSpecies] = 0.0;



		Collision_Freq_Heat_Transfer[iSpecies] = new double [nSpecies];
		Collision_Freq_Momentum_Tranfer[iSpecies] = new double [nSpecies];
		CollisionArea[iSpecies] = new double [nSpecies];
		CollisionVelo[iSpecies] = new double [nSpecies];

		dQT_dTemperature[iSpecies] = new double [nSpecies];

		rhoU[iSpecies] = new double [nDim];
		velocity[iSpecies] = new double [nDim];
		EMF[iSpecies] = new double [nDim];

		MomentumSource[iSpecies] = new double [nDim];
		dTemperature_drhoU[iSpecies] = new double [nDim];

		for (iDim = 0; iDim < nDim; iDim ++) {
			rhoU[iSpecies][iDim]  = 0.0;
			velocity[iSpecies][iDim]  = 0.0;
			EMF[iSpecies][iDim]  = 0.0;
			MomentumSource[iSpecies][iDim]  = 0.0;
			dTemperature_drhoU[iSpecies][iDim]  = 0.0;
		}
		for (jSpecies = 0; jSpecies < jSpecies; jSpecies ++) {
			Collision_Freq_Heat_Transfer[iSpecies][jSpecies]  = 0.0;
			Collision_Freq_Momentum_Tranfer[iSpecies][jSpecies]  = 0.0;
			CollisionArea[iSpecies][jSpecies]  = 0.0;
			CollisionVelo[iSpecies][jSpecies]  = 0.0;

			dQT_dTemperature[iSpecies][jSpecies]  = 0.0;
		}
	}

	n[nSpecies] = 1.0;
	MassSource[nSpecies] = 0.0;

	Species_Charge[9]  =  1.0;
	Species_Charge[10] =  1.0;
	Species_Charge[11] = -1.0;
	Species_Charge[16] =  1.0;
	Species_Charge[17] = -1.0;
	Species_Charge[18] = -1.0;

	SourceVector = new double[nVar];
	SourceVector_i = new double[nVar];
	SourceVector_j = new double[nVar];

	SourceJacobian = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar ++) {
		SourceVector[iVar] = 0.0;
		SourceVector_i[iVar] = 0.0;
		SourceVector_j[iVar] = 0.0;
		SourceJacobian[iVar] = new double [nVar];
		for (jVar = 0; jVar < nVar; jVar ++)
			SourceJacobian[iVar][jVar] = 0.0;
	}

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
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ )
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies ++ )
			CollisionArea[iSpecies][jSpecies] = PI_NUMBER*r12*r12; //TEMPORARY, CORRECT LATER

	/* The first three columns of Reactions has the indicies of the reactants and the last three columns has indicies of the products
	 * nS implies that there is no third reactant or product;
	 * */
	Reactions[0][0][0]=13;	Reactions[0][0][1]=0;	Reactions[0][0][2] =nS;		Reactions[0][1][0] =3;	Reactions[0][1][1]=12;	Reactions[0][1][2] =nS;
	Reactions[1][0][0]=13;	Reactions[1][0][1]=3;	Reactions[1][0][2] =nS;		Reactions[1][1][0] =2;	Reactions[1][1][1]=12;	Reactions[1][1][2] =nS;
	Reactions[2][0][0]=12;	Reactions[2][0][1]=12;	Reactions[2][0][2] =nS;		Reactions[2][1][0] =0;	Reactions[2][1][1]=nS;	Reactions[2][1][2] =nS;
	Reactions[3][0][0]=12;	Reactions[3][0][1]=13;	Reactions[3][0][2] =nS;		Reactions[3][1][0] =3;	Reactions[3][1][1]=nS;	Reactions[3][1][2] =nS;
	Reactions[4][0][0]=13;	Reactions[4][0][1]=13;	Reactions[4][0][2] =nS;		Reactions[4][1][0] =2;	Reactions[4][1][1]=nS;	Reactions[4][1][2] =nS;
	Reactions[5][0][0]=14;	Reactions[5][0][1]=nS;	Reactions[5][0][2] =nS;		Reactions[5][1][0] =12;	Reactions[5][1][1]=nS;	Reactions[5][1][2] =nS;
	Reactions[6][0][0]=14;	Reactions[6][0][1]=nS;	Reactions[6][0][2] =nS;		Reactions[6][1][0] =12;	Reactions[6][1][1]=nS;	Reactions[6][1][2] =nS;
	Reactions[7][0][0]=15;	Reactions[7][0][1]=0;	Reactions[7][0][2] =nS;		Reactions[7][1][0] =3;	Reactions[7][1][1]=14;	Reactions[7][1][2] =nS;
	Reactions[8][0][0]=15;	Reactions[8][0][1]=nS;	Reactions[8][0][2] =nS;		Reactions[8][1][0] =13;	Reactions[8][1][1]=nS;	Reactions[8][1][2] =nS;
	Reactions[9][0][0]=12;	Reactions[9][0][1]=0;	Reactions[9][0][2] =nS;		Reactions[9][1][0] =1;	Reactions[9][1][1]=nS;	Reactions[9][1][2] =nS;
	Reactions[10][0][0]=12;	Reactions[10][0][1]=1;	Reactions[10][0][2] =nS;	Reactions[10][1][0] =0;	Reactions[10][1][1]=0;	Reactions[10][1][2] =nS;
	Reactions[11][0][0]=1;	Reactions[11][0][1]=nS;	Reactions[11][0][2] =nS;	Reactions[11][1][0] =12;Reactions[11][1][1]=0;	Reactions[11][1][2] =nS;
	Reactions[12][0][0]=4;	Reactions[12][0][1]=0;	Reactions[12][0][2] =nS;	Reactions[12][1][0] =7;	Reactions[12][1][1]=12;	Reactions[12][1][2] =12;
	Reactions[13][0][0]=4;	Reactions[13][0][1]=12;	Reactions[13][0][2] =nS;	Reactions[13][1][0] =3;	Reactions[13][1][1]=15;	Reactions[13][1][2] =nS;
	Reactions[14][0][0]=4;	Reactions[14][0][1]=4;	Reactions[14][0][2] =nS;	Reactions[14][1][0] =6;	Reactions[14][1][1]=7;	Reactions[14][1][2] =nS;
	Reactions[15][0][0]=4;	Reactions[15][0][1]=nS;	Reactions[15][0][2] =nS;	Reactions[15][1][0] =7;	Reactions[15][1][1]=nS;	Reactions[15][1][2] =nS;
	Reactions[16][0][0]=4;	Reactions[16][0][1]=nS;	Reactions[16][0][2] =nS;	Reactions[16][1][0] =7;	Reactions[16][1][1]=nS;	Reactions[16][1][2] =nS;
	Reactions[17][0][0]=5;	Reactions[17][0][1]=nS;	Reactions[17][0][2] =nS;	Reactions[17][1][0] =4;	Reactions[17][1][1]=nS;	Reactions[17][1][2] =nS;
	Reactions[18][0][0]=5;	Reactions[18][0][1]=nS;	Reactions[18][0][2] =nS;	Reactions[18][1][0] =4;	Reactions[18][1][1]=nS;	Reactions[18][1][2] =nS;
	Reactions[19][0][0]=5;	Reactions[19][0][1]=nS;	Reactions[19][0][2] =nS;	Reactions[19][1][0] =4;	Reactions[19][1][1]=nS;	Reactions[19][1][2] =nS;
	Reactions[20][0][0]=5;	Reactions[20][0][1]=0;	Reactions[20][0][2] =nS;	Reactions[20][1][0] =7;	Reactions[20][1][1]=12;	Reactions[20][1][2] =12;
	Reactions[21][0][0]=8;	Reactions[21][0][1]=nS;	Reactions[21][0][2] =nS;	Reactions[21][1][0] =5;	Reactions[21][1][1]=nS;	Reactions[21][1][2] =nS;
	Reactions[22][0][0]=8;	Reactions[22][0][1]=0;	Reactions[22][0][2] =nS;	Reactions[22][1][0] =7;	Reactions[22][1][1]=12;	Reactions[22][1][2] =12;
	Reactions[23][0][0]=8;	Reactions[23][0][1]=3;	Reactions[23][0][2] =nS;	Reactions[23][1][0] =7;	Reactions[23][1][1]=13;	Reactions[23][1][2] =12;
	Reactions[24][0][0]=6;	Reactions[24][0][1]=nS;	Reactions[24][0][2] =nS;	Reactions[24][1][0] =5;	Reactions[24][1][1]=nS;	Reactions[24][1][2] =nS;
	Reactions[25][0][0]=6;	Reactions[25][0][1]=0;	Reactions[25][0][2] =nS;	Reactions[25][1][0] =7;	Reactions[25][1][1]=12;	Reactions[25][1][2] =12;
	Reactions[26][0][0]=13;	Reactions[26][0][1]=1;	Reactions[26][0][2] =nS;	Reactions[26][1][0] =3;	Reactions[26][1][1]=0;	Reactions[26][1][2] =nS;
	Reactions[27][0][0]=12;	Reactions[27][0][1]=0;	Reactions[27][0][2] =nS;	Reactions[27][1][0] =1;	Reactions[27][1][1]=nS;	Reactions[27][1][2] =nS;
	Reactions[28][0][0]=2;	Reactions[28][0][1]=nS;	Reactions[28][0][2] =nS;	Reactions[28][1][0] =4;	Reactions[28][1][1]=nS;	Reactions[28][1][2] =nS;
	Reactions[29][0][0]=2;	Reactions[29][0][1]=nS;	Reactions[29][0][2] =nS;	Reactions[29][1][0] =5;	Reactions[29][1][1]=nS;	Reactions[29][1][2] =nS;
	Reactions[30][0][0]=2;	Reactions[30][0][1]=nS;	Reactions[30][0][2] =nS;	Reactions[30][1][0] =6;	Reactions[30][1][1]=nS;	Reactions[30][1][2] =nS;
	Reactions[31][0][0]=2;	Reactions[31][0][1]=nS;	Reactions[31][0][2] =nS;	Reactions[31][1][0] =8;	Reactions[31][1][1]=nS;	Reactions[31][1][2] =nS;
	Reactions[32][0][0]=2;	Reactions[32][0][1]=nS;	Reactions[32][0][2] =nS;	Reactions[32][1][0] =13;Reactions[32][1][1]=13;	Reactions[32][1][2] =nS;
	Reactions[33][0][0]=0;	Reactions[33][0][1]=nS;	Reactions[33][0][2] =nS;	Reactions[33][1][0] =12;Reactions[33][1][1]=12;	Reactions[33][1][2] =nS;
	Reactions[34][0][0]=0;	Reactions[34][0][1]=nS;	Reactions[34][0][2] =nS;	Reactions[34][1][0] =12;Reactions[34][1][1]=14;	Reactions[34][1][2] =nS;
	Reactions[35][0][0]=10;	Reactions[35][0][1]=18;	Reactions[35][0][2] =nS;	Reactions[35][1][0] =13;Reactions[35][1][1]=13;	Reactions[35][1][2] =nS;
	Reactions[36][0][0]=10;	Reactions[36][0][1]=18;	Reactions[36][0][2] =nS;	Reactions[36][1][0] =13;Reactions[36][1][1]=15;	Reactions[36][1][2] =nS;
	Reactions[37][0][0]=9;	Reactions[37][0][1]=18;	Reactions[37][0][2] =nS;	Reactions[37][1][0] =12;Reactions[37][1][1]=12;	Reactions[37][1][2] =nS;
	Reactions[38][0][0]=18;	Reactions[38][0][1]=10;	Reactions[38][0][2] =nS;	Reactions[38][1][0] =2;	Reactions[38][1][1]=nS;	Reactions[38][1][2] =nS;
	Reactions[39][0][0]=18;	Reactions[39][0][1]=nS;	Reactions[39][0][2] =9;		Reactions[39][1][0] =nS;Reactions[39][1][1]=0;	Reactions[39][1][2] =nS;
	Reactions[40][0][0]=16;	Reactions[40][0][1]=18;	Reactions[40][0][2] =nS;	Reactions[40][1][0] =12;Reactions[40][1][1]=nS;	Reactions[40][1][2] =nS;
	Reactions[41][0][0]=18;	Reactions[41][0][1]=0;	Reactions[41][0][2] =nS;	Reactions[41][1][0] =11;Reactions[41][1][1]=nS;	Reactions[41][1][2] =nS;
	Reactions[42][0][0]=18;	Reactions[42][0][1]=0;	Reactions[42][0][2] =nS;	Reactions[42][1][0] =11;Reactions[42][1][1]=nS;	Reactions[42][1][2] =nS;
	Reactions[43][0][0]=18;	Reactions[43][0][1]=12;	Reactions[43][0][2] =nS;	Reactions[43][1][0] =17;Reactions[43][1][1]=nS;	Reactions[43][1][2] =nS;
	Reactions[44][0][0]=18;	Reactions[44][0][1]=nS;	Reactions[44][0][2] =0;		Reactions[44][1][0] =nS;Reactions[44][1][1]=11;	Reactions[44][1][2] =nS;
	Reactions[45][0][0]=18;	Reactions[45][0][1]=1;	Reactions[45][0][2] =nS;	Reactions[45][1][0] =11;Reactions[45][1][1]=12;	Reactions[45][1][2] =nS;
	Reactions[46][0][0]=18;	Reactions[46][0][1]=1;	Reactions[46][0][2] =nS;	Reactions[46][1][0] =17;Reactions[46][1][1]=0;	Reactions[46][1][2] =nS;
	Reactions[47][0][0]=11;	Reactions[47][0][1]=2;	Reactions[47][0][2] =nS;	Reactions[47][1][0] =0;	Reactions[47][1][1]=2;	Reactions[47][1][2] =18;
	Reactions[48][0][0]=11;	Reactions[48][0][1]=nS;	Reactions[48][0][2] =nS;	Reactions[48][1][0] =0;	Reactions[48][1][1]=nS;	Reactions[48][1][2] =18;
	Reactions[49][0][0]=11;	Reactions[49][0][1]=12;	Reactions[49][0][2] =nS;	Reactions[49][1][0] =1;	Reactions[49][1][1]=18;	Reactions[49][1][2] =nS;
	Reactions[50][0][0]=17;	Reactions[50][0][1]=12;	Reactions[50][0][2] =nS;	Reactions[50][1][0] =0;	Reactions[50][1][1]=18;	Reactions[50][1][2] =nS;
	Reactions[51][0][0]=17;	Reactions[51][0][1]=13;	Reactions[51][0][2] =nS;	Reactions[51][1][0] =3;	Reactions[51][1][1]=18;	Reactions[51][1][2] =nS;
	Reactions[52][0][0]=17;	Reactions[52][0][1]=0;	Reactions[52][0][2] =nS;	Reactions[52][1][0] =1;	Reactions[52][1][1]=18;	Reactions[52][1][2] =nS;
	Reactions[53][0][0]=11;	Reactions[53][0][1]=4;	Reactions[53][0][2] =nS;	Reactions[53][1][0] =0;	Reactions[53][1][1]=2;	Reactions[53][1][2] =18;
	Reactions[54][0][0]=11;	Reactions[54][0][1]=5;	Reactions[54][0][2] =nS;	Reactions[54][1][0] =0;	Reactions[54][1][1]=2;	Reactions[54][1][2] =18;
	Reactions[55][0][0]=17;	Reactions[55][0][1]=4;	Reactions[55][0][2] =nS;	Reactions[55][1][0] =12;Reactions[55][1][1]=2;	Reactions[55][1][2] =18;
	Reactions[56][0][0]=17;	Reactions[56][0][1]=5;	Reactions[56][0][2] =nS;	Reactions[56][1][0] =12;Reactions[56][1][1]=2;	Reactions[56][1][2] =18;
	Reactions[57][0][0]=2;	Reactions[57][0][1]=nS;	Reactions[57][0][2] =nS;	Reactions[57][1][0] =10;Reactions[57][1][1]=18;	Reactions[57][1][2] =nS;
	Reactions[58][0][0]=0;	Reactions[58][0][1]=nS;	Reactions[58][0][2] =nS;	Reactions[58][1][0] =9;	Reactions[58][1][1]=18;	Reactions[58][1][2] =nS;
	Reactions[59][0][0]=0;	Reactions[59][0][1]=nS;	Reactions[59][0][2] =nS;	Reactions[59][1][0] =12;Reactions[59][1][1]=16;	Reactions[59][1][2] =18;
	Reactions[60][0][0]=11;	Reactions[60][0][1]=12;	Reactions[60][0][2] =nS;	Reactions[60][1][0] =0;	Reactions[60][1][1]=17;	Reactions[60][1][2] =nS;
	Reactions[61][0][0]=16;	Reactions[61][0][1]=0;	Reactions[61][0][2] =nS;	Reactions[61][1][0] =9;	Reactions[61][1][1]=12;	Reactions[61][1][2] =nS;
	Reactions[62][0][0]=10;	Reactions[62][0][1]=0;	Reactions[62][0][2] =nS;	Reactions[62][1][0] =9;	Reactions[62][1][1]=2;	Reactions[62][1][2] =nS;
	Reactions[63][0][0]=10;	Reactions[63][0][1]=1;	Reactions[63][0][2] =nS;	Reactions[63][1][0] =9;	Reactions[63][1][1]=12;	Reactions[63][1][2] =2;

}

CSourcePieceWise_Plasma_Air::~CSourcePieceWise_Plasma_Air(void) {

}

void CSourcePieceWise_Plasma_Air::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	double rhoUsqr;

	//cout << " in source.cpp " <<" nSpecies = " << nSpecies << " nDiatomics = " << nDiatomics << endl;

	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		Density[iSpecies] = U_i[loc + 0];
		Energy[iSpecies] = U_i[loc + nDim+1];
		n[iSpecies] = Density[iSpecies]/Mass[iSpecies]; 			// m-3
		rhoUsqr = 0.0;
		for (iDim = 0; iDim < nDim; iDim ++) {
			rhoU[iSpecies][iDim]  = U_i[loc + iDim+1];
			rhoUsqr += rhoU[iSpecies][iDim]*rhoU[iSpecies][iDim];
			velocity[iSpecies][iDim] = rhoU[iSpecies][iDim]/ Density[iSpecies];
		}
		Pressure[iSpecies] = (Gamma_species[iSpecies]-1)*(Energy[iSpecies] - 0.5*(rhoUsqr)/Density[iSpecies]);
		Temperature[iSpecies] = Pressure[iSpecies]/(UNIVERSAL_GAS_CONSTANT/Molwt[iSpecies]*Density[iSpecies]);
	}

	n[nS] = 1.0;

	reactant_1 = Reactions[0][0][0];  reactant_2 = Reactions[0][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[0] = 1.5E-11*exp(-3600/T);
	reactant_1 = Reactions[1][0][0];  reactant_2 = Reactions[1][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[1] = 2.1E-11*exp(100/T);
	reactant_1 = Reactions[2][0][0];  reactant_2 = Reactions[2][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[2] = 4.5E-34*exp(630/T);
	reactant_1 = Reactions[3][0][0];  reactant_2 = Reactions[3][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[3] = 6.3E-33*exp(140/T);
	reactant_1 = Reactions[4][0][0];  reactant_2 = Reactions[4][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[4] = 8.3E-34*exp(500/T);
	reactant_1 = Reactions[5][0][0];  reactant_2 = Reactions[5][0][1];  T = Temperature[reactant_1]; 		RateofReaction[5] = 3.2E-11*exp(67/T);
	reactant_1 = Reactions[6][0][0];  reactant_2 = Reactions[6][0][1];  T = Temperature[reactant_1] ; 		RateofReaction[6] = 1.8E-11*exp(107/T);
	reactant_1 = Reactions[7][0][0];  reactant_2 = Reactions[7][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[7] = 1.0E-11*exp(-210/T);
	reactant_1 = Reactions[8][0][0];  reactant_2 = Reactions[8][0][1];  T = Temperature[reactant_1]; 		RateofReaction[8] = 5E-12*exp(-1620/T);

	reactant_1 = Reactions[9][0][0];  reactant_2 = Reactions[9][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[9] = 1.85E-35*exp(1057/T);
	reactant_1 = Reactions[10][0][0]; reactant_2 = Reactions[10][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[10] = 8E-12*exp(-2060/T);
	reactant_1 = Reactions[11][0][0]; reactant_2 = Reactions[11][0][1]; T = Temperature[reactant_1]; 		RateofReaction[11] = 7.26E-10*exp(-11400/T);

	RateofReaction[12] = 2.54E-12; 	RateofReaction[13] = 7E-12;	  RateofReaction[14] = 2E-12;
	RateofReaction[15] = 3E-18;		RateofReaction[16] = 7E-11;   RateofReaction[17] = 5E-11;
	RateofReaction[18] = 150000; 	RateofReaction[19] = 2.4E-10; RateofReaction[20] = 3E-10;
	RateofReaction[21] = 2E-13; 	RateofReaction[22] = 2.8E-11; RateofReaction[23] = 3.6E-11;


	Emag = 1E6;
	Nden = n[nSpecies-1];
	double EmagbyNden = 1E16*(Emag/Nden*1E4); // in units of 1E-16V/cm2

	RateofReaction[24] = 30000000; 	RateofReaction[25] = 3E-10;	  RateofReaction[26] = 2E-16;
	reactant_1 = Reactions[27][0][0]; reactant_2 = Reactions[27][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[27] = 6.2E-34*(300/T)*(300/T);

	constA = -8.57; constB = -11.3; RateofReaction[28] = pow(10,(constA+constB/(EmagbyNden)));
	//cout << " EmagbyNden = " << EmagbyNden  << endl;

	constA = -7.97; constB = -13.2; RateofReaction[29] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.82; constB = -22.0; RateofReaction[30] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -8.15; constB = -15.8; RateofReaction[31] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.96; constB = -22.7; RateofReaction[32] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -8.31; constB = -8.4; RateofReaction[33] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.86; constB = -17.2; RateofReaction[34] = pow(10,(constA+constB/(EmagbyNden)));

	Te = Temperature[nSpecies-1];

	RateofReaction[35] = 2.8E-7*sqrt(300/Te); 	RateofReaction[36] = 2E-7*sqrt(300/Te); 	RateofReaction[37] = 2E-7*(300/Te);
	RateofReaction[38] = 1E-19*pow((300/Te),4.5); 	RateofReaction[39] = 6E-27*pow((300/Te),1.5); 	RateofReaction[40] = 0.14E-7*pow(Te,-4.5);

	reactant_1 = Reactions[41][0][0]; reactant_2 = Reactions[41][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[41] = 	1.4E-29*(300/Te)*exp(-600/T)*exp(700*(Te-T)/(Te*T));
	reactant_1 = Reactions[42][0][0]; reactant_2 = Reactions[42][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[42] = 	1.07E-31*pow((300/Te),2)*exp(-70/T)*exp(1500*(Te-T)/(Te*T));

	RateofReaction[43] = 1E-31; 		RateofReaction[44] = 	1E-31; 		RateofReaction[45] = 	1E-9; 							 RateofReaction[46] = 	1E-11;

	reactant_1 = Reactions[47][0][0]; reactant_2 = Reactions[47][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[47] = 	1.9E-12*pow((T/300),.5)*exp(-4990/T);
	reactant_1 = Reactions[48][0][0]; reactant_2 = Reactions[48][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[48] = 	2.7E-10*pow((T/300),.5)*exp(-5590/T);

	RateofReaction[49] = 1.5E-10; 		RateofReaction[50] = 	2.6E-10; 	RateofReaction[51] = 	2.6E-10; 	 RateofReaction[52] = 	5E-15;
	RateofReaction[53] = 2.1E-9; 		RateofReaction[54] = 	2.5E-9; 	RateofReaction[55] = 	2.2E-9; 	 RateofReaction[56] = 	1.9E-9;
	constA = -7.76; constB = -37.0; 	RateofReaction[57] =  pow(10,(constA+constB/(EmagbyNden)));
	constA = -8.34; constB = -30.7; 	RateofReaction[58] =  pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.94; constB = -32.2; 	RateofReaction[59] =  pow(10,(constA+constB/(EmagbyNden)));				 RateofReaction[60] = 	3.3E-10;
	RateofReaction[61] = 0.199E-10; 	RateofReaction[62] = 	0.6E-10;	RateofReaction[63] = 	1E-10;

	/*for (iReactions = 0; iReactions < nReactions; iReactions ++)
		cout << "Reaction = " << iReactions <<  " Rate of Reaction = " << RateofReaction[iReactions] << endl;

	cout << " ************* " << endl;*/
	for (iSpecies = 0; iSpecies < nSpecies+1; iSpecies ++ )
		MassSource[iSpecies] = 0.0;

	for (iReactions = 0; iReactions < nReactions; iReactions ++) {
		numdensity_Reactants = 1.0;
		for (iReactants = 0; iReactants < 3; iReactants ++) {
			reactant = Reactions[iReactions][0][iReactants];
			numdensity_Reactants = numdensity_Reactants*(n[reactant]*1E-6);   // Number density in cm-3;
		}
		for (iProducts = 0; iProducts < 3; iProducts ++) {
			product = Reactions[iReactions][1][iProducts];
			MassSource[product] += RateofReaction[iReactions]*numdensity_Reactants*1E6; // Mass source in 1/(m3.s), later multiplied by species mass
		}

		for (iReactants = 0; iReactants < 3; iReactants ++) {
			reactant = Reactions[iReactions][0][iReactants];
			MassSource[reactant] -= RateofReaction[iReactions]*numdensity_Reactants*1E6; // Mass source in 1/(m3.s), later multiplied by species mass
			/*		if (product == 12) {
					cout << " iSpecies = " << product << " being produced in reaction = " << iReactions << " amount = " << Mass[product]*RateofReaction[iReactions]*numdensity_Reactants*1E6/Density[product] << endl;
				cout << " Rate of reaction = " << RateofReaction[iReactions] << " number densities = " << numdensity_Reactants*1E6 << endl;
				cout << endl;

			}*/
		}
	}

	/*	MagneticDipole[0] = 0.0;//-11.0612;
		MagneticDipole[1] = 0.0;
		MagneticDipole[2] = 0.0;
		for (iVar = 0; iVar < 3; iVar ++)
			MagneticField[iVar] = 0.0;


		mdotr = 0.0;
		distance = 0.0;

		for (iDim = 0; iDim < nDim; iDim ++ ) {
			mdotr +=MagneticDipole[iDim]*Coord_i[iDim];
			distance += Coord_i[iDim]*Coord_i[iDim];
		}


		distance = sqrt(distance);
		rpower5 = pow(distance,5);
		rcubed  = pow(distance,3);

		for (iDim = 0; iDim < nDim; iDim ++ )
			MagneticField[iDim] = MAGNETIC_CONSTANT/(4.0*PI_NUMBER) * (3.0 * Coord_i[iDim]* mdotr/rpower5 - MagneticDipole[iDim]/rcubed);

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
			if (nDim ==2)
				VcrossB[iSpecies][iDim] = 0.0;
			else {
				VcrossB[iSpecies][0] = velocity[iSpecies][1]*MagneticField[2] - velocity[iSpecies][2]*MagneticField[1];
				VcrossB[iSpecies][1] = velocity[iSpecies][2]*MagneticField[0] - velocity[iSpecies][0]*MagneticField[2];
				VcrossB[iSpecies][2] = velocity[iSpecies][0]*MagneticField[1] - velocity[iSpecies][1]*MagneticField[0];
			}
		}

		for (iDim = 0; iDim < nDim; iDim ++ ) {
			ElectricField[iDim] = - ConsVar_Grad[0][iDim];
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ )
				EMF[iSpecies][iDim] = ElectricField[iDim] + VcrossB[iSpecies][iDim];
		}
	 */
	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
		for (jSpecies =0; jSpecies < nSpecies; jSpecies ++) {
			CollisionVelo[iSpecies][jSpecies] = sqrt(8.0*BOLTZMANN_CONSTANT/PI_NUMBER* ( Temperature[iSpecies]/Mass[iSpecies]+Temperature[jSpecies]/Mass[jSpecies]));
			Collision_Freq_Momentum_Tranfer[iSpecies][jSpecies] = Density[jSpecies]/(Mass[jSpecies]+ Mass[iSpecies]) * CollisionArea[iSpecies][jSpecies] * CollisionVelo[iSpecies][jSpecies];
		}
		CollisionVelo[iSpecies][iSpecies] = 0.0;
		Collision_Freq_Momentum_Tranfer[iSpecies][iSpecies] = 0.0;
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		for (iDim = 0; iDim < nDim; iDim ++ ) {
			MomentumSource[iSpecies][iDim] = Density[iSpecies]*ELECTRON_CHARGE*Species_Charge[iSpecies]/Mass[iSpecies]*EMF[iSpecies][iDim];
			for (jSpecies = 0; jSpecies < nSpecies; jSpecies ++)
				MomentumSource[iSpecies][iDim] += Density[iSpecies]*Collision_Freq_Momentum_Tranfer[iSpecies][jSpecies] * (velocity[jSpecies][iDim]-velocity[iSpecies][iDim]);
		}
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
		EnergySource[iSpecies]  = 0.0;
		for (iDim = 0; iDim < nDim; iDim ++ )
			EnergySource[iSpecies] += MomentumSource[iSpecies][iDim]*velocity[iSpecies][iDim];
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		loc = iSpecies* (nDim+2);
		SourceVector_i[loc + 0] 			= Mass[iSpecies]*MassSource[iSpecies];
		for (iDim = 0; iDim < nDim; iDim++)
			SourceVector_i[loc + iDim+1]  = MomentumSource[iSpecies][iDim];
		SourceVector_i[loc + nDim+1] 		= EnergySource[iSpecies];
	}










	/*SOURCE AT POINT J */


	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		Density[iSpecies] = U_j[loc + 0];
		Energy[iSpecies] = U_j[loc + nDim+1];
		n[iSpecies] = Density[iSpecies]/Mass[iSpecies]; 			// m-3
		rhoUsqr = 0.0;
		for (iDim = 0; iDim < nDim; iDim ++) {
			rhoU[iSpecies][iDim]  = U_j[loc + iDim+1];
			rhoUsqr += rhoU[iSpecies][iDim]*rhoU[iSpecies][iDim];
			velocity[iSpecies][iDim] = rhoU[iSpecies][iDim]/ Density[iSpecies];
		}
		Pressure[iSpecies] = (Gamma_species[iSpecies]-1)*(Energy[iSpecies] - 0.5*(rhoUsqr)/Density[iSpecies]);
		Temperature[iSpecies] = Pressure[iSpecies]/(UNIVERSAL_GAS_CONSTANT/Molwt[iSpecies]*Density[iSpecies]);
	}

	n[nS] = 1.0;

	reactant_1 = Reactions[0][0][0];  reactant_2 = Reactions[0][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[0] = 1.5E-11*exp(-3600/T);
	reactant_1 = Reactions[1][0][0];  reactant_2 = Reactions[1][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[1] = 2.1E-11*exp(100/T);
	reactant_1 = Reactions[2][0][0];  reactant_2 = Reactions[2][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[2] = 4.5E-34*exp(630/T);
	reactant_1 = Reactions[3][0][0];  reactant_2 = Reactions[3][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[3] = 6.3E-33*exp(140/T);
	reactant_1 = Reactions[4][0][0];  reactant_2 = Reactions[4][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[4] = 8.3E-34*exp(500/T);
	reactant_1 = Reactions[5][0][0];  reactant_2 = Reactions[5][0][1];  T = Temperature[reactant_1]; 		RateofReaction[5] = 3.2E-11*exp(67/T);
	reactant_1 = Reactions[6][0][0];  reactant_2 = Reactions[6][0][1];  T = Temperature[reactant_1] ; 		RateofReaction[6] = 1.8E-11*exp(107/T);
	reactant_1 = Reactions[7][0][0];  reactant_2 = Reactions[7][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[7] = 1.0E-11*exp(-210/T);
	reactant_1 = Reactions[8][0][0];  reactant_2 = Reactions[8][0][1];  T = Temperature[reactant_1]; 		RateofReaction[8] = 5E-12*exp(-1620/T);

	reactant_1 = Reactions[9][0][0];  reactant_2 = Reactions[9][0][1];  T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[9] = 1.85E-35*exp(1057/T);
	reactant_1 = Reactions[10][0][0]; reactant_2 = Reactions[10][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[10] = 8E-12*exp(-2060/T);
	reactant_1 = Reactions[11][0][0]; reactant_2 = Reactions[11][0][1]; T = Temperature[reactant_1]; 		RateofReaction[11] = 7.26E-10*exp(-11400/T);

	RateofReaction[12] = 2.54E-12; 	RateofReaction[13] = 7E-12;	  RateofReaction[14] = 2E-12;
	RateofReaction[15] = 3E-18;		RateofReaction[16] = 7E-11;   RateofReaction[17] = 5E-11;
	RateofReaction[18] = 150000; 	RateofReaction[19] = 2.4E-10; RateofReaction[20] = 3E-10;
	RateofReaction[21] = 2E-13; 	RateofReaction[22] = 2.8E-11; RateofReaction[23] = 3.6E-11;


	Nden = n[nSpecies-1];
	EmagbyNden = 1E16*(Emag/Nden*1E4); // in units of 1E-16V/cm2

	RateofReaction[24] = 30000000; 	RateofReaction[25] = 3E-10;	  RateofReaction[26] = 2E-16;
	reactant_1 = Reactions[27][0][0]; reactant_2 = Reactions[27][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[27] = 6.2E-34*(300/T)*(300/T);

	constA = -8.57; constB = -11.3; RateofReaction[28] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.97; constB = -13.2; RateofReaction[29] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.82; constB = -22.0; RateofReaction[30] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -8.15; constB = -15.8; RateofReaction[31] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.96; constB = -22.7; RateofReaction[32] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -8.31; constB = -8.4; RateofReaction[33] = pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.86; constB = -17.2; RateofReaction[34] = pow(10,(constA+constB/(EmagbyNden)));

	Te = Temperature[nSpecies-1];

	RateofReaction[35] = 2.8E-7*sqrt(300/Te); 	RateofReaction[36] = 2E-7*sqrt(300/Te); 	RateofReaction[37] = 2E-7*(300/Te);
	RateofReaction[38] = 1E-19*pow((300/Te),4.5); 	RateofReaction[39] = 6E-27*pow((300/Te),1.5); 	RateofReaction[40] = 0.14E-7*pow(Te,-4.5);

	reactant_1 = Reactions[41][0][0]; reactant_2 = Reactions[41][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[41] = 	1.4E-29*(300/Te)*exp(-600/T)*exp(700*(Te-T)/(Te*T));
	reactant_1 = Reactions[42][0][0]; reactant_2 = Reactions[42][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[42] = 	1.07E-31*pow((300/Te),2)*exp(-70/T)*exp(1500*(Te-T)/(Te*T));

	RateofReaction[43] = 1E-31; 		RateofReaction[44] = 	1E-31; 		RateofReaction[45] = 	1E-9; 							 RateofReaction[46] = 	1E-11;

	reactant_1 = Reactions[47][0][0]; reactant_2 = Reactions[47][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[47] = 	1.9E-12*pow((T/300),.5)*exp(-4990/T);
	reactant_1 = Reactions[48][0][0]; reactant_2 = Reactions[48][0][1]; T = sqrt(Temperature[reactant_1] * Temperature[reactant_2]); RateofReaction[48] = 	2.7E-10*pow((T/300),.5)*exp(-5590/T);

	RateofReaction[49] = 1.5E-10; 		RateofReaction[50] = 	2.6E-10; 	RateofReaction[51] = 	2.6E-10; 	 RateofReaction[52] = 	5E-15;
	RateofReaction[53] = 2.1E-9; 		RateofReaction[54] = 	2.5E-9; 	RateofReaction[55] = 	2.2E-9; 	 RateofReaction[56] = 	1.9E-9;
	constA = -7.76; constB = -37.0; 	RateofReaction[57] =  pow(10,(constA+constB/(EmagbyNden)));
	constA = -8.34; constB = -30.7; 	RateofReaction[58] =  pow(10,(constA+constB/(EmagbyNden)));
	constA = -7.94; constB = -32.2; 	RateofReaction[59] =  pow(10,(constA+constB/(EmagbyNden)));				 RateofReaction[60] = 	3.3E-10;
	RateofReaction[61] = 0.199E-10; 	RateofReaction[62] = 	0.6E-10;	RateofReaction[63] = 	1E-10;

	for (iSpecies = 0; iSpecies < nSpecies+1; iSpecies ++ )
		MassSource[iSpecies] = 0.0;

	for (iVar = 0; iVar < nVar; iVar++)
		for (jVar = 0; jVar < nVar; jVar++)
			SourceJacobian[iVar][jVar] = 0.0;

	for (iReactions = 0; iReactions < nReactions; iReactions ++) {
		numdensity_Reactants = 1.0;
		for (iReactants = 0; iReactants < 3; iReactants ++) {
			reactant = Reactions[iReactions][0][iReactants];
			numdensity_Reactants = numdensity_Reactants*(n[reactant]*1E-6);   // Number density in cm-3;
		}
		for (iProducts = 0; iProducts < 3; iProducts ++) {
			product = Reactions[iReactions][1][iProducts];
			MassSource[product] += RateofReaction[iReactions]*numdensity_Reactants*1E6 * Mass[product]; // Mass source in 1/(m3.s), later multiplied by species mass
		}

		for (iReactants = 0; iReactants < 3; iReactants ++) {
			reactant = Reactions[iReactions][0][iReactants];
			MassSource[reactant] -= RateofReaction[iReactions]*numdensity_Reactants*1E6 * Mass[reactant]; // Mass source in 1/(m3.s), later multiplied by species mass
		}
	}






	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
		for (jSpecies =0; jSpecies < nSpecies; jSpecies ++) {
			CollisionVelo[iSpecies][jSpecies] = sqrt(8.0*BOLTZMANN_CONSTANT/PI_NUMBER* ( Temperature[iSpecies]/Mass[iSpecies]+Temperature[jSpecies]/Mass[jSpecies]));
			Collision_Freq_Momentum_Tranfer[iSpecies][jSpecies] = Density[jSpecies]/(Mass[jSpecies]+ Mass[iSpecies]) * CollisionArea[iSpecies][jSpecies] * CollisionVelo[iSpecies][jSpecies];
		}
		CollisionVelo[iSpecies][iSpecies] = 0.0;
		Collision_Freq_Momentum_Tranfer[iSpecies][iSpecies] = 0.0;
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		for (iDim = 0; iDim < nDim; iDim ++ ) {
			MomentumSource[iSpecies][iDim] = Density[iSpecies]*ELECTRON_CHARGE*Species_Charge[iSpecies]/Mass[iSpecies]*EMF[iSpecies][iDim];
			for (jSpecies = 0; jSpecies < nSpecies; jSpecies ++)
				MomentumSource[iSpecies][iDim] += Density[iSpecies]*Collision_Freq_Momentum_Tranfer[iSpecies][jSpecies] * (velocity[jSpecies][iDim]-velocity[iSpecies][iDim]);
		}
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
		EnergySource[iSpecies]  = 0.0;
		for (iDim = 0; iDim < nDim; iDim ++ )
			EnergySource[iSpecies] += MomentumSource[iSpecies][iDim]*velocity[iSpecies][iDim];
	}

//	double netmasssrc = 0.0;
//	double netmomensrc = 0.0;
//	double netenergysrc = 0.0;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		loc = iSpecies* (nDim+2);
		SourceVector_j[loc + 0] 			= MassSource[iSpecies];
		for (iDim = 0; iDim < nDim; iDim++)
			SourceVector_j[loc + iDim+1]    = MomentumSource[iSpecies][iDim];
		SourceVector_j[loc + nDim+1] 		= EnergySource[iSpecies];
		/*		netmasssrc += SourceVector[loc + 0];
		netmomensrc += SourceVector[loc + 1];
		netmomensrc += SourceVector[loc + 2];
		netenergysrc += SourceVector[loc + nDim+1];
		 */
	}

	/*	cout << " net mass source = " << netmasssrc;
	cout << " net momentum source = " << netmomensrc;
	cout << " net energy source = " << netenergysrc << endl;
	 */

	for (iVar = 0; iVar < nVar; iVar++) {
		SourceVector[iVar] = 0.5*(SourceVector_i[iVar] + SourceVector_j[iVar]);

		if (implicit) {
			if ( iSpecies < nDiatomics ) {
				loc = (nDim+3)*iSpecies;
				nVar_species = nDim+3;
			}
			else {
				loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				nVar_species = nDim+2;
			}
			diff_Source = (SourceVector_i[iVar] - SourceVector_j[iVar]);
			cout << "iVar = " << iVar;
			for (jVar = loc + 0; jVar < loc + nVar_species; jVar++) {
				SourceJacobian[iVar][jVar] = diff_Source/(U_i[jVar] - U_j[jVar] + EPS);
			}
			cout << " " << U_i[iVar]<<endl;
		}
	}

	for (iVar = 0; iVar < nVar; iVar ++) {
		val_residual[iVar] = SourceVector[iVar]*Volume;
		if (implicit) {
			for (jVar = 0; jVar < nVar; jVar ++) {
				val_Jacobian_i[iVar][jVar] = SourceJacobian[iVar][jVar]*Volume;
			}
		}
	}
}
CSourcePieceWise_PlasmaDiatomic::CSourcePieceWise_PlasmaDiatomic(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	unsigned short iSpecies, iReaction;

	nSpecies = config->GetnSpecies();
	nReactions = config->GetnReactions();
	nReactions = 24;
	GammaMonatomic = config->GetGammaMonatomic();
	GammaDiatomic = config->GetGammaDiatomic();

	Molar_Mass = new double[nSpecies];
	ChargeNumber = new double[nSpecies];
	w_dot = new double[nSpecies];
	Q_tv = new double[nSpecies];
	Q_elastic = new double [nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Molar_Mass[iSpecies] = config->GetMolar_Mass(iSpecies);
		ChargeNumber[iSpecies] = config->GetParticle_ChargeNumber(iSpecies);
	}

	ExtentOfReaction = new double[nReactions];
	ReactionRateFwd  = new double[nReactions];
	ReactionRateBkw  = new double[nReactions];
	fwdRxn					 = new double[nReactions];
	bkwRxn					 = new double[nReactions];
	Keq							 = new double[nReactions];
	Cf							 = new double[nReactions];
	eta							 = new double[nReactions];
	theta						 = new double[nReactions];
	RxnReactants		 = new int*[nReactions];
	RxnProducts 		 = new int*[nReactions];
	EqRxnConstants	 = new double*[nReactions];
	for (iReaction = 0; iReaction < nReactions; iReaction++) {
		//		Cf[iReaction] = config->GetArrheniusCoeff(iReaction);
		//		eta[iReaction] = config->GetArrheniusEta(iReaction);
		//		theta[iReaction] = config->GetArrheniusTheta(iReaction);
		RxnReactants[iReaction] = new int [nSpecies];
		RxnProducts[iReaction]  = new int [nSpecies];
		EqRxnConstants[iReaction] = new double[5];
		//		for (ii = 0; ii < 5; ii++)
		//			EqRxnConstants[iReaction][ii] = config->GetEqRxnConstants(iReaction,ii);
	}	
	Molecular_Mass = new double [nSpecies];
	Molecular_Diameter = new double [nSpecies];
	P = new double*[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		P[iSpecies] = new double [nDim];

	CharVibTemp = new double [nDiatomics];

	ElectricField = new double [nDim];

	Residual_New = new double[nVar];
	Residual_Baseline = new double[nVar];
	U_Baseline = new double[nVar];


	/*--- Define matrices for reactants and products ---*/
	//Comment: Should come up with a slicker way of doing this in the long run...
	RxnReactants[0][0]  = 2;	RxnReactants[0][1]  = 0;	RxnReactants[0][2]  = 0;	RxnReactants[0][3]  = 0;	RxnReactants[0][4]  = 0;	RxnReactants[0][5]  = 0;	RxnReactants[0][6]  = 0;	
	RxnReactants[1][0]  = 1;	RxnReactants[1][1]  = 1;	RxnReactants[1][2]  = 0;	RxnReactants[1][3]  = 0;	RxnReactants[1][4]  = 0;	RxnReactants[1][5]  = 0;	RxnReactants[1][6]  = 0;
	RxnReactants[2][0]  = 1;	RxnReactants[2][1]  = 0;	RxnReactants[2][2]  = 1;	RxnReactants[2][3]  = 0;	RxnReactants[2][4]  = 0;	RxnReactants[2][5]  = 0;	RxnReactants[2][6]  = 0;	
	RxnReactants[3][0]  = 1;	RxnReactants[3][1]  = 0;	RxnReactants[3][2]  = 0;	RxnReactants[3][3]  = 1;	RxnReactants[3][4]  = 0;	RxnReactants[3][5]  = 0;	RxnReactants[3][6]  = 0;
	RxnReactants[4][0]  = 1;	RxnReactants[4][1]  = 0;	RxnReactants[4][2]  = 0;	RxnReactants[4][3]  = 0;	RxnReactants[4][4]  = 1;	RxnReactants[4][5]  = 0;	RxnReactants[4][6]  = 0;	
	RxnReactants[5][0]  = 1;	RxnReactants[5][1]  = 0;	RxnReactants[5][2]  = 0;	RxnReactants[5][3]  = 0;	RxnReactants[5][4]  = 0;	RxnReactants[5][5]  = 1;	RxnReactants[5][6]  = 0;	
	RxnReactants[6][0]  = 1;	RxnReactants[6][1]  = 0;	RxnReactants[6][2]  = 0;	RxnReactants[6][3]  = 0;	RxnReactants[6][4]  = 0;	RxnReactants[6][5]  = 0;	RxnReactants[6][6]  = 1;	
	RxnReactants[7][0]  = 1;	RxnReactants[7][1]  = 1;	RxnReactants[7][2]  = 0;	RxnReactants[7][3]  = 0;	RxnReactants[7][4]  = 0;	RxnReactants[7][5]  = 0;	RxnReactants[7][6]  = 0;	
	RxnReactants[8][0]  = 0;	RxnReactants[8][1]  = 2;	RxnReactants[8][2]  = 0;	RxnReactants[8][3]  = 0;	RxnReactants[8][4]  = 0;	RxnReactants[8][5]  = 0;	RxnReactants[8][6]  = 0;	
	RxnReactants[9][0]  = 0;	RxnReactants[9][1]  = 1;	RxnReactants[9][2]  = 1;	RxnReactants[9][3]  = 0;	RxnReactants[9][4]  = 0;	RxnReactants[9][5]  = 0;	RxnReactants[9][6]  = 0;
	RxnReactants[10][0] = 0;	RxnReactants[10][1] = 1;	RxnReactants[10][2] = 0;	RxnReactants[10][3] = 1;	RxnReactants[10][4] = 0;	RxnReactants[10][5] = 0;	RxnReactants[10][6] = 0;	
	RxnReactants[11][0] = 0;	RxnReactants[11][1] = 1;	RxnReactants[11][2] = 0;	RxnReactants[11][3] = 0;	RxnReactants[11][4] = 1;	RxnReactants[11][5] = 0;	RxnReactants[11][6] = 0;
	RxnReactants[12][0] = 0;	RxnReactants[12][1] = 1;	RxnReactants[12][2] = 0;	RxnReactants[12][3] = 0;	RxnReactants[12][4] = 0;	RxnReactants[12][5] = 1;	RxnReactants[12][6] = 0;	
	RxnReactants[13][0] = 0;	RxnReactants[13][1] = 1;	RxnReactants[13][2] = 0;	RxnReactants[13][3] = 0;	RxnReactants[13][4] = 0;	RxnReactants[13][5] = 0;	RxnReactants[13][6] = 1;	
	RxnReactants[14][0] = 1;	RxnReactants[14][1] = 0;	RxnReactants[14][2] = 1;	RxnReactants[14][3] = 0;	RxnReactants[14][4] = 0;	RxnReactants[14][5] = 0;	RxnReactants[14][6] = 0;	
	RxnReactants[15][0] = 0;	RxnReactants[15][1] = 1;	RxnReactants[15][2] = 1;	RxnReactants[15][3] = 0;	RxnReactants[15][4] = 0;	RxnReactants[15][5] = 0;	RxnReactants[15][6] = 0;	
	RxnReactants[16][0] = 0;	RxnReactants[16][1] = 0;	RxnReactants[16][2] = 2;	RxnReactants[16][3] = 0;	RxnReactants[16][4] = 0;	RxnReactants[16][5] = 0;	RxnReactants[16][6] = 0;
	RxnReactants[17][0] = 0;	RxnReactants[17][1] = 0;	RxnReactants[17][2] = 1;	RxnReactants[17][3] = 1;	RxnReactants[17][4] = 0;	RxnReactants[17][5] = 0;	RxnReactants[17][6] = 0;
	RxnReactants[18][0] = 0;	RxnReactants[18][1] = 0;	RxnReactants[18][2] = 1;	RxnReactants[18][3] = 0;	RxnReactants[18][4] = 1;	RxnReactants[18][5] = 0;	RxnReactants[18][6] = 0;
	RxnReactants[19][0] = 0;	RxnReactants[19][1] = 0;	RxnReactants[19][2] = 1;	RxnReactants[19][3] = 0;	RxnReactants[19][4] = 0;	RxnReactants[19][5] = 1;	RxnReactants[19][6] = 0;
	RxnReactants[20][0] = 0;	RxnReactants[20][1] = 0;	RxnReactants[20][2] = 1;	RxnReactants[20][3] = 0;	RxnReactants[20][4] = 0;	RxnReactants[20][5] = 0;	RxnReactants[20][6] = 1;
	RxnReactants[21][0] = 1;	RxnReactants[21][1] = 0;	RxnReactants[21][2] = 0;	RxnReactants[21][3] = 0;	RxnReactants[21][4] = 0;	RxnReactants[21][5] = 1;	RxnReactants[21][6] = 0;
	RxnReactants[22][0] = 0;	RxnReactants[22][1] = 0;	RxnReactants[22][2] = 1;	RxnReactants[22][3] = 0;	RxnReactants[22][4] = 0;	RxnReactants[22][5] = 1;	RxnReactants[22][6] = 0;
	RxnReactants[23][0] = 0;	RxnReactants[23][1] = 0;	RxnReactants[23][2] = 0;	RxnReactants[23][3] = 0;	RxnReactants[23][4] = 1;	RxnReactants[23][5] = 1;	RxnReactants[23][6] = 0;

	RxnProducts[0][0]  = 1;	RxnProducts[0][1]  = 0;	RxnProducts[0][2]  = 0;	RxnProducts[0][3]  = 0;	RxnProducts[0][4]  = 2;	RxnProducts[0][5]  = 0;	RxnProducts[0][6]  = 0;	
	RxnProducts[1][0]  = 0;	RxnProducts[1][1]  = 1;	RxnProducts[1][2]  = 0;	RxnProducts[1][3]  = 0;	RxnProducts[1][4]  = 2;	RxnProducts[1][5]  = 0;	RxnProducts[1][6]  = 0;
	RxnProducts[2][0]  = 0;	RxnProducts[2][1]  = 0;	RxnProducts[2][2]  = 1;	RxnProducts[2][3]  = 0;	RxnProducts[2][4]  = 2;	RxnProducts[2][5]  = 0;	RxnProducts[2][6]  = 0;	
	RxnProducts[3][0]  = 0;	RxnProducts[3][1]  = 0;	RxnProducts[3][2]  = 0;	RxnProducts[3][3]  = 1;	RxnProducts[3][4]  = 2;	RxnProducts[3][5]  = 0;	RxnProducts[3][6]  = 0;
	RxnProducts[4][0]  = 0;	RxnProducts[4][1]  = 0;	RxnProducts[4][2]  = 0;	RxnProducts[4][3]  = 0;	RxnProducts[4][4]  = 3;	RxnProducts[4][5]  = 0;	RxnProducts[4][6]  = 0;	
	RxnProducts[5][0]  = 0;	RxnProducts[5][1]  = 0;	RxnProducts[5][2]  = 0;	RxnProducts[5][3]  = 0;	RxnProducts[5][4]  = 2;	RxnProducts[5][5]  = 1;	RxnProducts[5][6]  = 0;	
	RxnProducts[6][0]  = 0;	RxnProducts[6][1]  = 0;	RxnProducts[6][2]  = 0;	RxnProducts[6][3]  = 0;	RxnProducts[6][4]  = 2;	RxnProducts[6][5]  = 0;	RxnProducts[6][6]  = 1;	
	RxnProducts[7][0]  = 1;	RxnProducts[7][1]  = 0;	RxnProducts[7][2]  = 0;	RxnProducts[7][3]  = 0;	RxnProducts[7][4]  = 0;	RxnProducts[7][5]  = 2;	RxnProducts[7][6]  = 0;	
	RxnProducts[8][0]  = 0;	RxnProducts[8][1]  = 1;	RxnProducts[8][2]  = 0;	RxnProducts[8][3]  = 0;	RxnProducts[8][4]  = 0;	RxnProducts[8][5]  = 2;	RxnProducts[8][6]  = 0;	
	RxnProducts[9][0]  = 0;	RxnProducts[9][1]  = 0;	RxnProducts[9][2]  = 1;	RxnProducts[9][3]  = 0;	RxnProducts[9][4]  = 0;	RxnProducts[9][5]  = 2;	RxnProducts[9][6]  = 0;
	RxnProducts[10][0] = 0;	RxnProducts[10][1] = 0;	RxnProducts[10][2] = 0;	RxnProducts[10][3] = 1;	RxnProducts[10][4] = 0;	RxnProducts[10][5] = 2;	RxnProducts[10][6] = 0;	
	RxnProducts[11][0] = 0;	RxnProducts[11][1] = 0;	RxnProducts[11][2] = 0;	RxnProducts[11][3] = 0;	RxnProducts[11][4] = 1;	RxnProducts[11][5] = 2;	RxnProducts[11][6] = 0;
	RxnProducts[12][0] = 0;	RxnProducts[12][1] = 0;	RxnProducts[12][2] = 0;	RxnProducts[12][3] = 0;	RxnProducts[12][4] = 0;	RxnProducts[12][5] = 3;	RxnProducts[12][6] = 0;	
	RxnProducts[13][0] = 0;	RxnProducts[13][1] = 0;	RxnProducts[13][2] = 0;	RxnProducts[13][3] = 0;	RxnProducts[13][4] = 0;	RxnProducts[13][5] = 2;	RxnProducts[13][6] = 1;	
	RxnProducts[14][0] = 1;	RxnProducts[14][1] = 0;	RxnProducts[14][2] = 0;	RxnProducts[14][3] = 0;	RxnProducts[14][4] = 1;	RxnProducts[14][5] = 1;	RxnProducts[14][6] = 0;	
	RxnProducts[15][0] = 0;	RxnProducts[15][1] = 1;	RxnProducts[15][2] = 0;	RxnProducts[15][3] = 0;	RxnProducts[15][4] = 1;	RxnProducts[15][5] = 1;	RxnProducts[15][6] = 0;	
	RxnProducts[16][0] = 0;	RxnProducts[16][1] = 0;	RxnProducts[16][2] = 1;	RxnProducts[16][3] = 0;	RxnProducts[16][4] = 1;	RxnProducts[16][5] = 1;	RxnProducts[16][6] = 0;
	RxnProducts[17][0] = 0;	RxnProducts[17][1] = 0;	RxnProducts[17][2] = 0;	RxnProducts[17][3] = 1;	RxnProducts[17][4] = 1;	RxnProducts[17][5] = 1;	RxnProducts[17][6] = 0;
	RxnProducts[18][0] = 0;	RxnProducts[18][1] = 0;	RxnProducts[18][2] = 0;	RxnProducts[18][3] = 0;	RxnProducts[18][4] = 2;	RxnProducts[18][5] = 1;	RxnProducts[18][6] = 0;
	RxnProducts[19][0] = 0;	RxnProducts[19][1] = 0;	RxnProducts[19][2] = 0;	RxnProducts[19][3] = 0;	RxnProducts[19][4] = 1;	RxnProducts[19][5] = 2;	RxnProducts[19][6] = 0;
	RxnProducts[20][0] = 0;	RxnProducts[20][1] = 0;	RxnProducts[20][2] = 0;	RxnProducts[20][3] = 0;	RxnProducts[20][4] = 1;	RxnProducts[20][5] = 1;	RxnProducts[20][6] = 1;
	RxnProducts[21][0] = 0;	RxnProducts[21][1] = 0;	RxnProducts[21][2] = 1;	RxnProducts[21][3] = 0;	RxnProducts[21][4] = 1;	RxnProducts[21][5] = 0;	RxnProducts[21][6] = 0;
	RxnProducts[22][0] = 0;	RxnProducts[22][1] = 1;	RxnProducts[22][2] = 0;	RxnProducts[22][3] = 0;	RxnProducts[22][4] = 1;	RxnProducts[22][5] = 0;	RxnProducts[22][6] = 0;
	RxnProducts[23][0] = 0;	RxnProducts[23][1] = 0;	RxnProducts[23][2] = 0;	RxnProducts[23][3] = 1;	RxnProducts[23][4] = 0;	RxnProducts[23][5] = 0;	RxnProducts[23][6] = 1;

	/*--- Set Arrhenius constants ---*/
	//Comment: Already implemented in CConfig, but doing this for convenience.
	//Scalabrin Thesis (Park '90)
	Cf[0]  = 7.0E21;
	Cf[1]  = 7.0E21;
	Cf[2]  = 7.0E21;
	Cf[3]  = 7.0E21;
	Cf[4]  = 3.0E22;
	Cf[5]  = 3.0E22;
	Cf[6]  = 0.0;
	Cf[7]  = 2.0E21;
	Cf[8]  = 2.0E21;
	Cf[9]  = 2.0E21;
	Cf[10] = 2.0E21;
	Cf[11] = 1.0E22;
	Cf[12] = 1.0E22;
	Cf[13] = 0.0;
	Cf[14] = 5.0E15;
	Cf[15] = 5.0E15;
	Cf[16] = 5.0E15;
	Cf[17] = 5.0E15;
	Cf[18] = 1.1E17;
	Cf[19] = 1.1E17;
	Cf[20] = 0.0;
	Cf[21] = 6.4E17;
	Cf[22] = 8.4E12;
	Cf[23] = 5.3E12;

	//Candler Thesis
	/*	Cf[0]  = 3.700E18;
	Cf[1]  = 3.700E18;
	Cf[2]  = 3.700E18;
	Cf[3]  = 3.700E18;
	Cf[4]  = 1.110E19;
	Cf[5]  = 1.110E19;
	Cf[6]  = 1.110E21;
	Cf[7]  = 2.750E16;
	Cf[8]  = 2.750E16;
	Cf[9]  = 2.750E16;
	Cf[10] = 2.750E16;
	Cf[11] = 8.250E16;
	Cf[12] = 8.250E16;
	Cf[13] = 1.320E19;
	Cf[14] = 2.300E14;
	Cf[15] = 2.300E14;
	Cf[16] = 2.300E14;
	Cf[17] = 2.300E14;
	Cf[18] = 4.600E14;
	Cf[19] = 4.600E14;
	Cf[20] = 7.360E16;
	Cf[21] = 3.180E10;
	Cf[22] = 2.160E5;
	Cf[23] = 6.500E8; */
	//Comment: Candler Appendix A has a seventh parameter for another reaction that he doesn't explicitly state...

	//Scalabrin Thesis
	eta[0]  = -1.60;
	eta[1]  = -1.60;
	eta[2]  = -1.60;
	eta[3]  = -1.60;
	eta[4]  = -1.60;
	eta[5]  = -1.60;
	eta[6]  = -1.60;
	eta[7]  = -1.50;
	eta[8]  = -1.50;
	eta[9]  = -1.50;
	eta[10] = -1.50;
	eta[11] = -1.50;
	eta[12] = -1.50;
	eta[13] = -1.50;
	eta[14] = 0.0;
	eta[15] = 0.0;
	eta[16] = 0.0;
	eta[17] = 0.0;
	eta[18] = 0.0;
	eta[19] = 0.0;
	eta[20] = 0.0;
	eta[21] = -1.00;
	eta[22] = 0.0;
	eta[23] = 0.0;

	//Candler Thesis
	/*	eta[0]  = -1.600;
	eta[1]  = -1.600;
	eta[2]  = -1.600;
	eta[3]  = -1.600;
	eta[4]  = -1.600;
	eta[5]  = -1.600;
	eta[6]  = -1.600;
	eta[7]  = -1.000;
	eta[8]  = -1.000;
	eta[9]  = -1.000;
	eta[10] = -1.000;
	eta[11] = -1.000;
	eta[12] = -1.000;
	eta[13] = -1.000;
	eta[14] = -0.500;
	eta[15] = -0.500;
	eta[16] = -0.500;
	eta[17] = -0.500;
	eta[18] = -0.500;
	eta[19] = -0.500;
	eta[20] = -0.500;
	eta[21] = 0.100;
	eta[22] = 1.290;
	eta[23] = 0.000;*/
	//Comment: Candler Appendix A has a seventh parameter for another reaction that he doesn't explicitly state...

	//Scalabrin Thesis
	theta[0] = 113200.0;
	theta[1] = 113200.0;
	theta[2] = 113200.0;
	theta[3] = 113200.0;
	theta[4] = 113200.0;
	theta[5] = 113200.0;
	theta[6] = 113200.0;
	theta[7] = 59500.0;
	theta[8] = 59500.0;
	theta[9] = 59500.0;
	theta[10] = 59500.0;
	theta[11] = 59500.0;
	theta[12] = 59500.0;
	theta[13] = 59500.0;
	theta[14] = 75500.0;
	theta[15] = 75500.0;
	theta[16] = 75500.0;
	theta[17] = 75500.0;
	theta[18] = 75500.0;
	theta[19] = 75500.0;
	theta[20] = 75500.0;
	theta[21] = 38400.0; 
	theta[22] = 19450.0;
	theta[23] = 31900.0;

	//Candler Thesis
	/*	theta[0]  = 113200.0;
	theta[1]  = 113200.0;
	theta[2]  = 113200.0;
	theta[3]  = 113200.0;
	theta[4]  = 113200.0;
	theta[5]  = 113200.0;
	theta[6]  = 113200.0;
	theta[7]  = 59500.0;
	theta[8]  = 59500.0;
	theta[9]  = 59500.0;
	theta[10] = 59500.0;
	theta[11] = 59500.0;
	theta[12] = 59500.0;
	theta[13] = 59500.0;
	theta[14] = 75500.0;
	theta[15] = 75500.0;
	theta[16] = 75500.0;
	theta[17] = 75500.0;
	theta[18] = 75500.0;
	theta[19] = 75500.0;
	theta[20] = 75500.0;
	theta[21] = 37700.0;
	theta[22] = 19220.0;
	theta[23] = 32000.0;*/
	//Comment: Candler Appendix A has a seventh parameter for another reaction that he doesn't explicitly state...

	/*--- Set equilibrium reaction constants ---*/
	//Scalabrin Thesis (for 1E17 number density)
	EqRxnConstants[0][0] = 1.531;	EqRxnConstants[0][1] = 1.6061;	EqRxnConstants[0][2] = 1.2993;	EqRxnConstants[0][3] = -11.494;	EqRxnConstants[0][4] = -.00698;
	EqRxnConstants[1][0] = 1.531;	EqRxnConstants[1][1] = 1.6061;	EqRxnConstants[1][2] = 1.2993;	EqRxnConstants[1][3] = -11.494;	EqRxnConstants[1][4] = -.00698;
	EqRxnConstants[2][0] = 1.531;	EqRxnConstants[2][1] = 1.6061;	EqRxnConstants[2][2] = 1.2993;	EqRxnConstants[2][3] = -11.494;	EqRxnConstants[2][4] = -.00698;
	EqRxnConstants[3][0] = 1.531;	EqRxnConstants[3][1] = 1.6061;	EqRxnConstants[3][2] = 1.2993;	EqRxnConstants[3][3] = -11.494;	EqRxnConstants[3][4] = -.00698;
	EqRxnConstants[4][0] = 1.531;	EqRxnConstants[4][1] = 1.6061;	EqRxnConstants[4][2] = 1.2993;	EqRxnConstants[4][3] = -11.494;	EqRxnConstants[4][4] = -.00698;
	EqRxnConstants[5][0] = 1.531;	EqRxnConstants[5][1] = 1.6061;	EqRxnConstants[5][2] = 1.2993;	EqRxnConstants[5][3] = -11.494;	EqRxnConstants[5][4] = -.00698;
	EqRxnConstants[6][0] = 1.531;	EqRxnConstants[6][1] = 1.6061;	EqRxnConstants[6][2] = 1.2993;	EqRxnConstants[6][3] = -11.494;	EqRxnConstants[6][4] = -.00698;

	EqRxnConstants[7][0]  = 0.55388;	EqRxnConstants[7][1]  = 2.46;	EqRxnConstants[7][2]  = 1.7763;	EqRxnConstants[7][3]  = -6.572;	EqRxnConstants[7][4]  = 0.031445;
	EqRxnConstants[8][0]  = 0.55388;	EqRxnConstants[8][1]  = 2.46;	EqRxnConstants[8][2]  = 1.7763;	EqRxnConstants[8][3]  = -6.572;	EqRxnConstants[8][4]  = 0.031445;
	EqRxnConstants[9][0]  = 0.55388;	EqRxnConstants[9][1]  = 2.46;	EqRxnConstants[9][2]  = 1.7763;	EqRxnConstants[9][3]  = -6.572;	EqRxnConstants[9][4]  = 0.031445;
	EqRxnConstants[10][0] = 0.55388;	EqRxnConstants[10][1] = 2.46;	EqRxnConstants[10][2] = 1.7763;	EqRxnConstants[10][3] = -6.572;	EqRxnConstants[10][4] = 0.031445;
	EqRxnConstants[11][0] = 0.55388;	EqRxnConstants[11][1] = 2.46;	EqRxnConstants[11][2] = 1.7763;	EqRxnConstants[11][3] = -6.572;	EqRxnConstants[11][4] = 0.031445;
	EqRxnConstants[12][0] = 0.55388;	EqRxnConstants[12][1] = 2.46;	EqRxnConstants[12][2] = 1.7763;	EqRxnConstants[12][3] = -6.572;	EqRxnConstants[12][4] = 0.031445;
	EqRxnConstants[13][0] = 0.55388;	EqRxnConstants[13][1] = 2.46;	EqRxnConstants[13][2] = 1.7763;	EqRxnConstants[13][3] = -6.572;	EqRxnConstants[13][4] = 0.031445;

	EqRxnConstants[14][0] = 0.55889;	EqRxnConstants[14][1] = 0.71558;	EqRxnConstants[14][2] = 0.55396;	EqRxnConstants[14][3] = -7.5304;	EqRxnConstants[14][4] = -0.014089;
	EqRxnConstants[15][0] = 0.55889;	EqRxnConstants[15][1] = 0.71558;	EqRxnConstants[15][2] = 0.55396;	EqRxnConstants[15][3] = -7.5304;	EqRxnConstants[15][4] = -0.014089;
	EqRxnConstants[16][0] = 0.55889;	EqRxnConstants[16][1] = 0.71558;	EqRxnConstants[16][2] = 0.55396;	EqRxnConstants[16][3] = -7.5304;	EqRxnConstants[16][4] = -0.014089;
	EqRxnConstants[17][0] = 0.55889;	EqRxnConstants[17][1] = 0.71558;	EqRxnConstants[17][2] = 0.55396;	EqRxnConstants[17][3] = -7.5304;	EqRxnConstants[17][4] = -0.014089;
	EqRxnConstants[18][0] = 0.55889;	EqRxnConstants[18][1] = 0.71558;	EqRxnConstants[18][2] = 0.55396;	EqRxnConstants[18][3] = -7.5304;	EqRxnConstants[18][4] = -0.014089;
	EqRxnConstants[19][0] = 0.55889;	EqRxnConstants[19][1] = 0.71558;	EqRxnConstants[19][2] = 0.55396;	EqRxnConstants[19][3] = -7.5304;	EqRxnConstants[19][4] = -0.014089;
	EqRxnConstants[20][0] = 0.55889;	EqRxnConstants[20][1] = 0.71558;	EqRxnConstants[20][2] = 0.55396;	EqRxnConstants[20][3] = -7.5304;	EqRxnConstants[20][4] = -0.014089;

	EqRxnConstants[21][0] = 0.97646;	EqRxnConstants[21][1] = 0.89043;	EqRxnConstants[21][2] = 0.74572;	EqRxnConstants[21][3] = -3.9642;	EqRxnConstants[21][4] = 0.007123;	
	EqRxnConstants[22][0] = 0.004815;	EqRxnConstants[22][1] = -1.7443;	EqRxnConstants[22][2] = -1.2227;	EqRxnConstants[22][3] = -0.95824;	EqRxnConstants[22][4] = -0.045545;
	EqRxnConstants[23][0] = -0.57924;	EqRxnConstants[23][1] = -7.3079;	EqRxnConstants[23][2] = -1.9999;	EqRxnConstants[23][3] = -3.2294;	EqRxnConstants[23][4] = 0.016382;

	//Candler Thesis
	/*	EqRxnConstants[0][0]  = 3.898;	EqRxnConstants[0][1]  = -12.611;	EqRxnConstants[0][2]  = 0.683;	EqRxnConstants[0][3]  = -0.118;	EqRxnConstants[0][4] = 0.006;
	EqRxnConstants[1][0]  = 3.898;	EqRxnConstants[1][1]  = -12.611;	EqRxnConstants[1][2]  = 0.683;	EqRxnConstants[1][3]  = -0.118;	EqRxnConstants[1][4] = 0.006;	
	EqRxnConstants[2][0]  = 3.898;	EqRxnConstants[2][1]  = -12.611;	EqRxnConstants[2][2]  = 0.683;	EqRxnConstants[2][3]  = -0.118;	EqRxnConstants[2][4] = 0.006;	
	EqRxnConstants[3][0]  = 3.898;	EqRxnConstants[3][1]  = -12.611;	EqRxnConstants[3][2]  = 0.683;	EqRxnConstants[3][3]  = -0.118;	EqRxnConstants[3][4] = 0.006;	
	EqRxnConstants[4][0]  = 3.898;	EqRxnConstants[4][1]  = -12.611;	EqRxnConstants[4][2]  = 0.683;	EqRxnConstants[4][3]  = -0.118;	EqRxnConstants[4][4] = 0.006;	
	EqRxnConstants[5][0]  = 3.898;	EqRxnConstants[5][1]  = -12.611;	EqRxnConstants[5][2]  = 0.683;	EqRxnConstants[5][3]  = -0.118;	EqRxnConstants[5][4] = 0.006;	
	EqRxnConstants[6][0]  = 3.898;	EqRxnConstants[6][1]  = -12.611;	EqRxnConstants[6][2]  = 0.683;	EqRxnConstants[6][3]  = -0.118;	EqRxnConstants[6][4] = 0.006;	
	EqRxnConstants[7][0]  = 1.335;	EqRxnConstants[7][1]  = -4.127;		EqRxnConstants[7][2]  = -0.616;	EqRxnConstants[7][3]  = 0.093;	EqRxnConstants[7][4] = -0.005;	
	EqRxnConstants[8][0]  = 1.335;	EqRxnConstants[8][1]  = -4.127;		EqRxnConstants[8][2]  = -0.616;	EqRxnConstants[8][3]  = 0.093;	EqRxnConstants[8][4] = -0.005;	
	EqRxnConstants[9][0]  = 1.335;	EqRxnConstants[9][1]  = -4.127;		EqRxnConstants[9][2]  = -0.616;	EqRxnConstants[9][3]  = 0.093;	EqRxnConstants[9][4] = -0.005;	
	EqRxnConstants[10][0] = 1.335;	EqRxnConstants[10][1] = -4.127;		EqRxnConstants[10][2] = -0.616;	EqRxnConstants[10][3] = 0.093;	EqRxnConstants[10][4] = -0.005;	
	EqRxnConstants[11][0] = 1.335;	EqRxnConstants[11][1] = -4.127;		EqRxnConstants[11][2] = -0.616;	EqRxnConstants[11][3] = 0.093;	EqRxnConstants[11][4] = -0.005;	
	EqRxnConstants[12][0] = 1.335;	EqRxnConstants[12][1] = -4.127;		EqRxnConstants[12][2] = -0.616;	EqRxnConstants[12][3] = 0.093;	EqRxnConstants[12][4] = -0.005;	
	EqRxnConstants[13][0] = 1.335;	EqRxnConstants[13][1] = -4.127;		EqRxnConstants[13][2] = -0.616;	EqRxnConstants[13][3] = 0.093;	EqRxnConstants[13][4] = -0.005;	
	EqRxnConstants[14][0] = 1.549;	EqRxnConstants[14][1] = -7.784;		EqRxnConstants[14][2] = 0.228;	EqRxnConstants[14][3] = -0.043;	EqRxnConstants[14][4] = 0.002;	
	EqRxnConstants[15][0] = 1.549;	EqRxnConstants[15][1] = -7.784;		EqRxnConstants[15][2] = 0.228;	EqRxnConstants[15][3] = -0.043;	EqRxnConstants[15][4] = 0.002;	
	EqRxnConstants[16][0] = 1.549;	EqRxnConstants[16][1] = -7.784;		EqRxnConstants[16][2] = 0.228;	EqRxnConstants[16][3] = -0.043;	EqRxnConstants[16][4] = 0.002;	
	EqRxnConstants[17][0] = 1.549;	EqRxnConstants[17][1] = -7.784;		EqRxnConstants[17][2] = 0.228;	EqRxnConstants[17][3] = -0.043;	EqRxnConstants[17][4] = 0.002;	
	EqRxnConstants[18][0] = 1.549;	EqRxnConstants[18][1] = -7.784;		EqRxnConstants[18][2] = 0.228;	EqRxnConstants[18][3] = -0.043;	EqRxnConstants[18][4] = 0.002;	
	EqRxnConstants[19][0] = 1.549;	EqRxnConstants[19][1] = -7.784;		EqRxnConstants[19][2] = 0.228;	EqRxnConstants[19][3] = -0.043;	EqRxnConstants[19][4] = 0.002;	
	EqRxnConstants[20][0] = 1.549;	EqRxnConstants[20][1] = -7.784;		EqRxnConstants[20][2] = 0.228;	EqRxnConstants[20][3] = -0.043;	EqRxnConstants[20][4] = 0.002;	
	EqRxnConstants[21][0] = 2.349;	EqRxnConstants[21][1] = -4.828;		EqRxnConstants[21][2] = 0.455;	EqRxnConstants[21][3] = -0.075;	EqRxnConstants[21][4] = 0.004;	
	EqRxnConstants[22][0] = 0.215;	EqRxnConstants[22][1] = -3.652;		EqRxnConstants[22][2] = 0.843;	EqRxnConstants[22][3] = -0.136;	EqRxnConstants[22][4] = 0.007;	
	EqRxnConstants[23][0] = -6.234;	EqRxnConstants[23][1] = -5.536;		EqRxnConstants[23][2] = 0.494;	EqRxnConstants[23][3] = -0.058;	EqRxnConstants[23][4] = 0.003;	*/
	//Comment: Candler Appendix A has a seventh parameter for another reaction that he doesn't explicitly state...

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Molar_Mass[iSpecies] = config->GetMolar_Mass(iSpecies);
		Molecular_Mass[iSpecies] = Molar_Mass[iSpecies] / AVOGAD_CONSTANT;
		Molecular_Diameter[iSpecies] = config->GetMolecular_Diameter(iSpecies);
	}

	CharVibTemp[0] = 3395.0;
	CharVibTemp[1] = 2239.0;
	CharVibTemp[2] = 2817.0;
	CharVibTemp[3] = 2817.0;

	tol = 1E-15;
}

CSourcePieceWise_PlasmaDiatomic::~CSourcePieceWise_PlasmaDiatomic(void) {

	delete [] Molar_Mass;		delete [] Molecular_Mass;		delete [] Molecular_Diameter;
	delete [] w_dot;				delete [] Q_tv;							delete [] Q_elastic;
	delete [] ChargeNumber;

	delete [] ExtentOfReaction;
	delete [] ReactionRateFwd;
	delete [] ReactionRateBkw;
	delete [] fwdRxn;		delete [] bkwRxn;
	delete [] Cf;	delete [] eta;	delete [] theta;

	unsigned short iReaction;
	for (iReaction = 0; iReaction < nReactions; iReaction++) {
		delete [] RxnReactants[iReaction];
		delete [] RxnProducts[iReaction];
		delete [] EqRxnConstants[iReaction];
	}
	delete [] RxnReactants;	delete [] RxnProducts;	delete [] EqRxnConstants;

	unsigned short iSpecies;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		delete [] P[iSpecies];
	delete [] P;

	delete [] CharVibTemp;	
	delete [] ElectricField;

	delete [] Residual_New;
	delete [] Residual_Baseline;
	delete [] U_Baseline;
}

void CSourcePieceWise_PlasmaDiatomic::SetJacobian_Chemistry(double **val_Jacobian_i, CConfig *config) {

	unsigned short iVar, jVar;

	/*--- In a near future AD should be used instead of FD ---*/
	double FDEpsilon = 1E-8;

	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
		Residual_New[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++)
			val_Jacobian_i[iVar][jVar] = 0.0;
	}

	/*--- Compute the baseline line residual ---*/
	SetResidual_Chemistry(Residual_Baseline, config);

	/*--- Compute forward finite differences ---*/
	for (iVar = 0; iVar < nVar; iVar++) {

		/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
		U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
		for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
		SetResidual_Chemistry(Residual_New, config);

		/*--- Undo the change in U_i to keep constant the solution ---*/
		U_i[iVar] = U_i[iVar] - FDEpsilon;

		/*--- Save the new Jacobian (per column) ---*/
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;			
		}
	}
}


void CSourcePieceWise_PlasmaDiatomic::GetEq_Rxn_Coefficients(double **EqRxnConstants, CConfig *config) {
	unsigned short iSpecies, iReaction, loc, ii;
	double mixNumDensity, spNumDensity, interpFact;

	if (config->GetKind_GasModel() == AIR7) {
		// Park 90 Chemistry model implemented in Scalabrin.

		/*--- Allocate memory for the equilibrium coefficients ---*/
		RxnConstantTable = new double*[6];
		for (ii = 0; ii < 6; ii++)
			RxnConstantTable[ii] = new double[5];

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
}

void CSourcePieceWise_PlasmaDiatomic::SetResidual_Chemistry(double *val_residual, CConfig *config) {
	double *T_rxnf;
	double *T_rxnb;
	T_rxnf = new double[nReactions];
	T_rxnb = new double[nReactions];

	unsigned short iSpecies, iReaction, iVar, iDim, loc, counterFwd, counterBkw;
	nReactions = config->GetnReactions();

	tol = 1E-60;

	for (iVar = 0; iVar < nVar; iVar++) val_residual[iVar] = 0.0;

	/*--- Determine the extent of each chemical reaction ---*/
	for (iReaction = 0; iReaction < nReactions; iReaction++) {
		T_rxnf[iReaction] = 1.0;
		T_rxnb[iReaction] = 1.0;
		counterFwd = 0;
		counterBkw = 0;
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			if (RxnReactants[iReaction][iSpecies] != 0) {
				T_rxnf[iReaction] *= pow(Temp_tr_i[iSpecies], RxnReactants[iReaction][iSpecies]);
				counterFwd += RxnReactants[iReaction][iSpecies];
			}
			if (RxnProducts[iReaction][iSpecies] != 0) {
				T_rxnb[iReaction] *= pow(Temp_tr_i[iSpecies], RxnProducts[iReaction][iSpecies]);
				counterBkw += RxnProducts[iReaction][iSpecies];
			}
		}
		T_rxnf[iReaction] = pow(T_rxnf[iReaction], 1.0/double(counterFwd));
		T_rxnb[iReaction] = pow(T_rxnb[iReaction], 1.0/double(counterBkw));

		//		GetEq_Rxn_Coefficients(EqRxnConstants, config);


		//-----------------------  Scalabrin Implementation  -----------------------
		Keq[iReaction] = exp(EqRxnConstants[iReaction][0]*(T_rxnb[iReaction]/10000.0) + EqRxnConstants[iReaction][1]
		                                                                                                          + EqRxnConstants[iReaction][2]*log(10000.0/T_rxnb[iReaction]) + EqRxnConstants[iReaction][3]*(10000.0/T_rxnb[iReaction])
		                                                                                                          + EqRxnConstants[iReaction][4]*pow(10000.0/T_rxnb[iReaction], 2.0));

		ReactionRateFwd[iReaction] = Cf[iReaction] * pow(T_rxnf[iReaction], eta[iReaction]) * exp(-theta[iReaction]/T_rxnf[iReaction]);
		ReactionRateBkw[iReaction] = ReactionRateFwd[iReaction] / Keq[iReaction];

		fwdRxn[iReaction] = 1.0;
		bkwRxn[iReaction] = 1.0;
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			fwdRxn[iReaction] = fwdRxn[iReaction]*pow(0.001*U_i[loc+0]/Molar_Mass[iSpecies], RxnReactants[iReaction][iSpecies]);
			bkwRxn[iReaction] = bkwRxn[iReaction]*pow(0.001*U_i[loc+0]/Molar_Mass[iSpecies], RxnProducts[iReaction][iSpecies]);
		}
		fwdRxn[iReaction] = 1000.0 * ReactionRateFwd[iReaction] * fwdRxn[iReaction];
		bkwRxn[iReaction] = 1000.0 * ReactionRateBkw[iReaction] * bkwRxn[iReaction];
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		w_dot[iSpecies] = 0.0;
		for (iReaction = 0; iReaction < nReactions; iReaction++) {
			w_dot[iSpecies] =   (RxnProducts[iReaction][iSpecies] - RxnReactants[iReaction][iSpecies])
																								* (fwdRxn[iReaction] - bkwRxn[iReaction]);
		}
		w_dot[iSpecies] = Molar_Mass[iSpecies]*w_dot[iSpecies];
		val_residual[loc] = w_dot[iSpecies]*Volume;
		for (iDim = 0; iDim < nDim; iDim++)
			val_residual[loc+1+iDim] = w_dot[iSpecies] * U_i[loc+1+iDim]/U_i[loc] * Volume;
		val_residual[loc+1+nDim] = w_dot[iSpecies] * U_i[loc+1+nDim]/U_i[loc] * Volume;
		if (iSpecies < nDiatomics)
			val_residual[loc+2+nDim] = w_dot[iSpecies] * U_i[loc+2+nDim]/U_i[loc] * Volume;
	}








	//-----------------------  Candler Implementation  -----------------------
	/*		Keq = exp(EqRxnConstants[iReaction][0]
		                                    + EqRxnConstants[iReaction][1]*1.0E4/T_rxnTR
		                                    + EqRxnConstants[iReaction][2]*pow(1.0E4/T_rxnTR,2)
		                                    + EqRxnConstants[iReaction][3]*pow(1.0E4/T_rxnTR, 3)
		                                    + EqRxnConstants[iReaction][4]*pow(1.0E4/T_rxnTR, 4));		

		ReactionRateFwd[iReaction] = Cf[iReaction] * pow(T_rxnTR,eta[iReaction]) * exp(-theta[iReaction]/T_rxnTR);
		ReactionRateBkw[iReaction] = ReactionRateFwd[iReaction] / Keq;	*/

	/*--- Calculate reaction rate constants in log-space ---*/
	/*		Keq = EqRxnConstants[iReaction][0] + EqRxnConstants[iReaction][1]*1.0E4/T_rxnTR + EqRxnConstants[iReaction][2]*pow(1.0E4/T_rxnTR,2)
				+ EqRxnConstants[iReaction][3]*pow(1.0E4/T_rxnTR, 3) + EqRxnConstants[iReaction][4]*pow(1.0E4/T_rxnTR, 4);
		ReactionRateFwd[iReaction] = log(Cf[iReaction]) + eta[iSpecies]*log(T_rxnTR) - theta[iReaction]/T_rxnTR;
		ReactionRateBkw[iReaction] = ReactionRateFwd[iReaction] - Keq; */

	/*--- Convert back to normal-space ---*/
	/*		ReactionRateBkw[iReaction] = exp(ReactionRateBkw[iReaction]);
		ReactionRateFwd[iReaction] = exp(ReactionRateFwd[iReaction]);*/

	/*		fwdRxn = 1.0 * ReactionRateFwd[iReaction];
		bkwRxn = 1.0 * ReactionRateBkw[iReaction];*/

	//		cout << "Reaction: " << iReaction << " - Keq: " << Keq << " fwdRxn: " << fwdRxn << " bkwRxn: " << bkwRxn << endl;
	/*		double StoichCoeff_;
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			StoichCoeff_ = 0;
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			if (iReaction == 8)
				cout << "hi" << endl;
			if (RxnReactants[iReaction][iSpecies] != 0) {//Reactants
				StoichCoeff_ = RxnReactants[iReaction][iSpecies];
				fwdRxn = fwdRxn * pow((U_i[loc+0] / Molar_Mass[iSpecies]), fabs(StoichCoeff_));
			}
			if (RxnProducts[iReaction][iSpecies] != 0) { //Products
				StoichCoeff_ = RxnProducts[iReaction][iSpecies];
				bkwRxn = bkwRxn * pow((U_i[loc+0] / Molar_Mass[iSpecies]), fabs(StoichCoeff_));	
				if (iReaction == 8)
					if (iSpecies == 5)
						cout << "In Bkw: " << StoichCoeff_ << endl;
			} */
	/*if (iReaction == 8) {
				if (StoichCoeff_ !=0) {
					cout << "iSpecies: " << iSpecies << endl;
					cout << "Stoich: " << StoichCoeff_ << endl;
					cout << "U_i[0]: " << U_i[loc] << endl;
					cout << "Keq: " << Keq << endl;
					cout << "Exponent: " << EqRxnConstants[iReaction][0] + EqRxnConstants[iReaction][1]*1.0E4/T_rxnTR + EqRxnConstants[iReaction][2]*pow(1.0E4/T_rxnTR,2) + EqRxnConstants[iReaction][3]*pow(1.0E4/T_rxnTR, 3) + EqRxnConstants[iReaction][4]*pow(1.0E4/T_rxnTR, 4) << endl;
					cout << "TrxTR: " << T_rxnTR << endl;
					cout << "Kf: " << ReactionRateFwd[iReaction] << endl;
					cout << "fwdRxn: " << fwdRxn << endl;
					cout << "Kb: " << ReactionRateBkw[iReaction] << endl;
					cout << "bkwRxn: " << bkwRxn << endl;
					cout << "R: " << ExtentOfReaction[iSpecies] << endl;
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						if (RxnProducts[iReaction][iSpecies] != 0)
							cout << "products coeff[" << iSpecies <<  "]: " << RxnProducts[iReaction][iSpecies] << endl;
						if (RxnReactants[iReaction][iSpecies] != 0)
							cout << "reactants coeff[" << iSpecies << "]: " << RxnReactants[iReaction][iSpecies] << endl;
					}
					cin.get();
				}
			}*/
	/*		}
		ExtentOfReaction[iReaction] = -fwdRxn + bkwRxn;		
	}*/

	/*--- Check values for forward and backward reaction rates for numerical stability ---*/
	/*if (ReactionRateBkw[iReaction] < numeric_limits<double>::epsilon())
		ReactionRateBkw[iReaction] = numeric_limits<double>::epsilon();
	else if (ReactionRateBkw[iReaction] > 1.0E8)
		ReactionRateBkw[iReaction] = 1.0E8;
	if (ReactionRateFwd[iReaction] < numeric_limits<double>::epsilon())
		ReactionRateFwd[iReaction] = numeric_limits<double>::epsilon();
	else if (ReactionRateFwd[iReaction] > 1.0E8)
		ReactionRateFwd[iReaction] = 1.0E8;*/

	/*--- Solve for mass production based on extent of reaction ---*/
	/*	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);			
		w_dot[iSpecies] = 0.0;
		for (iReaction = 0; iReaction < nReactions; iReaction++) {
				w_dot[iSpecies] = w_dot[iSpecies] - RxnProducts[iReaction][iSpecies]*ExtentOfReaction[iReaction] + RxnReactants[iReaction][iSpecies]*ExtentOfReaction[iReaction];
		}
		w_dot[iSpecies] = Molar_Mass[iSpecies] * w_dot[iSpecies];
		val_residual[loc] = w_dot[iSpecies] * Volume;
		for (iDim = 0; iDim < nDim; iDim++)
			val_residual[loc+1+iDim] = w_dot[iSpecies] * (U_i[loc+1+iDim]/U_i[loc]) * Volume;
		val_residual[loc+1+nDim] = w_dot[iSpecies] * (U_i[loc+1+nDim]/U_i[loc]) * Volume;
		if (iSpecies < nDiatomics)
			val_residual[loc+2+nDim] = w_dot[iSpecies] * (U_i[loc+2+nDim]/U_i[loc]) * Volume;
	}*/
}

void CSourcePieceWise_PlasmaDiatomic::SetJacobian_ElecForce(double **val_Jacobian_i, CConfig *config) {

	unsigned short iVar, jVar;

	/*--- In a near future AD should be used instead of FD ---*/
	double *Residual_Baseline, *Residual_New, *U_Baseline;
	double FDEpsilon = 1E-8;

	Residual_Baseline = new double[nVar];
	Residual_New = new double[nVar];
	U_Baseline = new double [nVar];

	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
	}

	/*--- Compute the baseline line residual ---*/

	SetResidual_ElecForce(Residual_Baseline, config);

	/*--- Compute forward finite differences ---*/
	for (iVar = 0; iVar < nVar; iVar++) {

		/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
		U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
		for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
		SetResidual_ElecForce(Residual_New, config);

		/*--- Undo the change in U_i to keep constant the solution ---*/
		U_i[iVar] = U_i[iVar] - FDEpsilon;

		/*--- Save the new Jacobian (per column) ---*/
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;			
		}
	}
	delete [] Residual_New; delete [] Residual_Baseline; delete [] U_Baseline;
}

void CSourcePieceWise_PlasmaDiatomic::SetResidual_ElecForce(double *val_residual, CConfig *config) {
	unsigned short iDim, iSpecies, loc;
	double F_es;

	/*--- Solve for the electrostatic force on each constituent ---*/
	for (iDim = 0; iDim < nDim; iDim++) { 
		ElectricField[iDim] = -ConsVar_Grad[0][iDim];
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);			
			F_es = config->GetParticle_ChargeNumber(iSpecies) * (U_i[loc+0]/Molar_Mass[iSpecies]) * ELECTRON_CHARGE * ElectricField[iDim];
			val_residual[loc+1+iDim] = F_es * Volume;
		}
	}

	/*--- Compute jacobian of electrostatic force terms ---*/

}

void CSourcePieceWise_PlasmaDiatomic::SetJacobian_MomentumExch(double **val_Jacobian_i, CConfig *config) {

	unsigned short iVar, jVar;

	/*--- In a near future AD should be used instead of FD ---*/
	double *Residual_Baseline, *Residual_New, *U_Baseline;
	double FDEpsilon = 1E-8;

	Residual_Baseline = new double[nVar];
	Residual_New = new double[nVar];
	U_Baseline = new double [nVar];


	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
	}

	/*--- Compute the baseline line residual ---*/
	SetResidual_MomentumExch(Residual_Baseline, config);

	/*--- Compute forward finite differences ---*/
	for (iVar = 0; iVar < nVar; iVar++) {

		/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
		U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
		for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
		SetResidual_MomentumExch(Residual_New, config);

		/*--- Undo the change in U_i to keep constant the solution ---*/
		U_i[iVar] = U_i[iVar] - FDEpsilon;

		/*--- Save the new Jacobian (per column) ---*/
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
		}
	}

	delete [] Residual_Baseline; delete [] Residual_New; delete [] U_Baseline;
}

void CSourcePieceWise_PlasmaDiatomic::SetResidual_MomentumExch(double *val_residual, CConfig *config) {
	unsigned short iDim, iSpecies, jSpecies, iVar, iLoc, jLoc;
	double collisionFreq;

	/*--- Initialization ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		for (iDim = 0; iDim < nDim; iDim++)
			P[iSpecies][iDim] = 0.0;
	for (iVar = 0; iVar < nVar; iVar++) val_residual[iVar] = 0.0;


	/*--- Solve for momentum exchange between species due to collisions ---*/	
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);			
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
			if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
			else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);			
			if (iSpecies != jSpecies) {
				/*				collisionFreq = (U_i[iLoc+0]/Molecular_Mass[iSpecies]) * (U_i[jLoc+0]/Molecular_Mass[jSpecies]) * PI_NUMBER
				 * pow((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0,2.0)
				 * sqrt( 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies])
															 +8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]) );*/
				collisionFreq = (U_i[jLoc+0]/Molecular_Mass[jSpecies]) * PI_NUMBER
						* pow((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0,2.0)
				* sqrt( 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies])
						+8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]) );
				for (iDim = 0; iDim < nDim; iDim++) {
					P[iSpecies][iDim] +=  U_i[iLoc+0]/Molecular_Mass[iSpecies] * 
							(Molecular_Mass[iSpecies]*Molecular_Mass[jSpecies])/(Molecular_Mass[iSpecies] + Molecular_Mass[jSpecies])
							* collisionFreq * (U_i[jLoc+1+iDim]/U_i[jLoc+0] - U_i[iLoc+1+iDim]/U_i[iLoc+0]);
				}
			}
		}
		for (iDim = 0; iDim < nDim; iDim++)
			val_residual[iLoc+1+iDim] = P[iSpecies][iDim] * Volume;
	}
}

void CSourcePieceWise_PlasmaDiatomic::SetJacobian_EnergyExch(double **val_Jacobian_i, CConfig *config) {

	unsigned short iVar, jVar;

	/*--- In a near future AD should be used instead of FD ---*/
	double *Residual_Baseline, *Residual_New, *U_Baseline, *Residual_Baseline_ElecForce, *Residual_New_ElecForce;
	double FDEpsilon = 1E-8;

	Residual_Baseline = new double[nVar];
	Residual_New = new double[nVar];
	U_Baseline = new double [nVar];
	Residual_Baseline_ElecForce = new double[nVar];
	Residual_New_ElecForce = new double[nVar];

	/*--- Note that U_i is a pointer to the main solution, a change will affect the solution ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		U_Baseline[iVar] = U_i[iVar];
		Residual_Baseline[iVar] = 0.0;
		Residual_Baseline_ElecForce[iVar] = 0.0;
	}

	/*--- Compute the baseline line residual ---*/
	SetResidual_ElecForce(Residual_Baseline_ElecForce, config);
	SetResidual_EnergyExch(Residual_Baseline, Residual_Baseline_ElecForce, config);

	/*--- Compute forward finite differences ---*/
	for (iVar = 0; iVar < nVar; iVar++) {

		/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
		U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
		for (jVar = 0; jVar < nVar; jVar++) {
			Residual_New[jVar] = 0.0;
			Residual_New_ElecForce[jVar] = 0.0;
		}
		SetResidual_ElecForce(Residual_New_ElecForce, config);
		SetResidual_EnergyExch(Residual_New, Residual_New_ElecForce, config);

		/*--- Undo the change in U_i to keep constant the solution ---*/
		U_i[iVar] = U_i[iVar] - FDEpsilon;

		/*--- Save the new Jacobian (per column) ---*/
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;			
		}
	}
	delete [] Residual_Baseline; delete [] Residual_New; delete [] U_Baseline; delete [] Residual_Baseline_ElecForce; delete [] Residual_New_ElecForce;
}


void CSourcePieceWise_PlasmaDiatomic::SetResidual_EnergyExch(double *val_residual, double *val_residual_ElecForce, CConfig *config) {
	unsigned short iDim, iSpecies, jSpecies, iLoc, jLoc;
	double vel_dot_prod, collisionFreq, collisionXSection;	

	/*--- Energy transport via electrostatic work ---*/
	/*	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_residual[loc+1+nDim] = val_residual[loc+1+nDim] + (val_residual_ElecForce[loc+1+iDim]/Volume) * (U_i[loc+1+iDim]/U_i[loc]) * Volume;
		}																										
	}*/

	/*--- Energy transport via elastic collisions ---*/
	//Comment: From Lee 1985, and originally from Engineering Magnetohydrodynamics by Sutton (1965)
	/*Comment: Two terms, both from kinetic theory.  The first term accounts for changes in species temperature from 
						 elastic encounters between particles, the second term accounts for frictional heating between species*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) Q_elastic[iSpecies] = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){
			if (jSpecies != iSpecies) {
				if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
				else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				/* Original
				collisionFreq = (U_i[iLoc+0]/Molecular_Mass[iSpecies]) * (U_i[jLoc+0]/Molecular_Mass[jSpecies]) * PI_NUMBER
				 * (Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0
				 * sqrt( 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies])
														   +8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]) );*/				

				/* Reduced mass
				collisionFreq = (U_i[iLoc+0]/Molecular_Mass[iSpecies]) * (U_i[jLoc+0]/Molecular_Mass[jSpecies]) * PI_NUMBER
				 * pow((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0,2.0)
				 * sqrt( 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies])
																+8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]) );*/

				/* Varying collision cross section
				if (config->GetParticle_ChargeNumber(iSpecies) < 0 && config->GetParticle_ChargeNumber(jSpecies) == 0) {
					collisionXSection = 1.0E-20;
				} else if (config->GetParticle_ChargeNumber(iSpecies) < 0 && config->GetParticle_ChargeNumber(jSpecies) > 0) {
					cutoffParameter = sqrt(BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies] / (4.0*PI_NUMBER*(U_i[iLoc+0]/Molecular_Mass[iSpecies])*pow(ELECTRON_CHARGE,2.0)));
					collisionXSection = 4.0/3.0*2.0*PI_NUMBER*pow(ELECTRON_CHARGE,4.0)/pow(3.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies],2.0)
				 * log(1.0 + pow(cutoffParameter*3.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies],2.0)/pow(ELECTRON_CHARGE,4.0));
				} else if (config->GetParticle_ChargeNumber(iSpecies) > 0 && config->GetParticle_ChargeNumber(jSpecies) < 0) {
					cutoffParameter = sqrt(BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies] / (4.0*PI_NUMBER*(U_i[iLoc+0]/Molecular_Mass[jSpecies])*pow(ELECTRON_CHARGE,2.0)));
					collisionXSection = 4.0/3.0*2.0*PI_NUMBER*pow(ELECTRON_CHARGE,4.0)/pow(3.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies],2.0)
				 * log(1.0 + pow(cutoffParameter*3.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies],2.0)/pow(ELECTRON_CHARGE,4.0));
				} else {
					collisionXSection = PI_NUMBER * pow((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0,2.0);
				}*/


				collisionXSection = PI_NUMBER * pow((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0,2.0);
				collisionFreq = (U_i[jLoc+0]/Molecular_Mass[jSpecies]) * collisionXSection
						* sqrt( 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies])
								+8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]) );
				vel_dot_prod = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					vel_dot_prod += (U_i[iLoc+1+iDim] - U_i[jLoc+1+iDim])*U_i[iLoc+1+iDim];
				}

				/* Exchange from Lee and Sutton, heavy particles */
				Q_elastic[iSpecies] +=   (U_i[iLoc+0]*AVOGAD_CONSTANT/Molar_Mass[iSpecies]) * 2.0 
						* (Molecular_Mass[iSpecies]*Molecular_Mass[jSpecies])/(Molecular_Mass[iSpecies] + Molecular_Mass[jSpecies])
						/ Molecular_Mass[jSpecies] * (3.0/2.0*BOLTZMANN_CONSTANT*(Temp_tr_i[jSpecies] - Temp_tr_i[iSpecies]))
						+ (U_i[iLoc]*AVOGAD_CONSTANT/Molar_Mass[iSpecies])
						* (Molecular_Mass[iSpecies]*Molecular_Mass[jSpecies])/(Molecular_Mass[iSpecies] + Molecular_Mass[jSpecies])
						* collisionFreq * vel_dot_prod;

				/* Exchange from Lee and Sutton, electron simplification
				Q_elastic[iSpecies] += 2.0*U_i[iLoc+0]*collisionFreq/Molecular_Mass[jSpecies] 
				 * (3.0/2.0*BOLTZMANN_CONSTANT*(Temp_tr_i[jSpecies] - Temp_tr_i[iSpecies]));
															 + U_i[iLoc]*collisionFreq * vel_dot_prod; */
			}									
		}
		val_residual[iLoc+1+nDim] =  Q_elastic[iSpecies] * Volume;
	}

	/*--- Energy exchange: Tranlsation-vibration (Landau-Teller Form)---*/
	//Comment: L-T form good for low vib. temperature fields.  Start with this, but will need to improve
	/*	MixtureNumberDensity = 0.0;
	MixturePressure = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	
		MixtureNumberDensity += U_i[loc+0] / Molecular_Mass[iSpecies];
		MixturePressure += Pressure[iSpecies];
	}
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) {
			iLoc = (nDim+3)*iSpecies;
		//else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	
			E1 = CharVibTemp[iSpecies] * BOLTZMANN_CONSTANT;
			equilVibEnergy = E1 / (exp(E1/(BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies])) - 1.0);
			E_vib = U_i[iLoc+nDim+2];
			LTNumerator = 0.0;
			LTDenominator = 0.0;
			for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
				if (jSpecies < nDiatomics) jLoc = (nDim+3)*jSpecies;
				else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);	
				ReducedMass = Molecular_Mass[iSpecies]*Molecular_Mass[jSpecies] / (Molecular_Mass[iSpecies] + Molecular_Mass[jSpecies]);
				A_sr = 1.16*10E-3 * pow(ReducedMass,1.0/2.0) * pow(CharVibTemp[iSpecies],4.0/3.0);
				B_sr = 0.015 * pow(ReducedMass,1.0/4.0);
				RelaxationTime = 101325 / MixturePressure * exp(A_sr*(pow(T,-1.0/3.0)-B_sr)-18.42);
				LTNumerator += U_i[jLoc+0] / (Molecular_Mass[jSpecies] * MixtureNumberDensity);
				LTDenominator += (U_i[jLoc+0] / (Molecular_Mass[jSpecies] * MixtureNumberDensity)) / RelaxationTime;
			}
			LTRelaxationTime = LTNumerator / LTDenominator;
			LimitRelaxationTime = 1 / (LimitCrossSection*sqrt(8.0*UNIVERSAL_GAS_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies]))*MixtureNumberDensity);
			Q_tv[iSpecies] = U_i[iLoc+0]  * (equilVibEnergy - E_vib) / (LTRelaxationTime + LimitRelaxationTime);
		}
	} */

	/*--- Energy exchange: Translation-electron ---*/
	//Comment: Applied only to the electron gas
	//Comment: NOT FINISHED.  NEED INCLUDE TERM FOR THE COLLISION CROSS-SECTIONS!!!
	/*	double collisionXsection;
	unsigned int ELECTRON_SPECIES_NUM, loc_electron;
	ELECTRON_SPECIES_NUM = 7;
	loc_electron = (nDim+3)*nDiatomics + (nDim+2)*(ELECTRON_SPECIES_NUM-nDiatomics);
	Q_te = 3 * UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies] * U_i[loc_electron+0]
	 * (Temp_tr_i[iSpecies] - Temp_tr_i[ELECTRON_SPECIES_NUM])
	 * sqrt(8.0*(UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies]*Temp_tr_i[ELECTRON_SPECIES_NUM] / (PI_NUMBER * Molar_Mass[ELECTRON_SPECIES_NUM])));
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {			
		if (iSpecies != ELECTRON_SPECIES_NUM) {
			if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			if (ChargeNumber[iSpecies] == 0)
				collisionXsection = 1E-20;
			else if (ChargeNumber[iSpecies] > 0)
//				collisionXsection = 8*PI_NUMBER/27 * 
			Q_te = Q_te + U_i[loc+0]*AVOGAD_CONSTANT/(Molecular_Mass[iSpecies]*Molecular_Mass[iSpecies]) * collisionXsection;
		}												
	}
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		val_residual[loc+1+nDim] = 0.0; //Total energy residual
		val_residual[loc+2+nDim] = 0.0; //Vibrational energy residual
	} */
}

CSource_Template::CSource_Template(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSource_Template::~CSource_Template(void) { }

void CSource_Template::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) { }
