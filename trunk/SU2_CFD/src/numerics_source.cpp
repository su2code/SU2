/*!
 * \file numerics_source.cpp
 * \brief This file contains all the source term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.
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
	
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible = config->GetIncompressible();

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

//************************************************//
// Please do not delete //SU2_CPP2C comment lines //
//************************************************//

//SU2_CPP2C START CSourcePieceWise_TurbSA::SetResidual
//SU2_CPP2C CALL_LIST START
//SU2_CPP2C INVARS *U_i **PrimVar_Grad_i Laminar_Viscosity_i *TurbVar_i
//SU2_CPP2C OUTVARS *val_residual
//SU2_CPP2C VARS DOUBLE dist_i cv1_3 k2 cb1 Volume cw2 cw3_6 cb2 sigma
//SU2_CPP2C CALL_LIST END

//SU2_CPP2C DEFINE nDim

//SU2_CPP2C DECL_LIST START
//SU2_CPP2C VARS INT SCALAR iDim
//SU2_CPP2C VARS DOUBLE SCALAR Density_i DivVelocity Vorticity dist_0_2
//SU2_CPP2C VARS DOUBLE SCALAR nu Ji Ji_2 Ji_3 fv1 fv2 Omega Shat
//SU2_CPP2C VARS DOUBLE SCALAR r g g_6 glim fw norm2_Grad
//SU2_CPP2C DECL_LIST END

//SU2_CPP2C COMMENT START
	if (incompressible) Density_i = DensityInc_i;
	else {
//SU2_CPP2C COMMENT END
		Density_i = U_i[0];
//SU2_CPP2C COMMENT START
	}
//SU2_CPP2C COMMENT END
	
	val_residual[0] = 0.0;
//SU2_CPP2C COMMENT START
	if (implicit)
		val_Jacobian_i[0][0] = 0.0;
//SU2_CPP2C COMMENT END

	/*--- Computation of divergence of velocity and vorticity ---*/
	DivVelocity = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		DivVelocity += PrimVar_Grad_i[iDim+1][iDim];

	Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
	if (nDim == 3)
		Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) +
				(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );
	Vorticity = sqrt(Vorticity);

//SU2_CPP2C COMMENT START
	switch (config->GetKind_Turb_Model()) {
	case SA :
//SU2_CPP2C COMMENT END
		if (dist_i > 0.0) {

			/*--- Production term ---*/
			dist_0_2 = dist_i*dist_i;
			nu = Laminar_Viscosity_i/Density_i;
			Ji = TurbVar_i[0]/nu;
			Ji_2 = Ji*Ji;
			Ji_3 = Ji_2*Ji;
			fv1 = Ji_3/(Ji_3+cv1_3);
			fv2 = 1.0 - Ji/(1.0+Ji*fv1);
			Omega = Vorticity;
			Shat = max(Omega + TurbVar_i[0]*fv2/(k2*dist_0_2), TURB_EPS);
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
//SU2_CPP2C COMMENT START
			/*--- Implicit part ---*/
			if (implicit) {

				/*--- Production term ---*/
				dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
				dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
				if ( Shat <= TURB_EPS )
					dShat = 0.0;
				else
					dShat = (fv2+TurbVar_i[0]*dfv2)/(k2*dist_0_2);
//				val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;

				/*--- Destruction term ---*/
				dr = (Shat-TurbVar_i[0]*dShat)/(Shat*Shat*k2*dist_0_2);
				if (r == 10.) dr = 0.0;
				dg = dr*(1.+cw2*(6.*pow(r,5.)-1.));
				dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
//				val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] +	2.*fw)*TurbVar_i[0]/dist_0_2*Volume;
			}
//SU2_CPP2C COMMENT END
		}
//SU2_CPP2C COMMENT START
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
			Shat = max(Omega + TurbVar_i[0]*fv2/(Density_i*k2*dist_0_2),TURB_EPS);
			val_residual[0] += cb1*Shat*TurbVar_i[0]*Volume;

			/*--- Destruction term ---*/
			r = min(TurbVar_i[0]/(Density_i*Shat*k2*dist_0_2),10.);
			g = r + cw2*(pow(r,6.)-r);
			g_6 = pow(g,6.);
			glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
			fw = g*glim;
			val_residual[0] -= cw1*fw*TurbVar_i[0]*TurbVar_i[0]/(Density_i*dist_0_2)*Volume;

			/*--- Diffusion term ---*/
			nu_hat_i = TurbVar_i[0]/Density_i;
			norm2_Grad = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				norm2_Grad += pow(TurbVar_Grad_i[0][iDim]-ConsVar_Grad_i[0][iDim]*nu_hat_i,2.0);
			val_residual[0] += cb2/(sigma*Density_i)*norm2_Grad*Volume;

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
				dr = (Shat-TurbVar_i[0]*dShat)/(Shat*Shat*Density_i*k2*dist_0_2);
				if (r == 10.) dr = 0.0;
				dg = dr*(1.+cw2*(6.*pow(r,5.)-1.));
				dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
				val_Jacobian_i[0][0] -= cw1/Density_i*(dfw*TurbVar_i[0] + 2.*fw)*TurbVar_i[0]/dist_0_2*Volume;

				/*--- Diffusion term ---*/
				prod_grads = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					grad_nu_hat = TurbVar_Grad_i[0][iDim]-ConsVar_Grad_i[0][iDim]*nu_hat_i;
					prod_grads += grad_nu_hat*ConsVar_Grad_i[0][iDim]/Density_i;
				}
				val_Jacobian_i[0][0] -= 2.0*cb2/(sigma*Density_i)*prod_grads*Volume;
			}
		}
		break;
	}
//SU2_CPP2C COMMENT END

//SU2_CPP2C END CSourcePieceWise_TurbSA::SetResidual

}

CSourcePieceWise_TransLM::CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Spalart-Allmaras Closure constants ---*/
	cv1_3 = pow(7.1,3.0);
	k2 = pow(0.41,2.0);
	cb1 = 0.1355;
	cw2 = 0.3;
	cw3_6 = pow(2.0,6.0);
	sigma = 2./3.;
	cb2 = 0.622;
	cw1 = cb1/k2+(1+cb2)/sigma;

	/*-- gamma-theta closure constants --*/
	c_e1    = 1.0;
	c_a1    = 2.0;
	c_e2    = 50.0;
	c_a2    = 0.06;
	sigmaf  = 1.0;
	s1      = 2.0;
	c_theta = 0.03;
	sigmat  = 2.0;

	/*-- Correlation constants --*/
	flen_global  = 12.0;
	alpha_global = 0.85;
}

CSourcePieceWise_TransLM::~CSourcePieceWise_TransLM(void) { }

void CSourcePieceWise_TransLM::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	// TODO, Aniket

	/*-- Local intermediate variables --*/
	double rey_tc, flen, re_v, strain, f_onset1,f_onset2,f_onset3,f_onset,f_turb,tu;

	double prod, des;
	double f_lambda, re_theta, rey, re_theta_lim;
	double Velocity_Mag = 0.0, du_ds, theta, lambda, time_scale, delta_bl, delta, f_wake, var1, f_theta;
	double theta_bl;

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

	/*-- Note: no imcompressible for now! --*/

	if (dist_i > 0.0) {

		/*-- Spalart-Allmaras variable eq: --*/
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

		/*-- Intermittency eq.: --*/
		rey_tc = alpha_global*TurbVar_i[2];  // REth critical correlation
		flen   = flen_global/TurbVar_i[2];   // f_length correlation
		strain = Vorticity*sqrt(2.0);
		re_v   = U_i[0]*dist_0_2/Laminar_Viscosity_i*strain;  // Vorticity Reynolds number

		f_onset1 = re_v / (2.193*rey_tc);
		f_onset2 = min(max(f_onset1, pow(f_onset1,4)), 2.);
		f_onset3 = max(1. - pow(0.4*Ji,3),0.);
		f_onset  = max(f_onset2 - f_onset3, 0.);

		f_turb = exp(-pow(0.25*Ji,4));

		prod = flen*c_a1*strain*sqrt(f_onset*TurbVar_i[1]);
		prod = prod*(1. - c_e1*TurbVar_i[1]);

		des = c_a2*TurbVar_i[1]*Vorticity*f_turb;
		des = des*(c_e2*TurbVar_i[1] - 1.);

		// TODO: Find out what bjmat is all about in Vinod's implementation
		val_residual[1] = val_residual[1] + prod - des;

		/*--REtheta eq: --*/
		tu = config->GetTurbulenceIntensity_FreeStream();
		f_lambda = 1.;

		if (nDim==2) {
			Velocity_Mag = sqrt(U_i[1]*U_i[1]+U_i[2]*U_i[2])/U_i[0];
		} else if (nDim==3) {
			Velocity_Mag = sqrt(U_i[1]*U_i[1]+U_i[2]*U_i[2]+U_i[3]*U_i[3])/U_i[0];
		}

		du_ds = U_i[1]/(U_i[0]*Velocity_Mag) * PrimVar_Grad_i[1][0] +  // Streamwise velocity derivative
				U_i[2]/(U_i[0]*Velocity_Mag) * PrimVar_Grad_i[2][1];
		if (nDim==3)
			du_ds += U_i[3]/(U_i[0]*Velocity_Mag) * PrimVar_Grad_i[3][2];

		re_theta_lim = 20.;
		rey = config->GetReynolds();

		/*-- subiterations to solve REth correlation --*/
		for (int iter=0; iter<10; iter++) {
			if (tu <= 1.3) {
				re_theta = f_lambda * (1173.51-589.428*tu+0.2196/(tu*tu));
			} else {
				re_theta = 331.5 * f_lambda*pow(tu-0.5658,-0.671);
			}
			re_theta = max(re_theta/rey, re_theta_lim);

			theta  = re_theta * Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);

			lambda = U_i[0]*theta*theta/Laminar_Viscosity_i * du_ds;
			lambda = min(max(-0.1,lambda),0.1);

			if (lambda<=0.0) {
				f_lambda = 1. - (-12.986*lambda - 123.66*lambda*lambda -
						405.689*lambda*lambda*lambda)*exp(-pow(2./3*tu,1.5));
			} else {
				f_lambda = 1. + 0.275*(1.-exp(-35.*lambda))*exp(-2.*tu);
			}

			/*-- Calculate blending function f_theta --*/
			time_scale = 500.0*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag*Velocity_Mag);
			theta_bl   = TurbVar_i[2]*Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);
			delta_bl   = 7.5*theta_bl;
			delta      = 50*Vorticity*dist_i/Velocity_Mag*delta_bl;

			f_wake = 1.;

			var1 = (TurbVar_i[1]-1./c_e2)/(1.0-1./c_e2);
			var1 = 1. - pow(var1,2);
			f_theta = min(max(f_wake*exp(-pow(dist_i/delta,4)), var1),1.0);

			val_residual[2] = c_theta*U_i[0]/time_scale *  (1.-f_theta) * (re_theta-TurbVar_i[2]);

			// TODO: Left off here!
		}

		/*--- Implicit part ---*/
		if (implicit) {
			cout << "Implicit SAGT not implemented yet!!" << endl;
			int j;
			cin >> j;

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
}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double *constants,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible = config->GetIncompressible();

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Closure constants ---*/
	beta_star     = constants[6];
	sigma_omega_1 = constants[2];
	sigma_omega_2 = constants[3];
	beta_1        = constants[4];
	beta_2        = constants[5];
	alfa_1        = constants[8];
	alfa_2        = constants[9];
	a1            = constants[7];
}

CSourcePieceWise_TurbSST::~CSourcePieceWise_TurbSST(void) { }

void CSourcePieceWise_TurbSST::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

//************************************************//
// Please do not delete //SU2_CPP2C comment lines //
//************************************************//

//SU2_CPP2C START CSourcePieceWise_TurbSST::SetResidual
//SU2_CPP2C CALL_LIST START
//SU2_CPP2C INVARS *val_U_i
//SU2_CPP2C OUTVARS *val_laminar_viscosity_i
//SU2_CPP2C VARS DOUBLE Temperature_Ref Viscosity_Ref, Gamma_Minus_One
//SU2_CPP2C CALL_LIST END

//SU2_CPP2C DEFINE nDim NONE SA SST

//SU2_CPP2C DECL_LIST START
//SU2_CPP2C DECL_LIST END

	unsigned short iDim;
	double alfa_blended, beta_blended;
	double diverg, pk, pw, zeta;

	val_residual[0] = 0.0;
	val_residual[1] = 0.0;
//SU2_CPP2C COMMENT START
	if (implicit){
		val_Jacobian_i[0][0] = 0.0;		val_Jacobian_i[0][1] = 0.0;
		val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = 0.0;
	}
//SU2_CPP2C COMMENT END

	/*--- Computation of blended constants for the source terms---*/
	alfa_blended = F1_i*alfa_1 + (1.0 - F1_i)*alfa_2;
	beta_blended = F1_i*beta_1 + (1.0 - F1_i)*beta_2;

	if (dist_i > 0.0) {
		/*--- Production ---*/
		diverg = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			diverg += PrimVar_Grad_i[iDim+1][iDim];

		pk = Eddy_Viscosity_i*StrainMag*StrainMag - 2.0/3.0*U_i[0]*TurbVar_i[0]*diverg;
		pk = min(pk,20.0*beta_star*U_i[0]*TurbVar_i[1]*TurbVar_i[0]);
		pk = max(pk,0.0);

		zeta = max(TurbVar_i[1],StrainMag*F2_i/a1);
		pw = StrainMag*StrainMag - 2.0/3.0*zeta*diverg;
		pw = max(pw,0.0);

		val_residual[0] += pk*Volume;
		val_residual[1] += alfa_blended*U_i[0]*pw*Volume;

		/*--- Dissipation ---*/
		val_residual[0] -= beta_star*U_i[0]*TurbVar_i[1]*TurbVar_i[0]*Volume;
		val_residual[1] -= beta_blended*U_i[0]*TurbVar_i[1]*TurbVar_i[1]*Volume;

		/*--- Cross diffusion ---*/
		val_residual[1] += (1.0 - F1_i)*CDkw*Volume;

//SU2_CPP2C COMMENT START
		/*--- Implicit part ---*/
		if (implicit) {
			val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;		val_Jacobian_i[0][1] = 0.0;
			val_Jacobian_i[1][0] = 0.0;									val_Jacobian_i[1][1] = -2.0*beta_blended*TurbVar_i[1]*Volume;
		}
//SU2_CPP2C COMMENT END
	}

//SU2_CPP2C END CSourcePieceWise_TurbSST::SetResidual

}

CSourcePieceWise_FreeSurface::CSourcePieceWise_FreeSurface(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

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
	Ni_times_Nj = new double*[4];
	for (unsigned short iVar = 0; iVar < 4; iVar++) {
		Ni_times_Nj[iVar] = new double [4];
	}

	Ni_times_Nj[0][0] = 2.0/120*6;	Ni_times_Nj[0][1] = 1.0/120*6;	Ni_times_Nj[0][2] = 1.0/120*6;	Ni_times_Nj[0][3] = 1.0/120*6;
	Ni_times_Nj[1][0] = 1.0/120*6;	Ni_times_Nj[1][1] = 2.0/120*6;	Ni_times_Nj[1][2] = 1.0/120*6;	Ni_times_Nj[1][3] = 1.0/120*6;
	Ni_times_Nj[2][0] = 1.0/120*6;	Ni_times_Nj[2][1] = 1.0/120*6;	Ni_times_Nj[2][2] = 2.0/120*6;	Ni_times_Nj[2][3] = 1.0/120*6;
	Ni_times_Nj[3][0] = 1.0/120*6;	Ni_times_Nj[3][1] = 1.0/120*6;	Ni_times_Nj[3][2] = 1.0/120*6;	Ni_times_Nj[3][3] = 2.0/120*6;

}

CSourcePieceWise_Elec::~CSourcePieceWise_Elec(void) { }

void CSourcePieceWise_Elec::SetResidual_MacCormack(double *val_residual, CConfig *config) {

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

			rho_Pos  = 1.0/3.0*(U_0[0] + U_1[0] + U_2[0]) ;
			rho_Elec = 1.0/3.0*(U_0[1] + U_1[1] + U_2[1]) ;

			/*--- Source q ---*/
			mass_Elec = config->GetParticle_Mass(2);
			mass_Pos  = config->GetParticle_Mass(1);

			ne = rho_Elec / mass_Elec;
			np = rho_Pos  / mass_Pos;

			Kai_n = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne - np );
			alpha = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Elec / (mass_Elec*mass_Elec) + rho_Pos / (mass_Pos*mass_Pos));

			diff_ru_Pos_x  = 1.0/3.0*(ConsVar_Grad_0[1][0] + ConsVar_Grad_1[1][0] + ConsVar_Grad_2[1][0]);
			diff_ru_Elec_x = 1.0/3.0*(ConsVar_Grad_0[2][0] + ConsVar_Grad_1[2][0] + ConsVar_Grad_2[2][0]);
			diff_ru_Pos_y  = 1.0/3.0*(ConsVar_Grad_0[1][1] + ConsVar_Grad_1[1][1] + ConsVar_Grad_2[1][1]);
			diff_ru_Elec_y = 1.0/3.0*(ConsVar_Grad_0[2][1] + ConsVar_Grad_1[2][1] + ConsVar_Grad_2[2][1]);


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
	if (config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == O2 || config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {
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

void CSourcePieceWise_Elec::SetResidual(double *val_residual, CConfig *config) {

	if (config->GetKind_GasModel() == ARGON) {
		double Kai_n, ne, np;
		double rho_Pos = 0.0, rho_Elec = 0.0, mass_Pos;
		double a[4], b[4], Area = 0.0, f[4];
		mass_Pos  = config->GetParticle_Mass(1);
		if (nDim == 2) {
			for (unsigned short iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);			/*--- Norm of the normal component of area, area = 1/2*cross(a,b) ---*/

			rho_Pos  = 1.0/3.0*(U_0[0] + U_1[0] + U_2[0]) ;
			rho_Elec = 1.0/3.0*(U_0[1] + U_1[1] + U_2[1]) ;

			/*--- Source q ---*/

			ne = rho_Elec / ELECTRON_MASS;
			np = rho_Pos  / mass_Pos;

			Kai_n = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne - np );

			/*--- Residual = transpose(N) * source (q) * Area  ---*/
			val_residual[0] =-1.0/3.0 * Kai_n * Area;
			val_residual[1] =-1.0/3.0 * Kai_n * Area;
			val_residual[2] =-1.0/3.0 * Kai_n * Area;
		}
		if (nDim == 3) {

			f[0] = ELECTRON_CHARGE/FREE_PERMITTIVITY * ( U_0[1]/ELECTRON_MASS - U_0[0]/mass_Pos );
			f[1] = ELECTRON_CHARGE/FREE_PERMITTIVITY * ( U_1[1]/ELECTRON_MASS - U_1[0]/mass_Pos );
			f[2] = ELECTRON_CHARGE/FREE_PERMITTIVITY * ( U_2[1]/ELECTRON_MASS - U_2[0]/mass_Pos );
			f[3] = ELECTRON_CHARGE/FREE_PERMITTIVITY * ( U_3[1]/ELECTRON_MASS - U_3[0]/mass_Pos );

			/*--- Residual = transpose(N) * source (q) * Area  ---*/
			for (unsigned short iVar = 0; iVar < 4; iVar ++ ) {
				val_residual[iVar] = 0.0;
				for (unsigned short jVar = 0; jVar < 4; jVar ++ )
					val_residual[iVar] -= Ni_times_Nj[iVar][jVar] * f[jVar] * Volume;
			}
		}
	}
}

CSourceViscous_AdjFlow::CSourceViscous_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
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

CSourceViscous_AdjFlow::~CSourceViscous_AdjFlow(void) {
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

void CSourceViscous_AdjFlow::SetResidual (double *val_residual, CConfig *config) {

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

}

CSourcePieceWise_AdjDiscFlow::CSourcePieceWise_AdjDiscFlow() {

}

CSourcePieceWise_AdjDiscFlow::~CSourcePieceWise_AdjDiscFlow(void) {

}

void CSourcePieceWise_AdjDiscFlow::SetResidual () {

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

CSourcePieceWise_AdjDiscTurb::CSourcePieceWise_AdjDiscTurb() {

}

CSourcePieceWise_AdjDiscTurb::~CSourcePieceWise_AdjDiscTurb(void) {

}

void CSourcePieceWise_AdjDiscTurb::SetResidual() {

}

CSourcePieceWise_AdjElec::CSourcePieceWise_AdjElec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { 

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

CSourcePieceWise_AdjElec::~CSourcePieceWise_AdjElec(void) { }

void CSourcePieceWise_AdjElec::SetResidual(double *val_residual, CConfig *config) {
	val_residual[0] = Volume*sin(PI_NUMBER*Coord_i[0])*sin(PI_NUMBER*Coord_i[1]);
}

CSourcePieceWise_AdjDiscElec::CSourcePieceWise_AdjDiscElec() {

}

CSourcePieceWise_AdjDiscElec::~CSourcePieceWise_AdjDiscElec(void) { }

void CSourcePieceWise_AdjDiscElec::SetResidual() {

}

CSourcePieceWise_AdjLevelSet::CSourcePieceWise_AdjLevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { 
}

CSourcePieceWise_AdjLevelSet::~CSourcePieceWise_AdjLevelSet(void) { }

void CSourcePieceWise_AdjLevelSet::SetResidual(double *val_residual, CConfig *config) {}

CSourcePieceWise_AdjDiscLevelSet::CSourcePieceWise_AdjDiscLevelSet() {
}

CSourcePieceWise_AdjDiscLevelSet::~CSourcePieceWise_AdjDiscLevelSet(void) { }

void CSourcePieceWise_AdjDiscLevelSet::SetResidual() {}

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

CSourceConservative_AdjDiscFlow::CSourceConservative_AdjDiscFlow() {

}

CSourceConservative_AdjDiscFlow::~CSourceConservative_AdjDiscFlow(void) {

}

void CSourceConservative_AdjDiscFlow::SetResidual () {

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

CSourceConservative_AdjDiscTurb::CSourceConservative_AdjDiscTurb() {

}

CSourceConservative_AdjDiscTurb::~CSourceConservative_AdjDiscTurb(void) {
}

void CSourceConservative_AdjDiscTurb::SetResidual() {

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

CSourceRotationalFrame_AdjDiscFlow::CSourceRotationalFrame_AdjDiscFlow() { }

CSourceRotationalFrame_AdjDiscFlow::~CSourceRotationalFrame_AdjDiscFlow(void) { }

void CSourceRotationalFrame_AdjDiscFlow::SetResidual() {

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

	unsigned short iReaction;

	U_id = new double [nVar];
  Temp_tr_id = new double[nSpecies];

	nSpecies 		= config->GetnSpecies();
	nReactions 		= config->GetnReactions();
	GammaMonatomic 	= config->GetGammaMonatomic();
	GammaDiatomic 	= config->GetGammaDiatomic();
	nMonatomics 	= config->GetnMonatomics();
	nDiatomics 		= config->GetnDiatomics();

	Molar_Mass 		= new double[nSpecies];
	ChargeNumber 	= new double[nSpecies];
	w_dot 			= new double[nSpecies];					w_dotd = new double[nSpecies];
	Q_tv 			= new double[nSpecies];
	Q_elastic 		= new double [nSpecies];			Q_elasticd	  = new double [nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Molar_Mass[iSpecies] = config->GetMolar_Mass(iSpecies);
		ChargeNumber[iSpecies] = config->GetParticle_ChargeNumber(iSpecies);
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
			config->GetKind_GasModel() == AIR5 || config->GetKind_GasModel() == AIR7) { 
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
	/*	M1 = config->GetParticle_Mass(0);
		M3 = config->GetParticle_Mass(2);
		M2 = config->GetParticle_Mass(1);
	 */

	SourceVector = new double[nVar];
	SourceJacobian = new double*[nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar ++) {
		SourceVector[iVar] = 10.0;
		SourceJacobian[iVar] = new double [nVar];
		for (unsigned short jVar = 0; jVar < nVar; jVar ++)
			SourceJacobian[iVar][jVar] = 0.0;
	}


	MagneticField = new double [3];
	VcrossB = new double* [nSpecies];
	velocity = new double *[nSpecies];

	/* For the more accurate magnetic field model */
	VioncrossB = new double [nDim];
	Current_Density = new double [3];
	JcrossB = new double [3];
	Electric_Conductivity = 1250.0;

	for (iVar = 0; iVar < 3; iVar ++) {
		Current_Density[iVar] = 0.0;
		JcrossB[iVar] = 0.0;
	}

	/* Till here For the more accurate magnetic field model */

	EMF  = new double* [nSpecies];
	MagneticDipole = new double [3];
	ElectricField = new double [nDim];

	for (iVar = 0; iVar < 3; iVar ++) {
		MagneticField[iVar] = 0.0;
		MagneticDipole[iVar] = 0.0;
	}
	MagneticDipole[0] = config->GetMagneticDipole(0);;
	MagneticDipole[1] = config->GetMagneticDipole(1);
	MagneticDipole[2] = config->GetMagneticDipole(2);

	for (iDim = 0; iDim < nDim; iDim++)
		ElectricField[iDim] = 1.0;
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
	M1Avg = config->GetMolar_Mass(0);
	M2Avg = config->GetMolar_Mass(1);
	M3Avg = config->GetMolar_Mass(2);
	M1M2M3Avg3 = M1Avg * M2Avg* M3Avg;

	M1 = M1Avg/AVOGAD_CONSTANT;
	M3 = M3Avg/AVOGAD_CONSTANT;
	M2 = M1-M3;

	M1 = config->GetParticle_Mass(0);
	M3 = config->GetParticle_Mass(2);
	M2 = config->GetParticle_Mass(1);

	MagneticInteraction = config->GetMagneticInteractionParameter();
	ReferenceLength = config->GetReferenceLength();
	ReferenceDensity  = config->GetReferenceDensity();
	ReferenceSpeed = config->GetReferenceSpeed();

	Sigma_Bsqr = MagneticInteraction * ReferenceDensity * ReferenceSpeed / ReferenceLength;

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

	delete [] ElectricField;	delete[] MagneticField;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		delete [] velocity[iSpecies];
		delete [] EMF[iSpecies];
		delete [] VcrossB[iSpecies];
	}
	delete [] velocity;	delete [] EMF;	delete [] VcrossB;
	delete [] VioncrossB;	delete [] Current_Density; delete [] JcrossB;
	delete [] MagneticDipole;

}

void CSourcePieceWise_Plasma::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
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

	/*	double Mconstant = .15;
	double rB = 0.01825;

	for (iVar = 0; iVar < 3; iDim ++ )
		MagneticField[iVar] = Mconstant * pow(( rB/distance),3.3);
	 */

	double Bx = MagneticField[0];
	double By = MagneticField[1];
	double Bz = MagneticField[2];
	double s0 = Electric_Conductivity;


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
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ )
			EMF[iSpecies][iDim] = (ElectricField[iDim] + VcrossB[iSpecies][iDim]);
	}

	for (iDim = 0; iDim < nDim; iDim ++ )
		Current_Density[iDim] = Electric_Conductivity*(ElectricField[iDim] + VcrossB[1][iDim]);

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

	double delJcrossBdelm = s0/r2/r2 * ( Bx * (l2*Bz +   n2*By) + By * (n2*Bx - 2*m2*By) + Bz * (l2*Bx - 2*m2*Bz) );
	double delJcrossBdeln = s0/r2/r2 * ( Bx * (m2*By - 2*n2*Bx) + By * (m2*Bx +   l2*Bz) + Bz * (l2*By - 2*n2*Bz) );
	double delJcrossBdell = s0/r2/r2 * ( Bx * (m2*Bz - 2*l2*Bx) + By * (n2*Bz - 2*l2*By) + Bz * (m2*Bx +   n2*By) );

	double delJdotJdelm   =  1/s0*( Current_Density[0] * 0  	  - Bz*s0/r2*Current_Density[1] + By*s0/r2*Current_Density[2] );
	double delJdotJdeln   =  1/s0*( Current_Density[0] * Bz*s0/r2 - 0*Current_Density[0] 	    - Bx*s0/r2*Current_Density[2] );
	double delJdotJdell   =  1/s0*(-Current_Density[0] * By*s0/r2 + Bx*s0/r2*Current_Density[1] + 0*Current_Density[2] );

	/* Derivative of the first term in the source terms wrt all the conservative variables */
	unsigned short loc;

	loc = 0;
	SourceVector[loc+0] = M1*AvgNum*R;

	SourceJacobian[loc+0][loc+0] =  M1Avg*dR_dr1;
	SourceJacobian[loc+0][loc+1] = M1Avg*dR_dm1;
	SourceJacobian[loc+0][loc+2] = M1Avg*dR_dn1;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] = M1Avg*dR_dl1;

	SourceJacobian[loc+0][loc+nDim+1] = M1Avg*dR_de1;

	/* Derivative of the second term in the source terms wrt all the conservative variables */

	SourceVector[loc+1] = r1*nu12*(m2/r2 - m1/r1)+ r1*nu13*(m3/r3 - m1/r1);

	SourceJacobian[loc+1][loc+0] = nu12*m2/r2 + nu13*m3/r3;// + r1*(m2/r2-m1/r1)*dv12_dT1*dT1_dr1 + r1*(m3/r3 - m1/r1)*dv13_dT1*dT1_dr1;
	SourceJacobian[loc+1][loc+1] = -nu12 - nu13 + r1*(m2/r2-m1/r1)*dv12_dT1*dT1_dm1 + r1*(m3/r3 - m1/r1)*dv13_dT1*dT1_dm1;
	SourceJacobian[loc+1][loc+2] = 0.0;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] = 0;
	SourceJacobian[loc+1][loc+nDim+1] = r1*dv12_dT1*dT1_de1*(m2/r2 - m1/r1) + r1*dv13_dT1*dT1_de1*(m3/r3 - m1/r1);

	/* Derivative of the third term in the source terms wrt all the conservative variables */

	SourceVector[loc+2] = r1*nu12*(n2/r2 - n1/r1) + r1*nu13*(n3/r3 - n1/r1);

	SourceJacobian[loc+2][loc+0]  = nu12*n2/r2 + nu13*n3/r3;
	SourceJacobian[loc+2][loc+1]  = 0.0;
	SourceJacobian[loc+2][loc+2]  = -nu12 - nu13 + r1*(n2/r2-n1/r1)*dv12_dT1*dT1_dn1 + r1*(n3/r3 - n1/r1)*dv13_dT1*dT1_dn1;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] = 0;
	SourceJacobian[loc+2][loc+nDim+1] =  r1*(n2/r2 -n1/r1)*dv12_dT1*dT1_de1 + r1*(n3/r3 - n1/r1)*dv13_dT1*dT1_de1;


	if (nDim == 3) {
		SourceVector[loc+3] = r1*nu12*(l2/r2 - l1/r1) + r1*nu13*(l3/r3 - l1/r1);

		SourceJacobian[loc+3][loc+0]  = nu12*l2/r2 + nu13*l3/r3;
		SourceJacobian[loc+3][loc+1]  = 0.0;
		SourceJacobian[loc+3][loc+2]  = 0.0;
		SourceJacobian[loc+3][loc+3] = -nu12 - nu13 + r1*(l2/r2-l1/r1)*dv12_dT1*dT1_dl1 + r1*(l3/r3 - l1/r1)*dv13_dT1*dT1_dl1;
		SourceJacobian[loc+3][loc+nDim+1] =  r1*(l2/r2 -l1/r1)*dv12_dT1*dT1_de1 + r1*(l3/r3 - l1/r1)*dv13_dT1*dT1_de1;
	}
	/* Derivative of the Fourth term in the source terms wrt all the conservative variables */

	SourceVector[loc+nDim+1] = m1*nu12*(m2/r2 - m1/r1) + m1*nu13*(m3/r3 - m1/r1) +
			n1*nu12*(n2/r2 - n1/r1) + n1*nu13*(n3/r3 - n1/r1) +
			l1*nu12*(l2/r2 - l1/r1) + l1*nu13*(l3/r3 - l1/r1) + QT1;

	SourceJacobian[loc+nDim+1][loc+0]  = (m1*m1 + n1*n1 + l1*l1)/(r1*r1) * (nu12 + nu13);
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

	SourceVector[loc+1] = JcrossB[0] + r2*nu21*(m1/r1 - m2/r2) + r2*nu23*(m3/r3 - m2/r2);

	SourceJacobian[loc+1][loc+0] = nu21*m1/r1 + nu23*m3/r3 - JcrossB[0]/r2;
	SourceJacobian[loc+1][loc+1] = -nu21 - nu23 - s0 * ( Bz*Bz + By*By)/r2;
	SourceJacobian[loc+1][loc+2] = s0*By*Bx/r2;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] =  s0*Bx*Bz/r2;
	SourceJacobian[loc+1][loc+nDim+1]  = JcrossB[0]/e2 + r2*(m1/r1 - m2/r2)*dv21_dT2*dT2_de2 + r2*(m3/r3 - m2/r2)*dv23_dT2*dT2_de2;


	/* Derivative of the Seventh term in the source terms wrt all the conservative variables */

	SourceVector[loc+2]  =  JcrossB[1]  + r2*nu21*(n1/r1 - n2/r2) + r2*nu23*(n3/r3 - n2/r2);

	SourceJacobian[loc+2][loc+0] = nu21*n1/r1 + nu23*n3/r3 - JcrossB[1]/r2;
	SourceJacobian[loc+2][loc+1] = s0*Bx*By/r2;
	SourceJacobian[loc+2][loc+2] = -nu21 - nu23 -s0*(Bx*Bx + Bz*Bz)/r2;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] = s0*Bz*By/r2 ;
	SourceJacobian[loc+2][loc+nDim+1] = JcrossB[1]/e2 + r2*(n1/r1 - n2/r2)*dv21_dT2*dT2_de2 + r2*(n3/r3 - n2/r2)*dv23_dT2*dT2_de2;

	if (nDim == 3) {
		SourceVector[loc+3] = JcrossB[2]  + r2*nu21*(l1/r1 - l2/r2) + r2*nu23*(l3/r3 - l2/r2);

		SourceJacobian[loc+3][loc+0]  = nu21*l1/r1 + nu23*l3/r3 - JcrossB[2]/r2;
		SourceJacobian[loc+3][loc+1]  = s0*Bx*Bz/r2;
		SourceJacobian[loc+3][loc+2]  = s0*By*Bz/r2;
		SourceJacobian[loc+3][loc+3] = -nu21 - nu23 - s0*(By*By + Bx*Bx)/r2 ;
		SourceJacobian[loc+3][loc+nDim+1] = JcrossB[2]/e2 + r2*(l1/r1 -l2/r2)*dv21_dT2*dT2_de2 + r2*(l3/r3 - l2/r2)*dv23_dT2*dT2_de2;
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
	SourceJacobian[loc+1][loc+nDim+1] =  JcrossB[0]/e3 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_de3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_de3;

	/* Derivative of the Eleventh term in the source terms wrt all the conservative variables */
	SourceVector[loc+2]  = JcrossB[1]  + r3*nu31*(n1/r1 - n3/r3) + r3*nu32*(n2/r2 - n3/r3);

	SourceJacobian[loc+2][loc+0] = nu31*n1/r1 + nu32*n2/r2 - JcrossB[1]/r3;
	SourceJacobian[loc+2][loc+1] = s0*Bx*By/r3;
	SourceJacobian[loc+2][loc+2] =  -nu31 - nu32 -s0*(Bx*Bx + Bz*Bz)/r3;
	if (nDim ==3)  SourceJacobian[loc+2][loc+3] =  s0*Bz*By/r3 ;
	SourceJacobian[loc+2][loc+nDim+1] = JcrossB[1]/e3 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_de3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_de3;

	/* Derivative of the Twelfth term in the source terms wrt all the conservative variables */
	if (nDim ==3) {
		SourceVector[loc+3]  = JcrossB[2] + r3*nu31*(l1/r1 - l3/r3) + r3*nu32*(l2/r2 - l3/r3);
		SourceJacobian[loc+3][loc+0] = nu31*l1/r1 + nu32*l2/r2 - JcrossB[2]/r3;
		SourceJacobian[loc+3][loc+1] =  s0*Bx*Bz/r3;
		SourceJacobian[loc+3][loc+2] =  s0*By*Bz/r3;
		SourceJacobian[loc+3][loc+3] =  -nu31 - nu32 - s0*(By*By + Bx*Bx)/r3 ;
		SourceJacobian[loc+3][loc+nDim+1] =  JcrossB[2]/e3 + r3*(l1/r1 - l3/r3)*dv31_dT3*dT3_de3 + r3*(l2/r2 - l3/r3)*dv32_dT3*dT3_de3;
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



#ifdef EMFisEplusVcrossB
	/* Derivative of the fifth term in the source terms wrt all the conservative variables */
	loc = (nDim+2);

	SourceVector[loc+0]= mdotr;//-M2*AvgNum*R;

	SourceJacobian[loc+0][loc+0] = -M2*AvgNum*dR_dr2;
	SourceJacobian[loc+0][loc+1] = -M2*AvgNum*dR_dm2;
	SourceJacobian[loc+0][loc+2] = -M2*AvgNum*dR_dn2;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] =  -M2*AvgNum*dR_dl2;
	SourceJacobian[loc+0][loc+nDim+1] = -M2*AvgNum*dR_de2;

	/* Derivative of the Sixth term in the source terms with all the conservative variables */

	SourceVector[loc+1] = r2*ec*EMF[1][0]/M2 + r2*nu21*(m1/r1 - m2/r2) + r2*nu23*(m3/r3 - m2/r2);

	SourceJacobian[loc+1][loc+0] = nu21*m1/r1 + nu23*m3/r3 ;//+ ec*ElectricField[0]/M2;
	SourceJacobian[loc+1][loc+1] = -nu21 - nu23;
	SourceJacobian[loc+1][loc+2] = ec*MagneticField[2]/M2;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] =  -ec*MagneticField[1]/M2;
	SourceJacobian[loc+1][loc+nDim+1]  = 0 + r2*(m1/r1 - m2/r2)*dv21_dT2*dT2_de2 + r2*(m3/r3 - m2/r2)*dv23_dT2*dT2_de2;


	/* Derivative of the Seventh term in the source terms wrt all the conservative variables */

	SourceVector[loc+2]  =  r2*ec*EMF[1][1]/M2 + r2*nu21*(n1/r1 - n2/r2) + r2*nu23*(n3/r3 - n2/r2);

	SourceJacobian[loc+2][loc+0] = nu21*n1/r1 + nu23*n3/r3 ;//+ ec*ElectricField[1]/M2;
	SourceJacobian[loc+2][loc+1] = -ec*MagneticField[2]/M2;
	SourceJacobian[loc+2][loc+2] = -nu21 - nu23 ;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] = ec*MagneticField[0]/M2;
	SourceJacobian[loc+2][loc+nDim+1] = 0 + r2*(n1/r1 - n2/r2)*dv21_dT2*dT2_de2 + r2*(n3/r3 - n2/r2)*dv23_dT2*dT2_de2;

	if (nDim == 3) {
		SourceVector[loc+3] = r2*ec*EMF[1][2]/M2 + r2*nu21*(l1/r1 - l2/r2) + r2*nu23*(l3/r3 - l2/r2);

		SourceJacobian[loc+3][loc+0]  = nu21*l1/r1 + nu23*l3/r3 ;//+ ec*ElectricField[2]/M2;
		SourceJacobian[loc+3][loc+1]  = ec*MagneticField[1]/M2;
		SourceJacobian[loc+3][loc+2]  = -ec*MagneticField[0]/M2;
		SourceJacobian[loc+3][loc+3] = -nu21 - nu23;
		SourceJacobian[loc+3][loc+nDim+1] = 0 + r2*(l1/r1 -l2/r2)*dv21_dT2*dT2_de2 + r2*(l3/r3 - l2/r2)*dv23_dT2*dT2_de2;
	}

	/* Derivative of the Eight term in the source terms wrt all the conservative variables */

	SourceVector[loc+nDim+1] = m2*ec*EMF[1][0]/M2 + n2*ec*EMF[1][1]/M2  + l2*ec*EMF[1][2]/M2 +
			m2*nu21*(m1/r1 - m2/r2) + m2*nu23*(m3/r3 - m2/r2)+
			n2*nu21*(n1/r1 - n2/r2) + n2*nu23*(n3/r3 - n2/r2) +
			l2*nu21*(l1/r1 - l2/r2) + l2*nu23*(l3/r3 - l2/r2) + QT2;

	SourceJacobian[loc+nDim+1][loc+0] = -1/r2*(m2*ec*EMF[1][0]/M2 + n2*ec*EMF[1][1]/M2  + l2*ec*EMF[1][2]/M2) + (nu21+nu23)*(m2*m2 + n2*n2 + l2*l2)/(r2*r2) ;// +  dQT2_dT2*dT2_dr2;
	SourceJacobian[loc+nDim+1][loc+1] = nu21*(m1/r1 - 2*m2/r2) + nu23*(m3/r3 - 2*m2/r2) + ec*EMF[1][0]/M2  + dQT2_dT2*dT2_dm2;
	SourceJacobian[loc+nDim+1][loc+2] = nu21*(n1/r1 - 2*n2/r2) + nu23*(n3/r3 - 2*n2/r2) + ec*EMF[1][1]/M2  + dQT2_dT2*dT2_dn2;
	if (nDim ==3) SourceJacobian[loc+nDim+1][loc+3] = nu21*(l1/r1 - 2*l2/r2) + nu23*(l3/r3 - 2*l2/r2) + ec*EMF[1][2]/M2  + dQT2_dT2*dT2_dl2;
	SourceJacobian[loc+nDim+1][loc+nDim+1] = 0 +
			(n2*(n1/r1 - n2/r2) + m2*(m1/r1 - m2/r2) + l2*(l1/r1 - l2/r2))*dv21_dT2*dT2_de2 +
			(n2*(n3/r3 - n2/r2) + m2*(m3/r3 - m2/r2) + l2*(l3/r3 - l2/r2))*dv23_dT2*dT2_de2 + dQT2_dT2*dT2_de2;



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

	SourceJacobian[loc+1][loc+0] = nu31*m1/r1 + nu32*m2/r2;// - ec*ElectricField[0]/M3;
	SourceJacobian[loc+1][loc+1] =  -nu31 - nu32;
	SourceJacobian[loc+1][loc+2] =  -ec*MagneticField[2]/M3;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] = ec*MagneticField[1]/M3;
	SourceJacobian[loc+1][loc+nDim+1] =  0 + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_de3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_de3;

	/* Derivative of the Eleventh term in the source terms wrt all the conservative variables */
	SourceVector[loc+2]  = -r3*ec*EMF[2][1]/M3 + r3*nu31*(n1/r1 - n3/r3) + r3*nu32*(n2/r2 - n3/r3);

	SourceJacobian[loc+2][loc+0] = nu31*n1/r1 + nu32*n2/r2 ;//- ec*ElectricField[1]/M3 ;
	SourceJacobian[loc+2][loc+1] =  ec*MagneticField[2]/M3;
	SourceJacobian[loc+2][loc+2] =  -nu31 - nu32 ;
	if (nDim ==3)  SourceJacobian[loc+2][loc+3] = -ec*MagneticField[0]/M3;
	SourceJacobian[loc+2][loc+nDim+1] =  0 + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_de3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_de3;

	/* Derivative of the Twelfth term in the source terms wrt all the conservative variables */
	if (nDim ==3) {
		SourceVector[loc+3]  = -r3*ec*EMF[2][2]/M3 + r3*nu31*(l1/r1 - l3/r3) + r3*nu32*(l2/r2 - l3/r3);
		SourceJacobian[loc+3][loc+0] = nu31*l1/r1 + nu32*l2/r2 ;//- ec*ElectricField[2]/M3 ;
		SourceJacobian[loc+3][loc+1] =  -ec*MagneticField[1]/M3;
		SourceJacobian[loc+3][loc+2] =   ec*MagneticField[0]/M3;
		SourceJacobian[loc+3][loc+3] =  -nu31 - nu32 ;
		SourceJacobian[loc+3][loc+nDim+1] =  0 + r3*(l1/r1 - l3/r3)*dv31_dT3*dT3_de3 + r3*(l2/r2 - l3/r3)*dv32_dT3*dT3_de3;
	}

	SourceVector[loc+nDim+1]  = -m3*ec*EMF[2][0]/M3 - n3*ec*EMF[2][1]/M3 - l3*ec*EMF[2][2]/M3 +
			m3*nu31*(m1/r1 - m3/r3) + m3*nu32*(m2/r2 - m3/r3) +
			n3*nu31*(n1/r1 - n3/r3) + n3*nu32*(n2/r2 - n3/r3) +
			l3*nu31*(l1/r1 - l3/r3) + l3*nu32*(l2/r2 - l3/r3) + QT3;

	SourceJacobian[loc+nDim+1][loc+0] =  -1/r3*(-m3*ec*EMF[2][0]/M3 - n3*ec*EMF[2][1]/M3 - l3*ec*EMF[2][2]/M3) + (nu31+nu32)*(m3*m3 + n3*n3 +l3*l3)/(r3*r3);// + dQT3_dT3*dT3_dr3
	SourceJacobian[loc+nDim+1][loc+1]  = nu31*(m1/r1 - 2*m3/r3) + nu32*(m2/r2 - 2*m3/r3) -ec*EMF[2][0]/M3 + dQT3_dT3*dT3_dm3;
	SourceJacobian[loc+nDim+1][loc+2]  = nu31*(n1/r1 - 2*n3/r3) + nu32*(n2/r2 - 2*n3/r3) -ec*EMF[2][1]/M3 + dQT3_dT3*dT3_dn3;
	if (nDim ==3)SourceJacobian[loc+nDim+1][loc+3]  = nu31*(l1/r1 - 2*l3/r3) + nu32*(l2/r2 - 2*l3/r3) -ec*EMF[2][2]/M3 + dQT3_dT3*dT3_dl3;
	SourceJacobian[loc+nDim+1][loc+nDim+1]  = 0 +
			(m3*(m1/r1 - m3/r3) +  n3*(n1/r1 - n3/r3) +  l3*(l1/r1 - l3/r3))*dv31_dT3*dT3_de3 +
			(m3*(m2/r2 - m3/r3)+  n3*(n2/r2 - n3/r3) + l3*(l2/r2 - l3/r3))*dv32_dT3*dT3_de3 + dQT3_dT3*dT3_de3;

#endif
	for (unsigned short iVar = 0; iVar < nVar; iVar ++) {
		val_residual[iVar] = SourceVector[iVar]*Volume;
		if (implicit) {
			for (unsigned short jVar = 0; jVar < nVar; jVar ++) {
				val_Jacobian_i[iVar][jVar] = SourceJacobian[iVar][jVar]*Volume;
			}
		}
	}
}

void CSourcePieceWise_Plasma::SetResidual_Axisymmetric(double *val_residual, CConfig *config) {
//************************************************//
// Please do not delete //SU2_CPP2C comment lines //
//************************************************//

//SU2_CPP2C START CSourcePieceWise_Plasma::SetResidual_Axisymmetric
//SU2_CPP2C CALL_LIST START
//SU2_CPP2C INVARS *U_i
//SU2_CPP2C OUTVARS *val_residual
//SU2_CPP2C VARS DOUBLE *Coord_i Volume
//SU2_CPP2C CALL_LIST END

//SU2_CPP2C DEFINE nDim nVar nSpecies nDiatomics

//SU2_CPP2C DECL_LIST START

//SU2_CPP2C DECL_LIST END


	double yinv, Gamma, Density_i, Energy_i, Energy_vib_i, Energy_el_i, Pressure_i, Enthalpy_i, Enthalpy_formation_i, Velocity_i, sq_vel;
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

		Gamma = config->GetSpecies_Gamma(iSpecies);
		Enthalpy_formation_i = config->GetEnthalpy_Formation(iSpecies);
		Density_i    = U_i[loc+0];
		Energy_i     = U_i[loc+nDim+1] / Density_i;
		Energy_vib_i = 0.0;
		if (iSpecies < nDiatomics)
			Energy_vib_i = U_i[loc+nDim+2] / Density_i;
		Energy_el_i  = 0.0;
		Pressure_i = (Gamma-1.0) * Density_i * (Energy_i - 1.0/2.0*sq_vel - Enthalpy_formation_i - Energy_vib_i - Energy_el_i);
		Enthalpy_i   = (U_i[loc+nDim+1] + Pressure_i) / Density_i;

		val_residual[loc+0] = yinv*Volume*U_i[loc+2];
		val_residual[loc+1] = yinv*Volume*U_i[loc+1]*U_i[loc+2]/U_i[loc+0];
		val_residual[loc+2] = yinv*Volume*U_i[loc+2]*U_i[loc+2]/U_i[loc+0];
		val_residual[loc+3] = yinv*Volume*Enthalpy_i*U_i[loc+2];
		if (iSpecies < nDiatomics)
			val_residual[loc+4] = yinv*Volume*U_i[loc+4]*U_i[loc+2]/U_i[loc+0];
	}
//SU2_CPP2C END CSourcePieceWise_Plasma::SetResidual_Axisymmetric
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
			SetResidual_Axisymmetric_ad(Residual_Baseline, Residual_New, config);
			for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[jVar][iVar] = Residual_New[jVar];
		}

/*		cout << "++++ AD Jacobian (Axisymmetric) ++++" << endl;
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        cout << "\t" << val_Jacobian_i[iVar][jVar];
      }
      cout << endl;
    }
    cin.get();*/

	}

	/*--- Calculate jacobian using finite differencing ---*/
	else if (config->GetKind_SourJac_Plasma() == FINITE_DIFF) {
		double FDEpsilon;

		/*--- Compute the baseline line residual ---*/
		SetResidual_Axisymmetric(Residual_Baseline, config);

		/*--- Compute forward finite differences ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			FDEpsilon = 1E-6*U_Baseline[iVar];
			if (FDEpsilon < 1E-14)
				FDEpsilon = 1E-14;

			/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
			U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
			for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
			SetResidual_Axisymmetric(Residual_New, config);

			/*--- Undo the change in U_i to keep constant the solution ---*/
			U_i[iVar] = U_i[iVar] - FDEpsilon;

			/*--- Save the new Jacobian (per column) ---*/
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
				/*				cout << "Val Jacobian [" << jVar << "][" << iVar << "]: " << val_Jacobian_i[jVar][iVar] << endl;
				cout << "Residual_New[" << jVar << "]: " << Residual_New[jVar] << endl;
				cout << "Residual_Baseline[" << jVar << "]: " << Residual_Baseline[jVar] << endl;
				cout << U_i[jVar]<< endl;
				cout << FDEpsilon << endl;
				cin.get();*/
			}
		}
/*		cout << "++++ FD Jacobian (Axisymmetric) ++++" << endl;
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << val_Jacobian_i[iVar][jVar];
			}
			cout << endl;
		}
		cin.get();*/
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
			SetResidual_Chemistry_ad(Residual_Baseline, Residual_New, config);
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_i[jVar][iVar] = Residual_New[jVar];
		}

/*		cout << "++++ AD Jacobian (Chem) ++++" << endl;
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << val_Jacobian_i[iVar][jVar];
			}
			cout << endl;
		}
		cin.get();*/

	}

	/*--- Calculate jacobian using finite differencing ---*/
	else if (config->GetKind_SourJac_Plasma() == FINITE_DIFF) {
		double FDEpsilon;

		/*--- Compute the baseline line residual ---*/
		SetResidual_Chemistry(Residual_Baseline, config);

		/*--- Compute forward finite differences ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			FDEpsilon = 1E-6*U_Baseline[iVar];
			if (FDEpsilon < 1E-14)
				FDEpsilon = 1E-14;

			/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
			U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
			for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
			SetResidual_Chemistry(Residual_New, config);

			/*--- Undo the change in U_i to keep constant the solution ---*/
			U_i[iVar] = U_i[iVar] - FDEpsilon;

			/*--- Save the new Jacobian (per column) ---*/
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
				/*cout << "Val Jacobian [" << jVar << "][" << iVar << "] " << val_Jacobian_i[jVar][iVar] << endl;
        cout << "Residual_New[" << jVar << "]" << Residual_New[jVar] << endl;
        cout << "Residual_Baseline[" << jVar << "]" << Residual_Baseline[jVar] << endl;
        cout << U_i[jVar];
        cout << FDEpsilon << endl; 
        cin.get();*/
			}
		}
/*		cout << "++++ FD Jacobian (Chem) ++++" << endl;
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << val_Jacobian_i[iVar][jVar];
			}
			cout << endl;
		}
		cin.get();*/
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

void CSourcePieceWise_Plasma::SetResidual_Chemistry(double *val_residual, CConfig *config) {
//************************************************//
// Please do not delete //SU2_CPP2C comment lines //
//************************************************//

//SU2_CPP2C START CSourcePieceWise_Plasma::SetResidual_Chemistry
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
  double Energy_vib, Energy_el, Vel2, Gamma, Enthalpy_formation, Gas_constant;
  
	T_min = 800;
	epsilon = 80;
  
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);			
		Energy_vib = 0.0;
		Energy_el  = 0.0;
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += U_i[iLoc+iDim+1] / U_i[iLoc+0] * U_i[iLoc+iDim+1] / U_i[iLoc+0];
		if (iSpecies < nDiatomics) {
			Energy_vib = U_i[iLoc+nDim+2] / U_i[iLoc+0];
		}
    Gamma = config->GetSpecies_Gamma(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);    
		Temp_tr_i[iSpecies] = (Gamma-1.0)/Gas_constant * (U_i[iLoc+nDim+1]/U_i[iLoc] - 0.5*Vel2 - Enthalpy_formation - Energy_vib - Energy_el); 
  }

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
				T_rxnf[iReaction] *= Temp_tr_i[iSpecies];
				counterFwd++;
			}
			/*--- Products ---*/
			if (jSpecies != nSpecies) {
				T_rxnb[iReaction] *= Temp_tr_i[jSpecies];
				counterBkw++;
			}
		}
		T_rxnf[iReaction] = exp(1.0/counterFwd*log(T_rxnf[iReaction]));
		T_rxnb[iReaction] = exp(1.0/counterBkw*log(T_rxnb[iReaction]));

		/*--- Apply a modified temperature to ease the stiffness at low temperatures ---*/
		T_rxnf[iReaction] = 0.5 * (T_rxnf[iReaction]+T_min + sqrt((T_rxnf[iReaction]-T_min)*(T_rxnf[iReaction]-T_min) + epsilon*epsilon));
		T_rxnb[iReaction] = 0.5 * (T_rxnb[iReaction]+T_min + sqrt((T_rxnb[iReaction]-T_min)*(T_rxnb[iReaction]-T_min) + epsilon*epsilon));

//		GetEq_Rxn_Coefficients(EqRxnConstants, config);

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
		//		ReactionRateBkw[iReaction] = ReactionRateFwd[iReaction] / Keq[iReaction];

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

	/*	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		w_dot[iSpecies] = 0.0;
		for (iReaction = 0; iReaction < nReactions; iReaction++) {
			w_dot[iSpecies] =   (RxnProducts[iReaction][iSpecies] - RxnReactants[iReaction][iSpecies])
	 * (fwdRxn[iReaction] - bkwRxn[iReaction]);
		}
		w_dot[iSpecies] = Molar_Mass[iSpecies]*w_dot[iSpecies];
		val_residual[iLoc] = w_dot[iSpecies]*Volume;
		for (iDim = 0; iDim < nDim; iDim++)
			val_residual[iLoc+1+iDim] = w_dot[iSpecies] * U_i[iLoc+1+iDim]/U_i[iLoc] * Volume;
		val_residual[iLoc+1+nDim] = w_dot[iSpecies] * U_i[iLoc+1+nDim]/U_i[iLoc] * Volume;
		if (iSpecies < nDiatomics)
			val_residual[iLoc+2+nDim] = w_dot[iSpecies] * U_i[iLoc+2+nDim]/U_i[iLoc] * Volume;
	}	*/
//SU2_CPP2C END CSourcePieceWise_Plasma::SetResidual_Chemistry
}

void CSourcePieceWise_Plasma::SetResidual_Chemistry(double *val_residual, double **val_Jacobian, CConfig *config) {
	int ***Reactions;			// Array of chemical constituents participating in iReaction.
	unsigned short iSpecies, jSpecies, iReaction, iVar, jVar, iDim, ii;
	unsigned short loc, iLoc, jLoc;
	unsigned short counterFwd, counterBkw;

	tol = 1E-60;
	Reactions = config->GetReaction_Map();
	nReactions = config->GetnReactions();

	/*--- Initialize all components of the residual and Jacobian to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++)
			val_Jacobian[iVar][jVar] = 0.0;
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
				T_rxnf[iReaction] *= Temp_tr_i[iSpecies];
				counterFwd++;
			}
			/*--- Products ---*/
			if (jSpecies != nSpecies) {
				T_rxnb[iReaction] *= Temp_tr_i[jSpecies];
				counterBkw++;
			}
		}
		T_rxnf[iReaction] = pow(T_rxnf[iReaction], 1.0/double(counterFwd));
		T_rxnb[iReaction] = pow(T_rxnb[iReaction], 1.0/double(counterBkw));

		//		GetEq_Rxn_Coefficients(EqRxnConstants, config);

		/*--- Calculate equilibrium extent of reaction ---*/
		//NOTE: Scalabrin implementation
		Keq[iReaction] = exp(EqRxnConstants[iReaction][0]*(T_rxnb[iReaction]/10000.0)
				+ EqRxnConstants[iReaction][1]
				                            + EqRxnConstants[iReaction][2]*log(10000.0/T_rxnb[iReaction])
				                            + EqRxnConstants[iReaction][3]*(10000.0/T_rxnb[iReaction])
				                            + EqRxnConstants[iReaction][4]*pow(10000.0/T_rxnb[iReaction], 2.0));

		/*--- Calculate reaction rate coefficients ---*/
		//NOTE: Park implementation
		ReactionRateFwd[iReaction] = Cf[iReaction] * pow(T_rxnf[iReaction], eta[iReaction]) * exp(-theta[iReaction]/T_rxnf[iReaction]);
		ReactionRateBkw[iReaction] = Cf[iReaction] * pow(T_rxnb[iReaction], eta[iReaction]) * exp(-theta[iReaction]/T_rxnb[iReaction]) / Keq[iReaction];

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

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			fwdRxn[iReaction] = fwdRxn[iReaction]*pow(0.001*U_i[loc+0]/Molar_Mass[iSpecies], RxnReactants[iReaction][iSpecies]);
			bkwRxn[iReaction] = bkwRxn[iReaction]*pow(0.001*U_i[loc+0]/Molar_Mass[iSpecies], RxnProducts[iReaction][iSpecies]);
			/*			cout << "Molar_Mass[iSpecies] " << Molar_Mass[iSpecies] << endl;
			 cout << "RxnReactants[" << iReaction << "][" << iSpecies << "] " << RxnReactants[iReaction][iSpecies] << endl;
			 cout << "RxnProducts[" << iReaction << "][" << iSpecies << "] " << RxnProducts[iReaction][iSpecies] << endl;
			 cout << "U_i[" << loc+0 << "] " << U_i[loc+0] << endl;
			 cout << "fwdRxn[ " << iReaction << "] " << fwdRxn[iReaction] << endl;
			 cout << "bkwRxn[ " << iReaction << "] " << bkwRxn[iReaction] << endl;
			 cin.get(); */
		}
		fwdRxn[iReaction] = 1000.0 * ReactionRateFwd[iReaction] * fwdRxn[iReaction];
		bkwRxn[iReaction] = 1000.0 * ReactionRateBkw[iReaction] * bkwRxn[iReaction];
		/*		cout << "ReactionRateBkw[" << iReaction << "] " << ReactionRateBkw[iReaction] << endl;
		 cout << "ReactionRateFwd[" << iReaction << "] " << ReactionRateFwd[iReaction] << endl;*/
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

	SetResidual_ElecForce(Residual_Baseline, config);

	/*--- Compute forward finite differences ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		FDEpsilon = 0.01*U_Baseline[iVar];
		if (FDEpsilon == 0)
			FDEpsilon = 1E-9;

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
}

void CSourcePieceWise_Plasma::SetResidual_ElecForce(double *val_residual, CConfig *config) {
	//	unsigned short iDim, iSpecies, loc;
	//	double F_es;

	/*--- Solve for the electrostatic force on each constituent ---*/
	/*	for (iDim = 0; iDim < nDim; iDim++) {
	 ElectricField[iDim] = -ConsVar_Grad[0][iDim];
	 for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
	 if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	 else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
	 F_es = config->GetParticle_ChargeNumber(iSpecies) * (U_i[loc+0]/Molar_Mass[iSpecies]) * ELECTRON_CHARGE * ElectricField[iDim];
	 val_residual[loc+1+iDim] = F_es * Volume;
	 }
	 } */

	/*--- Compute jacobian of electrostatic force terms ---*/

}

void CSourcePieceWise_Plasma::SetResidual_ElecForce(double *val_residual, double **val_Jacobian, CConfig *config) {

	unsigned short jVar;
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.0;
		for(jVar = 0; jVar < nVar; jVar++)
			val_Jacobian[iVar][jVar] = 0.0;
	}
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
			SetResidual_MomentumExch_ad(Residual_Baseline, Residual_New, config);
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_i[jVar][iVar] = Residual_New[jVar];
		}
    
/*    cout << "++++ AD Jacobian (Momentum) ++++" << endl;
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << val_Jacobian_i[iVar][jVar];
			}
			cout << endl;
		}
		cin.get();*/
	}

	/*--- Calculate jacobian using finite differencing ---*/
	else if (config->GetKind_SourJac_Plasma() == FINITE_DIFF) {

		double FDEpsilon;

		/*--- Compute the baseline line residual ---*/
		SetResidual_MomentumExch(Residual_Baseline, config);

		/*--- Compute forward finite differences ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			FDEpsilon = 1E-6*U_Baseline[iVar];
			if (FDEpsilon < 1E-14)
				FDEpsilon = 1E-14;

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
/*    cout << "++++ FD Jacobian (Momentum) ++++" << endl;
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << val_Jacobian_i[iVar][jVar];
			}
			cout << endl;
		}
		cin.get();*/
	}
}

void CSourcePieceWise_Plasma::SetResidual_MomentumExch(double *val_residual, CConfig *config) {
//************************************************//
// Please do not delete //SU2_CPP2C comment lines //
//************************************************//

//SU2_CPP2C START CSourcePieceWise_Plasma::SetResidual_MomentumExch
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
  double Energy_vib, Energy_el, Vel2, Gamma, Enthalpy_formation, Gas_constant;
  
	/*--- Initialization ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		for (iDim = 0; iDim < nDim; iDim++)
			P[iSpecies][iDim] = 0.0;
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = 0.0;
  
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);			
		Energy_vib = 0.0;
		Energy_el  = 0.0;
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += U_i[iLoc+iDim+1] / U_i[iLoc+0] * U_i[iLoc+iDim+1] / U_i[iLoc+0];
		if (iSpecies < nDiatomics) {
			Energy_vib = U_i[iLoc+nDim+2] / U_i[iLoc+0];
		}
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);    
		Temp_tr_i[iSpecies] = (Gamma-1.0)/Gas_constant * (U_i[iLoc+nDim+1]/U_i[iLoc] - 0.5*Vel2 - Enthalpy_formation - Energy_vib - Energy_el); 
  }


	electron_temperature = Temp_tr_i[nSpecies -1];

	/*--- Solve for momentum exchange between species due to collisions ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
			if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
			else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);

			if (iSpecies != jSpecies) {
/*				collisionArea = PI_NUMBER * ((Molecular_Diameter[iSpecies] + Molecular_Diameter[jSpecies])/2.0)
																																	* ((Molecular_Diameter[iSpecies] + Molecular_Diameter[jSpecies])/2.0); */
        T_control = sqrt(Temp_tr_i[iSpecies]*Temp_tr_i[jSpecies]);
        collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] 
                      * pow(T_control, Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) 
                            + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2]);
/*        collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] 
        * exp((Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) 
              + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2])*log(T_control));*/

				/*----- An if-condition for the case when a positive charge collides with an electron -----*/
				if( (ChargeNumber[iSpecies] == 1 && jSpecies == nSpecies -1) || (iSpecies == nSpecies-1 && ChargeNumber[jSpecies] == 1)) {
					radius_electronIonCollision = ELECTRON_CHARGE*ELECTRON_CHARGE/(32.0*FREE_PERMITTIVITY*BOLTZMANN_CONSTANT*electron_temperature);
					collisionArea = PI_NUMBER * radius_electronIonCollision * radius_electronIonCollision;
				}

				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]));
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
//SU2_CPP2C END CSourcePieceWise_Plasma::SetResidual_MomentumExch
}

void CSourcePieceWise_Plasma::SetResidual_MomentumExch(double *val_residual, double **val_Jacobian, CConfig *config) {
	unsigned short iDim, iSpecies, jSpecies, iVar, jVar, iLoc, jLoc;
	double collisionFreq, collisionArea, velocity_Species_i, velocity_Species_j, collisionVelocity;
	double radius_electronIonCollision, electron_temperature;
	double dTs_drhos, *dTs_drhous, dTs_dEs;
	double dcij_drhos, *dcij_drhous, dcij_dEs;

	dTs_drhous = new double[nDim];
	dcij_drhous = new double[nDim];

	/*--- Initialization ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		for (iDim = 0; iDim < nDim; iDim++)
			P[iSpecies][iDim] = 0.0;
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar ++)
			val_Jacobian[iVar][jVar] = 0.0;
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ )
		Temp_tr_i[iSpecies] = max(Temp_tr_i[iSpecies], 1000.0);

	electron_temperature = Temp_tr_i[nSpecies -1];

	/*--- Solve for momentum exchange between species due to collisions ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
			if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
			else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);

			if (iSpecies != jSpecies) {
				collisionArea = PI_NUMBER * pow((Molecular_Diameter[iSpecies] + Molecular_Diameter[jSpecies])/2.0,2.0);

				/*----- An if-condition for the case when a positive charge collides with an electron -----*/
				if( (ChargeNumber[iSpecies] == 1 && jSpecies == nSpecies -1) || (iSpecies == nSpecies-1 && ChargeNumber[jSpecies] == 1)) {
					radius_electronIonCollision = ELECTRON_CHARGE*ELECTRON_CHARGE/(32.0*FREE_PERMITTIVITY*BOLTZMANN_CONSTANT*electron_temperature);
					collisionArea = PI_NUMBER * radius_electronIonCollision * radius_electronIonCollision;
				}

				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]));
				collisionVelocity  = sqrt(velocity_Species_i*velocity_Species_i + velocity_Species_j*velocity_Species_j);
				collisionFreq = (U_i[jLoc+0] * collisionArea * collisionVelocity)/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);
				/*				cout << "Momentum Collision Freq: " << collisionFreq << endl;
				cout << "Momentum XSection: " << collisionArea << endl;*/

				if (config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == O2 || config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {
					for (iDim = 0; iDim < nDim; iDim++)
						P[iSpecies][iDim] 						 +=  U_i[iLoc+0]* collisionFreq * (U_i[jLoc+1+iDim]/U_i[jLoc+0] - U_i[iLoc+1+iDim]/U_i[iLoc+0]);

					dTs_drhos = -U_i[iLoc+nDim+1]/(U_i[iLoc+0]*U_i[iLoc+0]);
					for (iDim = 0; iDim < nDim; iDim++)
						dTs_drhos += U_i[iLoc+iDim+1]*U_i[iLoc+iDim+1]/(pow(U_i[iLoc+0],3.0));
					if (iSpecies < nDiatomics)
						dTs_drhos += U_i[iLoc+nDim+2]/(U_i[iLoc+0]*U_i[iLoc+0]);
					//Need to include electronic energy!!! E_el
					for (iDim = 0; iDim < nDim; iDim++)
						dTs_drhous[iDim] = -U_i[iLoc+iDim+1]/(U_i[iLoc+0]*U_i[iLoc+0]);
					dTs_dEs = 1.0/U_i[iLoc+0];

					dcij_drhos = 0.5*pow(collisionVelocity,-0.5)*8.0*BOLTZMANN_CONSTANT/(PI_NUMBER*Molecular_Mass[iSpecies])*dTs_drhos;
					for (iDim = 0; iDim < nDim; iDim++)
						dcij_drhous[iDim] = 0.5*pow(collisionVelocity,-0.5)*8.0*BOLTZMANN_CONSTANT/(PI_NUMBER*Molecular_Mass[iSpecies])*dTs_drhous[iDim];
					dcij_dEs = 0.5*pow(collisionVelocity,-0.5)*8.0*BOLTZMANN_CONSTANT/(PI_NUMBER*Molecular_Mass[iSpecies])*dTs_dEs;

					for (iDim = 0; iDim < nDim; iDim++) {
						val_Jacobian[iLoc+iDim+1][iLoc+0]			 += collisionArea/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies])
																																																													*(/*dcij_drhos*(U_i[iLoc+0]*U_i[jLoc+iDim+1] - U_i[jLoc+0]*U_i[iLoc+iDim+1])*/
																																																															+ collisionVelocity*U_i[jLoc+iDim+1])*Volume;
						val_Jacobian[iLoc+iDim+1][iLoc+iDim+1] += collisionArea/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies])
																																																													*(/*dcij_drhous[iDim]*(U_i[iLoc+0]*U_i[jLoc+iDim+1] - U_i[jLoc+0]*U_i[iLoc+iDim+1])*/
																																																															- collisionVelocity*U_i[jLoc+0])*Volume;
						/*val_Jacobian[iLoc+iDim+1][iLoc+nDim+1] += collisionArea/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies])
						 *dcij_dEs*(U_i[iLoc+0]*U_i[jLoc+iDim+1] - U_i[jLoc+0]*U_i[iLoc+iDim+1])*Volume;*/
						val_Jacobian[iLoc+nDim+1][iLoc+0]     += collisionArea/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies])
																								    																																		 * collisionVelocity * U_i[jLoc+0]/(U_i[iLoc+0]*U_i[iLoc+0])*U_i[iLoc+iDim+1]*U_i[iLoc+iDim+1]*Volume;
						val_Jacobian[iLoc+nDim+1][iLoc+iDim+1] = collisionArea/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies])
																																																				 * collisionVelocity * (U_i[jLoc+iDim+1] - U_i[jLoc+0]/U_i[iLoc+0]*2.0*U_i[iLoc+iDim+1])*Volume;

						/*						val_Jacobian[iLoc+iDim+1][iLoc + 0] 	 += collisionFreq * U_i[jLoc+1+iDim]/U_i[jLoc+0] * Volume ;
						val_Jacobian[iLoc+iDim+1][iLoc + iDim+1] -= collisionFreq * Volume;
						val_Jacobian[iLoc+iDim+1][iLoc + nDim+1]  = 0.0;*/
					}
				} else {
					for (iDim = 0; iDim < nDim; iDim++) {
						P[iSpecies][iDim] 						 +=  U_i[iLoc+0]* collisionFreq * (U_i[jLoc+1+iDim]/U_i[jLoc+0] - U_i[iLoc+1+iDim]/U_i[iLoc+0]);
						val_Jacobian[iLoc+iDim+1][iLoc + 0] 	 += collisionFreq * U_i[jLoc+1+iDim]/U_i[jLoc+0] * Volume ;
						val_Jacobian[iLoc+iDim+1][iLoc + iDim+1] -= collisionFreq * Volume;
						val_Jacobian[iLoc+iDim+1][iLoc + nDim+1]  = 0.0;
					}
				}
			}
		}
		for (iDim = 0; iDim < nDim; iDim++) {
			val_residual[iLoc+1+iDim] = P[iSpecies][iDim] * Volume;
			val_residual[iLoc+nDim + 1] += P[iSpecies][iDim]* U_i[iLoc+1+iDim]/U_i[iLoc+0] * Volume;

		}

	}
	delete[] dTs_drhous;
	delete[] dcij_drhous;
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
			SetResidual_EnergyExch_ad(Residual_Baseline, Residual_New, config);
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_i[jVar][iVar] = Residual_New[jVar];
		}
    
/*    cout << "++++ AD Jacobian (Energy) ++++" << endl;
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << val_Jacobian_i[iVar][jVar];
			}
			cout << endl;
		}
		cin.get();*/
	}

	/*--- Calculate jacobian using finite differencing ---*/
	else if (config->GetKind_SourJac_Plasma() == FINITE_DIFF) {
		double FDEpsilon;

		/*--- Compute the baseline line residual ---*/
		SetResidual_EnergyExch(Residual_Baseline, config);

		/*--- Compute forward finite differences ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			FDEpsilon = 1E-6*U_Baseline[iVar];
			if (FDEpsilon < 1E-14)
				FDEpsilon = 1E-14;

			/*--- Recompute the residual, perturbation in the iVar component of the solution ---*/
			U_i[iVar] = U_Baseline[iVar] + FDEpsilon;
			for (jVar = 0; jVar < nVar; jVar++) Residual_New[jVar] = 0.0;
			SetResidual_EnergyExch(Residual_New, config);

			/*--- Undo the change in U_i to keep constant the solution ---*/
			U_i[iVar] = U_i[iVar] - FDEpsilon;

			/*--- Save the new Jacobian (per column) ---*/
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_i[jVar][iVar] = (Residual_New[jVar] - Residual_Baseline[jVar])/FDEpsilon;
			}
		}
/*    cout << "++++ FD Jacobian (Energy) ++++" << endl;
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				cout << "\t" << val_Jacobian_i[iVar][jVar];
			}
			cout << endl;
		}
		cin.get();*/
	}
}


void CSourcePieceWise_Plasma::SetResidual_EnergyExch(double *val_residual, double **val_Jacobian, CConfig *config) {
	unsigned short iDim, iSpecies, jSpecies, iLoc, jLoc, iVar, jVar;
	double collisionFreq, collisionArea, velocity_Species_i, velocity_Species_j, collisionVelocity, vel_dot_prod;
	double dTs_drhos, *dTs_drhous, dTs_dEs;

	dTs_drhous = new double[nDim];

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
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Q_elastic[iSpecies] = 0.0;
	}
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian[iVar][jVar] = 0.0;
		}
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){
			if (jSpecies != iSpecies) {
				if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
				else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
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


				//				collisionXSection = PI_NUMBER * pow((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0,2.0);
				collisionArea = PI_NUMBER * pow((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0,2.0);
				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]));
				collisionVelocity = sqrt(velocity_Species_i*velocity_Species_i + velocity_Species_j*velocity_Species_j);
				collisionFreq = (U_i[jLoc+0] * collisionArea * collisionVelocity)/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);
				//				collisionFreq = (U_i[jLoc+0]/Molecular_Mass[jSpecies] * collisionArea * collisionVelocity);

				/*				collisionFreq = (U_i[jLoc+0]/Molecular_Mass[jSpecies]) * collisionXSection
				 * sqrt( 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies])
								+8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]) );*/
				vel_dot_prod = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					vel_dot_prod += (U_i[iLoc+1+iDim] - U_i[jLoc+1+iDim])*U_i[iLoc+1+iDim];
				}

				/* Exchange from Lee and Sutton, heavy particles */
				/*				Q_elastic[iSpecies] +=   (U_i[iLoc+0]*AVOGAD_CONSTANT/Molar_Mass[iSpecies]) * 2.0
				 * (Molecular_Mass[iSpecies]*Molecular_Mass[jSpecies])/(Molecular_Mass[iSpecies] + Molecular_Mass[jSpecies])
															 / Molecular_Mass[jSpecies] * collisionFreq * (3.0/2.0*BOLTZMANN_CONSTANT*(Temp_tr_i[jSpecies] - Temp_tr_i[iSpecies]));
						+ (U_i[iLoc]*AVOGAD_CONSTANT/Molar_Mass[iSpecies])
				 * (Molecular_Mass[iSpecies]*Molecular_Mass[jSpecies])/(Molecular_Mass[iSpecies] + Molecular_Mass[jSpecies])
				 * collisionFreq * vel_dot_prod;*/
				double coefficient;
				coefficient = 2.0*U_i[iLoc+0] / (Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]) * collisionFreq * 3.0/2.0*BOLTZMANN_CONSTANT;
				Q_elastic[iSpecies] += coefficient*(Temp_tr_i[jSpecies] - Temp_tr_i[iSpecies]);

				dTs_drhos = -U_i[iLoc+nDim+1]/(U_i[iLoc+0]*U_i[iLoc+0]);
				for (iDim = 0; iDim < nDim; iDim++)
					dTs_drhos += U_i[iLoc+iDim+1]*U_i[iLoc+iDim+1]/(pow(U_i[iLoc+0],3.0));
				if (iSpecies < nDiatomics)
					dTs_drhos += U_i[iLoc+nDim+2]/(U_i[iLoc+0]*U_i[iLoc+0]);
				//Need to include electronic energy!!! E_el
				for (iDim = 0; iDim < nDim; iDim++)
					dTs_drhous[iDim] = -U_i[iLoc+iDim+1]/(U_i[iLoc+0]*U_i[iLoc+0]);
				dTs_dEs = 1.0/U_i[iLoc+0];

				val_Jacobian[iLoc+nDim+1][iLoc+0] += coefficient/U_i[iLoc+0]*(Temp_tr_i[jSpecies]-Temp_tr_i[iSpecies]-U_i[iLoc+0]*dTs_drhos)*Volume;
				for (iDim = 0; iDim < nDim; iDim++)
					val_Jacobian[iLoc+nDim+1][iLoc+iDim+1] -= coefficient*dTs_drhous[iDim]*Volume;
				val_Jacobian[iLoc+nDim+1][iLoc+nDim+1] -= coefficient*dTs_dEs*Volume;
				if (iSpecies < nDiatomics)
					cout << "NEED TO IMPLEMENT THE CONTRIBUTION TO ENERGY JACOBIAN FROM VIBRATIONAL ENERGY.  File: SetResidual_EnergyExch..." << endl;

				/*				cout << "NumDensity_i: " << (U_i[iLoc+0]*AVOGAD_CONSTANT/Molar_Mass[iSpecies]) << endl;
				cout << "2*Reduced Mass: " << 2.0 * (Molecular_Mass[iSpecies]*Molecular_Mass[jSpecies])/(Molecular_Mass[iSpecies] + Molecular_Mass[jSpecies]) << endl;
				cout << "Molecular Mass_j" << Molecular_Mass[jSpecies] << endl;
				cout << "Energy Collision Frequency: " << collisionFreq << endl;
				cout << "Boltzmann constant: " << BOLTZMANN_CONSTANT << endl;
				cout << "Overall coefficient: " << (U_i[iLoc+0]*AVOGAD_CONSTANT/Molar_Mass[iSpecies]) * 2.0
				 * (Molecular_Mass[iSpecies]*Molecular_Mass[jSpecies])/(Molecular_Mass[iSpecies] + Molecular_Mass[jSpecies])
				/ Molecular_Mass[jSpecies] * collisionFreq * BOLTZMANN_CONSTANT << endl;

				cout << "Q_elastic[" << iSpecies << "]: " << Q_elastic[iSpecies] << endl;
				cout << "Energy[" << iSpecies << "]: " << U_i[iLoc+nDim+1] << endl;
				cin.get();*/

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
	/*	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] *= 0.1;
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian[iVar][jVar] *= 0.1;
		}
	}*/

	delete[] dTs_drhous;
}

void CSourcePieceWise_Plasma::SetResidual_EnergyExch(double *val_residual, CConfig *config) {
//************************************************//
// Please do not delete //SU2_CPP2C comment lines //
//************************************************//

//SU2_CPP2C START CSourcePieceWise_Plasma::SetResidual_EnergyExch
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
  double Energy_vib, Energy_el, Vel2, Gamma, Enthalpy_formation, Gas_constant;
  
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
		Energy_vib = 0.0;
		Energy_el  = 0.0;
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += U_i[iLoc+iDim+1] / U_i[iLoc+0] * U_i[iLoc+iDim+1] / U_i[iLoc+0];
		if (iSpecies < nDiatomics) {
			Energy_vib = U_i[iLoc+nDim+2] / U_i[iLoc+0];
		}
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);    
		Temp_tr_i[iSpecies] = (Gamma-1.0)/Gas_constant * (U_i[iLoc+nDim+1]/U_i[iLoc] - 0.5*Vel2 - Enthalpy_formation - Energy_vib - Energy_el); 
    SpeciesPressure_i[iSpecies] = (Gamma-1.0) * U_i[iLoc] * (U_i[iLoc+nDim+1]/U_i[iLoc] - 0.5*Vel2 - Enthalpy_formation - Energy_vib - Energy_el);
  }
  

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){
			if (jSpecies != iSpecies) {
				if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
				else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);

/*				collisionArea = PI_NUMBER * ((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0)
																																	* ((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0);*/
        T_control = sqrt(Temp_tr_i[iSpecies]*Temp_tr_i[jSpecies]);
        collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] 
                        * pow(T_control, Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) 
                              + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2]);
/*        collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] 
        * exp((Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) 
               + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2])*log(T_control));*/
        
				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]));
				collisionVelocity = sqrt(velocity_Species_i*velocity_Species_i + velocity_Species_j*velocity_Species_j);
				collisionFreq = (U_i[jLoc+0] * collisionArea * collisionVelocity)/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);
				vel_dot_prod = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					vel_dot_prod += (U_i[iLoc+1+iDim] - U_i[jLoc+1+iDim])*U_i[iLoc+1+iDim];
				}

				/* Exchange from Lee and Sutton, heavy particles */
				coefficient = 2.0*U_i[iLoc+0] / (Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]) * collisionFreq * 3.0/2.0*BOLTZMANN_CONSTANT;
				Q_elastic[iSpecies] += coefficient*(Temp_tr_i[jSpecies] - Temp_tr_i[iSpecies]);
			}
		}
		val_residual[iLoc+1+nDim] =  Q_elastic[iSpecies] * Volume;
	}


	/*--- Translational-rotational & vibrational energy exchange via inelastic collisions ---*/
	//Comment: Landau-Teller formulation
	//Comment: Based on Scalabrin and Candler.  May need to re-visit with much more detail and derive for my formulation

	double tau_sr, tau_ps, LimitingXSection, AvgMolecularSpeed, ReducedMass, A_sr, estar_vs, e_vs, q_tr_vs;
	double MixturePressure, MixtureNumDensity;

	/*--- Calculate mixture quantities ---*/
	MixturePressure   = 0.0;
	MixtureNumDensity = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) {
      iLoc = (nDim+3)*iSpecies; 
      CharVibTemp[iSpecies] = config->GetCharVibTemp(iSpecies);
    } else 
      iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		MixturePressure += SpeciesPressure_i[iSpecies];
		MixtureNumDensity += U_i[iLoc+0]/Molecular_Mass[iSpecies];
	}

	for (iSpecies = 0; iSpecies < nDiatomics; iSpecies++) {
		iLoc = (nDim+3)*iSpecies;
		Q_tv[iSpecies] = 0.0;
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){
			if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
			else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);

			/*--- Calculate Landau-Teller relaxation time ---*/
			ReducedMass = Molar_Mass[iSpecies]*Molar_Mass[jSpecies]/(Molar_Mass[iSpecies]+Molar_Mass[jSpecies]);
			LimitingXSection = 1E-20*pow(50000/Temp_tr_i[iSpecies], 2.0);
//      LimitingXSection = 1E-20*exp(2.0*log(50000/Temp_tr_i[iSpecies]));
      
			AvgMolecularSpeed = sqrt(8.0*UNIVERSAL_GAS_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molar_Mass[iSpecies]));
      
			A_sr = 1.16 * 1E-3 * sqrt(ReducedMass) * pow(CharVibTemp[iSpecies], 4.0/3.0);
//      A_sr = 1.16 * 1E-3 * sqrt(ReducedMass) * exp(4.0/3.0*log(CharVibTemp[iSpecies]));
      
			tau_sr =   101325/(SpeciesPressure_i[iSpecies]+SpeciesPressure_i[jSpecies])
																			* exp(A_sr*(pow(sqrt(Temp_tr_i[iSpecies]*Temp_tr_i[jSpecies]),-1.0/3.0) - 0.015*pow(ReducedMass,0.25)) - 18.42);
/*      tau_sr =   101325/(SpeciesPressure_i[iSpecies]+SpeciesPressure_i[jSpecies])
                  * exp(A_sr*(exp(-1.0/3.0 * log(sqrt(Temp_tr_i[iSpecies]*Temp_tr_i[jSpecies]))) - 0.015*exp(0.25*log(ReducedMass))) - 18.42);*/
			tau_ps = 1.0/(LimitingXSection*AvgMolecularSpeed*MixtureNumDensity);

			estar_vs = UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies] * CharVibTemp[iSpecies]/(exp(CharVibTemp[iSpecies]/Temp_tr_i[jSpecies])-1.0);
			e_vs = U_i[iLoc+nDim+2]/U_i[iLoc+0];

			/*--- Energy transferred per-molecule from species r to species s ---*/
			q_tr_vs = Molecular_Mass[iSpecies] * (estar_vs - e_vs)/(tau_sr + tau_ps);
			//			Q_tv[iSpecies] += U_i[iLoc+0] * (estar_vs - e_vs)/(tau_sr + tau_ps);

			/*--- Convert to energy per volume for r and s species and multiply by volume for residual ---*/
			val_residual[iLoc+nDim+2] += q_tr_vs*U_i[iLoc+0]/Molecular_Mass[iSpecies]*Volume;
			val_residual[iLoc+nDim+1] += q_tr_vs*U_i[iLoc+0]/Molecular_Mass[iSpecies]*Volume;
			val_residual[jLoc+nDim+1] -= q_tr_vs*U_i[jLoc+0]/Molecular_Mass[jSpecies]*Volume;
		}
		//		val_residual[iLoc+nDim+2] = Q_tv[iSpecies]*Volume;
	}

//SU2_CPP2C END CSourcePieceWise_Plasma::SetResidual_EnergyExch
}


void CSourcePieceWise_Plasma::SetResidual_EnergyExch(double *val_residual, double *val_residual_ElecForce, double **val_Jacobian, CConfig *config) {
	unsigned short iDim, iSpecies, jSpecies, iLoc, jLoc, jVar;
	double collisionFreq, collisionXSection;
	double radius_electronIonCollision, velocity_Species_i, velocity_Species_j,collisionVelocity;
	double electron_temperature, dQ_elastic_dT;
	double Gamma, Gas_Constant, Cv, dT_drhoU, dT_de;

	/*--- Energy transport via electrostatic work ---*/
	/*	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
	 if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	 else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
	 for (iDim = 0; iDim < nDim; iDim++) {
	 val_residual[loc+1+nDim] = val_residual[loc+1+nDim] + (val_residual_ElecForce[loc+1+iDim]/Volume) * (U_i[loc+1+iDim]/U_i[loc]) * Volume;
	 }
	 }*/

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Q_elastic[iSpecies] = 0.0;
	}
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = 0.0;
		for (jVar = 0; jVar < nVar; jVar ++)
			val_Jacobian[iVar][jVar] = 0.0;
	}
	electron_temperature = Temp_tr_i[nSpecies-1];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		dQ_elastic_dT = 0.0;
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){

			if (jSpecies != iSpecies) {
				if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
				else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);

				collisionXSection = PI_NUMBER * pow((Molecular_Diameter[iSpecies] + Molecular_Diameter[jSpecies])/2.0,2.0);

				/*----- An if-condition for the case when a positive charge collides with an electron -----*/
				if( (ChargeNumber[iSpecies] == 1 && jSpecies == nSpecies -1) || (iSpecies == nSpecies-1 && ChargeNumber[jSpecies] == 1)) {
					radius_electronIonCollision = ELECTRON_CHARGE*ELECTRON_CHARGE/(32.0*FREE_PERMITTIVITY*BOLTZMANN_CONSTANT*electron_temperature);
					collisionXSection = PI_NUMBER * radius_electronIonCollision * radius_electronIonCollision;
				}

				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass[jSpecies]));
				collisionVelocity  = sqrt(velocity_Species_i*velocity_Species_i + velocity_Species_j*velocity_Species_j);
				collisionFreq = (U_i[jLoc+0] * collisionXSection * collisionVelocity)/(Molecular_Mass[jSpecies]+Molecular_Mass[iSpecies]);

				/* Exchange from Lee and Sutton, heavy particles */
				Gas_Constant = config->GetSpecies_Gas_Constant(jSpecies);
				Gamma = config->GetSpecies_Gamma(jSpecies);
				Cv = Gas_Constant/(Gamma-1.0);
				Q_elastic[iSpecies] += 2.0*U_i[iLoc+0] * Cv*collisionFreq*Molecular_Mass[jSpecies]/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies])* (Temp_tr_i[jSpecies] - Temp_tr_i[iSpecies]);
				dQ_elastic_dT += -2.0  * U_i[iLoc+0] *Cv* collisionFreq*Molecular_Mass[jSpecies]/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);

				val_Jacobian[iLoc+nDim+1][iLoc+0] = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					dT_drhoU = -U_i[iLoc+iDim+1]/U_i[iLoc+0]/(U_i[iLoc+0]*Cv);
					val_Jacobian[iLoc+nDim+1][iLoc+iDim+1] += dQ_elastic_dT*dT_drhoU * Volume;
				}
				dT_de = 	1.0/(Cv*U_i[iLoc + 0 ]);
				val_Jacobian[iLoc+nDim+1][iLoc+nDim+1] += dQ_elastic_dT*dT_de * Volume;
			}
		}
		val_residual[iLoc+1+nDim] =  Q_elastic[iSpecies] * Volume;
	}
}

CSourcePieceWise_Plasma_Air::CSourcePieceWise_Plasma_Air(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics,
		unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar, val_nSpecies, val_nDiatomics, val_nMonatomics, config) {

	Kb = BOLTZMANN_CONSTANT;

	zero = 1E-15;
	nS = nSpecies;
	Molar_Mass = new double [nSpecies];
	Species_Charge = new double [nSpecies];
	Cv_heatcap = new double [nSpecies];
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

	//	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
	//		Mass[iSpecies] = config->GetMolar_Mass(iSpecies) / AVOGAD_CONSTANT;


	dTemperature_drhoU = new double*[nSpecies];
	dTemperature_Energy = new double [nSpecies];
	QT =  new double [nSpecies];
	Collision_Freq_Heat_Transfer = new double* [nSpecies];
	Collision_Freq_Momentum_Tranfer = new double* [nSpecies];
	CollisionArea = new double* [nSpecies];
	CollisionVelo = new double* [nSpecies];
	dQT_dTemperature = new double *[nSpecies];


	nReactions = config->GetnReactions();
	Reactants = 0;
	Products = 1;
	RateofReaction = new double[nReactions];

	Reactions = new int**[nReactions];
	for (iReactions = 0; iReactions < nReactions; iReactions ++){
		Reactions[iReactions] = new int*[2];
		Reactions[iReactions][Reactants] = new int [3];
		Reactions[iReactions][Products] = new int [3];
	}

	for (iReactions = 0; iReactions < nReactions; iReactions ++){
		for (iVar = 0; iVar < 3; iVar ++ ) {
			Reactions[iReactions][Reactants][iVar] = 0;
			Reactions[iReactions][Products][iVar] = 0;
		}
		RateofReaction[iReactions] = 0.0;
	}

	Reactions = config->GetReaction_Map();

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		Molar_Mass[iSpecies] = config->GetMolar_Mass(iSpecies) ;
		Species_Charge[iSpecies] = config->GetCharge_Number(iSpecies);
		Gamma  =  config->GetSpecies_Gamma(iSpecies);
		n[iSpecies] = 1.0;
		Cv_heatcap[iSpecies] = 1/(Gamma-1.0)*config->GetSpecies_Gas_Constant(iSpecies);

		MassSource[iSpecies] = 0.0;
		EnergySource[iSpecies] = 0.0;
		Density[iSpecies] = 0.0;
		Pressure[iSpecies] = 0.0;
		Temperature[iSpecies] = 0.0;
		Energy[iSpecies] = 0.0;
		Energy_el = 0.0;
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
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies ++) {
			Collision_Freq_Heat_Transfer[iSpecies][jSpecies]  = 0.0;
			Collision_Freq_Momentum_Tranfer[iSpecies][jSpecies]  = 0.0;
			CollisionArea[iSpecies][jSpecies]  = 0.0;
			CollisionVelo[iSpecies][jSpecies]  = 0.0;

			dQT_dTemperature[iSpecies][jSpecies]  = 0.0;
		}
	}

	n[nSpecies] = 1.0;
	MassSource[nSpecies] = 0.0;

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

	ec = ELECTRON_CHARGE;    eps0 = FREE_PERMITTIVITY;
	Te = 10000;

	double r12 = 2E-10; //m

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ )
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies ++ )
			CollisionArea[iSpecies][jSpecies] = PI_NUMBER*r12*r12; //TEMPORARY, CORRECT LATER

}

CSourcePieceWise_Plasma_Air::~CSourcePieceWise_Plasma_Air(void) {

}

void CSourcePieceWise_Plasma_Air::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	double Gas_Constant;

	//cout << " in source.cpp " <<" nSpecies = " << nSpecies << " nDiatomics = " << nDiatomics << endl;



	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		Density[iSpecies] = U_i[loc + 0];
		Energy[iSpecies] = U_i[loc + nDim+1];
		Energy_vib = 0.0;
		if ( iSpecies < nDiatomics )
			Energy_vib = U_i[loc+nDim+2] / 	Density[iSpecies];
		Energy_el = 0.0;
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += U_i[loc+iDim+1]/U_i[loc+0] * U_i[loc+iDim+1]/U_i[loc+0];
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Pressure[iSpecies] = (Gamma-1.0) * Density[iSpecies] * (Energy[iSpecies] - 1.0/2.0*Vel2 - Enthalpy_formation - Energy_vib - Energy_el);
		n[iSpecies] = Density[iSpecies]/(Molar_Mass[iSpecies]*AVOGAD_CONSTANT); 			// m-3
		for (iDim = 0; iDim < nDim; iDim ++) {
			rhoU[iSpecies][iDim]  = U_i[loc + iDim+1];
			velocity[iSpecies][iDim] = rhoU[iSpecies][iDim]/ Density[iSpecies];
		}
		Gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
		Temperature[iSpecies] = Pressure[iSpecies]/(Gas_Constant*Density[iSpecies]);
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
	double Mass_Species_i, Mass_Species_j;
	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
		for (jSpecies =0; jSpecies < nSpecies; jSpecies ++) {
			Mass_Species_i = Molar_Mass[iSpecies]/AVOGAD_CONSTANT;
			Mass_Species_j = Molar_Mass[jSpecies]/AVOGAD_CONSTANT;

			CollisionVelo[iSpecies][jSpecies] = sqrt(8.0*BOLTZMANN_CONSTANT/PI_NUMBER* ( Temperature[iSpecies]/Mass_Species_i+Temperature[jSpecies]/Mass_Species_j));
			Collision_Freq_Momentum_Tranfer[iSpecies][jSpecies] = Density[jSpecies]/(Mass_Species_j+ Mass_Species_i) * CollisionArea[iSpecies][jSpecies] * CollisionVelo[iSpecies][jSpecies];
		}
		CollisionVelo[iSpecies][iSpecies] = 0.0;
		Collision_Freq_Momentum_Tranfer[iSpecies][iSpecies] = 0.0;
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		for (iDim = 0; iDim < nDim; iDim ++ ) {
			Mass_Species_i = Molar_Mass[iSpecies]/AVOGAD_CONSTANT;
			MomentumSource[iSpecies][iDim] = Density[iSpecies]*ELECTRON_CHARGE*Species_Charge[iSpecies]/Mass_Species_i*EMF[iSpecies][iDim];
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
		Mass_Species_i = Molar_Mass[iSpecies]/AVOGAD_CONSTANT;
		SourceVector_i[loc + 0] 			= Mass_Species_i*MassSource[iSpecies];
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
		Energy_vib = 0.0;
		if ( iSpecies < nDiatomics )
			Energy_vib = U_i[loc+nDim+2] / 	Density[iSpecies];
		Energy_el = 0.0;
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += U_i[loc+iDim+1]/U_i[loc+0] * U_i[loc+iDim+1]/U_i[loc+0];
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Pressure[iSpecies] = (Gamma-1.0) * Density[iSpecies] * (Energy[iSpecies] - 1.0/2.0*Vel2 - Enthalpy_formation - Energy_vib - Energy_el);
		n[iSpecies] = Density[iSpecies]/(Molar_Mass[iSpecies]*AVOGAD_CONSTANT); 			// m-3
		for (iDim = 0; iDim < nDim; iDim ++) {
			rhoU[iSpecies][iDim]  = U_i[loc + iDim+1];
			velocity[iSpecies][iDim] = rhoU[iSpecies][iDim]/ Density[iSpecies];
		}
		Gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
		Temperature[iSpecies] = Pressure[iSpecies]/(Gas_Constant*Density[iSpecies]);
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
			Mass_Species_i = Molar_Mass[product]/AVOGAD_CONSTANT;
			product = Reactions[iReactions][1][iProducts];
			MassSource[product] += RateofReaction[iReactions]*numdensity_Reactants*1E6 * Mass_Species_i; // Mass source in 1/(m3.s), later multiplied by species mass
		}

		for (iReactants = 0; iReactants < 3; iReactants ++) {
			Mass_Species_i = Molar_Mass[reactant]/AVOGAD_CONSTANT;
			reactant = Reactions[iReactions][0][iReactants];
			MassSource[reactant] -= RateofReaction[iReactions]*numdensity_Reactants*1E6 * Mass_Species_i; // Mass source in 1/(m3.s), later multiplied by species mass
		}
	}






	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
		Mass_Species_i = Molar_Mass[iSpecies]/AVOGAD_CONSTANT;

		for (jSpecies =0; jSpecies < nSpecies; jSpecies ++) {
			Mass_Species_j = Molar_Mass[jSpecies]/AVOGAD_CONSTANT;

			CollisionVelo[iSpecies][jSpecies] = sqrt(8.0*BOLTZMANN_CONSTANT/PI_NUMBER* ( Temperature[iSpecies]/Mass_Species_i+Temperature[jSpecies]/Mass_Species_j));
			Collision_Freq_Momentum_Tranfer[iSpecies][jSpecies] = Density[jSpecies]/(Mass_Species_j+ Mass_Species_i) * CollisionArea[iSpecies][jSpecies] * CollisionVelo[iSpecies][jSpecies];
		}
		CollisionVelo[iSpecies][iSpecies] = 0.0;
		Collision_Freq_Momentum_Tranfer[iSpecies][iSpecies] = 0.0;
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		for (iDim = 0; iDim < nDim; iDim ++ ) {
			Mass_Species_i = Molar_Mass[iSpecies]/AVOGAD_CONSTANT;
			MomentumSource[iSpecies][iDim] = Density[iSpecies]*ELECTRON_CHARGE*Species_Charge[iSpecies]/Mass_Species_i*EMF[iSpecies][iDim];
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


CSource_Template::CSource_Template(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSource_Template::~CSource_Template(void) { }

void CSource_Template::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) { }
