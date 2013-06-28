/*!
 * \file numerics_source.cpp
 * \brief This file contains all the source term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
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

CSourceNothing::CSourceNothing(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceNothing::~CSourceNothing(void) { }

CSourcePieceWise_TurbSA::CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
		CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	incompressible = config->GetIncompressible();
	//transition     = (config->GetKind_Trans_Model() == LM);
  transition = false; // Debugging, -AA
  rotating_frame = config->GetRotating_Frame();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Spalart-Allmaras closure constants ---*/
	cv1_3 = pow(7.1,3.0);
	k2 = pow(0.41,2.0);
	cb1 = 0.1355;
	cw2 = 0.3;
	cw3_6 = pow(2.0,6.0);
	sigma = 2./3.;
	cb2 = 0.622;
	cw1 = cb1/k2+(1+cb2)/sigma;

	/*--- LM transition model constants ---*/
	beta = 0.5;
	s1   = 2.0;

}

CSourcePieceWise_TurbSA::~CSourcePieceWise_TurbSA(void) { }

void CSourcePieceWise_TurbSA::SetResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CSourcePieceWise_TurbSA::SetResidual
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i **PrimVar_Grad_i Laminar_Viscosity_i *TurbVar_i **TurbVar_Grad_i
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE dist_i cv1_3 k2 cb1 cb2 cw1 cw2 cw3_6 Volume sigma TURB_EPS
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS INT SCALAR iDim
	//SU2_CPP2C VARS DOUBLE SCALAR Density_i DivVelocity Vorticity dist_i_2
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
  val_Jacobian_i[0][0] = 0.0;
	//SU2_CPP2C COMMENT END

	/*--- Computation of vorticity ---*/
	Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
	if (nDim == 3) Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) +
                               (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );
	Omega = max(sqrt(Vorticity), 1.0e-10);
	dist_i = max(dist_i, 1.0e-10);

  /*--- Rotational correction term ---*/
  if (rotating_frame) {
    div = PrimVar_Grad_i[1][0] + PrimVar_Grad_i[2][1];
    if (nDim == 3) div += PrimVar_Grad_i[3][2];
    StrainMag = 0.0;
    // add diagonals
    StrainMag += pow(PrimVar_Grad_i[1][0] - 1.0/3.0*div,2.0);
    StrainMag += pow(PrimVar_Grad_i[2][1] - 1.0/3.0*div,2.0);
    if (nDim == 3) StrainMag += pow(PrimVar_Grad_i[3][2] - 1.0/3.0*div,2.0);
    // add off diagonals
    StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[1][1]+PrimVar_Grad_i[2][0]),2.0);
    if (nDim == 3) {
      StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[1][2]+PrimVar_Grad_i[3][0]),2.0);
      StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[2][2]+PrimVar_Grad_i[3][1]),2.0);
    }
    StrainMag = sqrt(2.0*StrainMag);
    Omega += 2.0*min(0.0,StrainMag-Omega);
  }
  
	if (dist_i > 0.0) {
    
		/*--- Production term ---*/
		dist_i_2 = dist_i*dist_i;
		nu = Laminar_Viscosity_i/Density_i;
		Ji = TurbVar_i[0]/nu;
		Ji_2 = Ji*Ji;
		Ji_3 = Ji_2*Ji;
		fv1 = Ji_3/(Ji_3+cv1_3);
		fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);
    
		Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    inv_Shat = 1.0/max(Shat, 1.0e-10);
    
    /*--- Production term ---*/
		if (!transition) val_residual[0] += cb1*Shat*TurbVar_i[0]*Volume;
    else val_residual[0] += cb1*Shat*TurbVar_i[0]*Volume*intermittency;
    
		/*--- Destruction term ---*/
		r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
		g = r + cw2*(pow(r,6.0)-r);
		g_6 =	pow(g,6.0);
		glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
		fw = g*glim;
    
		if (!transition) val_residual[0] -= cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
		else val_residual[0] -= cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume*min(max(intermittency,0.1),1.0);
    
		/*--- Diffusion term ---*/
		norm2_Grad = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
		val_residual[0] += cb2/sigma*norm2_Grad*Volume;
    
		//SU2_CPP2C COMMENT START
    
		/*--- Implicit part ---*/
    
    /*--- Production term ---*/
    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
    val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
    
    /*--- Destruction term ---*/
    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    if (r == 10.0) dr = 0.0;
    dg = dr*(1.+cw2*(6.*pow(r,5.)-1.));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] +	2.*fw)*TurbVar_i[0]/dist_i_2*Volume;
    
		//SU2_CPP2C COMMENT END
	}
	//SU2_CPP2C COMMENT START
	//SU2_CPP2C COMMENT END
  
	//SU2_CPP2C END CSourcePieceWise_TurbSA::SetResidual
  
}

CSourcePieceWise_TransLM::CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
    
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
    
	/*--- Spalart-Allmaras closure constants ---*/
	cv1_3 = pow(7.1,3.0);
	k2 = pow(0.41,2.0);
	cb1 = 0.1355;
	cw2 = 0.3;
	cw3_6 = pow(2.0,6.0);
	sigma = 2./3.;
	cb2 = 0.622;
	cw1 = cb1/k2+(1+cb2)/sigma;
    
	/*-- Gamma-theta closure constants --*/
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
  
  /*-- For debugging -AA --*/
  debugme = 0;
}

CSourcePieceWise_TransLM::~CSourcePieceWise_TransLM(void) { }

void CSourcePieceWise_TransLM::SetResidual_TransLM(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config, double &gamma_sep) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CSourcePieceWise_TransLM::SetResidual_TransLM
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *TransVar_i
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *U_i **PrimVar_Grad_i Laminar_Viscosity_i Eddy_Viscosity_i dist_i
  //SU2_CPP2C VARS DOUBLE SCALAR c_a1 c_e1 c_a2 c_e2 c_theta alpha_global flen_global
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS DOUBLE SCALAR Vorticity 
	//SU2_CPP2C DECL_LIST END

	/*-- Local intermediate variables --*/
	double rey_tc, flen, re_v, strain, f_onset1,f_onset2,f_onset3,f_onset,f_turb,tu;
    
	double prod, des;
	double f_lambda, re_theta, rey, re_theta_lim, r_t, rey_t, mach;
	double Velocity_Mag = 0.0, du_ds, theta, lambda, time_scale, delta_bl, delta, f_wake, var1, f_theta;
	double theta_bl, f_reattach;
	double dU_dx, dU_dy, dU_dz;

	//SU2_CPP2C COMMENT START
  double val_residuald[2], TransVar_id[2];

	//SU2_CPP2C COMMENT END
    
	val_residual[0] = 0.0;
	val_residual[1] = 0.0;
  
	//SU2_CPP2C COMMENT START
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	if (implicit) {
		val_Jacobian_i[0][0] = 0.0;
		val_Jacobian_i[1][0] = 0.0;
		val_Jacobian_i[0][1] = 0.0;
		val_Jacobian_i[1][1] = 0.0;
	}
	//SU2_CPP2C COMMENT END
  
  /* -- These lines included just so Tapenade doesn't complain --*/
  rey  = 0.0;
  mach = 0.0;
  tu   = 0.0;
	//SU2_CPP2C COMMENT START
  /* -- These lines must be manually reinserted into the differentiated routine! --*/
  rey  = config->GetReynolds();
  mach = config->GetMach_FreeStreamND();
	tu   = config->GetTurbulenceIntensity_FreeStream();
	//SU2_CPP2C COMMENT END

	/*--- Compute vorticity and strain (TODO: Update for 3D) ---*/
  Vorticity = fabs(PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0]);
  
  /*-- Strain = sqrt(2*Sij*Sij) --*/
	strain = sqrt(2.*(    PrimVar_Grad_i[1][0]*PrimVar_Grad_i[1][0]
                     +  0.5*pow(PrimVar_Grad_i[1][1]+PrimVar_Grad_i[2][0],2)
                     +  PrimVar_Grad_i[2][1]*PrimVar_Grad_i[2][1]  ));
    
	/*-- Note: no incompressible for now! --*/
    
	if (dist_i > 0.0) {   // Only operate away from wall
        
		/*-- Intermittency eq.: --*/
    
    rey_tc = (4.45*pow(tu,3) - 5.7*pow(tu,2) + 1.37*tu + 0.585)*TransVar_i[1];
    flen   = 0.171*pow(tu,2) - 0.0083*tu + 0.0306;

		re_v   = U_i[0]*pow(dist_i,2.)/Laminar_Viscosity_i*strain;  // Vorticity Reynolds number
        
    /*-- f_onset controls transition onset location --*/
		r_t      = Eddy_Viscosity_i/Laminar_Viscosity_i;
		f_onset1 = re_v / (2.193*rey_tc);
		f_onset2 = min(max(f_onset1, pow(f_onset1,4.)), 2.);
		f_onset3 = max(1. - pow(0.4*r_t,3),0.);
		f_onset  = max(f_onset2 - f_onset3, 0.);
        
		f_turb = exp(-pow(0.25*r_t,4));  // Medida eq. 10
        
		prod = flen*c_a1*U_i[0]*strain*sqrt(f_onset*TransVar_i[0]);
		prod = prod*(1. - c_e1*TransVar_i[0]);
        
		des = c_a2*U_i[0]*Vorticity*TransVar_i[0]*f_turb;
		des = des*(c_e2*TransVar_i[0] - 1.);
        
		val_residual[0] = prod - des;
        
		/*-- REtheta eq: --*/
		if (nDim==2) {
			Velocity_Mag = sqrt(U_i[1]*U_i[1]+U_i[2]*U_i[2])/U_i[0];
		} else if (nDim==3) {
			Velocity_Mag = sqrt(U_i[1]*U_i[1]+U_i[2]*U_i[2]+U_i[3]*U_i[3])/U_i[0];
		}
        
		/*-- Gradient of velocity magnitude ---*/
		dU_dx = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][0]
                              +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][0]);
		if (nDim==3)
			dU_dx += 0.5*Velocity_Mag*( 2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][0]);
        
		dU_dy = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][1]
                              +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][1]);
		if (nDim==3)
			dU_dy += 0.5*Velocity_Mag*( 2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][1]);
        
		if (nDim==3)
			dU_dz = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][2]
                                      +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][2]
                                      +2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][2]);
        
		du_ds = U_i[1]/(U_i[0]*Velocity_Mag) * dU_dx +  // Streamwise velocity derivative
        U_i[2]/(U_i[0]*Velocity_Mag) * dU_dy;
		if (nDim==3)
			du_ds += U_i[3]/(U_i[0]*Velocity_Mag) * dU_dz;
        
		re_theta_lim = 20.;
        
		/*-- Fixed-point iterations to solve REth correlation --*/
		f_lambda = 1.;
		for (int iter=0; iter<10; iter++) {
			if (tu <= 1.3) {
				re_theta = f_lambda * (1173.51-589.428*tu+0.2196/(tu*tu));
			} else {
				re_theta = 331.5 * f_lambda*pow(tu-0.5658,-0.671);
			}
			re_theta = max(re_theta, re_theta_lim);
            
			theta  = re_theta * Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);
            
			lambda = U_i[0]*theta*theta*du_ds / Laminar_Viscosity_i;
			lambda = min(max(-0.1,lambda),0.1);
            
			if (lambda<=0.0) {
				f_lambda = 1. - (-12.986*lambda - 123.66*lambda*lambda -
                                 405.689*lambda*lambda*lambda)*exp(-pow(2./3*tu,1.5));
			} else {
				f_lambda = 1. + 0.275*(1.-exp(-35.*lambda))*exp(-2.*tu);
			}
		}

		/*-- Calculate blending function f_theta --*/
		time_scale = 500.0*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag*Velocity_Mag);

    // Deactivated the f_wake parameter... 
		//theta_bl   = TransVar_i[1]*Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);
		//delta_bl   = 7.5*theta_bl;
		//delta      = 50.0*Vorticity*dist_i/Velocity_Mag*delta_bl + 1e-20;
    //    
		//f_wake = 1.;
        
		var1 = (TransVar_i[0]-1./c_e2)/(1.0-1./c_e2);
		var1 = 1. - pow(var1,2);
    
		//f_theta = min(max(f_wake*exp(-pow(dist_i/delta,4)), var1),1.0);
		f_theta = min(var1,1.0);
        
		val_residual[1] = c_theta*U_i[0]/time_scale *  (1.-f_theta) * (re_theta-TransVar_i[1]);
	
	  //SU2_CPP2C COMMENT START
    cout << "val_res0: "  << val_residual[0]      << endl;
    cout << "val_res1: "  << val_residual[1]      << endl;
    cout << "dist_i: "    << dist_i               << endl;
    cout << "re_v: "      << re_v                 << endl;
    cout << "c_a1: "      << c_a1                 << endl;
    cout << "strain: "    << strain               << endl;
    cout << "primgrad10: "<< PrimVar_Grad_i[1][0] << endl;
    cout << "primgrad11: "<< PrimVar_Grad_i[1][1] << endl;
    cout << "primgrad20: "<< PrimVar_Grad_i[2][0] << endl;
    cout << "primgrad21: "<< PrimVar_Grad_i[2][1] << endl;
    cout << "f_onset: "   << f_onset              << endl;
    cout << "TransVar0: " << TransVar_i[0]        << endl;
    cout << "prod: "      << prod                 << endl;
    cout << "c_a2: "      << c_a2                 << endl;
    cout << "Vorticity: " << Vorticity            << endl;
    cout << "f_turb: "    << f_turb               << endl;
    cout << "des: "       << des                  << endl;
    cout << "du_ds: "     << du_ds                << endl;
    cout << "r_t:    "    << r_t                  << endl;
    cout << "rey_tc: "    << rey_tc               << endl;
    cout << "re_theta: "  << re_theta             << endl;
        
		/*-- Calculate term for separation correction --*/
		f_reattach = exp(-pow(0.05*r_t,4));
		gamma_sep = s1*max(0.,re_v/(3.235*rey_tc)-1.)*f_reattach;
		gamma_sep = min(gamma_sep,2.0)*f_theta;
        
		/*--- Implicit part ---*/
    TransVar_id[0] = 1.0; TransVar_id[1] = 0.0;
    CSourcePieceWise_TransLM__SetResidual_TransLM_d(TransVar_i, TransVar_id, val_residual, val_residuald, config);
    val_Jacobian_i[0][0] = val_residuald[0];
    val_Jacobian_i[1][0] = val_residuald[1];

    TransVar_id[0] = 0.0; TransVar_id[1] = 1.0;
    CSourcePieceWise_TransLM__SetResidual_TransLM_d(TransVar_i, TransVar_id, val_residual, val_residuald, config);
    val_Jacobian_i[0][1] = val_residuald[0];
    val_Jacobian_i[1][1] = val_residuald[1];

	  //SU2_CPP2C COMMENT END
	}
    //SU2_CPP2C END CSourcePieceWise_TransLM::SetResidual_TransLM
}


void CSourcePieceWise_TransLM::CSourcePieceWise_TransLM__SetResidual_TransLM_d(double *TransVar_i, double *TransVar_id, double *val_residual, double *val_residuald, CConfig *config)
{
    double rey_tc, flen, re_v, strain, f_onset1, f_onset2, f_onset3, f_onset, 
    f_turb, tu;
    double rey_tcd, f_onset1d, f_onset2d, f_onsetd;
    double prod, des;
    double prodd, desd;
    double f_lambda, re_theta, rey, re_theta_lim, r_t, rey_t, mach;
    double Velocity_Mag = 0.0, du_ds, theta, lambda, time_scale, delta_bl, 
    delta, f_wake, var1, f_theta;
    double var1d, f_thetad;
    double theta_bl, f_reattach;
    double dU_dx, dU_dy, dU_dz;
    double result1;
    double result1d;
    double arg1;
    double arg1d;
    double result2;
    double x2;
    double x1;
    double x1d;
    double y1;
    double y1d;
    val_residuald[0] = 0.0;
    val_residual[0] = 0.0;
    val_residuald[1] = 0.0;
    val_residual[1] = 0.0;
    /* -- These lines included just so Tapenade doesn't complain --*/
    rey = 0.0;
    mach = 0.0;
    tu = 0.0;
    rey  = config->GetReynolds();
    mach = config->GetMach_FreeStreamND();
  	tu   = config->GetTurbulenceIntensity_FreeStream();
    /*--- Compute vorticity and strain (TODO: Update for 3D) ---*/
    Vorticity = fabs(PrimVar_Grad_i[1][1] - PrimVar_Grad_i[2][0]);
    /*-- Strain = sqrt(2*Sij*Sij) --*/
    result1 = pow(PrimVar_Grad_i[1][1] + PrimVar_Grad_i[2][0], 2);
    arg1 = 2.*(PrimVar_Grad_i[1][0]*PrimVar_Grad_i[1][0]+0.5*result1+
        PrimVar_Grad_i[2][1]*PrimVar_Grad_i[2][1]);
    strain = sqrt(arg1);
    /*-- Note: no incompressible for now! --*/
    if (dist_i > 0.0) {
        /*-- Intermittency eq.: --*/
        // Only operate away from wall
        result1 = pow(tu, 3);
        result2 = pow(tu, 2);
        rey_tcd = (4.45*result1-5.7*result2+1.37*tu+0.585)*TransVar_id[1];
        rey_tc = (4.45*result1-5.7*result2+1.37*tu+0.585)*TransVar_i[1];
        result1 = pow(tu, 2);
        flen = 0.171*result1 - 0.0083*tu + 0.0306;
        result1 = pow(dist_i, 2.);
        re_v = U_i[0]*result1/Laminar_Viscosity_i*strain;
        /*-- f_onset controls transition onset location --*/
        // Vorticity Reynolds number
        r_t = Eddy_Viscosity_i/Laminar_Viscosity_i;
        f_onset1d = -(re_v*2.193*rey_tcd/(2.193*rey_tc*(2.193*rey_tc)));
        f_onset1 = re_v/(2.193*rey_tc);
        y1d = pow_d(f_onset1, f_onset1d, 4., &y1);
        if (f_onset1 < y1) {
            x1d = y1d;
            x1 = y1;
        } else {
            x1d = f_onset1d;
            x1 = f_onset1;
        }
        if (x1 > 2.) {
            f_onset2 = 2.;
            f_onset2d = 0.0;
        } else {
            f_onset2d = x1d;
            f_onset2 = x1;
        }
        result1 = pow(0.4*r_t, 3);
        x2 = 1. - result1;
        if (x2 < 0.)
            f_onset3 = 0.;
        else
            f_onset3 = x2;
        if (f_onset2 - f_onset3 < 0.) {
            f_onset = 0.;
            f_onsetd = 0.0;
        } else {
            f_onsetd = f_onset2d;
            f_onset = f_onset2 - f_onset3;
        }
        result1 = pow(0.25*r_t, 4);
        f_turb = exp(-result1);
        // Medida eq. 10
        arg1d = f_onsetd*TransVar_i[0] + f_onset*TransVar_id[0];
        arg1 = f_onset*TransVar_i[0];
        result1d = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
        result1 = sqrt(arg1);
        prodd = flen*c_a1*U_i[0]*strain*result1d;
        prod = flen*c_a1*U_i[0]*strain*result1;
        prodd = prodd*(1.-c_e1*TransVar_i[0]) - prod*c_e1*TransVar_id[0];
        prod = prod*(1.-c_e1*TransVar_i[0]);
        desd = c_a2*U_i[0]*Vorticity*f_turb*TransVar_id[0];
        des = c_a2*U_i[0]*Vorticity*TransVar_i[0]*f_turb;
        desd = desd*(c_e2*TransVar_i[0]-1.) + des*c_e2*TransVar_id[0];
        des = des*(c_e2*TransVar_i[0]-1.);
        val_residuald[0] = prodd - desd;
        val_residual[0] = prod - des;
        /*-- REtheta eq: --*/
        if (nDim == 2) {
            arg1 = U_i[1]*U_i[1] + U_i[2]*U_i[2];
            result1 = sqrt(arg1);
            Velocity_Mag = result1/U_i[0];
        } else
            if (nDim == 3) {
                arg1 = U_i[1]*U_i[1] + U_i[2]*U_i[2] + U_i[3]*U_i[3];
                result1 = sqrt(arg1);
                Velocity_Mag = result1/U_i[0];
            }
        /*-- Gradient of velocity magnitude ---*/
        dU_dx = 0.5*Velocity_Mag*(2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][0]+2*U_i[2
            ]/U_i[0]*PrimVar_Grad_i[2][0]);
        if (nDim == 3)
            dU_dx += 0.5*Velocity_Mag*(2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][0]);
        dU_dy = 0.5*Velocity_Mag*(2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][1]+2*U_i[2
            ]/U_i[0]*PrimVar_Grad_i[2][1]);
        if (nDim == 3)
            dU_dy += 0.5*Velocity_Mag*(2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][1]);
        if (nDim == 3)
            dU_dz = 0.5*Velocity_Mag*(2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][2]+2*
                U_i[2]/U_i[0]*PrimVar_Grad_i[2][2]+2*U_i[3]/U_i[0]*
                PrimVar_Grad_i[3][2]);
        du_ds = U_i[1]/(U_i[0]*Velocity_Mag)*dU_dx + U_i[2]/(U_i[0]*
            Velocity_Mag)*dU_dy;
        // Streamwise velocity derivative
        if (nDim == 3)
            du_ds += U_i[3]/(U_i[0]*Velocity_Mag)*dU_dz;
        re_theta_lim = 20.;
        /*-- Fixed-point iterations to solve REth correlation --*/
        f_lambda = 1.;
        {
          double x3;
          for (int iter = 0; iter < 10; ++iter) {
              if (tu <= 1.3)
                  re_theta = f_lambda*(1173.51-589.428*tu+0.2196/(tu*tu));
              else {
                  result1 = pow(tu - 0.5658, -0.671);
                  re_theta = 331.5*f_lambda*result1;
              }
              if (re_theta < re_theta_lim)
                  re_theta = re_theta_lim;
              else
                  re_theta = re_theta;
              theta = re_theta*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag);
              lambda = U_i[0]*theta*theta*du_ds/Laminar_Viscosity_i;
              if (-0.1 < lambda)
                  x3 = lambda;
              else
                  x3 = -0.1;
              if (x3 > 0.1)
                  lambda = 0.1;
              else
                  lambda = x3;
              if (lambda <= 0.0) {
                  result1 = pow(2./3*tu, 1.5);
                  f_lambda = 1. - (-12.986*lambda-123.66*lambda*lambda-405.689
                      *lambda*lambda*lambda)*exp(-result1);
              } else
                  f_lambda = 1. + 0.275*(1.-exp(-35.*lambda))*exp(-2.*tu);
          }
        }
        /*-- Calculate blending function f_theta --*/
        time_scale = 500.0*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag*
            Velocity_Mag);
        // Deactivated the f_wake parameter... 
        //theta_bl   = TransVar_i[1]*Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);
        //delta_bl   = 7.5*theta_bl;
        //delta      = 50.0*Vorticity*dist_i/Velocity_Mag*delta_bl + 1e-20;
        //    
        //f_wake = 1.;
        var1d = TransVar_id[0]/(1.0-1./c_e2);
        var1 = (TransVar_i[0]-1./c_e2)/(1.0-1./c_e2);
        result1d = pow_d(var1, var1d, 2, &result1);
        var1d = -result1d;
        var1 = 1. - result1;
        if (var1 > 1.0) {
            f_theta = 1.0;
            f_thetad = 0.0;
        } else {
            f_thetad = var1d;
            f_theta = var1;
        }
        val_residuald[1] = c_theta*U_i[0]*(-(f_thetad*(re_theta-TransVar_i[1])
            )-(1.-f_theta)*TransVar_id[1])/time_scale;
        val_residual[1] = c_theta*U_i[0]/time_scale*(1.-f_theta)*(re_theta-
            TransVar_i[1]);
    } else
        *val_residuald = 0.0;

}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double *constants,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
    
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
  val_Jacobian_i[0][0] = 0.0;		val_Jacobian_i[0][1] = 0.0;
  val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = 0.0;
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
    val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;		val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;                               val_Jacobian_i[1][1] = -2.0*beta_blended*TurbVar_i[1]*Volume;
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
  else {
     
		/*--- Evaluate the source term  ---*/
		val_residual[nDim] = Volume * U_i[0] * STANDART_GRAVITY;
    
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
            
			rho_Pos  = 1.0/3.0*(U_0[0] + U_1[0] + U_2[0]);
			rho_Elec = 1.0/3.0*(U_0[1] + U_1[1] + U_2[1]);
            
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
            
//            if (fabs(diff_ru_Elec_x) > 1E-5) {
//                cout << " alpha = " << alpha << " dt = " << dt << " Kai_np1 = " << Kai_np1 << " Kai_n = " << Kai_n;
//                cout << " n_rho = " << rho_Elec/mass_Elec - rho_Pos/mass_Pos << endl;
//                cout << " diff_ru_Elec_x = " << diff_ru_Elec_x;
//                cout << " diff_ru_Elec_y = " << diff_ru_Elec_y << endl;
//                cout << " source = " << val_residual[0] << ", " << val_residual[1] << ", " << val_residual[2] << endl;
//            }

		}
	}
	if (config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == ARGON_SID || config->GetKind_GasModel() == O2 || config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {
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
	double Gas_Constant = config->GetGas_ConstantND();

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

	/*--- Turn on laminar viscosity sensitivity for NS ---*/
	if ((!config->GetFrozen_Visc()) && (config->GetKind_Solver() != ADJ_RANS)) {    

		double Temperature_Ref = config->GetTemperature_Ref();
		double Temperature_Dim = Temp_i*Temperature_Ref;
		double dVisc_T = ((Laminar_Viscosity_i)/(2.0*Temperature_Dim*(Temperature_Dim + 110.3)))*(Temperature_Dim + 3.0*110.3)*Temperature_Ref;

		double Cp = (Gamma/Gamma_Minus_One)*Gas_Constant;
		double kappa_psi = (sigma_gradpsi + vel_sigma_gradpsi5)/mu_tot_1;
		double theta = (kappa_psi + Cp/PRANDTL*gradT_gradpsi5)*dVisc_T*Gamma_Minus_One/(Gas_Constant*Density);


        val_residual[0] += (theta*(sq_vel-Energy))*Volume;
        for (iDim = 0; iDim < nDim; iDim++)
            val_residual[iDim+1] -= theta*Velocity[iDim]*Volume;
        val_residual[nVar-1] += theta*Volume;


	}
    
    /*--- Turn on laminar/eddy viscosity sensitivity for Hybrid RANS ---*/
	if ((config->GetKind_Adjoint() == HYBRID) && (config->GetKind_Solver() == ADJ_RANS)) {

		double Temperature_Ref = config->GetTemperature_Ref();
		double Temperature_Dim = Temp_i*Temperature_Ref;
		double dVisc_T = ((Laminar_Viscosity_i)/(2.0*Temperature_Dim*(Temperature_Dim + 110.3)))*(Temperature_Dim + 3.0*110.3)*Temperature_Ref;
        
		double Cp = (Gamma/Gamma_Minus_One)*Gas_Constant;
		double kappa_psi = (sigma_gradpsi + vel_sigma_gradpsi5)/mu_tot_1 + Cp/PRANDTL_TURB*gradT_gradpsi5;
		double theta = (kappa_psi + Cp/PRANDTL*gradT_gradpsi5)*dVisc_T*Gamma_Minus_One/(Gas_Constant*Density);
        
        // If frozen hybrid, this doesn't get added
        if (!config->GetFrozen_Visc()) {
            val_residual[0] += (theta*(sq_vel-Energy))*Volume;
            for (iDim = 0; iDim < nDim; iDim++)
                val_residual[iDim+1] -= theta*Velocity[iDim]*Volume;
            val_residual[nVar-1] += theta*Volume;
        }
        
		// store this value for coupling
        kappapsi_Volume = kappa_psi*Volume;
            
        SetKappaPsiVolume(kappapsi_Volume);
        
	}

	/*--- Coupling terms coming from the continuous adjoint turbulent equations ---*/
	if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc()) && (config->GetKind_Adjoint() == CONTINUOUS)) {
    
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
		double dVisc_T;
		dVisc_T = 0.0;
        
		if (!config->GetFrozen_Visc()) {
      
			double Temperature_Ref = config->GetTemperature_Ref();
			double Temperature_Dim = Temp_i*Temperature_Ref;
			dVisc_T = ((Laminar_Viscosity_i)/(2.0*Temperature_Dim*(Temperature_Dim + 110.3)))*(Temperature_Dim + 3.0*110.3)*Temperature_Ref;
      
		}
    
		//		double mu1 = 1.404/config->GetReynolds();
		//		double mu2 = 0.404;
		//		double dVisc_T = Laminar_Viscosity_i*(Temp_i+3.0*mu2)/(2.0*Temp_i*(Temp_i+mu2));

		double Cp = (Gamma/Gamma_Minus_One)*Gas_Constant;
		double kappa_psi = (sigma_gradpsi + vel_sigma_gradpsi5)/mu_tot_1 + Cp/PRANDTL_TURB*gradT_gradpsi5;
		double cv1_const = 3.0*cv1_3/(Ji_3+cv1_3);
		double theta = (kappa_psi*(1.0-Eddy_Viscosity_i/Laminar_Viscosity_i*cv1_const) -
				Cp/PRANDTL_TURB*gradT_gradpsi5*(1.0-PRANDTL_TURB/PRANDTL))*dVisc_T*Gamma_Minus_One/(Gas_Constant*Density);
		double xi = kappa_psi*(1.0+cv1_const)*Eddy_Viscosity_i/Density;

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

			double alpha_coeff = Gamma_Minus_One/(Gas_Constant*Density)*dVisc_T;
			double beta_coeff = alpha_coeff*(sq_vel-Energy)-Laminar_Viscosity_i/Density;
			double Fs_coeff = TurbPsi_i[0]*(cb1*TurbVar_i[0]-cw1*TurbVar_i[0]*TurbVar_i[0]/dist_0_2*dfw_g*dg_r*dr_Shat)*
					dShat_fv2*(dfv2_Ji+dfv2_fv1*dfv1_Ji)*dJi_nu;
			double Gamma = Fs_coeff - gradTurbVar_gradTurbPsi/sigma;

			val_residual[0] -= (Gamma*beta_coeff - TurbVar_i[0]*vel_gradTurbPsi)/Density*Volume;
			for (iDim = 0; iDim < nDim; iDim++)
				val_residual[iDim+1] += (Gamma*alpha_coeff*Velocity[iDim] - TurbVar_i[0]*TurbPsi_Grad_i[0][iDim])/Density*Volume;
			val_residual[nVar-1] -= (Gamma*alpha_coeff)/Density*Volume;

            // this should improve stability (when commented):
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

CSourcePieceWise_AdjDiscTurb::CSourcePieceWise_AdjDiscTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
    
}

CSourcePieceWise_AdjDiscTurb::~CSourcePieceWise_AdjDiscTurb(void) {
    
}

void CSourcePieceWise_AdjDiscTurb::SetResidual() {
    
}

CSourcePieceWise_AdjDiscTurbSA::CSourcePieceWise_AdjDiscTurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	/*--- Closure constants ---*/
	cv1_3 = pow(7.1,3.0);
	k2 = pow(0.41,2.0);
	cb1 = 0.1355;
	cw2 = 0.3;
	cw3_6 = pow(2.0,6.0);
	sigma = 2./3.;
	cb2 = 0.622;
	cw1 = cb1/k2+(1+cb2)/sigma;
    
    TurbVar_id = new double[nVar];
    val_residual = new double [nVar];
    val_residuald = new double [nVar];
    
    unsigned short nFlowVar = (nDim+2);
    PrimVar_Grad_id = new double *[nFlowVar];
    for (unsigned short iPos = 0; iPos < nFlowVar; iPos++)
        PrimVar_Grad_id[iPos] = new double [nDim];
    
    TurbVar_Grad_id = new double *[nVar];
    for (unsigned short iPos = 0; iPos < nVar; iPos++)
        TurbVar_Grad_id[iPos] = new double [nDim];
    
    U_id = new double [nFlowVar];
	U_jd = new double [nFlowVar];


}

CSourcePieceWise_AdjDiscTurbSA::~CSourcePieceWise_AdjDiscTurbSA(void) {
    
    delete [] TurbVar_i;
    delete [] val_residual;
    delete [] val_residuald;
    
    unsigned short nFlowVar = (nDim+2);
    for (unsigned short iPos = 0; iPos < nFlowVar; iPos++)
        delete [] PrimVar_Grad_id[iPos];
    delete [] PrimVar_Grad_id;
    
    for (unsigned short iPos = 0; iPos < nVar; iPos++)
        delete [] TurbVar_Grad_id[iPos];
    delete [] TurbVar_Grad_id;
    
    delete [] U_id;
	delete [] U_jd;
    

}

void CSourcePieceWise_AdjDiscTurbSA::SetResidual(double **val_Jacobian_i, double *val_Jacobian_mui, double ***val_Jacobian_gradi, CConfig *config) {
	//unsigned short iVar, jVar;
	unsigned short iPos, jPos, kPos, lPos; // iVar, jVar are global and used in SetResidual_ad, so cannot be used here

	unsigned short nTotalVar, nFlowVar;
	nFlowVar = nDim + 2;
	nTotalVar = nVar + nFlowVar;

	// U_i sensitivity:

	for (iPos = 0; iPos < nFlowVar; iPos++){
// zero things first
		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		Laminar_Viscosity_id = 0.0;

		for (jPos = 0; jPos < nFlowVar; jPos++){
			for (kPos = 0; kPos < nDim; kPos++){
				PrimVar_Grad_id[jPos][kPos] = 0.0;
			}
		}

		for (jPos = 0; jPos < nVar; jPos++){
			for (kPos = 0; kPos < nDim; kPos++){
				TurbVar_Grad_id[jPos][kPos] = 0.0;
			}
		}

		U_id[iPos] = 1.0;

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			// Transpose each block: [jPos][iPos] -> [iPos][jPos]
			val_Jacobian_i[iPos][jPos] = val_residuald[jPos];

		}
        
//        // TEST AD *****************
//        for (jPos = 0; jPos < nVar; jPos++) {
//            //cout << "Dist: " << dist_i << endl;
//            //cout << "Direct: " << val_residual[jPos] << endl;
//            cout << "AD: " << val_residuald[jPos] << endl;
//        }
//       // cout << "--" << endl;
//        // AD *****************


	}
    
//    // TEST FD *****************
//    double *temp_U_i;
//    double *temp_TurbVar_i;
//    double temp_Laminar_Viscosity_i;
//    temp_U_i = new double[nFlowVar];
//    temp_TurbVar_i  = new double[nVar];
//    double *temp1_val_residual, *temp2_val_residual;
//    temp1_val_residual = new double[nVar];
//    temp2_val_residual = new double[nVar];
//    double **temp_PrimVar_Grad_i, **temp_TurbVar_Grad_i;
//    temp_PrimVar_Grad_i = new double*[nFlowVar];
//    temp_TurbVar_Grad_i = new double*[nVar];
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_PrimVar_Grad_i[jPos] = new double[nDim];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_Grad_i[jPos] = new double[nDim];
//    }
//    
//    double delta;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_U_i[jPos] = U_i[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//    }
//    
//    temp_Laminar_Viscosity_i = Laminar_Viscosity_i;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_PrimVar_Grad_i[jPos][kPos] = PrimVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_TurbVar_Grad_i[jPos][kPos] = TurbVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (iPos = 0; iPos < nFlowVar; iPos++){
//        
//		for (jPos = 0; jPos < nFlowVar; jPos++){
//			U_i[jPos] = temp_U_i[jPos];
//		}
//        
//		for (jPos = 0; jPos < nVar; jPos++){
//			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//		}
//        
//        Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//            for (kPos = 0; kPos < nDim; kPos++){
//                PrimVar_Grad_i[jPos][kPos] = temp_PrimVar_Grad_i[jPos][kPos];
//            }
//        }
//        
//        for (jPos = 0; jPos < nVar; jPos++){
//            for (kPos = 0; kPos < nDim; kPos++){
//                TurbVar_Grad_i[jPos][kPos] = temp_TurbVar_Grad_i[jPos][kPos];
//            }
//        }
//        
//        if (fabs(temp_U_i[iPos]) > 1e-15)
//            delta = 0.01*temp_U_i[iPos];
//        else
//            delta = 1e-15;
//        
//        U_i[iPos] = temp_U_i[iPos] - delta;
//        
//		this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//			U_i[jPos] = temp_U_i[jPos];
//		}
//        
//		for (jPos = 0; jPos < nVar; jPos++){
//			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//		}
//        
//        Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//        
//        U_i[iPos] = temp_U_i[iPos] + delta;
//        
//		this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//        
//        //cout << "U_i: " << temp_U_i[iPos] << endl;
//        
//        for (jPos = 0; jPos < nVar; jPos++) {
//            cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//        }
//	}
//    delete [] temp_U_i;
//    delete [] temp_TurbVar_i;
//    delete [] temp1_val_residual;
//    delete [] temp2_val_residual;
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        delete [] temp_PrimVar_Grad_i[jPos];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        delete [] temp_TurbVar_Grad_i[jPos];
//    }
//    delete [] temp_PrimVar_Grad_i;
//    delete [] temp_TurbVar_Grad_i;
//    cin.get();
//    // FD *****************


	// TurbVar_i sensitivity

	for (iPos = 0; iPos < nVar; iPos++){
		// zero things first
		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		Laminar_Viscosity_id = 0.0;

		for (jPos = 0; jPos < nFlowVar; jPos++){
			for (kPos = 0; kPos < nDim; kPos++){
				PrimVar_Grad_id[jPos][kPos] = 0.0;
			}
		}

		for (jPos = 0; jPos < nVar; jPos++){
			for (kPos = 0; kPos < nDim; kPos++){
				TurbVar_Grad_id[jPos][kPos] = 0.0;
			}
		}

		TurbVar_id[iPos] = 1.0;

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			// Transpose each block: [jPos][iPos] -> [iPos][jPos]
			val_Jacobian_i[iPos+nFlowVar][jPos] = val_residuald[jPos];

		}
        
//        // TEST AD *****************
//        for (jPos = 0; jPos < nVar; jPos++) {
//            //cout << "Dist: " << dist_i << endl;
//            //cout << "Direct: " << val_residual[jPos] << endl;
//            cout << "AD: " << val_residuald[jPos] << endl;
//        }
//        // cout << "--" << endl;
//        // AD *****************

	}
    
    
//    // TEST FD *****************
//    double *temp_U_i;
//    double *temp_TurbVar_i;
//    double temp_Laminar_Viscosity_i;
//    temp_U_i = new double[nFlowVar];
//    temp_TurbVar_i  = new double[nVar];
//    double *temp1_val_residual, *temp2_val_residual;
//    temp1_val_residual = new double[nVar];
//    temp2_val_residual = new double[nVar];
//    double **temp_PrimVar_Grad_i, **temp_TurbVar_Grad_i;
//    temp_PrimVar_Grad_i = new double*[nFlowVar];
//    temp_TurbVar_Grad_i = new double*[nVar];
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_PrimVar_Grad_i[jPos] = new double[nDim];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_Grad_i[jPos] = new double[nDim];
//    }
//    
//    double delta;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_U_i[jPos] = U_i[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//    }
//    
//    temp_Laminar_Viscosity_i = Laminar_Viscosity_i;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_PrimVar_Grad_i[jPos][kPos] = PrimVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_TurbVar_Grad_i[jPos][kPos] = TurbVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (iPos = 0; iPos < nVar; iPos++){
//        
//		for (jPos = 0; jPos < nFlowVar; jPos++){
//			U_i[jPos] = temp_U_i[jPos];
//		}
//        
//		for (jPos = 0; jPos < nVar; jPos++){
//			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//		}
//        
//        Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//            for (kPos = 0; kPos < nDim; kPos++){
//                PrimVar_Grad_i[jPos][kPos] = temp_PrimVar_Grad_i[jPos][kPos];
//            }
//        }
//        
//        for (jPos = 0; jPos < nVar; jPos++){
//            for (kPos = 0; kPos < nDim; kPos++){
//                TurbVar_Grad_i[jPos][kPos] = temp_TurbVar_Grad_i[jPos][kPos];
//            }
//        }
//        
//        if (fabs(TurbVar_i[iPos]) > 1e-15)
//            delta = 0.01*TurbVar_i[iPos];
//        else
//            delta = 1e-15;
//        
//        TurbVar_i[iPos] = temp_TurbVar_i[iPos] - delta;
//        
//		this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//			U_i[jPos] = temp_U_i[jPos];
//		}
//        
//		for (jPos = 0; jPos < nVar; jPos++){
//			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//		}
//        
//        Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//        
//        TurbVar_i[iPos] = temp_TurbVar_i[iPos] + delta;
//        
//		this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//        
//        //cout << "U_i: " << temp_U_i[iPos] << endl;
//        
//        for (jPos = 0; jPos < nVar; jPos++) {
//            cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//        }
//	}
//    delete [] temp_U_i;
//    delete [] temp_TurbVar_i;
//    delete [] temp1_val_residual;
//    delete [] temp2_val_residual;
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        delete [] temp_PrimVar_Grad_i[jPos];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        delete [] temp_TurbVar_Grad_i[jPos];
//    }
//    delete [] temp_PrimVar_Grad_i;
//    delete [] temp_TurbVar_Grad_i;
//    cin.get();
//    // FD *****************

	// Laminar_Viscosity_i sensitivity
		// zero things first
		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
			val_residuald[jPos] = 0.0;
		}

		Laminar_Viscosity_id = 0.0;

		for (jPos = 0; jPos < nFlowVar; jPos++){
			for (kPos = 0; kPos < nDim; kPos++){
				PrimVar_Grad_id[jPos][kPos] = 0.0;
			}
		}

		for (jPos = 0; jPos < nVar; jPos++){
			for (kPos = 0; kPos < nDim; kPos++){
				TurbVar_Grad_id[jPos][kPos] = 0.0;
			}
		}

		Laminar_Viscosity_id = 1.0;

		this->SetDirectResidual_ad(val_residual, val_residuald, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			val_Jacobian_mui[jPos] = val_residuald[jPos];
            //cout << val_residuald[jPos] << endl;
		}
    
    
    
//    // TEST AD *****************
//    for (jPos = 0; jPos < nVar; jPos++) {
//        //cout << "Dist: " << dist_i << endl;
//        //cout << "Direct: " << val_residual[jPos] << endl;
//        cout << "AD: " << val_residuald[jPos] << endl;
//    }
//    // cout << "--" << endl;
//    // AD *****************
    
//    // TEST FD *****************
//    double *temp_U_i;
//    double *temp_TurbVar_i;
//    double temp_Laminar_Viscosity_i;
//    temp_U_i = new double[nFlowVar];
//    temp_TurbVar_i  = new double[nVar];
//    double *temp1_val_residual, *temp2_val_residual;
//    temp1_val_residual = new double[nVar];
//    temp2_val_residual = new double[nVar];
//    double **temp_PrimVar_Grad_i, **temp_TurbVar_Grad_i;
//    temp_PrimVar_Grad_i = new double*[nFlowVar];
//    temp_TurbVar_Grad_i = new double*[nVar];
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_PrimVar_Grad_i[jPos] = new double[nDim];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_Grad_i[jPos] = new double[nDim];
//    }
//    
//    double delta;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_U_i[jPos] = U_i[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//    }
//    
//    temp_Laminar_Viscosity_i = Laminar_Viscosity_i;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_PrimVar_Grad_i[jPos][kPos] = PrimVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_TurbVar_Grad_i[jPos][kPos] = TurbVar_Grad_i[jPos][kPos];
//        }
//    }
//            
//		for (jPos = 0; jPos < nFlowVar; jPos++){
//			U_i[jPos] = temp_U_i[jPos];
//		}
//        
//		for (jPos = 0; jPos < nVar; jPos++){
//			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//		}
//        
//        Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//            for (kPos = 0; kPos < nDim; kPos++){
//                PrimVar_Grad_i[jPos][kPos] = temp_PrimVar_Grad_i[jPos][kPos];
//            }
//        }
//        
//        for (jPos = 0; jPos < nVar; jPos++){
//            for (kPos = 0; kPos < nDim; kPos++){
//                TurbVar_Grad_i[jPos][kPos] = temp_TurbVar_Grad_i[jPos][kPos];
//            }
//        }
//        
//        if (fabs(Laminar_Viscosity_i) > 1e-15)
//            delta = 0.01*Laminar_Viscosity_i;
//        else
//            delta = 1e-15;
//        
//        Laminar_Viscosity_i = temp_Laminar_Viscosity_i - delta;
//        
//		this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//        
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//			U_i[jPos] = temp_U_i[jPos];
//		}
//        
//		for (jPos = 0; jPos < nVar; jPos++){
//			TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//		}
//        
//        Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//        
//        Laminar_Viscosity_i = temp_Laminar_Viscosity_i + delta;
//        
//		this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//        
//        //cout << "U_i: " << temp_U_i[iPos] << endl;
//        
//        for (jPos = 0; jPos < nVar; jPos++) {
//            cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//        }
//    delete [] temp_U_i;
//    delete [] temp_TurbVar_i;
//    delete [] temp1_val_residual;
//    delete [] temp2_val_residual;
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        delete [] temp_PrimVar_Grad_i[jPos];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        delete [] temp_TurbVar_Grad_i[jPos];
//    }
//    delete [] temp_PrimVar_Grad_i;
//    delete [] temp_TurbVar_Grad_i;
//    cin.get();
//    // FD *****************
    

	// PrimVar_Grad_i sensitivity

		for (iPos = 0; iPos < nFlowVar; iPos++){
			for (lPos = 0; lPos < nDim; lPos++){
				// zero things first
				for (jPos = 0; jPos < nFlowVar; jPos++){
					U_id[jPos] = 0.0;
				}

				for (jPos = 0; jPos < nVar; jPos++){
					TurbVar_id[jPos] = 0.0;
					val_residuald[jPos] = 0.0;
				}

				Laminar_Viscosity_id = 0.0;

				for (jPos = 0; jPos < nFlowVar; jPos++){
					for (kPos = 0; kPos < nDim; kPos++){
						PrimVar_Grad_id[jPos][kPos] = 0.0;
					}
				}

				for (jPos = 0; jPos < nVar; jPos++){
					for (kPos = 0; kPos < nDim; kPos++){
						TurbVar_Grad_id[jPos][kPos] = 0.0;
					}
				}

				PrimVar_Grad_id[iPos][lPos] = 1.0;

				this->SetDirectResidual_ad(val_residual, val_residuald, config);

				for (jPos = 0; jPos < nVar; jPos++) {
					// Transpose each block: [jPos][iPos] -> [iPos][jPos]
					val_Jacobian_gradi[iPos][lPos][jPos] = val_residuald[jPos];

				}
                
//                // TEST AD *****************
//                for (jPos = 0; jPos < nVar; jPos++) {
//                    //cout << "Dist: " << dist_i << endl;
//                    //cout << "Direct: " << val_residual[jPos] << endl;
//                    cout << "AD: " << val_residuald[jPos] << endl;
//                }
//                // cout << "--" << endl;
//                // AD *****************

                
			}
		}
    
    
//    // TEST FD *****************
//    double *temp_U_i;
//    double *temp_TurbVar_i;
//    double temp_Laminar_Viscosity_i;
//    temp_U_i = new double[nFlowVar];
//    temp_TurbVar_i  = new double[nVar];
//    double *temp1_val_residual, *temp2_val_residual;
//    temp1_val_residual = new double[nVar];
//    temp2_val_residual = new double[nVar];
//    double **temp_PrimVar_Grad_i, **temp_TurbVar_Grad_i;
//    temp_PrimVar_Grad_i = new double*[nFlowVar];
//    temp_TurbVar_Grad_i = new double*[nVar];
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_PrimVar_Grad_i[jPos] = new double[nDim];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_Grad_i[jPos] = new double[nDim];
//    }
//    
//    double delta;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_U_i[jPos] = U_i[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//    }
//    
//    temp_Laminar_Viscosity_i = Laminar_Viscosity_i;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_PrimVar_Grad_i[jPos][kPos] = PrimVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_TurbVar_Grad_i[jPos][kPos] = TurbVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (iPos = 0; iPos < nFlowVar; iPos++){
//        for (lPos = 0; lPos < nDim; lPos++){
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        U_i[jPos] = temp_U_i[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//    }
//    
//    Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            PrimVar_Grad_i[jPos][kPos] = temp_PrimVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            TurbVar_Grad_i[jPos][kPos] = temp_TurbVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    if (fabs(PrimVar_Grad_i[iPos][lPos]) > 1e-15)
//        delta = 0.01*PrimVar_Grad_i[iPos][lPos];
//    else
//        delta = 1e-15;
//    
//    PrimVar_Grad_i[iPos][lPos] = temp_PrimVar_Grad_i[iPos][lPos] - delta;
//    
//    this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        U_i[jPos] = temp_U_i[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//    }
//    
//    Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//            
//            for (jPos = 0; jPos < nFlowVar; jPos++){
//                for (kPos = 0; kPos < nDim; kPos++){
//                    PrimVar_Grad_i[jPos][kPos] = temp_PrimVar_Grad_i[jPos][kPos];
//                }
//            }
//            
//            for (jPos = 0; jPos < nVar; jPos++){
//                for (kPos = 0; kPos < nDim; kPos++){
//                    TurbVar_Grad_i[jPos][kPos] = temp_TurbVar_Grad_i[jPos][kPos];
//                }
//            }
//    
//    PrimVar_Grad_i[iPos][lPos] = temp_PrimVar_Grad_i[iPos][lPos] + delta;
//    
//    this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//    
//    //cout << "U_i: " << temp_U_i[iPos] << endl;
//    
//    for (jPos = 0; jPos < nVar; jPos++) {
//        cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//    }
//            
//        }
//    }
//    
//    delete [] temp_U_i;
//    delete [] temp_TurbVar_i;
//    delete [] temp1_val_residual;
//    delete [] temp2_val_residual;
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        delete [] temp_PrimVar_Grad_i[jPos];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        delete [] temp_TurbVar_Grad_i[jPos];
//    }
//    delete [] temp_PrimVar_Grad_i;
//    delete [] temp_TurbVar_Grad_i;
//    cin.get();
//    // FD *****************

    

	// TurbVar_Grad_i sensitivity

		for (iPos = 0; iPos < nVar; iPos++){
			for (lPos = 0; lPos < nDim; lPos++){
				// zero things first
				for (jPos = 0; jPos < nFlowVar; jPos++){
					U_id[jPos] = 0.0;
				}

				for (jPos = 0; jPos < nVar; jPos++){
					TurbVar_id[jPos] = 0.0;
					val_residuald[jPos] = 0.0;
				}

				Laminar_Viscosity_id = 0.0;

				for (jPos = 0; jPos < nFlowVar; jPos++){
					for (kPos = 0; kPos < nDim; kPos++){
						PrimVar_Grad_id[jPos][kPos] = 0.0;
					}
				}

				for (jPos = 0; jPos < nVar; jPos++){
					for (kPos = 0; kPos < nDim; kPos++){
						TurbVar_Grad_id[jPos][kPos] = 0.0;
					}
				}

				TurbVar_Grad_id[iPos][lPos] = 1.0;

				this->SetDirectResidual_ad(val_residual, val_residuald, config);

				for (jPos = 0; jPos < nVar; jPos++) {
					// Transpose each block: [jPos][iPos] -> [iPos][jPos]
					val_Jacobian_gradi[iPos+nFlowVar][lPos][jPos] = val_residuald[jPos];

				}
                
//                // TEST AD *****************
//                for (jPos = 0; jPos < nVar; jPos++) {
//                    //cout << "Dist: " << dist_i << endl;
//                    //cout << "Direct: " << val_residual[jPos] << endl;
//                    cout << "AD: " << val_residuald[jPos] << endl;
//                }
//                // cout << "--" << endl;
//                // AD *****************

                
			}
		}
    
//    // TEST FD *****************
//    double *temp_U_i;
//    double *temp_TurbVar_i;
//    double temp_Laminar_Viscosity_i;
//    temp_U_i = new double[nFlowVar];
//    temp_TurbVar_i  = new double[nVar];
//    double *temp1_val_residual, *temp2_val_residual;
//    temp1_val_residual = new double[nVar];
//    temp2_val_residual = new double[nVar];
//    double **temp_PrimVar_Grad_i, **temp_TurbVar_Grad_i;
//    temp_PrimVar_Grad_i = new double*[nFlowVar];
//    temp_TurbVar_Grad_i = new double*[nVar];
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_PrimVar_Grad_i[jPos] = new double[nDim];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_Grad_i[jPos] = new double[nDim];
//    }
//    
//    double delta;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        temp_U_i[jPos] = U_i[jPos];
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        temp_TurbVar_i[jPos] = TurbVar_i[jPos];
//    }
//    
//    temp_Laminar_Viscosity_i = Laminar_Viscosity_i;
//    
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_PrimVar_Grad_i[jPos][kPos] = PrimVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (jPos = 0; jPos < nVar; jPos++){
//        for (kPos = 0; kPos < nDim; kPos++){
//            temp_TurbVar_Grad_i[jPos][kPos] = TurbVar_Grad_i[jPos][kPos];
//        }
//    }
//    
//    for (iPos = 0; iPos < nVar; iPos++){
//        for (lPos = 0; lPos < nDim; lPos++){
//            
//            for (jPos = 0; jPos < nFlowVar; jPos++){
//                U_i[jPos] = temp_U_i[jPos];
//            }
//            
//            for (jPos = 0; jPos < nVar; jPos++){
//                TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//            }
//            
//            Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//            
//            for (jPos = 0; jPos < nFlowVar; jPos++){
//                for (kPos = 0; kPos < nDim; kPos++){
//                    PrimVar_Grad_i[jPos][kPos] = temp_PrimVar_Grad_i[jPos][kPos];
//                }
//            }
//            
//            for (jPos = 0; jPos < nVar; jPos++){
//                for (kPos = 0; kPos < nDim; kPos++){
//                    TurbVar_Grad_i[jPos][kPos] = temp_TurbVar_Grad_i[jPos][kPos];
//                }
//            }
//            
//            if (fabs(TurbVar_Grad_i[iPos][lPos]) > 1e-15)
//                delta = 0.01*TurbVar_Grad_i[iPos][lPos];
//            else
//                delta = 1e-15;
//            
//            TurbVar_Grad_i[iPos][lPos] = temp_TurbVar_Grad_i[iPos][lPos] - delta;
//            
//            this->SetDirectResidual_ad(temp1_val_residual, val_residuald, config);
//            
//            for (jPos = 0; jPos < nFlowVar; jPos++){
//                U_i[jPos] = temp_U_i[jPos];
//            }
//            
//            for (jPos = 0; jPos < nVar; jPos++){
//                TurbVar_i[jPos] = temp_TurbVar_i[jPos];
//            }
//            
//            Laminar_Viscosity_i = temp_Laminar_Viscosity_i;
//            
//            for (jPos = 0; jPos < nFlowVar; jPos++){
//                for (kPos = 0; kPos < nDim; kPos++){
//                    PrimVar_Grad_i[jPos][kPos] = temp_PrimVar_Grad_i[jPos][kPos];
//                }
//            }
//            
//            for (jPos = 0; jPos < nVar; jPos++){
//                for (kPos = 0; kPos < nDim; kPos++){
//                    TurbVar_Grad_i[jPos][kPos] = temp_TurbVar_Grad_i[jPos][kPos];
//                }
//            }
//            
//            TurbVar_Grad_i[iPos][lPos] = temp_TurbVar_Grad_i[iPos][lPos] + delta;
//            
//            this->SetDirectResidual_ad(temp2_val_residual, val_residuald, config);
//            
//            //cout << "U_i: " << temp_U_i[iPos] << endl;
//            
//            for (jPos = 0; jPos < nVar; jPos++) {
//                cout << "FD: " << (temp2_val_residual[jPos] - temp1_val_residual[jPos])/(2*delta) << endl;
//            }
//            
//        }
//    }
//    
//    delete [] temp_U_i;
//    delete [] temp_TurbVar_i;
//    delete [] temp1_val_residual;
//    delete [] temp2_val_residual;
//    for (jPos = 0; jPos < nFlowVar; jPos++){
//        delete [] temp_PrimVar_Grad_i[jPos];
//    }
//    for (jPos = 0; jPos < nVar; jPos++){
//        delete [] temp_TurbVar_Grad_i[jPos];
//    }
//    delete [] temp_PrimVar_Grad_i;
//    delete [] temp_TurbVar_Grad_i;
//    cin.get();
//    // FD *****************


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

CSourceRotatingFrame_Flow::CSourceRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
    
}

CSourceRotatingFrame_Flow::~CSourceRotatingFrame_Flow(void) { }

void CSourceRotatingFrame_Flow::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {

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

CSourceRotatingFrame_AdjFlow::CSourceRotatingFrame_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceRotatingFrame_AdjFlow::~CSourceRotatingFrame_AdjFlow(void) { }

void CSourceRotatingFrame_AdjFlow::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
    
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

CSourceRotatingFrame_AdjDiscFlow::CSourceRotatingFrame_AdjDiscFlow() { }

CSourceRotatingFrame_AdjDiscFlow::~CSourceRotatingFrame_AdjDiscFlow(void) { }

void CSourceRotatingFrame_AdjDiscFlow::SetResidual() {
    
}

CSourceAxisymmetric_Flow::CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
    
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
    
}

CSourceAxisymmetric_Flow::~CSourceAxisymmetric_Flow(void) { }

void CSourceAxisymmetric_Flow::SetResidual(double *val_residual, double **Jacobian_i, CConfig *config) {
    
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

CSourceAxisymmetric_AdjFlow::CSourceAxisymmetric_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceAxisymmetric_AdjFlow::~CSourceAxisymmetric_AdjFlow(void) { }

void CSourceAxisymmetric_AdjFlow::SetResidual(double *val_residual, double **Jacobian_ii, CConfig *config) {
    
    double yinv;
    double Jacobian_Axisymmetric[4][4];
    
	if (Coord_i[1] > 0.0) yinv = 1.0/Coord_i[1];
	else yinv = 0.0;
    
    Jacobian_Axisymmetric[0][0] = 0;
    Jacobian_Axisymmetric[0][1] = 0;
    Jacobian_Axisymmetric[0][2] = 1.;
    Jacobian_Axisymmetric[0][3] = 0;
    
    Jacobian_Axisymmetric[1][0] = -U_i[1]*U_i[2]/(U_i[0]*U_i[0]);
    Jacobian_Axisymmetric[1][1] = U_i[2]/U_i[0];
    Jacobian_Axisymmetric[1][2] = U_i[1]/U_i[0];
    Jacobian_Axisymmetric[1][3] = 0;
    
    Jacobian_Axisymmetric[2][0] = -U_i[2]*U_i[2]/(U_i[0]*U_i[0]);
    Jacobian_Axisymmetric[2][1] = 0;
    Jacobian_Axisymmetric[2][2] = 2*U_i[2]/U_i[0];
    Jacobian_Axisymmetric[2][3] = 0;
    
    Jacobian_Axisymmetric[3][0] = -Gamma*U_i[2]*U_i[3]/(U_i[0]*U_i[0]) + (Gamma-1)*U_i[2]*(U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]*U_i[0]);
    Jacobian_Axisymmetric[3][1] = -(Gamma-1)*U_i[2]*U_i[1]/(U_i[0]*U_i[0]);
    Jacobian_Axisymmetric[3][2] = Gamma*U_i[3]/U_i[0] - 1/2*(Gamma-1)*( (U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]) + 2*U_i[2]*U_i[2]/(U_i[0]*U_i[0]) );
    Jacobian_Axisymmetric[3][3] = Gamma*U_i[2]/U_i[0];
    
    for (int iVar=0; iVar<4; iVar++)
        for (int jVar=0; jVar<4; jVar++)
            Jacobian_Axisymmetric[iVar][jVar] *= yinv*Volume;
    
    /* -- Residual = transpose(Jacobian) * psi --*/
    for (int iVar = 0; iVar < nVar; iVar++) {
    	val_residual[iVar] = 0.0;
        for (int jVar = 0; jVar < nVar; jVar++) {
            val_residual[iVar] += Jacobian_Axisymmetric[jVar][iVar]*Psi_i[jVar];
            Jacobian_ii[iVar][jVar] = Jacobian_Axisymmetric[jVar][iVar];
        }
    }
}

CSourcePieceWise_Plasma::CSourcePieceWise_Plasma(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nSpecies, unsigned short val_nDiatomics,
                                                 unsigned short val_nMonatomics, CConfig *config) : CNumerics(val_nDim, val_nVar, val_nSpecies, val_nDiatomics, val_nMonatomics, config) {
    
	unsigned short iReaction;
	implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	U_id = new double [nVar];
	Temp_tr_id = new double[nSpecies];
    
	nSpecies 		= config->GetnSpecies();
	nReactions 		= config->GetnReactions();
	//	GammaMonatomic 	= config->GetGammaMonatomic();
	//	GammaDiatomic 	= config->GetGammaDiatomic();
	nMonatomics 	= config->GetnMonatomics();
	nDiatomics 		= config->GetnDiatomics();
    
	Molar_Mass 		= new double[nSpecies];
	ChargeNumber 	= new double[nSpecies];
	w_dot 			= new double[nSpecies];					w_dotd = new double[nSpecies];
	Q_tv 			= new double[nSpecies];
	Q_elastic 		= new double [nSpecies];			Q_elasticd	  = new double [nSpecies];
	Enthalpy_formation = new double [nSpecies];
    
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Molar_Mass[iSpecies] = config->GetMolar_Mass(iSpecies);
		ChargeNumber[iSpecies] = config->GetParticle_ChargeNumber(iSpecies);
		Enthalpy_formation[iSpecies] = config->GetEnthalpy_Formation(iSpecies);
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
        config->GetKind_GasModel() == AIR5 || config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == ARGON_SID) {
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
	Energy_vib = 0.0; Energy_el = 0.0;
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
    
	SourceVector = new double[nVar];
	SourceJacobian = new double*[nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar ++) {
		SourceVector[iVar] = 10.0;
		SourceJacobian[iVar] = new double [nVar];
		for (unsigned short jVar = 0; jVar < nVar; jVar ++)
			SourceJacobian[iVar][jVar] = 0.0;
	}
    
    
	MagneticField = new double [3];
	Mag_Force 	  = new double* [nSpecies];
	VcrossB = new double* [nSpecies];
	velocity = new double *[nSpecies];
    
	/* For the more accurate magnetic field model */
	VioncrossB = new double [nDim];
	Current_Density = new double [3];
	JcrossB = new double [3];
	Electric_Conductivity = config->GetElec_Conductivity();
	vector_r = new double [nDim];
	dpcenter = new double [nDim];
    
	for (iVar = 0; iVar < 3; iVar ++) {
		Current_Density[iVar] = 0.0;
		JcrossB[iVar] = 0.0;
	}
    
	/* Till here For the more accurate magnetic field model */
    
	MagneticDipole = new double [3];
	ElectricField = new double [nDim];
    
	for (iVar = 0; iVar < 3; iVar ++) {
		MagneticField[iVar] = 0.0;
		MagneticDipole[iVar] = 0.0;
	}
	Stagnation_B = config->GetStagnation_B();
	dpcenter[0] = config->GetDipoleDist();
	dpcenter[1] = 0.0;
	dpcenter[2] = 0.0;
    
	MagneticDipole[0] = Stagnation_B*2*PI_NUMBER*pow((0.01825+dpcenter[0]),3)/MAGNETIC_CONSTANT;
	MagneticDipole[1] = 0.0;
	MagneticDipole[2] = 0.0;
    
	for (iDim = 0; iDim < nDim; iDim++)
		ElectricField[iDim] = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		velocity[iSpecies] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim ++) velocity[iSpecies][iDim]  = 0.0;
		VcrossB[iSpecies] = new double [3];
		for (iVar = 0; iVar < 3; iVar ++) VcrossB[iSpecies][iVar] = 0.0;
		Mag_Force[iSpecies] = new double [3];
		for (iVar = 0; iVar < 3; iVar ++) Mag_Force[iSpecies][iVar] = 0.0;
        
	}
    
	AvgNum = AVOGAD_CONSTANT;
	M1Avg = config->GetMolar_Mass(0);
	M2Avg = config->GetMolar_Mass(1);
	M3Avg = config->GetMolar_Mass(2);
	M1M2M3Avg3 = M1Avg * M2Avg* M3Avg;
    
	M1 = M1Avg/AVOGAD_CONSTANT;
	M3 = M3Avg/AVOGAD_CONSTANT;
	M2 = M1-M3;
    
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
	delete [] Enthalpy_formation;
    
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
    
	//delete [] ElectricField;
	delete[] MagneticField;
    
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		delete [] velocity[iSpecies];
		delete [] VcrossB[iSpecies];
	}
	delete [] velocity;		delete [] VcrossB;
	delete [] VioncrossB;	delete [] Current_Density; delete [] JcrossB;
	delete [] MagneticDipole;
	delete [] dpcenter;
	delete [] vector_r;
    
}

void CSourcePieceWise_Plasma::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
    
	tol = 1E-60;
	double zero = 1E-15;
	double dens_ratio;
    
    //	cout << " Electric Field = " << ElectricField[0] << " , " << ElectricField[1] << endl;
    //    cin.get();
    
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
    
	P1 = (gam1-1)*(e1 - 0.5*(m1*m1 + n1*n1 + l1*l1)/r1 - (Enthalpy_formation[0] + Energy_vib + Energy_el)*r1 );	// Partial Pressure of species 1
	P2 = (gam2-1)*(e2 - 0.5*(m2*m2 + n2*n2 + l2*l2)/r2 - (Enthalpy_formation[1] + Energy_vib + Energy_el)*r2 );	// Partial Pressure of species 1
	P3 = (gam3-1)*(e3 - 0.5*(m3*m3 + n3*n3 + l3*l3)/r3 - (Enthalpy_formation[2] + Energy_vib + Energy_el)*r3 );	// Partial Pressure of species 1
    
	//	P1 = (gam1-1)*(e1 - 0.5*(m1*m1 + n1*n1 + l1*l1)/r1);	// Partial Pressure of species 2
	//	P2 = (gam2-1)*(e2 - 0.5*(m2*m2 + n2*n2 + l2*l2)/r2);	// Partial Pressure of species 2
	//	P3 = (gam3-1)*(e3 - 0.5*(m3*m3 + n3*n3 + l3*l3)/r3);	// Partial Pressure of species 3
    
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
    
	kb1 = kf1/ke1; 	kb2 = kf2/ke2; 	kb3 = kf3/ke3;
    
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
    
	mdotr = 0.0;
	distance = 0.0;
    
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		vector_r[iDim] = Coord_i[iDim] - dpcenter[iDim];
		mdotr +=MagneticDipole[iDim]*vector_r[iDim];
		distance += vector_r[iDim]*vector_r[iDim];
	}
	distance = sqrt(distance);
	rpower5 = pow(distance,5);
	rcubed  = pow(distance,3);
	for (iDim = 0; iDim < nDim; iDim ++ )
		MagneticField[iDim] = MAGNETIC_CONSTANT/(4.0*PI_NUMBER) * (3.0 * vector_r[iDim]* mdotr/rpower5 - MagneticDipole[iDim]/rcubed);
    
    
	//
	//
	//	mdotr = 0.0;
	//	distance = 0.0;
	//	for (iDim = 0; iDim < nDim; iDim ++ ) {
	//		mdotr +=MagneticDipole[iDim]*Coord_i[iDim];
	//		distance += Coord_i[iDim]*Coord_i[iDim];
	//	}
	//
	//	distance = sqrt(distance);
	//	rpower5 = pow(distance,5);
	//	rcubed  = pow(distance,3);
	//
	//	for (iDim = 0; iDim < nDim; iDim ++ )
	//		MagneticField[iDim] = MAGNETIC_CONSTANT/(4.0*PI_NUMBER) * (3.0 * Coord_i[iDim]* mdotr/rpower5 - MagneticDipole[iDim]/rcubed);
    
    
    
    
    
	double Bx = MagneticField[0];
	double By = MagneticField[1];
	double Bz = MagneticField[2];
	double s0 = Electric_Conductivity;
    
    
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		if (nDim ==2) {
			VcrossB[iSpecies][0] = 0.0;	VcrossB[iSpecies][1] = 0.0;
			Mag_Force[iSpecies][0] = 0.0; Mag_Force[iSpecies][1] = 0.0;
		}
		else {
			VcrossB[iSpecies][0] = velocity[iSpecies][1]*MagneticField[2] - velocity[iSpecies][2]*MagneticField[1];
			VcrossB[iSpecies][1] = velocity[iSpecies][2]*MagneticField[0] - velocity[iSpecies][0]*MagneticField[2];
			VcrossB[iSpecies][2] = velocity[iSpecies][0]*MagneticField[1] - velocity[iSpecies][1]*MagneticField[0];
			Mag_Force[iSpecies][0] = 0.0; Mag_Force[iSpecies][1] = 0.0;	Mag_Force[iSpecies][2] = 0.0;
            
		}
	}
    
	for (iDim = 0; iDim < nDim; iDim ++ )
		Current_Density[iDim] = Electric_Conductivity*(ElectricField[iDim] + VcrossB[0][iDim]);
    
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
    
	double delJcrossBdelm = s0/r1/r1 * ( Bx * (l1*Bz +   n1*By) + By * (n1*Bx - 2*m1*By) + Bz * (l1*Bx - 2*m1*Bz) );
	double delJcrossBdeln = s0/r1/r1 * ( Bx * (m1*By - 2*n1*Bx) + By * (m1*Bx +   l1*Bz) + Bz * (l1*By - 2*n1*Bz) );
	double delJcrossBdell = s0/r1/r1 * ( Bx * (m1*Bz - 2*l1*Bx) + By * (n1*Bz - 2*l1*By) + Bz * (m1*Bx +   n1*By) );
    
	double delJdotJdelm   =  1.0/s0*( Current_Density[0] * 0  	  - Bz*s0/r1*Current_Density[1] + By*s0/r1*Current_Density[2] );
	double delJdotJdeln   =  1.0/s0*( Current_Density[0] * Bz*s0/r1 - 0*Current_Density[0] 	    - Bx*s0/r1*Current_Density[2] );
	double delJdotJdell   =  1.0/s0*(-Current_Density[0] * By*s0/r1 + Bx*s0/r1*Current_Density[1] + 0*Current_Density[2] );
    
	/* Derivative of the first term in the source terms wrt all the conservative variables */
	unsigned short loc;
    
	//	JcrossB[0] = 0.0;	JcrossB[1] = 0.0;	JcrossB[2] = 0.0;	JdotJoversgima = 0.0;
	//	delJcrossBdelm = 0.0; delJcrossBdeln = 0.0; delJcrossBdell = 0.0;
	//	delJdotJdelm = 0.0;   delJdotJdeln = 0.0;   delJdotJdell = 0.0;
    
    
	loc = 0;
	dens_ratio = r1/(r1+r2+r3);
	JcrossB[0] = dens_ratio*JcrossB[0]; JcrossB[1] = dens_ratio*JcrossB[1];	JcrossB[2] = dens_ratio*JcrossB[2];
	JdotJoversgima = dens_ratio*JdotJoversgima;
	delJcrossBdelm = dens_ratio*delJcrossBdelm; delJcrossBdeln = dens_ratio*delJcrossBdeln;	delJcrossBdell = dens_ratio*delJcrossBdell;
	delJdotJdelm   = dens_ratio*delJdotJdelm;	delJdotJdeln   = dens_ratio*delJdotJdeln;	delJdotJdell   = dens_ratio*delJdotJdell;
    
    
    
	//	for (iDim = 0; iDim < nDim; iDim ++)
	//		Mag_Force[0][iDim] = dens_ratio*JcrossB[iDim];
    
	SourceVector[loc+0] = M1*AvgNum*R;
    
	SourceJacobian[loc+0][loc+0] =  M1Avg*dR_dr1;
	SourceJacobian[loc+0][loc+1] = M1Avg*dR_dm1;
	SourceJacobian[loc+0][loc+2] = M1Avg*dR_dn1;
	if (nDim ==3) SourceJacobian[loc+0][loc+3] = M1Avg*dR_dl1;
    
	SourceJacobian[loc+0][loc+nDim+1] = M1Avg*dR_de1;
    
	/* Derivative of the second term in the source terms wrt all the conservative variables */
    
	SourceVector[loc+1] = JcrossB[0] + r1*nu12*(m2/r2 - m1/r1) + r1*nu13*(m3/r3 - m1/r1);
    
	SourceJacobian[loc+1][loc+0] = nu12*m2/r2 + nu13*m3/r3 - JcrossB[0]/r1;
	SourceJacobian[loc+1][loc+1] = -nu12 - nu13  - s0 * ( Bz*Bz + By*By)/r1;
	SourceJacobian[loc+1][loc+2] = s0*By*Bx/r1;
	if (nDim ==3) SourceJacobian[loc+1][loc+3] =  s0*Bx*Bz/r1;
	SourceJacobian[loc+1][loc+nDim+1] = JcrossB[0] + r1*dv12_dT1*dT1_de1*(m2/r2 - m1/r1) + r1*dv13_dT1*dT1_de1*(m3/r3 - m1/r1);
    
	/* Derivative of the third term in the source terms wrt all the conservative variables */
    
	SourceVector[loc+2] = JcrossB[1] + r1*nu12*(n2/r2 - n1/r1) + r1*nu13*(n3/r3 - n1/r1);
    
	SourceJacobian[loc+2][loc+0]  = nu12*n2/r2 + nu13*n3/r3 - JcrossB[1]/r1;
	SourceJacobian[loc+2][loc+1]  = s0*Bx*By/r1;
	SourceJacobian[loc+2][loc+2]  = -nu12 - nu13 -s0*(Bx*Bx + Bz*Bz)/r1;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] =  s0*Bz*By/r1;
	SourceJacobian[loc+2][loc+nDim+1] = JcrossB[1] + r1*(n2/r2 -n1/r1)*dv12_dT1*dT1_de1 + r1*(n3/r3 - n1/r1)*dv13_dT1*dT1_de1;
    
    
	if (nDim == 3) {
		SourceVector[loc+3] = JcrossB[2] + r1*nu12*(l2/r2 - l1/r1) + r1*nu13*(l3/r3 - l1/r1);
        
		SourceJacobian[loc+3][loc+0]  = nu12*l2/r2 + nu13*l3/r3 - JcrossB[2]/r1;
		SourceJacobian[loc+3][loc+1]  = s0*Bx*Bz/r1;
		SourceJacobian[loc+3][loc+2]  = s0*By*Bz/r1;
		SourceJacobian[loc+3][loc+3] = -nu12 - nu13 - s0*(By*By + Bx*Bx)/r1;
		SourceJacobian[loc+3][loc+nDim+1] = JcrossB[2] + r1*(l2/r2 -l1/r1)*dv12_dT1*dT1_de1 + r1*(l3/r3 - l1/r1)*dv13_dT1*dT1_de1;
	}
	/* Derivative of the Fourth term in the source terms wrt all the conservative variables */
    
	SourceVector[loc+nDim+1] = JcrossB[0]*m1/r1 + JcrossB[1]*n1/r1 + JcrossB[2]*l1/r1 + JdotJoversgima +
    m1*nu12*(m2/r2 - m1/r1) + m1*nu13*(m3/r3 - m1/r1) +
    n1*nu12*(n2/r2 - n1/r1) + n1*nu13*(n3/r3 - n1/r1) +
    l1*nu12*(l2/r2 - l1/r1) + l1*nu13*(l3/r3 - l1/r1) + QT1;
    
	SourceJacobian[loc+nDim+1][loc+0]  =  -2/r1*(JcrossB[0]*m1/r1 + JcrossB[1]*n1/r1 + JcrossB[2]*l1/r1) - 2*JdotJoversgima/r1 + (m1*m1 + n1*n1 + l1*l1)/(r1*r1) * (nu12 + nu13);
	SourceJacobian[loc+nDim+1][loc+1]  = nu12*(m2/r2 -2*m1/r1) + nu13*(m3/r3 - 2*m1/r1) + delJcrossBdelm + delJdotJdelm + dQT1_dT1*dT1_dm1;
	SourceJacobian[loc+nDim+1][loc+2]  = nu12*(n2/r2 -2*n1/r1) + nu13*(n3/r3 - 2*n1/r1) + delJcrossBdeln + delJdotJdeln + dQT1_dT1*dT1_dn1;
	if (nDim == 3) SourceJacobian[loc+nDim+1][loc+3]  = nu12*(l2/r2 -2*l1/r1) + nu13*(l3/r3 - 2*l1/r1)   + delJcrossBdell + delJdotJdell + dQT1_dT1*dT1_dl1;
    
	SourceJacobian[loc+nDim+1][loc+nDim+1]  = JdotJoversgima/e1 +
    (m1*(m2/r2 - m1/r1) + n1*(n2/r2 - n1/r1) + l1*(l2/r2 - l1/r1))*dv12_dT1*dT1_de1 +
    (m1*(m3/r3 - m1/r1) + n1*(n3/r3 - n1/r1) + l1*(l3/r3 - l1/r1))*dv13_dT1*dT1_de1 + dQT1_dT1*dT1_de1;
    
    
    
    
    
    
    
	/* Derivative of the fifth term in the source terms wrt all the conservative variables */
	loc = (nDim+2);
    
	for (iDim = 0; iDim < nDim; iDim ++ )
		Current_Density[iDim] = Electric_Conductivity*(ElectricField[iDim] + VcrossB[1][iDim]);
    
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		JcrossB[0] = Current_Density[1]*MagneticField[2] - Current_Density[2]*MagneticField[1];
		JcrossB[1] = Current_Density[2]*MagneticField[0] - Current_Density[0]*MagneticField[2];
		JcrossB[2] = Current_Density[0]*MagneticField[1] - Current_Density[1]*MagneticField[0];
	}
	JdotJoversgima = 0.0;
	for (iDim = 0; iDim < nDim; iDim ++ )
		JdotJoversgima += (Current_Density[iDim]*Current_Density[iDim]);
    
	JdotJoversgima = JdotJoversgima/s0;
	delJcrossBdelm = s0/r2/r2 * ( Bx * (l2*Bz +   n2*By) + By * (n2*Bx - 2*m2*By) + Bz * (l2*Bx - 2*m2*Bz) );
	delJcrossBdeln = s0/r2/r2 * ( Bx * (m2*By - 2*n2*Bx) + By * (m2*Bx +   l2*Bz) + Bz * (l2*By - 2*n2*Bz) );
	delJcrossBdell = s0/r2/r2 * ( Bx * (m2*Bz - 2*l2*Bx) + By * (n2*Bz - 2*l2*By) + Bz * (m2*Bx +   n2*By) );
    
	delJdotJdelm   =  1/s0*( Current_Density[0] * 0  	  - Bz*s0/r2*Current_Density[1] + By*s0/r2*Current_Density[2] );
	delJdotJdeln   =  1/s0*( Current_Density[0] * Bz*s0/r2 - 0*Current_Density[0] 	    - Bx*s0/r2*Current_Density[2] );
	delJdotJdell   =  1/s0*(-Current_Density[0] * By*s0/r2 + Bx*s0/r2*Current_Density[1] + 0*Current_Density[2] );
    
	//	JcrossB[0] = 0.0;	JcrossB[1] = 0.0;	JcrossB[2] = 0.0;	JdotJoversgima = 0.0;
	//	delJcrossBdelm = 0.0; delJcrossBdeln = 0.0; delJcrossBdell = 0.0;
	//	delJdotJdelm = 0.0;   delJdotJdeln = 0.0;   delJdotJdell = 0.0;
    
	dens_ratio = r2/(r1+r2+r3);
	JcrossB[0] = dens_ratio*JcrossB[0]; JcrossB[1] = dens_ratio*JcrossB[1];	JcrossB[2] = dens_ratio*JcrossB[2];
	JdotJoversgima = dens_ratio*JdotJoversgima;
	delJcrossBdelm = dens_ratio*delJcrossBdelm; delJcrossBdeln = dens_ratio*delJcrossBdeln;	delJcrossBdell = dens_ratio*delJcrossBdell;
	delJdotJdelm   = dens_ratio*delJdotJdelm;	delJdotJdeln   = dens_ratio*delJdotJdeln;	delJdotJdell   = dens_ratio*delJdotJdell;
    
    
    
	//	for (iDim = 0; iDim < nDim; iDim ++)
	//		Mag_Force[1][iDim] = dens_ratio*JcrossB[iDim];
    
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
	SourceJacobian[loc+1][loc+nDim+1]  = JcrossB[0] + r2*(m1/r1 - m2/r2)*dv21_dT2*dT2_de2 + r2*(m3/r3 - m2/r2)*dv23_dT2*dT2_de2;
    
    
	/* Derivative of the Seventh term in the source terms wrt all the conservative variables */
    
	SourceVector[loc+2]  =  JcrossB[1]  + r2*nu21*(n1/r1 - n2/r2) + r2*nu23*(n3/r3 - n2/r2);
    
	SourceJacobian[loc+2][loc+0] = nu21*n1/r1 + nu23*n3/r3 - JcrossB[1]/r2;
	SourceJacobian[loc+2][loc+1] = s0*Bx*By/r2;
	SourceJacobian[loc+2][loc+2] = -nu21 - nu23 -s0*(Bx*Bx + Bz*Bz)/r2;
	if (nDim ==3) SourceJacobian[loc+2][loc+3] = s0*Bz*By/r2;
	SourceJacobian[loc+2][loc+nDim+1] = JcrossB[1] + r2*(n1/r1 - n2/r2)*dv21_dT2*dT2_de2 + r2*(n3/r3 - n2/r2)*dv23_dT2*dT2_de2;
    
	if (nDim == 3) {
		SourceVector[loc+3] = JcrossB[2]  + r2*nu21*(l1/r1 - l2/r2) + r2*nu23*(l3/r3 - l2/r2);
        
		SourceJacobian[loc+3][loc+0]  = nu21*l1/r1 + nu23*l3/r3 - JcrossB[2]/r2;
		SourceJacobian[loc+3][loc+1]  = s0*Bx*Bz/r2;
		SourceJacobian[loc+3][loc+2]  = s0*By*Bz/r2;
		SourceJacobian[loc+3][loc+3] = -nu21 - nu23 - s0*(By*By + Bx*Bx)/r2;
		SourceJacobian[loc+3][loc+nDim+1] = JcrossB[2]+ r2*(l1/r1 -l2/r2)*dv21_dT2*dT2_de2 + r2*(l3/r3 - l2/r2)*dv23_dT2*dT2_de2;
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
    
	for (iDim = 0; iDim < nDim; iDim ++)
		Mag_Force[2][iDim] = JcrossB[iDim];
    
	//	JcrossB[0] = 0.0;	JcrossB[1] = 0.0;	JcrossB[2] = 0.0;	JdotJoversgima = 0.0;
	//	delJcrossBdelm = 0.0; delJcrossBdeln = 0.0; delJcrossBdell = 0.0;
	//	delJdotJdelm = 0.0;   delJdotJdeln = 0.0;   delJdotJdell = 0.0;
    
	dens_ratio = r3/(r1+r2+r3);
	JcrossB[0] = dens_ratio*JcrossB[0]; JcrossB[1] = dens_ratio*JcrossB[1];	JcrossB[2] = dens_ratio*JcrossB[2];
	JdotJoversgima = dens_ratio*JdotJoversgima;
	delJcrossBdelm = dens_ratio*delJcrossBdelm; delJcrossBdeln = dens_ratio*delJcrossBdeln;	delJcrossBdell = dens_ratio*delJcrossBdell;
	delJdotJdelm   = dens_ratio*delJdotJdelm;	delJdotJdeln   = dens_ratio*delJdotJdeln;	delJdotJdell   = dens_ratio*delJdotJdell;
    
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
	SourceJacobian[loc+1][loc+nDim+1] =  JcrossB[0] + r3*(m1/r1 - m3/r3)*dv31_dT3*dT3_de3 + r3*(m2/r2 - m3/r3)*dv32_dT3*dT3_de3;
    
	/* Derivative of the Eleventh term in the source terms wrt all the conservative variables */
	SourceVector[loc+2]  = JcrossB[1]  + r3*nu31*(n1/r1 - n3/r3) + r3*nu32*(n2/r2 - n3/r3);
    
	SourceJacobian[loc+2][loc+0] = nu31*n1/r1 + nu32*n2/r2 - JcrossB[1]/r3;
	SourceJacobian[loc+2][loc+1] = s0*Bx*By/r3;
	SourceJacobian[loc+2][loc+2] =  -nu31 - nu32 -s0*(Bx*Bx + Bz*Bz)/r3;
	if (nDim ==3)  SourceJacobian[loc+2][loc+3] =  s0*Bz*By/r3 ;
	SourceJacobian[loc+2][loc+nDim+1] = JcrossB[1] + r3*(n1/r1 - n3/r3)*dv31_dT3*dT3_de3 + r3*(n2/r2 - n3/r3)*dv32_dT3*dT3_de3;
    
	/* Derivative of the Twelfth term in the source terms wrt all the conservative variables */
	if (nDim ==3) {
		SourceVector[loc+3]  = JcrossB[2] + r3*nu31*(l1/r1 - l3/r3) + r3*nu32*(l2/r2 - l3/r3);
		SourceJacobian[loc+3][loc+0] = nu31*l1/r1 + nu32*l2/r2 - JcrossB[2]/r3;
		SourceJacobian[loc+3][loc+1] =  s0*Bx*Bz/r3;
		SourceJacobian[loc+3][loc+2] =  s0*By*Bz/r3;
		SourceJacobian[loc+3][loc+3] =  -nu31 - nu32 - s0*(By*By + Bx*Bx)/r3 ;
		SourceJacobian[loc+3][loc+nDim+1] =  JcrossB[2] + r3*(l1/r1 - l3/r3)*dv31_dT3*dT3_de3 + r3*(l2/r2 - l3/r3)*dv32_dT3*dT3_de3;
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
    
    
    
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		JcrossB[0] = Current_Density[1]*MagneticField[2] - Current_Density[2]*MagneticField[1];
		JcrossB[1] = Current_Density[2]*MagneticField[0] - Current_Density[0]*MagneticField[2];
		JcrossB[2] = Current_Density[0]*MagneticField[1] - Current_Density[1]*MagneticField[0];
	}
    
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
    
    
	double yinv, Enthalpy_i, Velocity_i, sq_vel;
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
        
		/*--- Set the primitive variables (T_tr,    vel,  T_vib,      p,    rho,      h,      c) and compute the mean ---*/
		/*--- Indices:                    (   0, iDim+1, nDim+1, nDim+2, nDim+3, nDim+4, nDim+5)                      ---*/
		Enthalpy_i = Varray_i[iSpecies][nDim+4];
        
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
			}
		}
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
			}
		}
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
	double T_tr;
    
	T_min = 800;
	epsilon = 80;
    
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
				T_tr = Varray_i[iSpecies][0];
				T_rxnf[iReaction] *= T_tr;
				counterFwd++;
			}
			/*--- Products ---*/
			if (jSpecies != nSpecies) {
				T_tr = Varray_i[jSpecies][0];
				T_rxnb[iReaction] *= T_tr;
				counterBkw++;
			}
		}
		T_rxnf[iReaction] = exp(1.0/counterFwd*log(T_rxnf[iReaction]));
		T_rxnb[iReaction] = exp(1.0/counterBkw*log(T_rxnb[iReaction]));
        
		/*--- Apply a modified temperature to ease the stiffness at low temperatures ---*/
		T_rxnf[iReaction] = 0.5 * (T_rxnf[iReaction]+T_min + sqrt((T_rxnf[iReaction]-T_min)*(T_rxnf[iReaction]-T_min) + epsilon*epsilon));
		T_rxnb[iReaction] = 0.5 * (T_rxnb[iReaction]+T_min + sqrt((T_rxnb[iReaction]-T_min)*(T_rxnb[iReaction]-T_min) + epsilon*epsilon));
        
		GetEq_Rxn_Coefficients(EqRxnConstants, config);
        
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
	//SU2_CPP2C END CSourcePieceWise_Plasma::SetResidual_Chemistry
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
	double T_tr_i, T_tr_j;
    
	/*--- Initialization ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		for (iDim = 0; iDim < nDim; iDim++)
			P[iSpecies][iDim] = 0.0;
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = 0.0;
    
	electron_temperature = Varray_i[nSpecies-1][0];
    
	/*--- Solve for momentum exchange between species due to collisions ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
		else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
			if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
			else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
            
			if (iSpecies != jSpecies) {
				T_tr_i = Varray_i[iSpecies][0];
				T_tr_j = Varray_i[jSpecies][0];
				T_control = sqrt(T_tr_i*T_tr_j);
				collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] * pow(T_control, Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2]);
                
				/*----- An if-condition for the case when a positive charge collides with an electron -----*/
				if( (ChargeNumber[iSpecies] == 1 && jSpecies == nSpecies -1) || (iSpecies == nSpecies-1 && ChargeNumber[jSpecies] == 1)) {
					radius_electronIonCollision = ELECTRON_CHARGE*ELECTRON_CHARGE/(32.0*FREE_PERMITTIVITY*BOLTZMANN_CONSTANT*electron_temperature);
					collisionArea = PI_NUMBER * radius_electronIonCollision * radius_electronIonCollision;
				}
                
				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*T_tr_i/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*T_tr_j/(PI_NUMBER*Molecular_Mass[jSpecies]));
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
	}
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
	double T_tr_i, T_tr_j;
    
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
        
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){
			if (jSpecies != iSpecies) {
				if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
				else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
                
				T_tr_i = Varray_i[iSpecies][0];
				T_tr_j = Varray_i[jSpecies][0];
				T_control = sqrt(T_tr_i*T_tr_j);
				collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] * pow(T_control, Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2]);
                
				velocity_Species_i = sqrt(8.0*BOLTZMANN_CONSTANT*T_tr_i/(PI_NUMBER*Molecular_Mass[iSpecies]));
				velocity_Species_j = sqrt(8.0*BOLTZMANN_CONSTANT*T_tr_j/(PI_NUMBER*Molecular_Mass[jSpecies]));
				collisionVelocity = sqrt(velocity_Species_i*velocity_Species_i + velocity_Species_j*velocity_Species_j);
				collisionFreq = (U_i[jLoc+0] * collisionArea * collisionVelocity)/(Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);
				vel_dot_prod = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					vel_dot_prod += (U_i[iLoc+1+iDim] - U_i[jLoc+1+iDim])*U_i[iLoc+1+iDim];
				}
                
				/* Exchange from Lee and Sutton, heavy particles */
				coefficient = 2.0*U_i[iLoc+0] / (Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]) * collisionFreq * 3.0/2.0*BOLTZMANN_CONSTANT;
				Q_elastic[iSpecies] += coefficient*(T_tr_j - T_tr_i);
			}
		}
		val_residual[iLoc+1+nDim] =  Q_elastic[iSpecies] * Volume;
	}
    
	/*--- Translational-rotational & vibrational energy exchange via inelastic collisions ---*/
	//Comment: Landau-Teller formulation
	//Comment: Based on Scalabrin and Candler.  May need to re-visit with much more detail and derive for my formulation
    
	double tau_sr, tau_ps, LimitingXSection, AvgMolecularSpeed, ReducedMass, A_sr, estar_vs, e_vs, q_tr_vs;
	double MixturePressure, MixtureNumDensity;
	double P_i, P_j;
    
	/*--- Calculate mixture quantities ---*/
	MixturePressure   = 0.0;
	MixtureNumDensity = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) {
			iLoc = (nDim+3)*iSpecies;
			CharVibTemp[iSpecies] = config->GetCharVibTemp(iSpecies);
		} else
			iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		MixturePressure += Varray_i[iSpecies][nDim+2];//SpeciesPressure_i[iSpecies];
		MixtureNumDensity += U_i[iLoc+0]/Molecular_Mass[iSpecies];
	}
    
	for (iSpecies = 0; iSpecies < nDiatomics; iSpecies++) {
		iLoc = (nDim+3)*iSpecies;
		Q_tv[iSpecies] = 0.0;
		for (jSpecies = 0; jSpecies < nSpecies; jSpecies++){
			if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
			else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
            
			T_tr_i = Varray_i[iSpecies][0];
			T_tr_j = Varray_i[jSpecies][0];
			P_i = Varray_i[iSpecies][nDim+2];
			P_j = Varray_i[jSpecies][nDim+2];
            
			/*--- Calculate Landau-Teller relaxation time ---*/
			ReducedMass = Molar_Mass[iSpecies]*Molar_Mass[jSpecies]/(Molar_Mass[iSpecies]+Molar_Mass[jSpecies]);
			LimitingXSection = 1E-20*pow(50000/T_tr_i, 2.0);
            
			AvgMolecularSpeed = sqrt(8.0*UNIVERSAL_GAS_CONSTANT*T_tr_i/(PI_NUMBER*Molar_Mass[iSpecies]));
            
			A_sr = 1.16 * 1E-3 * sqrt(ReducedMass) * pow(CharVibTemp[iSpecies], 4.0/3.0);
			tau_sr =   101325/(P_i+P_j)* exp(A_sr*(pow(sqrt(T_tr_i*T_tr_j),-1.0/3.0) - 0.015*pow(ReducedMass,0.25)) - 18.42);
			tau_ps = 1.0/(LimitingXSection*AvgMolecularSpeed*MixtureNumDensity);
            
			estar_vs = UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies] * CharVibTemp[iSpecies]/(exp(CharVibTemp[iSpecies]/T_tr_j)-1.0);
			e_vs = U_i[iLoc+nDim+2]/U_i[iLoc+0];
            
			/*--- Energy transferred per-molecule from species r to species s ---*/
			q_tr_vs = Molecular_Mass[iSpecies] * (estar_vs - e_vs)/(tau_sr + tau_ps);
            
			/*--- Convert to energy per volume for r and s species and multiply by volume for residual ---*/
			val_residual[iLoc+nDim+2] += q_tr_vs*U_i[iLoc+0]/Molecular_Mass[iSpecies]*Volume;
			val_residual[iLoc+nDim+1] += q_tr_vs*U_i[iLoc+0]/Molecular_Mass[iSpecies]*Volume;
			val_residual[jLoc+nDim+1] -= q_tr_vs*U_i[jLoc+0]/Molecular_Mass[jSpecies]*Volume;
		}
	}
	//SU2_CPP2C END CSourcePieceWise_Plasma::SetResidual_EnergyExch
}


CSource_Magnet::CSource_Magnet(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
    
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    
	MagneticField = new double [3];
	MagneticDipole = new double [3];
	Current_Density = new double [3];
	JcrossB = new double [3];
	VcrossB = new double [3];
	Electric_Conductivity = config->GetElec_Conductivity();
	vector_r = new double [nDim];
	dpcenter = new double [nDim];
	velocity = new double[nDim];
    
	for (iVar = 0; iVar < 3; iVar ++) {
		MagneticField[iVar] = 0.0;
		MagneticDipole[iVar] = 0.0;
		Current_Density[iVar] = 0.0;
		JcrossB[iVar] = 0.0;
		VcrossB[iVar] = 0.0;
	}
    
	Stagnation_B = config->GetStagnation_B();
	dpcenter[0] = config->GetDipoleDist();
	dpcenter[1] = 0.0;
	dpcenter[2] = 0.0;
    
	MagneticDipole[0] = Stagnation_B*2*PI_NUMBER*pow((0.01825+dpcenter[0]),3)/MAGNETIC_CONSTANT;
	MagneticDipole[1] = 0.0;
	MagneticDipole[2] = 0.0;
    
	if (nDim ==2) VcrossB[iDim] = 0.0;
    
    
}

CSource_Magnet::~CSource_Magnet(void) {
	delete []	MagneticField;
	delete []	MagneticDipole;
	delete []	Current_Density;
	delete []	JcrossB;
	delete []   VcrossB;
	delete []   velocity;
}

void CSource_Magnet::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
    
	double mdotr, distance, rpower5, rcubed;
	double r1,m1,n1,l1,e1;
	r1 = U_i[0];	// density of each species
	m1 = U_i[1];	// density*horizontal velocity
	n1 = U_i[2];	// density*vertical velocity
	e1 = U_i[3];	// energy per unit volume
	l1 = 0.0;			// density*vertical velocity of each species
    
	if (nDim ==3 ) {
		l1 = U_i[3];	// density*vertical velocity
		e1 = U_i[4];
	}
    
	velocity[0] = m1/r1; 	velocity[1] = n1/r1;
	if (nDim ==3) velocity[2] = l1/r1;
    
    
	mdotr = 0.0;
	distance = 0.0;
    
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		vector_r[iDim] = Coord_i[iDim] - dpcenter[iDim];
		mdotr +=MagneticDipole[iDim]*vector_r[iDim];
		distance += vector_r[iDim]*vector_r[iDim];
	}
	distance = sqrt(distance);
	rpower5 = pow(distance,5);
	rcubed  = pow(distance,3);
	for (iDim = 0; iDim < nDim; iDim ++ )
		MagneticField[iDim] = MAGNETIC_CONSTANT/(4.0*PI_NUMBER) * (3.0 * vector_r[iDim]* mdotr/rpower5 - MagneticDipole[iDim]/rcubed);
    
	double Bx = MagneticField[0];
	double By = MagneticField[1];
	double Bz = MagneticField[2];
	double s0 = Electric_Conductivity;
    
	if (nDim ==3) {
		VcrossB[0] = velocity[1]*MagneticField[2] - velocity[2]*MagneticField[1];
		VcrossB[1] = velocity[2]*MagneticField[0] - velocity[0]*MagneticField[2];
		VcrossB[2] = velocity[0]*MagneticField[1] - velocity[1]*MagneticField[0];
	}
    
	for (iDim = 0; iDim < nDim; iDim ++ )
		Current_Density[iDim] = Electric_Conductivity*VcrossB[iDim];
    
	for (iDim = 0; iDim < nDim; iDim ++ ) {
		JcrossB[0] = Current_Density[1]*MagneticField[2] - Current_Density[2]*MagneticField[1];
		JcrossB[1] = Current_Density[2]*MagneticField[0] - Current_Density[0]*MagneticField[2];
		JcrossB[2] = Current_Density[0]*MagneticField[1] - Current_Density[1]*MagneticField[0];
	}
    
	double JdotJoversgima = 0.0;
    
	for (iDim = 0; iDim < nDim; iDim ++ )
		JdotJoversgima += (Current_Density[iDim]*Current_Density[iDim]);
    
	JdotJoversgima = JdotJoversgima/s0 ;
    
	double delJcrossBdelm = s0/r1/r1 * ( Bx * (l1*Bz +   n1*By) + By * (n1*Bx - 2*m1*By) + Bz * (l1*Bx - 2*m1*Bz) );
	double delJcrossBdeln = s0/r1/r1 * ( Bx * (m1*By - 2*n1*Bx) + By * (m1*Bx +   l1*Bz) + Bz * (l1*By - 2*n1*Bz) );
	double delJcrossBdell = s0/r1/r1 * ( Bx * (m1*Bz - 2*l1*Bx) + By * (n1*Bz - 2*l1*By) + Bz * (m1*Bx +   n1*By) );
    
	double delJdotJdelm   =  1.0/s0*( Current_Density[0] * 0  	  - Bz*s0/r1*Current_Density[1] + By*s0/r1*Current_Density[2] );
	double delJdotJdeln   =  1.0/s0*( Current_Density[0] * Bz*s0/r1 - 0*Current_Density[0] 	    - Bx*s0/r1*Current_Density[2] );
	double delJdotJdell   =  1.0/s0*(-Current_Density[0] * By*s0/r1 + Bx*s0/r1*Current_Density[1] + 0*Current_Density[2] );
    
	/* Derivative of the second term in the source terms wrt all the conservative variables */
    
	val_residual[1] = JcrossB[0]*Volume;
	val_residual[2] = JcrossB[1]*Volume;
	if (nDim == 3)	val_residual[3] = JcrossB[2]*Volume;
	val_residual[nDim+1] = (JcrossB[0]*m1/r1 + JcrossB[1]*n1/r1 + JcrossB[2]*l1/r1 + JdotJoversgima)*Volume;
    
	if (implicit) {
		val_Jacobian_i[1][0] = -JcrossB[0]/r1*Volume;
		val_Jacobian_i[1][1] =- s0 * ( Bz*Bz + By*By)/r1*Volume;
		val_Jacobian_i[1][2] = s0*By*Bx/r1*Volume;
		if (nDim ==3) val_Jacobian_i[1][3] =  s0*Bx*Bz/r1*Volume;
		val_Jacobian_i[1][nDim+1] = JcrossB[0]*Volume;
        
		val_Jacobian_i[2][0]  = - JcrossB[1]/r1*Volume;
		val_Jacobian_i[2][1]  = s0*Bx*By/r1*Volume;
		val_Jacobian_i[2][2]  = -s0*(Bx*Bx + Bz*Bz)/r1*Volume;
		if (nDim ==3) val_Jacobian_i[2][3] =  s0*Bz*By/r1*Volume;
		val_Jacobian_i[2][nDim+1] = JcrossB[1]*Volume;
        
		if (nDim == 3) {
			val_Jacobian_i[3][0]  = - JcrossB[2]/r1*Volume;
			val_Jacobian_i[3][1]  = s0*Bx*Bz/r1*Volume;
			val_Jacobian_i[3][2]  = s0*By*Bz/r1*Volume;
			val_Jacobian_i[3][3] = - s0*(By*By + Bx*Bx)/r1*Volume;
			val_Jacobian_i[3][nDim+1] = JcrossB[2]*Volume;
		}
        
		val_Jacobian_i[nDim+1][0]  =  (-2/r1*(JcrossB[0]*m1/r1 + JcrossB[1]*n1/r1 + JcrossB[2]*l1/r1) - 2*JdotJoversgima/r1)*Volume;
		val_Jacobian_i[nDim+1][1]  = (delJcrossBdelm + delJdotJdelm)*Volume;
		val_Jacobian_i[nDim+1][2]  = (delJcrossBdeln + delJdotJdeln)*Volume;
		if (nDim == 3) val_Jacobian_i[nDim+1][3]  = (delJcrossBdell + delJdotJdell)*Volume;
		val_Jacobian_i[nDim+1][nDim+1]  = JdotJoversgima/e1*Volume;
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

void CSource_JouleHeating::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
    
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

CSource_Template::CSource_Template(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {}

CSource_Template::~CSource_Template(void) {
    
}

void CSource_Template::SetResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {}



















