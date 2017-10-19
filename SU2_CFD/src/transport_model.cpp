/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/transport_model.hpp"
#include "../../../../FluidProp/api/c/fluidprop.h"


//#ifdef HAVE_FluidProp
//#include "fluidprop.h"
//#endif

/*-------------------------------------------------*/
/*----------- Dynamic Viscosity Models ------------*/
/*-------------------------------------------------*/

CViscosityModel::CViscosityModel(void) {

  /*--- Attributes initialization ---*/

  Mu = 0.0;
  dmudrho_T = 0.0;
  dmudT_rho = 0.0;

}

CViscosityModel::~CViscosityModel(void) { }


CConstantViscosity::CConstantViscosity(void) : CViscosityModel() { }

CConstantViscosity::CConstantViscosity(su2double mu_const) : CViscosityModel() {

  /*--- Attributes initialization ---*/

  Mu = mu_const;
  dmudrho_T = 0.0;
  dmudT_rho = 0.0;

}

CConstantViscosity::~CConstantViscosity(void) { }




CSutherland::CSutherland(void) : CViscosityModel() {
  Mu_ref = 0.0;
  T_ref = 0.0;
  S = 0.0;

}

CSutherland::CSutherland(su2double mu_ref, su2double t_ref, su2double s) : CViscosityModel() {

  Mu_ref = mu_ref;
  T_ref = t_ref;
  S = s;
}

CSutherland::~CSutherland(void) { }


void CSutherland::SetViscosity(su2double T, su2double rho) {

  Mu = Mu_ref*pow((T/T_ref),(3.0/2.0))*((T_ref + S)/(T + S));

}

void CSutherland::SetDerViscosity(su2double T, su2double rho) {

  dmudrho_T = 0.0;
  dmudT_rho = Mu_ref*( (3.0/2.0)*pow( (T/T_ref),(1.0/2.0) )*( (T_ref + S)/(T + S) )
          -pow( (T/T_ref),(3.0/2.0) )*(T_ref + S)/(T + S)/(T + S) );

}



CFluidPropViscosity::CFluidPropViscosity(void) : CViscosityModel() {

}

CFluidPropViscosity::~CFluidPropViscosity(void) { }


void CFluidPropViscosity::SetViscosity(double T, double rho) {

	   double T2 = T, rho2 = rho, mu;

        double deta_dT, deta_drho, lambda, dlambda_dT, dlambda_drho, sigma; 
        fluidprop_alltransprops( "Td", T2, rho2, &mu, &deta_dT, &deta_drho, &lambda, &dlambda_dT, &dlambda_drho, &sigma);

        Mu = mu;

        if (strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "T = %f, rho = %f, mu = %f\n", T, rho, Mu);
        }

}

void CFluidPropViscosity::SetDerViscosity(double T, double rho)  {

	double T2 = T, rho2 = rho, dmu1, dmu2;

        double eta, lambda, dlambda_dT, dlambda_drho, sigma; 
        fluidprop_alltransprops( "Td", T2, rho2, &eta, &dmu1, &dmu2, &lambda, &dlambda_dT, &dlambda_drho, &sigma);

        dmudT_rho = dmu1;
		dmudrho_T = dmu2;

        if (strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "T = %f, rho = %f, dmudT_rho = %f, dmudrho_dT = %f\n", T, rho, dmudT_rho, dmudT_rho);
        }

}



CCO2_Viscosity :: CCO2_Viscosity(void) : CViscosityModel() {
  Mu_ref = 0.0;
  T_ref = 0.0;
}

CCO2_Viscosity::CCO2_Viscosity(su2double mu_ref, su2double t_ref, su2double t_crit) : CViscosityModel() {

  Mu_ref = mu_ref;
  T_crit = t_crit;
  T_ref = t_ref;

  visc_coeff = new su2double [6];

  visc_coeff[0] = -0.0662;
  visc_coeff[1] = -5.7641;
  visc_coeff[2] = +71.441;
  visc_coeff[3] = -467.67;
  visc_coeff[4] = +1480.7;
  visc_coeff[5] = -1777.9;

}

CCO2_Viscosity::~CCO2_Viscosity(void) {
	delete [] visc_coeff;
}


void CCO2_Viscosity::SetViscosity(su2double T, su2double rho) {

  unsigned int i;

  Mu = 0;

  for (i=0; i<6 ; i++)
	Mu += visc_coeff[i] * pow(1- T*T_ref/T_crit, i);

  Mu = pow(10, Mu*T_crit/T/T_ref + log(30)/log(10));
  Mu = Mu * 1e-6 / Mu_ref;

}

void CCO2_Viscosity::SetDerViscosity(su2double T, su2double rho) {

  unsigned int i;
  su2double der_mu_partial;

  dmudrho_T = 0;

  der_mu_partial = 0;

  for (i=0; i<6 ; i++)
	  der_mu_partial += visc_coeff[i] * pow(1- T*T_ref/T_crit, i);
  der_mu_partial = -der_mu_partial * T_crit/pow(T*T_ref, 2);

  for (i=1; i<6 ; i++)
	  der_mu_partial += visc_coeff[i] * i * pow(1- T*T_ref/T_crit, i-1) * T_crit/T/T_ref;

  dmudT_rho = Mu*Mu_ref * log(10) * der_mu_partial;
  dmudT_rho = dmudT_rho/ Mu_ref * T_ref;


}


CR22_Viscosity :: CR22_Viscosity(void) : CViscosityModel() {
  Mu_ref = 0.0;
  T_ref = 0.0;

}

CR22_Viscosity::CR22_Viscosity(su2double mu_ref, su2double t_ref) : CViscosityModel() {

  Mu_ref = mu_ref;
  T_ref = t_ref;

}

CR22_Viscosity::~CR22_Viscosity(void) {}


void CR22_Viscosity::SetViscosity(su2double T, su2double rho) {

  unsigned int i;

  Mu = -52.658 + 11.509*log(T*T_ref);
  Mu = Mu*1e-6/Mu_ref;

}

void CR22_Viscosity::SetDerViscosity(su2double T, su2double rho) {

  dmudrho_T = 0.0;
  dmudT_rho = 11.509*1/T/T_ref*1e-6;
  dmudT_rho = dmudT_rho/Mu_ref*T_ref;

}


CR12_Viscosity :: CR12_Viscosity(void) : CViscosityModel() {
  Mu_ref = 0.0;
  T_ref = 0.0;

}

CR12_Viscosity::CR12_Viscosity(su2double mu_ref, su2double t_ref) : CViscosityModel() {

  Mu_ref = mu_ref;
  T_ref = t_ref;

}

CR12_Viscosity::~CR12_Viscosity(void) {}


void CR12_Viscosity::SetViscosity(su2double T, su2double rho) {

  unsigned int i;

  Mu = -35.324 + 8.4054*log(T*T_ref);
  Mu = Mu*1e-6/Mu_ref;

}

void CR12_Viscosity::SetDerViscosity(su2double T, su2double rho) {

  dmudrho_T = 0.0;
  dmudT_rho = 8.4054*1/T/T_ref*1e-6;
  dmudT_rho = dmudT_rho/Mu_ref*T_ref;

}


/* ------------------------------------------------- */
/* ---------- Thermal Conductivity Models ---------- */
/* ------------------------------------------------- */

CConductivityModel::CConductivityModel(void) {

  /*--- Attributes initialization ---*/

  Kt = 0.0;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;

}

CConductivityModel::~CConductivityModel(void) { }


CConstantConductivity::CConstantConductivity(void) : CConductivityModel() { }

CConstantConductivity::CConstantConductivity(su2double kt_const) : CConductivityModel() {

  /*--- Attributes initialization ---*/

  Kt = kt_const;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;

}

CConstantConductivity::~CConstantConductivity(void) { }


CConstantPrandtl::CConstantPrandtl(void) : CConductivityModel() { }

CConstantPrandtl::CConstantPrandtl(su2double pr_const) : CConductivityModel() {

  /*--- Attributes initialization ---*/

  Pr_const = pr_const;

}

void CConstantPrandtl::SetConductivity(su2double T, su2double rho, su2double mu, su2double cp) {

  Kt = mu*cp/Pr_const;

}

void CConstantPrandtl::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {

  dktdrho_T = dmudrho_T*cp/Pr_const;
  dktdT_rho = dmudT_rho*cp/Pr_const;

}

CConstantPrandtl::~CConstantPrandtl(void) { }




CFluidPropConductivity::CFluidPropConductivity(void) : CConductivityModel() { }

CFluidPropConductivity::~CFluidPropConductivity(void) { }

CFluidPropConductivity::CFluidPropConductivity(double pr_const) : CConductivityModel() {

  //--- Attributes initialization ---

	Pr_const = pr_const;

}

void CFluidPropConductivity::SetConductivity(double T, double rho, double mu, double cp) {

	    double T2 = T, rho2 = rho, kt;
        double eta, deta_dT, deta_drho, dlambda_dT, dlambda_drho, sigma; 
        fluidprop_alltransprops( "Td", T2, rho2, &eta, &deta_dT, &deta_drho, &kt, &dlambda_dT, &dlambda_drho, &sigma);

        Kt = kt;

        if (strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "T = %f, rho = %f, Kt = %f\n", T, rho, Kt);
        }

}

void CFluidPropConductivity::SetDerConductivity(double T, double rho, double dmudrho_T, double dmudT_rho, double cp) {


	    double T2 = T, rho2 = rho, dk1, dk2;
        double eta, deta_dT, deta_drho, lambda, sigma; 
        fluidprop_alltransprops( "Td", T2, rho2, &eta, &deta_dT, &deta_drho, &lambda, &dk1, &dk2, &sigma);

        dktdT_rho = dk1;
        dktdrho_T = dk2;


        if (strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "T = %f, rho = %f, dktdT_rho = %f, dktdrho_dT = %f\n", T, rho, dktdT_rho, dktdT_rho);
        }

}


CCO2_Conductivity :: CCO2_Conductivity(void) : CConductivityModel() {
  K_ref = 0.0;
  T_ref = 0.0;

}

CCO2_Conductivity::CCO2_Conductivity(su2double k_ref, su2double t_ref, su2double t_crit) : CConductivityModel() {

  K_ref = k_ref;
  T_crit = t_crit;
  T_ref = t_ref;

  cond_coeff = new su2double [6];

  cond_coeff[0] = -0.0055;
  cond_coeff[1] = -8.1303;
  cond_coeff[2] = +86.575;
  cond_coeff[3] = -536.02;
  cond_coeff[4] = +1662.3;
  cond_coeff[5] = -1981.4;


}

CCO2_Conductivity::~CCO2_Conductivity(void) {
	delete [] cond_coeff;
}


void CCO2_Conductivity::SetConductivity(su2double T, su2double rho, su2double mu, su2double cp) {

  unsigned int i;

  Kt = 0;

  for (i=0; i<6 ; i++)
	Kt += cond_coeff[i] * pow(1- T*T_ref/T_crit, i);

  Kt = pow(10, (Kt*T_crit/T/T_ref + log(50)/log(10)));

  Kt = Kt * 1e-3 / K_ref;

}

void CCO2_Conductivity::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {

  unsigned int i;
  su2double der_k_partial;

  dktdrho_T = 0;

  der_k_partial = 0;

  for (i=0; i<6 ; i++)
	  der_k_partial += cond_coeff[i] * pow(1- T*T_ref/T_crit, i);
  der_k_partial = -der_k_partial * T_crit/pow(T*T_ref, 2);

  for (i=1; i<6 ; i++)
	  der_k_partial += cond_coeff[i] * i * pow(1- T*T_ref/T_crit, i-1) * T_crit/T/T_ref;

  dktdT_rho = Kt*K_ref * log(10) * der_k_partial;
  dktdT_rho = dktdT_rho/ K_ref * T_ref;


}




CR22_Conductivity :: CR22_Conductivity(void) : CConductivityModel() {
  K_ref = 0.0;
  T_ref = 0.0;

}

CR22_Conductivity::CR22_Conductivity(su2double k_ref, su2double t_ref) : CConductivityModel() {

  K_ref = k_ref;
  T_ref = t_ref;

}

CR22_Conductivity::~CR22_Conductivity(void) {}


void CR22_Conductivity::SetConductivity(su2double T, su2double rho, su2double mu, su2double cp) {


  Kt = 0.048703*T*T_ref-3.0435;
  Kt = Kt*1e-3/K_ref;

}

void CR22_Conductivity::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {

  dktdrho_T = 0;
  dktdT_rho = 0.048703;
  dktdT_rho = dktdT_rho *1e-3/K_ref * T_ref;
}

CR12_Conductivity :: CR12_Conductivity(void) : CConductivityModel() {
  K_ref = 0.0;
  T_ref = 0.0;

}

CR12_Conductivity::CR12_Conductivity(su2double k_ref, su2double t_ref) : CConductivityModel() {

  K_ref = k_ref;
  T_ref = t_ref;

}

CR12_Conductivity::~CR12_Conductivity(void) {}


void CR12_Conductivity::SetConductivity(su2double T, su2double rho, su2double mu, su2double cp) {


  Kt = 0.53316e-3*pow(T*T_ref, 1.7216);
  Kt = Kt*1e-3/K_ref;

}

void CR12_Conductivity::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {

  dktdrho_T = 0;
  dktdT_rho = 0.53316e-3 * 1.7216 * pow(T*T_ref, 0.7216);
  dktdT_rho = dktdT_rho *1e-3/K_ref * T_ref;
}
