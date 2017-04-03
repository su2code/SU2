/*!
 * \file variable_direct_turbulent.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 4.3.0 "Cardinal"
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/variable_structure.hpp"

CTurbVariable::CTurbVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	TS_Source = NULL;
  
}

CTurbVariable::CTurbVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
  
  unsigned short iVar;

  /*--- Array initialization ---*/
  
	TS_Source = NULL;
  
	/*--- Allocate space for the time spectral source terms ---*/
  
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		TS_Source = new su2double[nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			TS_Source[iVar] = 0.0;
	}
  
	/*--- Allocate space for the limiter ---*/
  
  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new su2double [nVar];
  Solution_Min = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
}

CTurbVariable::~CTurbVariable(void) { }

su2double CTurbVariable::GetmuT() { return muT; }

void CTurbVariable::SetmuT(su2double val_muT) { muT = val_muT; }

CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() { }

CTurbSAVariable::CTurbSAVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, config) {
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	/*--- Initialization of S-A variables ---*/
	Solution[0] = val_nu_tilde;		Solution_Old[0] = val_nu_tilde;
  
	/*--- Initialization of the eddy viscosity ---*/
	muT = val_muT;
  
	/*--- Allocate and initialize solution for the dual time strategy ---*/
	if (dual_time) {
		Solution_time_n[0]  = val_nu_tilde;
		Solution_time_n1[0] = val_nu_tilde;
	}

}

CTurbSAVariable::~CTurbSAVariable(void) {
  
  if (TS_Source != NULL) delete [] TS_Source;
  
}

CTurbMLVariable::CTurbMLVariable(void) : CTurbVariable() { }

CTurbMLVariable::CTurbMLVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, config) {
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	/*--- Initialization of S-A variables ---*/
	Solution[0] = val_nu_tilde;		Solution_Old[0] = val_nu_tilde;
  
	/*--- Initialization of the eddy viscosity ---*/
	muT = val_muT;
  
	/*--- Allocate and initialize solution for the dual time strategy ---*/
	if (dual_time) {
		Solution_time_n[0]  = val_nu_tilde;
		Solution_time_n1[0] = val_nu_tilde;
	}
  
}

CTurbMLVariable::~CTurbMLVariable(void) {
  
  if (TS_Source != NULL) delete [] TS_Source;
  
}

CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() { }

CTurbSSTVariable::CTurbSSTVariable(su2double val_kine, su2double val_omega, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar,
                                   su2double *constants, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	/*--- Initialization of variables ---*/
  
	Solution[0] = val_kine;     Solution_Old[0] = val_kine;
	Solution[1] = val_omega;	Solution_Old[1] = val_omega;
  
	sigma_om2 = constants[3];
	beta_star = constants[6];
  
	F1   = 1.0;
	F2   = 0.0;
	CDkw = 0.0;
  
	/*--- Initialization of eddy viscosity ---*/
  
	muT = val_muT;
  
	/*--- Allocate and initialize solution for the dual time strategy ---*/
  
	if (dual_time) {
		Solution_time_n[0]  = val_kine; Solution_time_n[1]  = val_omega;
		Solution_time_n1[0]  = val_kine; Solution_time_n1[1]  = val_omega;
	}
    
}

CTurbSSTVariable::~CTurbSSTVariable(void) {

  if (TS_Source != NULL) delete [] TS_Source;
  
}

void CTurbSSTVariable::SetBlendingFunc(su2double val_viscosity, su2double val_dist, su2double val_density) {
	unsigned short iDim;
	su2double arg2, arg2A, arg2B, arg1;
  
	/*--- Cross diffusion ---*/
  
	CDkw = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		CDkw += Gradient[0][iDim]*Gradient[1][iDim];
	CDkw *= 2.0*val_density*sigma_om2/Solution[1];
	CDkw = max(CDkw, pow(10.0, -20.0));
  
	/*--- F1 ---*/
  
        arg2A = sqrt(Solution[0])/(beta_star*Solution[1]*val_dist+EPS*EPS);
        arg2B = 500.0*val_viscosity / (val_density*val_dist*val_dist*Solution[1]+EPS*EPS);
	arg2 = max(arg2A, arg2B);
        arg1 = min(arg2, 4.0*val_density*sigma_om2*Solution[0] / (CDkw*val_dist*val_dist+EPS*EPS));
	F1 = tanh(pow(arg1, 4.0));
  
	/*--- F2 ---*/
  
	arg2 = max(2.0*arg2A, arg2B);
	F2 = tanh(pow(arg2, 2.0));
  
}


// swh
CTurbKEVariable::CTurbKEVariable(void) : CTurbVariable() { }

CTurbKEVariable::CTurbKEVariable(su2double val_kine, su2double val_epsi, su2double val_zeta, su2double val_f, su2double val_muT, su2double val_Tm, su2double val_Lm, unsigned short val_nDim, unsigned short val_nvar, su2double *constants, CConfig *config): CTurbVariable(val_nDim, val_nvar, config) {

//CTurbKEVariable::CTurbKEVariable(su2double val_kine, su2double val_epsi, su2double val_zeta, su2double val_muT, su2double val_Tm, su2double val_Lm, unsigned short val_nDim, unsigned short val_nvar, su2double *constants, CConfig *config): CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	/*--- Initialization of variables ---*/
	Solution[0] = val_kine; Solution_Old[0] = val_kine;
	Solution[1] = val_epsi;	Solution_Old[1] = val_epsi;
	Solution[2] = val_zeta; Solution_Old[2] = val_zeta;
	Solution[3] = val_f;	Solution_Old[3] = val_f;
	Tm  = val_Tm; //val_kine/val_epsi;
	Lm  = val_Lm; //pow(val_kine,1.5)/val_epsi;
  
	/*--- Initialization of eddy viscosity ---*/  
	muT = val_muT;
  
	/*--- Allocate and initialize solution for the dual time strategy ---*/
	if (dual_time) {
		Solution_time_n[0]  = val_kine; Solution_time_n[1]  = val_epsi;
		Solution_time_n1[0] = val_kine; Solution_time_n1[1] = val_epsi;
		Solution_time_n[2]  = val_zeta; Solution_time_n[3]  = val_f;
		Solution_time_n1[2] = val_zeta; Solution_time_n1[3] = val_f;
	}
    
}

CTurbKEVariable::~CTurbKEVariable(void) {

  if (TS_Source != NULL) delete [] TS_Source;
  
}

//void CTurbKEVariable::SetTLFunc(su2double val_viscosity, su2double val_dist, su2double val_density, su2double val_kine,su2double val_epsi,su2double val_zeta, su2double val_f, su2double StrainMag) {
void CTurbKEVariable::SetTLFunc(su2double val_viscosity, su2double val_dist, su2double val_density, su2double val_kine,su2double val_epsi,su2double val_zeta, su2double StrainMag, su2double VelInf, su2double L_Inf) {

	unsigned short iDim;
	su2double C_mu, C_T, C_L, C_eta;
	//        su2double Tm, Lm, nu, temp;
        su2double nu, temp, scalar_min;
  
        /*--- Molecular kinematic viscosity ---*/
        nu = val_viscosity/val_density;

        // move these, already in constants...
        //--- zeta-f ---//
	/*
        C_mu = 0.22;
        C_T = 6.0;
        C_L = 0.36;
        C_eta = 85.0;

        //--- Model time scale ---//
        temp = min(val_kine/val_epsi,0.6/(sqrt(6.0)*C_mu*StrainMag*val_zeta));
        Tm = max(temp,C_T*sqrt(nu/val_epsi));

        //--- Model length scale ---//
	temp = min(pow(val_kine,1.5)/val_epsi, sqrt(val_kine)/(sqrt(6.0)*C_mu*StrainMag*val_zeta));
        Lm = C_L * max(temp,C_eta*pow(pow(nu,3.0)/val_epsi,0.25));
*/

        //--- v2-f ---//
        C_mu = 0.22;
        C_T = 6.0;
        C_L = 0.23;
        C_eta = 70.0;

	//        scalar_min = 1.0E-12;
	scalar_min = 1.0E-8/(VelInf*VelInf); // setting based on tke min being 1e-
        su2double tke = max(val_kine, scalar_min*VelInf*VelInf);
        su2double tdr = max(val_epsi, scalar_min*VelInf*VelInf*VelInf/L_Inf);
        su2double v2 = max(val_zeta, 2.0/3.0*scalar_min*VelInf*VelInf);
	//        su2double S = max(StrainMag,scalar_min*VelInf/L_Inf);
        su2double S = max(StrainMag,1.0E-14);

        //--- Model time scale ---//
	//        temp = min(val_kine/val_epsi,0.6*val_kine/(sqrt(6.0)*C_mu*StrainMag*val_zeta));
	//        Tm = max(temp,C_T*sqrt(nu/val_epsi));
        temp = max(tke/tdr,C_T*sqrt(nu/tdr));
        Tm = min(temp,0.6*tke/(sqrt(6.0)*C_mu*S*v2));

        //--- Model length scale ---//
	temp = min(pow(tke,1.5)/tdr, pow(tke,1.5)/(sqrt(6.0)*C_mu*S*v2));
        Lm = C_L * max(temp,C_eta*pow(pow(nu,3.0)/tdr,0.25));


	// good here... cout<<" Lm in variable_direct: "<<Lm<<"\n";
  
}
