/*!
 * \file variable_direct_turbulent.cpp
 * \brief Definition of the solution fields.
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

#include "../include/variable_structure.hpp"

CTurbVariable::CTurbVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	TS_Source = NULL;
  
}

CTurbVariable::CTurbVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
  
  /*--- Array initialization ---*/
	TS_Source = NULL;
  
}

CTurbVariable::~CTurbVariable(void) {
  
	if (TS_Source != NULL) delete [] TS_Source;
  
}

double CTurbVariable::GetmuT(){ return muT; }

void CTurbVariable::SetmuT(double val_muT){ muT = val_muT; }

CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() { }

CTurbSAVariable::CTurbSAVariable(double val_nu_tilde, double val_muT, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;
  
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
  
	/*--- Allocate space for the time spectral source terms ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		TS_Source = new double[nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			TS_Source[iVar] = 0.0;
	}
  
	/*--- Allocate space for the limiter ---*/
	if (config->GetKind_SlopeLimit_Turb() != NONE) {
		Limiter = new double [nVar];
		Solution_Max = new double [nVar];
		Solution_Min = new double [nVar];
	}
  
}

CTurbSAVariable::~CTurbSAVariable(void) {
  
  if (TS_Source != NULL) delete [] TS_Source;
  
}

CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() { }

CTurbSSTVariable::CTurbSSTVariable(double val_kine, double val_omega, double val_muT, unsigned short val_ndim, unsigned short val_nvar,
                                   double *constants, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {
  
	// Initialization of variables
	Solution[0] = val_kine;     Solution_Old[0] = val_kine;
	Solution[1] = val_omega;	Solution_Old[1] = val_omega;
  
	sigma_om2 = constants[3];
	beta_star = constants[6];
  
	F1   = 1.0;
	F2   = 0.0;
	CDkw = 0.0;
  
	// Initialization of eddy viscosity
	muT = val_muT;
}

CTurbSSTVariable::~CTurbSSTVariable(void) { }

void CTurbSSTVariable::SetBlendingFunc(double val_viscosity, double val_dist, double val_density){
	unsigned short iDim;
	double arg2, arg2A, arg2B, arg1;
  
	// Cross diffusion
	CDkw = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		CDkw += Gradient[0][iDim]*Gradient[1][iDim];
	CDkw *= 2.0*val_density*sigma_om2/Solution[1];
	CDkw = max(CDkw,pow(10.0,-20.0));
  
	// F1
	arg2A = sqrt(Solution[0])/(beta_star*Solution[1]*val_dist);
	arg2B = 500.0*val_viscosity/(val_density*val_dist*val_dist*Solution[1]);
	arg2 = max(arg2A,arg2B);
	arg1 = min(arg2,4.0*val_density*sigma_om2*Solution[0]/(CDkw*val_dist*val_dist));
	F1 = tanh(pow(arg1,4.0));
  
	// F2
	arg2 = max(2.0*arg2A,arg2B);
	F2 = tanh(pow(arg2,2.0));
}

double CTurbSSTVariable::GetF1blending(){
	return F1;
}

double CTurbSSTVariable::GetF2blending(){
	return F2;
}

double CTurbSSTVariable::GetCrossDiff(){
	return CDkw;
}