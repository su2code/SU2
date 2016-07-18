/*!
 * \file variable_direct_turbulent.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 4.2.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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
  
	if (config->GetUnsteady_Simulation() == SPECTRAL_METHOD) {
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
