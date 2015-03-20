/*!
 * \file variable_linearized_mean.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

CLinEulerVariable::CLinEulerVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	DeltaU = NULL;
  ForceProj_Vector = NULL;
  
}

CLinEulerVariable::CLinEulerVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
	
  /*--- Array initialization ---*/
	DeltaU = NULL;
  ForceProj_Vector = NULL;
  
	/*--- Allocate structures ---*/
	Residual_Sum = new double [nVar];
	Residual_Old = new double [nVar];
	
	Undivided_Laplacian = new double [nVar];
	Limiter = new double [nVar];
	
	Grad_AuxVar = new double [nDim];
	Solution_time_n = new double [nVar];
	Solution_time_n1 = new double [nVar];
	
	Res_TruncError = new double [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_TruncError[iVar] = 0.0;
	
	/*--- Initialization of variables ---*/
	ForceProj_Vector = new double [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
	
	for (unsigned short iVar = 0; iVar < nDim+2; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}
}

CLinEulerVariable::CLinEulerVariable(double val_deltarho, double *val_deltavel,
                                     double val_deltae, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
	
  /*--- Array initialization ---*/
	DeltaU = NULL;
  ForceProj_Vector = NULL;
  
	/*--- Allocate structures ---*/
	Residual_Sum = new double [nVar];
	Residual_Old = new double [nVar];
	
	Undivided_Laplacian = new double [nVar];
	Limiter = new double [nVar];
	
	Grad_AuxVar = new double [nDim];
	Solution_time_n = new double [nVar];
	Solution_time_n1 = new double [nVar];
	
	Res_TruncError = new double [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_TruncError[iVar] = 0.0;
	
	/*--- Initialization of variables ---*/
	ForceProj_Vector = new double [nDim];
	Solution[0] = val_deltarho; 	Solution_Old[0] = val_deltarho;
	Solution[nVar-1] = val_deltae; Solution_Old[nVar-1] = val_deltae;
	
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		Solution[iDim+1] = val_deltavel[iDim];
		Solution_Old[iDim+1] = val_deltavel[iDim];
	}
}

CLinEulerVariable::~CLinEulerVariable(void) {
  
	if (DeltaU != NULL) delete [] DeltaU;
	if (ForceProj_Vector != NULL) delete [] ForceProj_Vector;
  
}

void CLinEulerVariable::SetDeltaPressure(double *val_velocity, double Gamma) {
	
	double Mod_Vel = 0.0;
	double Vel_dot_DeltaRhoVel = 0.0;
	for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		Mod_Vel += val_velocity[iDim] * val_velocity[iDim];
		Vel_dot_DeltaRhoVel += val_velocity[iDim] * Solution[iDim+1];
	}
	
	DeltaPressure = 0.5*Solution[0]*Mod_Vel+(Gamma-1.0)*(Solution[nVar-1]-Vel_dot_DeltaRhoVel);
}
