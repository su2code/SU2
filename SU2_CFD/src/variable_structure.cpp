/*!
 * \file variable_structure.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.3
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

unsigned short CVariable::nDim = 0;

CVariable::CVariable(void) { }

CVariable::CVariable(unsigned short val_nvar, CConfig *config) {

  /*--- Initialize the number of solution variables. This version
   of the constructor will be used primarily for converting the
   restart files into solution files (SU2_SOL). ---*/
	nVar = val_nvar;
  
	/*--- Allocate the solution array - here it is also possible
	 to allocate some extra flow variables that do not participate
	 in the simulation ---*/
	Solution = new double [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Solution[iVar] = 0.0;
  
}

CVariable::CVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config) {
	unsigned short iVar, iDim;
	
	/*--- Initializate the number of dimension and number of variables ---*/
	nDim = val_ndim;
	nVar = val_nvar;
  
	/*--- Allocate solution, solution old, residual and gradient 
	 which is common for all the problems, here it is also possible 
	 to allocate some extra flow variables that do not participate 
	 in the simulation ---*/
	Solution = new double [nVar];
	
	for (iVar = 0; iVar < nVar; iVar++)
		Solution[iVar] = 0.0;

	Solution_Old = new double [nVar];

	Residual = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Residual[iVar] = 0.0;
	
	Gradient = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Gradient[iVar] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim ++)
			Gradient[iVar][iDim] = 0.0;
	}
	
	if (config->GetUnsteady_Simulation() != NO) {
		Solution_time_n = new double [nVar];
		Solution_time_n1 = new double [nVar];
	}
	
}

CVariable::~CVariable(void) {
	unsigned short iVar;
	
	delete [] Solution;
	delete [] Residual;
	delete [] Solution_Old;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Gradient[iVar];
	delete [] Gradient;
  
/*
  if (Solution != NULL) delete [] Solution;
	if (Solution_Old != NULL) delete [] Solution_Old;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;
  
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_Max != NULL) delete [] Solution_Max;
	if (Solution_Min != NULL) delete [] Solution_Min;
	if (Grad_AuxVar != NULL) delete [] Grad_AuxVar;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Residual != NULL) delete [] Residual;
	if (Res_Conv != NULL) delete [] Res_Conv;
	if (Res_Visc != NULL) delete [] Res_Visc;
	if (Res_Sour != NULL) delete [] Res_Sour;
	if (Res_TruncError != NULL) delete [] Res_TruncError;
	if (Residual_Old != NULL) delete [] Residual_Old;
	if (Residual_Sum != NULL) delete [] Residual_Sum;
  
  if (Gradient != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Gradient[iVar];
    delete [] Gradient;
  }
  
  if (Res_Visc_RK != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Res_Visc_RK[iVar];
    delete [] Res_Visc_RK;
  }
  
  //  if (LimiterPrimitive != NULL) {
  //    for (iVar = 0; iVar < nVar; iVar++)
  //      delete LimiterPrimitive[iVar];
  //    delete [] LimiterPrimitive;
  //  }
*/
  
}

void CVariable::AddUnd_Lapl(double *val_und_lapl) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Undivided_Laplacian[iVar] += val_und_lapl[iVar];
}

void CVariable::SubtractUnd_Lapl(double *val_und_lapl) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Undivided_Laplacian[iVar] -= val_und_lapl[iVar];
}

void CVariable::SubtractUnd_Lapl(unsigned short val_var, double val_und_lapl) {
	Undivided_Laplacian[val_var] -= val_und_lapl;
}

void CVariable::SetUnd_LaplZero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Undivided_Laplacian[iVar] = 0.0;
}

void CVariable::SetUnd_Lapl(unsigned short val_var, double val_und_lapl) {
		Undivided_Laplacian[val_var] = val_und_lapl;
}

void CVariable::GetResidual(double *val_residual) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Residual[iVar];
}

void CVariable::SetSolution(double *val_solution) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Solution[iVar] = val_solution[iVar];
}

void CVariable::Set_OldSolution(void) {
	unsigned short iVar;
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution_Old[iVar] = Solution[iVar];
	}
}

void CVariable::AddConservativeSolution(unsigned short val_var, double val_solution,
		double val_density, double val_density_old, double lowerlimit, double upperlimit) {
	Solution[val_var] = min(max((Solution_Old[val_var]*val_density_old + val_solution)/val_density,
			lowerlimit),upperlimit);
}

void CVariable::Set_Solution(void) {
	unsigned short iVar;
	for (iVar = 0; iVar < nVar; iVar++)
		 Solution[iVar] = Solution_Old[iVar];
}

void CVariable::Set_Solution_time_n(void) {
	unsigned short iVar;
	for (iVar = 0; iVar < nVar; iVar++)
		Solution_time_n[iVar] = Solution[iVar];
}

void CVariable::Set_Solution_time_n1(void) {
	unsigned short iVar;
	for (iVar = 0; iVar < nVar; iVar++)
		Solution_time_n1[iVar] = Solution_time_n[iVar];
}

void CVariable::AddRes_TruncError(double *val_truncation_error) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_TruncError[iVar] += val_truncation_error[iVar];
}

void CVariable::SubtractRes_TruncError(double *val_truncation_error) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_TruncError[iVar] -= val_truncation_error[iVar];
}

double CVariable::GetProjRes(double *val_vector) {
	double ProjRes = 0;
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		ProjRes += Residual[iDim+1]*val_vector[iDim];
	return ProjRes;
}

void CVariable::SetResidual(double *val_residual) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Residual[iVar] = val_residual[iVar];
}

void CVariable::SetRes_Conv(double *val_residual) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_Conv[iVar] = val_residual[iVar];
}

void CVariable::SetRes_Visc(double *val_residual) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_Visc[iVar] = val_residual[iVar];
}

void CVariable::SetRes_Sour(double *val_residual) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_Sour[iVar] = val_residual[iVar];
}

void CVariable::SetRes_Visc_RK(double *val_residual, unsigned short iRKStep) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar][iRKStep] = val_residual[iVar];
}

void CVariable::SetResidual_Old(double *val_residual_old) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Residual_Old[iVar] = val_residual_old[iVar];
}

void CVariable::SetSolution_Old(double *val_solution_old) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Solution_Old[iVar] = val_solution_old[iVar];
}

void CVariable::AddResidual_Sum(double *val_residual) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Residual_Sum[iVar] += val_residual[iVar];
}

void CVariable::SetResidualZero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Residual[iVar] = 0.0;
}

void CVariable::Set_ResConv_Zero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_Conv[iVar] = 0.0;
}

void CVariable::Set_ResVisc_Zero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_Visc[iVar] = 0.0;
}

void CVariable::Set_ResSour_Zero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_Sour[iVar] = 0.0;
}

void CVariable::SetVelResidualZero(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Residual[iDim+1] = 0.0;
}

void CVariable::SetVel_ResVisc_Zero(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Res_Visc[iDim+1] = 0.0;
}

void CVariable::SetVel_ResSour_Zero(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Res_Sour[iDim+1] = 0.0;
}

void CVariable::SetVel_ResConv_Zero(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Res_Conv[iDim+1] = 0.0;
}

void CVariable::SetVelRes_TruncErrorZero(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Res_TruncError[iDim+1] = 0.0;
}

void CVariable::SetVelSolutionZero(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Solution[iDim+1] = 0.0;
}

void CVariable::SetVelSolutionOldZero(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Solution_Old[iDim+1] = 0.0;
}

void CVariable::SetVelSolutionOldDVector(void) { }

void CVariable::SetVelSolutionDVector(void) { }

void CVariable::SetSolutionZero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Solution[iVar] = 0.0;
}

void CVariable::SetResidualSumZero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Residual_Sum[iVar] = 0.0;
}

void CVariable::SetGradientZero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Gradient[iVar][iDim] = 0.0;
}

void CVariable::SetAuxVarGradientZero(void) {
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Grad_AuxVar[iDim] = 0.0;
}

void CVariable::SetGradient(double **val_gradient) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Gradient[iVar][iDim] = val_gradient[iVar][iDim];
}

void CVariable::SetRes_TruncErrorZero(void) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Res_TruncError[iVar] = 0.0;
}

void CVariable::GetResidual_Sum(double *val_residual) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = Residual_Sum[iVar];
}

void CVariable::GetRes_TruncError(double *val_trunc_error) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		val_trunc_error[iVar] = Res_TruncError[iVar];
}

CBaselineVariable::CBaselineVariable(void) : CVariable() { }

CBaselineVariable::CBaselineVariable(double *val_solution, unsigned short val_nvar, CConfig *config) : CVariable(val_nvar, config) {
  
	/*--- Solution initialization ---*/
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Solution[iVar] = val_solution[iVar];
  
}

CBaselineVariable::~CBaselineVariable(void) { }
