/*!
 * \file variable_adjoint.cpp
 * \brief Definition of the solution fields.
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

#include "../include/variable_structure.hpp"

CAdjPotentialVariable::CAdjPotentialVariable(void) : CVariable() {}

CAdjPotentialVariable::~CAdjPotentialVariable(void) { }

CAdjPotentialVariable::CAdjPotentialVariable(double val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) 
: CVariable(val_ndim, val_nvar, config) {

}

CAdjEulerVariable::CAdjEulerVariable(void) : CVariable() { }

CAdjEulerVariable::CAdjEulerVariable(double val_psirho, double *val_phi, double val_psie, unsigned short val_ndim, 
																		 unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim;
	
	incompressible = config->GetIncompressible();
	
	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];
	
	/*--- Allocate undivided laplacian, limiter, and auxiliar gradient ---*/
	Undivided_Laplacian = new double [nVar];
	Limiter = new double [nVar];
	Grad_AuxVar = new double [nDim];
	
	/*--- Allocate and initializate tructation  error ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;
	
	/*--- Allocate and initializate projection vector for wall boundary condition ---*/
	ForceProj_Vector = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
	
	/*--- Allocate and initializate interior boundary jump vector for near field boundary condition ---*/
	IntBoundary_Jump = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		IntBoundary_Jump[iVar] = 0.0;
	
	/*--- Allocate and initializate solution ---*/
	if (incompressible) {
		Solution[0] = 0.0; 	Solution_Old[0] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = 0.0;
			Solution_Old[iDim+1] = 0.0;
		}
	}
	else {
		Solution[0] = val_psirho; 	Solution_Old[0] = val_psirho;
		Solution[nVar-1] = val_psie; Solution_Old[nVar-1] = val_psie;	
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = val_phi[iDim];
			Solution_Old[iDim+1] = val_phi[iDim];
		}
	}
}

CAdjEulerVariable::CAdjEulerVariable(double *val_solution, unsigned short val_ndim, 
					 unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim;
	
	incompressible = config->GetIncompressible();

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];  Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];
	
	/*--- Allocate undivided laplacian, limiter and auxiliar gradient ---*/
	Undivided_Laplacian = new double [nVar];
	Limiter = new double [nVar];
	Grad_AuxVar = new double [nDim];
	
	/*--- Allocate and initializate tructation  error ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;
	
	/*--- Allocate and initializate projection vector for wall boundary condition ---*/
	ForceProj_Vector = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
	
	/*--- Allocate and initializate interior boundary jump vector for near field boundary condition ---*/
	IntBoundary_Jump = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		IntBoundary_Jump[iVar] = 0.0;

	/*--- Allocate and initializate solution (including dual time strategy) ---*/
	Solution_time_n = new double [nVar];
	Solution_time_n1 = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}
}

CAdjEulerVariable::~CAdjEulerVariable(void) {
	unsigned short iVar;
	
	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	
	delete [] Undivided_Laplacian;
	delete [] Limiter;
	delete [] Grad_AuxVar;
	delete [] TruncationError;
	delete [] ForceProj_Vector;
	delete [] IntBoundary_Jump;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;
}

void CAdjEulerVariable::SetTheta(double val_density, double *val_velocity, double val_enthalpy) {
	unsigned short iDim;
	
	Theta = val_density*Solution[0];
	Theta += val_density*val_enthalpy*Solution[nDim+1];
	
	for (iDim = 0; iDim < nDim; iDim++)
		Theta += val_density*val_velocity[iDim]*Solution[iDim+1];
}

CAdjNSVariable::CAdjNSVariable(void) : CAdjEulerVariable() { }

CAdjNSVariable::CAdjNSVariable(double *val_solution, unsigned short val_ndim, 
								   unsigned short val_nvar, CConfig *config) : CAdjEulerVariable(val_solution, val_ndim, val_nvar, config) {

}

CAdjNSVariable::CAdjNSVariable(double val_psirho, double *val_phi, double val_psie, 
								   unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CAdjEulerVariable(val_psirho, val_phi, val_psie, val_ndim, val_nvar, config) {


}

CAdjNSVariable::~CAdjNSVariable(void) {
	unsigned short iVar;
	
	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	
	delete [] Undivided_Laplacian;
	delete [] Limiter;
	delete [] Grad_AuxVar;
	delete [] TruncationError;
	delete [] ForceProj_Vector;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;
}

void CAdjNSVariable::SetTheta(double val_density, double *val_velocity, double val_enthalpy) {
	Theta = val_density*Solution[0];
	Theta += val_density*val_enthalpy*Solution[nDim+1];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Theta += val_density*val_velocity[iDim]*Solution[iDim+1];
}

CAdjTurbVariable::CAdjTurbVariable(void) : CVariable() {}

CAdjTurbVariable::~CAdjTurbVariable(void) {}

CAdjTurbVariable::CAdjTurbVariable(double val_psinu_inf, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) 
			: CVariable(val_ndim, val_nvar, config) {

	// Initialization of variables
	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_psinu_inf; 
		Solution_Old[iVar] = val_psinu_inf;
	}

	Residual_Old = new double [nVar];

	if (config->GetKind_SlopeLimit() != NONE) Limiter = new double [nVar];
}

CAdjPlasmaVariable::CAdjPlasmaVariable(void) : CVariable() { }

CAdjPlasmaVariable::CAdjPlasmaVariable(double val_psirho, double *val_phi, double val_psie, unsigned short val_ndim, 
																		 unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim;
		
	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];
	
	/*--- Allocate undivided laplacian, limiter, and auxiliar gradient ---*/
	Limiter = new double [nVar];
	Grad_AuxVar = new double [nDim];
	
	/*--- Allocate and initializate tructation  error ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;
	
	/*--- Allocate and initializate projection vector for wall boundary condition ---*/
	ForceProj_Vector = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
	
	/*--- Allocate and initializate solution ---*/
	unsigned short iSpecies, nSpecies = 7, nDiatomics = 0, loc;
	
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
		
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		
		Solution[loc+0] = val_psirho; 	Solution_Old[loc+0] = val_psirho;
		Solution[loc+nDim+1] = val_psie; Solution_Old[loc+nDim+1] = val_psie;	
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[loc+iDim+1] = val_phi[iDim];
			Solution_Old[loc+iDim+1] = val_phi[iDim];
		}
	}
	
}

CAdjPlasmaVariable::CAdjPlasmaVariable(double *val_solution, unsigned short val_ndim, 
																		 unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim;
	
	incompressible = config->GetIncompressible();
	
	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];
	
	/*--- Allocate undivided laplacian, limiter and auxiliar gradient ---*/
	Limiter = new double [nVar];
	Grad_AuxVar = new double [nDim];
	
	/*--- Allocate and initializate tructation  error ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;
	
	/*--- Allocate and initializate projection vector for wall boundary condition ---*/
	ForceProj_Vector = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
		
	/*--- Allocate and initializate solution (including dual time strategy) ---*/
	Solution_time_n = new double [nVar];
	Solution_time_n1 = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}
}

CAdjPlasmaVariable::~CAdjPlasmaVariable(void) {
	unsigned short iVar;
	
	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	
	delete [] Undivided_Laplacian;
	delete [] Limiter;
	delete [] Grad_AuxVar;
	delete [] TruncationError;
	delete [] ForceProj_Vector;
	delete [] IntBoundary_Jump;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;
}

void CAdjPlasmaVariable::SetTheta(double val_density, double *val_velocity, double val_enthalpy) {
	unsigned short iDim;
	
	Theta = val_density*Solution[0];
	Theta += val_density*val_enthalpy*Solution[nDim+1];
	
	for (iDim = 0; iDim < nDim; iDim++)
		Theta += val_density*val_velocity[iDim]*Solution[iDim+1];
}

CAdjLevelSetVariable::CAdjLevelSetVariable(void) : CVariable() {}

CAdjLevelSetVariable::CAdjLevelSetVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;
	
	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];
	
	/*--- Allocate limiter (upwind)---*/
	if (config->GetKind_SlopeLimit() != NONE) Limiter = new double [nVar];
	
}

CAdjLevelSetVariable::CAdjLevelSetVariable(double val_levelset, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar,config) {
	unsigned short iVar;
	
	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];
	
	/*--- Allocate limiter (upwind)---*/
	if (config->GetKind_SlopeLimit() != NONE) Limiter = new double [nVar];
	
	/*--- Solution and old solution initialization ---*/
	Solution[0] = val_levelset;		Solution_Old[0] = val_levelset;
	
	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) 
			|| (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
		Solution_time_n[0] = val_levelset;
		Solution_time_n1[0] = val_levelset;
	}
	
}

CAdjLevelSetVariable::~CAdjLevelSetVariable(void) {
	unsigned short iVar;
	
	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Limiter != NULL) delete [] Limiter;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;
}
