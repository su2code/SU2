/*!
 * \file variable_adjoint_plasma.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.8
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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

CAdjPlasmaVariable::CAdjPlasmaVariable(void) : CVariable() { }

CAdjPlasmaVariable::CAdjPlasmaVariable(double val_psirho, double *val_phi, double val_psie, double val_psievib, unsigned short val_ndim,
                                       unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim;
  
	/*--- Allocate residual structures ---*/
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	
	/*--- Allocate undivided laplacian, limiter, and auxiliar gradient ---*/
	Limiter = new double [nVar];
	Grad_AuxVar = new double [nDim];
	
	/*--- Allocate and initializate tructation  error ---*/
	Res_TruncError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_TruncError[iVar] = 0.0;
	
	/*--- Allocate and initializate projection vector for wall boundary condition ---*/
	ForceProj_Vector = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
	
	/*--- Allocate and initializate solution ---*/
	unsigned short iSpecies, nSpecies, nDiatomics, loc;
  
  nSpecies = config->GetnSpecies();
  nDiatomics = config->GetnDiatomics();
	
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
		
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		
		Solution[loc+0] = val_psirho; 	Solution_Old[loc+0] = val_psirho;
		Solution[loc+nDim+1] = val_psie; Solution_Old[loc+nDim+1] = val_psie;
    if (iSpecies < nDiatomics)
      Solution[loc+nDim+2] = val_psievib; Solution_Old[loc+nDim+2] = val_psievib;
    
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[loc+iDim+1] = val_phi[iDim];
			Solution_Old[loc+iDim+1] = val_phi[iDim];
		}
	}
	
}

CAdjPlasmaVariable::CAdjPlasmaVariable(double *val_solution, unsigned short val_ndim,
                                       unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim;
  
	/*--- Allocate residual structures ---*/
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	
	/*--- Allocate undivided laplacian, limiter and auxiliar gradient ---*/
	Limiter = new double [nVar];
	Grad_AuxVar = new double [nDim];
	
	/*--- Allocate and initializate tructation  error ---*/
	Res_TruncError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_TruncError[iVar] = 0.0;
	
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

CAdjPlasmaVariable::~CAdjPlasmaVariable(void) { }
