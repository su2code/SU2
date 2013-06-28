/*!
 * \file variable_direct.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
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

CVariable::CVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config) {
	unsigned short iVar;
	
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
	
}

CVariable::~CVariable(void) {
	
	delete [] Solution;
  
}

void CVariable::SetSolution(unsigned short val_var, double val_solution) { Solution[val_var] = val_solution; }

double *CVariable::GetSolution(void) { return Solution; }

double CVariable::GetSolution(unsigned short val_var) { return Solution[val_var]; }

void CVariable::SetSolution(double *val_solution) {
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		Solution[iVar] = val_solution[iVar];
}

CBaselineVariable::CBaselineVariable(void) : CVariable() { }

CBaselineVariable::CBaselineVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Solution[iVar] = val_solution[iVar];

}

CBaselineVariable::~CBaselineVariable(void) { }