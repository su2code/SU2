/*!
 * \file variable_direct_heat.cpp
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

CHeatVariable::CHeatVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
}

CHeatVariable::CHeatVariable(double *val_heat, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;
  
  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
	/*--- Allocate residual structures ---*/
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
  
	/*--- Allocate direct solution container for adjoint problem ---*/
	Solution_Direct = new double[nVar];
  
	/*--- Allocate aux gradient vector ---*/
	Grad_AuxVar = new double [nDim];
  
	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_heat[iVar];
		Solution_Old[iVar] = val_heat[iVar];
		Solution_Direct[iVar] = 0.0;
	}
  
}

CHeatVariable::~CHeatVariable(void) {
  
  if (Solution_Direct != NULL) delete [] Solution_Direct;
  
}