/*!
 * \file variable_direct_elasticity.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

CFEAVariable::CFEAVariable(void) : CVariable() { }

CFEAVariable::CFEAVariable(double *val_fea, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
	unsigned short iVar, iDim;
  
	/*--- Allocate residual structures ---*/
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
  
	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_fea[iVar];
		Solution_Old[iVar] = val_fea[iVar];
	}
  
  /*--- Allocate stress tensor ---*/
	Stress = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Stress[iDim] = new double [nDim];
  
}

CFEAVariable::~CFEAVariable(void) {
  unsigned short iDim;

	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Stress[iDim];
	delete [] Stress;
  delete [] Residual_Sum;
  delete [] Residual_Old;
  
}
