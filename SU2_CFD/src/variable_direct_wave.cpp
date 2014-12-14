/*!
 * \file variable_direct_wave.cpp
 * \brief Definition of the solution fields.
 * \author T. Economon
 * \version 3.2.5 "eagle"
 *
 * Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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

CWaveVariable::CWaveVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
}

CWaveVariable::CWaveVariable(double *val_wave, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
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
		Solution[iVar] = val_wave[iVar];
		Solution_Old[iVar] = val_wave[iVar];
		Solution_Direct[iVar] = 0.0;
	}
  
}

CWaveVariable::~CWaveVariable(void) {
  
	if (Solution_Direct != NULL) delete [] Solution_Direct;
  
}
