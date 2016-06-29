/*!
 * \file variable_direct_wave.cpp
 * \brief Definition of the solution fields.
 * \author T. Economon, F. Palacios
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

CWaveVariable::CWaveVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
}

CWaveVariable::CWaveVariable(su2double *val_wave, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
	unsigned short iVar;
  
  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
	/*--- Allocate residual structures ---*/
	Residual_Sum = new su2double [nVar]; Residual_Old = new su2double [nVar];
  
	/*--- Allocate direct solution container for adjoint problem ---*/
	Solution_Direct = new su2double[nVar];
  
	/*--- Allocate aux gradient vector ---*/
	Grad_AuxVar = new su2double [nDim];
  
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
