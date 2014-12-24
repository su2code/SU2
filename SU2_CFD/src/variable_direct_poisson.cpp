/*!
 * \file variable_direct_poisson.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios
 * \version 3.2.6 "eagle"
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

CPotentialVariable::CPotentialVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	Charge_Density = NULL;
  
}

CPotentialVariable::CPotentialVariable(double val_potential,
                                       unsigned short val_nDim,
                                       unsigned short val_nvar,
                                       CConfig *config) : CVariable(val_nDim,
                                                                    val_nvar,
                                                                    config) {
	unsigned short iVar;
  
	Residual_Old = new double [nVar];
	Residual_Sum = new double [nVar];
  
	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar< nVar; iVar++) {
		Solution[iVar] = val_potential;
		Solution_Old[iVar] = val_potential;
	}
	Charge_Density = new double [2];
  
}

CPotentialVariable::~CPotentialVariable(void) {
  
	if (Charge_Density != NULL) delete [] Charge_Density;
  
}
