/*!
 * \file variable_direct_electric.cpp
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

CPotentialVariable::CPotentialVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	Charge_Density = NULL;
	PlasmaRhoUGradient = NULL;
  
}

CPotentialVariable::CPotentialVariable(double val_potential, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;
  
	Residual_Old = new double [nVar];
	Residual_Sum = new double [nVar];
  
	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar< nVar; iVar++) {
		Solution[iVar] = val_potential;
		Solution_Old[iVar] = val_potential;
	}
	Charge_Density = new double [2];
  
	PlasmaRhoUGradient = new double* [3];
	for (iVar = 0; iVar < 3; iVar++)
		PlasmaRhoUGradient[iVar] = new double [nDim];
  
}

CPotentialVariable::~CPotentialVariable(void) {
  unsigned short iVar;
  
	if (Charge_Density != NULL) delete [] Charge_Density;
  if (PlasmaRhoUGradient != NULL) {
    for (iVar = 0; iVar < 3; iVar++)
      delete PlasmaRhoUGradient[iVar];
    delete [] PlasmaRhoUGradient;
  }
  
}