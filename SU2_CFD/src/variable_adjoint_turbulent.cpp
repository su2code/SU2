/*!
 * \file variable_adjoint_turbulent.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
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

CAdjTurbVariable::CAdjTurbVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	dmuT_dUTvar = NULL;
	dRTstar_dUTvar = NULL;
	dFT_dUTvar = NULL;
	EddyViscSens = NULL;
  
}

CAdjTurbVariable::CAdjTurbVariable(double val_psinu_inf, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
  
  /*--- Array initialization ---*/
	dmuT_dUTvar = NULL;
	dRTstar_dUTvar = NULL;
	dFT_dUTvar = NULL;
	EddyViscSens = NULL;
  
	/*--- Initialization of variables ---*/
	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_psinu_inf;
		Solution_Old[iVar] = val_psinu_inf;
	}
  
	Residual_Old = new double [nVar];
  
	if (config->GetKind_SlopeLimit() != NONE) Limiter = new double [nVar];
  
	// Hybrid adjoint
  //	if (config->GetKind_Adjoint() == HYBRID) {
  //		unsigned short totalnVar;
  //		totalnVar = nDim + 2 + nVar;
  //	// Allocate
  //		dmuT_dUTvar = new double [totalnVar];
  //		dRTstar_dUTvar = new double* [nNeigh+1];
  //		for (unsigned short iNeigh = 0; iNeigh < (nNeigh+1); iNeigh++) {
  //			dRTstar_dUTvar[iNeigh] = new double [totalnVar];
  //		}
  //
  //		dFT_dUTvar = new double* [nNeigh+1]; // In theory only needed at boundaries
  //		for (unsigned short iNeigh = 0; iNeigh < (nNeigh+1); iNeigh++) {
  //			dFT_dUTvar[iNeigh] = new double [totalnVar];
  //		}
  //	// Initialise
  //		for (unsigned short iVar = 0; iVar < totalnVar; iVar++) {
  //			dmuT_dUTvar[iVar] = 0.0;
  //		}
  //		for (unsigned short iNeigh = 0; iNeigh < (nNeigh+1); iNeigh++) {
  //			for (unsigned short iVar = 0; iVar < totalnVar; iVar++) {
  //				dRTstar_dUTvar[iNeigh][iVar] = 0.0;
  //				dFT_dUTvar[iNeigh][iVar] = 0.0;
  //			}
  //		}
  //	}
  
	if (config->GetKind_Adjoint() == HYBRID) {
		unsigned short nTotalVar;
		nTotalVar = nDim + 2 + nVar;
		EddyViscSens = new double[nTotalVar];
    for (unsigned short iVar = 0; iVar < nTotalVar; iVar++)
      EddyViscSens[iVar] = 0.0;
	}
  
}

CAdjTurbVariable::~CAdjTurbVariable(void) {
  
	if (dmuT_dUTvar   != NULL) delete [] dmuT_dUTvar;
	if (EddyViscSens  != NULL) delete [] EddyViscSens;
  
  //	if (dRTstar_dUTvar != NULL) {
  //    delete [] dRTstar_dUTvar;
  //  }
  //	if (dFT_dUTvar != NULL) {
  //    delete [] dFT_dUTvar;
  //  }
  
}