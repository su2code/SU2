/*!
 * \file variable_direct_transition.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.2 "eagle"
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

CTransLMVariable::CTransLMVariable(void) : CTurbVariable() {}

CTransLMVariable::CTransLMVariable(double val_nu_tilde, double val_intermittency, double val_REth,  unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_nDim, val_nvar,config) {
  
	// Initialization of variables
	Solution[0] = val_intermittency; Solution_Old[0] = val_intermittency;
	Solution[1] = val_REth;          Solution_Old[1] = val_REth;
  
}

CTransLMVariable::~CTransLMVariable(void) { }

void CTransLMVariable::SetGammaEff() {
  
	/* -- Correction for separation-induced transition -- */
	Solution[0] = max(Solution[0],gamma_sep);
  
}
