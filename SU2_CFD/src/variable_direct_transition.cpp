/*!
 * \file variable_direct_transition.cpp
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

CTransLMVariable::CTransLMVariable(void) : CTurbVariable() {}

CTransLMVariable::CTransLMVariable(double val_nu_tilde, double val_intermittency, double val_REth,  unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {
  
	// Initialization of variables
	Solution[0] = val_intermittency; Solution_Old[0] = val_intermittency;
	Solution[1] = val_REth;          Solution_Old[1] = val_REth;
  
}

CTransLMVariable::~CTransLMVariable(void) { }

void CTransLMVariable::SetGammaEff() {
  
	/* -- Correction for separation-induced transition -- */
	Solution[0] = max(Solution[0],gamma_sep);
  
}