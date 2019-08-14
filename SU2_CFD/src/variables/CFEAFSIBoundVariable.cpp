/*!
 * \file CFEABoundVariable.cpp
 * \brief Definition of the variables for FEM elastic structural problems.
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/variables/CFEAFSIBoundVariable.hpp"

CFEAFSIBoundVariable::CFEAFSIBoundVariable(void) : CFEABoundVariable() {

  FlowTraction          = NULL;    // Nodal traction due to the fluid (fsi)
  FlowTraction_n        = NULL;    // Nodal traction due to the fluid (fsi) at time n (for gen-alpha methods)

}

CFEAFSIBoundVariable::CFEAFSIBoundVariable(su2double *val_fea, unsigned short val_nDim, unsigned short val_nvar,
                                          CConfig *config) : CFEABoundVariable(val_fea, val_nDim, val_nvar, config) {

  unsigned short iVar;
  bool gen_alpha = (config->GetKind_TimeIntScheme_FEA() == GENERALIZED_ALPHA);

  /*--- Flow traction ---*/
  FlowTraction = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    FlowTraction[iVar]   = 0.0;
  }

  /*--- Generalized alpha integration method requires storing the old residuals ---*/
  FlowTraction_n = NULL;
  if (gen_alpha) {
    FlowTraction_n = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      FlowTraction_n[iVar] = 0.0;
    }
  }

}

CFEAFSIBoundVariable::~CFEAFSIBoundVariable(void) {

  if (FlowTraction          != NULL) delete [] FlowTraction;
  if (FlowTraction_n         != NULL) delete [] FlowTraction_n;

}
