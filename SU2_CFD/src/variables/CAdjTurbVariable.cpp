/*!
 * \file CAdjTurbVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
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

#include "../../include/variables/CAdjTurbVariable.hpp"

CAdjTurbVariable::CAdjTurbVariable(void) : CVariable() {

  /*--- Array initialization ---*/

  dmuT_dUTvar = NULL;
  dRTstar_dUTvar = NULL;
  dFT_dUTvar = NULL;
  EddyViscSens = NULL;

}

CAdjTurbVariable::CAdjTurbVariable(su2double val_psinu_inf, unsigned short val_nDim, unsigned short val_nvar,
                                   CConfig *config) : CVariable(val_nDim, val_nvar, config) {
  unsigned short iVar;

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

  Residual_Old = new su2double [nVar];

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;

}

CAdjTurbVariable::~CAdjTurbVariable(void) {

  if (dmuT_dUTvar   != NULL) delete [] dmuT_dUTvar;
  if (EddyViscSens  != NULL) delete [] EddyViscSens;

}
