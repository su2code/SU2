/*!
 * \file variable_direct_heat.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CHeatVariable::CHeatVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
  Solution_Direct = NULL;
  
}

CHeatVariable::CHeatVariable(su2double val_Heat, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {

  unsigned short iVar, iMesh, nMGSmooth = 0;
  bool low_fidelity = config->GetLowFidelitySim();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Initialization of heat variable ---*/
  Solution[0] = val_Heat;		Solution_Old[0] = val_Heat;

  /*--- Allocate residual structures ---*/

  Res_TruncError = new su2double [nVar];

  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }

  /*--- Only for residual smoothing (multigrid) ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

  if ((nMGSmooth > 0) || low_fidelity) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }

  /*--- Allocate and initialize solution for dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[0]  = val_Heat;
    Solution_time_n1[0] = val_Heat;
  }
}

CHeatVariable::~CHeatVariable(void) {  }
