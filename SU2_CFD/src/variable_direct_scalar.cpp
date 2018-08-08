/*!
 * \file variable_direct_scalar.cpp
 * \brief Definition of the scalar equation variables at each vertex.
 * \author T. Economon
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/variable_structure.hpp"

CScalarVariable::CScalarVariable(void) : CVariable() { }

CScalarVariable::CScalarVariable(su2double *val_scalar,
                                 unsigned short val_nDim,
                                 unsigned short val_nvar,
                                 CConfig *config) : CVariable(val_nDim,
                                                              val_nvar,
                                                              config) {
  
  unsigned short iVar, iMesh, nMGSmooth = 0;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Initialization of variables ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar]     = val_scalar[iVar];
    Solution_Old[iVar] = val_scalar[iVar];
  }
  
  /*--- Allocate and initialize solution for the dual time strategy. ---*/
  
  if (dual_time) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution_time_n[iVar]     = val_scalar[iVar];
      Solution_time_n1[iVar] = val_scalar[iVar];
    }
  }
  
  /*--- Always allocate the slope limiter and necessary aux. variables. ---*/
  
  Limiter = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new su2double[nVar];
  Solution_Min = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
  /*--- Allocate residual structures in case of multigrid. ---*/
  
  Res_TruncError = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }
  
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new su2double[nVar];
    Residual_Old = new su2double[nVar];
  }
  
}

CScalarVariable::~CScalarVariable(void) { }

