/*!
 * \file variable_adjoint_discrete.cpp
 * \brief Main subroutines for the discrete adjoint variable structure.
 * \author T. Albring
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

CDiscAdjVariable::CDiscAdjVariable() : CVariable(){

}

CDiscAdjVariable::CDiscAdjVariable(su2double* val_solution, unsigned short val_ndim,
                               unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config){

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  if (dual_time){
    DualTime_Derivative = new su2double[nVar];
    DualTime_Derivative_n = new su2double[nVar];
  }

  Solution_Direct = new su2double[nVar];

  Sensitivity = new su2double[nDim];

  unsigned short iVar,iDim;

  for (iDim = 0; iDim < nDim; iDim++){
    Sensitivity[iDim] = 0.0;
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution[iVar] = val_solution[iVar];
  }


  if (dual_time){
    for (iVar = 0; iVar < nVar; iVar++){
      Solution_time_n[iVar]  = 0.0;
      Solution_time_n1[iVar] = 0.0;
      DualTime_Derivative[iVar] = 0.0;
      DualTime_Derivative_n[iVar] = 0.0;
    }
  }
}


CDiscAdjVariable::~CDiscAdjVariable(){

}
