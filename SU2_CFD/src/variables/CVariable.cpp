/*!
 * \file CVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
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

#include "../../include/variables/CVariable.hpp"

unsigned short CVariable::nDim = 0;

CVariable::CVariable(void) {

  /*--- Array initialization ---*/
  Solution = NULL;
  Solution_Old = NULL;
  Solution_time_n = NULL;
  Solution_time_n1 = NULL;
  Gradient = NULL;
  Limiter = NULL;
  Solution_Max = NULL;
  Solution_Min = NULL;
  Grad_AuxVar = NULL;
  Undivided_Laplacian = NULL;
  Res_TruncError = NULL;
  Residual_Old = NULL;
  Residual_Sum = NULL;
  Solution_Adj_Old = NULL;

}

CVariable::CVariable(unsigned short val_nvar, CConfig *config) {

  /*--- Array initialization ---*/
  Solution = NULL;
  Solution_Old = NULL;
  Solution_time_n = NULL;
  Solution_time_n1 = NULL;
  Gradient = NULL;
  Rmatrix = NULL;
  Limiter = NULL;
  Solution_Max = NULL;
  Solution_Min = NULL;
  Grad_AuxVar = NULL;
  Undivided_Laplacian = NULL;
  Res_TruncError = NULL;
  Residual_Old = NULL;
  Residual_Sum = NULL;
  Solution_Adj_Old = NULL;

  /*--- Initialize the number of solution variables. This version
   of the constructor will be used primarily for converting the
   restart files into solution files (SU2_SOL). ---*/
  nVar = val_nvar;

  /*--- Allocate the solution array - here it is also possible
   to allocate some extra flow variables that do not participate
   in the simulation ---*/
  Solution = new su2double [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution[iVar] = 0.0;

}

CVariable::CVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config) {

  unsigned short iVar, iDim, jDim;

  /*--- Array initialization ---*/
  Solution = NULL;
  Solution_Old = NULL;
  Solution_time_n = NULL;
  Solution_time_n1 = NULL;
  Gradient = NULL;
  Rmatrix = NULL;
  Limiter = NULL;
  Solution_Max = NULL;
  Solution_Min = NULL;
  Grad_AuxVar = NULL;
  Undivided_Laplacian = NULL;
  Res_TruncError = NULL;
  Residual_Old = NULL;
  Residual_Sum = NULL;
  Solution_Adj_Old = NULL;

  /*--- Initializate the number of dimension and number of variables ---*/
  nDim = val_nDim;
  nVar = val_nvar;

  /*--- Allocate solution, solution old, residual and gradient
   which is common for all the problems, here it is also possible
   to allocate some extra flow variables that do not participate
   in the simulation ---*/
  Solution = new su2double [nVar];

  for (iVar = 0; iVar < nVar; iVar++)
    Solution[iVar] = 0.0;

  Solution_Old = new su2double [nVar];

  Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Gradient[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim ++)
      Gradient[iVar][iDim] = 0.0;
  }

  if (config->GetUnsteady_Simulation() != NO) {
    Solution_time_n = new su2double [nVar];
    Solution_time_n1 = new su2double [nVar];
  }
  else if (config->GetDynamic_Analysis() == DYNAMIC) {
    Solution_time_n = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) Solution_time_n[iVar] = 0.0;
  }

	if (config->GetFSI_Simulation() && config->GetDiscrete_Adjoint()){
	  Solution_Adj_Old = new su2double [nVar];
	}

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    Rmatrix = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Rmatrix[iDim] = new su2double[nDim];
      for (jDim = 0; jDim < nDim; jDim++)
        Rmatrix[iDim][jDim] = 0.0;
    }
  }

}

CVariable::~CVariable(void) {
  unsigned short iVar, iDim;

  if (Solution            != NULL) delete [] Solution;
  if (Solution_Old        != NULL) delete [] Solution_Old;
  if (Solution_time_n     != NULL) delete [] Solution_time_n;
  if (Solution_time_n1    != NULL) delete [] Solution_time_n1;
  if (Limiter             != NULL) delete [] Limiter;
  if (Solution_Max        != NULL) delete [] Solution_Max;
  if (Solution_Min        != NULL) delete [] Solution_Min;
  if (Grad_AuxVar         != NULL) delete [] Grad_AuxVar;
  if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
  if (Res_TruncError      != NULL) delete [] Res_TruncError;
  if (Residual_Old        != NULL) delete [] Residual_Old;
  if (Residual_Sum        != NULL) delete [] Residual_Sum;
  if (Solution_Adj_Old    != NULL) delete [] Solution_Adj_Old;

  if (Gradient != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Gradient[iVar];
    delete [] Gradient;
  }

  if (Rmatrix != NULL) {
    for (iDim = 0; iDim < nDim; iDim++)
      delete [] Rmatrix[iDim];
    delete [] Rmatrix;
  }

}
