/*!
 * \file variable_structure.cpp
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

#include "../include/variable_structure.hpp"

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
  
  unsigned short iVar, iDim;
  
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
  
}

CVariable::~CVariable(void) {
  unsigned short iVar;

  if (Solution            != NULL) delete [] Solution;
  if (Solution_Old        != NULL) delete [] Solution_Old;
  if (Solution_time_n     != NULL) delete [] Solution_time_n;
  if (Solution_time_n1    != NULL) delete [] Solution_time_n1;
  if (Limiter             != NULL) delete [] Limiter;
  if (Solution_Max        != NULL) delete [] Solution_Max;
  if (Solution_Min        != NULL) delete [] Solution_Min;
  if (Grad_AuxVar         != NULL) delete [] Grad_AuxVar;
  //if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian; // Need to break pointer dependence btwn CNumerics and CVariable
  if (Res_TruncError      != NULL) delete [] Res_TruncError;
  if (Residual_Old        != NULL) delete [] Residual_Old;
  if (Residual_Sum        != NULL) delete [] Residual_Sum;
  if (Solution_Adj_Old    != NULL) delete [] Solution_Adj_Old;
  
  if (Gradient != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Gradient[iVar];
    delete [] Gradient;
  }

}

void CVariable::AddUnd_Lapl(su2double *val_und_lapl) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Undivided_Laplacian[iVar] += val_und_lapl[iVar];
}

void CVariable::SubtractUnd_Lapl(su2double *val_und_lapl) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Undivided_Laplacian[iVar] -= val_und_lapl[iVar];
}

void CVariable::SubtractUnd_Lapl(unsigned short val_var, su2double val_und_lapl) {
  Undivided_Laplacian[val_var] -= val_und_lapl;
}

void CVariable::SetUnd_LaplZero(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Undivided_Laplacian[iVar] = 0.0;
  
}

void CVariable::SetUnd_Lapl(unsigned short val_var, su2double val_und_lapl) {
  
    Undivided_Laplacian[val_var] = val_und_lapl;
  
}

void CVariable::SetSolution(su2double *val_solution) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution[iVar] = val_solution[iVar];
  
}

void CVariable::Set_OldSolution(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_Old[iVar] = Solution[iVar];
  
}

void CVariable::Set_OldSolution_Adj(void) {

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_Adj_Old[iVar] = Solution[iVar];

}


void CVariable::AddSolution(unsigned short val_var, su2double val_solution) {
  
  Solution[val_var] = Solution_Old[val_var] + val_solution;
  
}

void CVariable::AddClippedSolution(unsigned short val_var, su2double val_solution,
                                   su2double lowerlimit, su2double upperlimit) {
  
  Solution[val_var] = min(max((Solution_Old[val_var] + val_solution), lowerlimit), upperlimit);
  
}

void CVariable::AddConservativeSolution(unsigned short val_var, su2double val_solution,
    su2double val_density, su2double val_density_old, su2double lowerlimit, su2double upperlimit) {
  
  Solution[val_var] = min(max((Solution_Old[val_var]*val_density_old + val_solution)/val_density,
      lowerlimit), upperlimit);
  
}

void CVariable::Set_Solution(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
     Solution[iVar] = Solution_Old[iVar];
  
}

void CVariable::Set_Solution_time_n(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_time_n[iVar] = Solution[iVar];
  
}

void CVariable::Set_Solution_time_n1(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_time_n1[iVar] = Solution_time_n[iVar];
  
}

void CVariable::Set_Solution_time_n(su2double *val_sol) {

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_time_n[iVar] = val_sol[iVar];

}

void CVariable::Set_Solution_time_n1(su2double *val_sol) {

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_time_n1[iVar] = val_sol[iVar];

}

void CVariable::AddRes_TruncError(su2double *val_truncation_error) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Res_TruncError[iVar] += val_truncation_error[iVar];
  
}

void CVariable::SubtractRes_TruncError(su2double *val_truncation_error) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Res_TruncError[iVar] -= val_truncation_error[iVar];
  
}

void CVariable::SetResidual_Old(su2double *val_residual_old) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Residual_Old[iVar] = val_residual_old[iVar];
  
}

void CVariable::SetSolution_Old(su2double *val_solution_old) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_Old[iVar] = val_solution_old[iVar];
  
}

void CVariable::SetSolution_time_n(su2double *val_solution_time_n) {

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution_time_n[iVar] = val_solution_time_n[iVar];

}

void CVariable::AddResidual_Sum(su2double *val_residual) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Residual_Sum[iVar] += val_residual[iVar];
  
}

void CVariable::SetVel_ResTruncError_Zero(void) {
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Res_TruncError[iDim+1] = 0.0;
  
}

void CVariable::SetEnergy_ResTruncError_Zero(void) {
  
  Res_TruncError[nDim+1] = 0.0;
  
}

void CVariable::SetVelSolutionZero(void) {
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution[iDim+1] = 0.0;
  
}

void CVariable::SetVelSolutionVector(su2double *val_vector) {
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution[iDim+1] = val_vector[iDim];
  
}

void CVariable::SetVelSolutionOldZero(void) {
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_Old[iDim+1] = 0.0;
  
}

void CVariable::SetVelSolutionOldVector(su2double *val_vector) {
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_Old[iDim+1] = val_vector[iDim];
  
}

void CVariable::SetSolutionZero(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution[iVar] = 0.0;
  
}

void CVariable::SetSolutionZero(unsigned short val_var) {
  
    Solution[val_var] = 0.0;
  
}

void CVariable::SetResidualSumZero(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Residual_Sum[iVar] = 0.0;
  
}

void CVariable::SetGradientZero(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Gradient[iVar][iDim] = 0.0;
  
}

void CVariable::SetAuxVarGradientZero(void) {
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Grad_AuxVar[iDim] = 0.0;
  
}

void CVariable::SetGradient(su2double **val_gradient) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Gradient[iVar][iDim] = val_gradient[iVar][iDim];
  
}

void CVariable::SetRes_TruncErrorZero(void) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Res_TruncError[iVar] = 0.0;
  
}

void CVariable::SetVal_ResTruncError_Zero(unsigned short val_var) {
  
    Res_TruncError[val_var] = 0.0;
  
}

void CVariable::GetResidual_Sum(su2double *val_residual) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Residual_Sum[iVar];
  
}

void CVariable::GetResTruncError(su2double *val_trunc_error) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    val_trunc_error[iVar] = Res_TruncError[iVar];
  
}

CBaselineVariable::CBaselineVariable(void) : CVariable() { }

CBaselineVariable::CBaselineVariable(su2double *val_solution, unsigned short val_nvar, CConfig *config) : CVariable(val_nvar, config) {
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Solution[iVar] = val_solution[iVar];
  
}

CBaselineVariable::~CBaselineVariable(void) { }
