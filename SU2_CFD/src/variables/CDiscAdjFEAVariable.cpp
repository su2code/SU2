/*!
 * \file CDiscAdjFEAVariable.cpp
 * \brief Definition of the variables for FEM adjoint elastic structural problems.
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

#include "../../include/variables/CDiscAdjFEAVariable.hpp"

CDiscAdjFEAVariable::CDiscAdjFEAVariable() : CVariable(){

  Sensitivity           = NULL;
  Solution_Direct       = NULL;

  Dynamic_Derivative          = NULL;
  Dynamic_Derivative_n        = NULL;
  Dynamic_Derivative_Vel      = NULL;
  Dynamic_Derivative_Vel_n    = NULL;
  Dynamic_Derivative_Accel    = NULL;
  Dynamic_Derivative_Accel_n  = NULL;

  Solution_Direct_Vel   = NULL;
  Solution_Direct_Accel = NULL;

  Solution_Vel   = NULL;
  Solution_Accel = NULL;

  Solution_Old_Vel      = NULL;
  Solution_Old_Accel    = NULL;

  Solution_Vel_time_n   = NULL;
  Solution_Accel_time_n = NULL;

  Cross_Term_Derivative = NULL;
  Geometry_CrossTerm_Derivative = NULL;

  Solution_BGS          = NULL;
  Solution_BGS_k        = NULL;

}

CDiscAdjFEAVariable::CDiscAdjFEAVariable(su2double* val_solution, unsigned short val_ndim, unsigned short val_nvar,
                                         CConfig *config) : CVariable(val_ndim, val_nvar, config){

  bool fsi = config->GetFSI_Simulation();

  Dynamic_Derivative          = NULL;
  Dynamic_Derivative_n        = NULL;
  Dynamic_Derivative_Vel      = NULL;
  Dynamic_Derivative_Vel_n    = NULL;
  Dynamic_Derivative_Accel    = NULL;
  Dynamic_Derivative_Accel_n  = NULL;

  Solution_Direct_Vel   = NULL;
  Solution_Direct_Accel = NULL;

  Solution_Vel          = NULL;
  Solution_Accel        = NULL;

  Solution_Old_Vel      = NULL;
  Solution_Old_Accel    = NULL;

  Solution_Vel_time_n   = NULL;
  Solution_Accel_time_n = NULL;

  Solution_Direct = new su2double[nVar];

  Sensitivity = new su2double[nDim];

  unsigned short iVar,iDim;

  for (iDim = 0; iDim < nDim; iDim++){
    Sensitivity[iDim] = 0.0;
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution[iVar] = val_solution[iVar];
  }

  Solution_BGS          = NULL;
  Solution_BGS_k        = NULL;
  Cross_Term_Derivative = NULL;
  Geometry_CrossTerm_Derivative = NULL;
  if (fsi){
    Cross_Term_Derivative = new su2double[nDim];
    Geometry_CrossTerm_Derivative = new su2double[nDim];
    Solution_BGS          = new su2double[nDim];
    Solution_BGS_k        = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Geometry_CrossTerm_Derivative [iDim] = 0.0;
      Cross_Term_Derivative[iDim] = 0.0;
      Solution_BGS[iDim]          = 0.0;
      Solution_BGS_k[iDim]        = 0.0;
    }
  }

}

CDiscAdjFEAVariable::CDiscAdjFEAVariable(su2double* val_solution, su2double* val_solution_accel, su2double* val_solution_vel,
                                         unsigned short val_ndim, unsigned short val_nvar, CConfig *config) :
                                         CVariable(val_ndim, val_nvar, config){

  bool fsi = config->GetFSI_Simulation();

  Dynamic_Derivative          = new su2double[nVar];
  Dynamic_Derivative_n        = new su2double[nVar];
  Dynamic_Derivative_Vel      = new su2double[nVar];
  Dynamic_Derivative_Vel_n    = new su2double[nVar];
  Dynamic_Derivative_Accel    = new su2double[nVar];
  Dynamic_Derivative_Accel_n  = new su2double[nVar];

  Solution_Direct_Vel         = new su2double[nVar];
  Solution_Direct_Accel       = new su2double[nVar];

  Solution_Vel                = new su2double[nVar];
  Solution_Accel              = new su2double[nVar];

  Solution_Old_Vel            = new su2double[nVar];
  Solution_Old_Accel          = new su2double[nVar];

  Solution_Vel_time_n         = new su2double[nVar];
  Solution_Accel_time_n       = new su2double[nVar];

  Solution_Direct = new su2double[nVar];

  Sensitivity = new su2double[nDim];

  unsigned short iVar,iDim;

  for (iDim = 0; iDim < nDim; iDim++){
    Sensitivity[iDim] = 0.0;
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution[iVar] = val_solution[iVar];
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution_Accel[iVar] = val_solution_accel[iVar];
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution_Vel[iVar] = val_solution_vel[iVar];
  }

  /*--- Initialize the rest to 0 ---*/

  for (iVar = 0; iVar < nVar; iVar++){
    Dynamic_Derivative[iVar]      = 0.0;
    Dynamic_Derivative_n[iVar]    = 0.0;
    Dynamic_Derivative_Vel[iVar]    = 0.0;
    Dynamic_Derivative_Vel_n[iVar]    = 0.0;
    Dynamic_Derivative_Accel[iVar]    = 0.0;
    Dynamic_Derivative_Accel_n[iVar]  = 0.0;

    Solution_Direct_Vel[iVar]     = 0.0;
    Solution_Direct_Accel[iVar]   = 0.0;

    Solution_Vel_time_n[iVar]     = 0.0;
    Solution_Accel_time_n[iVar]   = 0.0;

    Solution_Old_Vel[iVar]        = 0.0;
    Solution_Old_Accel[iVar]      = 0.0;

  }

  Solution_BGS          = NULL;
  Solution_BGS_k        = NULL;
  Cross_Term_Derivative = NULL;
  Geometry_CrossTerm_Derivative = NULL;
  if (fsi){
    Cross_Term_Derivative = new su2double[nDim];
    Geometry_CrossTerm_Derivative = new su2double[nDim];
    Solution_BGS_k        = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Geometry_CrossTerm_Derivative [iDim] = 0.0;
      Cross_Term_Derivative[iDim] = 0.0;
      Solution_BGS[iDim]          = 0.0;
      Solution_BGS_k[iDim]        = 0.0;
    }
  }

}


CDiscAdjFEAVariable::~CDiscAdjFEAVariable(){

  if (Sensitivity           != NULL) delete [] Sensitivity;
  if (Solution_Direct       != NULL) delete [] Solution_Direct;

  if (Dynamic_Derivative         != NULL) delete [] Dynamic_Derivative;
  if (Dynamic_Derivative_n       != NULL) delete [] Dynamic_Derivative_n;
  if (Dynamic_Derivative_Vel     != NULL) delete [] Dynamic_Derivative_Vel;
  if (Dynamic_Derivative_Vel_n   != NULL) delete [] Dynamic_Derivative_Vel_n;
  if (Dynamic_Derivative_Accel   != NULL) delete [] Dynamic_Derivative_Accel;
  if (Dynamic_Derivative_Accel_n != NULL) delete [] Dynamic_Derivative_Accel_n;

  if (Solution_Direct_Vel   != NULL) delete [] Solution_Direct_Vel;
  if (Solution_Direct_Accel != NULL) delete [] Solution_Direct_Accel;

  if (Solution_Vel          != NULL) delete [] Solution_Vel;
  if (Solution_Accel        != NULL) delete [] Solution_Accel;

  if (Solution_Vel_time_n   != NULL) delete [] Solution_Vel_time_n;
  if (Solution_Accel_time_n != NULL) delete [] Solution_Accel_time_n;

  if (Solution_Old_Vel      != NULL) delete [] Solution_Old_Vel;
  if (Solution_Old_Accel    != NULL) delete [] Solution_Old_Accel;

  if (Cross_Term_Derivative    != NULL) delete [] Cross_Term_Derivative;
  if (Geometry_CrossTerm_Derivative    != NULL) delete [] Geometry_CrossTerm_Derivative;

  if (Solution_BGS             != NULL) delete [] Solution_BGS;
  if (Solution_BGS_k           != NULL) delete [] Solution_BGS_k;

}
