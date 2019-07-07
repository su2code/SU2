/*!
 * \file CDiscAdjVariable.cpp
 * \brief Main subroutines for the discrete adjoint variable structure.
 * \author T. Albring
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

#include "../../include/variables/CDiscAdjVariable.hpp"

CDiscAdjVariable::CDiscAdjVariable() : CVariable() {

  /*--- Initialize arrays to NULL ---*/

  Solution_Direct = NULL;
  Sensitivity     = NULL;

  DualTime_Derivative   = NULL;
  DualTime_Derivative_n = NULL;

  Geometry_Direct       = NULL;
  Solution_Geometry     = NULL;
  Solution_Geometry_Old = NULL;
  Cross_Term_Derivative = NULL;

  Solution_BGS            = NULL;
  Solution_BGS_k          = NULL;
  Solution_Geometry_BGS_k = NULL;

  Geometry_CrossTerm_Derivative      = NULL;
  Geometry_CrossTerm_Derivative_Flow = NULL;

}

CDiscAdjVariable::CDiscAdjVariable(su2double* val_solution, unsigned short val_ndim, unsigned short val_nvar,
                                   CConfig *config) : CVariable(val_ndim, val_nvar, config) {

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  bool fsi = config->GetFSI_Simulation();

  /*--- Initialize arrays to NULL ---*/

  Solution_Direct = NULL;
  Sensitivity     = NULL;

  DualTime_Derivative   = NULL;
  DualTime_Derivative_n = NULL;

  Geometry_Direct       = NULL;
  Solution_Geometry     = NULL;
  Solution_Geometry_Old = NULL;
  Cross_Term_Derivative = NULL;

  Solution_BGS            = NULL;
  Solution_BGS_k          = NULL;
  Solution_Geometry_BGS_k = NULL;

  Geometry_CrossTerm_Derivative      = NULL;
  Geometry_CrossTerm_Derivative_Flow = NULL;

  if (dual_time) {
    DualTime_Derivative = new su2double[nVar];
    DualTime_Derivative_n = new su2double[nVar];
  }

  Solution_Direct = new su2double[nVar];

  Sensitivity = new su2double[nDim];

  unsigned short iVar,iDim;

  for (iDim = 0; iDim < nDim; iDim++) {
    Sensitivity[iDim] = 0.0;
  }

  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar] = val_solution[iVar];
  }


  if (dual_time) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution_time_n[iVar]  = 0.0;
      Solution_time_n1[iVar] = 0.0;
      DualTime_Derivative[iVar] = 0.0;
      DualTime_Derivative_n[iVar] = 0.0;
    }
  }

  if (fsi){
    Solution_Geometry       = new su2double[nDim];
    Geometry_Direct         = new su2double[nDim];
    Solution_Geometry_Old   = new su2double[nDim];
    Geometry_CrossTerm_Derivative = new su2double[nDim];
    Geometry_CrossTerm_Derivative_Flow = new su2double[nDim];
    Cross_Term_Derivative   = new su2double[nVar];
    Solution_BGS            = new su2double[nVar];
    Solution_BGS_k          = new su2double[nVar];
    Solution_Geometry_BGS_k = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Geometry_Direct[iDim]       = 0.0;
      Solution_Geometry[iDim]     = 1e-16;
      Solution_Geometry_Old[iDim] = 0.0;
      Solution_Geometry_BGS_k[iDim] = 0.0;
      Geometry_CrossTerm_Derivative[iDim] = 0.0;
      Geometry_CrossTerm_Derivative_Flow[iDim] = 0.0;
    }
    for (iVar = 0; iVar < nVar; iVar++) {
      Cross_Term_Derivative[iVar] = 0.0;
      Solution_BGS[iVar]          = 0.0;
      Solution_BGS_k[iVar]        = 0.0;
    }
  }

}

CDiscAdjVariable::~CDiscAdjVariable() {

  if (Geometry_Direct       != NULL) delete [] Geometry_Direct;
  if (Solution_Geometry     != NULL) delete [] Solution_Geometry;
  if (Solution_Geometry_Old != NULL) delete [] Solution_Geometry_Old;
  if (Cross_Term_Derivative != NULL) delete [] Cross_Term_Derivative;
  if (Geometry_CrossTerm_Derivative != NULL) delete [] Geometry_CrossTerm_Derivative;
  if (Geometry_CrossTerm_Derivative_Flow != NULL) delete [] Geometry_CrossTerm_Derivative_Flow;
  if (Solution_BGS          != NULL) delete [] Solution_BGS;
  if (Solution_BGS_k        != NULL) delete [] Solution_BGS_k;
  if (Solution_Geometry_BGS_k != NULL) delete [] Solution_Geometry_BGS_k;

  if (Solution_Direct != NULL) delete [] Solution_Direct;
  if (Sensitivity     != NULL) delete [] Sensitivity;

  if (DualTime_Derivative   != NULL) delete [] DualTime_Derivative;
  if (DualTime_Derivative_n != NULL) delete [] DualTime_Derivative_n;

}
