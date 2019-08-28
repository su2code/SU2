/*!
 * \file CDiscAdjMeshVariable.cpp
 * \brief Main subroutines for the discrete adjoint mesh variable structure.
 * \author Ruben Sanchez
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


#include "../../include/variables/CDiscAdjMeshBoundVariable.hpp"

CDiscAdjMeshBoundVariable::CDiscAdjMeshBoundVariable(su2double *val_coor, unsigned short val_nDim, CConfig *config) : CDiscAdjMeshVariable(val_coor, val_nDim, config) {

  unsigned short iDim;

  bool fsi = false;

  /*--- Initialize Boundary Displacement container to 0.0 ---*/
  Bound_Disp_Sens   = new su2double [nDim];
  Bound_Disp_Direct = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++){
    Bound_Disp_Sens[iDim]   = 0.0;
    Bound_Disp_Direct[iDim] = 0.0;
  }

  /*--- Container for the BGS solution at the previous iteration ---*/
  Solution_BGS_k        = NULL;
  if (fsi){
    Solution_BGS_k        = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Solution_BGS_k[iDim]        = 0.0;
    }
  }

}

CDiscAdjMeshBoundVariable::~CDiscAdjMeshBoundVariable(void) {

  if (Bound_Disp_Sens != NULL)   delete [] Bound_Disp_Sens;
  if (Bound_Disp_Direct != NULL) delete [] Bound_Disp_Direct;

}
