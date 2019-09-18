/*!
 * \file CMeshVariable.cpp
 * \brief Definition of the variables for mesh motion using a pseudo-elastic approach.
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

#include "../../include/variables/CMeshVariable.hpp"

CMeshVariable::CMeshVariable(const su2double *val_coor, unsigned short val_nDim, CConfig *config) : CVariable(val_nDim, config) {

  unsigned short iDim;

  /*--- Initialize pointers to NULL ---*/
  Mesh_Coord       = nullptr;

  Solution_Old     = nullptr;

  Solution_time_n  = nullptr;
  Solution_time_n1 = nullptr;

  /*--- Booleans that determine the kind of problems ---*/
  bool time_domain = config->GetTime_Domain();
  bool multizone = config->GetMultizone_Problem();

  /*--- Store the dimensionality of the problem ---*/
  nDim = val_nDim;

  /*--- Initalize the variables that will always be there in a problem with moving mesh ---*/
  Mesh_Coord    = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++){
    Mesh_Coord[iDim]    = val_coor[iDim];
  }

  /*--- Initialize the variables necessary when the problem is multizone ---*/
  if (multizone){
    Solution_Old    = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++){
      Solution_Old[iDim]    = 0.0;
    }
  }

  /*--- Initialize the variables necessary when the problem is time domain ---*/
  if (time_domain){
    Solution_time_n    = new su2double [nDim];
    Solution_time_n1   = new su2double [nDim];

    for (iDim = 0; iDim < nDim; iDim++){
      Solution_time_n[iDim]    = 0.0;
      Solution_time_n1[iDim]   = 0.0;
    }
  }

}

CMeshVariable::~CMeshVariable(void) {

  if (Mesh_Coord != nullptr) delete [] Mesh_Coord;

}
