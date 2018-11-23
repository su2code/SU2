/*!
 * \file variable_direct_mesh.cpp
 * \brief Definition of the variables for mesh motion using a pseudo-elastic approach.
 * \author R. Sanchez
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

CMeshVariable::CMeshVariable(su2double *val_coor, unsigned short val_nDim, CConfig *config) {

  unsigned short iDim;

  /*--- Initialize pointers to NULL ---*/
  Ref_Coord    = NULL;
  Curr_Coord   = NULL;
  Displacement = NULL;

  Displacement_Old = NULL;

  Displacement_n  = NULL;
  Displacement_n1 = NULL;

  Velocity    = NULL;
  Velocity_n  = NULL;
  Velocity_n1 = NULL;

  /*--- Booleans that determine the kind of problems ---*/
  time_domain = config->GetTime_Domain();
  multizone = config->GetMultizone_Problem();

  /*--- Store the dimensionality of the problem ---*/
  nDim = val_nDim;

  /*--- Initalize the variables that will always be there in a problem with moving mesh ---*/
  Ref_Coord    = new su2double [nDim];
  Curr_Coord   = new su2double [nDim];
  Displacement = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++){
    Ref_Coord[iDim]    = val_coor[iDim];
    Curr_Coord[iDim]   = val_coor[iDim];
    Displacement[iDim] = 0.0;
  }

  /*--- Initialize the variables necessary when the problem is multizone ---*/
  if (multizone){
    Displacement_Old    = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++){
      Displacement_Old[iDim]    = 0.0;
    }
  }

  /*--- Initialize the variables necessary when the problem is time domain ---*/
  if (time_domain){
    Displacement_n    = new su2double [nDim];
    Displacement_n1   = new su2double [nDim];

    Velocity    = new su2double [nDim];
    Velocity_n  = new su2double [nDim];
    Velocity_n1 = new su2double [nDim];

    for (iDim = 0; iDim < nDim; iDim++){
      Displacement_n[iDim]    = 0.0;
      Displacement_n1[iDim]   = 0.0;

      Velocity[iDim]    = 0.0;
      Velocity_n[iDim]  = 0.0;
      Velocity_n1[iDim] = 0.0;
    }
  }

}

CMeshVariable::~CMeshVariable(void) {

  if (Ref_Coord    != NULL) delete [] Ref_Coord;
  if (Curr_Coord   != NULL) delete [] Curr_Coord;
  if (Displacement != NULL) delete [] Displacement;

  if (Displacement_Old != NULL) delete [] Displacement_Old;

  if (Displacement_n  != NULL) delete [] Displacement_n;
  if (Displacement_n1 != NULL) delete [] Displacement_n1;

  if (Velocity    != NULL) delete [] Velocity;
  if (Velocity_n  != NULL) delete [] Velocity_n;
  if (Velocity_n1 != NULL) delete [] Velocity_n1;

}
