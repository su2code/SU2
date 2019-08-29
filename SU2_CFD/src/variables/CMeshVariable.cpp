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

CMeshVariable::CMeshVariable(Idx_t npoint, Idx_t ndim, CConfig *config) :
  CVariable(npoint, ndim, config) {

  /*--- Booleans that determine the kind of problems ---*/
  bool time_domain = config->GetTime_Domain();
  bool multizone = config->GetMultizone_Problem();

  /*--- Store the dimensionality of the problem ---*/
  nDim = ndim;

  /*--- Initalize the variables that will always be there in a problem with moving mesh ---*/
  Mesh_Coord.resize(nPoint,nDim) = su2double(0.0);

  /*--- Initialize the variables necessary when the problem is multizone ---*/
  if (multizone)
    Solution_Old.resize(nPoint,nDim) = su2double(0.0);

  /*--- Initialize the variables necessary when the problem is time domain ---*/
  if (time_domain) {
    Solution_time_n.resize(nPoint,nDim) = su2double(0.0);
    Solution_time_n1.resize(nPoint,nDim) = su2double(0.0);
  }
}
