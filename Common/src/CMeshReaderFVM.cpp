/*!
 * \file CMeshReaderFVM.cpp
 * \brief Helper class that provides the counts for each rank in a linear
 *        partitioning given the global count as input.
 * \author T. Economon
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

#include "../include/CMeshReaderFVM.hpp"

CMeshReaderFVM::CMeshReaderFVM(CConfig        *val_config,
                               unsigned short val_iZone,
                               unsigned short val_nZone) {
  
  /*--- Store MPI size ---*/
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  this->config = val_config;

  dimension = 0;
  
  numberOfLocalPoints = 0;
  numberOfGlobalPoints = 0;
  localPointCoordinates.clear();
  
  numberOfLocalElements = 0;
  numberOfGlobalElements = 0;
  localVolumeElementConnectivity.clear();
  
  numberOfMarkers = 0;
  markerNames.clear();
  surfaceElementConnectivity.clear();
  
}

CMeshReaderFVM::~CMeshReaderFVM(void) { }
