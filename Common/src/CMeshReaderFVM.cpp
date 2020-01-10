/*!
 * \file CMeshReaderFVM.cpp
 * \brief Helper class that provides the counts for each rank in a linear
 *        partitioning given the global count as input.
 * \author T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
