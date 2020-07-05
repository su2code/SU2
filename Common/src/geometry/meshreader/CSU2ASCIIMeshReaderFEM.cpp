/*!
 * \file CSU2ASCIIMeshReaderFVM.cpp
 * \brief Reads a native SU2 ASCII grid into linear partitions for the
 *        finite element solver (FEM).
 * \author T. Economon, E. van der Weide
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/geometry/meshreader/CSU2ASCIIMeshReaderFEM.hpp"

CSU2ASCIIMeshReaderFEM::CSU2ASCIIMeshReaderFEM(CConfig        *val_config,
                                               unsigned short val_iZone,
                                               unsigned short val_nZone)
: CSU2ASCIIMeshReaderBase(val_config, val_iZone, val_nZone) {
  
  /* Read the basic metadata and perform some basic error checks. */
  ReadMetadata();
  
  /* Read and store the interior elements, points, and surface elements. */
  ReadVolumeElementConnectivity();
  ReadPointCoordinates();
  ReadSurfaceElementConnectivity();
  
}

CSU2ASCIIMeshReaderFEM::~CSU2ASCIIMeshReaderFEM(void) { }

void CSU2ASCIIMeshReaderFEM::ReadPointCoordinates() {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CSU2ASCIIMeshReaderFEM::ReadVolumeElementConnectivity() {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CSU2ASCIIMeshReaderFEM::ReadSurfaceElementConnectivity() {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
