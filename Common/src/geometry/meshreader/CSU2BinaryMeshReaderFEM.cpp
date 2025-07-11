/*!
 * \file CSU2BinaryMeshReaderFEM.cpp
 * \brief Reads a native SU2 binary grid into linear partitions for the
 *        finite element solver (FEM).
 * \author T. Economon, E. van der Weide
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CSU2BinaryMeshReaderFEM.hpp"
#include "../../../include/fem/fem_standard_element.hpp"

CSU2BinaryMeshReaderFEM::CSU2BinaryMeshReaderFEM(CConfig* val_config, unsigned short val_iZone,
                                                 unsigned short val_nZone)
    : CSU2BinaryMeshReaderBase(val_config, val_iZone, val_nZone) {
  /* Read the basic metadata and perform some basic error checks. */
  ReadMetadata(val_config);

  /*--- Open the file with the mesh and go to the place where the data
        of the current zone is stored. ---*/
  mesh_file = fopen(meshFilename.c_str(), "rb");
  FastForwardToMyZone();

  /*--- Read the volume connectivity and distribute it
        linearly over the MPI ranks. ---*/
  ReadVolumeElementConnectivity();

  /*--- Read the coordinates of the points that are needed
        on this MPI rank. ---*/
  ReadPointCoordinates();

  /*--- Read the surface connectivity and store the surface elements whose
        corresponding volume element is stored on this MPI rank. ---*/
  ReadSurfaceElementConnectivity();

  fclose(mesh_file);
}

CSU2BinaryMeshReaderFEM::~CSU2BinaryMeshReaderFEM() = default;

void CSU2BinaryMeshReaderFEM::ReadPointCoordinates() { SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION); }

void CSU2BinaryMeshReaderFEM::ReadVolumeElementConnectivity() {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

void CSU2BinaryMeshReaderFEM::ReadSurfaceElementConnectivity() {
  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}
