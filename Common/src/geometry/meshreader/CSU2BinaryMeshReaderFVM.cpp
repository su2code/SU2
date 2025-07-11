/*!
 * \file CSU2BinaryMeshReaderFVM.cpp
 * \brief Reads a native SU2 binary grid into linear partitions for the
 *        finite volume solver (FVM).
 * \author T. Economon
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

#include "../../../include/geometry/meshreader/CSU2BinaryMeshReaderFVM.hpp"

CSU2BinaryMeshReaderFVM::CSU2BinaryMeshReaderFVM(CConfig* val_config, unsigned short val_iZone,
                                                 unsigned short val_nZone)
    : CSU2BinaryMeshReaderBase(val_config, val_iZone, val_nZone) {
  actuator_disk = (((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0)) &&
                   ((config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD) ||
                    ((config->GetKind_SU2() == SU2_COMPONENT::SU2_DEF) && (config->GetActDisk_SU2_DEF()))));
  if (config->GetActDisk_DoubleSurface()) actuator_disk = false;

  /* Read the basic metadata and perform some basic error checks. */
  ReadMetadata(val_config);

  /* If the mesh contains an actuator disk as a single surface,
   we need to first split the surface into repeated points and update
   the connectivity for each element touching the surface. */
  if (actuator_disk) SplitActuatorDiskSurface();

  /* Read and store the points, interior elements, and surface elements.
   We store only the points and interior elements on our rank's linear
   partition, but the master stores the entire set of surface connectivity. */
  mesh_file = fopen(meshFilename.c_str(), "rb");

  FastForwardToMyZone();
  ReadVolumeElementConnectivity();
  ReadPointCoordinates();
  ReadSurfaceElementConnectivity();

  fclose(mesh_file);
}

CSU2BinaryMeshReaderFVM::~CSU2BinaryMeshReaderFVM() = default;

void CSU2BinaryMeshReaderFVM::SplitActuatorDiskSurface() {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
