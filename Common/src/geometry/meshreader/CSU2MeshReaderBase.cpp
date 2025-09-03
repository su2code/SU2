/*!
 * \file CSU2MeshReaderBase.cpp
 * \brief Helper class for the reading of a native SU2 grid file.
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CSU2MeshReaderBase.hpp"

CSU2MeshReaderBase::CSU2MeshReaderBase(CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone)
    : CMeshReaderBase(val_config, val_iZone, val_nZone),
      myZone(val_iZone),
      nZones(val_nZone),
      meshFilename(config->GetMesh_FileName()) {}

CSU2MeshReaderBase::~CSU2MeshReaderBase(void) = default;
