/*!
 * \file CDummyMesh_DG.cpp
 * \brief Implementations of the member functions of CDummyMesh_DG.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/geometry/fem_grid/CDummyMeshFEM_DG.hpp"

CDummyMeshFEM_DG::CDummyMeshFEM_DG(CConfig *config): CMeshFEM_DG() {

  nZone = config->GetnZone();

  nPoint_P2PSend = new int[size] ();
  nPoint_P2PRecv = new int[size] ();

  nVertex = new unsigned long[config->GetnMarker_All()] ();

  Tag_to_Marker = new string[config->GetnMarker_All()];

  for (unsigned short i=0; i <= config->GetnLevels_TimeAccurateLTS(); i++){
    nMatchingFacesWithHaloElem.push_back(0);
  }

  boundaries.resize(config->GetnMarker_All());

  nDim = CConfig::GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());

}
