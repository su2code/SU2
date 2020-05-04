/*!
 * \file CDummyGeometry.hpp
 * \brief Implementation of the dummy geometry class used in "dry run" mode.
 * \author T. Albring
 * \version 7.0.4 "Blackbird"
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

#include "../../include/geometry/CDummyGeometry.hpp"


CDummyGeometry::CDummyGeometry(CConfig *config){

  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  nEdge      = 0;
  nPoint     = 0;
  nPointDomain = 0;
  nPointNode = 0;
  nElem      = 0;
  nMarker    = 0;
  nZone = config->GetnZone();

  nElem_Bound         = nullptr;
  Tag_to_Marker       = nullptr;
  elem                = nullptr;
  face                = nullptr;
  bound               = nullptr;
  node                = nullptr;
  edges               = nullptr;
  vertex              = nullptr;
  nVertex             = nullptr;
  newBound            = nullptr;
  nNewElem_Bound      = nullptr;
  Marker_All_SendRecv = nullptr;

  XCoordList.clear();
  Xcoord_plane.clear();
  Ycoord_plane.clear();
  Zcoord_plane.clear();
  FaceArea_plane.clear();
  Plane_points.clear();

  /*--- Arrays for defining the linear partitioning ---*/

  beg_node = nullptr;
  end_node = nullptr;

  nPointLinear = nullptr;
  nPointCumulative = nullptr;

  /*--- Containers for customized boundary conditions ---*/

  CustomBoundaryHeatFlux = nullptr;      //Customized heat flux wall
  CustomBoundaryTemperature = nullptr;   //Customized temperature wall

  /*--- MPI point-to-point data structures ---*/

  nP2PSend = 0;
  nP2PRecv = 0;

  countPerPoint = 0;

  bufD_P2PSend = nullptr;
  bufD_P2PRecv = nullptr;

  bufS_P2PSend = nullptr;
  bufS_P2PRecv = nullptr;

  req_P2PSend = nullptr;
  req_P2PRecv = nullptr;

  nPoint_P2PSend = new int[size];
  nPoint_P2PRecv = new int[size];

  Neighbors_P2PSend = nullptr;
  Neighbors_P2PRecv = nullptr;

  Local_Point_P2PSend = nullptr;
  Local_Point_P2PRecv = nullptr;

  /*--- MPI periodic data structures ---*/

  nPeriodicSend = 0;
  nPeriodicRecv = 0;

  countPerPeriodicPoint = 0;

  bufD_PeriodicSend = nullptr;
  bufD_PeriodicRecv = nullptr;

  bufS_PeriodicSend = nullptr;
  bufS_PeriodicRecv = nullptr;

  req_PeriodicSend = nullptr;
  req_PeriodicRecv = nullptr;

  nPoint_PeriodicSend = nullptr;
  nPoint_PeriodicRecv = nullptr;

  Neighbors_PeriodicSend = nullptr;
  Neighbors_PeriodicRecv = nullptr;

  Local_Point_PeriodicSend = nullptr;
  Local_Point_PeriodicRecv = nullptr;

  Local_Marker_PeriodicSend = nullptr;
  Local_Marker_PeriodicRecv = nullptr;

  nVertex = new unsigned long[config->GetnMarker_All()];

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    nVertex[iMarker] = 0;
  }

  Tag_to_Marker = new string[config->GetnMarker_All()];

  for (unsigned short iRank = 0; iRank < size; iRank++){
    nPoint_P2PRecv[iRank] = 0;
    nPoint_P2PSend[iRank] = 0;
  }

  nDim = CConfig::GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());

  config->SetnSpanWiseSections(0);
}

CDummyGeometry::~CDummyGeometry(){}
