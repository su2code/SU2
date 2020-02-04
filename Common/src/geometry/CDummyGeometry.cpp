/*!
 * \file CDummyGeometry.hpp
 * \brief Implementation of the dummy geometry class used in "dry run" mode.
 * \author T. Albring
 * \version 7.0.1 "Blackbird"
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

  nElem_Bound         = NULL;
  Tag_to_Marker       = NULL;
  elem                = NULL;
  face                = NULL;
  bound               = NULL;
  node                = NULL;
  edge                = NULL;
  vertex              = NULL;
  nVertex             = NULL;
  newBound            = NULL;
  nNewElem_Bound      = NULL;
  Marker_All_SendRecv = NULL;

  XCoordList.clear();
  Xcoord_plane.clear();
  Ycoord_plane.clear();
  Zcoord_plane.clear();
  FaceArea_plane.clear();
  Plane_points.clear();

  /*--- Arrays for defining the linear partitioning ---*/

  beg_node = NULL;
  end_node = NULL;

  nPointLinear = NULL;
  nPointCumulative = NULL;

  /*--- Containers for customized boundary conditions ---*/

  CustomBoundaryHeatFlux = NULL;      //Customized heat flux wall
  CustomBoundaryTemperature = NULL;   //Customized temperature wall

  /*--- MPI point-to-point data structures ---*/

  nP2PSend = 0;
  nP2PRecv = 0;

  countPerPoint = 0;

  bufD_P2PSend = NULL;
  bufD_P2PRecv = NULL;

  bufS_P2PSend = NULL;
  bufS_P2PRecv = NULL;

  req_P2PSend = NULL;
  req_P2PRecv = NULL;

  nPoint_P2PSend = new int[size];
  nPoint_P2PRecv = new int[size];

  Neighbors_P2PSend = NULL;
  Neighbors_P2PRecv = NULL;

  Local_Point_P2PSend = NULL;
  Local_Point_P2PRecv = NULL;

  /*--- MPI periodic data structures ---*/

  nPeriodicSend = 0;
  nPeriodicRecv = 0;

  countPerPeriodicPoint = 0;

  bufD_PeriodicSend = NULL;
  bufD_PeriodicRecv = NULL;

  bufS_PeriodicSend = NULL;
  bufS_PeriodicRecv = NULL;

  req_PeriodicSend = NULL;
  req_PeriodicRecv = NULL;

  nPoint_PeriodicSend = NULL;
  nPoint_PeriodicRecv = NULL;

  Neighbors_PeriodicSend = NULL;
  Neighbors_PeriodicRecv = NULL;

  Local_Point_PeriodicSend = NULL;
  Local_Point_PeriodicRecv = NULL;

  Local_Marker_PeriodicSend = NULL;
  Local_Marker_PeriodicRecv = NULL;

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
