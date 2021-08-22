/*!
 * \file CTetrahedron.cpp
 * \brief Main classes for defining the primal grid elements
 * \author F. Palacios
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/geometry/primal_grid/CTetrahedron.hpp"
#include "../../../include/option_structure.hpp"

constexpr unsigned short CTetrahedronConnectivity::Faces[4][3];
constexpr unsigned short CTetrahedronConnectivity::Neighbor_Nodes[4][3];
constexpr unsigned short CTetrahedronConnectivity::nNodesFace[4];
constexpr unsigned short CTetrahedronConnectivity::nNeighbor_Nodes[4];
constexpr unsigned short CTetrahedronConnectivity::nFaces;
constexpr unsigned short CTetrahedronConnectivity::nNodes;
constexpr unsigned short CTetrahedronConnectivity::nNeighbor_Elements;
constexpr unsigned short CTetrahedronConnectivity::VTK_Type;
constexpr unsigned short CTetrahedronConnectivity::maxNodesFace;

CTetrahedron::CTetrahedron(unsigned long val_point_0, unsigned long val_point_1,
               unsigned long val_point_2, unsigned long val_point_3) {
  unsigned short iNeighbor_Elements;

  /*--- Allocate and define face structure of the element ---*/
  Nodes = new unsigned long[GetnNodes()];
  Nodes[0] = val_point_0;
  Nodes[1] = val_point_1;
  Nodes[2] = val_point_2;
  Nodes[3] = val_point_3;

  /*--- Allocate and define neighbor elements to a element ---*/
  Neighbor_Elements = new long[GetnNeighbor_Elements()];
  for (iNeighbor_Elements = 0; iNeighbor_Elements<GetnNeighbor_Elements(); iNeighbor_Elements++) {
    Neighbor_Elements[iNeighbor_Elements]=-1;
  }

}

CTetrahedron::~CTetrahedron() {}

void CTetrahedron::Change_Orientation(void) {
  swap(Nodes[0],Nodes[1]);
}
