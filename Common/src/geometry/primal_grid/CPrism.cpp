/*!
 * \file CPrism.cpp
 * \brief Main classes for defining the primal grid elements
 * \author F. Palacios
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/geometry/primal_grid/CPrism.hpp"
#include "../../../include/option_structure.hpp"

unsigned short CPrism::Faces[5][4] = {{3,4,1,0},{5,2,1,4},{2,5,3,0},{0,1,2,2},{5,4,3,3}};

unsigned short CPrism::Neighbor_Nodes[6][3] = {{1,2,3},{0,2,4},{1,0,5},{0,4,5},{3,5,1},{4,3,2}};

unsigned short CPrism::nNodesFace[5] = {4,4,4,3,3};

unsigned short CPrism::nNeighbor_Nodes[6] = {3,3,3,3,3,3};

unsigned short CPrism::nFaces = N_FACES_PRISM;

unsigned short CPrism::nNodes = N_POINTS_PRISM;

unsigned short CPrism::nNeighbor_Elements = 5;

unsigned short CPrism::VTK_Type = 13;

unsigned short CPrism::maxNodesFace = 4;

CPrism::CPrism(unsigned long val_point_0, unsigned long val_point_1,
         unsigned long val_point_2, unsigned long val_point_3,
         unsigned long val_point_4, unsigned long val_point_5) : CPrimalGrid() {
  unsigned short iNeighbor_Elements;

  nDim = 3;

  /*--- Allocate and define face structure of the element ---*/
  Nodes = new unsigned long[nNodes];
  Nodes[0] = val_point_0;
  Nodes[1] = val_point_1;
  Nodes[2] = val_point_2;
  Nodes[3] = val_point_3;
  Nodes[4] = val_point_4;
  Nodes[5] = val_point_5;

  /*--- Allocate and define neighbor elements to a element ---*/
  nNeighbor_Elements = nFaces;
  Neighbor_Elements = new long[nNeighbor_Elements];
  for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
    Neighbor_Elements[iNeighbor_Elements]=-1;
  }

}

CPrism::~CPrism() {}

void CPrism::Change_Orientation(void) {
  swap(Nodes[0], Nodes[1]);
  swap(Nodes[3], Nodes[4]);
}
