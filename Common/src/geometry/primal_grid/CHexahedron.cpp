/*!
 * \file CHexahedron.cpp
 * \brief Main classes for defining the primal grid elements
 * \author F. Palacios
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

#include "../../../include/geometry/primal_grid/CHexahedron.hpp"

unsigned short CHexahedron::Faces[6][4] = {{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{0,3,2,1},{4,5,6,7}};

unsigned short CHexahedron::Neighbor_Nodes[8][3] = {{1,3,4},{0,2,5},{1,3,6},{0,2,7},{0,5,7},{4,6,1},{2,5,7},{4,3,6}};

unsigned short CHexahedron::nNodesFace[6] = {4,4,4,4,4,4};

unsigned short CHexahedron::nNeighbor_Nodes[8] = {3,3,3,3,3,3,3,3};

unsigned short CHexahedron::nFaces = 6;

unsigned short CHexahedron::nNodes = 8;

unsigned short CHexahedron::nNeighbor_Elements = 6;

unsigned short CHexahedron::VTK_Type = 12;

unsigned short CHexahedron::maxNodesFace = 4;

CHexahedron::CHexahedron(unsigned long val_point_0, unsigned long val_point_1,
             unsigned long val_point_2, unsigned long val_point_3,
             unsigned long val_point_4, unsigned long val_point_5,
             unsigned long val_point_6, unsigned long val_point_7) : CPrimalGrid() {
  unsigned short iDim, iFace, iNeighbor_Elements;

  /*--- Allocate center-of-gravity coordinates ---*/
  nDim = 3;
  Coord_CG = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord_CG[iDim] = 0.0;
  Coord_FaceElems_CG = new su2double* [nFaces];
  for (iFace = 0; iFace < nFaces; iFace++) {
    Coord_FaceElems_CG[iFace] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Coord_FaceElems_CG[iFace][iDim] = 0.0;
  }

  /*--- Allocate and define face structure of the element ---*/
  Nodes = new unsigned long[nNodes];
  Nodes[0] = val_point_0; Nodes[1] = val_point_1;
  Nodes[2] = val_point_2; Nodes[3] = val_point_3;
  Nodes[4] = val_point_4; Nodes[5] = val_point_5;
  Nodes[6] = val_point_6; Nodes[7] = val_point_7;

  /*--- Allocate and define neighbor elements to a element ---*/
  nNeighbor_Elements = nFaces;
  Neighbor_Elements = new long[nNeighbor_Elements];
  for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
    Neighbor_Elements[iNeighbor_Elements]=-1;
  }

}

CHexahedron::~CHexahedron() {
  unsigned short iFaces;

  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;

}

void CHexahedron::Change_Orientation(void) {
  unsigned long Point_1, Point_3, Point_5, Point_7;

  Point_1 = Nodes[1];
  Point_3 = Nodes[3];
  Point_5 = Nodes[5];
  Point_7 = Nodes[7];
  Nodes[1] = Point_3;
  Nodes[3] = Point_1;
  Nodes[5] = Point_7;
  Nodes[7] = Point_5;

}

