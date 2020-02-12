/*!
 * \file CQuadrilateral.cpp
 * \brief Main classes for defining the primal grid elements
 * \author F. Palacios
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

#include "../../../include/geometry/primal_grid/CQuadrilateral.hpp"

unsigned short CQuadrilateral::Faces[4][2] = {{0,1},{1,2},{2,3},{3,0}};

unsigned short CQuadrilateral::Neighbor_Nodes[4][2] = {{1,3},{2,0},{3,1},{0,2}};

unsigned short CQuadrilateral::nNodesFace[4] = {2,2,2,2};

unsigned short CQuadrilateral::nNeighbor_Nodes[4] = {2,2,2,2};

unsigned short CQuadrilateral::nFaces = 4;

unsigned short CQuadrilateral::nNodes = 4;

unsigned short CQuadrilateral::nNeighbor_Elements = 4;

unsigned short CQuadrilateral::VTK_Type = 9;

unsigned short CQuadrilateral::maxNodesFace = 2;

CQuadrilateral::CQuadrilateral(unsigned long val_point_0, unsigned long val_point_1,
             unsigned long val_point_2, unsigned long val_point_3, unsigned short val_nDim)
: CPrimalGrid() {
  unsigned short iDim, iFace, iNeighbor_Elements;

  /*--- Allocate CG coordinates ---*/
  nDim = val_nDim;
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
  Nodes[0] = val_point_0;
  Nodes[1] = val_point_1;
  Nodes[2] = val_point_2;
  Nodes[3] = val_point_3;


  nNeighbor_Elements = nFaces;
  Neighbor_Elements = new long[nNeighbor_Elements];
  for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
    Neighbor_Elements[iNeighbor_Elements]=-1;
  }

}

CQuadrilateral::~CQuadrilateral() {
  unsigned short iFaces;

  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;

}

void CQuadrilateral::Change_Orientation(void) {
  unsigned long Point_1, Point_3;

  Point_1 = Nodes[1];
  Point_3 = Nodes[3];
  Nodes[1] = Point_3;
  Nodes[3] = Point_1;

}
