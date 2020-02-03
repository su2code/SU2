/*!
 * \file CPrimalGrid.cpp
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

#include "../../../include/geometry/primal_grid/CPrimalGrid.hpp"

unsigned short CPrimalGrid::nDim;

CPrimalGrid::CPrimalGrid(void) {

  /*--- Set the default values for the pointers ---*/
  Nodes = NULL;
  Neighbor_Elements = NULL;
  ElementOwnsFace = NULL;
  PeriodIndexNeighbors = NULL;
  Coord_CG = NULL;
  Coord_FaceElems_CG = NULL;
  JacobianFaceIsConstant = NULL;
  GlobalIndex = 0;

}

CPrimalGrid::~CPrimalGrid() {

 if (Nodes != NULL) delete[] Nodes;
 if (Coord_CG != NULL) delete[] Coord_CG;
 if (Neighbor_Elements != NULL) delete[] Neighbor_Elements;
 if (ElementOwnsFace != NULL) delete[] ElementOwnsFace;
 if (PeriodIndexNeighbors != NULL) delete[] PeriodIndexNeighbors;
 if (JacobianFaceIsConstant != NULL) delete[] JacobianFaceIsConstant;
}

void CPrimalGrid::SetCoord_CG(su2double **val_coord) {
  unsigned short iDim, iNode, NodeFace, iFace;

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord, GetnNodes(), nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_CG[iDim] = 0.0;
    for (iNode = 0; iNode < GetnNodes();  iNode++)
      Coord_CG[iDim] += val_coord[iNode][iDim]/su2double(GetnNodes());
  }

  for (iFace = 0; iFace < GetnFaces();  iFace++)
    for (iDim = 0; iDim < nDim; iDim++) {
      Coord_FaceElems_CG[iFace][iDim] = 0.0;
      for (iNode = 0; iNode < GetnNodesFace(iFace); iNode++) {
        NodeFace = GetFaces(iFace, iNode);
        Coord_FaceElems_CG[iFace][iDim] += val_coord[NodeFace][iDim]/su2double(GetnNodesFace(iFace));
      }
    }

  AD::SetPreaccOut(Coord_CG, nDim);
  AD::SetPreaccOut(Coord_FaceElems_CG, GetnFaces(), nDim);
  AD::EndPreacc();

}

void CPrimalGrid::GetAllNeighbor_Elements() {
  cout << "( ";
  for (unsigned short iFace = 0; iFace < GetnFaces(); iFace++)
  {
    cout << GetNeighbor_Elements(iFace) << ", ";
  }
  cout << ")"  << endl;
}

void CPrimalGrid::InitializeJacobianConstantFaces(unsigned short val_nFaces) {

  /*--- Allocate the memory for JacobianFaceIsConstant and initialize
        its values to false.     ---*/
  JacobianFaceIsConstant = new bool[val_nFaces];
  for(unsigned short i=0; i<val_nFaces; ++i)
    JacobianFaceIsConstant[i] = false;
}

void CPrimalGrid::InitializeNeighbors(unsigned short val_nFaces) {

  /*--- Allocate the memory for Neighbor_Elements and PeriodIndexNeighbors and
        initialize the arrays to -1 to indicate that no neighbor is present and
        that no periodic transformation is needed to the neighbor. ---*/
  Neighbor_Elements    = new long[val_nFaces];
  ElementOwnsFace      = new bool[val_nFaces];
  PeriodIndexNeighbors = new short[val_nFaces];

  for(unsigned short i=0; i<val_nFaces; ++i) {
    Neighbor_Elements[i]    = -1;
    ElementOwnsFace[i]      =  false;
    PeriodIndexNeighbors[i] = -1;
  }
}
