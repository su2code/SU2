/*!
 * \file CPrimalGrid.cpp
 * \brief Main classes for defining the primal grid elements
 * \author F. Palacios
 * \version 7.0.7 "Blackbird"
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

#include "../../../include/geometry/primal_grid/CPrimalGrid.hpp"

unsigned short CPrimalGrid::nDim;

CPrimalGrid::CPrimalGrid(void) {

  /*--- Set the default values for the pointers ---*/
  Nodes              = nullptr;
  Neighbor_Elements  = nullptr;
  Coord_CG           = nullptr;
  Coord_FaceElems_CG = nullptr;
  GlobalIndex        = 0;

}

CPrimalGrid::~CPrimalGrid() {

 delete[] Nodes;
 delete[] Coord_CG;
 delete[] Neighbor_Elements;
}

void CPrimalGrid::SetCoord_CG(const su2double* const* val_coord) {
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
