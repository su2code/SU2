/*!
 * \file CPrimalGrid.cpp
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

#include "../../../include/geometry/primal_grid/CPrimalGrid.hpp"

unsigned short CPrimalGrid::nDim;

CPrimalGrid::CPrimalGrid(void) {

  /*--- Set the default values for the pointers ---*/
  Nodes = nullptr;
  Neighbor_Elements = nullptr;
  ElementOwnsFace = nullptr;
  PeriodIndexNeighbors = nullptr;
  JacobianFaceIsConstant = nullptr;
  GlobalIndex = 0;

}

CPrimalGrid::~CPrimalGrid() {

 delete[] Nodes;
 delete[] Neighbor_Elements;
 delete[] ElementOwnsFace;
 delete[] PeriodIndexNeighbors;
 delete[] JacobianFaceIsConstant;
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
