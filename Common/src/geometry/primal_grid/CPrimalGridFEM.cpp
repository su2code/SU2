/*!
 * \file CPrimalGridFEM.cpp
 * \brief Main classes for defining the primal grid elements
 * \author F. Palacios
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/geometry/primal_grid/CPrimalGridFEM.hpp"

CPrimalGridFEM::CPrimalGridFEM(unsigned long val_elemGlobalID, unsigned short val_VTK_Type,
                               unsigned short val_nPolyGrid, unsigned short val_nPolySol, unsigned short val_nDOFsGrid,
                               unsigned short val_nDOFsSol, unsigned long val_offDOfsSol, std::istringstream& elem_line)
    : CPrimalGrid(true, val_nDOFsGrid, nFacesOfElementType(val_VTK_Type)) {
  /*--- Store the integer data in the member variables of this object. ---*/
  VTK_Type = val_VTK_Type;
  nFaces = nFacesOfElementType(VTK_Type);

  nPolyGrid = val_nPolyGrid;
  nPolySol = val_nPolySol;
  nDOFsGrid = val_nDOFsGrid;
  nDOFsSol = val_nDOFsSol;

  elemIDGlobal = val_elemGlobalID;
  offsetDOFsSolGlobal = val_offDOfsSol;

  /*--- Read face structure of the element from elem_line. ---*/

  for (unsigned short i = 0; i < nDOFsGrid; i++) elem_line >> Nodes[i];

  /*--- If a linear element is used, the node numbering for non-simplices
        must be adapted. The reason is that compatability with the original
        SU2 format is maintained for linear elements, but for the FEM solver
        the nodes of the elements are stored row-wise.                       ---*/
  if (nPolyGrid == 1) {
    switch (VTK_Type) {
      case QUADRILATERAL:
        std::swap(Nodes[2], Nodes[3]);
        break;

      case HEXAHEDRON:
        std::swap(Nodes[2], Nodes[3]);
        std::swap(Nodes[6], Nodes[7]);
        break;

      case PYRAMID:
        std::swap(Nodes[2], Nodes[3]);
        break;
    }
  }
}

CPrimalGridFEM::CPrimalGridFEM(unsigned long val_elemGlobalID, unsigned short val_VTK_Type,
                               unsigned short val_nPolyGrid, unsigned short val_nPolySol, unsigned short val_nDOFsGrid,
                               unsigned short val_nDOFsSol, unsigned long val_offDOfsSol, const unsigned long* connGrid)
    : CPrimalGrid(true, val_nDOFsGrid, nFacesOfElementType(val_VTK_Type)) {
  /*--- Store the integer data in the member variables of this object. ---*/
  VTK_Type = val_VTK_Type;
  nFaces = nFacesOfElementType(VTK_Type);

  nPolyGrid = val_nPolyGrid;
  nPolySol = val_nPolySol;
  nDOFsGrid = val_nDOFsGrid;
  nDOFsSol = val_nDOFsSol;

  elemIDGlobal = val_elemGlobalID;
  offsetDOFsSolGlobal = val_offDOfsSol;

  /*--- Copy face structure of the element from connGrid. ---*/
  for (unsigned short i = 0; i < nDOFsGrid; i++) Nodes[i] = connGrid[i];
}

void CPrimalGridFEM::GetLocalCornerPointsAllFaces(unsigned short elementType, unsigned short nPoly,
                                                  unsigned short nDOFs, unsigned short& numFaces,
                                                  unsigned short nPointsPerFace[], unsigned long faceConn[6][4]) {
  /*--- Determine the element type and set the face data accordingly.
        The faceConn values are local to the element.                 ---*/

  numFaces = nFacesOfElementType(elementType);
  unsigned short nn2, nn3, nn4;
  switch (elementType) {
    case TRIANGLE:
      nPointsPerFace[0] = 2;
      faceConn[0][0] = 0;
      faceConn[0][1] = nPoly;
      nPointsPerFace[1] = 2;
      faceConn[1][0] = nPoly;
      faceConn[1][1] = nDOFs - 1;
      nPointsPerFace[2] = 2;
      faceConn[2][0] = nDOFs - 1;
      faceConn[2][1] = 0;
      break;

    case QUADRILATERAL:
      nn2 = nPoly * (nPoly + 1);
      nPointsPerFace[0] = 2;
      faceConn[0][0] = 0;
      faceConn[0][1] = nPoly;
      nPointsPerFace[1] = 2;
      faceConn[1][0] = nPoly;
      faceConn[1][1] = nDOFs - 1;
      nPointsPerFace[2] = 2;
      faceConn[2][0] = nDOFs - 1;
      faceConn[2][1] = nn2;
      nPointsPerFace[3] = 2;
      faceConn[3][0] = nn2;
      faceConn[3][1] = 0;
      break;

    case TETRAHEDRON:
      nn2 = (nPoly + 1) * (nPoly + 2) / 2 - 1;
      nn3 = nDOFs - 1;
      nPointsPerFace[0] = 3;
      faceConn[0][0] = 0;
      faceConn[0][1] = nPoly;
      faceConn[0][2] = nn2;
      nPointsPerFace[1] = 3;
      faceConn[1][0] = 0;
      faceConn[1][1] = nn3;
      faceConn[1][2] = nPoly;
      nPointsPerFace[2] = 3;
      faceConn[2][0] = 0;
      faceConn[2][1] = nn2;
      faceConn[2][2] = nn3;
      nPointsPerFace[3] = 3;
      faceConn[3][0] = nPoly;
      faceConn[3][1] = nn3;
      faceConn[3][2] = nn2;
      break;

    case PYRAMID:
      nn2 = (nPoly + 1) * (nPoly + 1) - 1;
      nn3 = nn2 - nPoly;
      nPointsPerFace[0] = 4;
      faceConn[0][0] = 0;
      faceConn[0][1] = nPoly;
      faceConn[0][2] = nn2;
      faceConn[0][3] = nn3;
      nPointsPerFace[1] = 3;
      faceConn[1][0] = 0;
      faceConn[1][1] = nDOFs - 1;
      faceConn[1][2] = nPoly;
      nPointsPerFace[2] = 3;
      faceConn[2][0] = nn3;
      faceConn[2][1] = nn2;
      faceConn[2][2] = nDOFs - 1;
      nPointsPerFace[3] = 3;
      faceConn[3][0] = 0;
      faceConn[3][1] = nn3;
      faceConn[3][2] = nDOFs - 1;
      nPointsPerFace[4] = 3;
      faceConn[4][0] = nPoly;
      faceConn[4][1] = nDOFs - 1;
      faceConn[4][2] = nn2;
      break;

    case PRISM:
      nn2 = (nPoly + 1) * (nPoly + 2) / 2;
      nn3 = nPoly * nn2;
      --nn2;
      nPointsPerFace[0] = 3;
      faceConn[0][0] = 0;
      faceConn[0][1] = nPoly;
      faceConn[0][2] = nn2;
      nPointsPerFace[1] = 3;
      faceConn[1][0] = nn3;
      faceConn[1][1] = nn2 + nn3;
      faceConn[1][2] = nPoly + nn3;
      nPointsPerFace[2] = 4;
      faceConn[2][0] = 0;
      faceConn[2][1] = nn3;
      faceConn[2][2] = nPoly + nn3;
      faceConn[2][3] = nPoly;
      nPointsPerFace[3] = 4;
      faceConn[3][0] = 0;
      faceConn[3][1] = nn2;
      faceConn[3][2] = nn2 + nn3;
      faceConn[3][3] = nn3;
      nPointsPerFace[4] = 4;
      faceConn[4][0] = nPoly;
      faceConn[4][1] = nPoly + nn3;
      faceConn[4][2] = nn2 + nn3;
      faceConn[4][3] = nn2;
      break;

    case HEXAHEDRON:
      nn2 = (nPoly + 1) * (nPoly + 1);
      nn4 = nPoly * nn2;
      --nn2;
      nn3 = nn2 - nPoly;
      nPointsPerFace[0] = 4;
      faceConn[0][0] = 0;
      faceConn[0][1] = nPoly;
      faceConn[0][2] = nn2;
      faceConn[0][3] = nn3;
      nPointsPerFace[1] = 4;
      faceConn[1][0] = nn4;
      faceConn[1][1] = nn3 + nn4;
      faceConn[1][2] = nn2 + nn4;
      faceConn[1][3] = nPoly + nn4;
      nPointsPerFace[2] = 4;
      faceConn[2][0] = 0;
      faceConn[2][1] = nn4;
      faceConn[2][2] = nPoly + nn4;
      faceConn[2][3] = nPoly;
      nPointsPerFace[3] = 4;
      faceConn[3][0] = nn3;
      faceConn[3][1] = nn2;
      faceConn[3][2] = nn2 + nn4;
      faceConn[3][3] = nn3 + nn4;
      nPointsPerFace[4] = 4;
      faceConn[4][0] = 0;
      faceConn[4][1] = nn3;
      faceConn[4][2] = nn3 + nn4;
      faceConn[4][3] = nn4;
      nPointsPerFace[5] = 4;
      faceConn[5][0] = nPoly;
      faceConn[5][1] = nPoly + nn4;
      faceConn[5][2] = nn2 + nn4;
      faceConn[5][3] = nn2;
      break;
  }
}

void CPrimalGridFEM::GetCornerPointsAllFaces(unsigned short& numFaces, unsigned short nPointsPerFace[],
                                             unsigned long faceConn[6][4]) const {
  /*--- Get the corner points of the faces local to the element. ---*/

  GetLocalCornerPointsAllFaces(VTK_Type, nPolyGrid, nDOFsGrid, numFaces, nPointsPerFace, faceConn);

  /*--- Convert the local values of faceConn to global values. ---*/

  for (unsigned short i = 0; i < numFaces; ++i) {
    for (unsigned short j = 0; j < nPointsPerFace[i]; ++j) {
      unsigned long nn = faceConn[i][j];
      faceConn[i][j] = Nodes[nn];
    }
  }
}
