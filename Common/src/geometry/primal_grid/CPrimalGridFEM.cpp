/*!
 * \file CPrimalGridFEM.cpp
 * \brief Main classes for defining the primal grid elements
 * \author F. Palacios
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../include/fem/fem_standard_element.hpp"

CPrimalGridFEM::CPrimalGridFEM(const unsigned long* dataElem, unsigned long& offsetSolDOFs)
    : CPrimalGrid(true, dataElem[3], nFacesOfElementType(dataElem[0])) {
  /*--- Store the meta data for this element. ---*/
  VTK_Type = static_cast<unsigned short>(dataElem[0]);
  nPolyGrid = static_cast<unsigned short>(dataElem[1]);
  nPolySol = static_cast<unsigned short>(dataElem[2]);
  nDOFsGrid = static_cast<unsigned short>(dataElem[3]);
  nDOFsSol = CFEMStandardElementBase::GetNDOFsStatic(VTK_Type, nPolySol);
  elemIDGlobal = dataElem[4];

  offsetDOFsSolGlobal = offsetSolDOFs;
  offsetSolDOFs += nDOFsSol;

  nFaces = nFacesOfElementType(VTK_Type);

  /*--- Allocate the memory for the global nodes of the element to define
        the geometry and copy the data from dataElem. ---*/

  for (unsigned short i = 0; i < nDOFsGrid; ++i) Nodes[i] = dataElem[i + 5];
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
