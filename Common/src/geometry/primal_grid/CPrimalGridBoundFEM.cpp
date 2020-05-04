/*!
 * \file CPrimalGridBoundFEM.cpp
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

#include "../../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"

CPrimalGridBoundFEM::CPrimalGridBoundFEM(unsigned long         val_elemGlobalID,
                                         unsigned long         val_domainElementID,
                                         unsigned short        val_VTK_Type,
                                         unsigned short        val_nPolyGrid,
                                         unsigned short        val_nDOFsGrid,
                                         vector<unsigned long> &val_nodes)
{
  /*--- Store the integer data in the member variables of this object. ---*/

  VTK_Type = val_VTK_Type;
  nDim = (VTK_Type == LINE) ? 1 : 2;

  nPolyGrid = val_nPolyGrid;
  nDOFsGrid = val_nDOFsGrid;

  boundElemIDGlobal = val_elemGlobalID;
  DomainElement     = val_domainElementID;

  /*--- Allocate the memory for the global nodes of the element to define
        the geometry and copy them from val_nodes.                        ---*/

  Nodes = new unsigned long[nDOFsGrid];
  for(unsigned short i=0; i<nDOFsGrid; i++)
    Nodes[i] = val_nodes[i];

  /*--- For a linear quadrilateral the two last node numbers must be swapped,
        such that the element numbering is consistent with the FEM solver.    ---*/

  if(nPolyGrid == 1 && VTK_Type == QUADRILATERAL) swap(Nodes[2], Nodes[3]);
}

CPrimalGridBoundFEM::~CPrimalGridBoundFEM(){}

void CPrimalGridBoundFEM::GetLocalCornerPointsFace(unsigned short elementType,
                                                   unsigned short nPoly,
                                                   unsigned short nDOFs,
                                                   unsigned short &nPointsPerFace,
                                                   unsigned long  faceConn[]) {

  /*--- Determine the element type and set the face data accordingly.
        The faceConn values are local to the element.                 ---*/

  switch( elementType ) {
    case LINE:
      nPointsPerFace = 2;
      faceConn[0] = 0; faceConn[1] = nPoly;
      break;

    case TRIANGLE:
      nPointsPerFace = 3;
      faceConn[0] = 0; faceConn[1] = nPoly; faceConn[2] = nDOFs -1;
      break;

    case QUADRILATERAL:
      unsigned short nn2 = nPoly*(nPoly+1);
      nPointsPerFace = 4;
      faceConn[0] = 0; faceConn[1] = nPoly; faceConn[2] = nDOFs -1; faceConn[3] = nn2;
      break;
  }
}

void CPrimalGridBoundFEM::GetCornerPointsAllFaces(unsigned short &nFaces,
                                                  unsigned short nPointsPerFace[],
                                                  unsigned long  faceConn[6][4]) {

  /*--- For compatibility reasons with the base class, nFaces is also an
        argument to this function. However, it is always 1, because the
        boundary element is the face.                                   ---*/
  nFaces = 1;

  /*--- Get the corner points of the face local to the element. ---*/

  unsigned long thisFaceConn[4] = {0, 0, 0, 0};
  GetLocalCornerPointsFace(VTK_Type, nPolyGrid, nDOFsGrid,
                           nPointsPerFace[0], thisFaceConn);


  /*--- Convert the local values of thisFaceConn to global values. ---*/

  for(unsigned short j=0; j<nPointsPerFace[0]; ++j)
    faceConn[0][j] = Nodes[thisFaceConn[j]];
}

void CPrimalGridBoundFEM::RemoveMultipleDonorsWallFunctions(void) {

  /* Sort donorElementsWallFunctions in increasing order and remove the
     the double entities. */
  sort(donorElementsWallFunctions.begin(), donorElementsWallFunctions.end());
  vector<unsigned long>::iterator lastEntry;
  lastEntry = unique(donorElementsWallFunctions.begin(),
                     donorElementsWallFunctions.end());
  donorElementsWallFunctions.erase(lastEntry, donorElementsWallFunctions.end());
}
