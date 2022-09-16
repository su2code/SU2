/*!
 * \file CMeshReader.cpp
 * \brief Helper class that provides the counts for each rank in a linear
 *        partitioning given the global count as input.
 * \author T. Economon
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/geometry/meshreader/CMeshReaderFVM.hpp"
#include "../../../include/geometry/primal_grid/CPrimalGridFEM.hpp"

CMeshReaderFVM::CMeshReaderFVM(CConfig  *val_config,
                              unsigned short val_iZone,
                              unsigned short val_nZone) :
  rank(SU2_MPI::GetRank()),
  size(SU2_MPI::GetSize()),
  config(val_config) {
}

void CMeshReaderFVM::DetermineFacesVolumeElements(vector<CFaceOfElement> &localFaces) {

  /*--- Loop over the locally stored volume elements. ---*/
  unsigned long ind = 0;
  for(unsigned long k=0; k<numberOfLocalElements; ++k) {

    /*--- Set the pointer where the information of this element
          is stored and determine the number of faces as well
          as the corner points of these faces. Note that 6 faces
          is the maximum for the elements considered. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long  faceConn[6][4];

    const unsigned long *elemInfo = localVolumeElementConnectivity.data() + ind;

    GetCornerPointsAllFaces(elemInfo, nFaces, nPointsPerFace, faceConn);

    /*--- Loop over the faces and add them to localFaces. ---*/
    for(unsigned short i=0; i<nFaces; ++i) {
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for(unsigned short j=0; j<nPointsPerFace[i]; ++j)
        thisFace.cornerPoints[j] = faceConn[i][j];
      thisFace.elemID0 = elemInfo[4];

      thisFace.CreateUniqueNumbering();
      localFaces.push_back(thisFace);
    }

    /*--- Update the index for the next element. ---*/
    const unsigned long nDOFsGrid = localVolumeElementConnectivity[ind+3];
    ind += nDOFsGrid+5;
  }

  /*--- Sort localFaces in increasing order and remove the double entities,
        such that the binary search later on is a bit more efficient. ---*/
  sort(localFaces.begin(), localFaces.end());
  vector<CFaceOfElement>::iterator lastFace;
  lastFace = unique(localFaces.begin(), localFaces.end());
  localFaces.erase(lastFace, localFaces.end());
}

void CMeshReaderFVM::GetCornerPointsAllFaces(const unsigned long *elemInfo,
                                          unsigned short      &numFaces,
                                          unsigned short      nPointsPerFace[],
                                          unsigned long       faceConn[6][4]) {

  /*--- Retrieve the element type, polynomial degree of the grid and
        number of DOFs for this element and set the pointer for the
        connectivity information. ---*/
  const unsigned short VTK_Type  = (const unsigned short) elemInfo[0];
  const unsigned short nPolyGrid = (const unsigned short) elemInfo[1];
  const unsigned short nDOFsGrid = (const unsigned short) elemInfo[3];
  const unsigned long  *conn     = elemInfo + 5;

  /*--- Call the static function GetLocalCornerPointsAllFaces of CPrimalGridFEM
        to determine the local numbering of the corner points of the faces. ---*/
  CPrimalGridFEM::GetLocalCornerPointsAllFaces(VTK_Type, nPolyGrid, nDOFsGrid,
                                               numFaces, nPointsPerFace, faceConn);

  /*--- Convert the local values of faceConn to global values. ---*/
  for(unsigned short i=0; i<numFaces; ++i) {
    for(unsigned short j=0; j<nPointsPerFace[i]; ++j) {
      unsigned long nn = faceConn[i][j];
      faceConn[i][j] = conn[nn];
    }
  }
}
