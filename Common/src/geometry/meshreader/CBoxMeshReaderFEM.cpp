/*!
 * \file CBoxMeshReaderFEM.cpp
 * \brief Reads a 3D box grid into linear partitions for the
 *        finite element solver (FEM).
 * \author T. Economon, E. van der Weide
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CBoxMeshReaderFEM.hpp"

CBoxMeshReaderFEM::CBoxMeshReaderFEM(const CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone)
    : CMeshReaderBase(val_config, val_iZone, val_nZone) {
  /* The box mesh is always 3D. */
  dimension = 3;

  /* Set the VTK type for the interior elements and the boundary elements. */
  KindElem = HEXAHEDRON;
  KindBound = QUADRILATERAL;

  /* The number of grid nodes in the i, j and k directions. */
  nNode = config->GetMeshBoxSize(0);
  mNode = config->GetMeshBoxSize(1);
  pNode = config->GetMeshBoxSize(2);

  /* Lengths for non-square domains. */
  Lx = config->GetMeshBoxLength(0);
  Ly = config->GetMeshBoxLength(1);
  Lz = config->GetMeshBoxLength(2);

  /* Offsets in x, y and z directions from 0.0. */
  Ox = config->GetMeshBoxOffset(0);
  Oy = config->GetMeshBoxOffset(1);
  Oz = config->GetMeshBoxOffset(2);

  /* Polynomial degree of the solution. */
  nPolySol = config->GetMeshBoxPSolFEM();

  /*--- Compute and store the interior elements, points, and surface elements
        for this rank, for which simple analytic formulae can be used. ---*/
  ComputeBoxVolumeConnectivity();
  ComputeBoxPointCoordinates();
  ComputeBoxSurfaceConnectivity();
}

CBoxMeshReaderFEM::~CBoxMeshReaderFEM() = default;

void CBoxMeshReaderFEM::ComputeBoxPointCoordinates() {
  /*--- Set the global count of points based on the grid dimensions. ---*/
  numberOfGlobalPoints = nNode * mNode * pNode;

  /*--- Loop over the local elements to determine the global
        point IDs to be stored on this rank. --*/
  unsigned long ind = 0;
  for (unsigned long i = 0; i < numberOfLocalElements; ++i) {
    /*--- Store the number of grid DOFs for this element and
          skip the meta data for this element (5 entries). ---*/
    const unsigned long nDOFsGrid = localVolumeElementConnectivity[ind + 3];
    ind += 5;

    /*--- Copy the connectivity to globalPointIDs. ---*/
    unsigned long* conn = localVolumeElementConnectivity.data() + ind;
    ind += nDOFsGrid;
    globalPointIDs.insert(globalPointIDs.end(), conn, conn + nDOFsGrid);
  }

  /*--- Sort the globalPointIDs and remove the duplicate entries. ---*/
  sort(globalPointIDs.begin(), globalPointIDs.end());
  vector<unsigned long>::iterator lastNode;
  lastNode = unique(globalPointIDs.begin(), globalPointIDs.end());
  globalPointIDs.erase(lastNode, globalPointIDs.end());

  /*--- Determine the number of locally stored points. ---*/
  numberOfLocalPoints = globalPointIDs.size();

  /*--- Allocate the memory for the locally stored points. ---*/
  localPointCoordinates.resize(dimension);
  for (int k = 0; k < dimension; k++) localPointCoordinates[k].resize(numberOfLocalPoints);

  /*--- Loop over the locally stored points. ---*/
  const unsigned long nPointIJ = nNode * mNode;
  for (unsigned long i = 0; i < numberOfLocalPoints; ++i) {
    /*--- Convert the global index to i,j,k indices. ---*/
    const unsigned long kNode = globalPointIDs[i] / nPointIJ;
    const unsigned long rem = globalPointIDs[i] - kNode * nPointIJ;
    const unsigned long jNode = rem / nNode;
    const unsigned long iNode = rem - jNode * nNode;

    /*--- Store the coordinates of the point. ---*/
    localPointCoordinates[0][i] = SU2_TYPE::GetValue(Lx * ((su2double)iNode) / ((su2double)(nNode - 1)) + Ox);
    localPointCoordinates[1][i] = SU2_TYPE::GetValue(Ly * ((su2double)jNode) / ((su2double)(mNode - 1)) + Oy);
    localPointCoordinates[2][i] = SU2_TYPE::GetValue(Lz * ((su2double)kNode) / ((su2double)(pNode - 1)) + Oz);
  }
}

void CBoxMeshReaderFEM::ComputeBoxVolumeConnectivity() {
  /*--- Set the global count of elements based on the grid dimensions. ---*/
  const unsigned long nElemI = nNode - 1;
  const unsigned long nElemIJ = nElemI * (mNode - 1);
  numberOfGlobalElements = nElemIJ * (pNode - 1);

  /*--- Get a partitioner to help with linear partitioning. ---*/
  CLinearPartitioner elemPartitioner(numberOfGlobalElements, 0);

  /*--- Determine the index of the first and last element to be stored
        on this rank and the number of local elements. ---*/
  const unsigned long firstIndex = elemPartitioner.GetFirstIndexOnRank(rank);
  const unsigned long lastIndex = elemPartitioner.GetLastIndexOnRank(rank);
  numberOfLocalElements = elemPartitioner.GetSizeOnRank(rank);

  /*--- Loop over the owned element range. ---*/
  for (unsigned long elem = firstIndex; elem < lastIndex; ++elem) {
    /*--- Retrieve the i,j,k indices of this element, which are the indices
          of the lower, left, back point of the element. --*/
    const unsigned long kNode = elem / nElemIJ;
    const unsigned long rem = elem - kNode * nElemIJ;
    const unsigned long jNode = rem / nElemI;
    const unsigned long iNode = rem - jNode * nElemI;

    /*--- Store the meta data of this element. ---*/
    localVolumeElementConnectivity.push_back(KindElem);
    localVolumeElementConnectivity.push_back(1);  // Pol. degree grid.
    localVolumeElementConnectivity.push_back(nPolySol);
    localVolumeElementConnectivity.push_back(8);     // Number of grid DOFs.
    localVolumeElementConnectivity.push_back(elem);  // Global elem ID.

    /*--- Store the volume connectivity in the format used
          by the FEM solver. ---*/
    localVolumeElementConnectivity.push_back(kNode * mNode * nNode + jNode * nNode + iNode);
    localVolumeElementConnectivity.push_back(kNode * mNode * nNode + jNode * nNode + iNode + 1);
    localVolumeElementConnectivity.push_back(kNode * mNode * nNode + (jNode + 1) * nNode + iNode);
    localVolumeElementConnectivity.push_back(kNode * mNode * nNode + (jNode + 1) * nNode + iNode + 1);
    localVolumeElementConnectivity.push_back((kNode + 1) * mNode * nNode + jNode * nNode + iNode);
    localVolumeElementConnectivity.push_back((kNode + 1) * mNode * nNode + jNode * nNode + iNode + 1);
    localVolumeElementConnectivity.push_back((kNode + 1) * mNode * nNode + (jNode + 1) * nNode + iNode);
    localVolumeElementConnectivity.push_back((kNode + 1) * mNode * nNode + (jNode + 1) * nNode + iNode + 1);
  }
}

void CBoxMeshReaderFEM::ComputeBoxSurfaceConnectivity() {
  /*--- Determine the number of elements in i-direction and the
        number of elements on a k-face. ---*/
  const unsigned long nElemI = nNode - 1;
  const unsigned long nElemIJ = nElemI * (mNode - 1);

  /*--- Get a partitioner to help with linear partitioning. ---*/
  CLinearPartitioner elemPartitioner(numberOfGlobalElements, 0);

  /*--- The box always has 6 markers. Allocate the required memory. ---*/
  numberOfMarkers = 6;
  numberOfLocalSurfaceElements.resize(numberOfMarkers, 0);
  surfaceElementConnectivity.resize(numberOfMarkers);
  markerNames.resize(numberOfMarkers);

  /*--- Loop over all faces on the xMin (= iMin) boundary. ---*/
  markerNames[0] = "x_minus";

  unsigned long ind = 0;
  for (unsigned long kNode = 0; kNode < pNode - 1; ++kNode) {
    for (unsigned long jNode = 0; jNode < mNode - 1; ++jNode, ++ind) {
      /*--- Determine the corresponding global element ID and check
            if it is stored on this rank. ---*/
      const unsigned long globalElemID = kNode * nElemIJ + jNode * nElemI;
      if (elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {
        /*--- The corresponding volume element is stored on this rank,
              hence store the surface element as well. ---*/
        surfaceElementConnectivity[0].push_back(KindBound);     // VTK type.
        surfaceElementConnectivity[0].push_back(1);             // Poly degree grid.
        surfaceElementConnectivity[0].push_back(4);             // Number of grid DOFs.
        surfaceElementConnectivity[0].push_back(ind);           // Global surface element ID.
        surfaceElementConnectivity[0].push_back(globalElemID);  // Global volume element ID.

        surfaceElementConnectivity[0].push_back(kNode * mNode * nNode + jNode * nNode);
        surfaceElementConnectivity[0].push_back((kNode + 1) * mNode * nNode + jNode * nNode);
        surfaceElementConnectivity[0].push_back(kNode * mNode * nNode + (jNode + 1) * nNode);
        surfaceElementConnectivity[0].push_back((kNode + 1) * mNode * nNode + (jNode + 1) * nNode);

        /*--- Update the number of surface elements for this marker. ---*/
        ++numberOfLocalSurfaceElements[0];
      }
    }
  }

  /*--- Loop over all faces on the xMax (= iMax) boundary. ---*/
  markerNames[1] = "x_plus";

  ind = 0;
  for (unsigned long kNode = 0; kNode < pNode - 1; ++kNode) {
    for (unsigned long jNode = 0; jNode < mNode - 1; ++jNode, ++ind) {
      /*--- Determine the corresponding global element ID and check
            if it is stored on this rank. ---*/
      const unsigned long globalElemID = kNode * nElemIJ + jNode * nElemI + nElemI - 1;
      if (elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {
        /*--- The corresponding volume element is stored on this rank,
              hence store the surface element as well. ---*/
        surfaceElementConnectivity[1].push_back(KindBound);     // VTK type.
        surfaceElementConnectivity[1].push_back(1);             // Poly degree grid.
        surfaceElementConnectivity[1].push_back(4);             // Number of grid DOFs.
        surfaceElementConnectivity[1].push_back(ind);           // Global surface element ID.
        surfaceElementConnectivity[1].push_back(globalElemID);  // Global volume element ID.

        surfaceElementConnectivity[1].push_back(kNode * mNode * nNode + jNode * nNode + (nNode - 1));
        surfaceElementConnectivity[1].push_back(kNode * mNode * nNode + (jNode + 1) * nNode + (nNode - 1));
        surfaceElementConnectivity[1].push_back((kNode + 1) * mNode * nNode + jNode * nNode + (nNode - 1));
        surfaceElementConnectivity[1].push_back((kNode + 1) * mNode * nNode + (jNode + 1) * nNode + (nNode - 1));

        /*--- Update the number of surface elements for this marker. ---*/
        ++numberOfLocalSurfaceElements[1];
      }
    }
  }

  /*--- Loop over all faces on the yMin (= jMin) boundary. ---*/
  markerNames[2] = "y_minus";

  ind = 0;
  for (unsigned long kNode = 0; kNode < pNode - 1; ++kNode) {
    for (unsigned long iNode = 0; iNode < nNode - 1; ++iNode, ++ind) {
      /*--- Determine the corresponding global element ID and check
           if it is stored on this rank. ---*/
      const unsigned long globalElemID = kNode * nElemIJ + iNode;
      if (elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {
        /*--- The corresponding volume element is stored on this rank,
              hence store the surface element as well. ---*/
        surfaceElementConnectivity[2].push_back(KindBound);     // VTK type.
        surfaceElementConnectivity[2].push_back(1);             // Poly degree grid.
        surfaceElementConnectivity[2].push_back(4);             // Number of grid DOFs.
        surfaceElementConnectivity[2].push_back(ind);           // Global surface element ID.
        surfaceElementConnectivity[2].push_back(globalElemID);  // Global volume element ID.

        surfaceElementConnectivity[2].push_back(kNode * mNode * nNode + iNode);
        surfaceElementConnectivity[2].push_back(kNode * mNode * nNode + iNode + 1);
        surfaceElementConnectivity[2].push_back((kNode + 1) * mNode * nNode + iNode);
        surfaceElementConnectivity[2].push_back((kNode + 1) * mNode * nNode + iNode + 1);

        /*--- Update the number of surface elements for this marker. ---*/
        ++numberOfLocalSurfaceElements[2];
      }
    }
  }

  /*--- Loop over all faces on the yMax (= jMax) boundary. ---*/
  markerNames[3] = "y_plus";

  ind = 0;
  for (unsigned long kNode = 0; kNode < pNode - 1; ++kNode) {
    for (unsigned long iNode = 0; iNode < nNode - 1; ++iNode, ++ind) {
      /*--- Determine the corresponding global element ID and check
           if it is stored on this rank. ---*/
      const unsigned long globalElemID = kNode * nElemIJ + (mNode - 2) * nElemI + iNode;
      if (elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {
        /*--- The corresponding volume element is stored on this rank,
              hence store the surface element as well. ---*/
        surfaceElementConnectivity[3].push_back(KindBound);     // VTK type.
        surfaceElementConnectivity[3].push_back(1);             // Poly degree grid.
        surfaceElementConnectivity[3].push_back(4);             // Number of grid DOFs.
        surfaceElementConnectivity[3].push_back(ind);           // Global surface element ID.
        surfaceElementConnectivity[3].push_back(globalElemID);  // Global volume element ID.

        surfaceElementConnectivity[3].push_back(kNode * mNode * nNode + (mNode - 1) * nNode + iNode);
        surfaceElementConnectivity[3].push_back(kNode * mNode * nNode + (mNode - 1) * nNode + iNode + 1);
        surfaceElementConnectivity[3].push_back((kNode + 1) * mNode * nNode + (mNode - 1) * nNode + iNode);
        surfaceElementConnectivity[3].push_back((kNode + 1) * mNode * nNode + (mNode - 1) * nNode + iNode + 1);

        /*--- Update the number of surface elements for this marker. ---*/
        ++numberOfLocalSurfaceElements[3];
      }
    }
  }

  /*--- Loop over all faces on the zMin (= kMin) boundary. ---*/
  markerNames[4] = "z_minus";

  ind = 0;
  for (unsigned long jNode = 0; jNode < mNode - 1; ++jNode) {
    for (unsigned long iNode = 0; iNode < nNode - 1; ++iNode, ++ind) {
      /*--- Determine the corresponding global element ID and check
           if it is stored on this rank. ---*/
      const unsigned long globalElemID = jNode * nElemI + iNode;
      if (elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {
        /*--- The corresponding volume element is stored on this rank,
              hence store the surface element as well. ---*/
        surfaceElementConnectivity[4].push_back(KindBound);     // VTK type.
        surfaceElementConnectivity[4].push_back(1);             // Poly degree grid.
        surfaceElementConnectivity[4].push_back(4);             // Number of grid DOFs.
        surfaceElementConnectivity[4].push_back(ind);           // Global surface element ID.
        surfaceElementConnectivity[4].push_back(globalElemID);  // Global volume element ID.

        surfaceElementConnectivity[4].push_back(jNode * nNode + iNode);
        surfaceElementConnectivity[4].push_back(jNode * nNode + iNode + 1);
        surfaceElementConnectivity[4].push_back((jNode + 1) * nNode + iNode);
        surfaceElementConnectivity[4].push_back((jNode + 1) * nNode + iNode + 1);

        /*--- Update the number of surface elements for this marker. ---*/
        ++numberOfLocalSurfaceElements[4];
      }
    }
  }

  /*--- Loop over all faces on the zMax (= kMax) boundary. ---*/
  markerNames[5] = "z_plus";

  ind = 0;
  for (unsigned long jNode = 0; jNode < mNode - 1; ++jNode) {
    for (unsigned long iNode = 0; iNode < nNode - 1; ++iNode, ++ind) {
      /*--- Determine the corresponding global element ID and check
           if it is stored on this rank. ---*/
      const unsigned long globalElemID = (pNode - 2) * nElemIJ + jNode * nElemI + iNode;
      if (elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {
        /*--- The corresponding volume element is stored on this rank,
              hence store the surface element as well. ---*/
        surfaceElementConnectivity[5].push_back(KindBound);     // VTK type.
        surfaceElementConnectivity[5].push_back(1);             // Poly degree grid.
        surfaceElementConnectivity[5].push_back(4);             // Number of grid DOFs.
        surfaceElementConnectivity[5].push_back(ind);           // Global surface element ID.
        surfaceElementConnectivity[5].push_back(globalElemID);  // Global volume element ID.

        surfaceElementConnectivity[5].push_back((pNode - 1) * mNode * nNode + jNode * nNode + iNode);
        surfaceElementConnectivity[5].push_back((pNode - 1) * mNode * nNode + jNode * nNode + iNode + 1);
        surfaceElementConnectivity[5].push_back((pNode - 1) * mNode * nNode + (jNode + 1) * nNode + iNode);
        surfaceElementConnectivity[5].push_back((pNode - 1) * mNode * nNode + (jNode + 1) * nNode + iNode + 1);

        /*--- Update the number of surface elements for this marker. ---*/
        ++numberOfLocalSurfaceElements[5];
      }
    }
  }
}
