/*!
 * \file CRectangularMeshReaderFVM.cpp
 * \brief Reads a 2D rectangular grid into linear partitions for the
 *        finite element solver (FEM).
 * \author T. Economon, E. van der Weide
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CRectangularMeshReaderFEM.hpp"

CRectangularMeshReaderFEM::CRectangularMeshReaderFEM(CConfig        *val_config,
                                                     unsigned short val_iZone,
                                                     unsigned short val_nZone)
: CMeshReaderFVM(val_config, val_iZone, val_nZone) {
  
  /* The rectangular mesh is always 2D. */
  dimension = 2;

  /* Set the VTK type for the interior elements and the boundary elements. */
  KindElem  = QUADRILATERAL;
  KindBound = LINE;

  /* The number of nodes in the i and j directions. */
  nNode = config->GetMeshBoxSize(0);
  mNode = config->GetMeshBoxSize(1);
  
  /* Lengths for non-square domains. */
  Lx = config->GetMeshBoxLength(0);
  Ly = config->GetMeshBoxLength(1);
  
  /* Offsets in x and y directions from 0.0. */
  Ox = config->GetMeshBoxOffset(0);
  Oy = config->GetMeshBoxOffset(1);

  /* Polynomial degree of the solution. */
  nPolySol = config->GetMeshBoxPSolFEM();

  /*--- Compute and store the interior elements, points, and surface elements
        for this rank, for which simple analytic formulae can be used. ---*/
  ComputeRectangularVolumeConnectivity();
  ComputeRectangularPointCoordinates();
  ComputeRectangularSurfaceConnectivity();
}

CRectangularMeshReaderFEM::~CRectangularMeshReaderFEM(void) { }

void CRectangularMeshReaderFEM::ComputeRectangularPointCoordinates() {

  /*--- Set the global count of points based on the grid dimensions. ---*/
  numberOfGlobalPoints = nNode*mNode;

  /*--- Loop over the local elements to determine the global
        point IDs to be stored on this rank. --*/
  unsigned long ind = 0;
  for(unsigned long i=0; i<numberOfLocalElements; ++i) {

    /*--- Store the number of grid DOFs for this element and
          skip the meta data for this element (5 entries). ---*/
    const unsigned long nDOFsGrid = localVolumeElementConnectivity[ind+3];
    ind += 5;

    /*--- Copy the connectivity to globalPointIDs. ---*/
    unsigned long *conn = localVolumeElementConnectivity.data() + ind;
    ind += nDOFsGrid;
    globalPointIDs.insert(globalPointIDs.end(), conn, conn+nDOFsGrid);
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
  for (int k = 0; k < dimension; k++)
    localPointCoordinates[k].resize(numberOfLocalPoints);

  /*--- Loop over the locally stored points. ---*/
  for(unsigned long i=0; i<numberOfLocalPoints; ++i) {

    /*--- Convert the global index to i,j indices. ---*/
    const unsigned long jNode = globalPointIDs[i]/nNode;
    const unsigned long iNode = globalPointIDs[i] - jNode*nNode;

    /*--- Store the coordinates of the point. ---*/
    localPointCoordinates[0][i] = SU2_TYPE::GetValue(Lx*((su2double)iNode)/((su2double)(nNode-1))+Ox); 
    localPointCoordinates[1][i] = SU2_TYPE::GetValue(Ly*((su2double)jNode)/((su2double)(mNode-1))+Oy);
  }
}

void CRectangularMeshReaderFEM::ComputeRectangularVolumeConnectivity() {

  /*--- Set the global count of elements based on the grid dimensions. ---*/
  const unsigned long nElemI = nNode-1;
  numberOfGlobalElements = nElemI*(mNode-1);

  /*--- Get a partitioner to help with linear partitioning. ---*/
  CLinearPartitioner elemPartitioner(numberOfGlobalElements,0);

  /*--- Determine the index of the first and last element to be stored
        on this rank and the number of local elements. ---*/
  const unsigned long firstIndex = elemPartitioner.GetFirstIndexOnRank(rank);
  const unsigned long lastIndex  = elemPartitioner.GetLastIndexOnRank(rank);
  numberOfLocalElements = elemPartitioner.GetSizeOnRank(rank);

  /*--- Loop over the owned element range. ---*/
  for(unsigned long elem=firstIndex; elem<lastIndex; ++elem) {

    /*--- Retrieve the i,j indices of this element, which are the indices
          of the lower, left point of the element. --*/
    const unsigned long jNode = elem/nElemI;
    const unsigned long iNode = elem - jNode*nElemI;

    /*--- Store the meta data of this element. ---*/
    localVolumeElementConnectivity.push_back(KindElem);
    localVolumeElementConnectivity.push_back(1);         // Pol. degree grid.
    localVolumeElementConnectivity.push_back(nPolySol);
    localVolumeElementConnectivity.push_back(4);         // Number of grid DOFs.
    localVolumeElementConnectivity.push_back(elem);      // Global elem ID.

    /*--- Store the volume connectivity in the format used
          by the FEM solver. ---*/
    localVolumeElementConnectivity.push_back(jNode*nNode + iNode);
    localVolumeElementConnectivity.push_back(jNode*nNode + iNode+1);
    localVolumeElementConnectivity.push_back((jNode+1)*nNode + iNode);
    localVolumeElementConnectivity.push_back((jNode+1)*nNode + iNode+1);
  }
}

void CRectangularMeshReaderFEM::ComputeRectangularSurfaceConnectivity() {

  /*--- Determine the number of elements in i-direction. ---*/
  const unsigned long nElemI = nNode-1;

  /*--- Get a partitioner to help with linear partitioning. ---*/
  CLinearPartitioner elemPartitioner(numberOfGlobalElements,0);

  /*--- The rectangle always has 4 markers. Allocate the required memory. ---*/
  numberOfMarkers = 4;
  numberOfLocalSurfaceElements.resize(numberOfMarkers, 0);
  surfaceElementConnectivity.resize(numberOfMarkers);
  markerNames.resize(numberOfMarkers);

  /*--- Loop over all faces on the yMin (= jMin) boundary. ---*/  
  markerNames[0] = "y_minus";

  for (unsigned long iNode = 0; iNode < nNode-1; ++iNode) {

     /*--- Determine the corresponding global element ID and check
          if it is stored on this rank. ---*/
    const unsigned long globalElemID = iNode;
    if(elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {

      /*--- The corresponding volume element is stored on this rank,
            hence store the surface element as well. ---*/
      surfaceElementConnectivity[0].push_back(KindBound);     // VTK type.
      surfaceElementConnectivity[0].push_back(1);             // Poly degree grid.
      surfaceElementConnectivity[0].push_back(2);             // Number of grid DOFs.
      surfaceElementConnectivity[0].push_back(iNode);         // Global surface element ID.
      surfaceElementConnectivity[0].push_back(globalElemID);  // Global volume element ID.

      surfaceElementConnectivity[0].push_back(iNode);
      surfaceElementConnectivity[0].push_back(iNode+1);

      /*--- Update the number of surface elements for this marker. ---*/
      ++numberOfLocalSurfaceElements[0];
    }
  }

  /*--- Loop over all faces on the xMax (= iMax) boundary. ---*/  
  markerNames[1] = "x_plus";

  for(unsigned long jNode = 0; jNode < mNode-1; ++jNode) {

    /*--- Determine the corresponding global element ID and check
          if it is stored on this rank. ---*/
    const unsigned long globalElemID = jNode*nElemI + nElemI-1;
    if(elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {

      /*--- The corresponding volume element is stored on this rank,
            hence store the surface element as well. ---*/
      surfaceElementConnectivity[1].push_back(KindBound);     // VTK type.
      surfaceElementConnectivity[1].push_back(1);             // Poly degree grid.
      surfaceElementConnectivity[1].push_back(2);             // Number of grid DOFs.
      surfaceElementConnectivity[1].push_back(jNode);         // Global surface element ID.
      surfaceElementConnectivity[1].push_back(globalElemID);  // Global volume element ID.

      surfaceElementConnectivity[1].push_back(jNode*nNode + (nNode-1));
      surfaceElementConnectivity[1].push_back((jNode+1)*nNode + (nNode-1));

      /*--- Update the number of surface elements for this marker. ---*/
      ++numberOfLocalSurfaceElements[1];
    }
  }

  /*--- Loop over all faces on the yMax (= jMax) boundary. ---*/  
  markerNames[2] = "y_plus";

  for (unsigned long iNode = 0; iNode < nNode-1; ++iNode) {

     /*--- Determine the corresponding global element ID and check
          if it is stored on this rank. ---*/
    const unsigned long globalElemID = (mNode-2)*nElemI + iNode;
    if(elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {

      /*--- The corresponding volume element is stored on this rank,
            hence store the surface element as well. ---*/
      surfaceElementConnectivity[2].push_back(KindBound);     // VTK type.
      surfaceElementConnectivity[2].push_back(1);             // Poly degree grid.
      surfaceElementConnectivity[2].push_back(2);             // Number of grid DOFs.
      surfaceElementConnectivity[2].push_back(iNode);         // Global surface element ID.
      surfaceElementConnectivity[2].push_back(globalElemID);  // Global volume element ID.

      surfaceElementConnectivity[2].push_back((mNode-1)*nNode + iNode+1);
      surfaceElementConnectivity[2].push_back((mNode-1)*nNode + iNode);

      /*--- Update the number of surface elements for this marker. ---*/
      ++numberOfLocalSurfaceElements[2];
    }
  }

  /*--- Loop over all faces on the xMin (= iMin) boundary. ---*/  
  markerNames[3] = "x_minus";

  for(unsigned long jNode = 0; jNode < mNode-1; ++jNode) {

    /*--- Determine the corresponding global element ID and check
          if it is stored on this rank. ---*/
    const unsigned long globalElemID = jNode*nElemI;
    if(elemPartitioner.GetRankContainingIndex(globalElemID) == static_cast<unsigned long>(rank)) {

      /*--- The corresponding volume element is stored on this rank,
            hence store the surface element as well. ---*/
      surfaceElementConnectivity[3].push_back(KindBound);     // VTK type.
      surfaceElementConnectivity[3].push_back(1);             // Poly degree grid.
      surfaceElementConnectivity[3].push_back(2);             // Number of grid DOFs.
      surfaceElementConnectivity[3].push_back(jNode);         // Global surface element ID.
      surfaceElementConnectivity[3].push_back(globalElemID);  // Global volume element ID.

      surfaceElementConnectivity[3].push_back((jNode+1)*nNode);
      surfaceElementConnectivity[3].push_back(jNode*nNode);

      /*--- Update the number of surface elements for this marker. ---*/
      ++numberOfLocalSurfaceElements[3];
    }
  }
}
