/*!
 * \file CRectangularMeshReaderFVM.cpp
 * \brief Reads a 2D rectangular grid into linear partitions for the
 *        finite volume solver (FVM).
 * \author T. Economon
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CRectangularMeshReaderFVM.hpp"

CRectangularMeshReaderFVM::CRectangularMeshReaderFVM(const CConfig* val_config, unsigned short val_iZone,
                                                     unsigned short val_nZone)
    : CMeshReaderFVM(val_config, val_iZone, val_nZone) {
  /* The rectangular mesh is always 2D. */
  dimension = 2;

  /* Set the VTK type for the interior elements and the boundary elements. */
  KindElem = QUADRILATERAL;
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

  /* Compute and store the points, interior elements, and surface elements.
   In these routines, we use a simple analytic formula to compute the
   coordinates and the node numbering. We store only the points and interior
   elements on our rank's linear partition, but the master stores the entire
   set of surface connectivity. */
  ComputeRectangularPointCoordinates();
  ComputeRectangularVolumeConnectivity();
  ComputeRectangularSurfaceConnectivity();
}

void CRectangularMeshReaderFVM::ComputeRectangularPointCoordinates() {
  /* Set the global count of points based on the grid dimensions. */
  numberOfGlobalPoints = (nNode) * (mNode);

  /* Get a partitioner to help with linear partitioning. */
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints, 0);

  /* Determine number of local points */
  for (unsigned long globalIndex = 0; globalIndex < numberOfGlobalPoints; globalIndex++) {
    if ((int)pointPartitioner.GetRankContainingIndex(globalIndex) == rank) {
      numberOfLocalPoints++;
    }
  }

  /* Loop over our analytically defined of coordinates and store only
   those that contain a node within our linear partition of points. */
  localPointCoordinates.resize(dimension);
  for (int k = 0; k < dimension; k++) localPointCoordinates[k].reserve(numberOfLocalPoints);
  unsigned long globalIndex = 0;
  for (unsigned long jNode = 0; jNode < mNode; jNode++) {
    for (unsigned long iNode = 0; iNode < nNode; iNode++) {
      if ((int)pointPartitioner.GetRankContainingIndex(globalIndex) == rank) {
        /* Store the coordinates more clearly. */
        const passivedouble x = SU2_TYPE::GetValue(Lx * ((su2double)iNode) / ((su2double)(nNode - 1)) + Ox);
        const passivedouble y = SU2_TYPE::GetValue(Ly * ((su2double)jNode) / ((su2double)(mNode - 1)) + Oy);

        /* Load into the coordinate class data structure. */
        localPointCoordinates[0].push_back(x);
        localPointCoordinates[1].push_back(y);
      }
      globalIndex++;
    }
  }
}

void CRectangularMeshReaderFVM::ComputeRectangularVolumeConnectivity() {
  /* Set the global count of elements based on the grid dimensions. */
  numberOfGlobalElements = (nNode - 1) * (mNode - 1);

  /* Get a partitioner to help with linear partitioning. */
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints, 0);

  /* Loop over our analytically defined of elements and store only those
   that contain a node within our linear partition of points. */
  numberOfLocalElements = 0;
  vector<unsigned long> connectivity(N_POINTS_HEXAHEDRON, 0);
  unsigned long globalIndex = 0;
  for (unsigned long jNode = 0; jNode < mNode - 1; jNode++) {
    for (unsigned long iNode = 0; iNode < nNode - 1; iNode++) {
      /* Compute connectivity based on the i,j index. */
      connectivity[0] = jNode * nNode + iNode;
      connectivity[1] = jNode * nNode + iNode + 1;
      connectivity[2] = (jNode + 1) * nNode + (iNode + 1);
      connectivity[3] = (jNode + 1) * nNode + iNode;

      /* Check whether any of the points is in our linear partition. */
      bool isOwned = false;
      for (unsigned short i = 0; i < N_POINTS_QUADRILATERAL; i++) {
        if ((int)pointPartitioner.GetRankContainingIndex(connectivity[i]) == rank) {
          isOwned = true;
        }
      }

      /* If so, we need to store the element locally. */
      if (isOwned) {
        localVolumeElementConnectivity.push_back(globalIndex);
        localVolumeElementConnectivity.push_back(KindElem);
        for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
          localVolumeElementConnectivity.push_back(connectivity[i]);
        }
        numberOfLocalElements++;
      }
      globalIndex++;
    }
  }
}

void CRectangularMeshReaderFVM::ComputeRectangularSurfaceConnectivity() {
  /* The rectangle alays has 4 markers. */
  numberOfMarkers = 4;
  surfaceElementConnectivity.resize(numberOfMarkers);
  markerNames.resize(numberOfMarkers);

  /* Compute and store the 4 sets of connectivity. */
  markerNames[0] = "y_minus";
  if (rank == MASTER_NODE) {
    for (unsigned long iNode = 0; iNode < nNode - 1; iNode++) {
      surfaceElementConnectivity[0].push_back(0);
      surfaceElementConnectivity[0].push_back(KindBound);
      surfaceElementConnectivity[0].push_back(iNode);
      surfaceElementConnectivity[0].push_back(iNode + 1);
      for (unsigned short i = 0; i < 6; i++) surfaceElementConnectivity[0].push_back(0);
    }
  }

  markerNames[1] = "x_plus";
  if (rank == MASTER_NODE) {
    for (unsigned long jNode = 0; jNode < mNode - 1; jNode++) {
      surfaceElementConnectivity[1].push_back(0);
      surfaceElementConnectivity[1].push_back(KindBound);
      surfaceElementConnectivity[1].push_back(jNode * nNode + (nNode - 1));
      surfaceElementConnectivity[1].push_back((jNode + 1) * nNode + (nNode - 1));
      for (unsigned short i = 0; i < 6; i++) surfaceElementConnectivity[1].push_back(0);
    }
  }

  markerNames[2] = "y_plus";
  if (rank == MASTER_NODE) {
    for (unsigned long iNode = 0; iNode < nNode - 1; iNode++) {
      surfaceElementConnectivity[2].push_back(0);
      surfaceElementConnectivity[2].push_back(KindBound);
      surfaceElementConnectivity[2].push_back((nNode * mNode - 1) - iNode);
      surfaceElementConnectivity[2].push_back((nNode * mNode - 1) - (iNode + 1));
      for (unsigned short i = 0; i < 6; i++) surfaceElementConnectivity[2].push_back(0);
    }
  }

  markerNames[3] = "x_minus";
  if (rank == MASTER_NODE) {
    for (unsigned long jNode = 0; jNode < mNode - 1; jNode++) {
      surfaceElementConnectivity[3].push_back(0);
      surfaceElementConnectivity[3].push_back(KindBound);
      surfaceElementConnectivity[3].push_back((jNode + 1) * nNode);
      surfaceElementConnectivity[3].push_back(jNode * nNode);
      for (unsigned short i = 0; i < 6; i++) surfaceElementConnectivity[3].push_back(0);
    }
  }
}
