/*!
 * \file CADTPointsOnlyClass.cpp
 * \brief Class for storing an ADT of only points in an arbitrary number of dimensions.
 * \author E. van der Weide
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

#include "../../include/adt/CADTPointsOnlyClass.hpp"
#include "../../include/parallelization/mpi_structure.hpp"
#include "../../include/option_structure.hpp"

CADTPointsOnlyClass::CADTPointsOnlyClass(unsigned short nDim, unsigned long nPoints, const su2double* coor,
                                         const unsigned long* pointID, const bool globalTree) {
  /* Allocate some thread-safe working variables if required. */
#ifdef HAVE_OMP
  FrontLeaves.resize(omp_get_max_threads());
  FrontLeavesNew.resize(omp_get_max_threads());
#endif

  /*--- Make a distinction between parallel and sequential mode. ---*/

#ifdef HAVE_MPI

  /* Parallel mode. Check whether a global or a local tree must be built. */
  if (globalTree) {
    /*--- A global tree must be built. All points are gathered on all ranks.
          First determine the number of points per rank and store them in such
          a way that the info can be used directly in Allgatherv. ---*/
    int rank, size;
    SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
    SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);

    vector<int> recvCounts(size), displs(size);
    int sizeLocal = (int)nPoints;

    SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, SU2_MPI::GetComm());
    displs[0] = 0;
    for (int i = 1; i < size; ++i) displs[i] = displs[i - 1] + recvCounts[i - 1];

    int sizeGlobal = displs.back() + recvCounts.back();

    /*--- Gather the local pointID's and the ranks of the nodes on all ranks. ---*/
    localPointIDs.resize(sizeGlobal);
    SU2_MPI::Allgatherv(pointID, sizeLocal, MPI_UNSIGNED_LONG, localPointIDs.data(), recvCounts.data(), displs.data(),
                        MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

    ranksOfPoints.resize(sizeGlobal);
    vector<int> rankLocal(sizeLocal, rank);
    SU2_MPI::Allgatherv(rankLocal.data(), sizeLocal, MPI_INT, ranksOfPoints.data(), recvCounts.data(), displs.data(),
                        MPI_INT, SU2_MPI::GetComm());

    /*--- Gather the coordinates of the points on all ranks. ---*/
    for (int i = 0; i < size; ++i) {
      recvCounts[i] *= nDim;
      displs[i] *= nDim;
    }

    coorPoints.resize(nDim * sizeGlobal);
    SU2_MPI::Allgatherv(coor, nDim * sizeLocal, MPI_DOUBLE, coorPoints.data(), recvCounts.data(), displs.data(),
                        MPI_DOUBLE, SU2_MPI::GetComm());
  } else {
    /*--- A local tree must be built. Copy the coordinates and point IDs and
          set the ranks to the rank of this processor. ---*/
    int rank;
    SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

    coorPoints.assign(coor, coor + nDim * nPoints);
    localPointIDs.assign(pointID, pointID + nPoints);
    ranksOfPoints.assign(nPoints, rank);
  }

#else

  /*--- Sequential mode. Copy the coordinates and point IDs and
        set the ranks to MASTER_NODE. ---*/
  coorPoints.assign(coor, coor + nDim * nPoints);
  localPointIDs.assign(pointID, pointID + nPoints);
  ranksOfPoints.assign(nPoints, MASTER_NODE);

#endif

  /*--- Build the tree. ---*/
  BuildADT(nDim, localPointIDs.size(), coorPoints.data());

  /*--- Reserve the memory for frontLeaves and frontLeavesNew,
        which are needed during the tree search. ---*/
  for (auto& vec : FrontLeaves) vec.reserve(200);
  for (auto& vec : FrontLeavesNew) vec.reserve(200);
}

void CADTPointsOnlyClass::DetermineNearestNode_impl(vector<unsigned long>& frontLeaves,
                                                    vector<unsigned long>& frontLeavesNew, const su2double* coor,
                                                    su2double& dist, unsigned long& pointID, int& rankID) const {
  const bool wasActive = AD::BeginPassive();

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Initialize the nearest node to the central node of the     ---*/
  /*---         root leaf. Note that the distance is the distance squared  ---*/
  /*---         to avoid a sqrt.                                           ---*/
  /*--------------------------------------------------------------------------*/

  unsigned long kk = leaves[0].centralNodeID, minIndex;
  const su2double* coorTarget = coorPoints.data() + nDimADT * kk;

  pointID = localPointIDs[kk];
  rankID = ranksOfPoints[kk];
  minIndex = kk;
  dist = 0.0;
  for (unsigned short l = 0; l < nDimADT; ++l) {
    const su2double ds = coor[l] - coorTarget[l];
    dist += ds * ds;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Traverse the tree and search for the nearest node.         ---*/
  /*---         During the tree traversal the currently stored distance    ---*/
  /*---         squared is modified, because the guaranteed minimum        ---*/
  /*---         distance squared of the children could be smaller.         ---*/
  /*--------------------------------------------------------------------------*/

  /* Start at the root leaf of the ADT, i.e. initialize frontLeaves such that
     it only contains the root leaf. Make sure to wipe out any data from a
     previous search. */
  frontLeaves.clear();
  frontLeaves.push_back(0);

  /* Infinite loop of the tree traversal. */
  for (;;) {
    /* Initialize the new front, i.e. the front for the next round, to empty. */
    frontLeavesNew.clear();

    /* Loop over the leaves of the current front. */
    for (unsigned long i = 0; i < frontLeaves.size(); ++i) {
      /* Store the current leaf a bit easier in ll and loop over its children. */
      const unsigned long ll = frontLeaves[i];
      for (unsigned short mm = 0; mm < 2; ++mm) {
        /* Determine whether this child contains a node or a leaf
           of the next level of the ADT. */
        kk = leaves[ll].children[mm];
        if (leaves[ll].childrenAreTerminal[mm]) {
          /*--- Child contains a node. Compute the distance squared to this node
                and store it if this distance squared is less than the currently
                stored value. ---*/
          coorTarget = coorPoints.data() + nDimADT * kk;
          su2double distTarget = 0;
          for (unsigned short l = 0; l < nDimADT; ++l) {
            const su2double ds = coor[l] - coorTarget[l];
            distTarget += ds * ds;
          }

          if (distTarget < dist) {
            dist = distTarget;
            pointID = localPointIDs[kk];
            rankID = ranksOfPoints[kk];
            minIndex = kk;
          }
        } else {
          /*--- Child contains a leaf. Determine the possible minimum distance
                squared to that leaf. ---*/
          su2double posDist = 0.0;
          for (unsigned short l = 0; l < nDimADT; ++l) {
            su2double ds = 0.0;
            if (coor[l] < leaves[kk].xMin[l])
              ds = coor[l] - leaves[kk].xMin[l];
            else if (coor[l] > leaves[kk].xMax[l])
              ds = coor[l] - leaves[kk].xMax[l];

            posDist += ds * ds;
          }

          /*--- Check if the possible minimum distance is less than the currently
                stored minimum distance. If so this leaf must be stored for the
                next round. In that case the distance squared to the central node is
                determined, which is used to update the currently stored value. ---*/
          if (posDist < dist) {
            frontLeavesNew.push_back(kk);

            const unsigned long jj = leaves[kk].centralNodeID;

            coorTarget = coorPoints.data() + nDimADT * jj;
            su2double distTarget = 0;
            for (unsigned short l = 0; l < nDimADT; ++l) {
              const su2double ds = coor[l] - coorTarget[l];
              distTarget += ds * ds;
            }

            if (distTarget < dist) {
              dist = distTarget;
              pointID = localPointIDs[jj];
              rankID = ranksOfPoints[jj];
              minIndex = jj;
            }
          }
        }
      }
    }

    /*--- End of the loop over the current front. Copy the data from
          frontLeavesNew to frontLeaves for the next round. If the new front
          is empty the entire tree has been traversed and a break can be made
          from the infinite loop. ---*/
    frontLeaves = frontLeavesNew;
    if (frontLeaves.empty()) break;
  }

  AD::EndPassive(wasActive);

  /* Recompute the distance to get the correct dependency if we use AD */
  coorTarget = coorPoints.data() + nDimADT * minIndex;
  dist = 0.0;
  for (unsigned short l = 0; l < nDimADT; ++l) {
    const su2double ds = coor[l] - coorTarget[l];
    dist += ds * ds;
  }

  /* At the moment the distance squared to the nearest node is stored.
     Take the sqrt to obtain the correct value. */
  dist = sqrt(dist);
}
