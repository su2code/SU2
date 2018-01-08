/*!
 * \file adt_structure.cpp
 * \brief Main subroutines for for carrying out geometrical searches using an
 *        alternating digital tree (ADT).
 * \author E. van der Weide
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/adt_structure.hpp"

/* Define the tolerance to decide whether or not a point is inside an element. */
const su2double tolInsideElem   =  1.e-10;
const su2double paramLowerBound = -1.0 - tolInsideElem;
const su2double paramUpperBound =  1.0 + tolInsideElem;

su2_adtComparePointClass::su2_adtComparePointClass(const su2double      *coor,
                                                   const unsigned short splitDir,
                                                   const unsigned short nDimADT)
  : pointCoor(coor),
    splitDirection(splitDir),
    nDim(nDimADT) {}

su2_BBoxTargetClass::su2_BBoxTargetClass(const unsigned long val_BBoxID,
                                         const su2double     val_posDist2,
                                         const su2double     val_guarDist2) {
  boundingBoxID      = val_BBoxID;
  possibleMinDist2   = val_posDist2;
  guaranteedMinDist2 = val_guarDist2;
}

bool su2_BBoxTargetClass::operator <(const su2_BBoxTargetClass &other) const {

  /* Make sure that the bounding boxes with the smallest possible distances
     are stored first. */
  if(possibleMinDist2   < other.possibleMinDist2)   return true;
  if(possibleMinDist2   > other.possibleMinDist2)   return false;
  if(guaranteedMinDist2 < other.guaranteedMinDist2) return true;
  if(guaranteedMinDist2 > other.guaranteedMinDist2) return false;
  if(boundingBoxID      < other.boundingBoxID)      return true;
  return false;
}

void su2_BBoxTargetClass::Copy(const su2_BBoxTargetClass &other) {
  boundingBoxID      = other.boundingBoxID;
  possibleMinDist2   = other.possibleMinDist2;
  guaranteedMinDist2 = other.guaranteedMinDist2;
}

void su2_adtNodeClass::Copy(const su2_adtNodeClass &other) {

  childrenAreTerminal[0] = other.childrenAreTerminal[0];
  childrenAreTerminal[1] = other.childrenAreTerminal[1];

  children[0] = other.children[0];
  children[1] = other.children[1];

  centralNodeID = other.centralNodeID;

  xMin = other.xMin;
  xMax = other.xMax;
}

void su2_adtBaseClass::BuildADT(unsigned short  nDim,
                                unsigned long   nPoints,
                                const su2double *coor) {

  /*---  Determine the number of leaves. It can be proved that nLeaves equals
         nPoints-1 for an optimally balanced tree. Take the exceptional case of
         nPoints == 1 into account and return if the tree is empty. ---*/
  nDimADT = nDim;
  isEmpty = false;
  nLeaves = nPoints -1;
  if(nPoints <= 1) ++nLeaves;
  if(nLeaves == 0) {isEmpty = true; return;}

  /*--- Allocate the memory for the leaves of the ADT and the minimum and
        maximum coordinates of the leaves. Note that these coordinates are
        stored in one vector, rather than that memory is allocated for the
        individual leaves. ---*/
  leaves.resize(nLeaves);
  coorMinLeaves.resize(nDim*nLeaves);
  coorMaxLeaves.resize(nDim*nLeaves);

  /*--- Define the vectors, which control the subdivision of the leaves. ---*/
  unsigned long nn = (nPoints+1)/2;
  vector<unsigned long> pointIDs(nPoints), pointIDsNew(nPoints);
  vector<unsigned long> nPointIDs(nn+1),   nPointIDsNew(nn+1);
  vector<unsigned long> curLeaf(nn),       curLeafNew(nn);

  /*--------------------------------------------------------------------------*/
  /*---                 Building of the actual ADT                         ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Initialize the arrays pointIDs, nPointIDs and curLeaf such that all
        points belong to the root leaf. Also set the counters nLeavesToDivide
        and nLeavesTot.  ---*/
  nPointIDs[0] = 0; nPointIDs[1] = nPoints;
  curLeaf[0] = 0;

  for(unsigned long i=0; i<nPoints; ++i) pointIDs[i] = i;

  unsigned long nLeavesToDivide = 1, nLeavesTot = 1;

  /*--- Loop to subdivide the leaves. The division is such that the ADT is
        optimally balanced.  ---*/
  for(;;) {

    /* Criterion to exit the loop. */
    if(nLeavesToDivide == 0) break;

    /* Initializations for the next round of subdivisions. */
    unsigned long nLeavesToDivideNew = 0;
    nPointIDsNew[0] = 0;

    /*--- Loop over the current number of leaves to be divided. ---*/
    for(unsigned long i=0; i<nLeavesToDivide; ++i) {

      /* Store the number of points present in the leaf in nn and the
         current leaf number in mm. */
      nn = nPointIDs[i+1] - nPointIDs[i];
      unsigned long mm = curLeaf[i];

      /*--- Set the pointers for the coordinates of the leaf to the correct
            locations in the vectors coorMinLeaves and coorMaxLeaves and
            determine the bounding box coordinates of this leaf. ---*/
      leaves[mm].xMin = coorMinLeaves.data() + nDim*mm;
      leaves[mm].xMax = coorMaxLeaves.data() + nDim*mm;

      unsigned long ll = nDim*pointIDs[nPointIDs[i]];
      for(unsigned short l=0; l<nDim; ++l)
        leaves[mm].xMin[l] = leaves[mm].xMax[l] = coor[ll+l];

      for(unsigned long j=(nPointIDs[i]+1); j<nPointIDs[i+1]; ++j) {
        ll = nDim*pointIDs[j];
        for(unsigned short l=0; l<nDim; ++l) {
          leaves[mm].xMin[l] = min(leaves[mm].xMin[l], coor[ll+l]);
          leaves[mm].xMax[l] = max(leaves[mm].xMax[l], coor[ll+l]);
        }
      }

      /*--- Determine the split direction for this leaf. The splitting is done
            in such a way that isotropy is reached as quickly as possible.
            Hence the split direction is the direction of the largest dimension
            of the leaf. ---*/
      unsigned short splitDir= 0;
      su2double distMax = -1.0;
      for(unsigned short l=0; l<nDim; ++l) {
        const su2double dist = leaves[mm].xMax[l] - leaves[mm].xMin[l];
        if(dist > distMax) {distMax = dist; splitDir = l;}
      }

      /* Sort the points of the current leaf in increasing order. The sorting
         is based on the coordinate in the split direction, for which the
         functor su2_adtComparePointClass is used. */
      sort(pointIDs.data() + nPointIDs[i], pointIDs.data() + nPointIDs[i+1],
           su2_adtComparePointClass(coor, splitDir, nDim));

      /* Determine the index of the node, which is approximately central
         in this leave. */
      leaves[mm].centralNodeID = pointIDs[nPointIDs[i] + nn/2];

      /*--- Determine the situation of the leaf. It is either a terminal leaf
            or a leaf that must be subdivided. ---*/
      if(nn <= 2) {

        /* Terminal leaf. Store the ID's of the points as children and
           indicate that the children are terminal. */
        leaves[mm].children[0] = pointIDs[nPointIDs[i]];
        leaves[mm].children[1] = pointIDs[nPointIDs[i+1]-1];

        leaves[mm].childrenAreTerminal[0] = true;
        leaves[mm].childrenAreTerminal[1] = true;
      }
      else {

        /* The leaf must be divided. Determine the number of points in the
           left leaf. This number is at least 2. The actual number stored in kk
           is this number plus an offset. Also initialize the counter nfl, which
           is used to store the bounding boxes in the arrays for the new round. */
        unsigned long kk  = (nn+1)/2 + nPointIDs[i];
        unsigned long nfl = nPointIDsNew[nLeavesToDivideNew];

        /* Copy the ID's of the left points into pointIDsNew. Also update the
           corresponding entry in nPointIDsNew. */
        for(unsigned long k=nPointIDs[i]; k<kk; ++k)
          pointIDsNew[nfl++] = pointIDs[k];

        nPointIDsNew[nLeavesToDivideNew+1] = nfl;

        /* Store the leaf info in the tree and in the leafs for next round and
           update nLeavesToDivideNew and nLeavesTot. */
        leaves[mm].children[0]            = nLeavesTot;
        leaves[mm].childrenAreTerminal[0] = false;

        curLeafNew[nLeavesToDivideNew] = nLeavesTot;
        ++nLeavesToDivideNew;
        ++nLeavesTot;

        /*--- The right leaf will only be created if it has more than one point
              in it, i.e. if the original leaf has more than three points.
              If the new leaf only has one point in it, it is not created.
              Instead, the point is stored in the current leaf. ---*/
        if(nn == 3)
        {
          /* Only three points present in the current leaf. The right leaf is
             not created and the last point is stored as the second child of
             the current leaf. */
          leaves[mm].children[1]            = pointIDs[nPointIDs[i+1]-1];
          leaves[mm].childrenAreTerminal[1] = true;
        }
        else {

          /* More than 3 points are present and thus the right leaf is created.
             Copy the ID's from pointIDs into pointIDsNew and update the
             counters for the new round. */
          unsigned long nfr = nPointIDsNew[nLeavesToDivideNew];

          for(unsigned long k=kk; k<nPointIDs[i+1]; ++k)
            pointIDsNew[nfr++] = pointIDs[k];

          nPointIDsNew[nLeavesToDivideNew+1] = nfr;

          /* Store the leaf info in the tree and in the leaves for next round and
             update nLeavesToDivideNew and nLeavesTot. */
          leaves[mm].children[1]            = nLeavesTot;
          leaves[mm].childrenAreTerminal[1] = false;

          curLeafNew[nLeavesToDivideNew] = nLeavesTot;
          ++nLeavesToDivideNew;
          ++nLeavesTot;
        }
      }
    }

    /* Set the data for the next round. */
    nLeavesToDivide = nLeavesToDivideNew;

    for(unsigned long i=0; i<=nLeavesToDivide; ++i)           nPointIDs[i] = nPointIDsNew[i];
    for(unsigned long i=0; i< nLeavesToDivide; ++i)           curLeaf[i]   = curLeafNew[i];
    for(unsigned long i=0; i<nPointIDs[nLeavesToDivide]; ++i) pointIDs[i]  = pointIDsNew[i];
  }
}

su2_adtPointsOnlyClass::su2_adtPointsOnlyClass(unsigned short nDim,
                                               unsigned long  nPoints,
                                               su2double      *coor,
                                               unsigned long  *pointID,
                                               const bool     globalTree) {

  /*--- Make a distinction between parallel and sequential mode. ---*/

#ifdef HAVE_MPI

  /* Parallel mode. Check whether a global or a local tree must be built. */
  if( globalTree ) {

    /*--- A global tree must be built. All points are gathered on all ranks.
          First determine the number of points per rank and store them in such
          a way that the info can be used directly in Allgatherv. ---*/
    int rank, size;
    SU2_MPI::Comm_rank(MPI_COMM_WORLD, &rank);
    SU2_MPI::Comm_size(MPI_COMM_WORLD, &size);

    vector<int> recvCounts(size), displs(size);
    int sizeLocal = (int) nPoints;

    SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1,
                       MPI_INT, MPI_COMM_WORLD);
    displs[0] = 0;
    for(int i=1; i<size; ++i) displs[i] = displs[i-1] + recvCounts[i-1];

    int sizeGlobal = displs.back() + recvCounts.back();

    /*--- Gather the local pointID's and the ranks of the nodes on all ranks. ---*/
    localPointIDs.resize(sizeGlobal);
    SU2_MPI::Allgatherv(pointID, sizeLocal, MPI_UNSIGNED_LONG, localPointIDs.data(),
                        recvCounts.data(), displs.data(), MPI_UNSIGNED_LONG,
                        MPI_COMM_WORLD);

    ranksOfPoints.resize(sizeGlobal);
    vector<int> rankLocal(sizeLocal, rank);
    SU2_MPI::Allgatherv(rankLocal.data(), sizeLocal, MPI_INT, ranksOfPoints.data(),
                        recvCounts.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);

    /*--- Gather the coordinates of the points on all ranks. ---*/
    for(int i=0; i<size; ++i) {recvCounts[i] *= nDim; displs[i] *= nDim;}

    coorPoints.resize(nDim*sizeGlobal);
    SU2_MPI::Allgatherv(coor, nDim*sizeLocal, MPI_DOUBLE, coorPoints.data(),
                        recvCounts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  }
  else {

    /*--- A local tree must be built. Copy the coordinates and point IDs and
          set the ranks to the rank of this processor. ---*/
    int rank;
    SU2_MPI::Comm_rank(MPI_COMM_WORLD, &rank);

    coorPoints.assign(coor, coor + nDim*nPoints);
    localPointIDs.assign(pointID, pointID + nPoints);
    ranksOfPoints.assign(nPoints, rank);
  }
  
#else

  /*--- Sequential mode. Copy the coordinates and point IDs and
        set the ranks to MASTER_NODE. ---*/
  coorPoints.assign(coor, coor + nDim*nPoints);
  localPointIDs.assign(pointID, pointID + nPoints);
  ranksOfPoints.assign(nPoints, MASTER_NODE);

#endif

  /*--- Build the tree. ---*/
  BuildADT(nDim, localPointIDs.size(), coorPoints.data());

  /*--- Reserve the memory for frontLeaves and frontLeavesNew,
        which are needed during the tree search. ---*/
  frontLeaves.reserve(200);
  frontLeavesNew.reserve(200);
}

void su2_adtPointsOnlyClass::DetermineNearestNode(const su2double *coor,
                                                  su2double       &dist,
                                                  unsigned long   &pointID,
                                                  int             &rankID) {

  AD_BEGIN_PASSIVE

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Initialize the nearest node to the central node of the     ---*/
  /*---         root leaf. Note that the distance is the distance squared  ---*/
  /*---         to avoid a sqrt.                                           ---*/
  /*--------------------------------------------------------------------------*/

  unsigned long kk = leaves[0].centralNodeID, minIndex;
  const su2double *coorTarget = coorPoints.data() + nDimADT*kk;

  pointID  = localPointIDs[kk];
  rankID   = ranksOfPoints[kk];
  minIndex = kk;
  dist = 0.0;
  for(unsigned short l=0; l<nDimADT; ++l) {
    const su2double ds = coor[l] - coorTarget[l];
    dist += ds*ds;
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
  for(;;) {

    /* Initialize the new front, i.e. the front for the next round, to empty. */
    frontLeavesNew.clear();

    /* Loop over the leaves of the current front. */
    for(unsigned long i=0; i<frontLeaves.size(); ++i) {

      /* Store the current leaf a bit easier in ll and loop over its children. */
      const unsigned long ll = frontLeaves[i];
      for(unsigned short mm=0; mm<2; ++mm) {

        /* Determine whether this child contains a node or a leaf
           of the next level of the ADT. */
        kk = leaves[ll].children[mm];
        if( leaves[ll].childrenAreTerminal[mm] ) {

          /*--- Child contains a node. Compute the distance squared to this node
                and store it if this distance squared is less than the currently
                stored value. ---*/
          coorTarget = coorPoints.data() + nDimADT*kk;
          su2double distTarget = 0;
          for(unsigned short l=0; l<nDimADT; ++l) {
            const su2double ds = coor[l] - coorTarget[l];
            distTarget += ds*ds;
          }

          if(distTarget < dist) {
            dist     = distTarget;
            pointID  = localPointIDs[kk];
            rankID   = ranksOfPoints[kk];
            minIndex = kk;
          }
        }
        else {

          /*--- Child contains a leaf. Determine the possible minimum distance
                squared to that leaf. ---*/
          su2double posDist = 0.0;
          for(unsigned short l=0; l<nDimADT; ++l) {
            su2double ds = 0.0;
            if(     coor[l] < leaves[kk].xMin[l]) ds = coor[l] - leaves[kk].xMin[l];
            else if(coor[l] > leaves[kk].xMax[l]) ds = coor[l] - leaves[kk].xMax[l];

            posDist += ds*ds;
          }

          /*--- Check if the possible minimum distance is less than the currently
                stored minimum distance. If so this leaf must be stored for the
                next round. In that case the distance squared to the central node is
                determined, which is used to update the currently stored value. ---*/
          if(posDist < dist) {
            frontLeavesNew.push_back(kk);

            const unsigned long jj = leaves[kk].centralNodeID;

            coorTarget = coorPoints.data() + nDimADT*jj;
            su2double distTarget = 0;
            for(unsigned short l=0; l<nDimADT; ++l) {
              const su2double ds = coor[l] - coorTarget[l];
              distTarget += ds*ds;
            }

            if(distTarget < dist) {
              dist     = distTarget;
              pointID  = localPointIDs[jj];
              rankID   = ranksOfPoints[jj];
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
    if(frontLeaves.size() == 0) break;
  }

  AD_END_PASSIVE

  /* Recompute the distance to get the correct dependency if we use AD */
  coorTarget = coorPoints.data() + nDimADT*minIndex;
  dist = 0.0;
  for(unsigned short l=0; l<nDimADT; ++l) {
    const su2double ds = coor[l] - coorTarget[l];
    dist += ds*ds;
  }


  /* At the moment the distance squared to the nearest node is stored.
     Take the sqrt to obtain the correct value. */
  dist = sqrt(dist);

}

su2_adtElemClass::su2_adtElemClass(unsigned short         val_nDim,
                                   vector<su2double>      &val_coor,
                                   vector<unsigned long>  &val_connElem,
                                   vector<unsigned short> &val_VTKElem,
                                   vector<unsigned short> &val_markerID,
                                   vector<unsigned long>  &val_elemID) {

  /* Copy the dimension of the problem into nDim. */
  nDim = val_nDim;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Gather the local grids on all ranks, such that the entire  ---*/
  /*---         grid to be searched is known on all ranks. This leads to   ---*/
  /*---         the most efficient algorithm in terms of performance, but  ---*/
  /*---         it could also lead to a memory bottleneck.                 ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Make a distinction between parallel and sequential mode. ---*/

#ifdef HAVE_MPI

  /*--- Parallel mode. The local grids are gathered on all ranks. For very
        large cases this could become a serious memory bottleneck and a
        parallel version may be needed.  ---*/

  /*--- First determine the number of points per rank and make them
        available to all ranks. ---*/
  int rank, size;
  SU2_MPI::Comm_rank(MPI_COMM_WORLD, &rank);
  SU2_MPI::Comm_size(MPI_COMM_WORLD, &size);

  vector<int> recvCounts(size), displs(size);
  int sizeLocal = (int) val_coor.size();

  SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1,
                     MPI_INT, MPI_COMM_WORLD);
  displs[0] = 0;
  for(int i=1; i<size; ++i) displs[i] = displs[i-1] + recvCounts[i-1];

  /*--- Correct the local connectivities with the offset of this rank,
        such that they get the correct values when all the points from
        all ranks are gathered. ---*/
  const unsigned long offsetRank = displs[rank]/nDim;
  for(unsigned long i=0; i<val_connElem.size(); ++i)
    val_connElem[i] += offsetRank;

  /*--- Gather the coordinates of the points on all ranks. ---*/
  int sizeGlobal = displs.back() + recvCounts.back();

  coorPoints.resize(sizeGlobal);
  SU2_MPI::Allgatherv(val_coor.data(), sizeLocal, MPI_DOUBLE, coorPoints.data(),
                      recvCounts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  /*--- Determine the number of elements per rank and make them
        available to all ranks. ---*/
  sizeLocal = (int) val_VTKElem.size();

  SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1,
                     MPI_INT, MPI_COMM_WORLD);
  displs[0] = 0;
  for(int i=1; i<size; ++i) displs[i] = displs[i-1] + recvCounts[i-1];

  /*--- Gather the element type, the possible markers of these elements
        and the local element ID's on all ranks. ---*/
  sizeGlobal = displs.back() + recvCounts.back();

  elemVTK_Type.resize(sizeGlobal);
  localMarkers.resize(sizeGlobal);
  localElemIDs.resize(sizeGlobal);

  SU2_MPI::Allgatherv(val_VTKElem.data(), sizeLocal, MPI_UNSIGNED_SHORT, elemVTK_Type.data(),
                      recvCounts.data(), displs.data(), MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);

  SU2_MPI::Allgatherv(val_markerID.data(), sizeLocal, MPI_UNSIGNED_SHORT, localMarkers.data(),
                      recvCounts.data(), displs.data(), MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);

  SU2_MPI::Allgatherv(val_elemID.data(), sizeLocal, MPI_UNSIGNED_LONG, localElemIDs.data(),
                      recvCounts.data(), displs.data(), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /*--- Create the content of ranksOfElems, which stores the original ranks
        where the elements come from. ---*/
  ranksOfElems.resize(sizeGlobal);

  for(int i=0; i<size; ++i) {
    for(int j=0; j<recvCounts[i]; ++j)
      ranksOfElems[displs[i]+j] = i;
  }

  /*--- Determine the size of the local connectivity per rank and make them
        available to all ranks. ---*/
  sizeLocal = (int) val_connElem.size();

  SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1,
                     MPI_INT, MPI_COMM_WORLD);
  displs[0] = 0;
  for(int i=1; i<size; ++i) displs[i] = displs[i-1] + recvCounts[i-1];

  /*--- Gather the element connectivity on all ranks. ---*/
  sizeGlobal = displs.back() + recvCounts.back();

  elemConns.resize(sizeGlobal);

  SU2_MPI::Allgatherv(val_connElem.data(), sizeLocal, MPI_UNSIGNED_LONG, elemConns.data(),
                      recvCounts.data(), displs.data(), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else

  /*--- Sequential mode. Copy the data from the arguments into the member
        variables and set the ranks to MASTER_NODE. ---*/
  coorPoints   = val_coor;
  elemConns    = val_connElem;
  elemVTK_Type = val_VTKElem;
  localMarkers = val_markerID;
  localElemIDs = val_elemID;

  ranksOfElems.assign(elemVTK_Type.size(), MASTER_NODE);

#endif

    /*--- Determine the values of the vector nDOFsPerElem, which contains the
        number of DOFs per element in cumulative storage format. ---*/
  const unsigned long nElem = elemVTK_Type.size();
  nDOFsPerElem.resize(nElem+1);

  nDOFsPerElem[0] = 0;
  for(unsigned long i=0; i<nElem; ++i) {

    unsigned short nDOFs = 0;
    switch( elemVTK_Type[i] ) {
      case LINE:          nDOFs = 2; break;
      case TRIANGLE:      nDOFs = 3; break;
      case QUADRILATERAL: nDOFs = 4; break;
      case TETRAHEDRON:   nDOFs = 4; break;
      case PYRAMID:       nDOFs = 5; break;
      case PRISM:         nDOFs = 6; break;
      case HEXAHEDRON:    nDOFs = 8; break;
    }

    nDOFsPerElem[i+1] = nDOFsPerElem[i] + nDOFs;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Create the bounding boxes of the elements. The coordinates ---*/
  /*---         of these bounding boxes can be interpreted as a point in a ---*/
  /*---         higher dimensional space. The ADT is built as a tree of    ---*/
  /*---         these points in this higher dimensional space.             ---*/
  /*--------------------------------------------------------------------------*/

  /* Allocate the memory for the bounding boxes of the elements. */
  BBoxCoor.resize(2*nDim*nElem);

  /*--- Loop over the elements to determine the minimum and maximum coordinates
        of the elements. These coordinates completely define the bounding box,
        which can be represented as a point in 2*nDim dimensional space. ---*/
  for(unsigned long i=0; i<nElem; ++i) {

    /* Set the pointers for the minimum and maximum coordinates for this
       bounding box. Initialize these coordinates to the first coordinate
       of the corresponding element. */
    su2double *BBMin = BBoxCoor.data() + 2*nDim*i;
    su2double *BBMax = BBMin + nDim;

    unsigned long    j  = nDOFsPerElem[i];
    const su2double *xP = coorPoints.data() + nDim*elemConns[j];

    for(unsigned short k=0; k<nDim; ++k)
      BBMin[k] = BBMax[k] = xP[k];

    /* Loop over the other vertices of the element to determine the
       minimum and maximum coordinates of the bounding box. */
    for(j=(nDOFsPerElem[i]+1); j<nDOFsPerElem[i+1]; ++j) {
      xP = coorPoints.data() + nDim*elemConns[j];

      for(unsigned short k=0; k<nDim; ++k) {
        BBMin[k] = min(BBMin[k], xP[k]);
        BBMax[k] = max(BBMax[k], xP[k]);
      }
    }

    /* Add a tolerance to the bounding box size, such that the overlap check in
       the tree traversal does not go wrong due to round off error. */
    for(unsigned short k=0; k<nDim; ++k) {
      const su2double lenScale = BBMax[k] - BBMin[k];
      const su2double tol      = max(1.e-25, 1.e-6*lenScale);

      BBMin[k] -= tol;
      BBMax[k] += tol;
    }
  }

  /* Build the ADT of the bounding boxes. */
  BuildADT(2*nDim, nElem, BBoxCoor.data());

  /*--- Reserve the memory for frontLeaves, frontLeavesNew and BBoxTargets,
        which are needed during the tree search. ---*/
  frontLeaves.reserve(200);
  frontLeavesNew.reserve(200);
  BBoxTargets.reserve(200);
}

void su2_adtElemClass::DetermineNearestElement(const su2double *coor,
                                               su2double       &dist,
                                               unsigned short  &markerID,
                                               unsigned long   &elemID,
                                               int             &rankID) {

  /*----------------------------------------------------------------------------*/
  /*--- Step 1: Initialize the distance (squared) to the quaranteed distance ---*/
  /*---         of the central bounding box of the root element.             ---*/
  /*----------------------------------------------------------------------------*/

  unsigned long kk = leaves[0].centralNodeID;
  const su2double *coorBBMin = BBoxCoor.data() + nDimADT*kk;
  const su2double *coorBBMax = coorBBMin + nDim;

  dist = 0.0;
  for(unsigned short k=0; k<nDim; ++k) {
    const su2double dsMin = fabs(coor[k] - coorBBMin[k]);
    const su2double dsMax = fabs(coor[k] - coorBBMax[k]);
    const su2double ds    = max(dsMin, dsMax);

    dist += ds*ds;
  }

  /*----------------------------------------------------------------------------*/
  /*--- Step 2: Traverse the tree and store the bounding boxes for which the ---*/
  /*---         possible minimum distance is less than the currently stored  ---*/
  /*---         guaranteed minimum distance. During this traversal the value ---*/
  /*---         of the currently stored minimum distance is adapted.         ---*/
  /*----------------------------------------------------------------------------*/

  /* Start at the root leaf of the ADT, i.e. initialize frontLeaves such that
     it only contains the root leaf. Make sure to wipe out any data from a
     previous search. */
  BBoxTargets.clear();
  frontLeaves.clear();
  frontLeaves.push_back(0);

  /* Infinite loop of the tree traversal. */
  for(;;) {

    /* Initialize the new front, i.e. the front for the next round, to empty. */
    frontLeavesNew.clear();

    /* Loop over the leaves of the current front. */
    for(unsigned long i=0; i<frontLeaves.size(); ++i) {

      /* Store the current leaf a bit easier in ll and loop over its children. */
      const unsigned long ll = frontLeaves[i];
      for(unsigned short mm=0; mm<2; ++mm) {

        /* Determine whether this child contains a node or a leaf
           of the next level of the ADT. */
        kk = leaves[ll].children[mm];
        if( leaves[ll].childrenAreTerminal[mm] ) {

          /*--- Child contains a bounding. Compute the possible distance squared
                to this bounding box. ---*/
          coorBBMin = BBoxCoor.data() + nDimADT*kk;
          coorBBMax = coorBBMin + nDim;

          su2double posDist2 = 0.0;
          for(unsigned short k=0; k<nDim; ++k) {
            su2double ds = 0.0;
            if(     coor[k] < coorBBMin[k]) ds = coor[k] - coorBBMin[k];
            else if(coor[k] > coorBBMax[k]) ds = coor[k] - coorBBMax[k];

            posDist2 += ds*ds;
          }

          /* Check if the possible minimum distance is less than the currently
             stored distance. If so, this bounding box is a candidate for the
             actual minimum distance and must be stored in BBoxTargets. */
          if(posDist2 < dist) {

            /*--- Compute the guaranteed minimum distance for this bounding box. ---*/
            su2double guarDist2 = 0.0;
            for(unsigned short k=0; k<nDim; ++k) { 
              const su2double dsMin = fabs(coor[k] - coorBBMin[k]);
              const su2double dsMax = fabs(coor[k] - coorBBMax[k]);
              const su2double ds    = max(dsMin, dsMax);

              guarDist2 += ds*ds;
            }

            /* Store this bounding box in BBoxTargets and update the currently
               stored value of the distance squared. */
            BBoxTargets.push_back(su2_BBoxTargetClass(kk, posDist2, guarDist2));
            dist = min(dist, guarDist2);
          }
        }
        else {

          /*--- Child contains a leaf. Determine the possible minimum distance
                squared to that leaf. ---*/
          coorBBMin = leaves[kk].xMin;
          coorBBMax = leaves[kk].xMax + nDim;

          su2double posDist2 = 0.0;
          for(unsigned short k=0; k<nDim; ++k) {
            su2double ds = 0.0;
            if(     coor[k] < coorBBMin[k]) ds = coor[k] - coorBBMin[k];
            else if(coor[k] > coorBBMax[k]) ds = coor[k] - coorBBMax[k];

            posDist2 += ds*ds;
          }

          /* Check if the possible minimum distance is less than the currently
             stored distance. If so this leaf must be stored for the next round. */
          if(posDist2 < dist) {
            frontLeavesNew.push_back(kk); 

            /*--- Determine the guaranteed minimum distance squared to the central
                  bounding box of this leaf and update the currently stored
                  minimum wall distance. ---*/
            kk = leaves[kk].centralNodeID;
            coorBBMin = BBoxCoor.data() + nDimADT*kk;
            coorBBMax = coorBBMin + nDim;

            su2double guarDist2 = 0.0;
            for(unsigned short k=0; k<nDim; ++k) {
              const su2double dsMin = fabs(coor[k] - coorBBMin[k]);
              const su2double dsMax = fabs(coor[k] - coorBBMax[k]);
              const su2double ds    = max(dsMin, dsMax);

              guarDist2 += ds*ds;
            }

            dist = min(dist, guarDist2);
          }
        }
      }
    }

    /*--- End of the loop over the current front. Copy the data from
          frontLeavesNew to frontLeaves for the next round. If the new front
          is empty the entire tree has been traversed and a break can be made
          from the infinite loop. ---*/
    frontLeaves = frontLeavesNew;
    if(frontLeaves.size() == 0) break;
  }

  /*----------------------------------------------------------------------------*/
  /*--- Step 3: Loop over the possible target bounding boxes and check if    ---*/
  /*---         the corresponding element minimizes the distance to the      ---*/
  /*---         given coordinate.                                            ---*/
  /*----------------------------------------------------------------------------*/

  /* Sort the bounding boxes in increasing order, such that the most likely
     candidates are checked first. */
  sort(BBoxTargets.begin(), BBoxTargets.end());

  /* Loop over the candidate bounding boxes. */
  for(unsigned long i=0; i<BBoxTargets.size(); ++i) {

    /* Break the loop if the possible minimum distance is larger than or equal
       to the currently stored value. In that case it does not make sense to
       check the remainder of the bounding boxes, as they are sorted in
       increasing order (based on the possible minimum distance. */
    if(BBoxTargets[i].possibleMinDist2 >= dist) break;

    /*--- Compute the distance squared to the element that corresponds to the
          current bounding box. If this distance is less than the current
          value, overwrite the return information of this function. ---*/
    const unsigned long ii = BBoxTargets[i].boundingBoxID;

    su2double dist2Elem;
    Dist2ToElement(ii, coor, dist2Elem);
    if(dist2Elem < dist) {
      dist     = dist2Elem;
      markerID = localMarkers[ii];
      elemID   = localElemIDs[ii];
      rankID   = ranksOfElems[ii];
    }
  }

  /* At the moment the square of the distance is stored in dist. Compute
     the correct value. */
  dist = sqrt(dist);
}

void su2_adtElemClass::Dist2ToElement(const unsigned long elemID,
                                      const su2double     *coor,
                                      su2double           &dist2Elem) {

  /*--- Make a distinction between the element types. ---*/
  switch( elemVTK_Type[elemID] ) {

    case LINE: {

      /*--- Element is a line. The number of space dimensions can be either
            1, 2 or 3. Store the indices where the coordinates of the vertices
            are stored in i0 and i1. ---*/
      unsigned long i0 = nDOFsPerElem[elemID];
      unsigned long i1 = i0 + 1;
      i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1];

      /*--- Call the function Dist2ToLine to do the actual work. ---*/
      Dist2ToLine(i0, i1, coor, dist2Elem);

      break;
    }

    /*------------------------------------------------------------------------*/

    case TRIANGLE: {

      /*--- Element is a triangle. The number of space dimensions can be either
            2 or 3. Store the indices where the coordinates of the vertices
            are stored in i0, i1 and i2. ---*/
      unsigned long i0 = nDOFsPerElem[elemID];
      unsigned long i1 = i0 + 1, i2 = i0 + 2;
      i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1]; i2 = nDim*elemConns[i2];

      /*--- Call the function Dist2ToTriangle to compute the distance to the
            triangle if the projection is inside the triangle. In that case the
            function returns true. If the projection is not inside the triangle,
            false is returned and the distance to each of the lines of the
            triangle is computed and the minimum is taken. ---*/
      su2double r, s;
      if( !Dist2ToTriangle(i0, i1, i2, coor, dist2Elem, r, s) ) {
        Dist2ToLine(i0, i1, coor, dist2Elem);

        su2double dist2Line;
        Dist2ToLine(i1, i2, coor, dist2Line); dist2Elem = min(dist2Elem, dist2Line);
        Dist2ToLine(i2, i0, coor, dist2Line); dist2Elem = min(dist2Elem, dist2Line);
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case QUADRILATERAL: {

      /*--- Element is a quadrilateral. The number of space dimensions can be
            either 2 or 3. Store the indices where the coordinates of the
            vertices are stored in i0, i1, i2 and i3. ---*/
      unsigned long i0 = nDOFsPerElem[elemID];
      unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0 + 3;
      i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1];
      i2 = nDim*elemConns[i2]; i3 = nDim*elemConns[i3];

      /*--- As the mapping to the standard quadrilateral contains a nonlinear
            term, the computation of the minimum distance involves a Newton
            algorithm. Consequently this computation is quite costly and for
            wrapped quadrilaterals not very stable when the projection is far
            outside the the element. Therefore in a predictor step the
            quadrilateral is divived into two triangles, for which the mapping
            to the standard element is linear, and it is checked if the
            projection is inside one of this triangles. ---*/
      su2double r, s;
      if( Dist2ToTriangle(i0, i1, i3, coor, dist2Elem, r, s) ) {

        /* Projection is inside the triangle i0, i1, i3. Compute the true
           distance to the quadrilateral, where the parametric coordinates
           r and s are used as initial guess. */
        Dist2ToQuadrilateral(i0, i1, i2, i3, coor, r, s, dist2Elem);
      }
      else if( Dist2ToTriangle(i2, i3, i1, coor, dist2Elem, r, s) ) {

        /* Projection is inside the triangle i2, i3, i1. Compute the true
           distance to the quadrilateral, where the parametric coordinates
           r and s are used as initial guess. Note that r and s must be
           converted to the parameters of the quadrilateral. */
        r = 1.0 - r;
        s = 1.0 - s;

        Dist2ToQuadrilateral(i0, i1, i2, i3, coor, r, s, dist2Elem);
      }
      else {

        /* The projection is outside the quadrilatral. Hence it suffices
           to check the distance to the surrounding lines of the quad. */
        Dist2ToLine(i0, i1, coor, dist2Elem);

        su2double dist2Line;
        Dist2ToLine(i1, i2, coor, dist2Line); dist2Elem = min(dist2Elem, dist2Line);
        Dist2ToLine(i2, i3, coor, dist2Line); dist2Elem = min(dist2Elem, dist2Line);
        Dist2ToLine(i3, i0, coor, dist2Line); dist2Elem = min(dist2Elem, dist2Line);
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case TETRAHEDRON:
    case PYRAMID:
    case PRISM:
    case HEXAHEDRON: {
      SU2_MPI::Error("3D elements not implemented yet", CURRENT_FUNCTION);
    }
  }
}

void su2_adtElemClass::Dist2ToLine(const unsigned long i0, 
                                   const unsigned long i1, 
                                   const su2double     *coor, 
                                   su2double           &dist2Line) {

  /*--- The line is parametrized by X = X0 + (r+1)*(X1-X0)/2, -1 <= r <= 1.
        As a consequence the minimum distance is found where the expression
        |V0 - r*V1| has a minimum, where the vectors V0 and V1 are defined
        as: V0 = coor - (X1+X0)/2, V1 = (X1-X0)/2. First construct the
        vectors V0 and V1. ---*/
  su2double V0[3], V1[3];
  for(unsigned short k=0; k<nDim; ++k) {
    V0[k] = coor[k] - 0.5*(coorPoints[i1+k] + coorPoints[i0+k]);
    V1[k] =           0.5*(coorPoints[i1+k] - coorPoints[i0+k]);
  }

  /*--- Determine the value of r, for which the minimum occurs. This is the
        ratio of the dot product of V0.V1 and V1.V1. Make sure that r is
        limited between -1 and 1 to make sure that the projection is inside
        the line element. ---*/
  su2double dotV0V1 = 0.0, dotV1V1 = 0.0;
  for(unsigned short k=0; k<nDim; ++k) {
    dotV0V1 += V0[k]*V1[k];
    dotV1V1 += V1[k]*V1[k];
  }
  su2double r = dotV0V1/dotV1V1;
  r = max(-1.0,min(1.0,r));

  /*--- Determine the minimum distance squared. ---*/
  dist2Line = 0.0;
  for(unsigned short k=0; k<nDim; ++k) {
    const su2double ds = V0[k] - r*V1[k];
    dist2Line += ds*ds;
  }
}

bool su2_adtElemClass::Dist2ToTriangle(const unsigned long i0,
                                       const unsigned long i1,
                                       const unsigned long i2,
                                       const su2double     *coor,
                                       su2double           &dist2Tria,
                                       su2double           &r,
                                       su2double           &s) {

  /*--- The triangle is parametrized by X = X0 + (r+1)*(X1-X0)/2 + (s+1)*(X2-X0)/2,
        r, s >= -1, r+s <= 0. As a consequence the minimum distance is found where
        the expression |V0 - r*V1 - s*V2| has a minimum, where the vectors V0, V1
        and V2 are defined as: V0 = coor - (X1+X2)/2, V1 = (X1-X0)/2, V2 = (X2-X0)/2.
        First construct the vectors V0, V1 and V2. ---*/
  su2double V0[3], V1[3], V2[3];
  for(unsigned short k=0; k<nDim; ++k) {
    V0[k] = coor[k] - 0.5*(coorPoints[i1+k] + coorPoints[i2+k]);
    V1[k] =           0.5*(coorPoints[i1+k] - coorPoints[i0+k]);
    V2[k] =           0.5*(coorPoints[i2+k] - coorPoints[i0+k]);
  }

  /*--- For the values of r and s where the minimum occurs the dot products V0.V1,
        V0.V2, V1.V1, V1.V2 and V2.V2 are needed. Compute these values. ---*/
  su2double dotV0V1 = 0.0, dotV0V2 = 0, dotV1V1 = 0.0, dotV1V2 = 0.0, dotV2V2 = 0.0;
  for(unsigned short k=0; k<nDim; ++k) {
    dotV0V1 += V0[k]*V1[k];
    dotV0V2 += V0[k]*V2[k];
    dotV1V1 += V1[k]*V1[k];
    dotV1V2 += V1[k]*V2[k];
    dotV2V2 += V2[k]*V2[k];
  }

  /*--- Solve the linear system for r and s. ---*/
  const su2double detInv = 1.0/(dotV1V1*dotV2V2 - dotV1V2*dotV1V2);

  r = detInv*(dotV0V1*dotV2V2 - dotV0V2*dotV1V2);
  s = detInv*(dotV0V2*dotV1V1 - dotV0V1*dotV1V2);

  /*--- Check if the projection is inside the triangle. ---*/
  if((r >= paramLowerBound) && (s >= paramLowerBound) && ((r+s) <= tolInsideElem)) {

    /*--- The projection of the coordinate is inside the triangle. Compute the
          minimum distance squared and return true. ---*/
    dist2Tria = 0.0;
    for(unsigned short k=0; k<nDim; ++k) {
      const su2double ds = V0[k] - r*V1[k] - s*V2[k];
      dist2Tria += ds*ds;
    }

    return true;
  }

  /*--- The projection of the coordinate is outside the triangle.
        Return false. ---*/
  return false;
}

void su2_adtElemClass::Dist2ToQuadrilateral(const unsigned long i0,
                                            const unsigned long i1,
                                            const unsigned long i2,
                                            const unsigned long i3,
                                            const su2double     *coor,
                                            su2double           &r,
                                            su2double           &s,
                                            su2double           &dist2Quad) {

  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton  = 1.e-10;

  /*--- The quadrilateral is parametrized by
        X = {(1-r)*(1-s)*X0 + (1+r)*(1-s)*X1 + (1+r)*(1+s)*X2 + (1-r)*(1+s)*X3}/4,
        -1 <= r,s <= 1. As a consequence the minimum distance is found where the
        expression |V0 - r*V1 - s*V2 - r*s*V3| has a minimum, where the vectors
        V0, V1, V2 and V3 are defined as: V0 = coor - (X0+X1+X2+X3)/4,
        V1 = (X1+X2-X0-X3)/4, V2 = (X2+X3-X0+X1), V3 = (X0+X2-X1-X3)/4. First
        construct the vectors V0, V1, V2 and V3. ---*/
  su2double V0[3], V1[3], V2[3], V3[3];
  for(unsigned short k=0; k<nDim; ++k) {
    V0[k] = coor[k] - 0.25*(coorPoints[i0+k] + coorPoints[i1+k]
          +                 coorPoints[i2+k] + coorPoints[i3+k]);
    V1[k] =           0.25*(coorPoints[i1+k] + coorPoints[i2+k]
          -                 coorPoints[i0+k] - coorPoints[i3+k]);
    V2[k] =           0.25*(coorPoints[i2+k] + coorPoints[i3+k]
          -                 coorPoints[i0+k] - coorPoints[i1+k]);
    V3[k] =           0.25*(coorPoints[i0+k] + coorPoints[i2+k]
          -                 coorPoints[i1+k] - coorPoints[i3+k]);
  }

  /*--- Newton algorithm to solve the nonlinear equations
        (V0 - r*V1 - s*V2 - r*s*V3).(-V1 - s*V3) = 0
        (V0 - r*V1 - s*V2 - r*s*V3).(-V2 - r*V3) = 0.
        These equations are the gradients of the distance function squared w.r.t.
        the parametric coordinates r and s. ---*/
  unsigned short itCount;
  for(itCount=0; itCount<maxIt; ++itCount) {

    /* Compute the vectors needed in the equations to be solved. */
    su2double distVec[3], distVecDr[3], distVecDs[3];
    for(unsigned short k=0; k<nDim; ++k) {
      distVec[k]   =  V0[k] - r*V1[k] - s*V2[k] - r*s*V3[k];
      distVecDr[k] = -V1[k] - s*V3[k];
      distVecDs[k] = -V2[k] - r*V3[k];
    }

    /*--- The functional to be minimized is the L2 norm (squared) of the
          distance vector. Compute the gradients of this functional as well
          as the elements of the Hessian matrix. ---*/
    su2double dfdr = 0.0, dfds = 0.0, H00 = 0.0, H11 = 0.0, H01 = 0.0;
    for(unsigned short k=0; k<nDim; ++k) {
      dfdr += distVec[k]  *distVecDr[k];
      dfds += distVec[k]  *distVecDs[k];
      H00  += distVecDr[k]*distVecDr[k];
      H11  += distVecDs[k]*distVecDs[k];
      H01  += distVecDr[k]*distVecDs[k] - distVec[k]*V3[k];
    }
       
    /*--- Make sure that the Hessian is positive definite. This is accomplished
          by checking the smallest eigenvalue. If this is smaller than a certain
          cutoff, the Hessian is modified. ---*/
    const su2double DD    = sqrt((H00-H11)*(H00-H11) + 4.0*H01*H01);
    const su2double lam1  = 0.5*(H00+H11 +DD);
    su2double       lam2  = 0.5*(H00+H11 -DD);
    const su2double thres = 1.0e-8*lam1;

    if(lam2 < thres) {
      lam2 = thres;
      const su2double tmp = (H00 - H11)*(lam1 - lam2)/DD;

      H00  = 0.5*(lam1 + lam2 + tmp);
      H11  = 0.5*(lam1 + lam2 - tmp);
      H01  = H01*(lam1 - lam2)/DD;
    }

    /*--- Compute the updates for r and s. ---*/
    const su2double detInv = 1.0/(H00*H11 - H01*H01);
    const su2double dr     = (dfdr*H11 - dfds*H01)*detInv;
    const su2double ds     = (dfds*H00 - dfdr*H01)*detInv;

    /*--- Compute the new values of r and s. ---*/
    r -= dr;
    s -= ds;

    /*--- Check for convergence. ---*/
    if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton) break;
  }

  /*--- Terminate if the optimization process did not converge. ---*/
  if(itCount == maxIt)
    SU2_MPI::Error("Newton did not converge", CURRENT_FUNCTION);

  /*--- Check if the projection is inside the quadrilateral. If not, terminate. ---*/
  if(r < paramLowerBound || r > paramUpperBound ||
     s < paramLowerBound || s > paramUpperBound)
    SU2_MPI::Error("Projection not inside the quadrilateral.", CURRENT_FUNCTION);

  /*--- Determine the minimum distance squared. ---*/
  dist2Quad = 0.0;
  for(unsigned short k=0; k<nDim; ++k) {
    const su2double ds = V0[k] - r*V1[k] - s*V2[k] - r*s*V3[k];
    dist2Quad += ds*ds;
  }
}
