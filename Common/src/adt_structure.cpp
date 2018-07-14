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

#include <iomanip>

su2_adtComparePointClass::su2_adtComparePointClass(const su2double      *coor,
                                                   const unsigned short splitDir,
                                                   const unsigned short nDimADT)
  : pointCoor(coor),
    splitDirection(splitDir),
    nDim(nDimADT) {}

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

su2_adtPointsOnlyClass::su2_adtPointsOnlyClass(unsigned short      nDim,
                                               unsigned long       nPoints,
                                               const su2double     *coor,
                                               const unsigned long *pointID) {

  /*--- Make a distinction between parallel and sequential mode. ---*/

#ifdef HAVE_MPI

  /*--- Parallel mode. All points are gathered on all ranks. First determine the
        number of points per rank and store them in such a way that the info can
        be used directly in Allgatherv. For now, we will use the regular 
        Allgather until we add Allgatherv to the SU2_MPI wrapper. ---*/
  
  int rank, iProcessor, nProcessor;
  unsigned long  iVertex, nBuffer;
  unsigned short iDim;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  unsigned long nLocalVertex = nPoints, nGlobalVertex = 0, MaxLocalVertex = 0;
  
  unsigned long *Buffer_Send_nVertex    = new unsigned long [1];
  unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];

  Buffer_Send_nVertex[0] = nLocalVertex;
  
  SU2_MPI::Allreduce(&nLocalVertex, &nGlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /*--- Gather the local pointID's and the ranks of the nodes on all ranks. ---*/
  
  unsigned long *Buffer_Send = new unsigned long[MaxLocalVertex];
  unsigned long *Buffer_Recv = new unsigned long[nProcessor*MaxLocalVertex];
  
  for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
    Buffer_Send[iVertex] = pointID[iVertex];
  }
  
  SU2_MPI::Allgather(Buffer_Send, MaxLocalVertex, MPI_UNSIGNED_LONG, Buffer_Recv, MaxLocalVertex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  /*--- Unpack the buffer into the local point ID vector. ---*/
  
  localPointIDs.resize(nGlobalVertex);
  
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
    for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++)
      localPointIDs.push_back( Buffer_Recv[iProcessor*MaxLocalVertex + iVertex] );

  /*--- Now gather the ranks for all points ---*/
  
  for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
    Buffer_Send[iVertex] = (unsigned long) rank;
  }
  
  SU2_MPI::Allgather(Buffer_Send, MaxLocalVertex, MPI_UNSIGNED_LONG, Buffer_Recv, MaxLocalVertex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  /*--- Unpack the ranks into the vector and delete buffer memory. ---*/
  
  ranksOfPoints.resize(nGlobalVertex);

  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
    for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++)
      ranksOfPoints.push_back( Buffer_Recv[iProcessor*MaxLocalVertex + iVertex] );
  
  delete [] Buffer_Send;  delete [] Buffer_Recv;
  
  /*--- Gather the coordinates of the points on all ranks. ---*/
  
  su2double *Buffer_Send_Coord = new su2double [MaxLocalVertex*nDim];
  su2double *Buffer_Recv_Coord = new su2double [nProcessor*MaxLocalVertex*nDim];
  
  nBuffer = MaxLocalVertex*nDim;
  
  for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
    for (iDim = 0; iDim < nDim; iDim++)
    Buffer_Send_Coord[iVertex*nDim + iDim] = coor[iVertex*nDim + iDim];
  }
  
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Recv_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);
  
  /*--- Unpack the coordinates into the vector and delete buffer memory. ---*/
  
  coorPoints.resize(nDim*nGlobalVertex);
  
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
    for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++)
      for (iDim = 0; iDim < nDim; iDim++)
      coorPoints.push_back( Buffer_Recv_Coord[iProcessor*MaxLocalVertex*nDim + iVertex*nDim + iDim] );
  
  delete [] Buffer_Send_Coord;   delete [] Buffer_Recv_Coord;
  delete [] Buffer_Send_nVertex; delete [] Buffer_Receive_nVertex;
  
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

void su2_adtPointsOnlyClass::Determine_N_NearestNodes(int N,
		const su2double *coor, su2double * dist, unsigned long* pointID,
		int* rankID) {

	/*--------------------------------------------------------------------------*/
	/*--- Step 1: Initialize the nearest node to the central node of the     ---*/
	/*---         root leaf. Note that the distance is the distance squared  ---*/
	/*---         to avoid a sqrt.                                           ---*/
	/*--------------------------------------------------------------------------*/

	unsigned long kk = leaves[0].centralNodeID;
	const su2double *coorTarget = coorPoints.data() + nDimADT * kk;

	for (int i = 0; i < N; i++) {
		pointID[0] = 0;
		rankID[0] = HUGE_VAL;
		dist[i] = HUGE_VAL;
	}

	//Initialize the first distance
	pointID[0] = localPointIDs[kk];
	rankID[0] = ranksOfPoints[kk];
	dist[0] = 0;
	//Normalized distance per component (i.e. component should never be 0)
	for (unsigned short l = 0; l < nDimADT; ++l) {
		const su2double ds = (coor[l] - coorTarget[l]) / coor[l];
		dist[0] += ds * ds;
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
					 and store it if this distance squared is less than any of the currently
					 stored values. ---*/
					coorTarget = coorPoints.data() + nDimADT * kk;
					su2double distTarget = 0;
					for (unsigned short l = 0; l < nDimADT; ++l) {
						const su2double ds = (coor[l] - coorTarget[l]) / coor[l];
						distTarget += ds * ds;
					}
					//Generalization of the algorithm to handle as many points as desired.
					int dist_i = 0;
					while (dist_i < N) {
						if (distTarget <= dist[dist_i]){
							if(pointID[dist_i] != localPointIDs[kk]){
							for (int dist_j = N - 1; dist_j > dist_i; dist_j--) {
								dist[dist_j] = dist[dist_j - 1];
								pointID[dist_j] = pointID[dist_j - 1];
								rankID[dist_j] = rankID[dist_j - 1];
							}
							dist[dist_i] = distTarget;
							pointID[dist_i] = localPointIDs[kk];
							rankID[dist_i] = ranksOfPoints[kk];
							}
							dist_i = N + 1; //break the loop
						}
						dist_i++;
					}
				} else {

					/*--- Child contains a leaf. Determine the possible minimum distance
					 squared to that leaf. ---*/
					su2double posDist = 0.0;
					for (unsigned short l = 0; l < nDimADT; ++l) {
						su2double ds = 0.0;
						if (coor[l] < leaves[kk].xMin[l])
							ds = (coor[l] - leaves[kk].xMin[l]) / coor[l];
						else if (coor[l] > leaves[kk].xMax[l])
							ds = (coor[l] - leaves[kk].xMax[l]) / coor[l];

						posDist += ds * ds;
					}

					/*--- Check if the possible minimum distance is less than the currently
					 worst stored minimum distance. If so this leaf must be stored for the
					 next round. In that case the distance squared to the central node is
					 determined, which is used to update the currently stored value. ---*/
					if (posDist < dist[N - 1]) {
						frontLeavesNew.push_back(kk);

						const unsigned long jj = leaves[kk].centralNodeID;

						coorTarget = coorPoints.data() + nDimADT * jj;
						su2double distTarget = 0;
						for (unsigned short l = 0; l < nDimADT; ++l) {
							const su2double ds = (coor[l] - coorTarget[l]) / coor[l];
							distTarget += ds * ds;
						}

						//Generalization of the algorithm to handle as many points as desired.
						int dist_i = 0;
						while (dist_i < N) {
							if (distTarget <= dist[dist_i]) {
								if(pointID[dist_i] != localPointIDs[jj]){
								for (int dist_j = N - 1; dist_j > dist_i; dist_j--) {
									dist[dist_j] = dist[dist_j - 1];
									pointID[dist_j] = pointID[dist_j - 1];
									rankID[dist_j] = rankID[dist_j - 1];
								}
								dist[dist_i] = distTarget;
								pointID[dist_i] = localPointIDs[jj];
								rankID[dist_i] = ranksOfPoints[jj];
								}
								dist_i = N + 1; //break the loop
							}
							dist_i++;
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
		if (frontLeaves.size() == 0)
			break;
	}

	/* At the moment the distance squared to the nearest node is stored.
	 Take the sqrt to obtain the correct value. */
	for (int i; i < N; i++) {
		dist[i] = sqrt(dist[i]);
	}

}
