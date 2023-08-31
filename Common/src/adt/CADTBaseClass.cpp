/*!
 * \file CADTBaseClass.cpp
 * \brief Base class for storing an ADT in an arbitrary number of dimensions.
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

#include "../../include/adt/CADTBaseClass.hpp"
#include "../../include/adt/CADTPointsOnlyClass.hpp"
#include "../../include/adt/CADTComparePointClass.hpp"

#include <algorithm>

void CADTBaseClass::BuildADT(unsigned short nDim, unsigned long nPoints, const su2double* coor) {
  /*---  Determine the number of leaves. It can be proved that nLeaves equals
         nPoints-1 for an optimally balanced tree. Take the exceptional case of
         nPoints == 1 into account and return if the tree is empty. ---*/
  nDimADT = nDim;
  isEmpty = false;
  nLeaves = nPoints - 1;
  if (nPoints <= 1) ++nLeaves;
  if (nLeaves == 0) {
    isEmpty = true;
    return;
  }

  /*--- Allocate the memory for the leaves of the ADT and the minimum and
        maximum coordinates of the leaves. Note that these coordinates are
        stored in one vector, rather than that memory is allocated for the
        individual leaves. ---*/
  leaves.resize(nLeaves);
  coorMinLeaves.resize(nDim * nLeaves);
  coorMaxLeaves.resize(nDim * nLeaves);

  /*--- Define the vectors, which control the subdivision of the leaves. ---*/
  unsigned long nn = (nPoints + 1) / 2;
  vector<unsigned long> pointIDs(nPoints), pointIDsNew(nPoints);
  vector<unsigned long> nPointIDs(nn + 1), nPointIDsNew(nn + 1);
  vector<unsigned long> curLeaf(nn), curLeafNew(nn);

  /*--------------------------------------------------------------------------*/
  /*---                 Building of the actual ADT                         ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Initialize the arrays pointIDs, nPointIDs and curLeaf such that all
        points belong to the root leaf. Also set the counters nLeavesToDivide
        and nLeavesTot.  ---*/
  nPointIDs[0] = 0;
  nPointIDs[1] = nPoints;
  curLeaf[0] = 0;

  for (unsigned long i = 0; i < nPoints; ++i) pointIDs[i] = i;

  unsigned long nLeavesToDivide = 1, nLeavesTot = 1;

  /*--- Loop to subdivide the leaves. The division is such that the ADT is
        optimally balanced.  ---*/
  for (;;) {
    /* Criterion to exit the loop. */
    if (nLeavesToDivide == 0) break;

    /* Initializations for the next round of subdivisions. */
    unsigned long nLeavesToDivideNew = 0;
    nPointIDsNew[0] = 0;

    /*--- Loop over the current number of leaves to be divided. ---*/
    for (unsigned long i = 0; i < nLeavesToDivide; ++i) {
      /* Store the number of points present in the leaf in nn and the
         current leaf number in mm. */
      nn = nPointIDs[i + 1] - nPointIDs[i];
      unsigned long mm = curLeaf[i];

      /*--- Set the pointers for the coordinates of the leaf to the correct
            locations in the vectors coorMinLeaves and coorMaxLeaves and
            determine the bounding box coordinates of this leaf. ---*/
      leaves[mm].xMin = coorMinLeaves.data() + nDim * mm;
      leaves[mm].xMax = coorMaxLeaves.data() + nDim * mm;

      unsigned long ll = nDim * pointIDs[nPointIDs[i]];
      for (unsigned short l = 0; l < nDim; ++l) leaves[mm].xMin[l] = leaves[mm].xMax[l] = coor[ll + l];

      for (unsigned long j = (nPointIDs[i] + 1); j < nPointIDs[i + 1]; ++j) {
        ll = nDim * pointIDs[j];
        for (unsigned short l = 0; l < nDim; ++l) {
          leaves[mm].xMin[l] = min(leaves[mm].xMin[l], coor[ll + l]);
          leaves[mm].xMax[l] = max(leaves[mm].xMax[l], coor[ll + l]);
        }
      }

      /*--- Determine the split direction for this leaf. The splitting is done
            in such a way that isotropy is reached as quickly as possible.
            Hence the split direction is the direction of the largest dimension
            of the leaf. ---*/
      unsigned short splitDir = 0;
      su2double distMax = -1.0;
      for (unsigned short l = 0; l < nDim; ++l) {
        const su2double dist = leaves[mm].xMax[l] - leaves[mm].xMin[l];
        if (dist > distMax) {
          distMax = dist;
          splitDir = l;
        }
      }

      /* Sort the points of the current leaf in increasing order. The sorting
         is based on the coordinate in the split direction, for which the
         functor CADTComparePointClass is used. */
      sort(pointIDs.data() + nPointIDs[i], pointIDs.data() + nPointIDs[i + 1],
           CADTComparePointClass(coor, splitDir, nDim));

      /* Determine the index of the node, which is approximately central
         in this leave. */
      leaves[mm].centralNodeID = pointIDs[nPointIDs[i] + nn / 2];

      /*--- Determine the situation of the leaf. It is either a terminal leaf
            or a leaf that must be subdivided. ---*/
      if (nn <= 2) {
        /* Terminal leaf. Store the ID's of the points as children and
           indicate that the children are terminal. */
        leaves[mm].children[0] = pointIDs[nPointIDs[i]];
        leaves[mm].children[1] = pointIDs[nPointIDs[i + 1] - 1];

        leaves[mm].childrenAreTerminal[0] = true;
        leaves[mm].childrenAreTerminal[1] = true;
      } else {
        /* The leaf must be divided. Determine the number of points in the
           left leaf. This number is at least 2. The actual number stored in kk
           is this number plus an offset. Also initialize the counter nfl, which
           is used to store the bounding boxes in the arrays for the new round. */
        unsigned long kk = (nn + 1) / 2 + nPointIDs[i];
        unsigned long nfl = nPointIDsNew[nLeavesToDivideNew];

        /* Copy the ID's of the left points into pointIDsNew. Also update the
           corresponding entry in nPointIDsNew. */
        for (unsigned long k = nPointIDs[i]; k < kk; ++k) pointIDsNew[nfl++] = pointIDs[k];

        nPointIDsNew[nLeavesToDivideNew + 1] = nfl;

        /* Store the leaf info in the tree and in the leafs for next round and
           update nLeavesToDivideNew and nLeavesTot. */
        leaves[mm].children[0] = nLeavesTot;
        leaves[mm].childrenAreTerminal[0] = false;

        curLeafNew[nLeavesToDivideNew] = nLeavesTot;
        ++nLeavesToDivideNew;
        ++nLeavesTot;

        /*--- The right leaf will only be created if it has more than one point
              in it, i.e. if the original leaf has more than three points.
              If the new leaf only has one point in it, it is not created.
              Instead, the point is stored in the current leaf. ---*/
        if (nn == 3) {
          /* Only three points present in the current leaf. The right leaf is
             not created and the last point is stored as the second child of
             the current leaf. */
          leaves[mm].children[1] = pointIDs[nPointIDs[i + 1] - 1];
          leaves[mm].childrenAreTerminal[1] = true;
        } else {
          /* More than 3 points are present and thus the right leaf is created.
             Copy the ID's from pointIDs into pointIDsNew and update the
             counters for the new round. */
          unsigned long nfr = nPointIDsNew[nLeavesToDivideNew];

          for (unsigned long k = kk; k < nPointIDs[i + 1]; ++k) pointIDsNew[nfr++] = pointIDs[k];

          nPointIDsNew[nLeavesToDivideNew + 1] = nfr;

          /* Store the leaf info in the tree and in the leaves for next round and
             update nLeavesToDivideNew and nLeavesTot. */
          leaves[mm].children[1] = nLeavesTot;
          leaves[mm].childrenAreTerminal[1] = false;

          curLeafNew[nLeavesToDivideNew] = nLeavesTot;
          ++nLeavesToDivideNew;
          ++nLeavesTot;
        }
      }
    }

    /* Set the data for the next round. */
    nLeavesToDivide = nLeavesToDivideNew;

    for (unsigned long i = 0; i <= nLeavesToDivide; ++i) nPointIDs[i] = nPointIDsNew[i];
    for (unsigned long i = 0; i < nLeavesToDivide; ++i) curLeaf[i] = curLeafNew[i];
    for (unsigned long i = 0; i < nPointIDs[nLeavesToDivide]; ++i) pointIDs[i] = pointIDsNew[i];
  }
}
