/*!
 * \file CADTElemClass.cpp
 * \brief Class for storing an ADT of (linear) elements in an arbitrary number of dimensions.
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

#include "../../include/adt/CADTElemClass.hpp"
#include "../../include/parallelization/mpi_structure.hpp"
#include "../../include/option_structure.hpp"

/* Define the tolerance to decide whether or not a point is inside an element. */
const su2double tolInsideElem = 1.e-10;
const su2double paramLowerBound = -1.0 - tolInsideElem;
const su2double paramUpperBound = 1.0 + tolInsideElem;

CADTElemClass::CADTElemClass(unsigned short val_nDim, vector<su2double>& val_coor, vector<unsigned long>& val_connElem,
                             vector<unsigned short>& val_VTKElem, vector<unsigned short>& val_markerID,
                             vector<unsigned long>& val_elemID, const bool globalTree) {
  /* Copy the dimension of the problem into nDim. */
  nDim = val_nDim;

  /* Allocate some thread-safe working variables if required. */
#ifdef HAVE_OMP
  BBoxTargets.resize(omp_get_max_threads());
  FrontLeaves.resize(omp_get_max_threads());
  FrontLeavesNew.resize(omp_get_max_threads());
#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: If a global tree must be built, gather the local grids on  ---*/
  /*---         all ranks, such that the entire grid to be searched is     ---*/
  /*---         known on all ranks. This leads to the most efficient       ---*/
  /*---         algorithm in terms of performance, but it could also lead  ---*/
  /*---         to a memory bottleneck. If a local tree must be built,     ---*/
  /*---         this gathering is not needed.                              ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Make a distinction between parallel and sequential mode. ---*/

#ifdef HAVE_MPI

  /* Parallel mode. Check whether a global or a local tree must be built. */
  if (globalTree) {
    /*--- The local grids are gathered on all ranks. For very large cases this
          could become a serious memory bottleneck and a parallel version may
          be needed.  ---*/

    /*--- First determine the number of points per rank and make them
          available to all ranks. ---*/
    int rank, size;
    SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
    SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);

    vector<int> recvCounts(size), displs(size);
    int sizeLocal = (int)val_coor.size();

    SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, SU2_MPI::GetComm());
    displs[0] = 0;
    for (int i = 1; i < size; ++i) displs[i] = displs[i - 1] + recvCounts[i - 1];

    /*--- Correct the local connectivities with the offset of this rank,
          such that they get the correct values when all the points from
          all ranks are gathered. ---*/
    const unsigned long offsetRank = displs[rank] / nDim;
    for (unsigned long i = 0; i < val_connElem.size(); ++i) val_connElem[i] += offsetRank;

    /*--- Gather the coordinates of the points on all ranks. ---*/
    int sizeGlobal = displs.back() + recvCounts.back();

    coorPoints.resize(sizeGlobal);
    SU2_MPI::Allgatherv(val_coor.data(), sizeLocal, MPI_DOUBLE, coorPoints.data(), recvCounts.data(), displs.data(),
                        MPI_DOUBLE, SU2_MPI::GetComm());

    /*--- Determine the number of elements per rank and make them
          available to all ranks. ---*/
    sizeLocal = (int)val_VTKElem.size();

    SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, SU2_MPI::GetComm());
    displs[0] = 0;
    for (int i = 1; i < size; ++i) displs[i] = displs[i - 1] + recvCounts[i - 1];

    /*--- Gather the element type, the possible markers of these elements
          and the local element ID's on all ranks. ---*/
    sizeGlobal = displs.back() + recvCounts.back();

    elemVTK_Type.resize(sizeGlobal);
    localMarkers.resize(sizeGlobal);
    localElemIDs.resize(sizeGlobal);

    SU2_MPI::Allgatherv(val_VTKElem.data(), sizeLocal, MPI_UNSIGNED_SHORT, elemVTK_Type.data(), recvCounts.data(),
                        displs.data(), MPI_UNSIGNED_SHORT, SU2_MPI::GetComm());

    SU2_MPI::Allgatherv(val_markerID.data(), sizeLocal, MPI_UNSIGNED_SHORT, localMarkers.data(), recvCounts.data(),
                        displs.data(), MPI_UNSIGNED_SHORT, SU2_MPI::GetComm());

    SU2_MPI::Allgatherv(val_elemID.data(), sizeLocal, MPI_UNSIGNED_LONG, localElemIDs.data(), recvCounts.data(),
                        displs.data(), MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

    /*--- Create the content of ranksOfElems, which stores the original ranks
          where the elements come from. ---*/
    ranksOfElems.resize(sizeGlobal);

    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < recvCounts[i]; ++j) ranksOfElems[displs[i] + j] = i;
    }

    /*--- Determine the size of the local connectivity per rank and make them
          available to all ranks. ---*/
    sizeLocal = (int)val_connElem.size();

    SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, SU2_MPI::GetComm());
    displs[0] = 0;
    for (int i = 1; i < size; ++i) displs[i] = displs[i - 1] + recvCounts[i - 1];

    /*--- Gather the element connectivity on all ranks. ---*/
    sizeGlobal = displs.back() + recvCounts.back();

    elemConns.resize(sizeGlobal);

    SU2_MPI::Allgatherv(val_connElem.data(), sizeLocal, MPI_UNSIGNED_LONG, elemConns.data(), recvCounts.data(),
                        displs.data(), MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
  } else {
    /*--- A local tree must be built. Copy the data from the arguments into the
          member variables and set the ranks to the rank of this processor. ---*/
    int rank;
    SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

    coorPoints = val_coor;
    elemConns = val_connElem;
    elemVTK_Type = val_VTKElem;
    localMarkers = val_markerID;
    localElemIDs = val_elemID;

    ranksOfElems.assign(elemVTK_Type.size(), rank);
  }
#else

  /*--- Sequential mode. Copy the data from the arguments into the member
        variables and set the ranks to MASTER_NODE. ---*/
  coorPoints = val_coor;
  elemConns = val_connElem;
  elemVTK_Type = val_VTKElem;
  localMarkers = val_markerID;
  localElemIDs = val_elemID;

  ranksOfElems.assign(elemVTK_Type.size(), MASTER_NODE);

#endif

  /*--- Determine the values of the vector nDOFsPerElem, which contains the
      number of DOFs per element in cumulative storage format. ---*/
  const unsigned long nElem = elemVTK_Type.size();
  nDOFsPerElem.resize(nElem + 1);

  nDOFsPerElem[0] = 0;
  for (unsigned long i = 0; i < nElem; ++i) {
    unsigned short nDOFs = 0;
    switch (elemVTK_Type[i]) {
      case LINE:
        nDOFs = 2;
        break;
      case TRIANGLE:
        nDOFs = 3;
        break;
      case QUADRILATERAL:
        nDOFs = 4;
        break;
      case TETRAHEDRON:
        nDOFs = 4;
        break;
      case PYRAMID:
        nDOFs = 5;
        break;
      case PRISM:
        nDOFs = 6;
        break;
      case HEXAHEDRON:
        nDOFs = 8;
        break;
    }

    nDOFsPerElem[i + 1] = nDOFsPerElem[i] + nDOFs;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Create the bounding boxes of the elements. The coordinates ---*/
  /*---         of these bounding boxes can be interpreted as a point in a ---*/
  /*---         higher dimensional space. The ADT is built as a tree of    ---*/
  /*---         these points in this higher dimensional space.             ---*/
  /*--------------------------------------------------------------------------*/

  /* Allocate the memory for the bounding boxes of the elements. */
  BBoxCoor.resize(2 * nDim * nElem);

  /*--- Loop over the elements to determine the minimum and maximum coordinates
        of the elements. These coordinates completely define the bounding box,
        which can be represented as a point in 2*nDim dimensional space. ---*/
  for (unsigned long i = 0; i < nElem; ++i) {
    /* Set the pointers for the minimum and maximum coordinates for this
       bounding box. Initialize these coordinates to the first coordinate
       of the corresponding element. */
    su2double* BBMin = BBoxCoor.data() + 2 * nDim * i;
    su2double* BBMax = BBMin + nDim;

    unsigned long j = nDOFsPerElem[i];
    const su2double* xP = coorPoints.data() + nDim * elemConns[j];

    for (unsigned short k = 0; k < nDim; ++k) BBMin[k] = BBMax[k] = xP[k];

    /* Loop over the other vertices of the element to determine the
       minimum and maximum coordinates of the bounding box. */
    for (j = (nDOFsPerElem[i] + 1); j < nDOFsPerElem[i + 1]; ++j) {
      xP = coorPoints.data() + nDim * elemConns[j];

      for (unsigned short k = 0; k < nDim; ++k) {
        BBMin[k] = min(BBMin[k], xP[k]);
        BBMax[k] = max(BBMax[k], xP[k]);
      }
    }

    /* Add a tolerance to the bounding box size, such that the overlap check in
       the tree traversal does not go wrong due to round off error. */
    for (unsigned short k = 0; k < nDim; ++k) {
      const su2double lenScale = BBMax[k] - BBMin[k];
      const su2double tol = max(1.e-25, 1.e-6 * lenScale);

      BBMin[k] -= tol;
      BBMax[k] += tol;
    }
  }

  /* Build the ADT of the bounding boxes. */
  BuildADT(2 * nDim, nElem, BBoxCoor.data());

  /*--- Reserve the memory for frontLeaves, frontLeavesNew and BBoxTargets,
        which are needed during the tree search. ---*/
  for (auto& vec : BBoxTargets) vec.reserve(200);
  for (auto& vec : FrontLeaves) vec.reserve(200);
  for (auto& vec : FrontLeavesNew) vec.reserve(200);
}

bool CADTElemClass::DetermineContainingElement_impl(vector<unsigned long>& frontLeaves,
                                                    vector<unsigned long>& frontLeavesNew, const su2double* coor,
                                                    unsigned short& markerID, unsigned long& elemID, int& rankID,
                                                    su2double* parCoor, su2double* weightsInterpol) const {
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
        unsigned long kk = leaves[ll].children[mm];
        if (leaves[ll].childrenAreTerminal[mm]) {
          /* Child contains a bounding box. Check if the coordinate is
             inside the bounding box. */
          const su2double* coorBBMin = BBoxCoor.data() + nDimADT * kk;
          const su2double* coorBBMax = coorBBMin + nDim;

          bool coorIsInside = true;
          for (unsigned short k = 0; k < nDim; ++k) {
            if (coor[k] < coorBBMin[k]) coorIsInside = false;
            if (coor[k] > coorBBMax[k]) coorIsInside = false;
          }

          if (coorIsInside) {
            /* Coordinate is inside the bounding box. Check if it
               is also inside the corresponding element. If so,
               set the required information and return true. */
            if (CoorInElement(kk, coor, parCoor, weightsInterpol)) {
              markerID = localMarkers[kk];
              elemID = localElemIDs[kk];
              rankID = ranksOfElems[kk];
              return true;
            }
          }
        } else {
          /* Child contains a leaf. If the coordinate is inside the leaf
             store the leaf for the next round. */
          const su2double* coorBBMin = leaves[kk].xMin;
          const su2double* coorBBMax = leaves[kk].xMax + nDim;

          bool coorIsInside = true;
          for (unsigned short k = 0; k < nDim; ++k) {
            if (coor[k] < coorBBMin[k]) coorIsInside = false;
            if (coor[k] > coorBBMax[k]) coorIsInside = false;
          }

          if (coorIsInside) frontLeavesNew.push_back(kk);
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

  /* If this point is reached, no element is found that contains the coordinate
     and false must be returned. */
  return false;
}

void CADTElemClass::DetermineNearestElement_impl(vector<CBBoxTargetClass>& BBoxTargets,
                                                 vector<unsigned long>& frontLeaves,
                                                 vector<unsigned long>& frontLeavesNew, const su2double* coor,
                                                 su2double& dist, unsigned short& markerID, unsigned long& elemID,
                                                 int& rankID) const {
  const bool wasActive = AD::BeginPassive();

  /*----------------------------------------------------------------------------*/
  /*--- Step 1: Initialize the distance (squared) to the quaranteed distance ---*/
  /*---         of the central bounding box of the root element.             ---*/
  /*----------------------------------------------------------------------------*/

  unsigned long kk = leaves[0].centralNodeID;
  const su2double* coorBBMin = BBoxCoor.data() + nDimADT * kk;
  const su2double* coorBBMax = coorBBMin + nDim;
  unsigned long jj = 0;

  dist = 0.0;
  su2double ds;
  ds = max(fabs(coor[0] - coorBBMin[0]), fabs(coor[0] - coorBBMax[0]));
  dist += ds * ds;
  ds = max(fabs(coor[1] - coorBBMin[1]), fabs(coor[1] - coorBBMax[1]));
  dist += ds * ds;
  if (nDim == 3) {
    ds = max(fabs(coor[2] - coorBBMin[2]), fabs(coor[2] - coorBBMax[2]));
    dist += ds * ds;
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
          /*--- Child contains a bounding. Compute the possible distance squared
                to this bounding box. Tf a coordinate is inside the box ds = 0. ---*/
          coorBBMin = BBoxCoor.data() + nDimADT * kk;
          coorBBMax = coorBBMin + nDim;

          su2double posDist2 = 0.0, ds;
          ds = min(0.0, coor[0] - coorBBMin[0]) + max(0.0, coor[0] - coorBBMax[0]);
          posDist2 += ds * ds;
          ds = min(0.0, coor[1] - coorBBMin[1]) + max(0.0, coor[1] - coorBBMax[1]);
          posDist2 += ds * ds;
          if (nDim == 3) {
            ds = min(0.0, coor[2] - coorBBMin[2]) + max(0.0, coor[2] - coorBBMax[2]);
            posDist2 += ds * ds;
          }

          /* Check if the possible minimum distance is less than or equal to
             the currently stored distance. If so, this bounding box is a
             candidate for the actual minimum distance and must be stored
             in BBoxTargets. */
          if (posDist2 <= dist) {
            /*--- Compute the guaranteed minimum distance for this bounding box. ---*/
            su2double guarDist2 = 0.0, ds;
            ds = max(fabs(coor[0] - coorBBMin[0]), fabs(coor[0] - coorBBMax[0]));
            guarDist2 += ds * ds;
            ds = max(fabs(coor[1] - coorBBMin[1]), fabs(coor[1] - coorBBMax[1]));
            guarDist2 += ds * ds;
            if (nDim == 3) {
              ds = max(fabs(coor[2] - coorBBMin[2]), fabs(coor[2] - coorBBMax[2]));
              guarDist2 += ds * ds;
            }

            /* Store this bounding box in BBoxTargets and update the currently
               stored value of the distance squared. */
            BBoxTargets.emplace_back(kk, posDist2, guarDist2);
            dist = min(dist, guarDist2);
          }
        } else {
          /*--- Child contains a leaf. Determine the possible minimum distance
                squared to that leaf. Same manual unroling as before. ---*/
          coorBBMin = leaves[kk].xMin;
          coorBBMax = leaves[kk].xMax + nDim;

          su2double posDist2 = 0.0, ds;
          ds = min(0.0, coor[0] - coorBBMin[0]) + max(0.0, coor[0] - coorBBMax[0]);
          posDist2 += ds * ds;
          ds = min(0.0, coor[1] - coorBBMin[1]) + max(0.0, coor[1] - coorBBMax[1]);
          posDist2 += ds * ds;
          if (nDim == 3) {
            ds = min(0.0, coor[2] - coorBBMin[2]) + max(0.0, coor[2] - coorBBMax[2]);
            posDist2 += ds * ds;
          }

          /* Check if the possible minimum distance is less than or equal to the currently
             stored distance. If so this leaf must be stored for the next round. */
          if (posDist2 <= dist) {
            frontLeavesNew.push_back(kk);

            /*--- Determine the guaranteed minimum distance squared to the central
                  bounding box of this leaf and update the currently stored
                  minimum wall distance. Same manual unroling as before. ---*/
            kk = leaves[kk].centralNodeID;
            coorBBMin = BBoxCoor.data() + nDimADT * kk;
            coorBBMax = coorBBMin + nDim;

            su2double guarDist2 = 0.0, ds;
            ds = max(fabs(coor[0] - coorBBMin[0]), fabs(coor[0] - coorBBMax[0]));
            guarDist2 += ds * ds;
            ds = max(fabs(coor[1] - coorBBMin[1]), fabs(coor[1] - coorBBMax[1]));
            guarDist2 += ds * ds;
            if (nDim == 3) {
              ds = max(fabs(coor[2] - coorBBMin[2]), fabs(coor[2] - coorBBMax[2]));
              guarDist2 += ds * ds;
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
    if (frontLeaves.empty()) break;
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
  for (unsigned long i = 0; i < BBoxTargets.size(); ++i) {
    /* Break the loop if the possible minimum distance is larger than
       the currently stored value. In that case it does not make sense to
       check the remainder of the bounding boxes, as they are sorted in
       increasing order (based on the possible minimum distance.
       Make sure that at least one bounding box is checked. */
    if (BBoxTargets[i].possibleMinDist2 > dist) break;

    /*--- Compute the distance squared to the element that corresponds to the
          current bounding box. If this distance is less than or equal to
          the current value, overwrite the return information of this function.
          The equal is necessary to avoid problems for extreme situations. ---*/
    const unsigned long ii = BBoxTargets[i].boundingBoxID;

    su2double dist2Elem;
    Dist2ToElement(ii, coor, dist2Elem);
    if (dist2Elem <= dist) {
      jj = ii;
      dist = dist2Elem;
      markerID = localMarkers[ii];
      elemID = localElemIDs[ii];
      rankID = ranksOfElems[ii];
    }
  }

  AD::EndPassive(wasActive);

  /* At the moment the square of the distance is stored in dist. Compute
     the correct value. */
  Dist2ToElement(jj, coor, dist);
  dist = sqrt(dist);
}

bool CADTElemClass::CoorInElement(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                                  su2double* weightsInterpol) const {
  /*--- Make a distinction between the element types. ---*/
  switch (elemVTK_Type[elemID]) {
    case TRIANGLE:
      return CoorInTriangle(elemID, coor, parCoor, weightsInterpol);
    case QUADRILATERAL:
      return CoorInQuadrilateral(elemID, coor, parCoor, weightsInterpol);
    case TETRAHEDRON:
      return CoorInTetrahedron(elemID, coor, parCoor, weightsInterpol);
    case PYRAMID:
      return CoorInPyramid(elemID, coor, parCoor, weightsInterpol);
    case PRISM:
      return CoorInPrism(elemID, coor, parCoor, weightsInterpol);
    case HEXAHEDRON:
      return CoorInHexahedron(elemID, coor, parCoor, weightsInterpol);

    default:
      /* This should not happen. */
      SU2_MPI::Error("Element type not recognized.", CURRENT_FUNCTION);
      return false;
  }
}

void CADTElemClass::Dist2ToElement(const unsigned long elemID, const su2double* coor, su2double& dist2Elem) const {
  /*--- Make a distinction between the element types. ---*/
  switch (elemVTK_Type[elemID]) {
    case LINE: {
      /*--- Element is a line. The number of space dimensions can be either
            1, 2 or 3. Store the indices where the coordinates of the vertices
            are stored in i0 and i1. ---*/
      unsigned long i0 = nDOFsPerElem[elemID];
      unsigned long i1 = i0 + 1;
      i0 = nDim * elemConns[i0];
      i1 = nDim * elemConns[i1];

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
      i0 = nDim * elemConns[i0];
      i1 = nDim * elemConns[i1];
      i2 = nDim * elemConns[i2];

      /*--- Call the function Dist2ToTriangle to compute the distance to the
            triangle if the projection is inside the triangle. In that case the
            function returns true. If the projection is not inside the triangle,
            false is returned and the distance to each of the lines of the
            triangle is computed and the minimum is taken. ---*/
      su2double r, s;
      if (!Dist2ToTriangle(i0, i1, i2, coor, dist2Elem, r, s)) {
        Dist2ToLine(i0, i1, coor, dist2Elem);

        su2double dist2Line;
        Dist2ToLine(i1, i2, coor, dist2Line);
        dist2Elem = min(dist2Elem, dist2Line);
        Dist2ToLine(i2, i0, coor, dist2Line);
        dist2Elem = min(dist2Elem, dist2Line);
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
      i0 = nDim * elemConns[i0];
      i1 = nDim * elemConns[i1];
      i2 = nDim * elemConns[i2];
      i3 = nDim * elemConns[i3];

      /*--- Call the function Dist2ToQuadrilateral to compute the distance to the
            quadrilateral if the projection is inside the quadrilateral. In that
            case the function returns true. If the projection is not inside the
            quadrilateral, false is returned and the distance to each of the lines
            of the quadrilateral is computed and the minimum is taken. ---*/
      su2double r, s;
      if (!Dist2ToQuadrilateral(i0, i1, i2, i3, coor, r, s, dist2Elem)) {
        Dist2ToLine(i0, i1, coor, dist2Elem);

        su2double dist2Line;
        Dist2ToLine(i1, i2, coor, dist2Line);
        dist2Elem = min(dist2Elem, dist2Line);
        Dist2ToLine(i2, i3, coor, dist2Line);
        dist2Elem = min(dist2Elem, dist2Line);
        Dist2ToLine(i3, i0, coor, dist2Line);
        dist2Elem = min(dist2Elem, dist2Line);
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

bool CADTElemClass::CoorInTriangle(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                                   su2double* weightsInterpol) const {
  /* Determine the indices of the three vertices of the triangle,
     multiplied by nDim (which is 2). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2;
  i0 = nDim * elemConns[i0];
  i1 = nDim * elemConns[i1];
  i2 = nDim * elemConns[i2];

  /* Determine the coordinates relative to the vertex 0. */
  const su2double xc = coor[0] - coorPoints[i0];
  const su2double yc = coor[1] - coorPoints[i0 + 1];

  const su2double x1 = coorPoints[i1] - coorPoints[i0];
  const su2double y1 = coorPoints[i1 + 1] - coorPoints[i0 + 1];

  const su2double x2 = coorPoints[i2] - coorPoints[i0];
  const su2double y2 = coorPoints[i2 + 1] - coorPoints[i0 + 1];

  /* The triangle is parametrized by X-X0 = (r+1)*(X1-X0)/2 + (s+1)*(X2-X0)/2,
     r, s >= -1, r+s <= 0. As this is a containment search, the number of
     dimesions is 2. As a consequence, the parametric coordinates r and s
     can be solved easily. Note that X0 is 0 in the above expression, because
     the coordinates are relative to node 0. */
  const su2double detInv = 2.0 / (x1 * y2 - x2 * y1);
  parCoor[0] = detInv * (xc * y2 - yc * x2) - 1.0;
  parCoor[1] = detInv * (yc * x1 - xc * y1) - 1.0;

  /* Check if the point resides within the triangle and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
      ((parCoor[0] + parCoor[1]) <= tolInsideElem)) {
    coorIsInside = true;

    weightsInterpol[0] = -0.5 * (parCoor[0] + parCoor[1]);
    weightsInterpol[1] = 0.5 * (parCoor[0] + 1.0);
    weightsInterpol[2] = 0.5 * (parCoor[1] + 1.0);
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::CoorInQuadrilateral(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                                        su2double* weightsInterpol) const {
  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton = 1.e-10;

  /* Determine the indices of the four vertices of the quadrilatral,
     multiplied by nDim (which is 2). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0 + 3;
  i0 = nDim * elemConns[i0];
  i1 = nDim * elemConns[i1];
  i2 = nDim * elemConns[i2];
  i3 = nDim * elemConns[i3];

  /* Determine the coordinates relative to the vertex 0. */
  const su2double xc = coor[0] - coorPoints[i0];
  const su2double yc = coor[1] - coorPoints[i0 + 1];

  const su2double x1 = coorPoints[i1] - coorPoints[i0];
  const su2double y1 = coorPoints[i1 + 1] - coorPoints[i0 + 1];

  const su2double x2 = coorPoints[i2] - coorPoints[i0];
  const su2double y2 = coorPoints[i2 + 1] - coorPoints[i0 + 1];

  const su2double x3 = coorPoints[i3] - coorPoints[i0];
  const su2double y3 = coorPoints[i3 + 1] - coorPoints[i0 + 1];

  /* The parametrization of the quadrilatral is nonlinear, which requires an
     iterative algorithm. Especially for highly skewed quadrilaterals, this
     could lead to convergence problems. Hence the quadrilateral is split into
     linear triangles to check if the point is actually within the quad.
     First check the triangle i0-i1-i3. See CoorInTriangle for more details
     on this test. */
  su2double detInv = 2.0 / (x1 * y3 - x3 * y1);
  parCoor[0] = detInv * (xc * y3 - yc * x3) - 1.0;
  parCoor[1] = detInv * (yc * x1 - xc * y1) - 1.0;

  bool coorIsInside = false;
  if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
      ((parCoor[0] + parCoor[1]) <= tolInsideElem))
    coorIsInside = true;

  /* Check triangle i2-i3-i1 if the coordinate is not inside i0-i1-i3. */
  if (!coorIsInside) {
    /* Define the coordinates w.r.t. vertex 2 using the numbering i2-i3-i1. */
    const su2double xxc = xc - x2, yyc = yc - y2;
    const su2double xx1 = x3 - x2, yy1 = y3 - y2;
    const su2double xx3 = x1 - x2, yy3 = y1 - y2;

    /* Check if the coordinate is inside this triangle. */
    detInv = 2.0 / (xx1 * yy3 - xx3 * yy1);
    parCoor[0] = detInv * (xxc * yy3 - yyc * xx3) - 1.0;
    parCoor[1] = detInv * (yyc * xx1 - xxc * yy1) - 1.0;

    if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
        ((parCoor[0] + parCoor[1]) <= tolInsideElem))
      coorIsInside = true;

    /* Convert the parametric coordinates to the ones used by the
       quadrilatral i0-i1-i2-i3. They serve as initial guess below. */
    parCoor[0] = -parCoor[0];
    parCoor[1] = -parCoor[1];
  }

  /* If the coordinate is in neither triangle, return false. */
  if (!coorIsInside) return false;

  /* The coordinate is inside the quadrilatral and an initial guess has been
     obtained by splitting the quad into two triangles. Carry out a Newton
     algorithm to obtain the true parametric coordinates.
     The quadrilateral is parametrized by
     X = {(1-r)*(1-s)*X0 + (1+r)*(1-s)*X1 + (1+r)*(1+s)*X2 + (1-r)*(1+s)*X3}/4,
     -1 <= r,s <= 1. As the coordinates are relative to X0, the first term
     drops from this equation. The nonlinear set of equations can be written as
     V0 - r*V1 - s*V2 - r*s*V3 = 0, where V0 = xc - (x1+x2+x3)/4,
     V1 = (x1+x2-x3)/4, V2 = (x2+x3-X1)/4, V3 = (x2-x1-x3)/4. First
     construct the vectors V0, V1, V2 and V3. */
  const su2double V0x = xc - 0.25 * (x1 + x2 + x3), V0y = yc - 0.25 * (y1 + y2 + y3);
  const su2double V1x = 0.25 * (x1 + x2 - x3), V1y = 0.25 * (y1 + y2 - y3);
  const su2double V2x = 0.25 * (x2 + x3 - x1), V2y = 0.25 * (y2 + y3 - y1);
  const su2double V3x = 0.25 * (x2 - x1 - x3), V3y = 0.25 * (y2 - y1 - y3);

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for (itCount = 0; itCount < maxIt; ++itCount) {
    /* Compute the values of the nonlinear functions. */
    const su2double f0 = V0x - parCoor[0] * V1x - parCoor[1] * V2x - parCoor[0] * parCoor[1] * V3x;
    const su2double f1 = V0y - parCoor[0] * V1y - parCoor[1] * V2y - parCoor[0] * parCoor[1] * V3y;

    /* Compute the negative of the Jacobian matrix. */
    const su2double a00 = V1x + parCoor[1] * V3x, a01 = V2x + parCoor[0] * V3x;
    const su2double a10 = V1y + parCoor[1] * V3y, a11 = V2y + parCoor[0] * V3y;

    /* Compute the update of the parametric coordinates. */
    detInv = 1.0 / (a00 * a11 - a01 * a10);
    const su2double dr = detInv * (f0 * a11 - f1 * a01);
    const su2double ds = detInv * (f1 * a00 - f0 * a10);

    /* Update the parametric coordinates. Note that the negative of the
       Jacobian is used, so the update must be added. */
    parCoor[0] += dr;
    parCoor[1] += ds;

    /* Check for convergence. */
    if (fabs(dr) <= tolNewton && fabs(ds) <= tolNewton) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if (itCount == maxIt) SU2_MPI::Error("Newton did not converge", CURRENT_FUNCTION);

  /* Check if the parametric coordinates are inside the quadrilateral.
     If not, something is seriously wrong, because the inside test has been
     done already with the triangles. */
  if (parCoor[0] < paramLowerBound || parCoor[0] > paramUpperBound || parCoor[1] < paramLowerBound ||
      parCoor[1] > paramUpperBound)
    SU2_MPI::Error("Point not inside the quadrilateral.", CURRENT_FUNCTION);

  /* Compute the interpolation weights. */
  const su2double omr = 0.5 * (1.0 - parCoor[0]), opr = 0.5 * (1.0 + parCoor[0]);
  const su2double oms = 0.5 * (1.0 - parCoor[1]), ops = 0.5 * (1.0 + parCoor[1]);

  weightsInterpol[0] = omr * oms;
  weightsInterpol[1] = opr * oms;
  weightsInterpol[2] = opr * ops;
  weightsInterpol[3] = omr * ops;

  /* Return true, because the search was successful. */
  return true;
}

bool CADTElemClass::CoorInTetrahedron(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                                      su2double* weightsInterpol) const {
  /* Determine the indices of the four vertices of the tetrahedron,
     multiplied by nDim (which is 3). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0 + 3;
  i0 = nDim * elemConns[i0];
  i1 = nDim * elemConns[i1];
  i2 = nDim * elemConns[i2];
  i3 = nDim * elemConns[i3];

  /* Determine the coordinates relative to the vertex 0. */
  const su2double xc = coor[0] - coorPoints[i0];
  const su2double yc = coor[1] - coorPoints[i0 + 1];
  const su2double zc = coor[2] - coorPoints[i0 + 2];

  const su2double x1 = coorPoints[i1] - coorPoints[i0];
  const su2double y1 = coorPoints[i1 + 1] - coorPoints[i0 + 1];
  const su2double z1 = coorPoints[i1 + 2] - coorPoints[i0 + 2];

  const su2double x2 = coorPoints[i2] - coorPoints[i0];
  const su2double y2 = coorPoints[i2 + 1] - coorPoints[i0 + 1];
  const su2double z2 = coorPoints[i2 + 2] - coorPoints[i0 + 2];

  const su2double x3 = coorPoints[i3] - coorPoints[i0];
  const su2double y3 = coorPoints[i3 + 1] - coorPoints[i0 + 1];
  const su2double z3 = coorPoints[i3 + 2] - coorPoints[i0 + 2];

  /* The tetrahedron is parametrized by
     X-X0 = (r+1)*(X1-X0)/2 + (s+1)*(X2-X0)/2 + (t+1)*(X3-X0)/2,
     r, s, t >= -1, r+s+t <= -1. As a consequence, the parametric coordinates
     r, s and t can be solved easily. Note that X0 is 0 in the above expression,
     because the coordinates are relative to node 0. */
  const su2double detInv =
      2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  parCoor[0] = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  parCoor[1] =
      -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  parCoor[2] = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* Check if the point resides within the tetrahedron and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) && (parCoor[2] >= paramLowerBound) &&
      ((parCoor[0] + parCoor[1] + parCoor[2]) <= paramLowerBound)) {
    coorIsInside = true;

    weightsInterpol[0] = -0.5 * (parCoor[0] + parCoor[1] + parCoor[2] + 1.0);
    weightsInterpol[1] = 0.5 * (parCoor[0] + 1.0);
    weightsInterpol[2] = 0.5 * (parCoor[1] + 1.0);
    weightsInterpol[3] = 0.5 * (parCoor[2] + 1.0);
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::CoorInPyramid(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                                  su2double* weightsInterpol) const {
  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton = 1.e-10;

  /* Determine the indices of the five vertices of the pyramid,
     multiplied by nDim (which is 3). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0 + 3, i4 = i0 + 4;
  i0 = nDim * elemConns[i0];
  i1 = nDim * elemConns[i1];
  i2 = nDim * elemConns[i2];
  i3 = nDim * elemConns[i3];
  i4 = nDim * elemConns[i4];

  /* Determine the coordinates relative to the vertex 0. */
  su2double xRel[5][3], xc[3];
  xc[0] = coor[0] - coorPoints[i0];
  xc[1] = coor[1] - coorPoints[i0 + 1];
  xc[2] = coor[2] - coorPoints[i0 + 2];

  xRel[0][0] = xRel[0][1] = xRel[0][2] = 0.0;

  xRel[1][0] = coorPoints[i1] - coorPoints[i0];
  xRel[1][1] = coorPoints[i1 + 1] - coorPoints[i0 + 1];
  xRel[1][2] = coorPoints[i1 + 2] - coorPoints[i0 + 2];

  xRel[2][0] = coorPoints[i2] - coorPoints[i0];
  xRel[2][1] = coorPoints[i2 + 1] - coorPoints[i0 + 1];
  xRel[2][2] = coorPoints[i2 + 2] - coorPoints[i0 + 2];

  xRel[3][0] = coorPoints[i3] - coorPoints[i0];
  xRel[3][1] = coorPoints[i3 + 1] - coorPoints[i0 + 1];
  xRel[3][2] = coorPoints[i3 + 2] - coorPoints[i0 + 2];

  xRel[4][0] = coorPoints[i4] - coorPoints[i0];
  xRel[4][1] = coorPoints[i4 + 1] - coorPoints[i0 + 1];
  xRel[4][2] = coorPoints[i4 + 2] - coorPoints[i0 + 2];

  /* Obtain an initial guess of the parametric coordinates by splitting the
     pyramid into tetrahedra. If this approach is not successful, this means
     that the point is not inside the true pyramid and false can be returned. */
  if (!InitialGuessContainmentPyramid(xc, xRel, parCoor)) return false;

  /* The pyramid is parametrized by X-X0 = (Xi-X0)*li, where the sum runs over
     i = 0..4, although i = 0 does not give a contribution. The Lagrangian
     interpolation functions for the pyramid are given by:
     l0 = (1-t-2*r)*(1-t-2*s)/(8*(1-t)), l1 = (1-t+2*r)*(1-t-2*s)/(8*(1-t))
     l2 = (1-t+2*r)*(1-t+2*s)/(8*(1-t)), l3 = (1-t-2*r)*(1-t+2*s)/(8*(1-t)),
     l4 = (1+t)/2.
     The boundaries are -1 <= t <= 1, (t-1)/2 <= r,s <= (1-t)/2.
     As all coordinates are taken relative to vertex 0, X0 drops out and the
     nonlinear set of equations can be written as
     V0 - V1*r - V2*s - V3*t - V4*r*s/(1-t)  = 0, where
     V0 = xc - (4*x4+x1+x2+x3)/8, V1 = (x1+x2-x3)/4, V2 = (x2+x3-X1)/4,
     V3 = (4*X4-x1-x2-x3)/8, V4 = (x2-x1-x3)/2. First construct these vectors. */
  const su2double V0x = xc[0] - 0.5 * xRel[4][0] - 0.125 * (xRel[1][0] + xRel[2][0] + xRel[3][0]);
  const su2double V0y = xc[1] - 0.5 * xRel[4][1] - 0.125 * (xRel[1][1] + xRel[2][1] + xRel[3][1]);
  const su2double V0z = xc[2] - 0.5 * xRel[4][2] - 0.125 * (xRel[1][2] + xRel[2][2] + xRel[3][2]);

  const su2double V1x = 0.25 * (xRel[1][0] + xRel[2][0] - xRel[3][0]);
  const su2double V1y = 0.25 * (xRel[1][1] + xRel[2][1] - xRel[3][1]);
  const su2double V1z = 0.25 * (xRel[1][2] + xRel[2][2] - xRel[3][2]);

  const su2double V2x = 0.25 * (xRel[2][0] + xRel[3][0] - xRel[1][0]);
  const su2double V2y = 0.25 * (xRel[2][1] + xRel[3][1] - xRel[1][1]);
  const su2double V2z = 0.25 * (xRel[2][2] + xRel[3][2] - xRel[1][2]);

  const su2double V3x = 0.5 * xRel[4][0] - 0.125 * (xRel[1][0] + xRel[2][0] + xRel[3][0]);
  const su2double V3y = 0.5 * xRel[4][1] - 0.125 * (xRel[1][1] + xRel[2][1] + xRel[3][1]);
  const su2double V3z = 0.5 * xRel[4][2] - 0.125 * (xRel[1][2] + xRel[2][2] + xRel[3][2]);

  const su2double V4x = 0.5 * (xRel[2][0] - xRel[1][0] - xRel[3][0]);
  const su2double V4y = 0.5 * (xRel[2][1] - xRel[1][1] - xRel[3][1]);
  const su2double V4z = 0.5 * (xRel[2][2] - xRel[1][2] - xRel[3][2]);

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for (itCount = 0; itCount < maxIt; ++itCount) {
    /* Compute the values of the nonlinear functions. */
    su2double oneMinT = 1.0 - parCoor[2];
    if (fabs(oneMinT) < 1.e-10) {
      if (oneMinT < 0.0)
        oneMinT = -1.e-10;
      else
        oneMinT = 1.e-10;
    }
    const su2double oneMinTInv = 1.0 / oneMinT;

    const su2double f0 =
        V0x - parCoor[0] * V1x - parCoor[1] * V2x - parCoor[2] * V3x - parCoor[0] * parCoor[1] * V4x * oneMinTInv;
    const su2double f1 =
        V0y - parCoor[0] * V1y - parCoor[1] * V2y - parCoor[2] * V3y - parCoor[0] * parCoor[1] * V4y * oneMinTInv;
    const su2double f2 =
        V0z - parCoor[0] * V1z - parCoor[1] * V2z - parCoor[2] * V3z - parCoor[0] * parCoor[1] * V4z * oneMinTInv;

    /* Compute the negative of the Jacobian matrix. */
    const su2double a00 = V1x + parCoor[1] * V4x * oneMinTInv;
    const su2double a01 = V2x + parCoor[0] * V4x * oneMinTInv;
    const su2double a02 = V3x + parCoor[0] * parCoor[1] * V4x * oneMinTInv * oneMinTInv;

    const su2double a10 = V1y + parCoor[1] * V4y * oneMinTInv;
    const su2double a11 = V2y + parCoor[0] * V4y * oneMinTInv;
    const su2double a12 = V3y + parCoor[0] * parCoor[1] * V4y * oneMinTInv * oneMinTInv;

    const su2double a20 = V1z + parCoor[1] * V4z * oneMinTInv;
    const su2double a21 = V2z + parCoor[0] * V4z * oneMinTInv;
    const su2double a22 = V3z + parCoor[0] * parCoor[1] * V4z * oneMinTInv * oneMinTInv;

    /* Compute the update of the parametric coordinates. */
    const su2double detInv = 1.0 / (a00 * a11 * a22 - a00 * a12 * a21 - a01 * a10 * a22 + a01 * a12 * a20 +
                                    a02 * a10 * a21 - a02 * a11 * a20);
    const su2double dr =
        detInv * (a01 * a12 * f2 - a01 * a22 * f1 - a02 * a11 * f2 + a02 * a21 * f1 + a11 * a22 * f0 - a12 * a21 * f0);
    const su2double ds =
        -detInv * (a00 * a12 * f2 - a00 * a22 * f1 - a02 * a10 * f2 + a02 * a20 * f1 + a10 * a22 * f0 - a12 * a20 * f0);
    const su2double dt =
        detInv * (a00 * a11 * f2 - a00 * a21 * f1 - a01 * a10 * f2 + a01 * a20 * f1 + a10 * a21 * f0 - a11 * a20 * f0);

    /* Update the parametric coordinates. Note that the negative of the
       Jacobian is used, so the update must be added. */
    parCoor[0] += dr;
    parCoor[1] += ds;
    parCoor[2] += dt;

    /* Check for convergence. */
    if (fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if (itCount == maxIt) SU2_MPI::Error("Newton did not converge", CURRENT_FUNCTION);

  /* Check if the point resides within the pyramid and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if ((parCoor[2] >= paramLowerBound) && (parCoor[2] <= paramUpperBound)) {
    const su2double lowRSBound = 0.5 * (parCoor[2] - 1.0) - tolInsideElem;
    const su2double uppRSBound = -lowRSBound;

    if ((parCoor[0] >= lowRSBound) && (parCoor[0] <= uppRSBound) && (parCoor[1] >= lowRSBound) &&
        (parCoor[1] <= uppRSBound)) {
      coorIsInside = true;

      su2double oneMinT = 1.0 - parCoor[2];
      if (fabs(oneMinT) < 1.e-10) {
        if (oneMinT < 0.0)
          oneMinT = -1.e-10;
        else
          oneMinT = 1.e-10;
      }
      const su2double oneMinTInv = 1.0 / oneMinT;

      const su2double omr = (1.0 - parCoor[2] - 2.0 * parCoor[0]);
      const su2double opr = (1.0 - parCoor[2] + 2.0 * parCoor[0]);
      const su2double oms = (1.0 - parCoor[2] - 2.0 * parCoor[1]);
      const su2double ops = (1.0 - parCoor[2] + 2.0 * parCoor[1]);

      weightsInterpol[0] = 0.125 * oneMinTInv * omr * oms;
      weightsInterpol[1] = 0.125 * oneMinTInv * opr * oms;
      weightsInterpol[2] = 0.125 * oneMinTInv * opr * ops;
      weightsInterpol[3] = 0.125 * oneMinTInv * omr * ops;
      weightsInterpol[4] = 0.5 * (1.0 + parCoor[2]);
    }
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::InitialGuessContainmentPyramid(const su2double xRelC[3], const su2double xRel[5][3],
                                                   su2double* parCoor) const {
  /* Tetrahedron, 0-1-3-4.
     Create the coordinates of the tetrahedron and of the point. */
  su2double x1 = xRel[1][0], y1 = xRel[1][1], z1 = xRel[1][2];
  su2double x2 = xRel[3][0], y2 = xRel[3][1], z2 = xRel[3][2];
  su2double x3 = xRel[4][0], y3 = xRel[4][1], z3 = xRel[4][2];

  su2double xc = xRelC[0], yc = xRelC[1], zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  su2double detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  parCoor[0] = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  parCoor[1] =
      -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  parCoor[2] = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, return true. */
  if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) && (parCoor[2] >= paramLowerBound) &&
      ((parCoor[0] + parCoor[1] + parCoor[2]) <= paramLowerBound))
    return true;

  /* Tetrahedron, 2-3-1-4.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[3][0] - xRel[2][0];
  y1 = xRel[3][1] - xRel[2][1];
  z1 = xRel[3][2] - xRel[2][2];
  x2 = xRel[1][0] - xRel[2][0];
  y2 = xRel[1][1] - xRel[2][1];
  z2 = xRel[1][2] - xRel[2][2];
  x3 = xRel[4][0] - xRel[2][0];
  y3 = xRel[4][1] - xRel[2][1];
  z3 = xRel[4][2] - xRel[2][2];

  xc = xRelC[0] - xRel[2][0];
  yc = xRelC[1] - xRel[2][1];
  zc = xRelC[2] - xRel[2][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  parCoor[0] = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  parCoor[1] =
      -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  parCoor[2] = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, adapt the parametric coordinates
     to the real pyramid and return true. */
  if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) && (parCoor[2] >= paramLowerBound) &&
      ((parCoor[0] + parCoor[1] + parCoor[2]) <= paramLowerBound)) {
    parCoor[0] = 1.0 - parCoor[0];
    parCoor[1] = 1.0 - parCoor[1];
    return true;
  }

  /* Tetrahedron, 1-2-0-4.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[2][0] - xRel[1][0];
  y1 = xRel[2][1] - xRel[1][1];
  z1 = xRel[2][2] - xRel[1][2];
  x2 = xRel[0][0] - xRel[1][0];
  y2 = xRel[0][1] - xRel[1][1];
  z2 = xRel[0][2] - xRel[1][2];
  x3 = xRel[4][0] - xRel[1][0];
  y3 = xRel[4][1] - xRel[1][1];
  z3 = xRel[4][2] - xRel[1][2];

  xc = xRelC[0] - xRel[1][0];
  yc = xRelC[1] - xRel[1][1];
  zc = xRelC[2] - xRel[1][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  parCoor[0] = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  parCoor[1] =
      -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  parCoor[2] = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, adapt the parametric coordinates
     to the real pyramid and return true. */
  if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) && (parCoor[2] >= paramLowerBound) &&
      ((parCoor[0] + parCoor[1] + parCoor[2]) <= paramLowerBound)) {
    const su2double r = parCoor[0];
    parCoor[0] = 1.0 - parCoor[1];
    parCoor[1] = r;
    return true;
  }

  /* Tetrahedron, 3-0-2-4.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[0][0] - xRel[3][0];
  y1 = xRel[0][1] - xRel[3][1];
  z1 = xRel[0][2] - xRel[3][2];
  x2 = xRel[2][0] - xRel[3][0];
  y2 = xRel[2][1] - xRel[3][1];
  z2 = xRel[2][2] - xRel[3][2];
  x3 = xRel[4][0] - xRel[3][0];
  y3 = xRel[4][1] - xRel[3][1];
  z3 = xRel[4][2] - xRel[3][2];

  xc = xRelC[0] - xRel[3][0];
  yc = xRelC[1] - xRel[3][1];
  zc = xRelC[2] - xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  parCoor[0] = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  parCoor[1] =
      -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  parCoor[2] = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, adapt the parametric coordinates
     to the real pyramid and return true. */
  if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) && (parCoor[2] >= paramLowerBound) &&
      ((parCoor[0] + parCoor[1] + parCoor[2]) <= paramLowerBound)) {
    const su2double r = parCoor[0];
    parCoor[0] = parCoor[1];
    parCoor[1] = 1.0 - r;
    return true;
  }

  /* The coordinate is in none of the four sub-tetrahedra of the pyramid. This
     implies that the point is not inside the true pyramid either and hence
     false is returned. */
  return false;
}

bool CADTElemClass::CoorInPrism(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                                su2double* weightsInterpol) const {
  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton = 1.e-10;

  /* Determine the indices of the six vertices of the prism,
     multiplied by nDim (which is 3). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0 + 3, i4 = i0 + 4, i5 = i0 + 5;
  i0 = nDim * elemConns[i0];
  i1 = nDim * elemConns[i1];
  i2 = nDim * elemConns[i2];
  i3 = nDim * elemConns[i3];
  i4 = nDim * elemConns[i4];
  i5 = nDim * elemConns[i5];

  /* Determine the coordinates relative to the vertex 0. */
  su2double xRel[6][3], xc[3];
  xc[0] = coor[0] - coorPoints[i0];
  xc[1] = coor[1] - coorPoints[i0 + 1];
  xc[2] = coor[2] - coorPoints[i0 + 2];

  xRel[0][0] = xRel[0][1] = xRel[0][2] = 0.0;

  xRel[1][0] = coorPoints[i1] - coorPoints[i0];
  xRel[1][1] = coorPoints[i1 + 1] - coorPoints[i0 + 1];
  xRel[1][2] = coorPoints[i1 + 2] - coorPoints[i0 + 2];

  xRel[2][0] = coorPoints[i2] - coorPoints[i0];
  xRel[2][1] = coorPoints[i2 + 1] - coorPoints[i0 + 1];
  xRel[2][2] = coorPoints[i2 + 2] - coorPoints[i0 + 2];

  xRel[3][0] = coorPoints[i3] - coorPoints[i0];
  xRel[3][1] = coorPoints[i3 + 1] - coorPoints[i0 + 1];
  xRel[3][2] = coorPoints[i3 + 2] - coorPoints[i0 + 2];

  xRel[4][0] = coorPoints[i4] - coorPoints[i0];
  xRel[4][1] = coorPoints[i4 + 1] - coorPoints[i0 + 1];
  xRel[4][2] = coorPoints[i4 + 2] - coorPoints[i0 + 2];

  xRel[5][0] = coorPoints[i5] - coorPoints[i0];
  xRel[5][1] = coorPoints[i5 + 1] - coorPoints[i0 + 1];
  xRel[5][2] = coorPoints[i5 + 2] - coorPoints[i0 + 2];

  /* Obtain an initial guess of the parametric coordinates by splitting the
     prism into tetrahedra. If this approach is not successful, this means
     that the point is not inside the true prism and false can be returned. */
  if (!InitialGuessContainmentPrism(xc, xRel, parCoor)) return false;

  /* The prism is parametrized by X-X0 = (Xi-X0)*li, where the sum runs over
     i = 0..5, although i = 0 does not give a contribution. The Lagrangian
     interpolation functions for the prism are given by:
     l0 = -(r+s)*(1-t)/4, l1 = (r+1)*(1-t)/4, l2 = (s+1)*(1-t)/4,
     l3 = -(r+s)*(1+t)/4, l4 = (r+1)*(1+t)/4, l5 = (s+1)*(1+t)/4.
     The boundaries are r,s >= -1, r+s <= 0, -1 <= t <= 1.
     As all coordinates are taken relative to vertex 0, X0 drops out and the
     nonlinear set of equations can be written as
     V0 - V1*r - V2*s - V3*t - V4*r*t - V5*s*t  = 0, where
     V0 = xc - (x1+x2+x4+x5)/4, V1 = (x1+x4-x3)/4, V2 = (x2+x5-x3)/4,
     V3 = (x4+x5-x1-x2)/4, V4 = (x4-x1-x3)/4, V5 = (x5-x2-x3)/4.
     First construct these vectors. */
  const su2double V0x = xc[0] - 0.25 * (xRel[1][0] + xRel[2][0] + xRel[4][0] + xRel[5][0]);
  const su2double V0y = xc[1] - 0.25 * (xRel[1][1] + xRel[2][1] + xRel[4][1] + xRel[5][1]);
  const su2double V0z = xc[2] - 0.25 * (xRel[1][2] + xRel[2][2] + xRel[4][2] + xRel[5][2]);

  const su2double V1x = 0.25 * (xRel[1][0] + xRel[4][0] - xRel[3][0]);
  const su2double V1y = 0.25 * (xRel[1][1] + xRel[4][1] - xRel[3][1]);
  const su2double V1z = 0.25 * (xRel[1][2] + xRel[4][2] - xRel[3][2]);

  const su2double V2x = 0.25 * (xRel[2][0] + xRel[5][0] - xRel[3][0]);
  const su2double V2y = 0.25 * (xRel[2][1] + xRel[5][1] - xRel[3][1]);
  const su2double V2z = 0.25 * (xRel[2][2] + xRel[5][2] - xRel[3][2]);

  const su2double V3x = 0.25 * (xRel[4][0] + xRel[5][0] - xRel[1][0] - xRel[2][0]);
  const su2double V3y = 0.25 * (xRel[4][1] + xRel[5][1] - xRel[1][1] - xRel[2][1]);
  const su2double V3z = 0.25 * (xRel[4][2] + xRel[5][2] - xRel[1][2] - xRel[2][2]);

  const su2double V4x = 0.25 * (xRel[4][0] - xRel[1][0] - xRel[3][0]);
  const su2double V4y = 0.25 * (xRel[4][1] - xRel[1][1] - xRel[3][1]);
  const su2double V4z = 0.25 * (xRel[4][2] - xRel[1][2] - xRel[3][2]);

  const su2double V5x = 0.25 * (xRel[5][0] - xRel[2][0] - xRel[3][0]);
  const su2double V5y = 0.25 * (xRel[5][1] - xRel[2][1] - xRel[3][1]);
  const su2double V5z = 0.25 * (xRel[5][2] - xRel[2][2] - xRel[3][2]);

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for (itCount = 0; itCount < maxIt; ++itCount) {
    const su2double f0 = V0x - parCoor[0] * V1x - parCoor[1] * V2x - parCoor[2] * V3x - parCoor[0] * parCoor[2] * V4x -
                         parCoor[1] * parCoor[2] * V5x;
    const su2double f1 = V0y - parCoor[0] * V1y - parCoor[1] * V2y - parCoor[2] * V3y - parCoor[0] * parCoor[2] * V4y -
                         parCoor[1] * parCoor[2] * V5y;
    const su2double f2 = V0z - parCoor[0] * V1z - parCoor[1] * V2z - parCoor[2] * V3z - parCoor[0] * parCoor[2] * V4z -
                         parCoor[1] * parCoor[2] * V5z;

    /* Compute the negative of the Jacobian matrix. */
    const su2double a00 = V1x + parCoor[2] * V4x;
    const su2double a01 = V2x + parCoor[2] * V5x;
    const su2double a02 = V3x + parCoor[0] * V4x + parCoor[1] * V5x;

    const su2double a10 = V1y + parCoor[2] * V4y;
    const su2double a11 = V2y + parCoor[2] * V5y;
    const su2double a12 = V3y + parCoor[0] * V4y + parCoor[1] * V5y;

    const su2double a20 = V1z + parCoor[2] * V4z;
    const su2double a21 = V2z + parCoor[2] * V5z;
    const su2double a22 = V3z + parCoor[0] * V4z + parCoor[1] * V5z;

    /* Compute the update of the parametric coordinates. */
    const su2double detInv = 1.0 / (a00 * a11 * a22 - a00 * a12 * a21 - a01 * a10 * a22 + a01 * a12 * a20 +
                                    a02 * a10 * a21 - a02 * a11 * a20);
    const su2double dr =
        detInv * (a01 * a12 * f2 - a01 * a22 * f1 - a02 * a11 * f2 + a02 * a21 * f1 + a11 * a22 * f0 - a12 * a21 * f0);
    const su2double ds =
        -detInv * (a00 * a12 * f2 - a00 * a22 * f1 - a02 * a10 * f2 + a02 * a20 * f1 + a10 * a22 * f0 - a12 * a20 * f0);
    const su2double dt =
        detInv * (a00 * a11 * f2 - a00 * a21 * f1 - a01 * a10 * f2 + a01 * a20 * f1 + a10 * a21 * f0 - a11 * a20 * f0);

    /* Update the parametric coordinates. Note that the negative of the
       Jacobian is used, so the update must be added. */
    parCoor[0] += dr;
    parCoor[1] += ds;
    parCoor[2] += dt;

    /* Check for convergence. */
    if (fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if (itCount == maxIt) SU2_MPI::Error("Newton did not converge", CURRENT_FUNCTION);

  /* Check if the point resides within the prism and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if ((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
      ((parCoor[0] + parCoor[1]) <= tolInsideElem) && (parCoor[2] >= paramLowerBound) &&
      (parCoor[2] <= paramUpperBound)) {
    coorIsInside = true;

    const su2double omt = 0.25 * (1.0 - parCoor[2]), opt = 0.25 * (1.0 + parCoor[2]);

    weightsInterpol[0] = -omt * (parCoor[0] + parCoor[1]);
    weightsInterpol[1] = omt * (parCoor[0] + 1.0);
    weightsInterpol[2] = omt * (parCoor[1] + 1.0);
    weightsInterpol[3] = -opt * (parCoor[0] + parCoor[1]);
    weightsInterpol[4] = opt * (parCoor[0] + 1.0);
    weightsInterpol[5] = opt * (parCoor[1] + 1.0);
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::InitialGuessContainmentPrism(const su2double xRelC[3], const su2double xRel[6][3],
                                                 su2double* parCoor) const {
  /* Tetrahedron, 0-1-2-3.
     Create the coordinates of the tetrahedron and of the point. */
  su2double x1 = xRel[1][0], y1 = xRel[1][1], z1 = xRel[1][2];
  su2double x2 = xRel[2][0], y2 = xRel[2][1], z2 = xRel[2][2];
  su2double x3 = xRel[3][0], y3 = xRel[3][1], z3 = xRel[3][2];

  su2double xc = xRelC[0], yc = xRelC[1], zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  su2double detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  su2double r =
      detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  su2double s =
      -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  su2double t =
      detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = r;
    parCoor[1] = s;
    parCoor[2] = t;
    return true;
  }

  /* Tetrahedron, 4-1-3-2.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[1][0] - xRel[4][0];
  y1 = xRel[1][1] - xRel[4][1];
  z1 = xRel[1][2] - xRel[4][2];
  x2 = xRel[3][0] - xRel[4][0];
  y2 = xRel[3][1] - xRel[4][1];
  z2 = xRel[3][2] - xRel[4][2];
  x3 = xRel[2][0] - xRel[4][0];
  y3 = xRel[2][1] - xRel[4][1];
  z3 = xRel[2][2] - xRel[4][2];

  xc = xRelC[0] - xRel[4][0];
  yc = xRelC[1] - xRel[4][1];
  zc = xRelC[2] - xRel[4][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - s;
    parCoor[1] = t;
    parCoor[2] = 1.0 - r;
    return true;
  }

  /* Tetrahedron, 3-5-4-2.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[5][0] - xRel[3][0];
  y1 = xRel[5][1] - xRel[3][1];
  z1 = xRel[5][2] - xRel[3][2];
  x2 = xRel[4][0] - xRel[3][0];
  y2 = xRel[4][1] - xRel[3][1];
  z2 = xRel[4][2] - xRel[3][2];
  x3 = xRel[2][0] - xRel[3][0];
  y3 = xRel[2][1] - xRel[3][1];
  z3 = xRel[2][2] - xRel[3][2];

  xc = xRelC[0] - xRel[3][0];
  yc = xRelC[1] - xRel[3][1];
  zc = xRelC[2] - xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = s;
    parCoor[1] = r;
    parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 3-5-4-0.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[5][0] - xRel[3][0];
  y1 = xRel[5][1] - xRel[3][1];
  z1 = xRel[5][2] - xRel[3][2];
  x2 = xRel[4][0] - xRel[3][0];
  y2 = xRel[4][1] - xRel[3][1];
  z2 = xRel[4][2] - xRel[3][2];
  x3 = xRel[0][0] - xRel[3][0];
  y3 = xRel[0][1] - xRel[3][1];
  z3 = xRel[0][2] - xRel[3][2];

  xc = xRelC[0] - xRel[3][0];
  yc = xRelC[1] - xRel[3][1];
  zc = xRelC[2] - xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = s;
    parCoor[1] = r;
    parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 1-0-4-5.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[0][0] - xRel[1][0];
  y1 = xRel[0][1] - xRel[1][1];
  z1 = xRel[0][2] - xRel[1][2];
  x2 = xRel[4][0] - xRel[1][0];
  y2 = xRel[4][1] - xRel[1][1];
  z2 = xRel[4][2] - xRel[1][2];
  x3 = xRel[5][0] - xRel[1][0];
  y3 = xRel[5][1] - xRel[1][1];
  z3 = xRel[5][2] - xRel[1][2];

  xc = xRelC[0] - xRel[1][0];
  yc = xRelC[1] - xRel[1][1];
  zc = xRelC[2] - xRel[1][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - r;
    parCoor[1] = t;
    parCoor[2] = s;
    return true;
  }

  /* Tetrahedron, 0-1-2-5.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[1][0];
  y1 = xRel[1][1];
  z1 = xRel[1][2];
  x2 = xRel[2][0];
  y2 = xRel[2][1];
  z2 = xRel[2][2];
  x3 = xRel[5][0];
  y3 = xRel[5][1];
  z3 = xRel[5][2];

  xc = xRelC[0];
  yc = xRelC[1];
  zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = r;
    parCoor[1] = s;
    parCoor[2] = t;
    return true;
  }

  /* The coordinate is in none of the six sub-tetrahedra of the prism. This
     implies that the point is not inside the true prism either and hence
     false is returned. */
  return false;
}

bool CADTElemClass::CoorInHexahedron(const unsigned long elemID, const su2double* coor, su2double* parCoor,
                                     su2double* weightsInterpol) const {
  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton = 1.e-10;

  /* Determine the indices of the eight vertices of the hexahedron,
     multiplied by nDim (which is 3). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0 + 3, i4 = i0 + 4, i5 = i0 + 5, i6 = i0 + 6, i7 = i0 + 7;
  i0 = nDim * elemConns[i0];
  i1 = nDim * elemConns[i1];
  i2 = nDim * elemConns[i2];
  i3 = nDim * elemConns[i3];
  i4 = nDim * elemConns[i4];
  i5 = nDim * elemConns[i5];
  i6 = nDim * elemConns[i6];
  i7 = nDim * elemConns[i7];

  /* Determine the coordinates relative to the vertex 0. */
  su2double xRel[8][3], xc[3];
  xc[0] = coor[0] - coorPoints[i0];
  xc[1] = coor[1] - coorPoints[i0 + 1];
  xc[2] = coor[2] - coorPoints[i0 + 2];

  xRel[0][0] = xRel[0][1] = xRel[0][2] = 0.0;

  xRel[1][0] = coorPoints[i1] - coorPoints[i0];
  xRel[1][1] = coorPoints[i1 + 1] - coorPoints[i0 + 1];
  xRel[1][2] = coorPoints[i1 + 2] - coorPoints[i0 + 2];

  xRel[2][0] = coorPoints[i2] - coorPoints[i0];
  xRel[2][1] = coorPoints[i2 + 1] - coorPoints[i0 + 1];
  xRel[2][2] = coorPoints[i2 + 2] - coorPoints[i0 + 2];

  xRel[3][0] = coorPoints[i3] - coorPoints[i0];
  xRel[3][1] = coorPoints[i3 + 1] - coorPoints[i0 + 1];
  xRel[3][2] = coorPoints[i3 + 2] - coorPoints[i0 + 2];

  xRel[4][0] = coorPoints[i4] - coorPoints[i0];
  xRel[4][1] = coorPoints[i4 + 1] - coorPoints[i0 + 1];
  xRel[4][2] = coorPoints[i4 + 2] - coorPoints[i0 + 2];

  xRel[5][0] = coorPoints[i5] - coorPoints[i0];
  xRel[5][1] = coorPoints[i5 + 1] - coorPoints[i0 + 1];
  xRel[5][2] = coorPoints[i5 + 2] - coorPoints[i0 + 2];

  xRel[6][0] = coorPoints[i6] - coorPoints[i0];
  xRel[6][1] = coorPoints[i6 + 1] - coorPoints[i0 + 1];
  xRel[6][2] = coorPoints[i6 + 2] - coorPoints[i0 + 2];

  xRel[7][0] = coorPoints[i7] - coorPoints[i0];
  xRel[7][1] = coorPoints[i7 + 1] - coorPoints[i0 + 1];
  xRel[7][2] = coorPoints[i7 + 2] - coorPoints[i0 + 2];

  /* Obtain an initial guess of the parametric coordinates by splitting the
     hexahedron into tetrahedra. If this approach is not successful, this means
     that the point is not inside the true hexahedron and false can be returned. */
  if (!InitialGuessContainmentHexahedron(xc, xRel, parCoor)) return false;

  /* The hexahedron is parametrized by X-X0 = (Xi-X0)*li, where the sum runs
     over i = 0..7, although i = 0 does not give a contribution. The Lagrangian
     interpolation functions for the hexahedron are given by:
     l0 = (1-r)*(1-s)*(1-t)/8, l1 = (1+r)*(1-s)*(1-t)/8,
     l2 = (1+r)*(1+s)*(1-t)/8, l3 = (1-r)*(1+s)*(1-t)/8,
     l4 = (1-r)*(1-s)*(1+t)/8, l5 = (1+r)*(1-s)*(1+t)/8,
     l6 = (1+r)*(1+s)*(1+t)/8, l7 = (1-r)*(1+s)*(1+t)/8.
     The boundaries are -1 <= r,s,t <= 1.
     As all coordinates are taken relative to vertex 0, X0 drops out and the
     nonlinear set of equations can be written as
     V0 - V1*r - V2*s - V3*t - V4*r*s - V5*r*t - V6*s*t - V7*r*s*t = 0, where
     V0 = xc - (x1+x2+x3+x4+x5+x6+x7)/8, V1 = (x1+x2-x3-x4+x5+x6-x7)/8,
     V2 =      (x2+x3-x1-x4-x5+x6+x7)/8, V3 = (x4+x5+x6+x7-x1-x2-x3)/8,
     V4 =      (x2+x4+x6-x1-x3-x5-x7)/8, V5 = (x3+x5+x6-x1-x2-x4-x7)/8,
     V6 =      (x1+x6+x7-x2-x3-x4-x5)/8, V7 = (x1+x3+x4+x6-x2-x5-x7)/8.
     First construct these vectors. */
  const su2double V0x =
      xc[0] - 0.125 * (xRel[1][0] + xRel[2][0] + xRel[3][0] + xRel[4][0] + xRel[5][0] + xRel[6][0] + xRel[7][0]);
  const su2double V0y =
      xc[1] - 0.125 * (xRel[1][1] + xRel[2][1] + xRel[3][1] + xRel[4][1] + xRel[5][1] + xRel[6][1] + xRel[7][1]);
  const su2double V0z =
      xc[2] - 0.125 * (xRel[1][2] + xRel[2][2] + xRel[3][2] + xRel[4][2] + xRel[5][2] + xRel[6][2] + xRel[7][2]);

  const su2double V1x =
      0.125 * (xRel[1][0] + xRel[2][0] - xRel[3][0] - xRel[4][0] + xRel[5][0] + xRel[6][0] - xRel[7][0]);
  const su2double V1y =
      0.125 * (xRel[1][1] + xRel[2][1] - xRel[3][1] - xRel[4][1] + xRel[5][1] + xRel[6][1] - xRel[7][1]);
  const su2double V1z =
      0.125 * (xRel[1][2] + xRel[2][2] - xRel[3][2] - xRel[4][2] + xRel[5][2] + xRel[6][2] - xRel[7][2]);

  const su2double V2x =
      0.125 * (xRel[2][0] + xRel[3][0] - xRel[1][0] - xRel[4][0] - xRel[5][0] + xRel[6][0] + xRel[7][0]);
  const su2double V2y =
      0.125 * (xRel[2][1] + xRel[3][1] - xRel[1][1] - xRel[4][1] - xRel[5][1] + xRel[6][1] + xRel[7][1]);
  const su2double V2z =
      0.125 * (xRel[2][2] + xRel[3][2] - xRel[1][2] - xRel[4][2] - xRel[5][2] + xRel[6][2] + xRel[7][2]);

  const su2double V3x =
      0.125 * (xRel[4][0] + xRel[5][0] + xRel[6][0] + xRel[7][0] - xRel[1][0] - xRel[2][0] - xRel[3][0]);
  const su2double V3y =
      0.125 * (xRel[4][1] + xRel[5][1] + xRel[6][1] + xRel[7][1] - xRel[1][1] - xRel[2][1] - xRel[3][1]);
  const su2double V3z =
      0.125 * (xRel[4][2] + xRel[5][2] + xRel[6][2] + xRel[7][2] - xRel[1][2] - xRel[2][2] - xRel[3][2]);

  const su2double V4x =
      0.125 * (xRel[2][0] + xRel[4][0] + xRel[6][0] - xRel[1][0] - xRel[3][0] - xRel[5][0] - xRel[7][0]);
  const su2double V4y =
      0.125 * (xRel[2][1] + xRel[4][1] + xRel[6][1] - xRel[1][1] - xRel[3][1] - xRel[5][1] - xRel[7][1]);
  const su2double V4z =
      0.125 * (xRel[2][2] + xRel[4][2] + xRel[6][2] - xRel[1][2] - xRel[3][2] - xRel[5][2] - xRel[7][2]);

  const su2double V5x =
      0.125 * (xRel[3][0] + xRel[5][0] + xRel[6][0] - xRel[1][0] - xRel[2][0] - xRel[4][0] - xRel[7][0]);
  const su2double V5y =
      0.125 * (xRel[3][1] + xRel[5][1] + xRel[6][1] - xRel[1][1] - xRel[2][1] - xRel[4][1] - xRel[7][1]);
  const su2double V5z =
      0.125 * (xRel[3][2] + xRel[5][2] + xRel[6][2] - xRel[1][2] - xRel[2][2] - xRel[4][2] - xRel[7][2]);

  const su2double V6x =
      0.125 * (xRel[1][0] + xRel[6][0] + xRel[7][0] - xRel[2][0] - xRel[3][0] - xRel[4][0] - xRel[5][0]);
  const su2double V6y =
      0.125 * (xRel[1][1] + xRel[6][1] + xRel[7][1] - xRel[2][1] - xRel[3][1] - xRel[4][1] - xRel[5][1]);
  const su2double V6z =
      0.125 * (xRel[1][2] + xRel[6][2] + xRel[7][2] - xRel[2][2] - xRel[3][2] - xRel[4][2] - xRel[5][2]);

  const su2double V7x =
      0.125 * (xRel[1][0] + xRel[3][0] + xRel[4][0] + xRel[6][0] - xRel[2][0] - xRel[5][0] - xRel[7][0]);
  const su2double V7y =
      0.125 * (xRel[1][1] + xRel[3][1] + xRel[4][1] + xRel[6][1] - xRel[2][1] - xRel[5][1] - xRel[7][1]);
  const su2double V7z =
      0.125 * (xRel[1][2] + xRel[3][2] + xRel[4][2] + xRel[6][2] - xRel[2][2] - xRel[5][2] - xRel[7][2]);

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for (itCount = 0; itCount < maxIt; ++itCount) {
    const su2double f0 = V0x - parCoor[0] * V1x - parCoor[1] * V2x - parCoor[2] * V3x - parCoor[0] * parCoor[1] * V4x -
                         parCoor[0] * parCoor[2] * V5x - parCoor[1] * parCoor[2] * V6x -
                         parCoor[0] * parCoor[1] * parCoor[2] * V7x;
    const su2double f1 = V0y - parCoor[0] * V1y - parCoor[1] * V2y - parCoor[2] * V3y - parCoor[0] * parCoor[1] * V4y -
                         parCoor[0] * parCoor[2] * V5y - parCoor[1] * parCoor[2] * V6y -
                         parCoor[0] * parCoor[1] * parCoor[2] * V7y;
    const su2double f2 = V0z - parCoor[0] * V1z - parCoor[1] * V2z - parCoor[2] * V3z - parCoor[0] * parCoor[1] * V4z -
                         parCoor[0] * parCoor[2] * V5z - parCoor[1] * parCoor[2] * V6z -
                         parCoor[0] * parCoor[1] * parCoor[2] * V7z;

    /* Compute the negative of the Jacobian matrix. */
    const su2double a00 = V1x + parCoor[1] * V4x + parCoor[2] * V5x + parCoor[1] * parCoor[2] * V7x;
    const su2double a01 = V2x + parCoor[0] * V4x + parCoor[2] * V6x + parCoor[0] * parCoor[2] * V7x;
    const su2double a02 = V3x + parCoor[0] * V5x + parCoor[1] * V6x + parCoor[0] * parCoor[1] * V7x;

    const su2double a10 = V1y + parCoor[1] * V4y + parCoor[2] * V5y + parCoor[1] * parCoor[2] * V7y;
    const su2double a11 = V2y + parCoor[0] * V4y + parCoor[2] * V6y + parCoor[0] * parCoor[2] * V7y;
    const su2double a12 = V3y + parCoor[0] * V5y + parCoor[1] * V6y + parCoor[0] * parCoor[1] * V7y;

    const su2double a20 = V1z + parCoor[1] * V4z + parCoor[2] * V5z + parCoor[1] * parCoor[2] * V7z;
    const su2double a21 = V2z + parCoor[0] * V4z + parCoor[2] * V6z + parCoor[0] * parCoor[2] * V7z;
    const su2double a22 = V3z + parCoor[0] * V5z + parCoor[1] * V6z + parCoor[0] * parCoor[1] * V7z;

    /* Compute the update of the parametric coordinates. */
    const su2double detInv = 1.0 / (a00 * a11 * a22 - a00 * a12 * a21 - a01 * a10 * a22 + a01 * a12 * a20 +
                                    a02 * a10 * a21 - a02 * a11 * a20);
    const su2double dr =
        detInv * (a01 * a12 * f2 - a01 * a22 * f1 - a02 * a11 * f2 + a02 * a21 * f1 + a11 * a22 * f0 - a12 * a21 * f0);
    const su2double ds =
        -detInv * (a00 * a12 * f2 - a00 * a22 * f1 - a02 * a10 * f2 + a02 * a20 * f1 + a10 * a22 * f0 - a12 * a20 * f0);
    const su2double dt =
        detInv * (a00 * a11 * f2 - a00 * a21 * f1 - a01 * a10 * f2 + a01 * a20 * f1 + a10 * a21 * f0 - a11 * a20 * f0);

    /* Update the parametric coordinates. Note that the negative of the
       Jacobian is used, so the update must be added. */
    parCoor[0] += dr;
    parCoor[1] += ds;
    parCoor[2] += dt;

    /* Check for convergence. */
    if (fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if (itCount == maxIt) SU2_MPI::Error("Newton did not converge", CURRENT_FUNCTION);

  /* Check if the point resides within the hexcahedron and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if ((parCoor[0] >= paramLowerBound) && (parCoor[0] <= paramUpperBound) && (parCoor[1] >= paramLowerBound) &&
      (parCoor[1] <= paramUpperBound) && (parCoor[2] >= paramLowerBound) && (parCoor[2] <= paramUpperBound)) {
    coorIsInside = true;

    const su2double omr = 0.5 * (1.0 - parCoor[0]), opr = 0.5 * (1.0 + parCoor[0]);
    const su2double oms = 0.5 * (1.0 - parCoor[1]), ops = 0.5 * (1.0 + parCoor[1]);
    const su2double omt = 0.5 * (1.0 - parCoor[2]), opt = 0.5 * (1.0 + parCoor[2]);

    weightsInterpol[0] = omr * oms * omt;
    weightsInterpol[1] = opr * oms * omt;
    weightsInterpol[2] = opr * ops * omt;
    weightsInterpol[3] = omr * ops * omt;
    weightsInterpol[4] = omr * oms * opt;
    weightsInterpol[5] = opr * oms * opt;
    weightsInterpol[6] = opr * ops * opt;
    weightsInterpol[7] = omr * ops * opt;
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::InitialGuessContainmentHexahedron(const su2double xRelC[3], const su2double xRel[8][3],
                                                      su2double* parCoor) const {
  /* Tetrahedron, 0-1-2-5.
     Create the coordinates of the tetrahedron and of the point. */
  su2double x1 = xRel[1][0], y1 = xRel[1][1], z1 = xRel[1][2];
  su2double x2 = xRel[2][0], y2 = xRel[2][1], z2 = xRel[2][2];
  su2double x3 = xRel[5][0], y3 = xRel[5][1], z3 = xRel[5][2];

  su2double xc = xRelC[0], yc = xRelC[1], zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  su2double detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  su2double r =
      detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  su2double s =
      -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  su2double t =
      detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = r;
    parCoor[1] = s;
    parCoor[2] = t;
    return true;
  }

  /* Tetrahedron, 4-7-5-0.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[7][0] - xRel[4][0];
  y1 = xRel[7][1] - xRel[4][1];
  z1 = xRel[7][2] - xRel[4][2];
  x2 = xRel[5][0] - xRel[4][0];
  y2 = xRel[5][1] - xRel[4][1];
  z2 = xRel[5][2] - xRel[4][2];
  x3 = xRel[0][0] - xRel[4][0];
  y3 = xRel[0][1] - xRel[4][1];
  z3 = xRel[0][2] - xRel[4][2];

  xc = xRelC[0] - xRel[4][0];
  yc = xRelC[1] - xRel[4][1];
  zc = xRelC[2] - xRel[4][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = s;
    parCoor[1] = r;
    parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 6-7-5-2.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[7][0] - xRel[6][0];
  y1 = xRel[7][1] - xRel[6][1];
  z1 = xRel[7][2] - xRel[6][2];
  x2 = xRel[5][0] - xRel[6][0];
  y2 = xRel[5][1] - xRel[6][1];
  z2 = xRel[5][2] - xRel[6][2];
  x3 = xRel[2][0] - xRel[6][0];
  y3 = xRel[2][1] - xRel[6][1];
  z3 = xRel[2][2] - xRel[6][2];

  xc = xRelC[0] - xRel[6][0];
  yc = xRelC[1] - xRel[6][1];
  zc = xRelC[2] - xRel[6][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - s;
    parCoor[1] = 1.0 - r;
    parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 3-0-2-7.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[0][0] - xRel[3][0];
  y1 = xRel[0][1] - xRel[3][1];
  z1 = xRel[0][2] - xRel[3][2];
  x2 = xRel[2][0] - xRel[3][0];
  y2 = xRel[2][1] - xRel[3][1];
  z2 = xRel[2][2] - xRel[3][2];
  x3 = xRel[7][0] - xRel[3][0];
  y3 = xRel[7][1] - xRel[3][1];
  z3 = xRel[7][2] - xRel[3][2];

  xc = xRelC[0] - xRel[3][0];
  yc = xRelC[1] - xRel[3][1];
  zc = xRelC[2] - xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - s;
    parCoor[1] = r;
    parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 0-5-2-7.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[5][0];
  y1 = xRel[5][1];
  z1 = xRel[5][2];
  x2 = xRel[2][0];
  y2 = xRel[2][1];
  z2 = xRel[2][2];
  x3 = xRel[7][0];
  y3 = xRel[7][1];
  z3 = xRel[7][2];

  xc = xRelC[0];
  yc = xRelC[1];
  zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = 1.0 + r + s;
    parCoor[1] = 1.0 + s + t;
    parCoor[2] = 1.0 + r + t;
    return true;
  }

  /* Tetrahedron, 0-1-3-4.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[1][0];
  y1 = xRel[1][1];
  z1 = xRel[1][2];
  x2 = xRel[3][0];
  y2 = xRel[3][1];
  z2 = xRel[3][2];
  x3 = xRel[4][0];
  y3 = xRel[4][1];
  z3 = xRel[4][2];

  xc = xRelC[0];
  yc = xRelC[1];
  zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = r;
    parCoor[1] = s;
    parCoor[2] = t;
    return true;
  }

  /* Tetrahedron, 7-6-4-3.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[6][0] - xRel[7][0];
  y1 = xRel[6][1] - xRel[7][1];
  z1 = xRel[6][2] - xRel[7][2];
  x2 = xRel[4][0] - xRel[7][0];
  y2 = xRel[4][1] - xRel[7][1];
  z2 = xRel[4][2] - xRel[7][2];
  x3 = xRel[3][0] - xRel[7][0];
  y3 = xRel[3][1] - xRel[7][1];
  z3 = xRel[3][2] - xRel[7][2];

  xc = xRelC[0] - xRel[7][0];
  yc = xRelC[1] - xRel[7][1];
  zc = xRelC[2] - xRel[7][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = r;
    parCoor[1] = 1.0 - s;
    parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 5-4-6-1.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[4][0] - xRel[5][0];
  y1 = xRel[4][1] - xRel[5][1];
  z1 = xRel[4][2] - xRel[5][2];
  x2 = xRel[6][0] - xRel[5][0];
  y2 = xRel[6][1] - xRel[5][1];
  z2 = xRel[6][2] - xRel[5][2];
  x3 = xRel[1][0] - xRel[5][0];
  y3 = xRel[1][1] - xRel[5][1];
  z3 = xRel[1][2] - xRel[5][2];

  xc = xRelC[0] - xRel[5][0];
  yc = xRelC[1] - xRel[5][1];
  zc = xRelC[2] - xRel[5][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - r;
    parCoor[1] = s;
    parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 2-3-1-6.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[3][0] - xRel[2][0];
  y1 = xRel[3][1] - xRel[2][1];
  z1 = xRel[3][2] - xRel[2][2];
  x2 = xRel[1][0] - xRel[2][0];
  y2 = xRel[1][1] - xRel[2][1];
  z2 = xRel[1][2] - xRel[2][2];
  x3 = xRel[6][0] - xRel[2][0];
  y3 = xRel[6][1] - xRel[2][1];
  z3 = xRel[6][2] - xRel[2][2];

  xc = xRelC[0] - xRel[2][0];
  yc = xRelC[1] - xRel[2][1];
  zc = xRelC[2] - xRel[2][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - r;
    parCoor[1] = 1.0 - s;
    parCoor[2] = t;
    return true;
  }

  /* Tetrahedron, 3-4-1-6.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[4][0] - xRel[3][0];
  y1 = xRel[4][1] - xRel[3][1];
  z1 = xRel[4][2] - xRel[3][2];
  x2 = xRel[1][0] - xRel[3][0];
  y2 = xRel[1][1] - xRel[3][1];
  z2 = xRel[1][2] - xRel[3][2];
  x3 = xRel[6][0] - xRel[3][0];
  y3 = xRel[6][1] - xRel[3][1];
  z3 = xRel[6][2] - xRel[3][2];

  xc = xRelC[0] - xRel[3][0];
  yc = xRelC[1] - xRel[3][1];
  zc = xRelC[2] - xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0 / (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1);
  r = detInv * (x2 * y3 * zc - x2 * yc * z3 - x3 * y2 * zc + x3 * yc * z2 + xc * y2 * z3 - xc * y3 * z2) - 1.0;
  s = -detInv * (x1 * y3 * zc - x1 * yc * z3 - x3 * y1 * zc + x3 * yc * z1 + xc * y1 * z3 - xc * y3 * z1) - 1.0;
  t = detInv * (x1 * y2 * zc - x1 * yc * z2 - x2 * y1 * zc + x2 * yc * z1 + xc * y1 * z2 - xc * y2 * z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) && ((r + s + t) <= paramLowerBound)) {
    parCoor[0] = 1.0 + s + t;
    parCoor[1] = -1.0 - r - s;
    parCoor[2] = 1.0 + r + t;
    return true;
  }

  /* The coordinate is in none of the ten sub-tetrahedra of the hexahedron.
     This implies that the point is not inside the true hexahedron either
     and hence false is returned. */
  return false;
}

void CADTElemClass::Dist2ToLine(const unsigned long i0, const unsigned long i1, const su2double* coor,
                                su2double& dist2Line) const {
  /*--- The line is parametrized by X = X0 + (r+1)*(X1-X0)/2, -1 <= r <= 1.
        As a consequence the minimum distance is found where the expression
        |V0 - r*V1| has a minimum, where the vectors V0 and V1 are defined
        as: V0 = coor - (X1+X0)/2, V1 = (X1-X0)/2. First construct the
        vectors V0 and V1. ---*/
  su2double V0[3], V1[3];
  for (unsigned short k = 0; k < nDim; ++k) {
    V0[k] = coor[k] - 0.5 * (coorPoints[i1 + k] + coorPoints[i0 + k]);
    V1[k] = 0.5 * (coorPoints[i1 + k] - coorPoints[i0 + k]);
  }

  /*--- Determine the value of r, for which the minimum occurs. This is the
        ratio of the dot product of V0.V1 and V1.V1. Make sure that r is
        limited between -1 and 1 to make sure that the projection is inside
        the line element. ---*/
  su2double dotV0V1 = 0.0, dotV1V1 = 0.0;
  for (unsigned short k = 0; k < nDim; ++k) {
    dotV0V1 += V0[k] * V1[k];
    dotV1V1 += V1[k] * V1[k];
  }
  su2double r = dotV0V1 / dotV1V1;
  r = max(-1.0, min(1.0, r));

  /*--- Determine the minimum distance squared. ---*/
  dist2Line = 0.0;
  for (unsigned short k = 0; k < nDim; ++k) {
    const su2double ds = V0[k] - r * V1[k];
    dist2Line += ds * ds;
  }
}

bool CADTElemClass::Dist2ToTriangle(const unsigned long i0, const unsigned long i1, const unsigned long i2,
                                    const su2double* coor, su2double& dist2Tria, su2double& r, su2double& s) const {
  constexpr unsigned short nDim = 3;  // boundary triangles only exist in 3D

  /*--- The triangle is parametrized by X = X0 + (r+1)*(X1-X0)/2 + (s+1)*(X2-X0)/2,
        r, s >= -1, r+s <= 0. As a consequence the minimum distance is found where
        the expression |V0 - r*V1 - s*V2| has a minimum, where the vectors V0, V1
        and V2 are defined as: V0 = coor - (X1+X2)/2, V1 = (X1-X0)/2, V2 = (X2-X0)/2.
        First construct the vectors V0, V1 and V2. ---*/
  su2double V0[3], V1[3], V2[3];
  for (unsigned short k = 0; k < nDim; ++k) {
    V0[k] = coor[k] - 0.5 * (coorPoints[i1 + k] + coorPoints[i2 + k]);
    V1[k] = 0.5 * (coorPoints[i1 + k] - coorPoints[i0 + k]);
    V2[k] = 0.5 * (coorPoints[i2 + k] - coorPoints[i0 + k]);
  }

  /*--- For the values of r and s where the minimum occurs the dot products V0.V1,
        V0.V2, V1.V1, V1.V2 and V2.V2 are needed. Compute these values. ---*/
  su2double dotV0V1 = 0.0, dotV0V2 = 0, dotV1V1 = 0.0, dotV1V2 = 0.0, dotV2V2 = 0.0;
  for (unsigned short k = 0; k < nDim; ++k) {
    dotV0V1 += V0[k] * V1[k];
    dotV0V2 += V0[k] * V2[k];
    dotV1V1 += V1[k] * V1[k];
    dotV1V2 += V1[k] * V2[k];
    dotV2V2 += V2[k] * V2[k];
  }

  /*--- Solve the linear system for r and s. ---*/
  const su2double detInv = 1.0 / (dotV1V1 * dotV2V2 - dotV1V2 * dotV1V2);

  r = detInv * (dotV0V1 * dotV2V2 - dotV0V2 * dotV1V2);
  s = detInv * (dotV0V2 * dotV1V1 - dotV0V1 * dotV1V2);

  /*--- Check if the projection is inside the triangle. ---*/
  if ((r >= paramLowerBound) && (s >= paramLowerBound) && ((r + s) <= tolInsideElem)) {
    /*--- The projection of the coordinate is inside the triangle. Compute the
          minimum distance squared and return true. ---*/
    dist2Tria = 0.0;
    for (unsigned short k = 0; k < nDim; ++k) {
      const su2double ds = V0[k] - r * V1[k] - s * V2[k];
      dist2Tria += ds * ds;
    }

    return true;
  }

  /*--- The projection of the coordinate is outside the triangle.
        Return false. ---*/
  return false;
}

bool CADTElemClass::Dist2ToQuadrilateral(const unsigned long i0, const unsigned long i1, const unsigned long i2,
                                         const unsigned long i3, const su2double* coor, su2double& r, su2double& s,
                                         su2double& dist2Quad) const {
  constexpr unsigned short nDim = 3;  // boundary quadrilaterals only exist in 3D

  /* Definition of the maximum number of iterations in the iterative solver
     and the tolerance level. */
  const unsigned short maxIt = 10;
  const su2double tolIt = 1.e-10;

  /*--- The quadrilateral is parametrized by
        X = {(1-r)*(1-s)*X0 + (1+r)*(1-s)*X1 + (1+r)*(1+s)*X2 + (1-r)*(1+s)*X3}/4,
        -1 <= r,s <= 1. As a consequence the minimum distance is found where the
        expression |V0 - r*V1 - s*V2 - r*s*V3| has a minimum, where the vectors
        V0, V1, V2 and V3 are defined as: V0 = coor - (X0+X1+X2+X3)/4,
        V1 = (X1+X2-X0-X3)/4, V2 = (X2+X3-X0-X1), V3 = (X0+X2-X1-X3)/4. First
        construct the vectors V0, V1, V2 and V3. ---*/
  su2double V0[3], V1[3], V2[3], V3[3];
  for (unsigned short k = 0; k < nDim; ++k) {
    V0[k] = coor[k] - 0.25 * (coorPoints[i0 + k] + coorPoints[i1 + k] + coorPoints[i2 + k] + coorPoints[i3 + k]);
    V1[k] = 0.25 * (coorPoints[i1 + k] + coorPoints[i2 + k] - coorPoints[i0 + k] - coorPoints[i3 + k]);
    V2[k] = 0.25 * (coorPoints[i2 + k] + coorPoints[i3 + k] - coorPoints[i0 + k] - coorPoints[i1 + k]);
    V3[k] = 0.25 * (coorPoints[i0 + k] + coorPoints[i2 + k] - coorPoints[i1 + k] - coorPoints[i3 + k]);
  }

  /*--- The minimization problem results in the following
        set of nonlinear equations.
        (V0 - r*V1 - s*V2 - r*s*V3).(-V1 - s*V3) = 0
        (V0 - r*V1 - s*V2 - r*s*V3).(-V2 - r*V3) = 0.
        These equations are the gradients of the distance function squared w.r.t.
        the parametric coordinates r and s. The coefficients for these equations
        involve a couple of dot products, which are computed first.   ---*/
  su2double V0V1 = 0.0, V0V2 = 0.0, V0V3 = 0.0, V1V1 = 0.0, V1V2 = 0.0, V1V3 = 0.0, V2V2 = 0.0, V2V3 = 0.0, V3V3 = 0.0;
  for (unsigned short k = 0; k < nDim; ++k) {
    V0V1 += V0[k] * V1[k];
    V0V2 += V0[k] * V2[k];
    V0V3 += V0[k] * V3[k];
    V1V1 += V1[k] * V1[k];
    V1V2 += V1[k] * V2[k];
    V1V3 += V1[k] * V3[k];
    V2V2 += V2[k] * V2[k];
    V2V3 += V2[k] * V3[k];
    V3V3 += V3[k] * V3[k];
  }

  /*--- The first equation given above is linear in r and can be used to
        express r as a function of s, i.e.
        r = (V0V1 + (V0V3-V1V2) s - V2V3 s^2)/(V1V1 + 2 V1V3 s + V3V3 s^2).
        This is substituted in the second equation and after some rearranging
        a fifth order equation for s is obtained, i.e.
        a5 s^5 + a4 s^4 + a3 s^3 + a2 s^2 + a1 s + a0 = 0.
        The coefficients for this equation are computed below. ---*/
  const su2double t1 = V0V1 * V0V1;
  const su2double t3 = V0V3 - V1V2;
  const su2double t6 = V1V1 * V1V1;
  const su2double t11 = V0V1 * V2V3 * 2.0;
  const su2double t14 = V0V3 * V0V3;
  const su2double t16 = V0V3 * V1V2 * 2.0;
  const su2double t17 = V1V2 * V1V2;
  const su2double t22 = V1V3 * V1V3;
  const su2double t34 = V3V3 * V0V2;
  const su2double t50 = V2V3 * V2V3;
  const su2double t55 = V3V3 * V3V3;

  su2double a0 = -V1V1 * t3 * V0V1 - t6 * V0V2 + V1V3 * t1;
  su2double a1 = t6 * V2V2 + V1V1 * (-4.0 * V0V2 * V1V3 + t11 - t14 + t16 - t17) + t1 * V3V3;
  su2double a2 = -4.0 * t22 * V0V2 + V1V3 * (4.0 * V1V1 * V2V2 + t11 - t14 + t16 - t17) +
                 (V0V1 * V3V3 + 3.0 * V2V3 * V1V1) * V0V3 + V1V1 * (-3.0 * V1V2 * V2V3 - 2.0 * t34) -
                 V0V1 * V1V2 * V3V3;
  su2double a3 = 4.0 * V2V2 * t22 + V1V3 * (4.0 * V2V3 * t3 - 4.0 * t34) + 2.0 * (V2V2 * V3V3 - t50) * V1V1;
  su2double a4 = -t55 * V0V2 + V3V3 * (4.0 * V2V2 * V1V3 + V2V3 * t3) - 3.0 * t50 * V1V3;
  su2double a5 = t55 * V2V2 - V3V3 * t50;

  /* Determine the maximum of these coefficients and scale them. */
  su2double scaleFact = max(fabs(a0), fabs(a1));
  scaleFact = max(scaleFact, fabs(a2));
  scaleFact = max(scaleFact, fabs(a3));
  scaleFact = max(scaleFact, fabs(a4));
  scaleFact = max(scaleFact, fabs(a5));

  scaleFact = 1.0 / scaleFact;

  a0 *= scaleFact;
  a1 *= scaleFact;
  a2 *= scaleFact;
  a3 *= scaleFact;
  a4 *= scaleFact;
  a5 *= scaleFact;

  /* The coefficients that occur in the derivative of f, i.e.
     b4 s^4 + b3 s^3 + b2 s^2 + b1 s + b0. */
  const su2double b4 = 5.0 * a5;
  const su2double b3 = 4.0 * a4;
  const su2double b2 = 3.0 * a3;
  const su2double b1 = 2.0 * a2;
  const su2double b0 = a1;

  /* Initial guess for s. */
  s = 0.0;

  /*--- Newtons algorithm to solve the nonlinear equation
        a5 s^5 + a4 s^4 + a3 s^3 + a2 s^2 + a1 s + a0 = 0. ---*/
  unsigned short itCount;
  for (itCount = 0; itCount < maxIt; ++itCount) {
    /* Store the old value of s for determining the actual update. */
    const su2double sOld = s;

    /* Compute the values of s2 to s5. */
    const su2double s2 = s * s;
    const su2double s3 = s2 * s;
    const su2double s4 = s3 * s;
    const su2double s5 = s4 * s;

    /* Compute the value of the function and its first derivative. */
    const su2double f = a5 * s5 + a4 * s4 + a3 * s3 + a2 * s2 + a1 * s + a0;
    const su2double df = b4 * s4 + b3 * s3 + b2 * s2 + b1 * s + b0;

    /* Compute the value of the update and the new value of s. */
    su2double ds = f / df;
    s -= ds;

    /* Clipping, such that s is bounded to a bit outside the quadrilateral. */
    s = max(s, -1.5);
    s = min(s, 1.5);

    /* Set the actual value of ds. */
    ds = sOld - s;

    /*--- Check for convergence. ---*/
    if (fabs(ds) <= tolIt) break;
  }

  /* Check if Newtons algorithm did not converge. */
  if (itCount == maxIt) {
    /* Newtons algorithm did not converge. The most likely reason is that there is
       no root for s on the interval paramLowerBound to paramUpperBound. Check this.
       First compute the function value for paramLowerBound. */
    s = paramLowerBound;

    su2double s2 = s * s;
    su2double s3 = s2 * s;
    su2double s4 = s3 * s;
    su2double s5 = s4 * s;

    su2double f = a5 * s5 + a4 * s4 + a3 * s3 + a2 * s2 + a1 * s + a0;

    /* Loop over a number of sampling points to check if the sign of f changes.
       If it does, this means that there is a root between s-ds and s. Hence
       a break from this loop can be made. */
    const su2double ds = (paramUpperBound - paramLowerBound) * 0.1;
    for (itCount = 1; itCount <= 10; ++itCount) {
      const su2double fOld = f;
      s += ds;

      s2 = s * s;
      s3 = s2 * s;
      s4 = s3 * s;
      s5 = s4 * s;
      f = a5 * s5 + a4 * s4 + a3 * s3 + a2 * s2 + a1 * s + a0;

      if (f * fOld <= 0.0) break;
    }

    /* If the end of the loop is reached this means that there is no root for s
       between the lower and upper bound. Return false. */
    if (itCount > 10) return false;

    /* Newtons algorithm did not converge, but there is a root. Do a crude approximation
       of the root using bisection. */
    s -= 0.5 * ds;
  }

  /* Check if s is inside the quadrilateral. If not, return false. */
  if (s < paramLowerBound || s > paramUpperBound) return false;

  /* Compute the corresponding value of r and check if it is inside the
     quadrilateral. If not return false. */
  const su2double s2 = s * s;

  r = (V0V1 + (V0V3 - V1V2) * s - V2V3 * s2) / (V1V1 + 2.0 * V1V3 * s + V3V3 * s2);
  if (r < paramLowerBound || r > paramUpperBound) return false;

  /*--- The projection is inside the quadrilateral. Determine the minimum distance
        squared and return true to indicate that the projection is inside. ---*/
  dist2Quad = 0.0;
  for (unsigned short k = 0; k < nDim; ++k) {
    const su2double ds = V0[k] - r * V1[k] - s * V2[k] - r * s * V3[k];
    dist2Quad += ds * ds;
  }

  return true;
}
