/*!
 * \file fem_geometry_structure.cpp
 * \brief Functions for creating the primal grid for the FEM solver.
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

#include "../../include/fem/fem_geometry_structure.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridFEM.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"
#include "../../include/adt/CADTElemClass.hpp"
#include "../../include/adt/CADTPointsOnlyClass.hpp"

/* Prototypes for Lapack functions, if MKL or LAPACK is used. */
#if defined(HAVE_MKL) || defined(HAVE_LAPACK)
extern "C" void dpotrf_(char*, int*, passivedouble*, int*, int*);
extern "C" void dpotri_(char*, int*, passivedouble*, int*, int*);
#endif

bool CLong3T::operator<(const CLong3T& other) const {
  if (long0 != other.long0) return (long0 < other.long0);
  if (long1 != other.long1) return (long1 < other.long1);
  if (long2 != other.long2) return (long2 < other.long2);

  return false;
}

CReorderElements::CReorderElements(const unsigned long val_GlobalElemID, const unsigned short val_TimeLevel,
                                   const bool val_CommSolution, const unsigned short val_VTK_Type,
                                   const unsigned short val_nPolySol, const bool val_JacConstant) {
  /* Copy the global elment ID, time level and whether or not this element
     must be communicated. */
  globalElemID = val_GlobalElemID;
  timeLevel = val_TimeLevel;
  commSolution = val_CommSolution;

  /* Create the element type used in this class, which stores information of
     the VTK type, the polynomial degree of the solution and whether or not the
     Jacobian of the transformation is constant. As it is possible that the
     polynomial degree of the solution is zero, this convention is different
     from the convention used in the SU2 grid file. */
  elemType = val_VTK_Type + 100 * val_nPolySol;
  if (!val_JacConstant) elemType += 50;
}

bool CReorderElements::operator<(const CReorderElements& other) const {
  /* Elements with the lowest time level are stored first. */
  if (timeLevel != other.timeLevel) return timeLevel < other.timeLevel;

  /* Next comparison is whether or not the element must communicate its
     solution data to other ranks. Elements which do not need to do this
     are stored first. */
  if (commSolution != other.commSolution) return other.commSolution;

  /* Elements of the same element type must be stored as contiguously as
     possible to allow for the simultaneous treatment of elements in the
     matrix multiplications. */
  if (elemType != other.elemType) return elemType < other.elemType;

  /* The final comparison is based on the global element ID. */
  return globalElemID < other.globalElemID;
}

bool CSortFaces::operator()(const CFaceOfElement& f0, const CFaceOfElement& f1) {
  /*--- Comparison in case both faces are boundary faces. ---*/
  if (f0.faceIndicator >= 0 && f1.faceIndicator >= 0) {
    /* Both faces are boundary faces. The first comparison is the boundary
       marker, which is stored in faceIndicator. */
    if (f0.faceIndicator != f1.faceIndicator) return f0.faceIndicator < f1.faceIndicator;

    /* Both faces belong to the same boundary marker. The second comparison is
       based on the on the local volume ID's of the adjacent elements. As the
       volumes are sorted according to the time levels for time accurate
       local stepping, there is no need to do this check seperately here. */
    unsigned long ind0 = f0.elemID0 < nVolElemTot ? f0.elemID0 : f0.elemID1;
    unsigned long ind1 = f1.elemID0 < nVolElemTot ? f1.elemID0 : f1.elemID1;

    return ind0 < ind1;
  }

  /*--- Comparison in case both faces are internal faces. ---*/
  if (f0.faceIndicator == -1 && f1.faceIndicator == -1) {
    /* Both faces are internal faces. First determine the minimum and maximum
       ID of its adjacent elements.  */
    unsigned long elemIDMin0 = min(f0.elemID0, f0.elemID1);
    unsigned long elemIDMax0 = max(f0.elemID0, f0.elemID1);

    unsigned long elemIDMin1 = min(f1.elemID0, f1.elemID1);
    unsigned long elemIDMax1 = max(f1.elemID0, f1.elemID1);

    /* Determine the situation. */
    if (elemIDMax0 < nVolElemTot && elemIDMax1 < nVolElemTot) {
      /* Both faces are matching internal faces. Determine whether or not these
         faces are local faces, i.e. faces between locally owned elements. */
      const bool face0IsLocal = elemIDMax0 < nVolElemOwned;
      const bool face1IsLocal = elemIDMax1 < nVolElemOwned;

      /* Check if both faces have the same status, i.e. either local or
         not local. */
      if (face0IsLocal == face1IsLocal) {
        /* Both faces are either local or not local. Determine the time level
           of the faces, which is the minimum value of the adjacent volume
           elements. */
        const unsigned short timeLevel0 = min(volElem[elemIDMin0].timeLevel, volElem[elemIDMax0].timeLevel);
        const unsigned short timeLevel1 = min(volElem[elemIDMin1].timeLevel, volElem[elemIDMax1].timeLevel);

        /* Internal faces with the same status are first sorted according to
           their time level. Faces with the smallest time level are numbered
           first. Note this is only relevant for time accurate local time
           stepping. */
        if (timeLevel0 != timeLevel1) return timeLevel0 < timeLevel1;

        /* The faces belong to the same time level. They are sorted according
           to their element ID's in order to increase cache performance. */
        if (elemIDMin0 != elemIDMin1) return elemIDMin0 < elemIDMin1;
        return elemIDMax0 < elemIDMax1;
      } /* One face is a local face and the other is not. Make sure that
   the local faces are numbered first. */
      return face0IsLocal;

    } else if (elemIDMax0 >= nVolElemTot && elemIDMax1 >= nVolElemTot) {
      /* Both faces are non-matching internal faces. Sort them according to
         their relevant element ID. The time level is not taken into account
         yet, because non-matching faces are not possible at the moment with
         time accurate local time stepping. */
      return elemIDMin0 < elemIDMin1;
    } else {
      /* One face is a matching internal face and the other face is a
         non-matching internal face. Make sure that the non-matching face
         is numbered after the matching face. This is accomplished by comparing
         the maximum element ID's. */
      return elemIDMax0 < elemIDMax1;
    }
  }

  /*--- One face is a boundary face and the other face is an internal face.
        Make sure that the boundary face is numbered first. This can be
        accomplished by using the > operator for faceIndicator. ---*/
  return f0.faceIndicator > f1.faceIndicator;
}

bool CSortBoundaryFaces::operator()(const CSurfaceElementFEM& f0, const CSurfaceElementFEM& f1) {
  /* First sorting criterion is the index of the standard element. The
     boundary faces should be sorted per standard element. Note that the
     time level is not taken into account here, because it is assumed that
     the surface elements to be sorted belong to one time level. */
  if (f0.indStandardElement != f1.indStandardElement) return f0.indStandardElement < f1.indStandardElement;

  /* The standard elements are the same. The second criterion is the
     corresponding volume IDs of the surface elements. */
  return f0.volElemID < f1.volElemID;
}

bool CPointFEM::operator<(const CPointFEM& other) const {
  if (periodIndexToDonor != other.periodIndexToDonor) return periodIndexToDonor < other.periodIndexToDonor;
  return globalID < other.globalID;
}

bool CPointFEM::operator==(const CPointFEM& other) const {
  return (globalID == other.globalID && periodIndexToDonor == other.periodIndexToDonor);
}

void CVolumeElementFEM::GetCornerPointsAllFaces(unsigned short& numFaces, unsigned short nPointsPerFace[],
                                                unsigned long faceConn[6][4]) {
  /*--- Get the corner connectivities of the faces, local to the element. ---*/
  CPrimalGridFEM::GetLocalCornerPointsAllFaces(VTK_Type, nPolyGrid, nDOFsGrid, numFaces, nPointsPerFace, faceConn);

  /*--- Convert the local values of faceConn to global values. ---*/
  for (unsigned short i = 0; i < numFaces; ++i) {
    for (unsigned short j = 0; j < nPointsPerFace[i]; ++j) {
      unsigned long nn = faceConn[i][j];
      faceConn[i][j] = nodeIDsGrid[nn];
    }
  }
}

bool CInternalFaceElementFEM::operator<(const CInternalFaceElementFEM& other) const {
  /* First comparison is the standard element, such that elements with
     the same same standard elements are grouped together. */
  if (indStandardElement != other.indStandardElement) return indStandardElement < other.indStandardElement;

  /* The second comparison is the element ID on side 0. */
  if (elemID0 != other.elemID0) return elemID0 < other.elemID0;

  /* The final comparison is the element ID on side 1. */
  return elemID1 < other.elemID1;
}

void CSurfaceElementFEM::GetCornerPointsFace(unsigned short& nPointsPerFace, unsigned long faceConn[]) {
  /*--- Get the corner connectivities of the face, local to the element. ---*/
  CPrimalGridBoundFEM::GetLocalCornerPointsFace(VTK_Type, nPolyGrid, nDOFsGrid, nPointsPerFace, faceConn);

  /*--- Convert the local values of faceConn to global values. ---*/
  for (unsigned short j = 0; j < nPointsPerFace; ++j) {
    unsigned long nn = faceConn[j];
    faceConn[j] = nodeIDsGrid[nn];
  }
}

CMeshFEM::CMeshFEM(CGeometry* geometry, CConfig* config) {
  /*--- Allocate the memory for blasFunctions. ---*/
  blasFunctions = new CBlasStructure;

  /*--- The new FEM mesh class has the same problem dimension/zone. ---*/
  nDim = geometry->GetnDim();
  nZone = geometry->GetnZone();

  /*--- Determine the number of variables stored per DOF.              ---*/
  /*--- This is specifically taylored towards perfect gas computations ---*/
  /*--- and must be generalized later. As this number is also needed   ---*/
  /*--- in the solver, it is worthwhile to create a separate function  ---*/
  /*--- for this purpose.                                              ---*/
  MAIN_SOLVER Kind_Solver = config->GetKind_Solver();
  const bool compressible = (Kind_Solver == MAIN_SOLVER::FEM_EULER) ||
                            (Kind_Solver == MAIN_SOLVER::FEM_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::FEM_RANS) ||
                            (Kind_Solver == MAIN_SOLVER::FEM_LES);

  unsigned short nVar;
  if (compressible)
    nVar = nDim + 2;
  else
    nVar = nDim + 1;

  /*--- Determine a mapping from the global point ID to the local index
        of the points.            ---*/
  map<unsigned long, unsigned long> globalPointIDToLocalInd;
  for (unsigned long i = 0; i < geometry->GetnPoint(); ++i)
    globalPointIDToLocalInd[geometry->nodes->GetGlobalIndex(i)] = i;

  /*----------------------------------------------------------------------------*/
  /*--- Step 1: Communicate the elements and the boundary elements to the    ---*/
  /*---         ranks where they will be stored during the computation.      ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Determine the ranks to which I have to send my elements. ---*/
  vector<int> sendToRank(size, 0);

  for (unsigned long i = 0; i < geometry->GetnElem(); ++i) {
    sendToRank[geometry->elem[i]->GetColor()] = 1;
  }

  map<int, int> rankToIndCommBuf;
  for (int i = 0; i < size; ++i) {
    if (sendToRank[i]) {
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  /*--- Definition of the communication buffers, used to send the element data
        to the correct ranks.                ---*/
  int nRankSend = (int)rankToIndCommBuf.size();
  vector<vector<short> > shortSendBuf(nRankSend, vector<short>(0));
  vector<vector<long> > longSendBuf(nRankSend, vector<long>(0));
  vector<vector<su2double> > doubleSendBuf(nRankSend, vector<su2double>(0));

  /*--- The first element of longSendBuf will contain the number of elements, which
        are stored in the communication buffers. Initialize this value to 0. ---*/
  for (int i = 0; i < nRankSend; ++i) longSendBuf[i].push_back(0);

  /*--- Determine the number of ranks, from which this rank will receive elements. ---*/
  int nRankRecv = nRankSend;

#ifdef HAVE_MPI
  vector<int> sizeRecv(size, 1);

  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- Loop over the local elements to fill the communication buffers with element data. ---*/
  for (unsigned long i = 0; i < geometry->GetnElem(); ++i) {
    int ind = (int)geometry->elem[i]->GetColor();
    map<int, int>::const_iterator MI = rankToIndCommBuf.find(ind);
    ind = MI->second;

    ++longSendBuf[ind][0]; /* The number of elements in the buffers must be incremented. */

    shortSendBuf[ind].push_back(geometry->elem[i]->GetVTK_Type());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNPolyGrid());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNPolySol());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNDOFsGrid());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNDOFsSol());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetnFaces());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetTimeLevel());
    shortSendBuf[ind].push_back((short)geometry->elem[i]->GetJacobianConsideredConstant());

    longSendBuf[ind].push_back(geometry->elem[i]->GetGlobalElemID());
    longSendBuf[ind].push_back(geometry->elem[i]->GetGlobalOffsetDOFsSol());

    for (unsigned short j = 0; j < geometry->elem[i]->GetNDOFsGrid(); ++j)
      longSendBuf[ind].push_back(geometry->elem[i]->GetNode(j));

    for (unsigned short j = 0; j < geometry->elem[i]->GetnFaces(); ++j)
      longSendBuf[ind].push_back(geometry->elem[i]->GetNeighbor_Elements(j));

    for (unsigned short j = 0; j < geometry->elem[i]->GetnFaces(); ++j) {
      shortSendBuf[ind].push_back(geometry->elem[i]->GetPeriodicIndex(j));
      shortSendBuf[ind].push_back((short)geometry->elem[i]->GetJacobianConstantFace(j));
      shortSendBuf[ind].push_back((short)geometry->elem[i]->GetOwnerFace(j));
    }

    doubleSendBuf[ind].push_back(geometry->elem[i]->GetLengthScale());
  }

  /*--- Determine for each rank to which I have to send elements the data of
        the corresponding nodes.   ---*/
  for (int i = 0; i < nRankSend; ++i) {
    /*--- Determine the vector with node IDs in the connectivity
          of the elements for this rank.   ---*/
    vector<long> nodeIDs;

    unsigned long indL = 3;
    unsigned long indS = 3;
    for (long j = 0; j < longSendBuf[i][0]; ++j) {
      short nDOFsGrid = shortSendBuf[i][indS], nFaces = shortSendBuf[i][indS + 2];
      indS += 3 * nFaces + 8;

      for (short k = 0; k < nDOFsGrid; ++k, ++indL) nodeIDs.push_back(longSendBuf[i][indL]);
      indL += nFaces + 2;
    }

    /*--- Sort nodeIDs in increasing order and remove the double entities. ---*/
    sort(nodeIDs.begin(), nodeIDs.end());
    auto lastNodeID = unique(nodeIDs.begin(), nodeIDs.end());
    nodeIDs.erase(lastNodeID, nodeIDs.end());

    /*--- Add the number of node IDs and the node IDs itself to longSendBuf[i]. ---*/
    longSendBuf[i].push_back(nodeIDs.size());
    longSendBuf[i].insert(longSendBuf[i].end(), nodeIDs.begin(), nodeIDs.end());

    /*--- Copy the coordinates to doubleSendBuf. ---*/
    for (unsigned long j = 0; j < nodeIDs.size(); ++j) {
      map<unsigned long, unsigned long>::const_iterator LMI;
      LMI = globalPointIDToLocalInd.find(nodeIDs[j]);

      if (LMI == globalPointIDToLocalInd.end()) SU2_MPI::Error("Entry not found in map", CURRENT_FUNCTION);

      unsigned long ind = LMI->second;
      for (unsigned short l = 0; l < nDim; ++l) doubleSendBuf[i].push_back(geometry->nodes->GetCoord(ind, l));
    }
  }

  /*--- Loop over the boundaries to send the boundary data to the appropriate rank. ---*/
  nMarker = geometry->GetnMarker();
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    /* Store the current indices in the longSendBuf, which are used to store the
       number of boundary elements sent to this rank. Initialize this value to 0. */
    vector<long> indLongBuf(nRankSend);
    for (int i = 0; i < nRankSend; ++i) {
      indLongBuf[i] = longSendBuf[i].size();
      longSendBuf[i].push_back(0);
    }

    /* Loop over the local boundary elements in geometry for this marker. */
    for (unsigned long i = 0; i < geometry->GetnElem_Bound(iMarker); ++i) {
      /* Determine the local ID of the corresponding domain element. */
      unsigned long elemID = geometry->bound[iMarker][i]->GetDomainElement() - geometry->beg_node[rank];

      /* Determine to which rank this boundary element must be sent.
         That is the same as its corresponding domain element.
         Update the corresponding index in longSendBuf. */
      int ind = (int)geometry->elem[elemID]->GetColor();
      map<int, int>::const_iterator MI = rankToIndCommBuf.find(ind);
      ind = MI->second;

      ++longSendBuf[ind][indLongBuf[ind]];

      /* Get the donor information for the wall function treatment. */
      const unsigned short nDonors = geometry->bound[iMarker][i]->GetNDonorsWallFunctions();
      const unsigned long* donors = geometry->bound[iMarker][i]->GetDonorsWallFunctions();

      /* Store the data for this boundary element in the communication buffers. */
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetVTK_Type());
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNPolyGrid());
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNDOFsGrid());
      shortSendBuf[ind].push_back(nDonors);

      longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetDomainElement());
      longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetGlobalElemID());

      for (unsigned short j = 0; j < geometry->bound[iMarker][i]->GetNDOFsGrid(); ++j)
        longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNode(j));

      for (unsigned short j = 0; j < nDonors; ++j) longSendBuf[ind].push_back(donors[j]);
    }
  }

  /*--- Definition of the communication buffers, used to receive
        the element data from the other correct ranks.        ---*/
  vector<vector<short> > shortRecvBuf(nRankRecv, vector<short>(0));
  vector<vector<long> > longRecvBuf(nRankRecv, vector<long>(0));
  vector<vector<su2double> > doubleRecvBuf(nRankRecv, vector<su2double>(0));

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/
  map<int, int>::const_iterator MI;

#ifdef HAVE_MPI

  /*--- Parallel mode. Send all the data using non-blocking sends. ---*/
  vector<SU2_MPI::Request> commReqs(3 * nRankSend);
  MI = rankToIndCommBuf.begin();

  for (int i = 0; i < nRankSend; ++i, ++MI) {
    int dest = MI->first;
    SU2_MPI::Isend(shortSendBuf[i].data(), shortSendBuf[i].size(), MPI_SHORT, dest, dest, SU2_MPI::GetComm(),
                   &commReqs[3 * i]);
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG, dest, dest + 1, SU2_MPI::GetComm(),
                   &commReqs[3 * i + 1]);
    SU2_MPI::Isend(doubleSendBuf[i].data(), doubleSendBuf[i].size(), MPI_DOUBLE, dest, dest + 2, SU2_MPI::GetComm(),
                   &commReqs[3 * i + 2]);
  }

  /* Loop over the number of ranks from which I receive data. */
  for (int i = 0; i < nRankRecv; ++i) {
    /* Block until a message with shorts arrives from any processor.
       Determine the source and the size of the message.   */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_SHORT, &sizeMess);

    /* Allocate the memory for the short receive buffer and receive the message. */
    shortRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(shortRecvBuf[i].data(), sizeMess, MPI_SHORT, source, rank, SU2_MPI::GetComm(), &status);

    /* Block until the corresponding message with longs arrives, determine
       its size, allocate the memory and receive the message. */
    SU2_MPI::Probe(source, rank + 1, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);
    longRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(longRecvBuf[i].data(), sizeMess, MPI_LONG, source, rank + 1, SU2_MPI::GetComm(), &status);

    /* Idem for the message with doubles. */
    SU2_MPI::Probe(source, rank + 2, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_DOUBLE, &sizeMess);
    doubleRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(doubleRecvBuf[i].data(), sizeMess, MPI_DOUBLE, source, rank + 2, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(3 * nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffers. ---*/
  shortRecvBuf[0] = shortSendBuf[0];
  longRecvBuf[0] = longSendBuf[0];
  doubleRecvBuf[0] = doubleSendBuf[0];

#endif

  /*--- Release the memory of the send buffers. To make sure that all
        the memory is deleted, the swap function is used. ---*/
  for (int i = 0; i < nRankSend; ++i) {
    vector<short>().swap(shortSendBuf[i]);
    vector<long>().swap(longSendBuf[i]);
    vector<su2double>().swap(doubleSendBuf[i]);
  }

  /*--- Allocate the memory for the number of elements for every boundary
        marker and initialize them to zero.     ---*/
  nElem_Bound = new unsigned long[nMarker];
  for (unsigned short i = 0; i < nMarker; ++i) nElem_Bound[i] = 0;

  /*--- Determine the global element ID's of the elements stored on this rank.
        Sort them in increasing order, such that an easy search can be done.
        In the same loop determine the upper bound for the local nodes (without
        halos) and the number of boundary elements for every marker. ---*/
  nElem = nPoint = 0;
  for (int i = 0; i < nRankRecv; ++i) nElem += longRecvBuf[i][0];

  vector<unsigned long> globalElemID;
  globalElemID.reserve(nElem);

  for (int i = 0; i < nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for (long j = 0; j < longRecvBuf[i][0]; ++j) {
      globalElemID.push_back(longRecvBuf[i][indL]);

      unsigned short nDOFsGrid = shortRecvBuf[i][indS + 3];
      unsigned short nFaces = shortRecvBuf[i][indS + 5];
      indS += 3 * nFaces + 8;
      indL += nDOFsGrid + nFaces + 2;
    }

    long nNodesThisRank = longRecvBuf[i][indL];
    nPoint += nNodesThisRank;
    indL += nNodesThisRank + 1;

    for (unsigned iMarker = 0; iMarker < nMarker; ++iMarker) {
      long nBoundElemThisRank = longRecvBuf[i][indL];
      ++indL;
      nElem_Bound[iMarker] += nBoundElemThisRank;

      for (long j = 0; j < nBoundElemThisRank; ++j) {
        short nDOFsBoundElem = shortRecvBuf[i][indS + 2];
        short nDonorsWallFunction = shortRecvBuf[i][indS + 3];
        indS += 4;
        indL += nDOFsBoundElem + nDonorsWallFunction + 2;
      }
    }
  }

  sort(globalElemID.begin(), globalElemID.end());

  /*--- Determine the global element ID's of the halo elements. A vector of
        CUnsignedLong2T is used for this purpose, such that a possible periodic
        transformation can be taken into account. Neighbors with a periodic
        transformation will always become a halo element, even if the element
        is stored on this rank. Furthermore a vector of CReorderElements
        is created for the owned elements to be able to reorder the elements,
        see step 2 below. ---*/
  vector<CUnsignedLong2T> haloElements;
  vector<CReorderElements> ownedElements;
  unsigned short maxTimeLevelLoc = 0;

  for (int i = 0; i < nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for (long j = 0; j < longRecvBuf[i][0]; ++j) {
      unsigned long globalID = longRecvBuf[i][indL];
      unsigned short VTK_Type = shortRecvBuf[i][indS];
      unsigned short nPolySol = shortRecvBuf[i][indS + 2];
      unsigned short nDOFsGrid = shortRecvBuf[i][indS + 3];
      unsigned short nFaces = shortRecvBuf[i][indS + 5];
      unsigned short timeLevel = shortRecvBuf[i][indS + 6];
      bool JacConstant = (bool)shortRecvBuf[i][indS + 7];
      bool commSol = false;

      /* If the integration rule for constant and non-constant Jacobian
         elements is the same, simply set JacConstant to false, such that
         the elements with constant and non-constant Jacobians are
         considered the same. */
      if (JacConstant) {
        const auto orderExactStraight = (unsigned short)ceil(nPolySol * config->GetQuadrature_Factor_Straight());
        const auto orderExactCurved = (unsigned short)ceil(nPolySol * config->GetQuadrature_Factor_Curved());
        if (orderExactStraight == orderExactCurved) JacConstant = false;
      }

      /* Update the local value of the maximum time level. */
      maxTimeLevelLoc = max(maxTimeLevelLoc, timeLevel);

      /*--- Loop over the faces of this element to determine the halo elements
            and the information needed to reorder the owned elements. ---*/
      indS += 8;
      indL += nDOFsGrid + 2;
      for (unsigned short k = 0; k < nFaces; ++k, indS += 3, ++indL) {
        if (longRecvBuf[i][indL] != -1) { /* -1 indicates a boundary face. */

          /* Check if the neighbor of the face is also an owned element.
             Per definition an element for which a periodic transformation
             is carried out, is a halo element, even if the parent element
             is stored locally. */
          bool neighborIsInternal = false;
          if (shortRecvBuf[i][indS] == -1) { /* -1 indicates no periodic transformation. */
            neighborIsInternal = binary_search(globalElemID.begin(), globalElemID.end(), longRecvBuf[i][indL]);
          }

          /* Check if this neighbor is not internal and if the element owns the face. */
          if (!neighborIsInternal) {
            if (shortRecvBuf[i][indS + 2]) {
              /* The face is owned by this element. As the neighboring element
                 is not owned, this implies that a halo element must be created. */
              haloElements.emplace_back(longRecvBuf[i][indL],
                                        shortRecvBuf[i][indS] + 1); /* The +1, because haloElements */
            }                                                       /* are unsigned longs.          */
            else {
              /* The face is not owned by this element and therefore it is owned
                 by the neighboring element on a different rank. Consequently the
                 solution of this element must be communicated. */
              commSol = true;
            }
          }
        }
      }

      /* Store the required data for this element in ownedElements. */
      ownedElements.emplace_back(globalID, timeLevel, commSol, VTK_Type, nPolySol, JacConstant);
    }

    /* Skip the part with the node numbers. */
    long nNodesThisRank = longRecvBuf[i][indL];
    indL += nNodesThisRank + 1;

    /* Loop over the boundary markers. */
    for (unsigned iMarker = 0; iMarker < nMarker; ++iMarker) {
      long nBoundElemThisRank = longRecvBuf[i][indL];
      ++indL;

      /* Loop over the boundary elements coming from this rank. */
      for (long j = 0; j < nBoundElemThisRank; ++j) {
        short nDOFsBoundElem = shortRecvBuf[i][indS + 2];
        short nDonorsWallFunctions = shortRecvBuf[i][indS + 3];
        indS += 4;
        indL += nDOFsBoundElem + 2;

        /* Loop over the donors for the wall functions. */
        for (short k = 0; k < nDonorsWallFunctions; ++k, ++indL) {
          /* Check for an external donor. If external, store it in
             haloElements with no periodic transformation. */
          if (!binary_search(globalElemID.begin(), globalElemID.end(), longRecvBuf[i][indL]))
            haloElements.emplace_back(longRecvBuf[i][indL], 0);
        }
      }
    }
  }

  /* Sort the halo elements in increasing order and remove the double entities. */
  sort(haloElements.begin(), haloElements.end());
  auto lastHalo = unique(haloElements.begin(), haloElements.end());
  haloElements.erase(lastHalo, haloElements.end());

  /* Determine the maximum global time level and possibly reset the number
     of time levels in config. */
  unsigned short maxTimeLevelGlob = maxTimeLevelLoc;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&maxTimeLevelLoc, &maxTimeLevelGlob, 1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

  const unsigned short nTimeLevels = maxTimeLevelGlob + 1;
  config->SetnLevels_TimeAccurateLTS(nTimeLevels);

  /*----------------------------------------------------------------------------*/
  /*--- Step 2: Find out on which rank the halo elements are stored.         ---*/
  /*----------------------------------------------------------------------------*/

  /* Determine the number of owned elements and the total number of
     elements stored on this rank. */
  nVolElemOwned = globalElemID.size();
  nVolElemTot = nVolElemOwned + haloElements.size();

  /* Determine the map from the global element ID to the current storage
     sequence of ownedElements. */
  map<unsigned long, unsigned long> mapGlobalElemIDToInd;
  for (unsigned long i = 0; i < nVolElemOwned; ++i) mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = i;

  /* Determine the number of elements per rank of the originally partitioned grid
     stored in cumulative storage format. */
  vector<unsigned long> nElemPerRankOr(size + 1);

  for (int i = 0; i < size; ++i) nElemPerRankOr[i] = geometry->beg_node[i];
  nElemPerRankOr[size] = geometry->end_node[size - 1];

  /* Determine to which ranks I have to send messages to find out the information
     of the halos stored on this rank. */
  sendToRank.assign(size, 0);

  for (unsigned long i = 0; i < haloElements.size(); ++i) {
    /* Determine the rank where this halo element was originally stored. */
    vector<unsigned long>::iterator low;
    low = lower_bound(nElemPerRankOr.begin(), nElemPerRankOr.end(), haloElements[i].long0);
    unsigned long rankHalo = low - nElemPerRankOr.begin();
    if (*low > haloElements[i].long0) --rankHalo;

    sendToRank[rankHalo] = 1;
  }

  rankToIndCommBuf.clear();
  for (int i = 0; i < size; ++i) {
    if (sendToRank[i]) {
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  /* Resize the first index of the long send buffers for the communication of
     the halo data.        */
  nRankSend = (int)rankToIndCommBuf.size();
  longSendBuf.resize(nRankSend);

  /* Determine the number of ranks, from which this rank will receive elements. */
  nRankRecv = nRankSend;

#ifdef HAVE_MPI
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /* Loop over the local halo elements to fill the communication buffers. */
  for (unsigned long i = 0; i < haloElements.size(); ++i) {
    /* Determine the rank where this halo element was originally stored. */
    vector<unsigned long>::iterator low;
    low = lower_bound(nElemPerRankOr.begin(), nElemPerRankOr.end(), haloElements[i].long0);
    unsigned long ind = low - nElemPerRankOr.begin();
    if (*low > haloElements[i].long0) --ind;

    /* Convert this rank to the index in the send buffer. */
    MI = rankToIndCommBuf.find((int)ind);
    ind = MI->second;

    /* Store the global element ID and the periodic index in the long buffer.
       The subtraction of 1 is there to obtain the correct periodic index.
       In haloElements a +1 is added, because this variable is of unsigned long,
       which cannot handle negative numbers. */
    long perIndex = haloElements[i].long1 - 1;

    longSendBuf[ind].push_back(haloElements[i].long0);
    longSendBuf[ind].push_back(perIndex);
  }

  /* Define a second set of long receive buffers, because the information from
     the first set is still needed later on. Also define the vector to
     store the ranks from which the message came. */
  vector<vector<long> > longSecondRecvBuf(nRankRecv, vector<long>(0));
  vector<int> sourceRank(nRankRecv);

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/

#ifdef HAVE_MPI

  /* Parallel mode. Send all the data using non-blocking sends. */
  commReqs.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for (int i = 0; i < nRankSend; ++i, ++MI) {
    int dest = MI->first;
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG, dest, dest, SU2_MPI::GetComm(),
                   &commReqs[i]);
  }

  /* Loop over the number of ranks from which I receive data. */
  for (int i = 0; i < nRankRecv; ++i) {
    /* Block until a message with longs arrives from any processor.
       Determine the source and the size of the message and receive it. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG, sourceRank[i], rank, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for (int i = 0; i < nRankRecv; ++i) longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /*--- Release the memory of the send buffers. To make sure that all the memory
        is deleted, the swap function is used. Afterwards resize the first index
        of the send buffers to nRankRecv, because this number of messages must
        be sent back to the sending ranks with halo information. ---*/
  for (int i = 0; i < nRankSend; ++i) {
    vector<long>().swap(longSendBuf[i]);
  }

  longSendBuf.resize(nRankRecv);

#ifdef HAVE_MPI
  /* Resize the vector of the communication requests to the number of messages
     to be sent by this rank. Only in parallel node. */
  commReqs.resize(nRankRecv);
#endif

  /*--- Loop over the receive buffers to fill and send the send buffers again. ---*/
  for (int i = 0; i < nRankRecv; ++i) {
    /* Determine the number of elements present in longSecondRecvBuf[i] and
       reserve the memory for the send buffer. */
    const long nElemBuf = longSecondRecvBuf[i].size() / 2;
    longSendBuf[i].reserve(3 * nElemBuf);

    /* Loop over the elements stored in the receive buffer. */
    for (long j = 0; j < nElemBuf; ++j) {
      /* Get the global element ID and periodic index from the receive buffer. */
      const long globalID = longSecondRecvBuf[i][2 * j];
      const long perInd = longSecondRecvBuf[i][2 * j + 1];

      /* Determine the local index of the element in the original partitioning.
         Check if the index is valid. */
      const long localID = globalID - geometry->beg_node[rank];
      if (localID < 0 || localID >= (long)geometry->nPointLinear[rank]) {
        ostringstream message;
        message << localID << " " << geometry->nPointLinear[rank] << endl;
        message << "Invalid local element ID";
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /* Determine which rank owns this element and store everything in the
         send buffer. */
      longSendBuf[i].push_back(globalID);
      longSendBuf[i].push_back(perInd);
      longSendBuf[i].push_back(geometry->elem[localID]->GetColor());
    }

    /* Release the memory of this receive buffer. */
    vector<long>().swap(longSecondRecvBuf[i]);

    /*--- Send the send buffer back to the calling rank.
          Only in parallel mode of course.     ---*/
#ifdef HAVE_MPI
    int dest = sourceRank[i];
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG, dest, dest + 1, SU2_MPI::GetComm(),
                   &commReqs[i]);
#endif
  }

  /*--- Resize the first index of the receive buffers to nRankSend, such that
        the requested halo information can be received.     ---*/
  longSecondRecvBuf.resize(nRankSend);

  /*--- Receive the communication data from the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/
#ifdef HAVE_MPI

  /* Parallel mode. Loop over the number of ranks from which I receive data
     in the return communication, i.e. nRankSend. */
  for (int i = 0; i < nRankSend; ++i) {
    /* Block until a message with longs arrives from any processor.
       Determine the source and the size of the message.   */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank + 1, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    /* Allocate the memory for the long receive buffer and receive the message. */
    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG, source, rank + 1, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends and synchronize the ranks, because
     wild cards have been used. */
  SU2_MPI::Waitall(nRankRecv, commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for (int i = 0; i < nRankRecv; ++i) longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /* Release the memory of the send buffers. To make sure that all
     the memory is deleted, the swap function is used. */
  for (int i = 0; i < nRankRecv; ++i) vector<long>().swap(longSendBuf[i]);

  /* Copy the data from the receive buffers into a class of CLong3T, such that
     it can be sorted in increasing order. Note that the rank of the element
     is stored first, followed by its global ID and last the periodic index. */
  vector<CLong3T> haloData;
  for (int i = 0; i < nRankSend; ++i) {
    const long nElemBuf = longSecondRecvBuf[i].size() / 3;

    for (long j = 0; j < nElemBuf; ++j) {
      const long j3 = 3 * j;
      haloData.emplace_back(longSecondRecvBuf[i][j3 + 2], longSecondRecvBuf[i][j3], longSecondRecvBuf[i][j3 + 1]);
    }

    /* Release the memory of this receive buffer. */
    vector<long>().swap(longSecondRecvBuf[i]);
  }

  /* Sort halo data in increasing order. */
  sort(haloData.begin(), haloData.end());

  /* Determine the number of halo elements per rank in cumulative storage.
     The first element of this vector is nVolElemOwned, such that this vector
     contains the starting position in the vector volElem. */
  vector<unsigned long> nHaloElemPerRank(size + 1, 0);
  for (unsigned long i = 0; i < haloData.size(); ++i) ++nHaloElemPerRank[haloData[i].long0 + 1];

  nHaloElemPerRank[0] = nVolElemOwned;
  for (int i = 0; i < size; ++i) nHaloElemPerRank[i + 1] += nHaloElemPerRank[i];

  if (nHaloElemPerRank[size] != nVolElemTot)
    SU2_MPI::Error("Inconsistency in total number of volume elements", CURRENT_FUNCTION);

  /* Determine the number of ranks to which I have to send data in this cycle. */
  sendToRank.assign(size, 0);
  rankToIndCommBuf.clear();
  for (int i = 0; i < size; ++i) {
    if (nHaloElemPerRank[i + 1] > nHaloElemPerRank[i]) {
      sendToRank[i] = 1;
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  nRankSend = (int)rankToIndCommBuf.size();

  /* Store the value of nRankSend for later use. */
  const int nRankSendHaloInfo = nRankSend;

  /* Determine the number of ranks, from which this rank will receive elements. */
  nRankRecv = nRankSend;

#ifdef HAVE_MPI
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /* Copy the data to be sent to the send buffers. */
  longSendBuf.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for (int i = 0; i < nRankSend; ++i, ++MI) {
    int dest = MI->first;
    for (unsigned long j = nHaloElemPerRank[dest]; j < nHaloElemPerRank[dest + 1]; ++j) {
      const unsigned long jj = j - nVolElemOwned;
      longSendBuf[i].push_back(haloData[jj].long1);
      longSendBuf[i].push_back(haloData[jj].long2);
    }
  }

  /* Resize the first index of the long receive buffer. */
  longSecondRecvBuf.resize(nRankRecv);

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/
#ifdef HAVE_MPI

  /* Parallel mode. Send all the data using non-blocking sends. */
  commReqs.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for (int i = 0; i < nRankSend; ++i, ++MI) {
    int dest = MI->first;
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG, dest, dest, SU2_MPI::GetComm(),
                   &commReqs[i]);
  }

  /* Resize the vector to store the ranks from which the message came. */
  sourceRank.resize(nRankRecv);

  /* Loop over the number of ranks from which I receive data. */
  for (int i = 0; i < nRankRecv; ++i) {
    /* Block until a message with longs arrives from any processor.
       Determine the source and the size of the message and receive it. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG, sourceRank[i], rank, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends and synchronize the ranks,
     because wild cards have been used. */
  SU2_MPI::Waitall(nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for (int i = 0; i < nRankSend; ++i) longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /* Release the memory of the send buffers. To make sure that all the memory
     is deleted, the swap function is used. */
  for (int i = 0; i < nRankSend; ++i) {
    vector<long>().swap(longSendBuf[i]);
  }

  /*--- Loop over the receive buffers to flag the locally owned elements for
        communication. Although the face information has already been used to
        do this when ownedElements are constructed, elements that are donors
        for the wall function treatment and are not direct neighbors may
        have been missed. ---*/
  for (int i = 0; i < nRankRecv; ++i) {
    const unsigned long nElemBuf = longSecondRecvBuf[i].size() / 2;
    for (unsigned long j = 0; j < nElemBuf; ++j) {
      const unsigned long elemID = longSecondRecvBuf[i][2 * j];

      auto MMI = mapGlobalElemIDToInd.find(elemID);
      if (MMI == mapGlobalElemIDToInd.end())
        SU2_MPI::Error("Entry not found in mapGlobalElemIDToInd", CURRENT_FUNCTION);

      ownedElements[MMI->second].SetCommSolution(true);
    }
  }

  /*----------------------------------------------------------------------------*/
  /*--- Step 3: Determine the numbering of the owned elements. The following ---*/
  /*---         criteria are used for the owned elements.                    ---*/
  /*---         - Time level of the element: elements with the smallest time ---*/
  /*---           level are number first, etc.                               ---*/
  /*---         - For each time level the elements that do not need to send  ---*/
  /*---           their data are numbered first, followed by the elements    ---*/
  /*---           that must send their data to other ranks. Note that not    ---*/
  /*---           sending the solution does not mean that the residual can   ---*/
  /*---           be built without communication. It is possible that a face ---*/
  /*---           is owned by a local element, but it is adjacent to an      ---*/
  /*---           element owned by a different rank. In that case the data   ---*/
  /*---           from the neighboring element is communicated and stored in ---*/
  /*---           a halo element. However, the residual of these internal    ---*/
  /*---           elements do not receive a contribution computed on a       ---*/
  /*---           different rank.                                            ---*/
  /*---         - Some of the elements that do not need to send their data   ---*/
  /*---           to other ranks may still be put in that part to increase   ---*/
  /*---           the efficiency when treating element simultaneously, see   ---*/
  /*---           next point.                                                ---*/
  /*---         - Elements of the same type are put contiguously in memory   ---*/
  /*---           as much as possible to facilitate the simultaneous         ---*/
  /*---           treatment in the matrix multiplications to increase        ---*/
  /*---           performance.                                               ---*/
  /*---         - A reverse Cuthill McKee renumbering takes place to obtain  ---*/
  /*---           better cache performance for the face residuals.           ---*/
  /*----------------------------------------------------------------------------*/

  /* Determine the number of elements/faces that are treated simultaneously */
  /* in the matrix products to obtain good gemm performance.                */
  const unsigned short nElemSimul = config->GetSizeMatMulPadding() / nVar;

  /* Determine the number of different element types present. */
  map<unsigned short, unsigned short> mapElemTypeToInd;
  for (auto OEI = ownedElements.begin(); OEI != ownedElements.end(); ++OEI) {
    const unsigned short elType = OEI->GetElemType();

    if (mapElemTypeToInd.find(elType) == mapElemTypeToInd.end()) {
      const unsigned short ind = mapElemTypeToInd.size();
      mapElemTypeToInd[elType] = ind;
    }
  }

  /* Determine the number of elements per element type per time level.
     Make a distinction between internal elements (elements that are not
     communicated) and elements that are communicated. Later on these vectors
     will be put in cumulative storage format, which explains the +1 for
     the second index. */
  vector<vector<unsigned long> > nInternalElem(nTimeLevels, vector<unsigned long>(mapElemTypeToInd.size() + 1));
  vector<vector<unsigned long> > nCommElem(nTimeLevels, vector<unsigned long>(mapElemTypeToInd.size() + 1));
  for (unsigned short i = 0; i < nTimeLevels; ++i) {
    for (unsigned long j = 0; j <= mapElemTypeToInd.size(); ++j) {
      nInternalElem[i][j] = 0;
      nCommElem[i][j] = 0;
    }
  }

  for (auto OEI = ownedElements.begin(); OEI != ownedElements.end(); ++OEI) {
    const unsigned short elType = OEI->GetElemType();
    map<unsigned short, unsigned short>::const_iterator MI = mapElemTypeToInd.find(elType);
    unsigned short ind = MI->second + 1;

    if (OEI->GetCommSolution())
      ++nCommElem[OEI->GetTimeLevel()][ind];
    else
      ++nInternalElem[OEI->GetTimeLevel()][ind];
  }

  /* Loop again over the owned elements and check if elements, which do not
     have to send their solution, should be flagged as such in order to improve
     the gemm performance. */
  for (auto OEI = ownedElements.begin(); OEI != ownedElements.end(); ++OEI) {
    /* Check for an internal element, i.e. an element for which the solution
       does not need to be communicated. */
    if (!OEI->GetCommSolution()) {
      /* Determine the time level and the index of the element type. */
      const unsigned short tLev = OEI->GetTimeLevel();

      const unsigned short elType = OEI->GetElemType();
      map<unsigned short, unsigned short>::const_iterator MI = mapElemTypeToInd.find(elType);
      const unsigned short ind = MI->second + 1;

      /* Determine whether or not this element must be moved from internal to
         comm elements to improve performance. Change the appropriate data
         when this must happen. */
      if (nInternalElem[tLev][ind] % nElemSimul && nCommElem[tLev][ind] % nElemSimul) {
        OEI->SetCommSolution(true);
        --nInternalElem[tLev][ind];
        ++nCommElem[tLev][ind];
      }
    }
  }

  /* Check whether there are enough elements in every partition for every
     time level to guarantee a good gemm performance. */
  unsigned long nFullChunks = 0, nPartialChunks = 0;
  for (unsigned short tLev = 0; tLev < nTimeLevels; ++tLev) {
    for (unsigned long j = 1; j < nInternalElem[tLev].size(); ++j) {
      nFullChunks += nInternalElem[tLev][j] / nElemSimul;
      if (nInternalElem[tLev][j] % nElemSimul) ++nPartialChunks;
    }

    for (unsigned long j = 1; j < nCommElem[tLev].size(); ++j) {
      nFullChunks += nCommElem[tLev][j] / nElemSimul;
      if (nCommElem[tLev][j] % nElemSimul) ++nPartialChunks;
    }
  }

  unsigned long tooManyPartChunksLoc = 0;
  if (5 * nPartialChunks >= (nPartialChunks + nFullChunks)) tooManyPartChunksLoc = 1;

  /* Determine the number of ranks which contain too many partial chunks.
     The result only needs to be known on the master node. */
  unsigned long nRanksTooManyPartChunks = tooManyPartChunksLoc;
#ifdef HAVE_MPI
  SU2_MPI::Reduce(&tooManyPartChunksLoc, &nRanksTooManyPartChunks, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE,
                  SU2_MPI::GetComm());
#endif

  if ((rank == MASTER_NODE) && (nRanksTooManyPartChunks != 0) && (size > 1)) {
    cout << endl << "                 WARNING" << endl;
    cout << "There are " << nRanksTooManyPartChunks << " partitions for which "
         << " the simultaneous treatment of volume elements is not optimal." << endl;
    cout << "This could lead to a significant load imbalance." << endl;
  }

  /* Put nInternalElem and nCommElem in cumulative storage format. */
  for (unsigned short tLev = 0; tLev < nTimeLevels; ++tLev) {
    if (tLev) nInternalElem[tLev][0] = nCommElem[tLev - 1].back();
    for (unsigned long j = 1; j < nInternalElem[tLev].size(); ++j) nInternalElem[tLev][j] += nInternalElem[tLev][j - 1];

    nCommElem[tLev][0] = nInternalElem[tLev].back();
    for (unsigned long j = 1; j < nCommElem[tLev].size(); ++j) nCommElem[tLev][j] += nCommElem[tLev][j - 1];
  }

  /* In the solver only the number of elements per time level is needed.
     The number of owned elements, nVolElemOwnedPerTimeLevel, is stored in
     cumulative storage format and can be retrieved easily from nCommElem.
     The number of internal elements per time level is stored normally and
     can be retrieved from nInternalElem. */
  nVolElemOwnedPerTimeLevel.resize(nTimeLevels + 1);
  nVolElemInternalPerTimeLevel.resize(nTimeLevels);

  nVolElemOwnedPerTimeLevel[0] = 0;
  for (unsigned short tLev = 0; tLev < nTimeLevels; ++tLev) {
    nVolElemOwnedPerTimeLevel[tLev + 1] = nCommElem[tLev].back();
    nVolElemInternalPerTimeLevel[tLev] = nInternalElem[tLev].back() - nInternalElem[tLev][0];
  }

  /* Sort the elements of ownedElements in increasing order. */
  sort(ownedElements.begin(), ownedElements.end());

  /* At the moment the elements are sorted per time level and element type.  */
  /* This is in principle enough to make to code work with the simultaneous  */
  /* treatment for volume elements. Below follows a renumbering within the   */
  /* time level and element type (Reverse Cuthill-McKee), such that a better */
  /* cache performance is obtained for the surface elements. For high order  */
  /* elements this effect is probably marginal.                              */

  /* Determine the map from the global element ID to the current storage
     sequence of ownedElements. */
  mapGlobalElemIDToInd.clear();
  for (unsigned long i = 0; i < nVolElemOwned; ++i) mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = i;

  /*--- Create the graph of local elements. The halo elements are ignored. ---*/
  vector<vector<unsigned long> > neighElem(nVolElemOwned, vector<unsigned long>(0));

  nRankRecv = (int)longRecvBuf.size();
  for (int i = 0; i < nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for (long j = 0; j < longRecvBuf[i][0]; ++j) {
      unsigned long globalID = longRecvBuf[i][indL];
      unsigned short nDOFsGrid = shortRecvBuf[i][indS + 3];
      unsigned short nFaces = shortRecvBuf[i][indS + 5];

      auto MMI = mapGlobalElemIDToInd.find(globalID);
      unsigned long ind = MMI->second;

      indS += 8;
      indL += nDOFsGrid + 2;
      for (unsigned short k = 0; k < nFaces; ++k, indS += 3, ++indL) {
        if ((longRecvBuf[i][indL] != -1) && (shortRecvBuf[i][indS] == -1)) {  // Check for internal owned node.

          MMI = mapGlobalElemIDToInd.find(longRecvBuf[i][indL]);
          if (MMI != mapGlobalElemIDToInd.end()) neighElem[ind].push_back(MMI->second);
        }
      }
    }
  }

  /* Sort the neighbors of each element in increasing order. */
  for (unsigned long i = 0; i < nVolElemOwned; ++i) sort(neighElem[i].begin(), neighElem[i].end());

  /* Define the vector, which contains the new numbering of the owned elements
     w.r.t. to the numbering currently stored in ownedElements. Note that signed
     longs are used for this purpose, because the vector is initialized with -1
     to indicate that no new number has been assigned yet. */
  vector<long> oldElemToNewElem(nVolElemOwned, -1);

  /*--- While loop to carry out the renumbering. A while loop is present,
        because the local partition may not be contiguous. ---*/
  unsigned long nElemRenumbered = 0;
  while (nElemRenumbered < nVolElemOwned) {
    /* Determine the first element in the list that has not been renumbered. */
    unsigned long indBeg;
    for (indBeg = 0; indBeg < nVolElemOwned; ++indBeg)
      if (oldElemToNewElem[indBeg] == -1) break;

    /* Determine the time level and element type of the element indBeg and
       determine the element range in which the starting element for the
       renumbering must be sought. */
    unsigned short timeLevel = ownedElements[indBeg].GetTimeLevel();
    unsigned short elType = ownedElements[indBeg].GetElemType();

    map<unsigned short, unsigned short>::const_iterator MI = mapElemTypeToInd.find(elType);
    unsigned short ind = MI->second;

    unsigned long indEnd;
    if (ownedElements[indBeg].GetCommSolution())
      indEnd = nCommElem[timeLevel][ind + 1];
    else
      indEnd = nInternalElem[timeLevel][ind + 1];

    /* Determine the element in the range [indBeg,indEnd) with the least number
       of neighbors that has not been renumbered yet. This is the starting
       element for the current renumbering round. */
    for (unsigned long i = (indBeg + 1); i < indEnd; ++i) {
      if ((oldElemToNewElem[i] == -1) && (neighElem[i].size() < neighElem[indBeg].size())) indBeg = i;
    }

    /* Start of the Reverse Cuthil McKee renumbering. */
    vector<unsigned long> frontElements(1, indBeg);
    while (!frontElements.empty()) {
      /* Vector, which stores the front for the next round. */
      vector<unsigned long> frontElementsNew;

      /* Loop over the elements of the current front. */
      for (unsigned long i = 0; i < frontElements.size(); ++i) {
        /* Carry out the renumbering for this element. */
        const unsigned long iFront = frontElements[i];

        timeLevel = ownedElements[iFront].GetTimeLevel();
        elType = ownedElements[iFront].GetElemType();
        MI = mapElemTypeToInd.find(elType);
        ind = MI->second;

        if (ownedElements[iFront].GetCommSolution())
          oldElemToNewElem[iFront] = nCommElem[timeLevel][ind]++;
        else
          oldElemToNewElem[iFront] = nInternalElem[timeLevel][ind]++;

        /* Store the neighbors that have not been renumbered yet in the front
           for the next round. Set its index to -2 to indicate that the element
           is already on the new front. */
        for (unsigned long j = 0; j < neighElem[iFront].size(); ++j) {
          if (oldElemToNewElem[neighElem[iFront][j]] == -1) {
            frontElementsNew.push_back(neighElem[iFront][j]);
            oldElemToNewElem[neighElem[iFront][j]] = -2;
          }
        }
      }

      /* Update the counter nElemRenumbered. */
      nElemRenumbered += frontElements.size();

      /* Sort frontElementsNew in increasing order. */
      sort(frontElementsNew.begin(), frontElementsNew.end());

      /* Store the new front elements in frontElements for the next round. */
      frontElements = frontElementsNew;
    }
  }

  if (nElemRenumbered != nVolElemOwned) SU2_MPI::Error("Something went wrong in the renumbering", CURRENT_FUNCTION);

  /* Determine the final mapping from the global element number to the local
     entry for the owned elements. First clear mapGlobalElemIDToInd before
     it can be used to store its correct content. */
  mapGlobalElemIDToInd.clear();
  for (unsigned long i = 0; i < nVolElemOwned; ++i)
    mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = oldElemToNewElem[i];

    /*----------------------------------------------------------------------------*/
    /*--- Step 4: Store the elements, nodes and boundary elements in the data  ---*/
    /*---         structures used by the FEM solver.                           ---*/
    /*----------------------------------------------------------------------------*/

    /*--- Check in parallel mode for empty partitions. If present, print a warning.
          The solver is capable of handling empty partitions, but it may not be
          efficient. ---*/
#ifdef HAVE_MPI
  unsigned long thisPartitionEmpty = nVolElemOwned ? 0 : 1;
  unsigned long nEmptyPartitions = 0;

  SU2_MPI::Reduce(&thisPartitionEmpty, &nEmptyPartitions, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE,
                  SU2_MPI::GetComm());

  if (rank == MASTER_NODE && nEmptyPartitions) {
    cout << endl << "         WARNING" << endl;
    cout << "There are " << nEmptyPartitions << " empty partitions present." << endl;
    cout << "SU2 is able to handle this, but it may be inefficient." << endl << endl;
  }
#endif

  /*--- Allocate the memory for the volume elements, the nodes
        and the surface elements of the boundaries.    ---*/
  volElem.resize(nVolElemTot);
  meshPoints.reserve(nPoint);

  boundaries.resize(nMarker);
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    boundaries[iMarker].markerTag = config->GetMarker_All_TagBound(iMarker);
    boundaries[iMarker].periodicBoundary = config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY;
    boundaries[iMarker].surfElem.reserve(nElem_Bound[iMarker]);
  }

  /*--- Copy the data from the communication buffers. ---*/
  for (int i = 0; i < nRankRecv; ++i) {
    /* The data for the volume elements. Loop over these elements in the buffer. */
    unsigned long indL = 1, indS = 0, indD = 0;
    for (long j = 0; j < longRecvBuf[i][0]; ++j) {
      /* Determine the location in volElem where this data must be stored. */
      unsigned long elemID = longRecvBuf[i][indL++];
      auto MMI = mapGlobalElemIDToInd.find(elemID);
      unsigned long ind = MMI->second;

      /* Store the data. */
      volElem[ind].elemIsOwned = true;
      volElem[ind].rankOriginal = rank;
      volElem[ind].periodIndexToDonor = -1;

      volElem[ind].VTK_Type = shortRecvBuf[i][indS++];
      volElem[ind].nPolyGrid = shortRecvBuf[i][indS++];
      volElem[ind].nPolySol = shortRecvBuf[i][indS++];
      volElem[ind].nDOFsGrid = shortRecvBuf[i][indS++];
      volElem[ind].nDOFsSol = shortRecvBuf[i][indS++];
      volElem[ind].nFaces = shortRecvBuf[i][indS++];
      volElem[ind].timeLevel = shortRecvBuf[i][indS++];

      volElem[ind].JacIsConsideredConstant = (bool)shortRecvBuf[i][indS++];

      volElem[ind].elemIDGlobal = elemID;
      volElem[ind].offsetDOFsSolGlobal = longRecvBuf[i][indL++];

      volElem[ind].nodeIDsGrid.resize(volElem[ind].nDOFsGrid);
      volElem[ind].JacFacesIsConsideredConstant.resize(volElem[ind].nFaces);
      volElem[ind].ElementOwnsFaces.resize(volElem[ind].nFaces);

      for (unsigned short k = 0; k < volElem[ind].nDOFsGrid; ++k) volElem[ind].nodeIDsGrid[k] = longRecvBuf[i][indL++];

      for (unsigned short k = 0; k < volElem[ind].nFaces; ++k) {
        long neighBorID = longRecvBuf[i][indL++];

        ++indS;  // At this location the periodic index of the face is stored in
                 // shortRecvBuf, which is not stored in volElem.
        volElem[ind].JacFacesIsConsideredConstant[k] = (bool)shortRecvBuf[i][indS++];
        volElem[ind].ElementOwnsFaces[k] = (bool)shortRecvBuf[i][indS++];

        if (neighBorID == -1) volElem[ind].ElementOwnsFaces[k] = true;  // Boundary faces are always owned.
      }

      volElem[ind].lenScale = doubleRecvBuf[i][indD++];
    }

    /* The data for the nodes. Loop over these nodes in the buffer and store
       them in meshPoints. */
    unsigned long nNodesThisRank = longRecvBuf[i][indL++];
    for (unsigned long j = 0; j < nNodesThisRank; ++j) {
      CPointFEM thisPoint;
      thisPoint.globalID = longRecvBuf[i][indL++];
      thisPoint.periodIndexToDonor = -1;
      for (unsigned short k = 0; k < nDim; ++k) thisPoint.coor[k] = doubleRecvBuf[i][indD++];

      meshPoints.push_back(thisPoint);
    }

    /* The data for the boundary markers. Loop over them. */
    for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
      unsigned long nElemThisRank = longRecvBuf[i][indL++];
      for (unsigned long j = 0; j < nElemThisRank; ++j) {
        CSurfaceElementFEM thisSurfElem;

        thisSurfElem.VTK_Type = shortRecvBuf[i][indS++];
        thisSurfElem.nPolyGrid = shortRecvBuf[i][indS++];
        thisSurfElem.nDOFsGrid = shortRecvBuf[i][indS++];
        const short nDonors = shortRecvBuf[i][indS++];

        thisSurfElem.volElemID = longRecvBuf[i][indL++];
        thisSurfElem.boundElemIDGlobal = longRecvBuf[i][indL++];

        thisSurfElem.nodeIDsGrid.resize(thisSurfElem.nDOFsGrid);
        for (unsigned short k = 0; k < thisSurfElem.nDOFsGrid; ++k)
          thisSurfElem.nodeIDsGrid[k] = longRecvBuf[i][indL++];

        indL += nDonors;

        /* Convert the global volume element ID to the local one.
           It is essential to do this before the sorting. */
        map<unsigned long, unsigned long>::iterator MMI;
        MMI = mapGlobalElemIDToInd.find(thisSurfElem.volElemID);
        thisSurfElem.volElemID = MMI->second;

        /* Store the surface element in the data structure for this boundary. */
        boundaries[iMarker].surfElem.push_back(thisSurfElem);
      }
    }
  }

  /* Sort meshPoints in increasing order and remove the double entities. */
  sort(meshPoints.begin(), meshPoints.end());
  auto lastPoint = unique(meshPoints.begin(), meshPoints.end());
  meshPoints.erase(lastPoint, meshPoints.end());

  /* Clear the contents of the map globalPointIDToLocalInd and fill
     it with the information present in meshPoints. */
  globalPointIDToLocalInd.clear();
  for (unsigned long i = 0; i < meshPoints.size(); ++i) globalPointIDToLocalInd[meshPoints[i].globalID] = i;

  /*--- All the data from the receive buffers has been copied in the local
        data structures. Release the memory of the receive buffers. To make
        sure that all the memory is deleted, the swap function is used. ---*/
  for (int i = 0; i < nRankRecv; ++i) {
    vector<short>().swap(shortRecvBuf[i]);
    vector<long>().swap(longRecvBuf[i]);
    vector<su2double>().swap(doubleRecvBuf[i]);
  }

  /*--- Sort the surface elements of the boundaries in increasing order. ---*/
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker)
    sort(boundaries[iMarker].surfElem.begin(), boundaries[iMarker].surfElem.end());

  /*----------------------------------------------------------------------------*/
  /*--- Step 5: Obtain the information of the halo elements, which are       ---*/
  /*---         sorted per time level and afterwards per rank, where the     ---*/
  /*---         sequence is determined by the numbering on the sending rank. ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Resize the first index of the send buffers to nRankRecv, because this
        number of messages must be sent back to the sending ranks with halo
        information. ---*/
  nRankRecv = (int)longSecondRecvBuf.size();
  shortSendBuf.resize(nRankRecv);
  longSendBuf.resize(nRankRecv);
  doubleSendBuf.resize(nRankRecv);

#ifdef HAVE_MPI
  /* Resize the vector of the communication requests to the number of messages
     to be sent by this rank. Only in parallel node. */
  commReqs.resize(3 * nRankRecv);
#endif

  /*--- Loop over the receive buffers to fill the send buffers again. ---*/
  for (int i = 0; i < nRankRecv; ++i) {
    /* Loop over the elements in this receive buffer to determine the local
       index on this rank. Note that also the periodic index must be stored,
       hence use an CUnsignedLong2T for this purpose. As -1 cannot be stored
       for an unsigned long a 1 is added to the periodic transformation. */
    const unsigned long nElemBuf = longSecondRecvBuf[i].size() / 2;
    vector<CUnsignedLong2T> elemBuf(nElemBuf);

    for (unsigned long j = 0; j < nElemBuf; ++j) {
      const unsigned long j2 = 2 * j;

      const unsigned long elemID = longSecondRecvBuf[i][j2];
      auto MMI = mapGlobalElemIDToInd.find(elemID);
      if (MMI == mapGlobalElemIDToInd.end())
        SU2_MPI::Error("Entry not found in mapGlobalElemIDToInd", CURRENT_FUNCTION);

      elemBuf[j].long0 = MMI->second;
      elemBuf[j].long1 = longSecondRecvBuf[i][j2 + 1] + 1;
    }

    /* Release the memory of the long receive buffer via the swap function
       and sort elemBuf in increasing order. */
    vector<long>().swap(longSecondRecvBuf[i]);

    sort(elemBuf.begin(), elemBuf.end());

    /* Store the number of elements in the first element of the long send buffer. */
    longSendBuf[i].push_back(nElemBuf);

    /* Vector with node IDs that must be returned to this calling rank.
       Note that also the periodic index must be stored, hence use an
       CUnsignedLong2T for this purpose. */
    vector<CUnsignedLong2T> nodeIDs;

    /* Loop over the elements to fill the send buffers. */
    for (unsigned long j = 0; j < nElemBuf; ++j) {
      /* Store the global element ID in the long buffer,
         the periodic index in the short send buffer and
         the length scale in the double send buffer. */
      const unsigned long indV = elemBuf[j].long0;
      longSendBuf[i].push_back(volElem[indV].elemIDGlobal);

      const short perIndex = (short)elemBuf[j].long1 - 1;  // Note the -1.
      shortSendBuf[i].push_back(perIndex);

      doubleSendBuf[i].push_back(volElem[indV].lenScale);

      /* Store the other relevant information of this element in the short
         and long communication buffers. Also store the node IDs and the
         periodic transformation in nodeIDs. */
      shortSendBuf[i].push_back(volElem[indV].VTK_Type);
      shortSendBuf[i].push_back(volElem[indV].nPolyGrid);
      shortSendBuf[i].push_back(volElem[indV].nPolySol);
      shortSendBuf[i].push_back(volElem[indV].nDOFsGrid);
      shortSendBuf[i].push_back(volElem[indV].nDOFsSol);
      shortSendBuf[i].push_back(volElem[indV].nFaces);
      shortSendBuf[i].push_back(volElem[indV].timeLevel);

      for (unsigned short k = 0; k < volElem[indV].nDOFsGrid; ++k) {
        longSendBuf[i].push_back(volElem[indV].nodeIDsGrid[k]);
        nodeIDs.emplace_back(volElem[indV].nodeIDsGrid[k], elemBuf[j].long1);
      }

      for (unsigned short k = 0; k < volElem[indV].nFaces; ++k)
        shortSendBuf[i].push_back((short)volElem[indV].JacFacesIsConsideredConstant[k]);
    }

    /* Sort nodeIDs in increasing order and remove the double entities. */
    sort(nodeIDs.begin(), nodeIDs.end());
    auto lastNodeID = unique(nodeIDs.begin(), nodeIDs.end());
    nodeIDs.erase(lastNodeID, nodeIDs.end());

    /* Add the number of node IDs and the node IDs itself to longSendBuf[i]
       and the periodix index to shortSendBuf. Note again the -1 for the
       periodic index, because an unsigned long cannot represent -1, the
       value for the periodic index when no peridicity is present. */
    longSendBuf[i].push_back(nodeIDs.size());
    for (unsigned long j = 0; j < nodeIDs.size(); ++j) {
      longSendBuf[i].push_back(nodeIDs[j].long0);
      shortSendBuf[i].push_back((short)nodeIDs[j].long1 - 1);
    }

    /*--- Copy the coordinates to doubleSendBuf. ---*/
    for (unsigned long j = 0; j < nodeIDs.size(); ++j) {
      map<unsigned long, unsigned long>::const_iterator LMI;
      LMI = globalPointIDToLocalInd.find(nodeIDs[j].long0);

      if (LMI == globalPointIDToLocalInd.end()) SU2_MPI::Error("Entry not found in map", CURRENT_FUNCTION);

      unsigned long ind = LMI->second;
      for (unsigned short l = 0; l < nDim; ++l) doubleSendBuf[i].push_back(meshPoints[ind].coor[l]);
    }

    /*--- Send the communication buffers back to the calling rank.
          Only in parallel mode of course.     ---*/
#ifdef HAVE_MPI
    int dest = sourceRank[i];
    SU2_MPI::Isend(shortSendBuf[i].data(), shortSendBuf[i].size(), MPI_SHORT, dest, dest + 1, SU2_MPI::GetComm(),
                   &commReqs[3 * i]);
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG, dest, dest + 2, SU2_MPI::GetComm(),
                   &commReqs[3 * i + 1]);
    SU2_MPI::Isend(doubleSendBuf[i].data(), doubleSendBuf[i].size(), MPI_DOUBLE, dest, dest + 3, SU2_MPI::GetComm(),
                   &commReqs[3 * i + 2]);
#endif
  }

  /*--- Resize the first index of the receive buffers to nRankSendHaloInfo,
        such that the requested halo information can be received.     ---*/
  nRankSend = nRankSendHaloInfo;
  shortRecvBuf.resize(nRankSend);
  longRecvBuf.resize(nRankSend);
  doubleRecvBuf.resize(nRankSend);

  /* Resize the vector to store the ranks from which the message came. */
  sourceRank.resize(nRankSend);

  /*--- Receive the communication data from the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/
#ifdef HAVE_MPI

  /* Parallel mode. Loop over the number of ranks from which I receive data
     in the return communication, i.e. nRankSend. */
  for (int i = 0; i < nRankSend; ++i) {
    /* Block until a message with shorts arrives from any processor.
       Determine the source and the size of the message.   */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank + 1, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_SHORT, &sizeMess);

    /* Allocate the memory for the short receive buffer and receive the message. */
    shortRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(shortRecvBuf[i].data(), sizeMess, MPI_SHORT, sourceRank[i], rank + 1, SU2_MPI::GetComm(), &status);

    /* Block until the corresponding message with longs arrives, determine
       its size, allocate the memory and receive the message. */
    SU2_MPI::Probe(sourceRank[i], rank + 2, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);
    longRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(longRecvBuf[i].data(), sizeMess, MPI_LONG, sourceRank[i], rank + 2, SU2_MPI::GetComm(), &status);

    /* Idem for the message with doubles. */
    SU2_MPI::Probe(sourceRank[i], rank + 3, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_DOUBLE, &sizeMess);
    doubleRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(doubleRecvBuf[i].data(), sizeMess, MPI_DOUBLE, sourceRank[i], rank + 3, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(3 * nRankRecv, commReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /* Sequential mode. Simply copy the buffers. Note that nRankSend is at most 1. */
  for (int i = 0; i < nRankSend; ++i) {
    sourceRank[i] = i;
    shortRecvBuf[i] = shortSendBuf[i];
    longRecvBuf[i] = longSendBuf[i];
    doubleRecvBuf[i] = doubleSendBuf[i];
  }

#endif

  /*--- Release the memory of the send buffers. To make sure that all
        the memory is deleted, the swap function is used. ---*/
  for (int i = 0; i < nRankRecv; ++i) {
    vector<short>().swap(shortSendBuf[i]);
    vector<long>().swap(longSendBuf[i]);
    vector<su2double>().swap(doubleSendBuf[i]);
  }

  /*----------------------------------------------------------------------------*/
  /*--- Step 6: Build the layer of halo elements from the information in the ---*/
  /*---         receive buffers shortRecvBuf, longRecvBuf and doubleRecvBuf. ---*/
  /*----------------------------------------------------------------------------*/

  /*--- The halo elements must be sorted first based on the time level, followed
        by the rank and finally the index in the receive buffers (which represents
        the sequence on the sending rank. This sorting can be accomplished by
        using a vector of CLong3T. The contents of this vector is build below. ---*/
  vector<CLong3T> haloElemInfo;
  haloElemInfo.reserve(nVolElemTot - nVolElemOwned);

  for (int i = 0; i < nRankSend; ++i) {
    unsigned long indS = 0;

    for (long j = 0; j < longRecvBuf[i][0]; ++j) {
      const unsigned short nFaces = shortRecvBuf[i][indS + 6];
      haloElemInfo.emplace_back(shortRecvBuf[i][indS + 7], sourceRank[i], j);
      indS += nFaces + 8;
    }
  }

  sort(haloElemInfo.begin(), haloElemInfo.end());

  /*--- Loop over the receive buffers to store the information of the
        halo elements and the halo points. ---*/
  vector<CPointFEM> haloPoints;
  for (int i = 0; i < nRankSend; ++i) {
    /* Initialization of the indices in the communication buffers. */
    unsigned long indL = 1, indS = 0, indD = 0;

    /* Loop over the halo elements received from this rank. */
    for (long j = 0; j < longRecvBuf[i][0]; ++j) {
      /* Create an object of CLong3T and search for its position
         in haloElemInfo to determine the index in volElem, where
         this element is stored. */
      CLong3T thisElem(shortRecvBuf[i][indS + 7], sourceRank[i], j);
      vector<CLong3T>::iterator low;
      low = lower_bound(haloElemInfo.begin(), haloElemInfo.end(), thisElem);

      unsigned long indV = low - haloElemInfo.begin();
      indV += nVolElemOwned;

      /* Retrieve the data from the communication buffers. */
      volElem[indV].elemIDGlobal = longRecvBuf[i][indL++];
      volElem[indV].rankOriginal = sourceRank[i];

      volElem[indV].periodIndexToDonor = shortRecvBuf[i][indS++];
      volElem[indV].VTK_Type = shortRecvBuf[i][indS++];
      volElem[indV].nPolyGrid = shortRecvBuf[i][indS++];
      volElem[indV].nPolySol = shortRecvBuf[i][indS++];
      volElem[indV].nDOFsGrid = shortRecvBuf[i][indS++];
      volElem[indV].nDOFsSol = shortRecvBuf[i][indS++];
      volElem[indV].nFaces = shortRecvBuf[i][indS++];
      volElem[indV].timeLevel = shortRecvBuf[i][indS++];

      volElem[indV].nodeIDsGrid.resize(volElem[indV].nDOFsGrid);
      for (unsigned short k = 0; k < volElem[indV].nDOFsGrid; ++k)
        volElem[indV].nodeIDsGrid[k] = longRecvBuf[i][indL++];

      volElem[indV].JacFacesIsConsideredConstant.resize(volElem[indV].nFaces);
      for (unsigned short k = 0; k < volElem[indV].nFaces; ++k)
        volElem[indV].JacFacesIsConsideredConstant[k] = (bool)shortRecvBuf[i][indS++];

      /* Give the member variables that are not obtained via communication their
         values. Some of these variables are not used for halo elements. */
      volElem[indV].elemIsOwned = false;
      volElem[indV].JacIsConsideredConstant = false;
      volElem[indV].offsetDOFsSolGlobal = ULONG_MAX;

      /* Halo elements do not own a face per definition. */
      volElem[indV].ElementOwnsFaces.assign(volElem[indV].nFaces, false);

      /* Get the length scale from the double receive buffer.*/
      volElem[indV].lenScale = doubleRecvBuf[i][indD++];
    }

    /* Store the information of the points in haloPoints. */
    const long nPointsThisRank = longRecvBuf[i][indL++];
    for (long j = 0; j < nPointsThisRank; ++j) {
      CPointFEM thisPoint;
      thisPoint.globalID = longRecvBuf[i][indL++];
      thisPoint.periodIndexToDonor = shortRecvBuf[i][indS++];
      for (unsigned short l = 0; l < nDim; ++l) thisPoint.coor[l] = doubleRecvBuf[i][indD++];

      haloPoints.push_back(thisPoint);
    }

    /* The communication buffers from this rank are not needed anymore.
       Delete them using the swap function. */
    vector<short>().swap(shortRecvBuf[i]);
    vector<long>().swap(longRecvBuf[i]);
    vector<su2double>().swap(doubleRecvBuf[i]);
  }

  /* Remove the duplicate entries from haloPoints. */
  sort(haloPoints.begin(), haloPoints.end());
  lastPoint = unique(haloPoints.begin(), haloPoints.end());
  haloPoints.erase(lastPoint, haloPoints.end());

  /* Initialization of some variables to sort the halo points. */
  Global_nPoint = geometry->GetGlobal_nPoint();
  unsigned long InvalidPointID = Global_nPoint + 10;
  short InvalidPerInd = SHRT_MAX;

  /*--- Search for the nonperiodic halo points in the local points to see
        if these points are already stored on this rank. If this is the
        case invalidate this halo and decrease the number of halo points.
        Afterwards remove the invalid halos from the vector.       ---*/
  unsigned long nHaloPoints = haloPoints.size();
  for (unsigned long i = 0; i < haloPoints.size(); ++i) {
    if (haloPoints[i].periodIndexToDonor != -1) break;  // Test for nonperiodic.

    if (binary_search(meshPoints.begin(), meshPoints.end(), haloPoints[i])) {
      haloPoints[i].globalID = InvalidPointID;
      haloPoints[i].periodIndexToDonor = InvalidPerInd;
      --nHaloPoints;
    }
  }

  sort(haloPoints.begin(), haloPoints.end());
  haloPoints.resize(nHaloPoints);

  /* Increase the capacity of meshPoints, such that the halo points can be
     stored in there as well. Note that in case periodic points are present
     this is an upper bound. Add the non-periodic halo points to meshPoints. */
  meshPoints.reserve(meshPoints.size() + nHaloPoints);

  for (unsigned long i = 0; i < haloPoints.size(); ++i) {
    if (haloPoints[i].periodIndexToDonor != -1) break;  // Test for nonperiodic.

    meshPoints.push_back(haloPoints[i]);
  }

  /* Create a map from the global point ID and periodic index to the local
     index in the vector meshPoints. First store the points already present
     in meshPoints. */
  map<CUnsignedLong2T, unsigned long> mapGlobalPointIDToInd;
  for (unsigned long i = 0; i < meshPoints.size(); ++i) {
    CUnsignedLong2T globIndAndPer;
    globIndAndPer.long0 = meshPoints[i].globalID;
    globIndAndPer.long1 = meshPoints[i].periodIndexToDonor + 1;  // Note the +1 again.

    mapGlobalPointIDToInd[globIndAndPer] = i;
  }

  /*--- Convert the global indices in the boundary connectivities to local ones.
        Note that the volume ID's already contain the local number. ---*/
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    for (unsigned long i = 0; i < boundaries[iMarker].surfElem.size(); ++i) {
      /* Convert the global node ID's to local values. Note that for these node
         ID's no periodic transformation can be present. */
      for (unsigned short j = 0; j < boundaries[iMarker].surfElem[i].nDOFsGrid; ++j) {
        CUnsignedLong2T searchItem(boundaries[iMarker].surfElem[i].nodeIDsGrid[j], 0);
        map<CUnsignedLong2T, unsigned long>::const_iterator LLMI;
        LLMI = mapGlobalPointIDToInd.find(searchItem);
        boundaries[iMarker].surfElem[i].nodeIDsGrid[j] = LLMI->second;
      }
    }
  }

  /*--- The only halo points that must be added to meshPoints are the periodic
        halo points. It must be checked whether or not the periodic points in
        haloPoints match with points in meshPoints. This is done below. ---*/
  for (unsigned long iLow = 0; iLow < haloPoints.size();) {
    /* Determine the upper index for this periodic transformation. */
    unsigned long iUpp;
    for (iUpp = iLow + 1; iUpp < haloPoints.size(); ++iUpp)
      if (haloPoints[iUpp].periodIndexToDonor != haloPoints[iLow].periodIndexToDonor) break;

    /* Check for a true periodic index. */
    short perIndex = haloPoints[iLow].periodIndexToDonor;
    if (perIndex != -1) {
      /* Easier storage of the surface elements. */
      vector<CSurfaceElementFEM>& surfElem = boundaries[perIndex].surfElem;

      /*--- In the loop below the coordinates of the points of this local
            periodic boundary as well as a matching tolerance are determined.
            A vector of point ID's is also created, which is needed later on
            when it is checked whether or not a matching point is already
            stored in meshPoints. ---*/
      vector<long> indInPoints(meshPoints.size(), -1);
      vector<unsigned long> IDsPoints;
      vector<su2double> coordPoints;
      vector<su2double> tolPoints;

      for (unsigned long j = 0; j < surfElem.size(); ++j) {
        /* Determine the tolerance for equal points, which is a small value
           times the length scale of the adjacent volume element. */
        const su2double tolElem = 1.e-2 * volElem[surfElem[j].volElemID].lenScale;

        /* Loop over the nodes of this surface grid and update the points on
           this periodic boundary. */
        for (unsigned short k = 0; k < surfElem[j].nDOFsGrid; ++k) {
          unsigned long nn = surfElem[j].nodeIDsGrid[k];

          if (indInPoints[nn] == -1) {
            /* Point is not stored yet in pointsBoundary. Do so now. */
            indInPoints[nn] = IDsPoints.size();
            IDsPoints.push_back(nn);
            tolPoints.push_back(tolElem);

            for (unsigned short l = 0; l < nDim; ++l) coordPoints.push_back(meshPoints[nn].coor[l]);
          } else {
            /* Point is already stored. Update the tolerance. */
            nn = indInPoints[nn];
            tolPoints[nn] = min(tolPoints[nn], tolElem);
          }
        }
      }

      /* Create a local ADT of the points on the periodic boundary. */
      CADTPointsOnlyClass periodicADT(nDim, IDsPoints.size(), coordPoints.data(), IDsPoints.data(), false);

      /* Get the data for the periodic transformation to the donor. */
      auto center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(perIndex));
      auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(perIndex));
      auto trans = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(perIndex));

      /*--- Compute the rotation matrix and translation vector for the
            transformation from the donor. This is the transpose of the
            transformation to the donor. ---*/

      /* Store (center-trans) as it is constant and will be added on. */
      su2double translation[] = {center[0] - trans[0], center[1] - trans[1], center[2] - trans[2]};

      /* Store angles separately for clarity. Compute sines/cosines. */
      su2double theta = angles[0];
      su2double phi = angles[1];
      su2double psi = angles[2];

      su2double cosTheta = cos(theta), cosPhi = cos(phi), cosPsi = cos(psi);
      su2double sinTheta = sin(theta), sinPhi = sin(phi), sinPsi = sin(psi);

      /* Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis, then z-axis. */
      su2double rotMatrix[3][3];
      rotMatrix[0][0] = cosPhi * cosPsi;
      rotMatrix[0][1] = cosPhi * sinPsi;
      rotMatrix[0][2] = -sinPhi;

      rotMatrix[1][0] = sinTheta * sinPhi * cosPsi - cosTheta * sinPsi;
      rotMatrix[1][1] = sinTheta * sinPhi * sinPsi + cosTheta * cosPsi;
      rotMatrix[1][2] = sinTheta * cosPhi;

      rotMatrix[2][0] = cosTheta * sinPhi * cosPsi + sinTheta * sinPsi;
      rotMatrix[2][1] = cosTheta * sinPhi * sinPsi - sinTheta * cosPsi;
      rotMatrix[2][2] = cosTheta * cosPhi;

      /* Loop over the halo points for this periodic transformation. */
      for (unsigned long i = iLow; i < iUpp; ++i) {
        /* Apply the periodic transformation to the coordinates
           stored in this halo point. */
        su2double dx = haloPoints[i].coor[0] - center[0];
        su2double dy = haloPoints[i].coor[1] - center[1];
        su2double dz = nDim == 3 ? haloPoints[i].coor[2] - center[2] : su2double(0.0);

        haloPoints[i].coor[0] = rotMatrix[0][0] * dx + rotMatrix[0][1] * dy + rotMatrix[0][2] * dz + translation[0];
        haloPoints[i].coor[1] = rotMatrix[1][0] * dx + rotMatrix[1][1] * dy + rotMatrix[1][2] * dz + translation[1];
        haloPoints[i].coor[2] = rotMatrix[2][0] * dx + rotMatrix[2][1] * dy + rotMatrix[2][2] * dz + translation[2];

        /* Search for the nearest coordinate in the ADT. */
        su2double dist;
        unsigned long pointID;
        int rankID;

        periodicADT.DetermineNearestNode(haloPoints[i].coor, dist, pointID, rankID);

        /* Check whether the distance is less equal to the tolerance for
           a matching point. */
        const unsigned long nn = indInPoints[pointID];
        if (dist <= tolPoints[nn]) {
          /* The distance to the nearest point is less than the tolerance,
             hence this periodically transformed point is present on the
             boundary. Store it as such in the map mapGlobalPointIDToInd. */
          CUnsignedLong2T globIndAndPer;
          globIndAndPer.long0 = haloPoints[i].globalID;
          globIndAndPer.long1 = haloPoints[i].periodIndexToDonor + 1;  // Note the +1 again.

          mapGlobalPointIDToInd[globIndAndPer] = pointID;
        } else {
          /* The distance to the nearest point is larger than the tolerance,
             hence this periodically transformed point is not present yet on
             this rank. Store it in the mapping to the local points and
             create it in meshPoints. */
          CUnsignedLong2T globIndAndPer;
          globIndAndPer.long0 = haloPoints[i].globalID;
          globIndAndPer.long1 = haloPoints[i].periodIndexToDonor + 1;  // Note the +1 again.

          mapGlobalPointIDToInd[globIndAndPer] = meshPoints.size();
          meshPoints.push_back(haloPoints[i]);
        }
      }
    }

    /* Set iLow to iUpp for the next periodic transformation. */
    iLow = iUpp;
  }

  /*--- Convert the global node numbering in the elements to a local numbering and
        determine the value of factTimeLevel. This is the number of local time steps
        of the element relative to the largest time step of an element in the mesh.
        This value can only differ from 1 when time accurate local time stepping is
        used. ---*/
  for (unsigned long i = 0; i < nVolElemTot; ++i) {
    for (unsigned short j = 0; j < volElem[i].nDOFsGrid; ++j) {
      CUnsignedLong2T searchItem(volElem[i].nodeIDsGrid[j],
                                 volElem[i].periodIndexToDonor + 1);  // Again the +1.
      map<CUnsignedLong2T, unsigned long>::const_iterator LLMI;
      LLMI = mapGlobalPointIDToInd.find(searchItem);
      volElem[i].nodeIDsGrid[j] = LLMI->second;
    }

    const unsigned short diffTimeLevel = nTimeLevels - 1 - volElem[i].timeLevel;
    volElem[i].factTimeLevel = pow(2, diffTimeLevel);
  }

  /* Determine the number of halo elements per time level in cumulative
     storage format. */
  nVolElemHaloPerTimeLevel.assign(nTimeLevels + 1, 0);
  for (unsigned long i = nVolElemOwned; i < nVolElemTot; ++i) ++nVolElemHaloPerTimeLevel[volElem[i].timeLevel + 1];

  nVolElemHaloPerTimeLevel[0] = nVolElemOwned;
  for (unsigned short i = 0; i < nTimeLevels; ++i) nVolElemHaloPerTimeLevel[i + 1] += nVolElemHaloPerTimeLevel[i];
}

void CMeshFEM::ComputeGradientsCoordinatesFace(const unsigned short nIntegration, const unsigned short nDOFs,
                                               const su2double* matDerBasisInt, const unsigned long* DOFs,
                                               su2double* derivCoor, CConfig* config) {
  /* Allocate the memory to store the values of dxdr, dydr, etc. */
  vector<su2double> helpDxdrVec(nIntegration * nDim * nDim);
  su2double* dxdrVec = helpDxdrVec.data();

  /* Determine the gradients of the Cartesian coordinates w.r.t. the
     parametric coordinates. */
  ComputeGradientsCoorWRTParam(nIntegration, nDOFs, matDerBasisInt, DOFs, dxdrVec, config);

  /* Make a distinction between 2D and 3D to compute the derivatives drdx,
     drdy, etc. */
  switch (nDim) {
    case 2: {
      /* 2D computation. Store the offset between the r and s derivatives. */
      const unsigned short off = 2 * nIntegration;

      /* Loop over the integration points. */
      unsigned short ii = 0;
      for (unsigned short j = 0; j < nIntegration; ++j) {
        /* Retrieve the values of dxdr, dydr, dxds and dyds from dxdrVec
           in this integration point. */
        const unsigned short jx = 2 * j;
        const unsigned short jy = jx + 1;
        const su2double dxdr = dxdrVec[jx], dydr = dxdrVec[jy];
        const su2double dxds = dxdrVec[jx + off], dyds = dxdrVec[jy + off];

        /* Compute the inverse relations drdx, drdy, dsdx, dsdy. */
        const su2double Jinv = 1.0 / (dxdr * dyds - dxds * dydr);

        derivCoor[ii++] = dyds * Jinv;   // drdx
        derivCoor[ii++] = -dxds * Jinv;  // drdy
        derivCoor[ii++] = -dydr * Jinv;  // dsdx
        derivCoor[ii++] = dxdr * Jinv;   // dsdy
      }
      break;
    }

    case 3: {
      /* 3D computation. Store the offset between the r and s and r and t derivatives. */
      const unsigned short offS = 3 * nIntegration, offT = 6 * nIntegration;

      /* Loop over the integration points. */
      unsigned short ii = 0;
      for (unsigned short j = 0; j < nIntegration; ++j) {
        /* Retrieve the values of dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt
           and dzdt from dxdrVec in this integration point. */
        const unsigned short jx = 3 * j;
        const unsigned short jy = jx + 1, jz = jx + 2;
        const su2double dxdr = dxdrVec[jx], dydr = dxdrVec[jy], dzdr = dxdrVec[jz];
        const su2double dxds = dxdrVec[jx + offS], dyds = dxdrVec[jy + offS], dzds = dxdrVec[jz + offS];
        const su2double dxdt = dxdrVec[jx + offT], dydt = dxdrVec[jy + offT], dzdt = dxdrVec[jz + offT];

        /* Compute the inverse relations drdx, drdy, drdz, dsdx, dsdy, dsdz,
           dtdx, dtdy, dtdz. */
        const su2double Jinv = 1.0 / (dxdr * (dyds * dzdt - dzds * dydt) - dxds * (dydr * dzdt - dzdr * dydt) +
                                      dxdt * (dydr * dzds - dzdr * dyds));

        derivCoor[ii++] = (dyds * dzdt - dzds * dydt) * Jinv;  // drdx
        derivCoor[ii++] = (dzds * dxdt - dxds * dzdt) * Jinv;  // drdy
        derivCoor[ii++] = (dxds * dydt - dyds * dxdt) * Jinv;  // drdz

        derivCoor[ii++] = (dzdr * dydt - dydr * dzdt) * Jinv;  // dsdx
        derivCoor[ii++] = (dxdr * dzdt - dzdr * dxdt) * Jinv;  // dsdy
        derivCoor[ii++] = (dydr * dxdt - dxdr * dydt) * Jinv;  // dsdz

        derivCoor[ii++] = (dydr * dzds - dzdr * dyds) * Jinv;  // dtdx
        derivCoor[ii++] = (dzdr * dxds - dxdr * dzds) * Jinv;  // dtdy
        derivCoor[ii++] = (dxdr * dyds - dydr * dxds) * Jinv;  // dtdz
      }

      break;
    }
  }
}

void CMeshFEM::ComputeGradientsCoorWRTParam(const unsigned short nIntegration, const unsigned short nDOFs,
                                            const su2double* matDerBasisInt, const unsigned long* DOFs,
                                            su2double* derivCoor, CConfig* config) {
  /* Allocate the memory to store the coordinates as right hand side. */
  vector<su2double> vecRHS(nDOFs * nDim);

  /* Loop over the grid DOFs of the element and copy the coordinates in
     vecRHS in row major order. */
  unsigned long ic = 0;
  for (unsigned short j = 0; j < nDOFs; ++j) {
    for (unsigned short k = 0; k < nDim; ++k, ++ic) vecRHS[ic] = meshPoints[DOFs[j]].coor[k];
  }

  /* Carry out the matrix matrix product The last argument is NULL, such
     that this gemm call is ignored in the profiling. Replace by config if
     if should be included. */
  blasFunctions->gemm(nDim * nIntegration, nDim, nDOFs, matDerBasisInt, vecRHS.data(), derivCoor, nullptr);
}

void CMeshFEM::ComputeNormalsFace(const unsigned short nIntegration, const unsigned short nDOFs, const su2double* dr,
                                  const su2double* ds, const unsigned long* DOFs, su2double* normals) {
  /* Initialize the counter ii to 0. ii is the index in normals where the
     information is stored. */
  unsigned int ii = 0;

  /* Make a distinction between 2D and 3D. */
  switch (nDim) {
    case 2: {
      /* 2D computation. Loop over the integration points of the face. */
      for (unsigned short j = 0; j < nIntegration; ++j) {
        /*--- Loop over the number of DOFs of the face to compute dxdr
              and dydr. ---*/
        const su2double* drr = &dr[j * nDOFs];
        su2double dxdr = 0.0, dydr = 0.0;
        for (unsigned short k = 0; k < nDOFs; ++k) {
          dxdr += drr[k] * meshPoints[DOFs[k]].coor[0];
          dydr += drr[k] * meshPoints[DOFs[k]].coor[1];
        }

        /* Determine the length of the tangential vector (dxdr, dydr), which
           is also the length of the corresponding normal vector. Also compute
           the inverse of the length. Make sure that a division by zero is
           avoided, although this is most likely never active. */
        const su2double lenNorm = sqrt(dxdr * dxdr + dydr * dydr);
        const su2double invLenNorm = lenNorm < su2double(1.e-35) ? su2double(1.e+35) : 1.0 / lenNorm;

        /* Store the corresponding unit normal vector and its length. The
           direction of the normal vector is such that it is outward pointing
           for the element on side 0 of the face. */
        normals[ii++] = dydr * invLenNorm;
        normals[ii++] = -dxdr * invLenNorm;
        normals[ii++] = lenNorm;
      }

      break;
    }

    case 3: {
      /* 3D computation. Loop over the integration points of the face. */
      for (unsigned short j = 0; j < nIntegration; ++j) {
        /*--- Loop over the number of DOFs of the face to compute dxdr,
              dxds, dydr, dyds, dzdr and dzds. ---*/
        const su2double *drr = &dr[j * nDOFs], *dss = &ds[j * nDOFs];
        su2double dxdr = 0.0, dydr = 0.0, dzdr = 0.0, dxds = 0.0, dyds = 0.0, dzds = 0.0;
        for (unsigned short k = 0; k < nDOFs; ++k) {
          dxdr += drr[k] * meshPoints[DOFs[k]].coor[0];
          dydr += drr[k] * meshPoints[DOFs[k]].coor[1];
          dzdr += drr[k] * meshPoints[DOFs[k]].coor[2];

          dxds += dss[k] * meshPoints[DOFs[k]].coor[0];
          dyds += dss[k] * meshPoints[DOFs[k]].coor[1];
          dzds += dss[k] * meshPoints[DOFs[k]].coor[2];
        }

        /* Compute the vector product dxdr X dxds, where x is the coordinate
           vector (x,y,z). Compute the length of this vector, which is an area,
           as well as the inverse.  Make sure that a division by zero is
           avoided, although this is most likely never active. */
        const su2double nx = dydr * dzds - dyds * dzdr;
        const su2double ny = dxds * dzdr - dxdr * dzds;
        const su2double nz = dxdr * dyds - dxds * dydr;

        const su2double lenNorm = sqrt(nx * nx + ny * ny + nz * nz);
        const su2double invLenNorm = lenNorm < su2double(1.e-35) ? su2double(1.e+35) : 1.0 / lenNorm;

        /* Store the components of the unit normal as well as its length in the
           normals. Note that the current direction of the normal is pointing
           into the direction of the element on side 0 of the face. However,
           in the actual computation of the integral over the faces, it is
           assumed that the vector points in the opposite direction.
           Hence the normal vector must be reversed. */
        normals[ii++] = -nx * invLenNorm;
        normals[ii++] = -ny * invLenNorm;
        normals[ii++] = -nz * invLenNorm;
        normals[ii++] = lenNorm;
      }
      break;
    }
  }
}

void CMeshFEM::MetricTermsBoundaryFaces(CBoundaryFEM* boundary, CConfig* config) {
  /* Determine whether or not the viscous terms are needed. */
  const bool viscousTerms = (config->GetKind_Solver() != MAIN_SOLVER::FEM_EULER &&
                             config->GetKind_Solver() != MAIN_SOLVER::DISC_ADJ_FEM_EULER);

  /*--- Loop over the boundary faces stored on this rank. ---*/
  for (unsigned long i = 0; i < boundary->surfElem.size(); ++i) {
    /*--------------------------------------------------------------------------*/
    /*--- Step 1: Allocate the memory for the face metric terms.             ---*/
    /*---         Depending on the case, not all of this memory is needed.   ---*/
    /*---         - Unit normals + area (nDim+1 per integration point)       ---*/
    /*---         - drdx, dsdx, etc. (nDim*nDim per integration point)       ---*/
    /*---         - grid velocities in the integration points.               ---*/
    /*--------------------------------------------------------------------------*/

    /* Determine the corresponding standard face element and get the
       relevant information from it. */
    const unsigned short ind = boundary->surfElem[i].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Allocate the several metric terms. ---*/
    boundary->surfElem[i].metricNormalsFace.resize(nInt * (nDim + 1));

    if (viscousTerms) boundary->surfElem[i].metricCoorDerivFace.resize(nInt * nDim * nDim);

    /* Allocate the memory for the grid velocities and initialize them to the
       default value of zero. */
    boundary->surfElem[i].gridVelocities.assign(nInt * nDim, 0.0);

    /*--------------------------------------------------------------------------*/
    /*--- Step 2: Determine the actual metric data in the integration points ---*/
    /*---         of the faces.                                              ---*/
    /*--------------------------------------------------------------------------*/

    /* Call the function ComputeNormalsFace to compute the unit normals and
       its corresponding area in the integration points. */
    unsigned short nDOFs = standardBoundaryFacesGrid[ind].GetNDOFsFace();
    const su2double* dr = standardBoundaryFacesGrid[ind].GetDrBasisFaceIntegration();
    const su2double* ds = standardBoundaryFacesGrid[ind].GetDsBasisFaceIntegration();

    ComputeNormalsFace(nInt, nDOFs, dr, ds, boundary->surfElem[i].DOFsGridFace.data(),
                       boundary->surfElem[i].metricNormalsFace.data());

    /* Compute the derivatives of the parametric coordinates w.r.t. the
       Cartesian coordinates, i.e. drdx, drdy, etc. in the integration points
       of the face, if needed. */
    if (viscousTerms) {
      nDOFs = standardBoundaryFacesGrid[ind].GetNDOFsElem();
      dr = standardBoundaryFacesGrid[ind].GetMatDerBasisElemIntegration();

      ComputeGradientsCoordinatesFace(nInt, nDOFs, dr, boundary->surfElem[i].DOFsGridElement.data(),
                                      boundary->surfElem[i].metricCoorDerivFace.data(), config);
    }
  }
}

void CMeshFEM::SetPositive_ZArea(CConfig* config) {
  /*---------------------------------------------------------------------------*/
  /*--- Step 1: Determine the local contribution to the positive z area.    ---*/
  /*---------------------------------------------------------------------------*/

  /* Loop over the boundary markers. */
  su2double PositiveZArea = 0.0;
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    /* Determine whether or not this boundary contributes. */
    const unsigned short Boundary = config->GetMarker_All_KindBC(iMarker);
    const unsigned short Monitoring = config->GetMarker_All_Monitoring(iMarker);

    if (((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL) ||
         (Boundary == LOAD_BOUNDARY) || (Boundary == DISPLACEMENT_BOUNDARY)) &&
        (Monitoring == YES)) {
      /* Easier storage of the surface elements for this marker. */
      const vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;

      /* Loop over the surface elements. */
      for (unsigned long i = 0; i < surfElem.size(); ++i) {
        /* Determine the number of integration points and their weights via
           the corresponding standard element. */
        const unsigned short ind = surfElem[i].indStandardElement;
        const unsigned short nInt = standardBoundaryFacesGrid[ind].GetNIntegration();
        const su2double* weights = standardBoundaryFacesGrid[ind].GetWeightsIntegration();

        /* Loop over the integration points for this element and update PositiveZArea
           if the normal has a negative z-component. In that case it will give a
           positive contribution to PositiveZArea as this must take the normal pointing
           into the element into account. */
        for (unsigned short j = 0; j < nInt; ++j) {
          /* Store the normal data for this integration point a bit easier and update
             PositiveZArea, if needed. */
          const su2double* normal = &surfElem[i].metricNormalsFace[j * (nDim + 1)];
          if (normal[nDim - 1] < 0.0) PositiveZArea -= weights[j] * normal[nDim - 1] * normal[nDim];
        }
      }
    }
  }

  /*---------------------------------------------------------------------------*/
  /*--- Step 2: Perform an Allreduce such that the global value of          ---*/
  /*---         PositiveZArea is known on all ranks.                        ---*/
  /*---------------------------------------------------------------------------*/

#ifdef HAVE_MPI
  su2double locArea = PositiveZArea;
  SU2_MPI::Allreduce(&locArea, &PositiveZArea, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*---------------------------------------------------------------------------*/
  /*--- Step 3: Set the reference area, if this was not specified and write ---*/
  /*---         a message about the projected area if I am the master.      ---*/
  /*---------------------------------------------------------------------------*/

  if (config->GetRefArea() == 0.0) config->SetRefArea(PositiveZArea);

  if (rank == MASTER_NODE) {
    if (nDim == 2)
      cout << "Area projection in the y-plane = " << PositiveZArea << "." << endl;
    else
      cout << "Area projection in the z-plane = " << PositiveZArea << "." << endl;
  }
}

CMeshFEM_DG::CMeshFEM_DG(CGeometry* geometry, CConfig* config) : CMeshFEM(geometry, config) {}

void CMeshFEM_DG::SetGlobal_to_Local_Point() {
  Global_to_Local_Point.clear();
  unsigned long ii = 0;
  for (unsigned long i = 0; i < nVolElemOwned; ++i) {
    for (unsigned short j = 0; j < volElem[i].nDOFsSol; ++j, ++ii) {
      Global_to_Local_Point[volElem[i].offsetDOFsSolGlobal + j] = ii;
    }
  }
}

void CMeshFEM_DG::CoordinatesIntegrationPoints() {
  /*--------------------------------------------------------------------*/
  /*--- Step 1: The integration points of the owned volume elements. ---*/
  /*--------------------------------------------------------------------*/

  /*--- Loop over the owned elements to compute the coordinates
        in the integration points. ---*/
  for (unsigned long l = 0; l < nVolElemOwned; ++l) {
    /* Get the required data from the corresponding standard element. */
    const unsigned short ind = volElem[l].indStandardElement;
    const unsigned short nInt = standardElementsGrid[ind].GetNIntegration();
    const su2double* lag = standardElementsGrid[ind].GetBasisFunctionsIntegration();

    /* Allocate the memory for the coordinates of the integration points. */
    volElem[l].coorIntegrationPoints.resize(nDim * nInt);

    /* Store the grid DOFs of this element a bit easier. */
    const unsigned short nDOFs = volElem[l].nDOFsGrid;
    const unsigned long* DOFs = volElem[l].nodeIDsGrid.data();

    /*--- Loop over the integration points of this element and compute
          the coordinates of the integration points. */
    for (unsigned short i = 0; i < nInt; ++i) {
      su2double* coor = volElem[l].coorIntegrationPoints.data() + i * nDim;
      const unsigned short jj = i * nDOFs;
      for (unsigned short j = 0; j < nDim; ++j) {
        coor[j] = 0.0;
        for (unsigned short k = 0; k < nDOFs; ++k) coor[j] += lag[jj + k] * meshPoints[DOFs[k]].coor[j];
      }
    }
  }

  /*----------------------------------------------------------------------*/
  /*--- Step 2: The integration points of the internal matching faces. ---*/
  /*----------------------------------------------------------------------*/

  /*--- Loop over the internal matching faces to compute the coordinates
        in the integration points. ---*/
  for (unsigned long l = 0; l < matchingFaces.size(); ++l) {
    /* Get the required data from the corresponding standard element. */
    const unsigned short ind = matchingFaces[l].indStandardElement;
    const unsigned short nInt = standardMatchingFacesGrid[ind].GetNIntegration();
    const unsigned short nDOFs = standardMatchingFacesGrid[ind].GetNDOFsFaceSide0();
    const su2double* lag = standardMatchingFacesGrid[ind].GetBasisFaceIntegrationSide0();

    /* Allocate the memory for the coordinates of the integration points. */
    matchingFaces[l].coorIntegrationPoints.resize(nDim * nInt);

    /* Store the grid DOFs of this face a bit easier. */
    const unsigned long* DOFs = matchingFaces[l].DOFsGridFaceSide0.data();

    /*--- Loop over the integration points of this face and compute
          the coordinates of the integration points. */
    for (unsigned short i = 0; i < nInt; ++i) {
      su2double* coor = matchingFaces[l].coorIntegrationPoints.data() + i * nDim;
      const unsigned short jj = i * nDOFs;
      for (unsigned short j = 0; j < nDim; ++j) {
        coor[j] = 0.0;
        for (unsigned short k = 0; k < nDOFs; ++k) coor[j] += lag[jj + k] * meshPoints[DOFs[k]].coor[j];
      }
    }
  }

  /*-------------------------------------------------------------*/
  /*--- Step 3: The integration points of the boundary faces. ---*/
  /*-------------------------------------------------------------*/

  /*--- Loop over the boundary markers. Make sure to exclude the periodic
        boundaries, because these are not physical. ---*/
  for (unsigned short iMarker = 0; iMarker < boundaries.size(); ++iMarker) {
    if (!boundaries[iMarker].periodicBoundary) {
      /* Abbreviate the surface elements a bit easier. */
      vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;

      /* Loop over the boundary faces and determine the coordinates
         in the integration points. */
      for (unsigned long l = 0; l < surfElem.size(); ++l) {
        /* Get the required data from the corresponding standard element. */
        const unsigned short ind = surfElem[l].indStandardElement;
        const unsigned short nInt = standardBoundaryFacesGrid[ind].GetNIntegration();
        const unsigned short nDOFs = surfElem[l].nDOFsGrid;
        const su2double* lag = standardBoundaryFacesGrid[ind].GetBasisFaceIntegration();

        /* Allocate the memory for the coordinates of the integration points. */
        surfElem[l].coorIntegrationPoints.resize(nDim * nInt);

        /* Store the grid DOFs of this face a bit easier. */
        const unsigned long* DOFs = surfElem[l].DOFsGridFace.data();

        /*--- Loop over the integration points of this face and compute
              the coordinates of the integration points. */
        for (unsigned short i = 0; i < nInt; ++i) {
          su2double* coor = surfElem[l].coorIntegrationPoints.data() + i * nDim;
          const unsigned short jj = i * nDOFs;
          for (unsigned short j = 0; j < nDim; ++j) {
            coor[j] = 0.0;
            for (unsigned short k = 0; k < nDOFs; ++k) coor[j] += lag[jj + k] * meshPoints[DOFs[k]].coor[j];
          }
        }
      }
    }
  }
}

void CMeshFEM_DG::CoordinatesSolDOFs() {
  /*--- Loop over the owned elements to compute the coordinates
        of the solution DOFs. ---*/
  for (unsigned long l = 0; l < nVolElemOwned; ++l) {
    /* Store the grid DOFs and the number of solDOFS a bit easier. */
    const unsigned short nDOFsGrid = volElem[l].nDOFsGrid;
    const unsigned short nDOFsSol = volElem[l].nDOFsSol;
    const unsigned long* DOFsGrid = volElem[l].nodeIDsGrid.data();

    /* Get the required data from the corresponding standard element. */
    const unsigned short ind = volElem[l].indStandardElement;
    const su2double* lag = standardElementsGrid[ind].GetBasisFunctionsSolDOFs();

    /* Allocate the memory for the coordinates of the solution DOFs. */
    volElem[l].coorSolDOFs.resize(nDim * nDOFsSol);

    /*--- Loop over the solution DOFs of this element and compute the
          coordinates. Note that in case the polynomial degree of the grid
          and solution is the same the grid DOFs and solution DOFs coincide
          and the interpolation boils down to a copy of the coordinates. ---*/
    for (unsigned short i = 0; i < nDOFsSol; ++i) {
      su2double* coor = volElem[l].coorSolDOFs.data() + i * nDim;
      const unsigned short jj = i * nDOFsGrid;
      for (unsigned short j = 0; j < nDim; ++j) {
        coor[j] = 0.0;
        for (unsigned short k = 0; k < nDOFsGrid; ++k) coor[j] += lag[jj + k] * meshPoints[DOFsGrid[k]].coor[j];
      }
    }
  }
}

void CMeshFEM_DG::CreateFaces(CConfig* config) {
  /*---------------------------------------------------------------------------*/
  /*--- Step 1: Determine the faces of the locally stored part of the grid. ---*/
  /*---------------------------------------------------------------------------*/

  /*--- Loop over the volume elements stored on this rank, including the halos. ---*/
  vector<CFaceOfElement> localFaces;

  for (unsigned long k = 0; k < nVolElemTot; ++k) {
    /* Determine the corner points of all the faces of this element. */
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long faceConn[6][4];

    volElem[k].GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /* Loop over the faces of this element, set the appropriate information,
       create a unique numbering and add the faces to localFaces. */
    for (unsigned short i = 0; i < nFaces; ++i) {
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for (unsigned short j = 0; j < nPointsPerFace[i]; ++j) thisFace.cornerPoints[j] = faceConn[i][j];

      /* Copy the data from the volume element. */
      thisFace.elemID0 = k;
      thisFace.nPolyGrid0 = volElem[k].nPolyGrid;
      thisFace.nPolySol0 = volElem[k].nPolySol;
      thisFace.nDOFsElem0 = volElem[k].nDOFsSol;
      thisFace.elemType0 = volElem[k].VTK_Type;
      thisFace.faceID0 = i;
      thisFace.faceIndicator = -2;  // Initialized to an invalid face.

      thisFace.JacFaceIsConsideredConstant = volElem[k].JacFacesIsConsideredConstant[i];
      thisFace.elem0IsOwner = volElem[k].ElementOwnsFaces[i];

      /* Determine whether or not the integration rule for straight and curved
         faces is the same. If that is the case, the value of JacFaceIsConsideredConstant
         is set to false. Hence it is only needed to carry out this check for faces
         with a constant Jacobian. This is done to reduce the number of standard elements. */
      if (thisFace.JacFaceIsConsideredConstant) {
        auto orderExactStraight = (unsigned short)ceil(thisFace.nPolyGrid0 * config->GetQuadrature_Factor_Straight());
        auto orderExactCurved = (unsigned short)ceil(thisFace.nPolyGrid0 * config->GetQuadrature_Factor_Curved());

        if (orderExactStraight == orderExactCurved) {
          orderExactStraight = (unsigned short)ceil(thisFace.nPolySol0 * config->GetQuadrature_Factor_Straight());
          orderExactCurved = (unsigned short)ceil(thisFace.nPolySol0 * config->GetQuadrature_Factor_Curved());
          if (orderExactStraight == orderExactCurved) thisFace.JacFaceIsConsideredConstant = false;
        }
      }

      /* Renumber the corner points of the face, but keep the orientation.
         Afterwards, add it to localFaces. */
      thisFace.CreateUniqueNumberingWithOrientation();

      localFaces.push_back(thisFace);
    }
  }

  /*--- Sort the the local faces in increasing order. ---*/
  sort(localFaces.begin(), localFaces.end());

  /*--- Loop over the faces to merge the matching faces. ---*/
  for (unsigned long i = 1; i < localFaces.size(); ++i) {
    /* Check for a matching face with the previous face in the vector.
       Note that the == operator only checks the node IDs. */
    if (localFaces[i] == localFaces[i - 1]) {
      /* Faces are matching. Check if it should be kept, i.e. if one
         of the elements owns the face. */
      if (localFaces[i].elem0IsOwner || localFaces[i - 1].elem0IsOwner) {
        /* Store the data for this matching face in faces[i-1]. */
        if (localFaces[i].elemID0 < nVolElemTot) {
          localFaces[i - 1].elemID0 = localFaces[i].elemID0;
          localFaces[i - 1].nPolyGrid0 = localFaces[i].nPolyGrid0;
          localFaces[i - 1].nPolySol0 = localFaces[i].nPolySol0;
          localFaces[i - 1].nDOFsElem0 = localFaces[i].nDOFsElem0;
          localFaces[i - 1].elemType0 = localFaces[i].elemType0;
          localFaces[i - 1].faceID0 = localFaces[i].faceID0;
          localFaces[i - 1].elem0IsOwner = localFaces[i].elem0IsOwner;
        } else {
          localFaces[i - 1].elemID1 = localFaces[i].elemID1;
          localFaces[i - 1].nPolyGrid1 = localFaces[i].nPolyGrid1;
          localFaces[i - 1].nPolySol1 = localFaces[i].nPolySol1;
          localFaces[i - 1].nDOFsElem1 = localFaces[i].nDOFsElem1;
          localFaces[i - 1].elemType1 = localFaces[i].elemType1;
          localFaces[i - 1].faceID1 = localFaces[i].faceID1;
        }

        /* Adapt the boolean to indicate whether or not the face has a constant
           Jacobian of the transformation, although in principle this info
           should be the same for both faces. */
        if (!(localFaces[i - 1].JacFaceIsConsideredConstant && localFaces[i].JacFaceIsConsideredConstant)) {
          localFaces[i - 1].JacFaceIsConsideredConstant = false;
        }

        /* Set this face indicator to -1 to indicate an internal face
           and set elem0IsOwner for localFaces[i] to false. */
        localFaces[i - 1].faceIndicator = -1;
        localFaces[i].elem0IsOwner = false;
      }
    }
  }

  /*--- Loop over the boundary markers and its boundary elements to search for
        the corresponding faces in localFaces. These faces should be found.
        Note that periodic boundaries are skipped, because these are treated
        via the halo elements, which are already in place. ---*/
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    if (!boundaries[iMarker].periodicBoundary) {
      for (unsigned long k = 0; k < boundaries[iMarker].surfElem.size(); ++k) {
        /* Determine the corner points of the face of this element. */
        unsigned short nPointsPerFace;
        unsigned long faceConn[4];

        boundaries[iMarker].surfElem[k].GetCornerPointsFace(nPointsPerFace, faceConn);

        /* Create an object of CFaceOfElement to carry out the search. */
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace;
        for (unsigned short j = 0; j < nPointsPerFace; ++j) thisFace.cornerPoints[j] = faceConn[j];
        thisFace.CreateUniqueNumberingWithOrientation();

        /* Search for thisFace in localFaces. It must be found. */
        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        bool thisFaceFound = false;
        if (low != localFaces.end()) {
          if (!(thisFace < *low)) thisFaceFound = true;
        }

        if (thisFaceFound) {
          low->faceIndicator = iMarker;

          /* A few additional checks. */
          bool side0IsBoundary = low->elemID0 < nVolElemTot;
          unsigned long elemID = side0IsBoundary ? low->elemID0 : low->elemID1;
          unsigned short nPoly = side0IsBoundary ? low->nPolyGrid0 : low->nPolyGrid1;

          if (elemID != boundaries[iMarker].surfElem[k].volElemID || nPoly != boundaries[iMarker].surfElem[k].nPolyGrid)
            SU2_MPI::Error(string("Element ID and/or polynomial degree do not match ") +
                               string("for this boundary element. This should not happen."),
                           CURRENT_FUNCTION);
        } else
          SU2_MPI::Error("Boundary face not found in localFaces. This should not happen.", CURRENT_FUNCTION);
      }
    }
  }

  /*--- It is possible that owned non-matching faces are present in the list.
        These faces are indicated by an owned face and a faceIndicator of -2.
        To avoid that these faces are removed afterwards, set their
        faceIndicator to -1. ---*/
  for (unsigned long i = 0; i < localFaces.size(); ++i) {
    if (localFaces[i].faceIndicator == -2 && localFaces[i].elem0IsOwner) localFaces[i].faceIndicator = -1;
  }

  /*--- Remove the invalid faces. This is accomplished by giving the face four
        points and global node ID's that are larger than the largest local point
        ID in the grid. In this way the sorting operator puts these faces at the
        end of the vector, see also the < operator of CFaceOfElement. ---*/
  unsigned long nFacesLoc = localFaces.size();
  for (unsigned long i = 0; i < localFaces.size(); ++i) {
    if (localFaces[i].faceIndicator == -2) {
      unsigned long invalID = meshPoints.size();
      localFaces[i].nCornerPoints = 4;
      localFaces[i].cornerPoints[0] = invalID;
      localFaces[i].cornerPoints[1] = invalID;
      localFaces[i].cornerPoints[2] = invalID;
      localFaces[i].cornerPoints[3] = invalID;
      --nFacesLoc;
    }
  }

  sort(localFaces.begin(), localFaces.end());
  localFaces.resize(nFacesLoc);

  /*---------------------------------------------------------------------------*/
  /*--- Step 2: Preparation of the localFaces vector, such that the info    ---*/
  /*---         stored in this vector can be separated in a contribution    ---*/
  /*---         from the internal faces and a contribution from the faces   ---*/
  /*---         that belong to physical boundaries.                         ---*/
  /*---------------------------------------------------------------------------*/

  /* Sort localFaces again, but now such that the boundary faces are numbered
     first, followed by the matching faces and at the end of localFaces the
     non-matching faces are stored. In order to carry out this sorting the
     functor CSortFaces is used for comparison. Within the categories
     the sorting depends on the time levels of the adjacent volume elements
     of the face. */
  sort(localFaces.begin(), localFaces.end(), CSortFaces(nVolElemOwned, nVolElemTot, volElem.data()));

  /*--- Carry out a possible swap of side 0 and side 1 of the faces. This is
        done for the following reasons (in order of importance).
        1: The element on side 0 must always be an owned element.
        2: The VTK type of the element on side 0 must be less than or equal
           to the VTK type of the element on side 1. This is done to reduce
           the number of standard elements for the matching faces.
        3: For the same VTK types the element on side 0 must be the smallest.

        Furthermore, make sure that for boundary faces and non-matching faces
        the corresponding element is always on side 0. ---*/
  for (unsigned long i = 0; i < localFaces.size(); ++i) {
    /* Determine whether or not the ajacent elements must be swapped. */
    bool swapElements;
    if (localFaces[i].elemID0 < nVolElemTot && localFaces[i].elemID1 < nVolElemTot) {
      /* This is an internal matching face. First check if one of the adjacent
         elements is not an owned element. */
      if (localFaces[i].elemID0 >= nVolElemOwned || localFaces[i].elemID1 >= nVolElemOwned) {
        /* One of the element is not owned. Make sure that the owned element
           is on side 0. */
        swapElements = localFaces[i].elemID0 > localFaces[i].elemID1;
      } else if (localFaces[i].elemType0 == localFaces[i].elemType1) {
        /* The same element type on both sides. Make sure that the element
           with the smallest ID is stored on side 0 of the face. */
        swapElements = localFaces[i].elemID0 > localFaces[i].elemID1;
      } else {
        /* Different element types. Make sure that the lowest element
           type will be stored on side 0 of the face. */
        swapElements = localFaces[i].elemType0 > localFaces[i].elemType1;
      }
    } else {
      /* Either a boundary face or a non-matching face. It must be swapped
         if the element is currently on side 1 of the face. */
      swapElements = localFaces[i].elemID1 < nVolElemTot;
    }

    /* Swap the adjacent elements of the face, if needed. Note that
       also the sequence of the corner points must be altered in order
       to obey the right hand rule. */
    if (swapElements) {
      swap(localFaces[i].elemID0, localFaces[i].elemID1);
      swap(localFaces[i].nPolyGrid0, localFaces[i].nPolyGrid1);
      swap(localFaces[i].nPolySol0, localFaces[i].nPolySol1);
      swap(localFaces[i].nDOFsElem0, localFaces[i].nDOFsElem1);
      swap(localFaces[i].elemType0, localFaces[i].elemType1);
      swap(localFaces[i].faceID0, localFaces[i].faceID1);

      if (localFaces[i].nCornerPoints == 2)
        swap(localFaces[i].cornerPoints[0], localFaces[i].cornerPoints[1]);
      else
        swap(localFaces[i].cornerPoints[0], localFaces[i].cornerPoints[2]);

      localFaces[i].elem0IsOwner = !localFaces[i].elem0IsOwner;
    }
  }

  /*--- For triangular faces with a pyramid as an adjacent element, it must be
        made sure that the first corner point does not coincide with the top of
        the pyramid. Otherwise it is impossible to carry the transformation to
        the standard pyramid element. ---*/
  for (unsigned long i = 0; i < localFaces.size(); ++i) {
    /* Check for a triangular face. */
    if (localFaces[i].nCornerPoints == 3) {
      /* Determine if the corner points correspond to a corner point of a
         pyramid. A pyramid can in principle occur on both sides of the face. */
      bool cornerIsTopPyramid[] = {false, false, false};

      if (localFaces[i].elemType0 == PYRAMID) {
        unsigned long topPyramid = volElem[localFaces[i].elemID0].nodeIDsGrid.back();
        if (localFaces[i].cornerPoints[0] == topPyramid) cornerIsTopPyramid[0] = true;
        if (localFaces[i].cornerPoints[1] == topPyramid) cornerIsTopPyramid[1] = true;
        if (localFaces[i].cornerPoints[2] == topPyramid) cornerIsTopPyramid[2] = true;
      }

      if (localFaces[i].elemType1 == PYRAMID) {
        unsigned long topPyramid = volElem[localFaces[i].elemID1].nodeIDsGrid.back();
        if (localFaces[i].cornerPoints[0] == topPyramid) cornerIsTopPyramid[0] = true;
        if (localFaces[i].cornerPoints[1] == topPyramid) cornerIsTopPyramid[1] = true;
        if (localFaces[i].cornerPoints[2] == topPyramid) cornerIsTopPyramid[2] = true;
      }

      /* Check if the first element corresponds to the top of a pyramid. */
      if (cornerIsTopPyramid[0]) {
        /* The sequence of the points of the face must be altered. It is done in
           such a way that the orientation remains the same. Store the original
           point numbering. */
        unsigned long tmp[] = {localFaces[i].cornerPoints[0], localFaces[i].cornerPoints[1],
                               localFaces[i].cornerPoints[2]};

        /* Determine the situation. */
        if (!cornerIsTopPyramid[1]) {
          /* cornerPoint[1] is not a top of a pyramid. Hence this will become
             the point 0 of the triangle. */
          localFaces[i].cornerPoints[0] = tmp[1];
          localFaces[i].cornerPoints[1] = tmp[2];
          localFaces[i].cornerPoints[2] = tmp[0];
        } else {
          /* Only cornerPoint[2] is not a top of a pyramid. Hence this will
             become point 0 of the triangle. */
          localFaces[i].cornerPoints[0] = tmp[2];
          localFaces[i].cornerPoints[1] = tmp[0];
          localFaces[i].cornerPoints[2] = tmp[1];
        }
      } else if (cornerIsTopPyramid[1] && !cornerIsTopPyramid[2]) {
        /* cornerPoint[1] is a top of a pyramid and the other two corners
           are not. Change the sequence, such that cornerPoint[2] will become
           point 0 of the triangle. */
        unsigned long tmp[] = {localFaces[i].cornerPoints[0], localFaces[i].cornerPoints[1],
                               localFaces[i].cornerPoints[2]};
        localFaces[i].cornerPoints[0] = tmp[2];
        localFaces[i].cornerPoints[1] = tmp[0];
        localFaces[i].cornerPoints[2] = tmp[1];
      }
    }
  }

  /*--- Determine the number of matching and non-matching internal faces.
        For the matching faces, determine these numbers per time level
        and also make a distinction between internal faces and faces that
        involve a halo element. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
  nMatchingFacesInternal.assign(nTimeLevels + 1, 0);
  nMatchingFacesWithHaloElem.assign(nTimeLevels + 1, 0);

  unsigned long nNonMatchingFaces = 0;
  for (unsigned long i = 0; i < localFaces.size(); ++i) {
    if (localFaces[i].faceIndicator == -1) {
      const unsigned long e0 = localFaces[i].elemID0;
      const unsigned long e1 = localFaces[i].elemID1;

      if (e1 < nVolElemTot) {
        const unsigned short timeLevel = min(volElem[e0].timeLevel, volElem[e1].timeLevel);
        if (e1 < nVolElemOwned)
          ++nMatchingFacesInternal[timeLevel + 1];
        else
          ++nMatchingFacesWithHaloElem[timeLevel + 1];
      } else
        ++nNonMatchingFaces;
    }
  }

  if (nNonMatchingFaces) {
    ostringstream message;
    message << nNonMatchingFaces << " non-matching internal faces found.\n"
            << "This is not supported yet." << endl;
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /* Put nMatchingFacesInternal and nMatchingFacesWithHaloElem in
     cumulative storage format. */
  for (unsigned short i = 0; i < nTimeLevels; ++i) nMatchingFacesInternal[i + 1] += nMatchingFacesInternal[i];

  nMatchingFacesWithHaloElem[0] = nMatchingFacesInternal[nTimeLevels];
  for (unsigned short i = 0; i < nTimeLevels; ++i) nMatchingFacesWithHaloElem[i + 1] += nMatchingFacesWithHaloElem[i];

  /*---------------------------------------------------------------------------*/
  /*--- Step 3: Create the local face based data structure for the internal ---*/
  /*---         faces. These are needed for the computation of the surface  ---*/
  /*---         integral in DG-FEM.                                         ---*/
  /*---------------------------------------------------------------------------*/

  /* Allocate the memory for the matching faces. */
  matchingFaces.resize(nMatchingFacesWithHaloElem[nTimeLevels]);

  /*--- Loop over the volume elements to determine the maximum number
        of DOFs for the volume elements. Allocate the memory for the
        vector used to store the DOFs of the element.  ---*/
  unsigned short nDOFsVolMax = 0;
  for (unsigned long i = 0; i < nVolElemTot; ++i) {
    nDOFsVolMax = max(nDOFsVolMax, volElem[i].nDOFsGrid);
    nDOFsVolMax = max(nDOFsVolMax, volElem[i].nDOFsSol);
  }

  vector<unsigned long> DOFsElem(nDOFsVolMax);

  /*--- Loop over localFaces to create the connectivity information, which
        is stored in matchingFaces. ---*/
  unsigned long ii = 0;
  for (unsigned long i = 0; i < localFaces.size(); ++i) {
    if (localFaces[i].faceIndicator == -1 && localFaces[i].elemID1 < nVolElemTot) {
      /*--- Determine the number of DOFs on both sides of the face for both the
            grid and solution. This value depends on the face type. ---*/
      unsigned short sizeDOFsGridFaceSide0 = 0, sizeDOFsGridFaceSide1 = 0;
      unsigned short sizeDOFsSolFaceSide0 = 0, sizeDOFsSolFaceSide1 = 0;

      unsigned short VTK_Type = 0;  // To avoid a compiler warning.
      switch (localFaces[i].nCornerPoints) {
        case 2:
          /* Face is a line. */
          VTK_Type = LINE;
          sizeDOFsGridFaceSide0 = localFaces[i].nPolyGrid0 + 1;
          sizeDOFsGridFaceSide1 = localFaces[i].nPolyGrid1 + 1;
          sizeDOFsSolFaceSide0 = localFaces[i].nPolySol0 + 1;
          sizeDOFsSolFaceSide1 = localFaces[i].nPolySol1 + 1;
          break;

        case 3:
          /* Face is a triangle. */
          VTK_Type = TRIANGLE;
          sizeDOFsGridFaceSide0 = (localFaces[i].nPolyGrid0 + 1) * (localFaces[i].nPolyGrid0 + 2) / 2;
          sizeDOFsGridFaceSide1 = (localFaces[i].nPolyGrid1 + 1) * (localFaces[i].nPolyGrid1 + 2) / 2;
          sizeDOFsSolFaceSide0 = (localFaces[i].nPolySol0 + 1) * (localFaces[i].nPolySol0 + 2) / 2;
          sizeDOFsSolFaceSide1 = (localFaces[i].nPolySol1 + 1) * (localFaces[i].nPolySol1 + 2) / 2;
          break;

        case 4:
          /* Face is a quadrilateral. */
          VTK_Type = QUADRILATERAL;
          sizeDOFsGridFaceSide0 = (localFaces[i].nPolyGrid0 + 1) * (localFaces[i].nPolyGrid0 + 1);
          sizeDOFsGridFaceSide1 = (localFaces[i].nPolyGrid1 + 1) * (localFaces[i].nPolyGrid1 + 1);
          sizeDOFsSolFaceSide0 = (localFaces[i].nPolySol0 + 1) * (localFaces[i].nPolySol0 + 1);
          sizeDOFsSolFaceSide1 = (localFaces[i].nPolySol1 + 1) * (localFaces[i].nPolySol1 + 1);
          break;
      }

      /* Abbreviate the adjacent element ID's and store them in matchingFaces. */
      const unsigned long v0 = localFaces[i].elemID0;
      const unsigned long v1 = localFaces[i].elemID1;

      matchingFaces[ii].elemID0 = v0;
      matchingFaces[ii].elemID1 = v1;

      /* Allocate the memory for the connectivities of the face. */
      matchingFaces[ii].DOFsGridFaceSide0.resize(sizeDOFsGridFaceSide0);
      matchingFaces[ii].DOFsGridFaceSide1.resize(sizeDOFsGridFaceSide1);
      matchingFaces[ii].DOFsSolFaceSide0.resize(sizeDOFsSolFaceSide0);
      matchingFaces[ii].DOFsSolFaceSide1.resize(sizeDOFsSolFaceSide1);

      matchingFaces[ii].DOFsGridElementSide0.resize(volElem[v0].nDOFsGrid);
      matchingFaces[ii].DOFsGridElementSide1.resize(volElem[v1].nDOFsGrid);
      matchingFaces[ii].DOFsSolElementSide0.resize(volElem[v0].nDOFsSol);
      matchingFaces[ii].DOFsSolElementSide1.resize(volElem[v1].nDOFsSol);

      /*--- Create the connectivities of the adjacent elements in the correct
            sequence as well as the connectivities of the face. The
            connectivities of both sides of the face are determined for the
            grid and solution. First for side 0.   ---*/
      bool swapFaceInElementSide0;
      for (unsigned short j = 0; j < volElem[v0].nDOFsSol; ++j) DOFsElem[j] = volElem[v0].offsetDOFsSolLocal + j;

      CreateConnectivitiesFace(VTK_Type, localFaces[i].cornerPoints, volElem[v0].VTK_Type, volElem[v0].nPolyGrid,
                               volElem[v0].nodeIDsGrid, volElem[v0].nPolySol, DOFsElem.data(), swapFaceInElementSide0,
                               matchingFaces[ii].DOFsSolFaceSide0.data(), matchingFaces[ii].DOFsSolElementSide0.data());

      for (unsigned short j = 0; j < volElem[v0].nDOFsGrid; ++j) DOFsElem[j] = volElem[v0].nodeIDsGrid[j];

      CreateConnectivitiesFace(VTK_Type, localFaces[i].cornerPoints, volElem[v0].VTK_Type, volElem[v0].nPolyGrid,
                               volElem[v0].nodeIDsGrid, volElem[v0].nPolyGrid, DOFsElem.data(), swapFaceInElementSide0,
                               matchingFaces[ii].DOFsGridFaceSide0.data(),
                               matchingFaces[ii].DOFsGridElementSide0.data());

      /*--- And also for side 1 of the face. ---*/
      bool swapFaceInElementSide1;
      for (unsigned short j = 0; j < volElem[v1].nDOFsSol; ++j) DOFsElem[j] = volElem[v1].offsetDOFsSolLocal + j;

      CreateConnectivitiesFace(VTK_Type, localFaces[i].cornerPoints, volElem[v1].VTK_Type, volElem[v1].nPolyGrid,
                               volElem[v1].nodeIDsGrid, volElem[v1].nPolySol, DOFsElem.data(), swapFaceInElementSide1,
                               matchingFaces[ii].DOFsSolFaceSide1.data(), matchingFaces[ii].DOFsSolElementSide1.data());

      for (unsigned short j = 0; j < volElem[v1].nDOFsGrid; ++j) DOFsElem[j] = volElem[v1].nodeIDsGrid[j];

      CreateConnectivitiesFace(VTK_Type, localFaces[i].cornerPoints, volElem[v1].VTK_Type, volElem[v1].nPolyGrid,
                               volElem[v1].nodeIDsGrid, volElem[v1].nPolyGrid, DOFsElem.data(), swapFaceInElementSide1,
                               matchingFaces[ii].DOFsGridFaceSide1.data(),
                               matchingFaces[ii].DOFsGridElementSide1.data());

      /*--- Search in the standard elements for faces for a matching
            standard element. If not found, create a new standard element.
            Note that both the grid and the solution representation must
            match with the standard element. ---*/
      unsigned long j;
      for (j = 0; j < standardMatchingFacesSol.size(); ++j) {
        if (standardMatchingFacesSol[j].SameStandardMatchingFace(
                VTK_Type, localFaces[i].JacFaceIsConsideredConstant, localFaces[i].elemType0, localFaces[i].nPolySol0,
                localFaces[i].elemType1, localFaces[i].nPolySol1, swapFaceInElementSide0, swapFaceInElementSide1) &&
            standardMatchingFacesGrid[j].SameStandardMatchingFace(
                VTK_Type, localFaces[i].JacFaceIsConsideredConstant, localFaces[i].elemType0, localFaces[i].nPolyGrid0,
                localFaces[i].elemType1, localFaces[i].nPolyGrid1, swapFaceInElementSide0, swapFaceInElementSide1)) {
          matchingFaces[ii].indStandardElement = j;
          break;
        }
      }

      /* Create the new standard elements if no match was found. */
      if (j == standardMatchingFacesSol.size()) {
        standardMatchingFacesSol.emplace_back(VTK_Type, localFaces[i].elemType0, localFaces[i].nPolySol0,
                                              localFaces[i].elemType1, localFaces[i].nPolySol1,
                                              localFaces[i].JacFaceIsConsideredConstant, swapFaceInElementSide0,
                                              swapFaceInElementSide1, config);

        standardMatchingFacesGrid.emplace_back(
            VTK_Type, localFaces[i].elemType0, localFaces[i].nPolyGrid0, localFaces[i].elemType1,
            localFaces[i].nPolyGrid1, localFaces[i].JacFaceIsConsideredConstant, swapFaceInElementSide0,
            swapFaceInElementSide1, config, standardMatchingFacesSol[j].GetOrderExact());
        matchingFaces[ii].indStandardElement = j;
      }

      /* Update the counter ii for the next internal matching face. */
      ++ii;
    }
  }

  /* Carry out an additional sorting of matchingFaces, such that faces with the
     same standard face are contiguous in memory as much as possible. This allows
     for a simultaneous treatment of faces, which increases the gemm performance.
     Note that this sorting takes place for each time level, because the grouping
     of faces per time level is essential for the correct functioning of the code.
     The same holds for the distinction between internal faces and faces that
     are adjacent to a halo element. */
  for (unsigned short i = 0; i < nTimeLevels; ++i) {
    sort(matchingFaces.begin() + nMatchingFacesInternal[i], matchingFaces.begin() + nMatchingFacesInternal[i + 1]);

    sort(matchingFaces.begin() + nMatchingFacesWithHaloElem[i],
         matchingFaces.begin() + nMatchingFacesWithHaloElem[i + 1]);
  }

  /* Loop over the faces to flag the elements that share one or more faces
     with elements from a lower time level. */
  vector<bool> elemAdjLowTimeLevel(nVolElemTot, false);
  for (unsigned long i = 0; i < matchingFaces.size(); ++i) {
    const unsigned long e0 = matchingFaces[i].elemID0;
    const unsigned long e1 = matchingFaces[i].elemID1;

    if (volElem[e0].timeLevel > volElem[e1].timeLevel) elemAdjLowTimeLevel[e0] = true;
    if (volElem[e1].timeLevel > volElem[e0].timeLevel) elemAdjLowTimeLevel[e1] = true;
  }

  /* Determine the list of elements per time level, which share one or more
     faces with elements of the lower time level. Make a distinction between
     owned and halo elements. */
  ownedElemAdjLowTimeLevel.resize(nTimeLevels);
  haloElemAdjLowTimeLevel.resize(nTimeLevels);

  for (unsigned long i = 0; i < nVolElemOwned; ++i) {
    if (elemAdjLowTimeLevel[i]) ownedElemAdjLowTimeLevel[volElem[i].timeLevel].push_back(i);
  }

  for (unsigned long i = nVolElemOwned; i < nVolElemTot; ++i) {
    if (elemAdjLowTimeLevel[i]) haloElemAdjLowTimeLevel[volElem[i].timeLevel].push_back(i);
  }

  /* Determine the offsets of the solution DOFs in the solution working vectors
     of the solver for the different time levels. These are needed when time
     accurate local time stepping is employed. */
  vector<unsigned long> counterDOFs(nTimeLevels, 0);

  for (unsigned long i = 0; i < nVolElemTot; ++i) {
    const unsigned short lev = volElem[i].timeLevel;
    volElem[i].offsetDOFsSolThisTimeLevel = counterDOFs[lev];
    counterDOFs[lev] += volElem[i].nDOFsSol;

    if (elemAdjLowTimeLevel[i]) {
      volElem[i].offsetDOFsSolPrevTimeLevel = counterDOFs[lev - 1];
      counterDOFs[lev - 1] += volElem[i].nDOFsSol;
    } else {
      volElem[i].offsetDOFsSolPrevTimeLevel = ULONG_MAX;
    }
  }

  /*---------------------------------------------------------------------------*/
  /*--- Step 4: Create the local face based data structure for the faces    ---*/
  /*---         that belong to the physical boundaries.                     ---*/
  /*---------------------------------------------------------------------------*/

  /* Loop over the boundary markers. ---*/
  unsigned long indBegMarker = 0;
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    /* Initialize the number of surface elements per time level to zero. */
    boundaries[iMarker].nSurfElem.assign(nTimeLevels + 1, 0);

    /*--- The periodic boundaries are skipped, because these are not physical
          boundaries and are treated via the halo elements. These have already
          been created. ---*/
    if (!boundaries[iMarker].periodicBoundary) {
      /* Determine the end index for this marker in localFaces. Note that
         localFaces is sorted such that the boundary faces are first and
         also grouped per boundary, see the functor CSortFaces. */
      unsigned long indEndMarker = indBegMarker;
      for (; indEndMarker < localFaces.size(); ++indEndMarker) {
        if (localFaces[indEndMarker].faceIndicator != (short)iMarker) break;
      }

      /* Easier storage of the surface elements for this boundary. */
      vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;

      /*--- Loop over range in localFaces for this boundary marker and create
            the connectivity information, which is stored in surfElem. ---*/
      for (unsigned long i = indBegMarker; i < indEndMarker; ++i) {
        /* Determine the time level of the adjacent element and increment
           the number of surface elements for this time level. The +1 is there,
           because this vector will be put in cumulative storage format
           afterwards. */
        ii = i - indBegMarker;
        const unsigned short timeLevel = volElem[surfElem[ii].volElemID].timeLevel;
        ++boundaries[iMarker].nSurfElem[timeLevel + 1];

        /*--- Determine the number of DOFs of the face for both the grid
              and solution. This value depends on the face type. ---*/
        unsigned short sizeDOFsGridFace = 0, sizeDOFsSolFace = 0;

        unsigned short VTK_Type = 0;  // To avoid a compiler warning.
        switch (localFaces[i].nCornerPoints) {
          case 2:
            /* Face is a line. */
            VTK_Type = LINE;
            sizeDOFsGridFace = localFaces[i].nPolyGrid0 + 1;
            sizeDOFsSolFace = localFaces[i].nPolySol0 + 1;
            break;

          case 3:
            /* Face is a triangle. */
            VTK_Type = TRIANGLE;
            sizeDOFsGridFace = (localFaces[i].nPolyGrid0 + 1) * (localFaces[i].nPolyGrid0 + 2) / 2;
            sizeDOFsSolFace = (localFaces[i].nPolySol0 + 1) * (localFaces[i].nPolySol0 + 2) / 2;
            break;

          case 4:
            /* Face is a quadrilateral. */
            VTK_Type = QUADRILATERAL;
            sizeDOFsGridFace = (localFaces[i].nPolyGrid0 + 1) * (localFaces[i].nPolyGrid0 + 1);
            sizeDOFsSolFace = (localFaces[i].nPolySol0 + 1) * (localFaces[i].nPolySol0 + 1);
            break;
        }

        /* Abbreviate the adjacent element ID. */
        const unsigned long v0 = localFaces[i].elemID0;

        /* Allocate the memory for the connectivities of the face. */
        surfElem[ii].DOFsGridFace.resize(sizeDOFsGridFace);
        surfElem[ii].DOFsSolFace.resize(sizeDOFsSolFace);

        surfElem[ii].DOFsGridElement.resize(volElem[v0].nDOFsGrid);
        surfElem[ii].DOFsSolElement.resize(volElem[v0].nDOFsSol);

        /*--- Create the connectivities of the adjacent element in the correct
              sequence as well as the connectivities of the face. The
              connectivities are determined for the grid and solution. ---*/
        bool swapFaceInElement;
        for (unsigned short j = 0; j < volElem[v0].nDOFsSol; ++j) DOFsElem[j] = volElem[v0].offsetDOFsSolLocal + j;

        CreateConnectivitiesFace(VTK_Type, localFaces[i].cornerPoints, volElem[v0].VTK_Type, volElem[v0].nPolyGrid,
                                 volElem[v0].nodeIDsGrid, volElem[v0].nPolySol, DOFsElem.data(), swapFaceInElement,
                                 surfElem[ii].DOFsSolFace.data(), surfElem[ii].DOFsSolElement.data());

        for (unsigned short j = 0; j < volElem[v0].nDOFsGrid; ++j) DOFsElem[j] = volElem[v0].nodeIDsGrid[j];

        CreateConnectivitiesFace(VTK_Type, localFaces[i].cornerPoints, volElem[v0].VTK_Type, volElem[v0].nPolyGrid,
                                 volElem[v0].nodeIDsGrid, volElem[v0].nPolyGrid, DOFsElem.data(), swapFaceInElement,
                                 surfElem[ii].DOFsGridFace.data(), surfElem[ii].DOFsGridElement.data());

        /*--- Search in the standard elements for boundary faces for a matching
              standard element. If not found, create a new standard element.
              Note that both the grid and the solution representation must
              match with the standard element. ---*/
        unsigned long j;
        for (j = 0; j < standardBoundaryFacesSol.size(); ++j) {
          if (standardBoundaryFacesSol[j].SameStandardBoundaryFace(VTK_Type, localFaces[i].JacFaceIsConsideredConstant,
                                                                   localFaces[i].elemType0, localFaces[i].nPolySol0,
                                                                   swapFaceInElement) &&
              standardBoundaryFacesGrid[j].SameStandardBoundaryFace(VTK_Type, localFaces[i].JacFaceIsConsideredConstant,
                                                                    localFaces[i].elemType0, localFaces[i].nPolyGrid0,
                                                                    swapFaceInElement)) {
            surfElem[ii].indStandardElement = j;
            break;
          }
        }

        /* Create the new standard elements if no match was found. */
        if (j == standardBoundaryFacesSol.size()) {
          standardBoundaryFacesSol.emplace_back(VTK_Type, localFaces[i].elemType0, localFaces[i].nPolySol0,
                                                localFaces[i].JacFaceIsConsideredConstant, swapFaceInElement, config);

          standardBoundaryFacesGrid.emplace_back(VTK_Type, localFaces[i].elemType0, localFaces[i].nPolyGrid0,
                                                 localFaces[i].JacFaceIsConsideredConstant, swapFaceInElement, config,
                                                 standardBoundaryFacesSol[j].GetOrderExact());
          surfElem[ii].indStandardElement = j;
        }
      }

      /* Set indBegMarker to indEndMarker for the next marker. */
      indBegMarker = indEndMarker;

      /* Put boundaries[iMarker].nSurfElem in cumulative storage. */
      for (unsigned short i = 0; i < nTimeLevels; ++i)
        boundaries[iMarker].nSurfElem[i + 1] += boundaries[iMarker].nSurfElem[i];

      /* Carry out another sort of the surface elements, such that surface elements
         belonging to the same standard element are contiguous in memory. In this
         way a simultaneous treatment of several surface elements is possible, which
         increases the performance. Note that this renumbering must happen within
         each time level. */
      for (unsigned short i = 0; i < nTimeLevels; ++i)
        sort(surfElem.begin() + boundaries[iMarker].nSurfElem[i],
             surfElem.begin() + boundaries[iMarker].nSurfElem[i + 1], CSortBoundaryFaces());
    }
  }
}

void CMeshFEM_DG::CreateStandardVolumeElements(CConfig* config) {
  /*--- Loop over the volume elements and create new standard elements if needed.
        Note that a new standard element is created when either the solution
        element or the grid element does not match. Note further that for the
        standard element of the grid the  order for the integration of the
        solution is used, such that the metric terms are computed in the correct
        integration points in case the polynomial order of the solution differs
        from that of the grid. ---*/
  for (unsigned long i = 0; i < nVolElemTot; ++i) {
    /* Check the existing standard elements in the list. */
    unsigned long j;
    for (j = 0; j < standardElementsSol.size(); ++j) {
      if (standardElementsSol[j].SameStandardElement(volElem[i].VTK_Type, volElem[i].nPolySol,
                                                     volElem[i].JacIsConsideredConstant) &&
          standardElementsGrid[j].SameStandardElement(volElem[i].VTK_Type, volElem[i].nPolyGrid,
                                                      volElem[i].JacIsConsideredConstant)) {
        volElem[i].indStandardElement = j;
        break;
      }
    }

    /* Create the new standard elements if no match was found. */
    if (j == standardElementsSol.size()) {
      standardElementsSol.emplace_back(volElem[i].VTK_Type, volElem[i].nPolySol, volElem[i].JacIsConsideredConstant,
                                       config);

      standardElementsGrid.emplace_back(volElem[i].VTK_Type, volElem[i].nPolyGrid, volElem[i].JacIsConsideredConstant,
                                        config, standardElementsSol[j].GetOrderExact(),
                                        standardElementsSol[j].GetRDOFs(), standardElementsSol[j].GetSDOFs(),
                                        standardElementsSol[j].GetTDOFs());
      volElem[i].indStandardElement = j;
    }
  }
}

void CMeshFEM_DG::SetSendReceive(const CConfig* config) {
  /*----------------------------------------------------------------------------*/
  /*--- Step 1: Determine the ranks from which this rank has to receive data ---*/
  /*---         during the actual communication of halo data, as well as the ---*/
  /*            data that must be communicated.                              ---*/
  /*----------------------------------------------------------------------------*/

  /* Determine for every element the local offset of the solution DOFs.
     The if statement in the line below is to avoid problems when the
     partition is empty. */
  if (nVolElemTot) volElem[0].offsetDOFsSolLocal = 0;
  for (unsigned long i = 1; i < nVolElemTot; ++i)
    volElem[i].offsetDOFsSolLocal = volElem[i - 1].offsetDOFsSolLocal + volElem[i - 1].nDOFsSol;

  /* Determine the ranks from which this rank will receive halo data. */
  vector<int> recvFromRank(size, 0);
  for (unsigned long i = nVolElemOwned; i < nVolElemTot; ++i) recvFromRank[volElem[i].rankOriginal] = 1;

  map<int, int> rankToIndRecvBuf;
  for (int i = 0; i < size; ++i) {
    if (recvFromRank[i]) {
      int ind = (int)rankToIndRecvBuf.size();
      rankToIndRecvBuf[i] = ind;
    }
  }

  ranksRecv.resize(rankToIndRecvBuf.size());
  map<int, int>::const_iterator MI = rankToIndRecvBuf.begin();
  for (unsigned long i = 0; i < rankToIndRecvBuf.size(); ++i, ++MI) ranksRecv[i] = MI->first;

  /* Define and determine the buffers to send the global indices of my halo
     elements to the appropriate ranks and the vectors which store the
     elements that I will receive from these ranks. */
  vector<vector<unsigned long> > longBuf(rankToIndRecvBuf.size(), vector<unsigned long>(0));
  entitiesRecv.resize(rankToIndRecvBuf.size());

  for (unsigned long i = nVolElemOwned; i < nVolElemTot; ++i) {
    MI = rankToIndRecvBuf.find(volElem[i].rankOriginal);
    longBuf[MI->second].push_back(volElem[i].elemIDGlobal);

    entitiesRecv[MI->second].push_back(i);
  }

  /* Determine the mapping from global element ID to local owned element ID. */
  map<unsigned long, unsigned long> globalElemIDToLocalInd;
  for (unsigned long i = 0; i < nVolElemOwned; ++i) globalElemIDToLocalInd[volElem[i].elemIDGlobal] = i;

#ifdef HAVE_MPI

  /*--- Parallel mode. First determine the number of ranks to which this
        rank has to send halo data during the actual exchange. ---*/
  int nRankSend;
  vector<int> sizeReduce(size, 1);

  SU2_MPI::Reduce_scatter(recvFromRank.data(), &nRankSend, sizeReduce.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /* Resize ranksSend and the first index of entitiesSend to the number of
     ranks to which this rank has to send data. */
  ranksSend.resize(nRankSend);
  entitiesSend.resize(nRankSend);

  /*--- Send all the data using non-blocking sends. ---*/
  vector<SU2_MPI::Request> commReqs(ranksRecv.size());

  for (unsigned long i = 0; i < ranksRecv.size(); ++i) {
    int dest = ranksRecv[i];
    SU2_MPI::Isend(longBuf[i].data(), longBuf[i].size(), MPI_UNSIGNED_LONG, dest, dest, SU2_MPI::GetComm(),
                   &commReqs[i]);
  }

  /* Loop over the number of ranks from which I receive data about the
     global element ID's that I must send. */
  for (int i = 0; i < nRankSend; ++i) {
    /* Block until a message arrivesi and determine the source. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    ranksSend[i] = status.MPI_SOURCE;

    /* Determine the size of the message, allocate the memory for the
       receive buffer and receive the message. */
    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    entitiesSend[i].resize(sizeMess);
    SU2_MPI::Recv(entitiesSend[i].data(), sizeMess, MPI_UNSIGNED_LONG, ranksSend[i], rank, SU2_MPI::GetComm(), &status);

    /* Convert the global indices currently stored in entitiesSend[i]
       to local indices. */
    for (int j = 0; j < sizeMess; ++j) {
      map<unsigned long, unsigned long>::const_iterator LMI;
      LMI = globalElemIDToLocalInd.find(entitiesSend[i][j]);

      if (LMI == globalElemIDToLocalInd.end()) SU2_MPI::Error("This should not happen", CURRENT_FUNCTION);

      entitiesSend[i][j] = LMI->second;
    }
  }

  /* Complete the non-blocking sends and synchronize the rank, because
     wild cards have been used. */
  SU2_MPI::Waitall(ranksRecv.size(), commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else
  /* Sequential mode. Resize ranksSend and the first index of entitiesSend to
     the number of ranks to which this rank has to send data. This number is
     only non-zero when periodic boundaries are present in the grid. */
  ranksSend.resize(ranksRecv.size());
  entitiesSend.resize(ranksRecv.size());

  /* Convert the global element ID's of longBuf to local indices, which are
     stored in entitiesSend[0]. Note that an additional test for longBuf.size()
     is necessary to avoid problems. */
  if (longBuf.size()) {
    entitiesSend[0].resize(longBuf[0].size());

    for (unsigned long i = 0; i < longBuf[0].size(); ++i) {
      map<unsigned long, unsigned long>::const_iterator LMI;
      LMI = globalElemIDToLocalInd.find(longBuf[0][i]);

      if (LMI == globalElemIDToLocalInd.end()) SU2_MPI::Error("This should not happen", CURRENT_FUNCTION);

      entitiesSend[0][i] = LMI->second;
    }
  }

#endif

  /*----------------------------------------------------------------------------*/
  /*--- Step 2: Determine the rotational periodic transformations as well as ---*/
  /*---         the halo elements for which these must be applied.           ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Loop over the markers and determine the mapping for the rotationally
        periodic transformations. The mapping is from the marker to the first
        index in the vectors to store the rotationally periodic halo elements. ---*/
  map<short, unsigned short> mapRotationalPeriodicToInd;

  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
      if (fabs(angles[0]) > 1.e-5 || fabs(angles[1]) > 1.e-5 || fabs(angles[2]) > 1.e-5) {
        unsigned short curSize = mapRotationalPeriodicToInd.size();
        mapRotationalPeriodicToInd[iMarker] = curSize;
      }
    }
  }

  /* Store the rotationally periodic indices in rotPerMarkers. */
  rotPerMarkers.reserve(mapRotationalPeriodicToInd.size());
  for (auto SMI = mapRotationalPeriodicToInd.begin(); SMI != mapRotationalPeriodicToInd.end(); ++SMI)
    rotPerMarkers.push_back(SMI->first);

  /* Resize the first index of rotPerHalos to the correct size. */
  rotPerHalos.resize(mapRotationalPeriodicToInd.size());

  /*--- Loop over the volume elements and store the indices of the rotationally
        periodic halo elements in rotPerHalos.     ---*/
  for (unsigned long i = nVolElemOwned; i < nVolElemTot; ++i) {
    if (volElem[i].periodIndexToDonor > -1) {
      map<short, unsigned short>::const_iterator SMI;
      SMI = mapRotationalPeriodicToInd.find(volElem[i].periodIndexToDonor);

      if (SMI != mapRotationalPeriodicToInd.end()) rotPerHalos[SMI->second].push_back(i);
    }
  }
}

void CMeshFEM_DG::CreateConnectivitiesFace(const unsigned short VTK_TypeFace, const unsigned long* cornerPointsFace,
                                           const unsigned short VTK_TypeElem, const unsigned short nPolyGrid,
                                           const vector<unsigned long>& elemNodeIDsGrid, const unsigned short nPolyConn,
                                           const unsigned long* connElem, bool& swapFaceInElement,
                                           unsigned long* modConnFace, unsigned long* modConnElem) {
  /*--- Set swapFaceInElement to false. This variable is only relevant for
        triangular faces of a pyramid and quadrilateral faces of a prism.
        Only for these situations this variable will be passed to the
        appropriate functions. ---*/
  swapFaceInElement = false;

  /*--- Make a distinction between the types of the volume element and call
        the appropriate function to do the actual job. ---*/
  switch (VTK_TypeElem) {
    case TRIANGLE:
      CreateConnectivitiesLineAdjacentTriangle(cornerPointsFace, nPolyGrid, elemNodeIDsGrid, nPolyConn, connElem,
                                               modConnFace, modConnElem);
      break;

    case QUADRILATERAL:
      CreateConnectivitiesLineAdjacentQuadrilateral(cornerPointsFace, nPolyGrid, elemNodeIDsGrid, nPolyConn, connElem,
                                                    modConnFace, modConnElem);
      break;

    case TETRAHEDRON:
      CreateConnectivitiesTriangleAdjacentTetrahedron(cornerPointsFace, nPolyGrid, elemNodeIDsGrid, nPolyConn, connElem,
                                                      modConnFace, modConnElem);
      break;

    case PYRAMID: {
      switch (VTK_TypeFace) {
        case TRIANGLE:
          CreateConnectivitiesTriangleAdjacentPyramid(cornerPointsFace, nPolyGrid, elemNodeIDsGrid, nPolyConn, connElem,
                                                      swapFaceInElement, modConnFace, modConnElem);
          break;

        case QUADRILATERAL:
          CreateConnectivitiesQuadrilateralAdjacentPyramid(cornerPointsFace, nPolyGrid, elemNodeIDsGrid, nPolyConn,
                                                           connElem, modConnFace, modConnElem);
          break;
      }

      break;
    }

    case PRISM: {
      switch (VTK_TypeFace) {
        case TRIANGLE:
          CreateConnectivitiesTriangleAdjacentPrism(cornerPointsFace, nPolyGrid, elemNodeIDsGrid, nPolyConn, connElem,
                                                    modConnFace, modConnElem);
          break;

        case QUADRILATERAL:
          CreateConnectivitiesQuadrilateralAdjacentPrism(cornerPointsFace, nPolyGrid, elemNodeIDsGrid, nPolyConn,
                                                         connElem, swapFaceInElement, modConnFace, modConnElem);
          break;
      }

      break;
    }

    case HEXAHEDRON:
      CreateConnectivitiesQuadrilateralAdjacentHexahedron(cornerPointsFace, nPolyGrid, elemNodeIDsGrid, nPolyConn,
                                                          connElem, modConnFace, modConnElem);
      break;
  }
}

void CMeshFEM_DG::CreateConnectivitiesLineAdjacentQuadrilateral(
    const unsigned long* cornerPointsLine, const unsigned short nPolyGrid, const vector<unsigned long>& quadNodeIDsGrid,
    const unsigned short nPolyConn, const unsigned long* connQuad, unsigned long* modConnLine,
    unsigned long* modConnQuad) {
  /* Determine the indices of the four corner points of the quadrilateral. */
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPolyGrid;
  const unsigned short ind2 = (nPolyGrid + 1) * (nPolyGrid + 1) - 1;
  const unsigned short ind3 = ind2 - nPolyGrid;

  /* Easier storage of the two corner points of the line in the new numbering. */
  const unsigned long vert0 = cornerPointsLine[0];
  const unsigned long vert1 = cornerPointsLine[1];

  /*--- There exists a linear mapping from the indices of the numbering used in
        quadNodeIDsGrid to the indices of the target numbering. This mapping is of
        the form ii = a + b*i + c*j and jj = d + e*i + f*j, where ii,jj are the
        indices of the new numbering and i,j the indices of the numbering used
        in quadNodeIDsGrid (and connQuad). The values of the coefficients a,b,c,d,e,f
        depend on how the corner points of the line coincide with the corner points
        of the quad. This is determined below. The bool verticesDontMatch is
        there to check if vertices do not match. This should not happen, but
        it is checked for security. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
  bool verticesDontMatch = false;

  if (vert0 == quadNodeIDsGrid[ind0]) {
    /* Vert0 coincides with vertex 0 of the quad connectivity.
       Determine the situation for vert1. */
    if (vert1 == quadNodeIDsGrid[ind1]) {
      /* The new numbering is the same as the original numbering. */
      a = d = 0;
      b = f = 1;
      c = e = 0;
    } else if (vert1 == quadNodeIDsGrid[ind3]) {
      /* The i and j numbering are swapped. This is a left handed transformation. */
      a = d = 0;
      b = f = 0;
      c = e = 1;
    } else {
      verticesDontMatch = true;  // Vert1 does not match with a neigbor.
    }
  } else if (vert0 == quadNodeIDsGrid[ind1]) {
    /* Vert0 coincides with vertex 1 of the quad connectivity.
       Determine the situation for vert1. */
    if (vert1 == quadNodeIDsGrid[ind2]) {
      /* The i-direction of the new numbering corresponds to the j-direction
         of the original numbering, while the new j-direction is the negative
         i-direction of the original numbering. */
      a = 0;
      d = nPolyConn;
      b = f = 0;
      c = 1;
      e = -1;
    } else if (vert1 == quadNodeIDsGrid[ind0]) {
      /* The i-direction is negated, while the j-direction coincides.
         This is a left handed transformation. */
      a = nPolyConn;
      d = 0;
      b = -1;
      f = 1;
      c = e = 0;
    } else {
      verticesDontMatch = true;  // Vert1 does not match with a neigbor.
    }
  } else if (vert0 == quadNodeIDsGrid[ind2]) {
    /* Vert0 coincides with vertex 2 of the quad connectivity.
       Determine the situation for vert1. */
    if (vert1 == quadNodeIDsGrid[ind3]) {
      /* Both the i- and j-direction are negated. */
      a = d = nPolyConn;
      b = f = -1;
      c = e = 0;
    } else if (vert1 == quadNodeIDsGrid[ind1]) {
      /* The new i-direction is the original negative j-direction, while the
         new j-direction is the original negative i-direction. This is a
         left handed transformation. */
      a = d = nPolyConn;
      b = f = 0;
      c = e = -1;
    } else {
      verticesDontMatch = true;  // Vert1 does not match with a neigbor.
    }
  } else if (vert0 == quadNodeIDsGrid[ind3]) {
    /* Vert0 coincides with vertex 3 of the quad connectivity.
       Determine the situation for vert1. */
    if (vert1 == quadNodeIDsGrid[ind0]) {
      /* The new i-direction is the original negative j-direction, while the
         new j-direction is the original i-direction. */
      a = nPolyConn;
      d = 0;
      b = f = 0;
      c = -1;
      e = 1;
    } else if (vert1 == quadNodeIDsGrid[ind2]) {
      /* The i-directions coincide while the j-direction is negated.
         This is a left handed transformation. */
      a = 0;
      d = nPolyConn;
      b = 1;
      f = -1;
      c = e = 0;
    } else {
      verticesDontMatch = true;  // Vert1 does not match with a neigbor.
    }
  } else {
    /* Vert0 does not match with any of the corner vertices of the quad. */
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Loop over the DOFs of the original quad to create the connectivity
        of the quad that corresponds to the new numbering. ---*/
  unsigned short ind = 0;
  for (unsigned short j = 0; j <= nPolyConn; ++j) {
    for (unsigned short i = 0; i <= nPolyConn; ++i, ++ind) {
      /*--- Determine the ii and jj indices of the new numbering, convert it to
            a 1D index and shore the modified index in modConnQuad. ---*/
      unsigned short ii = a + i * b + j * c;
      unsigned short jj = d + i * e + j * f;
      unsigned short iind = jj * (nPolyConn + 1) + ii;

      modConnQuad[iind] = connQuad[ind];
    }
  }

  /*--- The line corresponds to face 0 of the quadrilateral. Hence the first
        nPolyConn+1 entries in modConnQuad are the DOFs of the line. Copy
        these entries from modConnQuad. ---*/
  for (unsigned short i = 0; i <= nPolyConn; ++i) modConnLine[i] = modConnQuad[i];
}

void CMeshFEM_DG::CreateConnectivitiesLineAdjacentTriangle(const unsigned long* cornerPointsLine,
                                                           const unsigned short nPolyGrid,
                                                           const vector<unsigned long>& triaNodeIDsGrid,
                                                           const unsigned short nPolyConn,
                                                           const unsigned long* connTria, unsigned long* modConnLine,
                                                           unsigned long* modConnTria) {
  /* Determine the indices of the 3 corner vertices of the triangle. */
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPolyGrid;
  const unsigned short ind2 = (nPolyGrid + 1) * (nPolyGrid + 2) / 2 - 1;

  /* Easier storage of the two corner points of the line in the new numbering. */
  const unsigned long vert0 = cornerPointsLine[0];
  const unsigned long vert1 = cornerPointsLine[1];

  /*--- There exists a linear mapping from the indices of the numbering used in
        triaNodeIDsGrid to the indices of the target numbering. This mapping is of
        the form ii = a + b*i + c*j and jj = d + e*i + f*j, where ii,jj are the
        indices of the new numbering and i,j the indices of the numbering used
        in triaNodeIDsGrid (and connTria). The values of the coefficients a,b,c,d,e,f
        depend on how the corner points of the line coincide with the corner points
        of the triangle. This is determined below. The bool verticesDontMatch is
        there to check if vertices do not match. This should not happen, but
        it is checked for security. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
  bool verticesDontMatch = false;

  if (vert0 == triaNodeIDsGrid[ind0]) {
    /* Vert0 coincides with vertex 0 of the triangle connectivity.
       Determine the situation for vert1. */
    if (vert1 == triaNodeIDsGrid[ind1]) {
      /* The new numbering is the same as the original numbering. */
      a = 0;
      b = 1;
      c = 0;
      d = 0;
      e = 0;
      f = 1;
    } else if (vert1 == triaNodeIDsGrid[ind2]) {
      /* The i and j numbering are swapped. This is a left handed transformation. */
      a = 0;
      b = 0;
      c = 1;
      d = 0;
      e = 1;
      f = 0;
    } else {
      verticesDontMatch = true;  // Vert1 does not match with a neigbor.
    }
  } else if (vert0 == triaNodeIDsGrid[ind1]) {
    /* Vert0 coincides with vertex 1 of the triangle connectivity.
       Determine the situation for vert1. */
    if (vert1 == triaNodeIDsGrid[ind2]) {
      /* The i-direction of the new numbering corresponds to the j-direction
         of the original numbering, while the new j-direction corresponds to
         a combination of the original i- and j-direction. */
      a = 0;
      b = 0;
      c = 1;
      d = nPolyConn;
      e = -1;
      f = -1;
    } else if (vert1 == triaNodeIDsGrid[ind0]) {
      /* The i-direction of the new numbering corresponds to a combination of
         the original i- and j-direction, while the new j-direction corresponds
         to the j-direction of the original numbering. This is a left
         handed transformation. */
      a = nPolyConn;
      b = -1;
      c = -1;
      d = 0;
      e = 0;
      f = 1;
    } else {
      verticesDontMatch = true;  // Vert1 does not match with a neigbor.
    }
  } else if (vert0 == triaNodeIDsGrid[ind2]) {
    /* Vert0 coincides with vertex 2 of the triangle connectivity.
       Determine the situation for vert1. */
    if (vert1 == triaNodeIDsGrid[ind0]) {
      /* The i-direction of the new numbering corresponds to a combination of
         the original i- and j-direction, while the new j-direction corresponds
         to the i-direction of the original numbering. */
      a = nPolyConn;
      b = -1;
      c = -1;
      d = 0;
      e = 1;
      f = 0;
    } else if (vert1 == triaNodeIDsGrid[ind1]) {
      /* The i-direction of the new numbering corresponds to the i-direction of
         the original numbering, while the new j-direction corresponds to a
         combination of the original i- and j-direction. This is a left handed
         transformation. */
      a = 0;
      b = 1;
      c = 0;
      d = nPolyConn;
      e = -1;
      f = -1;
    } else {
      verticesDontMatch = true;  // Vert1 does not match with a neigbor.
    }
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Loop over the DOFs of the original triangle to create the connectivity
        of the triangle that corresponds to the new numbering. ---*/
  unsigned short ind = 0;
  for (unsigned short j = 0; j <= nPolyConn; ++j) {
    for (unsigned short i = 0; i <= (nPolyConn - j); ++i, ++ind) {
      /*--- Determine the ii and jj indices of the new numbering, convert it to
            a 1D index and shore the modified index in modConnTria. ---*/
      unsigned short ii = a + i * b + j * c;
      unsigned short jj = d + i * e + j * f;

      unsigned short iind = jj * (nPolyConn + 1) + ii - jj * (jj - 1) / 2;

      modConnTria[iind] = connTria[ind];
    }
  }

  /*--- The line corresponds to face 0 of the triangle. Hence the first
        nPolyConn+1 entries in modConnTria are the DOFs of the line. Copy
        these entries from modConnTria. ---*/
  for (unsigned short i = 0; i <= nPolyConn; ++i) modConnLine[i] = modConnTria[i];
}

void CMeshFEM_DG::CreateConnectivitiesQuadrilateralAdjacentHexahedron(
    const unsigned long* cornerPointsQuad, const unsigned short nPolyGrid, const vector<unsigned long>& hexaNodeIDsGrid,
    const unsigned short nPolyConn, const unsigned long* connHexa, unsigned long* modConnQuad,
    unsigned long* modConnHexa) {
  /* Determine the indices of the eight corner points of the hexahedron. */
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPolyGrid;
  const unsigned short ind2 = (nPolyGrid + 1) * (nPolyGrid + 1) - 1;
  const unsigned short ind3 = ind2 - nPolyGrid;
  const unsigned short ind4 = (nPolyGrid + 1) * (nPolyGrid + 1) * nPolyGrid;
  const unsigned short ind5 = ind1 + ind4;
  const unsigned short ind6 = ind2 + ind4;
  const unsigned short ind7 = ind3 + ind4;

  /* Easier storage of the four corner points of the quad in the new numbering. */
  const unsigned long vert0 = cornerPointsQuad[0];
  const unsigned long vert1 = cornerPointsQuad[1];
  const unsigned long vert2 = cornerPointsQuad[2];
  const unsigned long vert3 = cornerPointsQuad[3];

  /*--- There exists a linear mapping from the indices of the numbering used in
        hexaNodeIDsGrid to the indices of the target numbering. This mapping is of
        the form ii = a + b*i + c*j + d*k, jj = e + f*i + g*j + h*k and
        kk = l + m*i + n*j + o*k, where ii,jj,kk are the indices of the new
        numbering and i,j,k the indices of the numbering used in in hexaNodeIDsGrid
        (and connHexa). The values of the coefficients a,b,c,d,e,f,g,h,l,m,n,o
        depend on how the corner points of the quad coincide with the corner points
        of the hexahedron. This is determined below. The bool verticesDontMatch is
        there to check if vertices do not match. This should not happen, but
        it is checked for security. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, l = 0, m = 0, n = 0, o = 0;
  bool verticesDontMatch = false;

  if (vert0 == hexaNodeIDsGrid[ind0] && vert1 == hexaNodeIDsGrid[ind1] &&         // ii = i.
      vert2 == hexaNodeIDsGrid[ind2] && vert3 == hexaNodeIDsGrid[ind3]) {         // jj = j.
    b = g = o = 1;                                                                // kk = k.
  } else if (vert0 == hexaNodeIDsGrid[ind0] && vert1 == hexaNodeIDsGrid[ind3] &&  // ii = j.
             vert2 == hexaNodeIDsGrid[ind2] && vert3 == hexaNodeIDsGrid[ind1]) {  // jj = i.
    c = f = o = 1;                                                                // kk = k.
  } else if (vert0 == hexaNodeIDsGrid[ind0] && vert1 == hexaNodeIDsGrid[ind1] &&  // ii = i.
             vert2 == hexaNodeIDsGrid[ind5] && vert3 == hexaNodeIDsGrid[ind4]) {  // jj = k.
    b = h = n = 1;                                                                // kk = j.
  } else if (vert0 == hexaNodeIDsGrid[ind0] && vert1 == hexaNodeIDsGrid[ind4] &&  // ii = k.
             vert2 == hexaNodeIDsGrid[ind5] && vert3 == hexaNodeIDsGrid[ind1]) {  // jj = i.
    d = f = n = 1;                                                                // kk = j.
  } else if (vert0 == hexaNodeIDsGrid[ind0] && vert1 == hexaNodeIDsGrid[ind3] &&  // ii = j.
             vert2 == hexaNodeIDsGrid[ind7] && vert3 == hexaNodeIDsGrid[ind4]) {  // jj = k.
    c = h = m = 1;                                                                // kk = i.
  } else if (vert0 == hexaNodeIDsGrid[ind0] && vert1 == hexaNodeIDsGrid[ind4] &&  // ii = k.
             vert2 == hexaNodeIDsGrid[ind7] && vert3 == hexaNodeIDsGrid[ind3]) {  // jj = j.
    d = g = m = 1;                                                                // kk = i.
  }

  else if (vert0 == hexaNodeIDsGrid[ind1] && vert1 == hexaNodeIDsGrid[ind0] &&  // ii = nPoly-i.
           vert2 == hexaNodeIDsGrid[ind3] && vert3 == hexaNodeIDsGrid[ind2]) {  // jj = j.
    a = nPolyConn;
    b = -1;
    g = o = 1;                                                                    // kk = k.
  } else if (vert0 == hexaNodeIDsGrid[ind1] && vert1 == hexaNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == hexaNodeIDsGrid[ind3] && vert3 == hexaNodeIDsGrid[ind0]) {  // jj = nPoly-i.
    e = nPolyConn;
    f = -1;
    c = o = 1;                                                                    // kk = k.
  } else if (vert0 == hexaNodeIDsGrid[ind1] && vert1 == hexaNodeIDsGrid[ind0] &&  // ii = nPoly-i.
             vert2 == hexaNodeIDsGrid[ind4] && vert3 == hexaNodeIDsGrid[ind5]) {  // jj = k.
    a = nPolyConn;
    b = -1;
    h = n = 1;                                                                    // kk = j.
  } else if (vert0 == hexaNodeIDsGrid[ind1] && vert1 == hexaNodeIDsGrid[ind5] &&  // ii = k.
             vert2 == hexaNodeIDsGrid[ind4] && vert3 == hexaNodeIDsGrid[ind0]) {  // jj = nPoly-i.
    e = nPolyConn;
    f = -1;
    d = n = 1;                                                                    // kk = j.
  } else if (vert0 == hexaNodeIDsGrid[ind1] && vert1 == hexaNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == hexaNodeIDsGrid[ind6] && vert3 == hexaNodeIDsGrid[ind5]) {  // jj = k.
    l = nPolyConn;
    m = -1;
    c = h = 1;                                                                    // kk = nPoly-i.
  } else if (vert0 == hexaNodeIDsGrid[ind1] && vert1 == hexaNodeIDsGrid[ind5] &&  // ii = k.
             vert2 == hexaNodeIDsGrid[ind6] && vert3 == hexaNodeIDsGrid[ind2]) {  // jj = j.
    l = nPolyConn;
    m = -1;
    d = g = 1;  // kk = nPoly-i.
  }

  else if (vert0 == hexaNodeIDsGrid[ind2] && vert1 == hexaNodeIDsGrid[ind1] &&  // ii = nPoly-j.
           vert2 == hexaNodeIDsGrid[ind0] && vert3 == hexaNodeIDsGrid[ind3]) {  // jj = nPoly-i.
    a = e = nPolyConn;
    c = f = -1;
    o = 1;                                                                        // kk = k.
  } else if (vert0 == hexaNodeIDsGrid[ind2] && vert1 == hexaNodeIDsGrid[ind3] &&  // ii = nPoly-i.
             vert2 == hexaNodeIDsGrid[ind0] && vert3 == hexaNodeIDsGrid[ind1]) {  // jj = nPoly-j.
    a = e = nPolyConn;
    b = g = -1;
    o = 1;                                                                        // kk = k.
  } else if (vert0 == hexaNodeIDsGrid[ind2] && vert1 == hexaNodeIDsGrid[ind1] &&  // ii = nPoly-j.
             vert2 == hexaNodeIDsGrid[ind5] && vert3 == hexaNodeIDsGrid[ind6]) {  // jj = k.
    a = l = nPolyConn;
    c = m = -1;
    h = 1;                                                                        // kk = nPoly-i.
  } else if (vert0 == hexaNodeIDsGrid[ind2] && vert1 == hexaNodeIDsGrid[ind6] &&  // ii = k.
             vert2 == hexaNodeIDsGrid[ind5] && vert3 == hexaNodeIDsGrid[ind1]) {  // jj = nPoly-j.
    e = l = nPolyConn;
    g = m = -1;
    d = 1;                                                                        // kk = nPoly-i.
  } else if (vert0 == hexaNodeIDsGrid[ind2] && vert1 == hexaNodeIDsGrid[ind3] &&  // ii = nPoly-i.
             vert2 == hexaNodeIDsGrid[ind7] && vert3 == hexaNodeIDsGrid[ind6]) {  // jj = k.
    a = l = nPolyConn;
    b = n = -1;
    h = 1;                                                                        // kk = nPoly-j.
  } else if (vert0 == hexaNodeIDsGrid[ind2] && vert1 == hexaNodeIDsGrid[ind6] &&  // ii = k.
             vert2 == hexaNodeIDsGrid[ind7] && vert3 == hexaNodeIDsGrid[ind3]) {  // jj = nPoly-i.
    e = l = nPolyConn;
    f = n = -1;
    d = 1;  // kk = nPoly-j.
  }

  else if (vert0 == hexaNodeIDsGrid[ind3] && vert1 == hexaNodeIDsGrid[ind0] &&  // ii = nPoly-j.
           vert2 == hexaNodeIDsGrid[ind1] && vert3 == hexaNodeIDsGrid[ind2]) {  // jj = i.
    a = nPolyConn;
    c = -1;
    f = o = 1;                                                                    // kk = k.
  } else if (vert0 == hexaNodeIDsGrid[ind3] && vert1 == hexaNodeIDsGrid[ind2] &&  // ii = i.
             vert2 == hexaNodeIDsGrid[ind1] && vert3 == hexaNodeIDsGrid[ind0]) {  // jj = nPoly-j.
    e = nPolyConn;
    g = -1;
    b = o = 1;                                                                    // kk = k.
  } else if (vert0 == hexaNodeIDsGrid[ind3] && vert1 == hexaNodeIDsGrid[ind0] &&  // ii = nPoly-j.
             vert2 == hexaNodeIDsGrid[ind4] && vert3 == hexaNodeIDsGrid[ind7]) {  // jj = k.
    a = nPolyConn;
    c = -1;
    h = m = 1;                                                                    // kk = i.
  } else if (vert0 == hexaNodeIDsGrid[ind3] && vert1 == hexaNodeIDsGrid[ind7] &&  // ii = k.
             vert2 == hexaNodeIDsGrid[ind4] && vert3 == hexaNodeIDsGrid[ind0]) {  // jj = nPoly-j.
    e = nPolyConn;
    g = -1;
    d = m = 1;                                                                    // kk = i.
  } else if (vert0 == hexaNodeIDsGrid[ind3] && vert1 == hexaNodeIDsGrid[ind2] &&  // ii = i.
             vert2 == hexaNodeIDsGrid[ind6] && vert3 == hexaNodeIDsGrid[ind7]) {  // jj = k.
    l = nPolyConn;
    n = -1;
    b = h = 1;                                                                    // kk = nPoly-j.
  } else if (vert0 == hexaNodeIDsGrid[ind3] && vert1 == hexaNodeIDsGrid[ind7] &&  // ii = k.
             vert2 == hexaNodeIDsGrid[ind6] && vert3 == hexaNodeIDsGrid[ind2]) {  // jj = i.
    l = nPolyConn;
    n = -1;
    d = f = 1;  // kk = nPoly-j.
  }

  else if (vert0 == hexaNodeIDsGrid[ind4] && vert1 == hexaNodeIDsGrid[ind5] &&  // ii = i.
           vert2 == hexaNodeIDsGrid[ind6] && vert3 == hexaNodeIDsGrid[ind7]) {  // jj = j.
    l = nPolyConn;
    o = -1;
    b = g = 1;                                                                    // kk = nPoly-k.
  } else if (vert0 == hexaNodeIDsGrid[ind4] && vert1 == hexaNodeIDsGrid[ind7] &&  // ii = j.
             vert2 == hexaNodeIDsGrid[ind6] && vert3 == hexaNodeIDsGrid[ind5]) {  // jj = i.
    l = nPolyConn;
    o = -1;
    c = f = 1;                                                                    // kk = nPoly-k.
  } else if (vert0 == hexaNodeIDsGrid[ind4] && vert1 == hexaNodeIDsGrid[ind5] &&  // ii = i.
             vert2 == hexaNodeIDsGrid[ind1] && vert3 == hexaNodeIDsGrid[ind0]) {  // jj = nPoly-k.
    e = nPolyConn;
    h = -1;
    b = n = 1;                                                                    // kk = j.
  } else if (vert0 == hexaNodeIDsGrid[ind4] && vert1 == hexaNodeIDsGrid[ind0] &&  // ii = nPoly-k.
             vert2 == hexaNodeIDsGrid[ind1] && vert3 == hexaNodeIDsGrid[ind5]) {  // jj = i.
    a = nPolyConn;
    d = -1;
    f = n = 1;                                                                    // kk = j.
  } else if (vert0 == hexaNodeIDsGrid[ind4] && vert1 == hexaNodeIDsGrid[ind7] &&  // ii = j.
             vert2 == hexaNodeIDsGrid[ind3] && vert3 == hexaNodeIDsGrid[ind0]) {  // jj = nPoly-k.
    e = nPolyConn;
    h = -1;
    c = m = 1;                                                                    // kk = i.
  } else if (vert0 == hexaNodeIDsGrid[ind4] && vert1 == hexaNodeIDsGrid[ind0] &&  // ii = nPoly-k.
             vert2 == hexaNodeIDsGrid[ind3] && vert3 == hexaNodeIDsGrid[ind7]) {  // jj = j.
    a = nPolyConn;
    d = -1;
    g = m = 1;  // kk = i.
  }

  else if (vert0 == hexaNodeIDsGrid[ind5] && vert1 == hexaNodeIDsGrid[ind6] &&  // ii = j.
           vert2 == hexaNodeIDsGrid[ind7] && vert3 == hexaNodeIDsGrid[ind4]) {  // jj = nPoly-i.
    e = l = nPolyConn;
    f = o = -1;
    c = 1;                                                                        // kk = nPoly-k.
  } else if (vert0 == hexaNodeIDsGrid[ind5] && vert1 == hexaNodeIDsGrid[ind4] &&  // ii = nPoly-i.
             vert2 == hexaNodeIDsGrid[ind7] && vert3 == hexaNodeIDsGrid[ind6]) {  // jj = j.
    a = l = nPolyConn;
    b = o = -1;
    g = 1;                                                                        // kk = nPoly-k.
  } else if (vert0 == hexaNodeIDsGrid[ind5] && vert1 == hexaNodeIDsGrid[ind6] &&  // ii = j.
             vert2 == hexaNodeIDsGrid[ind2] && vert3 == hexaNodeIDsGrid[ind1]) {  // jj = nPoly-k.
    e = l = nPolyConn;
    h = m = -1;
    c = 1;                                                                        // kk = nPoly-i.
  } else if (vert0 == hexaNodeIDsGrid[ind5] && vert1 == hexaNodeIDsGrid[ind1] &&  // ii = nPoly-k.
             vert2 == hexaNodeIDsGrid[ind2] && vert3 == hexaNodeIDsGrid[ind6]) {  // jj = j.
    a = l = nPolyConn;
    d = m = -1;
    g = 1;                                                                        // kk = nPoly-i.
  } else if (vert0 == hexaNodeIDsGrid[ind5] && vert1 == hexaNodeIDsGrid[ind1] &&  // ii = nPoly-k.
             vert2 == hexaNodeIDsGrid[ind0] && vert3 == hexaNodeIDsGrid[ind4]) {  // jj = nPoly-i.
    a = e = nPolyConn;
    d = f = -1;
    n = 1;                                                                        // kk = j.
  } else if (vert0 == hexaNodeIDsGrid[ind5] && vert1 == hexaNodeIDsGrid[ind4] &&  // ii = nPoly-i.
             vert2 == hexaNodeIDsGrid[ind0] && vert3 == hexaNodeIDsGrid[ind1]) {  // jj = nPoly-k.
    a = e = nPolyConn;
    b = h = -1;
    n = 1;  // kk = j.
  }

  else if (vert0 == hexaNodeIDsGrid[ind6] && vert1 == hexaNodeIDsGrid[ind7] &&  // ii = nPoly-i.
           vert2 == hexaNodeIDsGrid[ind4] && vert3 == hexaNodeIDsGrid[ind5]) {  // jj = nPoly-j.
    a = e = l = nPolyConn;
    b = g = o = -1;                                                               // kk = nPoly-k.
  } else if (vert0 == hexaNodeIDsGrid[ind6] && vert1 == hexaNodeIDsGrid[ind5] &&  // ii = nPoly-j.
             vert2 == hexaNodeIDsGrid[ind4] && vert3 == hexaNodeIDsGrid[ind7]) {  // jj = nPoly-i.
    a = e = l = nPolyConn;
    c = f = o = -1;                                                               // kk = nPoly-k.
  } else if (vert0 == hexaNodeIDsGrid[ind6] && vert1 == hexaNodeIDsGrid[ind7] &&  // ii = nPoly-i.
             vert2 == hexaNodeIDsGrid[ind3] && vert3 == hexaNodeIDsGrid[ind2]) {  // jj = nPoly-k.
    a = e = l = nPolyConn;
    b = h = n = -1;                                                               // kk = nPoly-j.
  } else if (vert0 == hexaNodeIDsGrid[ind6] && vert1 == hexaNodeIDsGrid[ind2] &&  // ii = nPoly-k.
             vert2 == hexaNodeIDsGrid[ind3] && vert3 == hexaNodeIDsGrid[ind7]) {  // jj = nPoly-i.
    a = e = l = nPolyConn;
    d = f = n = -1;                                                               // kk = nPoly-j.
  } else if (vert0 == hexaNodeIDsGrid[ind6] && vert1 == hexaNodeIDsGrid[ind2] &&  // ii = nPoly-k.
             vert2 == hexaNodeIDsGrid[ind1] && vert3 == hexaNodeIDsGrid[ind5]) {  // jj = nPoly-j.
    a = e = l = nPolyConn;
    d = g = m = -1;                                                               // kk = nPoly-i.
  } else if (vert0 == hexaNodeIDsGrid[ind6] && vert1 == hexaNodeIDsGrid[ind5] &&  // ii = nPoly-j.
             vert2 == hexaNodeIDsGrid[ind1] && vert3 == hexaNodeIDsGrid[ind2]) {  // jj = nPoly-k.
    a = e = l = nPolyConn;
    c = h = m = -1;  // kk = nPoly-i.
  }

  else if (vert0 == hexaNodeIDsGrid[ind7] && vert1 == hexaNodeIDsGrid[ind4] &&  // ii = nPoly-j.
           vert2 == hexaNodeIDsGrid[ind5] && vert3 == hexaNodeIDsGrid[ind6]) {  // jj = i.
    a = l = nPolyConn;
    c = o = -1;
    f = 1;                                                                        // kk = nPoly-k.
  } else if (vert0 == hexaNodeIDsGrid[ind7] && vert1 == hexaNodeIDsGrid[ind6] &&  // ii = i.
             vert2 == hexaNodeIDsGrid[ind5] && vert3 == hexaNodeIDsGrid[ind4]) {  // jj = nPoly-j.
    e = l = nPolyConn;
    g = o = -1;
    b = 1;                                                                        // kk = nPoly-k.
  } else if (vert0 == hexaNodeIDsGrid[ind7] && vert1 == hexaNodeIDsGrid[ind4] &&  // ii = nPoly-j.
             vert2 == hexaNodeIDsGrid[ind0] && vert3 == hexaNodeIDsGrid[ind3]) {  // jj = nPoly-k.
    a = e = nPolyConn;
    c = h = -1;
    m = 1;                                                                        // kk = i.
  } else if (vert0 == hexaNodeIDsGrid[ind7] && vert1 == hexaNodeIDsGrid[ind3] &&  // ii = nPoly-k.
             vert2 == hexaNodeIDsGrid[ind0] && vert3 == hexaNodeIDsGrid[ind4]) {  // jj = nPoly-j.
    a = e = nPolyConn;
    d = g = -1;
    m = 1;                                                                        // kk = i.
  } else if (vert0 == hexaNodeIDsGrid[ind7] && vert1 == hexaNodeIDsGrid[ind6] &&  // ii = i.
             vert2 == hexaNodeIDsGrid[ind2] && vert3 == hexaNodeIDsGrid[ind3]) {  // jj = nPoly-k.
    e = l = nPolyConn;
    h = n = -1;
    b = 1;                                                                        // kk = nPoly-j.
  } else if (vert0 == hexaNodeIDsGrid[ind7] && vert1 == hexaNodeIDsGrid[ind3] &&  // ii = nPoly-k.
             vert2 == hexaNodeIDsGrid[ind2] && vert3 == hexaNodeIDsGrid[ind6]) {  // jj = i.
    a = l = nPolyConn;
    d = n = -1;
    f = 1;  // kk = nPoly-j.
  }

  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Loop over the DOFs of the original hexahedron to create the connectivity
        of the hexahedron that corresponds to the new numbering. ---*/
  const unsigned short nn2 = (nPolyConn + 1) * (nPolyConn + 1);
  unsigned short ind = 0;
  for (unsigned short k = 0; k <= nPolyConn; ++k) {
    for (unsigned short j = 0; j <= nPolyConn; ++j) {
      for (unsigned short i = 0; i <= nPolyConn; ++i, ++ind) {
        /*--- Determine the ii, jj and kk indices of the new numbering, convert it to
              a 1D index and shore the modified index in modConnHexa. ---*/
        unsigned short ii = a + i * b + j * c + k * d;
        unsigned short jj = e + i * f + j * g + k * h;
        unsigned short kk = l + i * m + j * n + k * o;
        unsigned short iind = kk * nn2 + jj * (nPolyConn + 1) + ii;

        modConnHexa[iind] = connHexa[ind];
      }
    }
  }

  /*--- The quad corresponds to face 0 of the hexahedron. Hence the first
        nn2 entries in modConnHexa are the DOFs of the quad. Copy
        these entries from modConnHexa. ---*/
  for (unsigned short i = 0; i < nn2; ++i) modConnQuad[i] = modConnHexa[i];
}

void CMeshFEM_DG::CreateConnectivitiesQuadrilateralAdjacentPrism(
    const unsigned long* cornerPointsQuad, const unsigned short nPolyGrid,
    const vector<unsigned long>& prismNodeIDsGrid, const unsigned short nPolyConn, const unsigned long* connPrism,
    bool& swapFaceInElement, unsigned long* modConnQuad, unsigned long* modConnPrism) {
  /* Determine the indices of the six corner points of the prism. */
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPolyGrid;
  const unsigned short ind2 = (nPolyGrid + 1) * (nPolyGrid + 2) / 2 - 1;
  const unsigned short ind3 = (nPolyGrid + 1) * (nPolyGrid + 2) * nPolyGrid / 2;
  const unsigned short ind4 = ind1 + ind3;
  const unsigned short ind5 = ind2 + ind3;

  /* Easier storage of the four corner points of the quad in the new numbering. */
  const unsigned long vert0 = cornerPointsQuad[0];
  const unsigned long vert1 = cornerPointsQuad[1];
  const unsigned long vert2 = cornerPointsQuad[2];
  const unsigned long vert3 = cornerPointsQuad[3];

  /*--- There exists a linear mapping from the indices of the numbering used in
        prismNodeIDsGrid to the indices of the target numbering. This mapping is
        of the form ii = a + b*i + c*j, jj = d + e*i + f*j and kk = g + h*k,
        where ii,jj,kk are the indices of the new numbering and i,j,k the indices
        of the numbering used in in prismNodeIDsGrid (and connPrism). Note that
        the k-direction can at most be reversed compared to the original definition.
        This is due to the way the basis functions are defined for a prism. As a
        consequence also the variable swapFaceInElement is needed to cover all
        possible situations. The values of the coefficients a,b,c,d,e,f,g,h depend
        on how the corner points of the quad coincide with the corner points of the
        prism. This is determined below. The bool verticesDontMatch is there to
        check if vertices do not match. This should not happen, but it is checked
        for security. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;
  bool verticesDontMatch = false;

  if (vert0 == prismNodeIDsGrid[ind0] && vert1 == prismNodeIDsGrid[ind1] &&  // ii = i.
      vert2 == prismNodeIDsGrid[ind4] && vert3 == prismNodeIDsGrid[ind3]) {  // jj = j.
    b = f = h = 1;
    swapFaceInElement = false;                                                      // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind0] && vert1 == prismNodeIDsGrid[ind3] &&  // ii = i.
             vert2 == prismNodeIDsGrid[ind4] && vert3 == prismNodeIDsGrid[ind1]) {  // jj = j.
    b = f = h = 1;
    swapFaceInElement = true;                                                       // kk = k. Plus swap.
  } else if (vert0 == prismNodeIDsGrid[ind0] && vert1 == prismNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind5] && vert3 == prismNodeIDsGrid[ind3]) {  // jj = i.
    c = e = h = 1;
    swapFaceInElement = false;                                                      // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind0] && vert1 == prismNodeIDsGrid[ind3] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind5] && vert3 == prismNodeIDsGrid[ind2]) {  // jj = i.
    c = e = h = 1;
    swapFaceInElement = true;  // kk = k. Plus swap.
  }

  else if (vert0 == prismNodeIDsGrid[ind1] && vert1 == prismNodeIDsGrid[ind0] &&  // ii = nPoly-i-j.
           vert2 == prismNodeIDsGrid[ind3] && vert3 == prismNodeIDsGrid[ind4]) {  // jj = j.
    a = nPolyConn;
    b = c = -1;
    f = h = 1;
    swapFaceInElement = false;                                                      // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind1] && vert1 == prismNodeIDsGrid[ind4] &&  // ii = nPoly-i-j.
             vert2 == prismNodeIDsGrid[ind3] && vert3 == prismNodeIDsGrid[ind0]) {  // jj = j.
    a = nPolyConn;
    b = c = -1;
    f = h = 1;
    swapFaceInElement = true;                                                       // kk = k. Plus swap.
  } else if (vert0 == prismNodeIDsGrid[ind1] && vert1 == prismNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind5] && vert3 == prismNodeIDsGrid[ind4]) {  // jj = nPoly-i-j.
    d = nPolyConn;
    e = f = -1;
    c = h = 1;
    swapFaceInElement = false;                                                      // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind1] && vert1 == prismNodeIDsGrid[ind4] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind5] && vert3 == prismNodeIDsGrid[ind2]) {  // jj = nPoly-i-j.
    d = nPolyConn;
    e = f = -1;
    c = h = 1;
    swapFaceInElement = true;  // kk = k. Plus swap.
  }

  else if (vert0 == prismNodeIDsGrid[ind2] && vert1 == prismNodeIDsGrid[ind0] &&  // ii = nPoly-i-j.
           vert2 == prismNodeIDsGrid[ind3] && vert3 == prismNodeIDsGrid[ind5]) {  // jj = i.
    a = nPolyConn;
    b = c = -1;
    e = h = 1;
    swapFaceInElement = false;                                                      // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind2] && vert1 == prismNodeIDsGrid[ind5] &&  // ii = nPoly-i-j.
             vert2 == prismNodeIDsGrid[ind3] && vert3 == prismNodeIDsGrid[ind0]) {  // jj = i.
    a = nPolyConn;
    b = c = -1;
    e = h = 1;
    swapFaceInElement = true;                                                       // kk = k. Plus swap.
  } else if (vert0 == prismNodeIDsGrid[ind2] && vert1 == prismNodeIDsGrid[ind1] &&  // ii = i.
             vert2 == prismNodeIDsGrid[ind4] && vert3 == prismNodeIDsGrid[ind5]) {  // jj = nPoly-i-j.
    d = nPolyConn;
    e = f = -1;
    b = h = 1;
    swapFaceInElement = false;                                                      // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind2] && vert1 == prismNodeIDsGrid[ind5] &&  // ii = i.
             vert2 == prismNodeIDsGrid[ind4] && vert3 == prismNodeIDsGrid[ind1]) {  // jj = nPoly-i-j.
    d = nPolyConn;
    e = f = -1;
    b = h = 1;
    swapFaceInElement = true;  // kk = k. Plus swap.
  }

  else if (vert0 == prismNodeIDsGrid[ind3] && vert1 == prismNodeIDsGrid[ind4] &&  // ii = i.
           vert2 == prismNodeIDsGrid[ind1] && vert3 == prismNodeIDsGrid[ind0]) {  // jj = j.
    g = nPolyConn;
    b = f = 1;
    h = -1;
    swapFaceInElement = false;                                                      // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind3] && vert1 == prismNodeIDsGrid[ind0] &&  // ii = i.
             vert2 == prismNodeIDsGrid[ind1] && vert3 == prismNodeIDsGrid[ind4]) {  // jj = j.
    g = nPolyConn;
    b = f = 1;
    h = -1;
    swapFaceInElement = true;                                                       // kk = nPoly-k. Plus swap.
  } else if (vert0 == prismNodeIDsGrid[ind3] && vert1 == prismNodeIDsGrid[ind5] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind2] && vert3 == prismNodeIDsGrid[ind0]) {  // jj = i.
    g = nPolyConn;
    c = e = 1;
    h = -1;
    swapFaceInElement = false;                                                      // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind3] && vert1 == prismNodeIDsGrid[ind0] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind2] && vert3 == prismNodeIDsGrid[ind5]) {  // jj = i.
    g = nPolyConn;
    c = e = 1;
    h = -1;
    swapFaceInElement = true;  // kk = nPoly-k. Plus swap.
  }

  else if (vert0 == prismNodeIDsGrid[ind4] && vert1 == prismNodeIDsGrid[ind3] &&  // ii = nPoly-i-j.
           vert2 == prismNodeIDsGrid[ind0] && vert3 == prismNodeIDsGrid[ind1]) {  // jj = j.
    a = g = nPolyConn;
    b = c = h = -1;
    f = 1;
    swapFaceInElement = false;                                                      // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind4] && vert1 == prismNodeIDsGrid[ind1] &&  // ii = nPoly-i-j.
             vert2 == prismNodeIDsGrid[ind0] && vert3 == prismNodeIDsGrid[ind3]) {  // jj = j.
    a = g = nPolyConn;
    b = c = h = -1;
    f = 1;
    swapFaceInElement = true;                                                       // kk = nPoly-k. Plus swap.
  } else if (vert0 == prismNodeIDsGrid[ind4] && vert1 == prismNodeIDsGrid[ind5] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind2] && vert3 == prismNodeIDsGrid[ind1]) {  // jj = nPoly-i-j.
    d = g = nPolyConn;
    e = f = h = -1;
    c = 1;
    swapFaceInElement = false;                                                      // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind4] && vert1 == prismNodeIDsGrid[ind1] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind2] && vert3 == prismNodeIDsGrid[ind5]) {  // jj = nPoly-i-j.
    d = g = nPolyConn;
    e = f = h = -1;
    c = 1;
    swapFaceInElement = true;  // kk = nPoly-k. Plus swap.
  }

  else if (vert0 == prismNodeIDsGrid[ind5] && vert1 == prismNodeIDsGrid[ind3] &&  // ii = nPoly-i-j.
           vert2 == prismNodeIDsGrid[ind0] && vert3 == prismNodeIDsGrid[ind2]) {  // jj = i.
    a = g = nPolyConn;
    b = c = h = -1;
    e = 1;
    swapFaceInElement = false;                                                      // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind5] && vert1 == prismNodeIDsGrid[ind2] &&  // ii = nPoly-i-j.
             vert2 == prismNodeIDsGrid[ind0] && vert3 == prismNodeIDsGrid[ind3]) {  // jj = i.
    a = g = nPolyConn;
    b = c = h = -1;
    e = 1;
    swapFaceInElement = true;                                                       // kk = nPoly-k. Plus swap.
  } else if (vert0 == prismNodeIDsGrid[ind5] && vert1 == prismNodeIDsGrid[ind4] &&  // ii = i.
             vert2 == prismNodeIDsGrid[ind1] && vert3 == prismNodeIDsGrid[ind2]) {  // jj = nPoly-i-j.
    d = g = nPolyConn;
    e = f = h = -1;
    b = 1;
    swapFaceInElement = false;                                                      // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind5] && vert1 == prismNodeIDsGrid[ind2] &&  // ii = i.
             vert2 == prismNodeIDsGrid[ind1] && vert3 == prismNodeIDsGrid[ind4]) {  // jj = nPoly-i-j.
    d = g = nPolyConn;
    e = f = h = -1;
    b = 1;
    swapFaceInElement = true;  // kk = nPoly-k. Plus swap.
  }

  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Loop over the DOFs of the original prism to create the connectivity
        of the prism that corresponds to the new numbering. ---*/
  const unsigned short kOff = (nPolyConn + 1) * (nPolyConn + 2) / 2;
  unsigned short ind = 0;
  for (unsigned short k = 0; k <= nPolyConn; ++k) {
    for (unsigned short j = 0; j <= nPolyConn; ++j) {
      unsigned short uppBoundI = nPolyConn - j;
      for (unsigned short i = 0; i <= uppBoundI; ++i, ++ind) {
        /*--- Determine the ii, jj and kk indices of the new numbering, convert it to
              a 1D index and shore the modified index in modConnPrism. ---*/
        unsigned short ii = a + i * b + j * c;
        unsigned short jj = d + i * e + j * f;
        unsigned short kk = g + h * k;
        unsigned short iind = kk * kOff + jj * (nPolyConn + 1) + ii - jj * (jj - 1) / 2;

        modConnPrism[iind] = connPrism[ind];
      }
    }
  }

  /*--- Determine the connectivity of the quadrilateral face. ---*/
  if (swapFaceInElement) {
    /*--- The parametric coordinates r and s of the quad must be swapped
          w.r.t. to the parametric coordinates of the face of the prism.
          This means that the coordinate r of the quad runs from the base
          triangle to the top triangle of the prism. This corresponds to
          the k-direction of the prism. Hence the s-direction of the quad
          corresponds to the i-direction of the prism. ---*/
    for (unsigned short k = 0; k <= nPolyConn; ++k) {
      for (unsigned short i = 0; i <= nPolyConn; ++i) {
        const unsigned short iind = i * (nPolyConn + 1) + k;
        modConnQuad[iind] = modConnPrism[k * kOff + i];
      }
    }
  } else {
    /*--- The parametric coordinates r and s of the quad correspond to the
          parametric coordinates in i- and k-direction of the face of the
          prism. So an easy copy can be made. ---*/
    unsigned short iind = 0;
    for (unsigned short k = 0; k <= nPolyConn; ++k) {
      for (unsigned short i = 0; i <= nPolyConn; ++i, ++iind) {
        modConnQuad[iind] = modConnPrism[k * kOff + i];
      }
    }
  }
}

void CMeshFEM_DG::CreateConnectivitiesQuadrilateralAdjacentPyramid(
    const unsigned long* cornerPointsQuad, const unsigned short nPolyGrid, const vector<unsigned long>& pyraNodeIDsGrid,
    const unsigned short nPolyConn, const unsigned long* connPyra, unsigned long* modConnQuad,
    unsigned long* modConnPyra) {
  /* Determine the indices of the four corner points of the quadrilateral
     base of the pyramid. Note that the top of the pyramid is not needed
     in the comparison, because only triangular faces are attached to it. */
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPolyGrid;
  const unsigned short ind2 = (nPolyGrid + 1) * (nPolyGrid + 1) - 1;
  const unsigned short ind3 = ind2 - nPolyGrid;

  /* Easier storage of the four corner points of the quad in the new numbering. */
  const unsigned long vert0 = cornerPointsQuad[0];
  const unsigned long vert1 = cornerPointsQuad[1];
  const unsigned long vert2 = cornerPointsQuad[2];
  const unsigned long vert3 = cornerPointsQuad[3];

  /*--- There exists a linear mapping from the indices of the numbering used in
        pyraNodeIDsGrid to the indices of the target numbering. This mapping is of
        the form ii = a + b*i + c*j, jj = d + e*i + f*j and kk = k, where
        ii,jj,kk are the indices of the new numbering and i,j,k the indices of
        the numbering used in in pyraNodeIDsGrid (and connPyra). The values of
        the coefficients a,b,c,d,e,f depend on how the corner points of the quad
        coincide with the corner points of the quadrilateral base of the pyramid.
        This is determined below. The bool verticesDontMatch is there to check
        if vertices do not match. This should not happen, but it is checked for
        security. Note that this mapping is less complicated than for a prism or
        hexahedron, because a pyramid only has one quadrilateral face and hence
        the k-direction has to match. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
  bool verticesDontMatch = false;

  if (vert0 == pyraNodeIDsGrid[ind0] && vert1 == pyraNodeIDsGrid[ind1] &&  // ii = i.
      vert2 == pyraNodeIDsGrid[ind2] && vert3 == pyraNodeIDsGrid[ind3]) {  // jj = j.
    b = f = 1;
  } else if (vert0 == pyraNodeIDsGrid[ind0] && vert1 == pyraNodeIDsGrid[ind3] &&  // ii = j.
             vert2 == pyraNodeIDsGrid[ind2] && vert3 == pyraNodeIDsGrid[ind1]) {  // jj = i.
    c = e = 1;
  }

  else if (vert0 == pyraNodeIDsGrid[ind1] && vert1 == pyraNodeIDsGrid[ind2] &&  // ii = j.
           vert2 == pyraNodeIDsGrid[ind3] && vert3 == pyraNodeIDsGrid[ind0]) {  // jj = nPoly-i.
    d = nPolyConn;
    c = 1;
    e = -1;
  } else if (vert0 == pyraNodeIDsGrid[ind1] && vert1 == pyraNodeIDsGrid[ind0] &&  // ii = nPoly-i.
             vert2 == pyraNodeIDsGrid[ind3] && vert3 == pyraNodeIDsGrid[ind2]) {  // jj = j.
    a = nPolyConn;
    b = -1;
    f = 1;
  }

  else if (vert0 == pyraNodeIDsGrid[ind2] && vert1 == pyraNodeIDsGrid[ind3] &&  // ii = nPoly-i.
           vert2 == pyraNodeIDsGrid[ind0] && vert3 == pyraNodeIDsGrid[ind1]) {  // jj = nPoly-j.
    a = d = nPolyConn;
    b = f = -1;
  } else if (vert0 == pyraNodeIDsGrid[ind2] && vert1 == pyraNodeIDsGrid[ind1] &&  // ii = nPoly-j.
             vert2 == pyraNodeIDsGrid[ind0] && vert3 == pyraNodeIDsGrid[ind3]) {  // jj = nPoly-i.
    a = d = nPolyConn;
    c = e = -1;
  }

  else if (vert0 == pyraNodeIDsGrid[ind3] && vert1 == pyraNodeIDsGrid[ind0] &&  // ii = nPoly-j.
           vert2 == pyraNodeIDsGrid[ind1] && vert3 == pyraNodeIDsGrid[ind2]) {  // jj = i.
    a = nPolyConn;
    c = -1;
    e = 1;
  } else if (vert0 == pyraNodeIDsGrid[ind3] && vert1 == pyraNodeIDsGrid[ind2] &&  // ii = i.
             vert2 == pyraNodeIDsGrid[ind1] && vert3 == pyraNodeIDsGrid[ind0]) {  // jj = nPoly-j.
    d = nPolyConn;
    b = 1;
    f = -1;
  }

  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Loop over the DOFs of the original pyramid to create the connectivity
        of the pyramid that corresponds to the new numbering. Note that the
        outer k-loop is the same for both numberings. ---*/
  unsigned short mPoly = nPolyConn;
  unsigned short offLevel = 0;

  for (unsigned short k = 0; k <= nPolyConn; ++k, --mPoly) {
    unsigned short ind = offLevel;

    /* The variables of a and d in the transformation are actually flexible. */
    /* Create the correct values of a and d for the current k-level.         */
    const signed short aa = a ? mPoly : 0;
    const signed short dd = d ? mPoly : 0;

    /* Loop over the DOFs of the current quadrilateral. */
    for (unsigned short j = 0; j <= mPoly; ++j) {
      for (unsigned short i = 0; i <= mPoly; ++i, ++ind) {
        /*--- Determine the ii and jj indices of the new numbering, convert it
              to a 1D index and shore the modified index in modConnPyra. ---*/
        unsigned short ii = aa + i * b + j * c;
        unsigned short jj = dd + i * e + j * f;
        unsigned short iind = offLevel + jj * (mPoly + 1) + ii;

        modConnPyra[iind] = connPyra[ind];
      }
    }

    /* Update offLevel for the next quadrilateral level of the pyramid. */
    offLevel += (mPoly + 1) * (mPoly + 1);
  }

  /*--- The quad corresponds to face 0 of the pyramid. Hence the first
        nn2 entries in modConnPyra are the DOFs of the quad. Copy
        these entries from modConnPyra. ---*/
  const unsigned short nn2 = (nPolyConn + 1) * (nPolyConn + 1);
  for (unsigned short i = 0; i < nn2; ++i) modConnQuad[i] = modConnPyra[i];
}

void CMeshFEM_DG::CreateConnectivitiesTriangleAdjacentPrism(const unsigned long* cornerPointsTria,
                                                            const unsigned short nPolyGrid,
                                                            const vector<unsigned long>& prismNodeIDsGrid,
                                                            const unsigned short nPolyConn,
                                                            const unsigned long* connPrism, unsigned long* modConnTria,
                                                            unsigned long* modConnPrism) {
  /* Determine the indices of the six corner points of the prism. */
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPolyGrid;
  const unsigned short ind2 = (nPolyGrid + 1) * (nPolyGrid + 2) / 2 - 1;
  const unsigned short ind3 = (nPolyGrid + 1) * (nPolyGrid + 2) * nPolyGrid / 2;
  const unsigned short ind4 = ind1 + ind3;
  const unsigned short ind5 = ind2 + ind3;

  /* Easier storage of the three corner points of the triangle in the new numbering. */
  const unsigned long vert0 = cornerPointsTria[0];
  const unsigned long vert1 = cornerPointsTria[1];
  const unsigned long vert2 = cornerPointsTria[2];

  /*--- There exists a linear mapping from the indices of the numbering used in
        prismNodeIDsGrid to the indices of the target numbering. This mapping is
        of the form ii = a + b*i + c*j, jj = d + e*i + f*j and kk = g + h*k,
        where ii,jj,kk are the indices of the new numbering and i,j,k the indices
        of the numbering used in in prismNodeIDsGrid (and connPrism). Note that
        the k-direction can at most be reversed compared to the original definition.
        This is due to the way the basis functions are defined for a prism.
        The values of the coefficients a,b,c,d,e,f,g,h depend on how the corner
        points of the triangle coincide with the corner points of the prism. This
        is determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;
  bool verticesDontMatch = false;

  if (vert0 == prismNodeIDsGrid[ind0] && vert1 == prismNodeIDsGrid[ind1] &&         // ii = i.
      vert2 == prismNodeIDsGrid[ind2]) {                                            // jj = j.
    b = f = h = 1;                                                                  // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind0] && vert1 == prismNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind1]) {                                     // jj = i.
    c = e = h = 1;                                                                  // kk = k.
  }

  else if (vert0 == prismNodeIDsGrid[ind1] && vert1 == prismNodeIDsGrid[ind0] &&  // ii = nPoly-i-j.
           vert2 == prismNodeIDsGrid[ind2]) {                                     // jj = j.
    a = nPolyConn;
    b = c = -1;
    f = h = 1;                                                                      // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind1] && vert1 == prismNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind0]) {                                     // jj = nPoly-i-j.
    d = nPolyConn;
    e = f = -1;
    c = h = 1;  // kk = k.
  }

  else if (vert0 == prismNodeIDsGrid[ind2] && vert1 == prismNodeIDsGrid[ind0] &&  // ii = nPoly-i-j.
           vert2 == prismNodeIDsGrid[ind1]) {                                     // jj = i.
    a = nPolyConn;
    b = c = -1;
    e = h = 1;                                                                      // kk = k.
  } else if (vert0 == prismNodeIDsGrid[ind2] && vert1 == prismNodeIDsGrid[ind1] &&  // ii = i.
             vert2 == prismNodeIDsGrid[ind0]) {                                     // jj = nPoly-i-j.
    d = nPolyConn;
    e = f = -1;
    b = h = 1;  // kk = k.
  }

  else if (vert0 == prismNodeIDsGrid[ind3] && vert1 == prismNodeIDsGrid[ind4] &&  // ii = i.
           vert2 == prismNodeIDsGrid[ind5]) {                                     // jj = j.
    g = nPolyConn;
    b = f = 1;
    h = -1;                                                                         // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind3] && vert1 == prismNodeIDsGrid[ind5] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind4]) {                                     // jj = i.
    g = nPolyConn;
    c = e = 1;
    h = -1;  // kk = nPoly-k.
  }

  else if (vert0 == prismNodeIDsGrid[ind4] && vert1 == prismNodeIDsGrid[ind3] &&  // ii = nPoly-i-j.
           vert2 == prismNodeIDsGrid[ind5]) {                                     // jj = j.
    a = g = nPolyConn;
    b = c = h = -1;
    f = 1;                                                                          // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind4] && vert1 == prismNodeIDsGrid[ind5] &&  // ii = j.
             vert2 == prismNodeIDsGrid[ind3]) {                                     // jj = nPoly-i-j.
    d = g = nPolyConn;
    e = f = h = -1;
    c = 1;  // kk = nPoly-k.
  }

  else if (vert0 == prismNodeIDsGrid[ind5] && vert1 == prismNodeIDsGrid[ind3] &&  // ii = nPoly-i-j.
           vert2 == prismNodeIDsGrid[ind4]) {                                     // jj = i.
    a = g = nPolyConn;
    b = c = h = -1;
    e = 1;                                                                          // kk = nPoly-k.
  } else if (vert0 == prismNodeIDsGrid[ind5] && vert1 == prismNodeIDsGrid[ind4] &&  // ii = i.
             vert2 == prismNodeIDsGrid[ind3]) {                                     // jj = nPoly-i-j.
    d = g = nPolyConn;
    e = f = h = -1;
    b = 1;  // kk = nPoly-k.
  }

  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Loop over the DOFs of the original prism to create the connectivity
        of the prism that corresponds to the new numbering. ---*/
  const unsigned short kOff = (nPolyConn + 1) * (nPolyConn + 2) / 2;
  unsigned short ind = 0;
  for (unsigned short k = 0; k <= nPolyConn; ++k) {
    for (unsigned short j = 0; j <= nPolyConn; ++j) {
      unsigned short uppBoundI = nPolyConn - j;
      for (unsigned short i = 0; i <= uppBoundI; ++i, ++ind) {
        /*--- Determine the ii, jj and kk indices of the new numbering, convert it to
              a 1D index and shore the modified index in modConnPrism. ---*/
        unsigned short ii = a + i * b + j * c;
        unsigned short jj = d + i * e + j * f;
        unsigned short kk = g + h * k;
        unsigned short iind = kk * kOff + jj * (nPolyConn + 1) + ii - jj * (jj - 1) / 2;

        modConnPrism[iind] = connPrism[ind];
      }
    }
  }

  /*--- The triangle corresponds to face 0 of the prism. Hence the first
        kOff entries in modConnPrism are the DOFs of the triangle. Copy
        these entries from modConnPrism. ---*/
  for (unsigned short i = 0; i < kOff; ++i) modConnTria[i] = modConnPrism[i];
}

void CMeshFEM_DG::CreateConnectivitiesTriangleAdjacentPyramid(const unsigned long* cornerPointsTria,
                                                              const unsigned short nPolyGrid,
                                                              const vector<unsigned long>& pyraNodeIDsGrid,
                                                              const unsigned short nPolyConn,
                                                              const unsigned long* connPyra, bool& swapFaceInElement,
                                                              unsigned long* modConnTria, unsigned long* modConnPyra) {
  /* Determine the indices of the five corner points of the pyramid. */
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPolyGrid;
  const unsigned short ind2 = (nPolyGrid + 1) * (nPolyGrid + 1) - 1;
  const unsigned short ind3 = ind2 - nPolyGrid;
  const unsigned short ind4 = (nPolyGrid + 1) * (nPolyGrid + 2) * (2 * nPolyGrid + 3) / 6 - 1;

  /* Check if the top of the pyramid coincides with either corner point 1 or
     corner point 2 of the triangle. Set swapFaceInElement accordingly. */
  if (cornerPointsTria[1] == pyraNodeIDsGrid[ind4])
    swapFaceInElement = true;
  else if (cornerPointsTria[2] == pyraNodeIDsGrid[ind4])
    swapFaceInElement = false;
  else
    SU2_MPI::Error(string("Top of the pyramid does not coincide with either corner point 1 or 2.\n") +
                       string("This should not happen"),
                   CURRENT_FUNCTION);

  /* Easier storage of the two other corner points of the triangle in the new numbering.
     vert0 always corresponds to cornerPointsTria[0], while vert1 contains the other
     point that is not equal to the top of the pyramid. This depends on swapFaceInElement. */
  const unsigned long vert0 = cornerPointsTria[0];
  const unsigned long vert1 = swapFaceInElement ? cornerPointsTria[2] : cornerPointsTria[1];

  /*--- There exists a linear mapping from the indices of the numbering used in
        pyraNodeIDsGrid to the indices of the target numbering. This mapping is of
        the form ii = a + b*i + c*j, jj = d + e*i + f*j and kk = k, where
        ii,jj,kk are the indices of the new numbering and i,j,k the indices of
        the numbering used in in pyraNodeIDsGrid (and connPyra). The values of
        the coefficients a,b,c,d,e,f depend on how the corner points of the triangle
        coincide with the corner points of the triangular face of the pyramid.
        This is determined below. The bool verticesDontMatch is there to check
        if vertices do not match. This should not happen, but it is checked for
        security. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
  bool verticesDontMatch = false;

  if (vert0 == pyraNodeIDsGrid[ind0] && vert1 == pyraNodeIDsGrid[ind1]) {         // ii = i.
    b = f = 1;                                                                    // jj = j.
  } else if (vert0 == pyraNodeIDsGrid[ind0] && vert1 == pyraNodeIDsGrid[ind3]) {  // ii = j.
    c = e = 1;                                                                    // jj = i.
  }

  else if (vert0 == pyraNodeIDsGrid[ind1] && vert1 == pyraNodeIDsGrid[ind2]) {  // ii = j.
    d = nPolyConn;
    c = 1;
    e = -1;                                                                       // jj = nPoly-i.
  } else if (vert0 == pyraNodeIDsGrid[ind1] && vert1 == pyraNodeIDsGrid[ind0]) {  // ii = nPoly-i.
    a = nPolyConn;
    b = -1;
    f = 1;  // jj = j.
  }

  else if (vert0 == pyraNodeIDsGrid[ind2] && vert1 == pyraNodeIDsGrid[ind3]) {  // ii = nPoly-i.
    a = d = nPolyConn;
    b = f = -1;                                                                   // jj = nPoly-j.
  } else if (vert0 == pyraNodeIDsGrid[ind2] && vert1 == pyraNodeIDsGrid[ind1]) {  // ii = nPoly-j.
    a = d = nPolyConn;
    c = e = -1;  // jj = nPoly-i.
  }

  else if (vert0 == pyraNodeIDsGrid[ind3] && vert1 == pyraNodeIDsGrid[ind0]) {  // ii = nPoly-j.
    a = nPolyConn;
    c = -1;
    e = 1;                                                                        // jj = i.
  } else if (vert0 == pyraNodeIDsGrid[ind3] && vert1 == pyraNodeIDsGrid[ind2]) {  // ii = i.
    d = nPolyConn;
    b = 1;
    f = -1;  // jj = nPoly-j.
  }

  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Loop over the DOFs of the original pyramid to create the connectivity
        of the pyramid that corresponds to the new numbering. Note that the
        outer k-loop is the same for both numberings. ---*/
  unsigned short mPoly = nPolyConn;
  unsigned short offLevel = 0;

  for (unsigned short k = 0; k <= nPolyConn; ++k, --mPoly) {
    unsigned short ind = offLevel;

    /* The variables of a and d in the transformation are actually flexible. */
    /* Create the correct values of a and d for the current k-level.         */
    const signed short aa = a ? mPoly : 0;
    const signed short dd = d ? mPoly : 0;

    /* Loop over the DOFs of the current quadrilateral. */
    for (unsigned short j = 0; j <= mPoly; ++j) {
      for (unsigned short i = 0; i <= mPoly; ++i, ++ind) {
        /*--- Determine the ii and jj indices of the new numbering, convert it
              to a 1D index and shore the modified index in modConnPyra. ---*/
        unsigned short ii = aa + i * b + j * c;
        unsigned short jj = dd + i * e + j * f;
        unsigned short iind = offLevel + jj * (mPoly + 1) + ii;

        modConnPyra[iind] = connPyra[ind];
      }
    }

    /* Update offLevel for the next quadrilateral level of the pyramid. */
    offLevel += (mPoly + 1) * (mPoly + 1);
  }

  /*--- Determine the connectivity of the triangular face. ---*/
  if (swapFaceInElement) {
    /*--- The parametric coordinates r and s of the triangle must be swapped
          w.r.t. to the parametric coordinates of the face of the pyramid.
          This means that the coordinate r of the triangle runs from the base
          to the top of the pyramid. This corresponds to the k-direction of the
          pyramid. Hence the s-direction of the triangle corresponds to the
          i-direction of the pyramid. ---*/
    mPoly = nPolyConn;
    offLevel = 0;
    for (unsigned short k = 0; k <= nPolyConn; ++k, --mPoly) {
      for (unsigned short i = 0; i <= mPoly; ++i) {
        unsigned short iind = i * (nPolyConn + 1) + k - i * (i - 1) / 2;
        modConnTria[iind] = modConnPyra[offLevel + i];
      }
      offLevel += (mPoly + 1) * (mPoly + 1);
    }
  } else {
    /*--- The parametric coordinates r and s of the triangle correspond to the
          parametric coordinates in i- and k-direction of the face of the
          pyramid. So an easy copy can be made. ---*/
    mPoly = nPolyConn;
    offLevel = 0;
    unsigned short iind = 0;
    for (unsigned short k = 0; k <= nPolyConn; ++k, --mPoly) {
      for (unsigned short i = 0; i <= mPoly; ++i, ++iind) {
        modConnTria[iind] = modConnPyra[offLevel + i];
      }
      offLevel += (mPoly + 1) * (mPoly + 1);
    }
  }
}

void CMeshFEM_DG::CreateConnectivitiesTriangleAdjacentTetrahedron(
    const unsigned long* cornerPointsTria, const unsigned short nPolyGrid, const vector<unsigned long>& tetNodeIDsGrid,
    const unsigned short nPolyConn, const unsigned long* connTet, unsigned long* modConnTria,
    unsigned long* modConnTet) {
  /* Determine the indices of the four corner points of the tetrahedron. */
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPolyGrid;
  const unsigned short ind2 = (nPolyGrid + 1) * (nPolyGrid + 2) / 2 - 1;
  const unsigned short ind3 = (nPolyGrid + 1) * (nPolyGrid + 2) * (nPolyGrid + 3) / 6 - 1;

  /* Easier storage of the three corner points of the triangle in the new numbering. */
  const unsigned long vert0 = cornerPointsTria[0];
  const unsigned long vert1 = cornerPointsTria[1];
  const unsigned long vert2 = cornerPointsTria[2];

  /*--- There exists a linear mapping from the indices of the numbering used in
        tetNodeIDsGrid to the indices of the target numbering. This mapping is of
        the form ii = a + b*i + c*j + d*k, jj = e + f*i + g*j + h*k and
        kk = l + m*i + n*j + o*k, where ii,jj,kk are the indices of the new
        numbering and i,j,k the indices of the numbering used in in tetNodeIDsGrid
        (and connTet). The values of the coefficients a,b,c,d,e,f,g,h,l,m,n,o
        depend on how the corner points of the triangle coincide with the corner points
        of the tetrahedron. This is determined below. The bool verticesDontMatch is
        there to check if vertices do not match. This should not happen, but
        it is checked for security. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, l = 0, m = 0, n = 0, o = 0;
  bool verticesDontMatch = false;

  if (vert0 == tetNodeIDsGrid[ind0] && vert1 == tetNodeIDsGrid[ind1] &&         // ii = i.
      vert2 == tetNodeIDsGrid[ind2]) {                                          // jj = j.
    b = g = o = 1;                                                              // kk = k.
  } else if (vert0 == tetNodeIDsGrid[ind0] && vert1 == tetNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == tetNodeIDsGrid[ind1]) {                                   // jj = i.
    c = f = o = 1;                                                              // kk = k.
  } else if (vert0 == tetNodeIDsGrid[ind0] && vert1 == tetNodeIDsGrid[ind1] &&  // ii = i.
             vert2 == tetNodeIDsGrid[ind3]) {                                   // jj = k.
    b = h = n = 1;                                                              // kk = j.
  } else if (vert0 == tetNodeIDsGrid[ind0] && vert1 == tetNodeIDsGrid[ind3] &&  // ii = k.
             vert2 == tetNodeIDsGrid[ind1]) {                                   // jj = i.
    d = f = n = 1;                                                              // kk = j.
  } else if (vert0 == tetNodeIDsGrid[ind0] && vert1 == tetNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == tetNodeIDsGrid[ind3]) {                                   // jj = k.
    c = h = m = 1;                                                              // kk = i.
  } else if (vert0 == tetNodeIDsGrid[ind0] && vert1 == tetNodeIDsGrid[ind3] &&  // ii = k.
             vert2 == tetNodeIDsGrid[ind2]) {                                   // jj = j.
    d = g = m = 1;                                                              // kk = i.
  }

  else if (vert0 == tetNodeIDsGrid[ind1] && vert1 == tetNodeIDsGrid[ind0] &&  // ii = nPoly-i-j-k.
           vert2 == tetNodeIDsGrid[ind2]) {                                   // jj = j.
    a = nPolyConn;
    b = c = d = -1;
    g = o = 1;                                                                  // kk = k.
  } else if (vert0 == tetNodeIDsGrid[ind1] && vert1 == tetNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == tetNodeIDsGrid[ind0]) {                                   // jj = nPoly-i-j-k.
    e = nPolyConn;
    f = g = h = -1;
    c = o = 1;                                                                  // kk = k.
  } else if (vert0 == tetNodeIDsGrid[ind1] && vert1 == tetNodeIDsGrid[ind0] &&  // ii = nPoly-i-j-k.
             vert2 == tetNodeIDsGrid[ind3]) {                                   // jj = k.
    a = nPolyConn;
    b = c = d = -1;
    h = n = 1;                                                                  // kk = j.
  } else if (vert0 == tetNodeIDsGrid[ind1] && vert1 == tetNodeIDsGrid[ind3] &&  // ii = k.
             vert2 == tetNodeIDsGrid[ind0]) {                                   // jj = nPoly-i-j-k.
    e = nPolyConn;
    f = g = h = -1;
    d = n = 1;                                                                  // kk = j.
  } else if (vert0 == tetNodeIDsGrid[ind1] && vert1 == tetNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == tetNodeIDsGrid[ind3]) {                                   // jj = k.
    l = nPolyConn;
    m = n = o = -1;
    c = h = 1;                                                                  // kk = nPoly-i-j-k.
  } else if (vert0 == tetNodeIDsGrid[ind1] && vert1 == tetNodeIDsGrid[ind3] &&  // ii = k.
             vert2 == tetNodeIDsGrid[ind2]) {                                   // jj = j.
    l = nPolyConn;
    m = n = o = -1;
    d = g = 1;  // kk = nPoly-i-j-k.
  }

  else if (vert0 == tetNodeIDsGrid[ind2] && vert1 == tetNodeIDsGrid[ind0] &&  // ii = nPoly-i-j-k.
           vert2 == tetNodeIDsGrid[ind1]) {                                   // jj = i.
    a = nPolyConn;
    b = c = d = -1;
    f = o = 1;                                                                  // kk = k.
  } else if (vert0 == tetNodeIDsGrid[ind2] && vert1 == tetNodeIDsGrid[ind1] &&  // ii = i.
             vert2 == tetNodeIDsGrid[ind0]) {                                   // jj = nPoly-i-j-k.
    e = nPolyConn;
    f = g = h = -1;
    b = o = 1;                                                                  // kk = k.
  } else if (vert0 == tetNodeIDsGrid[ind2] && vert1 == tetNodeIDsGrid[ind0] &&  // ii = nPoly-i-j-k.
             vert2 == tetNodeIDsGrid[ind3]) {                                   // jj = k.
    a = nPolyConn;
    b = c = d = -1;
    h = m = 1;                                                                  // kk = i.
  } else if (vert0 == tetNodeIDsGrid[ind2] && vert1 == tetNodeIDsGrid[ind3] &&  // ii = k.
             vert2 == tetNodeIDsGrid[ind0]) {                                   // jj = nPoly-i-j-k.
    e = nPolyConn;
    f = g = h = -1;
    d = m = 1;                                                                  // kk = i.
  } else if (vert0 == tetNodeIDsGrid[ind2] && vert1 == tetNodeIDsGrid[ind0] &&  // ii = nPoly-i-j-k.
             vert2 == tetNodeIDsGrid[ind3]) {                                   // jj = k.
    a = nPolyConn;
    b = c = d = -1;
    h = m = 1;                                                                  // kk = i.
  } else if (vert0 == tetNodeIDsGrid[ind2] && vert1 == tetNodeIDsGrid[ind1] &&  // ii = i.
             vert2 == tetNodeIDsGrid[ind3]) {                                   // jj = k.
    l = nPolyConn;
    m = n = o = -1;
    b = h = 1;                                                                  // kk = nPoly-i-j-k.
  } else if (vert0 == tetNodeIDsGrid[ind2] && vert1 == tetNodeIDsGrid[ind3] &&  // ii = k.
             vert2 == tetNodeIDsGrid[ind1]) {                                   // jj = i.
    l = nPolyConn;
    m = n = o = -1;
    d = f = 1;  // kk = nPoly-i-j-k.
  }

  else if (vert0 == tetNodeIDsGrid[ind3] && vert1 == tetNodeIDsGrid[ind0] &&  // ii = nPoly-i-j-k.
           vert2 == tetNodeIDsGrid[ind1]) {                                   // jj = i.
    a = nPolyConn;
    b = c = d = -1;
    f = n = 1;                                                                  // kk = j.
  } else if (vert0 == tetNodeIDsGrid[ind3] && vert1 == tetNodeIDsGrid[ind1] &&  // ii = i.
             vert2 == tetNodeIDsGrid[ind0]) {                                   // jj = nPoly-i-j-k.
    e = nPolyConn;
    f = g = h = -1;
    b = n = 1;                                                                  // kk = j.
  } else if (vert0 == tetNodeIDsGrid[ind3] && vert1 == tetNodeIDsGrid[ind0] &&  // ii = nPoly-i-j-k.
             vert2 == tetNodeIDsGrid[ind2]) {                                   // jj = j.
    a = nPolyConn;
    b = c = d = -1;
    g = m = 1;                                                                  // kk = i.
  } else if (vert0 == tetNodeIDsGrid[ind3] && vert1 == tetNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == tetNodeIDsGrid[ind0]) {                                   // jj = nPoly-i-j-k.
    e = nPolyConn;
    f = g = h = -1;
    c = m = 1;                                                                  // kk = i.
  } else if (vert0 == tetNodeIDsGrid[ind3] && vert1 == tetNodeIDsGrid[ind1] &&  // ii = i.
             vert2 == tetNodeIDsGrid[ind2]) {                                   // jj = j.
    l = nPolyConn;
    m = n = o = -1;
    b = g = 1;                                                                  // kk = nPoly-i-j-k.
  } else if (vert0 == tetNodeIDsGrid[ind3] && vert1 == tetNodeIDsGrid[ind2] &&  // ii = j.
             vert2 == tetNodeIDsGrid[ind1]) {                                   // jj = i.
    l = nPolyConn;
    m = n = o = -1;
    c = f = 1;  // kk = nPoly-i-j-k.
  }

  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Some constants to convert the (ii,jj,kk) indices to a 1D index. ---*/
  const unsigned short abv1 = (11 + 12 * nPolyConn + 3 * nPolyConn * nPolyConn);
  const unsigned short abv2 = (2 * nPolyConn + 3) * 3;
  const unsigned short abv3 = (nPolyConn + 2) * 3;

  /*--- Loop over the DOFs of the original tetrahedron to create the connectivity
        of the tetrahedron that corresponds to the new numbering. ---*/
  unsigned short ind = 0;
  for (unsigned short k = 0; k <= nPolyConn; ++k) {
    unsigned short uppBoundJ = nPolyConn - k;
    for (unsigned short j = 0; j <= uppBoundJ; ++j) {
      unsigned short uppBoundI = nPolyConn - k - j;
      for (unsigned short i = 0; i <= uppBoundI; ++i, ++ind) {
        /*--- Determine the ii, jj and kk indices of the new numbering, convert it to
              a 1D index and shore the modified index in modConnTet. ---*/
        unsigned short ii = a + i * b + j * c + k * d;
        unsigned short jj = e + i * f + j * g + k * h;
        unsigned short kk = l + i * m + j * n + k * o;
        unsigned short iind =
            (abv1 * kk + abv2 * jj + 6 * ii - abv3 * kk * kk - 6 * kk * jj - 3 * jj * jj + kk * kk * kk) / 6;

        modConnTet[iind] = connTet[ind];
      }
    }
  }

  /*--- The triangle corresponds to face 0 of the tetrahedron. Hence the first
        nn2 entries in modConnTet are the DOFs of the triangle. Copy these
        entries from modConnTet. ---*/
  const unsigned short nn2 = (nPolyConn + 1) * (nPolyConn + 2) / 2;
  for (unsigned short i = 0; i < nn2; ++i) modConnTria[i] = modConnTet[i];
}

void CMeshFEM_DG::MetricTermsMatchingFaces(CConfig* config) {
  /* Determine whether or not the viscous terms are needed. */
  bool viscousTerms = (config->GetKind_Solver() != MAIN_SOLVER::FEM_EULER &&
                       config->GetKind_Solver() != MAIN_SOLVER::DISC_ADJ_FEM_EULER);

  /* Loop over the internal matching faces. */
  for (unsigned long i = 0; i < matchingFaces.size(); ++i) {
    /*--------------------------------------------------------------------------*/
    /*--- Step 1: Allocate the memory for the face metric terms.             ---*/
    /*---         Depending on the case, not all of this memory is needed.   ---*/
    /*---         - Unit normals + area (nDim+1 per integration point)       ---*/
    /*---         - drdx, dsdx, etc. for both sides (nDim*nDim per           ---*/
    /*---                                            integration point)      ---*/
    /*---         - grid velocities in the integration points.               ---*/
    /*--------------------------------------------------------------------------*/

    /* Determine the corresponding standard face element and get the
       relevant information from it. */
    const unsigned short ind = matchingFaces[i].indStandardElement;
    const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();

    /* Allocate the memory for the vectors to store the different face
       metric terms. */
    matchingFaces[i].metricNormalsFace.resize(nInt * (nDim + 1));

    if (viscousTerms) {
      matchingFaces[i].metricCoorDerivFace0.resize(nInt * nDim * nDim);
      matchingFaces[i].metricCoorDerivFace1.resize(nInt * nDim * nDim);
    }

    /* Allocate the memory for the grid velocities and initialize them to the
       default value of zero. */
    matchingFaces[i].gridVelocities.assign(nInt * nDim, 0.0);

    /*--------------------------------------------------------------------------*/
    /*--- Step 2: Determine the actual metric data in the integration points ---*/
    /*---         of the faces.                                              ---*/
    /*--------------------------------------------------------------------------*/

    /* Call the function ComputeNormalsFace to compute the unit normals and
       its corresponding area in the integration points. The data from side 0
       is used, but the data from side 1 should give the same result. */
    unsigned short nDOFs = standardMatchingFacesGrid[ind].GetNDOFsFaceSide0();
    const su2double* dr = standardMatchingFacesGrid[ind].GetDrBasisFaceIntegrationSide0();
    const su2double* ds = standardMatchingFacesGrid[ind].GetDsBasisFaceIntegrationSide0();

    ComputeNormalsFace(nInt, nDOFs, dr, ds, matchingFaces[i].DOFsGridFaceSide0.data(),
                       matchingFaces[i].metricNormalsFace.data());

    /* Check whether the viscous terms are needed in the discretization. */
    if (viscousTerms) {
      /* Compute the derivatives of the parametric coordinates w.r.t. the
         Cartesian coordinates, i.e. drdx, drdy, etc. in the integration points
         on side 0 of the face. */
      nDOFs = standardMatchingFacesGrid[ind].GetNDOFsElemSide0();
      dr = standardMatchingFacesGrid[ind].GetMatDerBasisElemIntegrationSide0();

      ComputeGradientsCoordinatesFace(nInt, nDOFs, dr, matchingFaces[i].DOFsGridElementSide0.data(),
                                      matchingFaces[i].metricCoorDerivFace0.data(), config);

      /* Compute the derivatives of the parametric coordinates w.r.t. the
         Cartesian coordinates, i.e. drdx, drdy, etc. in the integration points
         on side 1 of the face. */
      nDOFs = standardMatchingFacesGrid[ind].GetNDOFsElemSide1();
      dr = standardMatchingFacesGrid[ind].GetMatDerBasisElemIntegrationSide1();

      ComputeGradientsCoordinatesFace(nInt, nDOFs, dr, matchingFaces[i].DOFsGridElementSide1.data(),
                                      matchingFaces[i].metricCoorDerivFace1.data(), config);
    }
  }
}

void CMeshFEM_DG::LengthScaleVolumeElements() {
  /* Initialize the length scale of the elements to zero. */
  for (unsigned long i = 0; i < nVolElemTot; ++i) volElem[i].lenScale = 0.0;

  /*-------------------------------------------------------------------*/
  /*--- Step 1: Determine the maximum value of the Jacobians of the ---*/
  /*---         surrounding faces of the elements.                  ---*/
  /*-------------------------------------------------------------------*/

  /*--- First the physical boundaries. Loop over them. ---*/
  for (unsigned short j = 0; j < boundaries.size(); ++j) {
    if (!boundaries[j].periodicBoundary) {
      for (unsigned long i = 0; i < boundaries[j].surfElem.size(); ++i) {
        /* Easier storage of corresponding volume element and the metric data
           for the normals of this surface element, determine the number of
           integration points for this surface element and loop over them. */
        const unsigned long volInd = boundaries[j].surfElem[i].volElemID;
        const su2double* normals = boundaries[j].surfElem[i].metricNormalsFace.data();
        const unsigned short ind = boundaries[j].surfElem[i].indStandardElement;
        const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

        /* Loop over the integration points and determine the maximum area
           scale for this face. This is stored in the variable lenScale of
           the corresponding volume element. */
        for (unsigned short k = 0; k < nInt; ++k)
          volElem[volInd].lenScale = max(volElem[volInd].lenScale, normals[k * (nDim + 1) + nDim]);
      }
    }
  }

  /*--- Loop over the matching internal faces to update the maximum Jacobian
        of the surrounding faces. ---*/
  for (unsigned long i = 0; i < matchingFaces.size(); ++i) {
    /* Easier storage of the ID's of the two adjacent elements, determine the
       corresponding standard face element, the number of integration points
       and store the metric data of the face a bit easier. */
    const unsigned long v0 = matchingFaces[i].elemID0;
    const unsigned long v1 = matchingFaces[i].elemID1;
    const su2double* normals = matchingFaces[i].metricNormalsFace.data();
    const unsigned short ind = matchingFaces[i].indStandardElement;
    const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();

    /* Loop over the integration points and determine the maximum area
       scale for this face. This is stored in the variable lenScale of
       the neigbhoring volume elements. */
    for (unsigned short k = 0; k < nInt; ++k) {
      const su2double area = normals[k * (nDim + 1) + nDim];
      volElem[v0].lenScale = max(volElem[v0].lenScale, area);
      volElem[v1].lenScale = max(volElem[v1].lenScale, area);
    }
  }

  /*-------------------------------------------------------------------*/
  /*--- Step 2: Communicate the currently stored data of the halo   ---*/
  /*---         elements to the ranks where the original is stored. ---*/
  /*-------------------------------------------------------------------*/

#ifdef HAVE_MPI

  /*--- Define the buffers needed for the communication. ---*/
  vector<SU2_MPI::Request> sendRequests(ranksSend.size());
  vector<SU2_MPI::Request> recvRequests(ranksRecv.size());

  vector<vector<su2double> > sendBuf(ranksSend.size());
  vector<vector<su2double> > recvBuf(ranksRecv.size());

  /*--- Send the data using non-blocking sends. Exclude self communication.
        Note that the inverse communication pattern must be used for this
        purpose.  ---*/
  unsigned long nRecvRequests = 0;
  for (unsigned long i = 0; i < ranksRecv.size(); ++i) {
    if (ranksRecv[i] != rank) {
      /* Allocate the memory for communication buffer. */
      recvBuf[i].resize(entitiesRecv[i].size());

      /* Fill the buffer to be sent, which is the receive
         information of the original communication pattern. */
      for (unsigned long j = 0; j < entitiesRecv[i].size(); ++j) recvBuf[i][j] = volElem[entitiesRecv[i][j]].lenScale;

      /* Send the data. */
      SU2_MPI::Isend(recvBuf[i].data(), recvBuf[i].size(), MPI_DOUBLE, ranksRecv[i], ranksRecv[i] + 10,
                     SU2_MPI::GetComm(), &recvRequests[nRecvRequests++]);
    }
  }

  /*--- Post the receives using non-blocking receives. Exclude self communication.
        Also here the inverse communication pattern must be used.  ---*/
  unsigned long nSendRequests = 0;
  for (unsigned long i = 0; i < ranksSend.size(); ++i) {
    if (ranksSend[i] != rank) {
      /* Allocate the memory for communication buffer and post the receive. */
      sendBuf[i].resize(entitiesSend[i].size());

      SU2_MPI::Irecv(sendBuf[i].data(), sendBuf[i].size(), MPI_DOUBLE, ranksSend[i], rank + 10, SU2_MPI::GetComm(),
                     &sendRequests[nSendRequests++]);
    }
  }

#endif

  /*--- Carry out the self communication, if needed. Note that this part
        is outside the MPI part. First determine if self communication
        takes place at all. ---*/
  bool selfComm = false;
  unsigned long selfRecvInd, selfSendInd = 0;
  for (selfRecvInd = 0; selfRecvInd < ranksRecv.size(); ++selfRecvInd) {
    if (ranksRecv[selfRecvInd] == rank) {
      selfComm = true;
      break;
    }
  }

  /*--- Determine the send index for self communication and update the
        data of the owners with the halo data, if self communication
        takes place. ---*/
  if (selfComm) {
    for (selfSendInd = 0; selfSendInd < ranksSend.size(); ++selfSendInd) {
      if (ranksSend[selfSendInd] == rank) break;
    }

    for (unsigned long j = 0; j < entitiesSend[selfSendInd].size(); ++j) {
      const unsigned long v0 = entitiesSend[selfSendInd][j];
      const unsigned long v1 = entitiesRecv[selfRecvInd][j];

      volElem[v0].lenScale = max(volElem[v0].lenScale, volElem[v1].lenScale);
    }
  }

#ifdef HAVE_MPI

  /*--- Complete the nonblocking communication. ---*/
  SU2_MPI::Waitall(nRecvRequests, recvRequests.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Waitall(nSendRequests, sendRequests.data(), MPI_STATUSES_IGNORE);

  /*--- Update the variable currently stored in lenScale for the
        owned elements. ---*/
  for (unsigned long i = 0; i < ranksSend.size(); ++i) {
    if (ranksSend[i] != rank) {
      for (unsigned long j = 0; j < entitiesSend[i].size(); ++j) {
        const unsigned long ind = entitiesSend[i][j];
        volElem[ind].lenScale = max(volElem[ind].lenScale, sendBuf[i][j]);
      }
    }
  }

#endif

  /*-------------------------------------------------------------------*/
  /*--- Step 3: Determine the minimum value of the Jacobians of the ---*/
  /*---         elements and from that the length scale of the      ---*/
  /*---         owned elements.                                     ---*/
  /*-------------------------------------------------------------------*/

  /* Determine the number of metric terms per integration point.
     This depends on the number of spatial dimensions of the problem. */
  const unsigned short nMetricPerPoint = nDim * nDim + 1;

  /* Loop over the owned volume elements. */
  for (unsigned long i = 0; i < nVolElemOwned; ++i) {
    /* Easier storage of the metric terms and determine the number of
       integration points for this element. */
    const su2double* metric = volElem[i].metricTerms.data();
    const unsigned short ind = volElem[i].indStandardElement;
    const unsigned short nInt = standardElementsSol[ind].GetNIntegration();

    /* Loop over the integration points and determine the minimum Jacobian
       for this element. Note that the Jacobian is the first variable stored
       in the metric terms of the integration points. */
    su2double minJacElem = metric[0];
    for (unsigned short k = 1; k < nInt; ++k) minJacElem = min(minJacElem, metric[k * nMetricPerPoint]);

    /* Determine the length scale of the element, for which the length
       scale of the reference element, 2.0, must be taken into account. */
    const su2double maxJacFaces = volElem[i].lenScale;
    volElem[i].lenScale = 2.0 * minJacElem / maxJacFaces;
  }

  /*-------------------------------------------------------------------*/
  /*--- Step 4: Set the length scale for the halo elements.         ---*/
  /*-------------------------------------------------------------------*/

#ifdef HAVE_MPI

  /*--- Send the data using non-blocking sends. Exclude self communication. ---*/
  nSendRequests = 0;
  for (unsigned long i = 0; i < ranksSend.size(); ++i) {
    if (ranksSend[i] != rank) {
      /* Fill the send buffer. */
      for (unsigned long j = 0; j < entitiesSend[i].size(); ++j) sendBuf[i][j] = volElem[entitiesSend[i][j]].lenScale;

      /* Send the data. */
      SU2_MPI::Isend(sendBuf[i].data(), sendBuf[i].size(), MPI_DOUBLE, ranksSend[i], ranksSend[i] + 20,
                     SU2_MPI::GetComm(), &sendRequests[nSendRequests++]);
    }
  }

  /*--- Post the receives using non-blocking receives.
        Exclude self communication. ---*/
  nRecvRequests = 0;
  for (unsigned long i = 0; i < ranksRecv.size(); ++i) {
    if (ranksRecv[i] != rank) {
      /* Post the nonblocking receive. */
      SU2_MPI::Irecv(recvBuf[i].data(), recvBuf[i].size(), MPI_DOUBLE, ranksRecv[i], rank + 20, SU2_MPI::GetComm(),
                     &recvRequests[nRecvRequests++]);
    }
  }

#endif

  /*--- Carry out the self communication, if needed. Note that this part
        is outside the MPI part. ---*/
  if (selfComm) {
    for (unsigned long j = 0; j < entitiesSend[selfSendInd].size(); ++j) {
      const unsigned long v0 = entitiesSend[selfSendInd][j];
      const unsigned long v1 = entitiesRecv[selfRecvInd][j];
      volElem[v1].lenScale = volElem[v0].lenScale;
    }
  }

#ifdef HAVE_MPI

  /*--- Complete the nonblocking communication. ---*/
  SU2_MPI::Waitall(nSendRequests, sendRequests.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Waitall(nRecvRequests, recvRequests.data(), MPI_STATUSES_IGNORE);

  /*--- Copy the data from the receive buffers into the volElem. ---*/
  for (unsigned long i = 0; i < ranksRecv.size(); ++i) {
    if (ranksRecv[i] != rank) {
      for (unsigned long j = 0; j < entitiesRecv[i].size(); ++j) volElem[entitiesRecv[i][j]].lenScale = recvBuf[i][j];
    }
  }

#endif
}

void CMeshFEM_DG::MetricTermsSurfaceElements(CConfig* config) {
  /*--- Compute the metric terms of the internal matching faces. ---*/
  MetricTermsMatchingFaces(config);

  /*--- Loop over the physical boundaries and compute the metric
        terms of the boundary. ---*/
  for (unsigned short i = 0; i < boundaries.size(); ++i) {
    if (!boundaries[i].periodicBoundary) MetricTermsBoundaryFaces(&boundaries[i], config);
  }
}

void CMeshFEM_DG::MetricTermsVolumeElements(CConfig* config) {
  /* Find out whether or not the full, inverse of the full and the lumped mass
     matrices are needed. In principle always the full mass matrix is used,
     unless the user explicitly indices to use the lumped mass matrix for
     steady and dual time stepping. Note that for the DG_FEM the mass matrix
     is local to the elements. Finally, when the non-aliased predictor is used
     for ADER-DG for the Navier-Stokes equations, the derivatives of the
     metric terms are needed for the computation of the 2nd derivatives. */
  const bool UseLumpedMassMatrix = config->GetUse_Lumped_MassMatrix_DGFEM();

  bool FullMassMatrix = false, FullInverseMassMatrix = false;
  bool LumpedMassMatrix = false, DerMetricTerms = false;

  if (config->GetTime_Marching() == TIME_MARCHING::STEADY ||
      config->GetTime_Marching() == TIME_MARCHING::ROTATIONAL_FRAME) {
    if (UseLumpedMassMatrix)
      LumpedMassMatrix = true;
    else
      FullInverseMassMatrix = true;
  } else if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST ||
             config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND ||
             config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE) {
    if (UseLumpedMassMatrix)
      FullMassMatrix = LumpedMassMatrix = true;
    else
      FullInverseMassMatrix = true;
  } else {
    /* Time accurate explicit time integration scheme. */
    FullMassMatrix = LumpedMassMatrix = false;
    FullInverseMassMatrix = true;

    /* For ADER-DG, check if the derivative of the metric terms are needed. */
    if (config->GetKind_TimeIntScheme_Flow() == ADER_DG) {
      MAIN_SOLVER solver = config->GetKind_Solver();
      if (solver == MAIN_SOLVER::FEM_NAVIER_STOKES || solver == MAIN_SOLVER::FEM_RANS ||
          solver == MAIN_SOLVER::FEM_LES) {
        if (config->GetKind_ADER_Predictor() == ADER_NON_ALIASED_PREDICTOR) DerMetricTerms = true;
      }
    }
  }

  /* Determine the number of metric terms per integration point.
     This depends on the number of spatial dimensions of the problem.
     Also determine the additional number of metric terms per integration
     point for the computation of the second derivatives. */
  const unsigned short nMetricPerPoint = nDim * nDim + 1;
  const unsigned short nMetric2ndDerPerPoint = nDim * (nDim + nDim * (nDim - 1) / 2);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Determine the metric terms, drdx, drdy, drdz, dsdx, etc.   ---*/
  /*---         and the Jacobian in the integration points and solution    ---*/
  /*---         DOFs of the owned volume elements.                         ---*/
  /*--------------------------------------------------------------------------*/

  /* Loop over the owned volume elements. */
  for (unsigned long i = 0; i < nVolElemOwned; ++i) {
    /* Easier storage of the index of the corresponding standard element
       and determine the number of integration points as well as the number
       of DOFs for the grid and solution. */
    const unsigned short ind = volElem[i].indStandardElement;
    const unsigned short nInt = standardElementsGrid[ind].GetNIntegration();
    const unsigned short nDOFsGrid = volElem[i].nDOFsGrid;
    const unsigned short nDOFsSol = volElem[i].nDOFsSol;

    /* Allocate the memory for the grid velocities. They are initialized to
       zero as default value, i.e. not moving. */
    volElem[i].gridVelocities.assign(nDim * nInt, 0.0);
    volElem[i].gridVelocitiesSolDOFs.assign(nDim * nDOFsSol, 0.0);

    /* Allocate the memory for the metric terms of this element. */
    volElem[i].metricTerms.resize(nMetricPerPoint * nInt);
    volElem[i].metricTermsSolDOFs.resize(nMetricPerPoint * nDOFsSol);

    /* Get the pointer to the matrix storage of the basis functions
       and its derivatives. The first nDOFsGrid*nInt entries of this matrix
       correspond to the interpolation data to the integration points
       and are not needed. Hence this part is skipped. */
    const su2double* matBasisInt = standardElementsGrid[ind].GetMatBasisFunctionsIntegration();
    const su2double* matDerBasisInt = matBasisInt + nDOFsGrid * nInt;

    /* Also get the pointer to the matrix storage of the derivatives of the
       basis functions in the sol DOFs. */
    const su2double* matDerBasisSolDOFs = standardElementsGrid[ind].GetMatDerBasisFunctionsSolDOFs();

    /* Allocate the memory for the result vectors. */
    vector<su2double> helpVecResultInt(nInt * nDim * nDim);
    vector<su2double> helpVecResultDOFsSol(nDOFsSol * nDim * nDim);
    su2double* vecResultInt = helpVecResultInt.data();
    su2double* vecResultDOFsSol = helpVecResultDOFsSol.data();

    /* Compute the gradient of the coordinates w.r.t. the parametric
       coordinates for this element in the integration points. */
    ComputeGradientsCoorWRTParam(nInt, nDOFsGrid, matDerBasisInt, volElem[i].nodeIDsGrid.data(), vecResultInt, config);

    /* Compute the gradient of the coordinates w.r.t. the parametric
       coordinates for this element in the solution DOFs. */
    ComputeGradientsCoorWRTParam(nDOFsSol, nDOFsGrid, matDerBasisSolDOFs, volElem[i].nodeIDsGrid.data(),
                                 vecResultDOFsSol, config);

    /* Convert the values of dxdr, dydr, etc. to the required metric terms
       for both the integration points and the solution DOFs. */
    VolumeMetricTermsFromCoorGradients(nInt, vecResultInt, volElem[i].metricTerms);

    VolumeMetricTermsFromCoorGradients(nDOFsSol, vecResultDOFsSol, volElem[i].metricTermsSolDOFs);

    /* Check for negative Jacobians in the integrations points and at the
       location of the solution DOFs. */
    bool negJacobian = false;
    for (unsigned short j = 0; j < nInt; ++j) {
      if (volElem[i].metricTerms[nMetricPerPoint * j] <= 0.0) negJacobian = true;
    }

    for (unsigned short j = 0; j < nDOFsSol; ++j) {
      if (volElem[i].metricTermsSolDOFs[nMetricPerPoint * j] <= 0.0) negJacobian = true;
    }

    if (negJacobian) SU2_MPI::Error("Negative Jacobian found.", CURRENT_FUNCTION);
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Determine the derivatives of the metric terms in the       ---*/
  /*---         integration points, if needed.                             ---*/
  /*--------------------------------------------------------------------------*/

  if (DerMetricTerms) {
    /* Loop over the owned volume elements. */
    for (unsigned long i = 0; i < nVolElemOwned; ++i) {
      /* Easier storage of the index of the corresponding standard element
         and determine the number of integration points as well as the number
         of DOFs for the grid. */
      const unsigned short ind = volElem[i].indStandardElement;
      const unsigned short nInt = standardElementsGrid[ind].GetNIntegration();
      const unsigned short nDOFsGrid = volElem[i].nDOFsGrid;

      /* Get the pointer to the matrix storage of the basis functions
         and its derivatives. The first nDOFsGrid*nInt entries of this matrix
         correspond to the interpolation data to the integration points
         and are not needed. Hence this part is skipped. */
      const su2double* matBasisInt = standardElementsGrid[ind].GetMatBasisFunctionsIntegration();
      const su2double* matDerBasisInt = matBasisInt + nDOFsGrid * nInt;

      /* Also get the pointer to the matrix storage of the derivatives of the
         basis functions in the owned DOFs, i.e. the grid DOFs. */
      const su2double* matDerBasisGridDOFs = standardElementsGrid[ind].GetMatDerBasisFunctionsOwnDOFs();

      /* Compute the gradient of the coordinates w.r.t. the parametric
         coordinates for this element in the solution DOFs. */
      vector<su2double> helpVecResultGridDOFs(nDOFsGrid * nDim * nDim);
      su2double* vecResultGridDOFs = helpVecResultGridDOFs.data();

      ComputeGradientsCoorWRTParam(nDOFsGrid, nDOFsGrid, matDerBasisGridDOFs, volElem[i].nodeIDsGrid.data(),
                                   vecResultGridDOFs, config);

      /* Convert the values of dxdr, dydr, etc. to the required metric terms
         in the grid DOFs. */
      vector<su2double> metricGridDOFs(nDOFsGrid * nMetricPerPoint);
      VolumeMetricTermsFromCoorGradients(nDOFsGrid, vecResultGridDOFs, metricGridDOFs);

      /*--- The metric terms currently stored in metricGridDOFs are scaled
            with the Jacobian and also the Jacobian is part of the metric
            terms. For the derivatives of the metric terms, the Jacobian
            is not needed, but the original unscaled terms are needed.
            This is done in the loop below. ---*/
      for (unsigned short j = 0; j < nDOFsGrid; ++j) {
        su2double* metOld = metricGridDOFs.data() + j * nMetricPerPoint;
        su2double* metNew = metricGridDOFs.data() + j * (nMetricPerPoint - 1);

        const su2double JacInv = 1.0 / metOld[0];
        for (unsigned short k = 1; k < nMetricPerPoint; ++k) metNew[k - 1] = JacInv * metOld[k];
      }

      /* Compute the derivatives of the metric terms w.r.t. the
         parametric coordinates in the integration points.
         First allocate the memory to store the results. */
      vector<su2double> helpVecDerMetrics(nDim * nInt * (nMetricPerPoint - 1));
      su2double* vecDerMetrics = helpVecDerMetrics.data();

      /* Carry out the matrix multiplication. The last argument is NULL, such
         that this gemm call is ignored in the profiling. Replace by config if
         it should be included. */
      blasFunctions->gemm(nDim * nInt, nMetricPerPoint - 1, nDOFsGrid, matDerBasisInt, metricGridDOFs.data(),
                          vecDerMetrics, nullptr);

      /* Allocate the memory for the additional metric terms needed to
         compute the second derivatives. */
      volElem[i].metricTerms2ndDer.resize(nInt * nMetric2ndDerPerPoint);

      /*--- Loop over the integration points to compute the additional metric
            terms needed for the second derivatives. This is a combination of
            the original metric terms and derivatives of these terms.
            Make a distinction between 2D and 3D. ---*/
      switch (nDim) {
        case 2: {
          /* 2D computation. Loop over the integration points. */
          for (unsigned short j = 0; j < nInt; ++j) {
            /* Set the pointers where the data for this integration
               point starts. */
            const su2double* metric = volElem[i].metricTerms.data() + j * nMetricPerPoint;
            const su2double* rDerMetric = vecDerMetrics + j * (nMetricPerPoint - 1);
            const su2double* sDerMetric = rDerMetric + nInt * (nMetricPerPoint - 1);
            su2double* metric2ndDer = volElem[i].metricTerms2ndDer.data() + j * nMetric2ndDerPerPoint;

            /* More readable abbreviations for the metric terms
               and its derivatives. */
            const su2double JacInv = 1.0 / metric[0];
            const su2double drdx = JacInv * metric[1], drdy = JacInv * metric[2];
            const su2double dsdx = JacInv * metric[3], dsdy = JacInv * metric[4];

            const su2double ddrdx_dr = rDerMetric[0], ddrdy_dr = rDerMetric[1];
            const su2double ddsdx_dr = rDerMetric[2], ddsdy_dr = rDerMetric[3];

            const su2double ddrdx_ds = sDerMetric[0], ddrdy_ds = sDerMetric[1];
            const su2double ddsdx_ds = sDerMetric[2], ddsdy_ds = sDerMetric[3];

            /* Compute the metric terms needed to compute the Cartesian 2nd
               derivatives. Note that for the cross derivatives an average of
               the two possibilities is taken. */
            metric2ndDer[0] = drdx * ddrdx_dr + dsdx * ddrdx_ds;
            metric2ndDer[1] = drdx * ddsdx_dr + dsdx * ddsdx_ds;

            metric2ndDer[2] = 0.5 * (drdx * ddrdy_dr + dsdx * ddrdy_ds + drdy * ddrdx_dr + dsdy * ddrdx_ds);
            metric2ndDer[3] = 0.5 * (drdx * ddsdy_dr + dsdx * ddsdy_ds + drdy * ddsdx_dr + dsdy * ddsdx_ds);

            metric2ndDer[4] = drdy * ddrdy_dr + dsdy * ddrdy_ds;
            metric2ndDer[5] = drdy * ddsdy_dr + dsdy * ddsdy_ds;
          }

          break;
        }

        case 3: {
          /* 3D computation. Loop over the integration points. */
          for (unsigned short j = 0; j < nInt; ++j) {
            /* Set the pointers where the data for this integration
               point starts. */
            const su2double* metric = volElem[i].metricTerms.data() + j * nMetricPerPoint;
            const su2double* rDerMetric = vecDerMetrics + j * (nMetricPerPoint - 1);
            const su2double* sDerMetric = rDerMetric + nInt * (nMetricPerPoint - 1);
            const su2double* tDerMetric = sDerMetric + nInt * (nMetricPerPoint - 1);
            su2double* metric2ndDer = volElem[i].metricTerms2ndDer.data() + j * nMetric2ndDerPerPoint;

            /* More readable abbreviations for the metric terms
               and its derivatives. */
            const su2double JacInv = 1.0 / metric[0];
            const su2double drdx = JacInv * metric[1], drdy = JacInv * metric[2], drdz = JacInv * metric[3];
            const su2double dsdx = JacInv * metric[4], dsdy = JacInv * metric[5], dsdz = JacInv * metric[6];
            const su2double dtdx = JacInv * metric[7], dtdy = JacInv * metric[8], dtdz = JacInv * metric[9];

            const su2double ddrdx_dr = rDerMetric[0], ddrdy_dr = rDerMetric[1], ddrdz_dr = rDerMetric[2];
            const su2double ddsdx_dr = rDerMetric[3], ddsdy_dr = rDerMetric[4], ddsdz_dr = rDerMetric[5];
            const su2double ddtdx_dr = rDerMetric[6], ddtdy_dr = rDerMetric[7], ddtdz_dr = rDerMetric[8];

            const su2double ddrdx_ds = sDerMetric[0], ddrdy_ds = sDerMetric[1], ddrdz_ds = sDerMetric[2];
            const su2double ddsdx_ds = sDerMetric[3], ddsdy_ds = sDerMetric[4], ddsdz_ds = sDerMetric[5];
            const su2double ddtdx_ds = sDerMetric[6], ddtdy_ds = sDerMetric[7], ddtdz_ds = sDerMetric[8];

            const su2double ddrdx_dt = tDerMetric[0], ddrdy_dt = tDerMetric[1], ddrdz_dt = tDerMetric[2];
            const su2double ddsdx_dt = tDerMetric[3], ddsdy_dt = tDerMetric[4], ddsdz_dt = tDerMetric[5];
            const su2double ddtdx_dt = tDerMetric[6], ddtdy_dt = tDerMetric[7], ddtdz_dt = tDerMetric[8];

            /* Compute the metric terms needed to compute the Cartesian 2nd
               derivatives. Note that for the cross derivatives an average of
               the two possibilities is taken. */
            metric2ndDer[0] = drdx * ddrdx_dr + dsdx * ddrdx_ds + dtdx * ddrdx_dt;
            metric2ndDer[1] = drdx * ddsdx_dr + dsdx * ddsdx_ds + dtdx * ddsdx_dt;
            metric2ndDer[2] = drdx * ddtdx_dr + dsdx * ddtdx_ds + dtdx * ddtdx_dt;

            metric2ndDer[3] = 0.5 * (drdx * ddrdy_dr + dsdx * ddrdy_ds + dtdx * ddrdy_dt + drdy * ddrdx_dr +
                                     dsdy * ddrdx_ds + dtdy * ddrdx_dt);
            metric2ndDer[4] = 0.5 * (drdx * ddsdy_dr + dsdx * ddsdy_ds + dtdx * ddsdy_dt + drdy * ddsdx_dr +
                                     dsdy * ddsdx_ds + dtdy * ddsdx_dt);
            metric2ndDer[5] = 0.5 * (drdx * ddtdy_dr + dsdx * ddtdy_ds + dtdx * ddtdy_dt + drdy * ddtdx_dr +
                                     dsdy * ddtdx_ds + dtdy * ddtdx_dt);

            metric2ndDer[6] = drdy * ddrdy_dr + dsdy * ddrdy_ds + dtdy * ddrdy_dt;
            metric2ndDer[7] = drdy * ddsdy_dr + dsdy * ddsdy_ds + dtdy * ddsdy_dt;
            metric2ndDer[8] = drdy * ddtdy_dr + dsdy * ddtdy_ds + dtdy * ddtdy_dt;

            metric2ndDer[9] = 0.5 * (drdx * ddrdz_dr + dsdx * ddrdz_ds + dtdx * ddrdz_dt + drdz * ddrdx_dr +
                                     dsdz * ddrdx_ds + dtdz * ddrdx_dt);
            metric2ndDer[10] = 0.5 * (drdx * ddsdz_dr + dsdx * ddsdz_ds + dtdx * ddsdz_dt + drdz * ddsdx_dr +
                                      dsdz * ddsdx_ds + dtdz * ddsdx_dt);
            metric2ndDer[11] = 0.5 * (drdx * ddtdz_dr + dsdx * ddtdz_ds + dtdx * ddtdz_dt + drdz * ddtdx_dr +
                                      dsdz * ddtdx_ds + dtdz * ddtdx_dt);

            metric2ndDer[12] = 0.5 * (drdy * ddrdz_dr + dsdy * ddrdz_ds + dtdy * ddrdz_dt + drdz * ddrdy_dr +
                                      dsdz * ddrdy_ds + dtdz * ddrdy_dt);
            metric2ndDer[13] = 0.5 * (drdy * ddsdz_dr + dsdy * ddsdz_ds + dtdy * ddsdz_dt + drdz * ddsdy_dr +
                                      dsdz * ddsdy_ds + dtdz * ddsdy_dt);
            metric2ndDer[14] = 0.5 * (drdy * ddtdz_dr + dsdy * ddtdz_ds + dtdy * ddtdz_dt + drdz * ddtdy_dr +
                                      dsdz * ddtdy_ds + dtdz * ddtdy_dt);

            metric2ndDer[15] = drdz * ddrdz_dr + dsdz * ddrdz_ds + dtdz * ddrdz_dt;
            metric2ndDer[16] = drdz * ddsdz_dr + dsdz * ddsdz_ds + dtdz * ddsdz_dt;
            metric2ndDer[17] = drdz * ddtdz_dr + dsdz * ddtdz_ds + dtdz * ddtdz_dt;
          }

          break;
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Determine the mass matrix (or its inverse) and/or the      ---*/
  /*---         lumped mass matrix.                                        ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the time coefficients in the iteration matrix of the ADER-DG
     predictor step. */
  if (config->GetKind_TimeIntScheme_Flow() == ADER_DG) TimeCoefficientsPredictorADER_DG(config);

  /* Define the double vector to define the values of the mass matrix
     contributions of the standard element in the integration points. */
  vector<vector<su2double> > valMMInt;
  valMMInt.resize(standardElementsSol.size());

  /* Loop over the different standard elements. */
  for (unsigned long i = 0; i < standardElementsSol.size(); ++i) {
    /* Get the required data from the standard element. */
    const unsigned short nInt = standardElementsSol[i].GetNIntegration();
    const unsigned short nDOFs = standardElementsSol[i].GetNDOFs();
    const su2double* lag = standardElementsSol[i].GetBasisFunctionsIntegration();
    const su2double* w = standardElementsSol[i].GetWeightsIntegration();

    /* Allocate the memory for valMMInt[i]. */
    valMMInt[i].resize(nDOFs * nDOFs * nInt);

    /* Create the values of valMMInt[i], which is a matrix of size
       nDOFs2XnInt stored in row major order. */
    unsigned int ll = 0;
    for (unsigned short k = 0; k < nDOFs; ++k)
      for (unsigned short j = 0; j < nDOFs; ++j)
        for (unsigned short l = 0; l < nInt; ++l, ++ll)
          valMMInt[i][ll] = w[l] * lag[l * nDOFs + k] * lag[l * nDOFs + j];
  }

  /* Loop over the owned volume elements. */
  for (unsigned long i = 0; i < nVolElemOwned; ++i) {
    /* Easier storage of the index of the corresponding standard element and
       get the necessary information from it. */
    const unsigned short ind = volElem[i].indStandardElement;
    const unsigned short nInt = standardElementsSol[ind].GetNIntegration();
    const unsigned short nDOFs = volElem[i].nDOFsSol;
    const unsigned int nDOFs2 = nDOFs * nDOFs;
    const su2double* valInt = valMMInt[ind].data();

    /* Copy the data of the Jacobian in the integration points in
       a separate vector. */
    vector<su2double> JacVec(nInt);
    su2double* Jac = JacVec.data();

    for (unsigned short l = 0; l < nInt; ++l) Jac[l] = volElem[i].metricTerms[l * nMetricPerPoint];

    /* Allocate the memory for working vector for the construction
       of the mass matrix. */
    vector<su2double> massMat(nDOFs2);

    /* Create the mass matrix, which is a BLAS gemv operation. */
    blasFunctions->gemv(nDOFs2, nInt, valInt, Jac, massMat.data());

    /* Store the full mass matrix in volElem[i], if needed. */
    if (FullMassMatrix) volElem[i].massMatrix = massMat;

    /*--- Check if the lumped mass matrix is needed. ---*/
    if (LumpedMassMatrix) {
      /* Allocate the memory for the lumped mass matrix and initialize them to zero.. */
      volElem[i].lumpedMassMatrix.assign(nDOFs, 0.0);

      /* Loop over the DOFs to compute the elements of the local lumped mass
         matrix. It is the sum of the absolute values of the row. */
      for (unsigned short j = 0; j < nDOFs; ++j) {
        const su2double* val = massMat.data() + j * nDOFs;
        for (unsigned short k = 0; k < nDOFs; ++k) volElem[i].lumpedMassMatrix[j] += fabs(val[k]);
      }
    }

    /*--- Check if the inverse of mass matrix is needed. ---*/
    if (FullInverseMassMatrix) {
      /*--- Check if LAPACK/MKL can be used to compute the inverse. ---*/

#if (defined(HAVE_MKL) || defined(HAVE_LAPACK)) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))

      /* The inverse can be computed using the Lapack routines dpotrf
         and dpotri. As the mass matrix is positive definite, a
         Cholesky decomposition is used, which is much more efficient
         than a standard inverse. */
      char uplo = 'L';
      int NN = nDOFs, errorCode;
      dpotrf_(&uplo, &NN, massMat.data(), &NN, &errorCode);
      if (errorCode != 0) {
        ostringstream message;
        if (errorCode < 0) {
          message << "Something wrong when calling dpotrf. Error code: " << errorCode;
        } else {
          message << "Mass matrix not positive definite. " << endl;
          message << "This is most likely caused by a too low accuracy of the quadrature rule," << endl;
          message << "possibly combined with a low quality element." << endl;
          message << "Increase the accuracy of the quadrature rule.";
        }

        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      dpotri_(&uplo, &NN, massMat.data(), &NN, &errorCode);
      if (errorCode != 0) {
        ostringstream message;
        if (errorCode < 0) {
          message << "Something wrong when calling dpotri. Error code: " << errorCode;
        } else {
          message << "Mass matrix is singular. " << endl;
          message << "The is most likely caused by a too low accuracy of the quadrature rule, " << endl;
          message << "possibly combined with a low quality element." << endl;
          message << "Increase the accuracy of the quadrature rule.";
        }

        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /* The Lapack routines for a Cholesky decomposition only store the upper
         part of the matrix. Copy the data to the lower part. */
      for (unsigned short k = 0; k < nDOFs; ++k) {
        for (unsigned short j = (k + 1); j < nDOFs; ++j) {
          massMat[j * nDOFs + k] = massMat[k * nDOFs + j];
        }
      }
#else
      /* No support for Lapack. Hence an internal routine is used.
         This does not do all the checking the Lapack routines do. */
      CFEMStandardElementBase::InverseMatrix(nDOFs, massMat);
#endif
      /* Store the inverse of the mass matrix in volElem[i]. */
      volElem[i].invMassMatrix = massMat;
    }
  }
}

void CMeshFEM_DG::TimeCoefficientsPredictorADER_DG(CConfig* config) {
  /*--------------------------------------------------------------------------*/
  /*--- Determine the coefficients that appear in the iteration matrix of  ---*/
  /*--- the prediction step of ADER-DG.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the number of time DOFs in the predictor step of ADER as well
     as their location on the interval [-1..1]. */
  const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();
  const su2double* TimeDOFs = config->GetTimeDOFsADER_DG();

  /* Compute the Vandermonde matrix and its inverse in the time DOFs. */
  vector<su2double> rTimeDOFs(nTimeDOFs);
  for (unsigned short i = 0; i < nTimeDOFs; ++i) rTimeDOFs[i] = TimeDOFs[i];

  vector<su2double> V(nTimeDOFs * nTimeDOFs);
  CFEMStandardElementBase timeElement;
  timeElement.Vandermonde1D(nTimeDOFs, rTimeDOFs, V);

  vector<su2double> VInv = V;
  timeElement.InverseMatrix(nTimeDOFs, VInv);

  /* Compute the Vandermonde matrix for r = 1, i.e. the end of the interval. */
  vector<su2double> rEnd(1);
  rEnd[0] = 1.0;
  vector<su2double> VEnd(nTimeDOFs);
  timeElement.Vandermonde1D(nTimeDOFs, rEnd, VEnd);

  /*--- Determine the matrix products VEnd*VInv to get the correct expression
        for the values of the Lagrangian interpolation functions for r = 1. ---*/
  vector<su2double> lEnd(nTimeDOFs, 0.0);
  for (unsigned short j = 0; j < nTimeDOFs; ++j) {
    for (unsigned short k = 0; k < nTimeDOFs; ++k) lEnd[j] += VEnd[k] * VInv[j * nTimeDOFs + k];
  }

  /*--- To reduce the error due to round off, make sure that the row sum is 1.
        Also check if the difference is not too large to be solely caused by
        roundoff.  ---*/
  su2double val = 0.0;
  for (unsigned short j = 0; j < nTimeDOFs; ++j) val += lEnd[j];

  if (fabs(val - 1.0) > 1.e-6) SU2_MPI::Error("Difference is too large to be caused by roundoff", CURRENT_FUNCTION);

  val = 1.0 / val;
  for (unsigned short j = 0; j < nTimeDOFs; ++j) lEnd[j] *= val;

  /*--- Compute the mass matrix in time, which is the inverse of V^T V. Note that
        this definition is different from the one used in Hesthaven, because their
        V is our V^T. Furthermore note that the mass matrix should be a diagonal
        matrix, due to the choice of the DOFs (Gauss Legendre points). ---*/
  vector<su2double> MassTime(nTimeDOFs * nTimeDOFs, 0.0);
  for (unsigned short j = 0; j < nTimeDOFs; ++j) {
    for (unsigned short i = 0; i < nTimeDOFs; ++i) {
      const unsigned short ji = j * nTimeDOFs + i;
      for (unsigned short k = 0; k < nTimeDOFs; ++k) MassTime[ji] += V[k * nTimeDOFs + j] * V[k * nTimeDOFs + i];
    }
  }

  timeElement.InverseMatrix(nTimeDOFs, MassTime);

  /* Compute the gradient of the Vandermonde matrix in the time DOFs. */
  vector<su2double> VDr(nTimeDOFs * nTimeDOFs);
  timeElement.GradVandermonde1D(nTimeDOFs, rTimeDOFs, VDr);

  /* Compute the product VDr VInv. Store the result in V. Note that in both
     matrices the transpose is stored compared to the definition of Hesthaven. */
  for (unsigned short j = 0; j < nTimeDOFs; ++j) {
    for (unsigned short i = 0; i < nTimeDOFs; ++i) {
      const unsigned short ji = j * nTimeDOFs + i;
      V[ji] = 0.0;
      for (unsigned short k = 0; k < nTimeDOFs; ++k) V[ji] += VDr[k * nTimeDOFs + j] * VInv[i * nTimeDOFs + k];
    }
  }

  /* Compute the time stiffness matrix S = M Dr. */
  vector<su2double> S(nTimeDOFs * nTimeDOFs, 0.0);
  for (unsigned short j = 0; j < nTimeDOFs; ++j) {
    for (unsigned short i = 0; i < nTimeDOFs; ++i) {
      const unsigned short ji = j * nTimeDOFs + i;
      for (unsigned short k = 0; k < nTimeDOFs; ++k) S[ji] += MassTime[j * nTimeDOFs + k] * V[k * nTimeDOFs + i];
    }
  }

  /* Compute the time coefficients of the iteration matrix in the predictor step. */
  timeCoefADER_DG.assign(nTimeDOFs * nTimeDOFs, 0.0);
  for (unsigned short j = 0; j < nTimeDOFs; ++j) {
    for (unsigned short i = 0; i < nTimeDOFs; ++i)
      timeCoefADER_DG[j * nTimeDOFs + i] = lEnd[j] * lEnd[i] - S[i * nTimeDOFs + j];
  }

  /* Determine the inverse of timeCoefADER_DG, because this is needed in the
     iteration matrix of the predictor step of ADER-DG. */
  timeElement.InverseMatrix(nTimeDOFs, timeCoefADER_DG);

  /*--------------------------------------------------------------------------*/
  /*--- Determine the interpolation matrix from the time DOFs of ADER-DG   ---*/
  /*--- to the time integration points. Note that in many cases this       ---*/
  /*--- interpolation matrix is the identity matrix, i.e. the DOFs and the ---*/
  /*--- integration points coincide.                                       ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the number of time integration points of ADER as well
     as their location on the interval [-1..1]. */
  unsigned short nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();
  const su2double* TimeIntegrationPoints = config->GetTimeIntegrationADER_DG();

  /* Compute the Vandermonde matrix for the time integration points. */
  vector<su2double> rTimeIntPoints(nTimeIntegrationPoints);
  for (unsigned short i = 0; i < nTimeIntegrationPoints; ++i) rTimeIntPoints[i] = TimeIntegrationPoints[i];

  V.resize(nTimeIntegrationPoints * nTimeDOFs);
  timeElement.Vandermonde1D(nTimeDOFs, rTimeIntPoints, V);

  /* Allocate the memory for the time interpolation coefficients and initialize
     them to zero. */
  timeInterpolDOFToIntegrationADER_DG.assign(nTimeIntegrationPoints * nTimeDOFs, 0.0);

  /* Determine the matrix products V*VInv to get the correct expression for
     mTimeInterpolDOFToIntegration. Note that from mathematical point of view the
     transpose of V*VInv is stored, because in this way the interpolation data
     for a time integration point is contiguous in memory. */
  for (unsigned short j = 0; j < nTimeDOFs; ++j) {
    for (unsigned short i = 0; i < nTimeIntegrationPoints; ++i) {
      const unsigned short ii = i * nTimeDOFs + j;

      for (unsigned short k = 0; k < nTimeDOFs; ++k) {
        const unsigned short indV = k * nTimeIntegrationPoints + i;
        const unsigned short indVInv = j * nTimeDOFs + k;

        timeInterpolDOFToIntegrationADER_DG[ii] += V[indV] * VInv[indVInv];
      }
    }
  }

  /*--- To reduce the error due to round off, make sure that the row sum is 1.
        Also check if the difference is not too large to be solely caused by
        roundoff.  ---*/
  for (unsigned short j = 0; j < nTimeIntegrationPoints; ++j) {
    const unsigned short jj = j * nTimeDOFs;

    val = 0.0;
    for (unsigned short i = 0; i < nTimeDOFs; ++i) val += timeInterpolDOFToIntegrationADER_DG[jj + i];

    if (fabs(val - 1.0) > 1.e-6) SU2_MPI::Error("Difference is too large to be caused by roundoff", CURRENT_FUNCTION);

    val = 1.0 / val;
    for (unsigned short i = 0; i < nTimeDOFs; ++i) timeInterpolDOFToIntegrationADER_DG[jj + i] *= val;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Determine the interpolation matrix from the time DOFs of ADER-DG   ---*/
  /*--- of adjacent elements of a higher time level to the time            ---*/
  /*--- integration points.                                                ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the location of the time integration points from the
     perspective of an adjacent element with twice the time step. Note that
     the time step of this adjacent element covers two fine time steps and
     hence the number of integration points double and the position is
     determined via the transformation xxi = 0.5(xi-1) and xi varies between
     -1 and 3. */
  rTimeIntPoints.resize(2 * nTimeIntegrationPoints);
  for (unsigned short i = 0; i < nTimeIntegrationPoints; ++i)
    rTimeIntPoints[i] = 0.5 * (TimeIntegrationPoints[i] - 1.0);

  for (unsigned short i = 0; i < nTimeIntegrationPoints; ++i) {
    const unsigned short ii = i + nTimeIntegrationPoints;
    rTimeIntPoints[ii] = 0.5 * (TimeIntegrationPoints[i] + 1.0);
  }

  nTimeIntegrationPoints *= 2;

  /* Compute the Vandermonde matrix for these time integration points. */
  V.resize(nTimeIntegrationPoints * nTimeDOFs);
  timeElement.Vandermonde1D(nTimeDOFs, rTimeIntPoints, V);

  /* Allocate the memory for the time interpolation coefficients and initialize
     them to zero. */
  timeInterpolAdjDOFToIntegrationADER_DG.assign(nTimeIntegrationPoints * nTimeDOFs, 0.0);

  /* Determine the matrix products V*VInv to get the correct expression for
     mTimeInterpolAdjDOFToIntegration. Note that from mathematical point of view
     the transpose of V*VInv is stored, because in this way the interpolation
     data for a time integration point is contiguous in memory. */
  for (unsigned short j = 0; j < nTimeDOFs; ++j) {
    for (unsigned short i = 0; i < nTimeIntegrationPoints; ++i) {
      const unsigned short ii = i * nTimeDOFs + j;

      for (unsigned short k = 0; k < nTimeDOFs; ++k) {
        const unsigned short indV = k * nTimeIntegrationPoints + i;
        const unsigned short indVInv = j * nTimeDOFs + k;

        timeInterpolAdjDOFToIntegrationADER_DG[ii] += V[indV] * VInv[indVInv];
      }
    }
  }

  /*--- To reduce the error due to round off, make sure that the row sum is 1.
        Also check if the difference is not too large to be solely caused by
        roundoff.  ---*/
  for (unsigned short j = 0; j < nTimeIntegrationPoints; ++j) {
    const unsigned short jj = j * nTimeDOFs;

    val = 0.0;
    for (unsigned short i = 0; i < nTimeDOFs; ++i) val += timeInterpolAdjDOFToIntegrationADER_DG[jj + i];

    if (fabs(val - 1.0) > 1.e-6) SU2_MPI::Error("Difference is too large to be caused by roundoff", CURRENT_FUNCTION);

    val = 1.0 / val;
    for (unsigned short i = 0; i < nTimeDOFs; ++i) timeInterpolAdjDOFToIntegrationADER_DG[jj + i] *= val;
  }
}

void CMeshFEM_DG::VolumeMetricTermsFromCoorGradients(const unsigned short nEntities, const su2double* gradCoor,
                                                     vector<su2double>& metricTerms) {
  /*--- Convert the dxdr, dydr, etc., stored in coorGradients, to the
        required metric terms. Make a distinction between 2D and 3D. ---*/
  switch (nDim) {
    case 2: {
      /* 2D computation. Store the offset between the r and s derivatives. */
      const unsigned short off = 2 * nEntities;

      /* Loop over the entities and store the metric terms. */
      unsigned short ii = 0;
      for (unsigned short j = 0; j < nEntities; ++j) {
        const unsigned short jx = 2 * j;
        const unsigned short jy = jx + 1;
        const su2double dxdr = gradCoor[jx], dydr = gradCoor[jy];
        const su2double dxds = gradCoor[jx + off], dyds = gradCoor[jy + off];

        metricTerms[ii++] = dxdr * dyds - dxds * dydr;  // J
        metricTerms[ii++] = dyds;                       // J drdx
        metricTerms[ii++] = -dxds;                      // J drdy
        metricTerms[ii++] = -dydr;                      // J dsdx
        metricTerms[ii++] = dxdr;                       // J dsdy
      }

      break;
    }

    case 3: {
      /* 3D computation. Store the offset between the r and s and r and t derivatives. */
      unsigned short offS = 3 * nEntities, offT = 6 * nEntities;

      /* Loop over the entities and store the metric terms. */
      unsigned short ii = 0;
      for (unsigned short j = 0; j < nEntities; ++j) {
        const unsigned short jx = 3 * j;
        const unsigned short jy = jx + 1, jz = jx + 2;
        const su2double dxdr = gradCoor[jx], dydr = gradCoor[jy], dzdr = gradCoor[jz];
        const su2double dxds = gradCoor[jx + offS], dyds = gradCoor[jy + offS], dzds = gradCoor[jz + offS];
        const su2double dxdt = gradCoor[jx + offT], dydt = gradCoor[jy + offT], dzdt = gradCoor[jz + offT];

        metricTerms[ii++] = dxdr * (dyds * dzdt - dzds * dydt) - dxds * (dydr * dzdt - dzdr * dydt) +
                            dxdt * (dydr * dzds - dzdr * dyds);  // J

        metricTerms[ii++] = dyds * dzdt - dzds * dydt;  // J drdx
        metricTerms[ii++] = dzds * dxdt - dxds * dzdt;  // J drdy
        metricTerms[ii++] = dxds * dydt - dyds * dxdt;  // J drdz

        metricTerms[ii++] = dzdr * dydt - dydr * dzdt;  // J dsdx
        metricTerms[ii++] = dxdr * dzdt - dzdr * dxdt;  // J dsdy
        metricTerms[ii++] = dydr * dxdt - dxdr * dydt;  // J dsdz

        metricTerms[ii++] = dydr * dzds - dzdr * dyds;  // J dtdx
        metricTerms[ii++] = dzdr * dxds - dxdr * dzds;  // J dtdy
        metricTerms[ii++] = dxdr * dyds - dydr * dxds;  // J dtdz
      }

      break;
    }
  }
}

void CMeshFEM_DG::WallFunctionPreprocessing(CConfig* config) {
  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Check whether wall functions are used at all.              ---*/
  /*--------------------------------------------------------------------------*/

  bool wallFunctions = false;
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX: {
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) wallFunctions = true;
        break;
      }
      default: /* Just to avoid a compiler warning. */
        break;
    }
  }

  /* If no wall functions are used, nothing needs to be done and a
     return can be made. */
  if (!wallFunctions) return;

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Build the local ADT of the volume elements. The halo       ---*/
  /*---         elements are included such that also direct neighbors are  ---*/
  /*---         guaranteed to be found as donor element. Note that the ADT ---*/
  /*---         is built with the linear subelements. This is done to      ---*/
  /*---         avoid relatively many expensive Newton solves for high     ---*/
  /*---         order elements.                                            ---*/
  /*--------------------------------------------------------------------------*/

  /* Define the vectors, which store the mapping from the subelement to the
     parent element, subelement ID within the parent element, the element
     type and the connectivity of the subelements. */
  vector<unsigned long> parentElement;
  vector<unsigned short> subElementIDInParent;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned long> elemConn;

  /* Loop over the locally stored volume elements (including halo elements)
     to create the connectivity of the subelements. */
  for (unsigned long l = 0; l < nVolElemTot; ++l) {
    /* Determine the necessary data from the corresponding standard element. */
    const unsigned short ind = volElem[l].indStandardElement;

    unsigned short VTK_Type[] = {standardElementsGrid[ind].GetVTK_Type1(), standardElementsGrid[ind].GetVTK_Type2()};
    unsigned short nSubElems[] = {0, 0};
    unsigned short nDOFsPerSubElem[] = {0, 0};
    const unsigned short* connSubElems[] = {nullptr, nullptr};

    if (VTK_Type[0] != NONE) {
      nSubElems[0] = standardElementsGrid[ind].GetNSubElemsType1();
      nDOFsPerSubElem[0] = standardElementsGrid[ind].GetNDOFsPerSubElem(VTK_Type[0]);
      connSubElems[0] = standardElementsGrid[ind].GetSubConnType1();
    }

    if (VTK_Type[1] != NONE) {
      nSubElems[1] = standardElementsGrid[ind].GetNSubElemsType2();
      nDOFsPerSubElem[1] = standardElementsGrid[ind].GetNDOFsPerSubElem(VTK_Type[1]);
      connSubElems[1] = standardElementsGrid[ind].GetSubConnType2();
    }

    /* Abbreviate the grid DOFs of this element a bit easier. */
    const unsigned long* DOFs = volElem[l].nodeIDsGrid.data();

    /* Loop over the number of subelements and store the required data. */
    unsigned short jj = 0;
    for (unsigned short i = 0; i < 2; ++i) {
      unsigned short kk = 0;
      for (unsigned short j = 0; j < nSubElems[i]; ++j, ++jj) {
        parentElement.push_back(l);
        subElementIDInParent.push_back(jj);
        VTK_TypeElem.push_back(VTK_Type[i]);

        for (unsigned short k = 0; k < nDOFsPerSubElem[i]; ++k, ++kk) elemConn.push_back(DOFs[connSubElems[i][kk]]);
      }
    }
  }

  /* Copy the coordinates in a vector that can be used by the ADT. */
  vector<su2double> volCoor;
  volCoor.reserve(nDim * meshPoints.size());

  for (unsigned long l = 0; l < meshPoints.size(); ++l) {
    for (unsigned short k = 0; k < nDim; ++k) volCoor.push_back(meshPoints[l].coor[k]);
  }

  /* Build the local ADT. */
  CADTElemClass localVolumeADT(nDim, volCoor, elemConn, VTK_TypeElem, subElementIDInParent, parentElement, false);

  /* Release the memory of the vectors used to build the ADT. To make sure
     that all the memory is deleted, the swap function is used. */
  vector<unsigned short>().swap(subElementIDInParent);
  vector<unsigned short>().swap(VTK_TypeElem);
  vector<unsigned long>().swap(parentElement);
  vector<unsigned long>().swap(elemConn);
  vector<su2double>().swap(volCoor);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3. Search for donor elements at the exchange locations in     ---*/
  /*---         the local elements.                                        ---*/
  /*--------------------------------------------------------------------------*/

  /* Loop over the markers and select the ones for which a wall function
     treatment must be carried out. */
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX: {
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) {
          /* An LES wall model is used for this boundary marker. Determine
             which wall model and allocate the memory for the member variable. */
          switch (config->GetWallFunction_Treatment(Marker_Tag)) {
            case WALL_FUNCTIONS::EQUILIBRIUM_MODEL: {
              if (rank == MASTER_NODE) cout << "Marker " << Marker_Tag << " uses an Equilibrium Wall Model." << endl;

              boundaries[iMarker].wallModel = new CWallModel1DEQ(config, Marker_Tag);
              break;
            }
            case WALL_FUNCTIONS::LOGARITHMIC_MODEL: {
              if (rank == MASTER_NODE)
                cout << "Marker " << Marker_Tag << " uses the Reichardt and Kader analytical laws for the Wall Model."
                     << endl;

              boundaries[iMarker].wallModel = new CWallModelLogLaw(config, Marker_Tag);
              break;
            }
            default: {
              SU2_MPI::Error("Wall function not present yet", CURRENT_FUNCTION);
            }
          }

          /* Retrieve the double information for this wall model. The height
             of the exchange location is the first element of this array. */
          const su2double* doubleInfo = config->GetWallFunction_DoubleInfo(Marker_Tag);

          /* Easier storage of the surface elements and loop over them. */
          vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;
          for (unsigned long l = 0; l < surfElem.size(); ++l) {
            /* Determine the corresponding standard face element and get the
               relevant information from it. Note that the standard element
               of the solution must be taken and not of the grid. */
            const unsigned short ind = surfElem[l].indStandardElement;
            const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

            /* Allocate the memory for the memory to store the donors and
               the parametric weights. The donor elements are stored in an
               CUnsignedLong2T, such that they can be sorted. */
            vector<CUnsignedLong2T> donorElements(nInt);
            vector<su2double> parCoorInDonor(nInt * nDim);

            /* Loop over the integration points. */
            for (unsigned short i = 0; i < nInt; ++i) {
              /* Easier storage of the normals and coordinates of this
                 integration point. */
              const su2double* normals = surfElem[l].metricNormalsFace.data() + i * (nDim + 1);
              const su2double* Coord = surfElem[l].coorIntegrationPoints.data() + i * nDim;

              /* Determine the coordinates of the exchange location. Note that the normals
                 point out of the domain, so the normals must be subtracted. */
              su2double coorExchange[] = {0.0, 0.0, 0.0};  // To avoid a compiler warning.
              for (unsigned short iDim = 0; iDim < nDim; ++iDim)
                coorExchange[iDim] = Coord[iDim] - doubleInfo[0] * normals[iDim];

              /* Search for the element, which contains the exchange location. */
              unsigned short subElem;
              unsigned long parElem;
              int rank;
              su2double parCoor[3], weightsInterpol[8];
              if (localVolumeADT.DetermineContainingElement(coorExchange, subElem, parElem, rank, parCoor,
                                                            weightsInterpol)) {
                /* Subelement found that contains the exchange location. However,
                   what is needed is the location in the high order parent element.
                   Determine this. */
                HighOrderContainmentSearch(coorExchange, parElem, subElem, weightsInterpol, parCoor);

                /* Store the info in donorElements and parCoorInDonor. */
                donorElements[i].long0 = parElem;
                donorElements[i].long1 = i;
                for (unsigned short iDim = 0; iDim < nDim; ++iDim) parCoorInDonor[nDim * i + iDim] = parCoor[iDim];
              } else {
                /* No subelement found that contains the exchange location.
                   The partitioning is done such that this should not happen.
                   Print an error message and exit. */
                SU2_MPI::Error("Exchange location not found in ADT", CURRENT_FUNCTION);
              }
            }

            /* Sort donorElements in increasing order, such that the
               integration points with the same donor are grouped together. */
            sort(donorElements.begin(), donorElements.end());

            /* Store the donor information in the member variables of surfElem. */
            surfElem[l].donorsWallFunction.push_back(donorElements[0].long0);
            surfElem[l].nIntPerWallFunctionDonor.push_back(0);
            surfElem[l].intPerWallFunctionDonor.resize(nInt);
            surfElem[l].intPerWallFunctionDonor[0] = donorElements[0].long1;

            for (unsigned short i = 1; i < nInt; ++i) {
              if (donorElements[i].long0 != surfElem[l].donorsWallFunction.back()) {
                surfElem[l].donorsWallFunction.push_back(donorElements[i].long0);
                surfElem[l].nIntPerWallFunctionDonor.push_back(i);
              }
              surfElem[l].intPerWallFunctionDonor[i] = donorElements[i].long1;
            }

            surfElem[l].nIntPerWallFunctionDonor.push_back(nInt);

            /* Determine whether or not halo information is needed to apply
               the boundary conditions for this surface. If needed, set
               haloInfoNeededForBC of the boundary to true. As DonorsWallFunction
               is sorted in increasing order, it suffices to check the last element. */
            if (surfElem[l].donorsWallFunction.back() >= nVolElemOwned) boundaries[iMarker].haloInfoNeededForBC = true;

            /* Allocate the memory of the first index of the interpolation
               matrices for the donordata. */
            surfElem[l].matWallFunctionDonor.resize(surfElem[l].donorsWallFunction.size());

            /* Loop over the different donors for the wall function data. */
            for (unsigned long j = 0; j < surfElem[l].donorsWallFunction.size(); ++j) {
              /* Easier storage of the donor and the vector for the
                 interpolation data for this donor. */
              const unsigned long donor = surfElem[l].donorsWallFunction[j];
              vector<su2double>& matDonor = surfElem[l].matWallFunctionDonor[j];

              /* Determine the number of DOFs in the donor element. Note that
                 the standard element of the solution must be used for this
                 purpose. */
              const unsigned short ind = volElem[donor].indStandardElement;
              const unsigned short nDOFs = standardElementsSol[ind].GetNDOFs();

              /* Allocate the memory for the vector used to compute the
                 Lagrangian interpolation functions in the parametric coordinates,
                 i.e. the interpolation weights, and for matDonor, which stores
                 all these weights. */
              const unsigned short nIntThisDonor =
                  surfElem[l].nIntPerWallFunctionDonor[j + 1] - surfElem[l].nIntPerWallFunctionDonor[j];
              vector<su2double> lagBasis(nDOFs);
              matDonor.reserve(nIntThisDonor * nDOFs);

              /* Loop over the integration points for this donor element. */
              for (unsigned short i = surfElem[l].nIntPerWallFunctionDonor[j];
                   i < surfElem[l].nIntPerWallFunctionDonor[j + 1]; ++i) {
                /* Determine the pointer to the parametric coordinates of the
                   donor point of this integration point. */
                const unsigned short ii = surfElem[l].intPerWallFunctionDonor[i];
                const su2double* parCoor = parCoorInDonor.data() + ii * nDim;

                /* Determine the interpolation weights for this point and
                   store it in matDonor. */
                standardElementsSol[ind].BasisFunctionsInPoint(parCoor, lagBasis);
                matDonor.insert(matDonor.end(), lagBasis.begin(), lagBasis.end());
              }
            }
          }
        }

        break;
      }
      default: /* Just to avoid a compiler warning. */
        break;
    }
  }
}

void CMeshFEM_DG::HighOrderContainmentSearch(const su2double* coor, const unsigned long parElem,
                                             const unsigned short subElem, const su2double* weightsSubElem,
                                             su2double* parCoor) {
  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton = 1.e-10;

  /*--------------------------------------------------------------------------*/
  /* Step 1: Create an initial guess for the parametric coordinates from      */
  /*         interpolation in the linear sub-element of the parent element.   */
  /*--------------------------------------------------------------------------*/

  /* Abbreviate the index in the standard element a bit easier. */
  const unsigned short ind = volElem[parElem].indStandardElement;

  /* Define the variables to store the number of DOFs and the connectivity
     of the sub element in which the given coordinate resides. */
  unsigned short nDOFsPerSubElem = 0;
  const unsigned short* connSubElems;

  /* Check if the sub element is of the first sub-element type. */
  const unsigned short nSubElemType1 = standardElementsGrid[ind].GetNSubElemsType1();
  if (subElem < nSubElemType1) {
    /* Determine the element type and set nDOFsPerSubElem. */
    switch (standardElementsGrid[ind].GetVTK_Type1()) {
      case TRIANGLE:
        nDOFsPerSubElem = 3;
        break;
      case QUADRILATERAL:
        nDOFsPerSubElem = 4;
        break;
      case TETRAHEDRON:
        nDOFsPerSubElem = 4;
        break;
      case PYRAMID:
        nDOFsPerSubElem = 5;
        break;
      case PRISM:
        nDOFsPerSubElem = 6;
        break;
      case HEXAHEDRON:
        nDOFsPerSubElem = 8;
        break;
      default:
        break; /* Just to avoid a compiler warning. */
    }

    /* Set the connectivity for the correct subelement. */
    connSubElems = standardElementsGrid[ind].GetSubConnType1() + subElem * nDOFsPerSubElem;
  } else {
    /* The sub-element is of the second sub-element type. Determine the
       element type and set nDOFsPerSubElem. */
    switch (standardElementsGrid[ind].GetVTK_Type2()) {
      case TRIANGLE:
        nDOFsPerSubElem = 3;
        break;
      case QUADRILATERAL:
        nDOFsPerSubElem = 4;
        break;
      case TETRAHEDRON:
        nDOFsPerSubElem = 4;
        break;
      case PYRAMID:
        nDOFsPerSubElem = 5;
        break;
      case PRISM:
        nDOFsPerSubElem = 6;
        break;
      case HEXAHEDRON:
        nDOFsPerSubElem = 8;
        break;
      default:
        break; /* Just to avoid a compiler warning. */
    }

    /* Set the connectivity for the correct subelement. */
    connSubElems = standardElementsGrid[ind].GetSubConnType2() + (subElem - nSubElemType1) * nDOFsPerSubElem;
  }

  /* Get the parametric coordinates of the DOFs from the standard element. */
  const vector<su2double>* locDOFs[] = {standardElementsGrid[ind].GetRDOFs(), standardElementsGrid[ind].GetSDOFs(),
                                        standardElementsGrid[ind].GetTDOFs()};

  /* Create the initial guess of the parametric coordinates by interpolation
     in the sub-element. */
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
    parCoor[iDim] = 0.0;
    const su2double* coorDOFs = locDOFs[iDim]->data();
    for (unsigned short i = 0; i < nDOFsPerSubElem; ++i) parCoor[iDim] += weightsSubElem[i] * coorDOFs[connSubElems[i]];
  }

  /*--------------------------------------------------------------------------*/
  /* Step 2: The Newton algorithm to compute the parametric coordinates in    */
  /*         the high order element.                                          */
  /*--------------------------------------------------------------------------*/

  /* Determine the number of DOFs of the standard element and allocate the
     memory to store the basis functions and its derivatives. */
  const unsigned short nDOFs = standardElementsGrid[ind].GetNDOFs();

  vector<su2double> lagBasis(nDOFs);
  vector<vector<su2double> > dLagBasis(nDim, vector<su2double>(nDOFs));

  /* Abbreviate the grid DOFs of this element a bit easier. */
  const unsigned long* DOFs = volElem[parElem].nodeIDsGrid.data();

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for (itCount = 0; itCount < maxIt; ++itCount) {
    /* Compute the Lagrangian basis functions and its derivatives in
       the current parametric coordinate. */
    standardElementsGrid[ind].BasisFunctionsAndDerivativesInPoint(parCoor, lagBasis, dLagBasis);

    /* Make a distinction between 2D and 3D in order to have the most
       efficient code. */
    bool converged = false;
    switch (nDim) {
      case 2: {
        /* Two dimensional computation. Compute the values of the function
           and minus the Jacobian matrix. */
        su2double f0 = coor[0], f1 = coor[1];
        su2double a00 = 0.0, a01 = 0.0, a10 = 0.0, a11 = 0.0;
        for (unsigned short i = 0; i < nDOFs; ++i) {
          const su2double x = meshPoints[DOFs[i]].coor[0];
          const su2double y = meshPoints[DOFs[i]].coor[1];

          f0 -= x * lagBasis[i];
          f1 -= y * lagBasis[i];

          a00 += x * dLagBasis[0][i];
          a01 += x * dLagBasis[1][i];
          a10 += y * dLagBasis[0][i];
          a11 += y * dLagBasis[1][i];
        }

        /* Compute the updates of the parametric values. As minus the
           Jacobian is computed, the updates should be added to parCoor. */
        const su2double detInv = 1.0 / (a00 * a11 - a01 * a10);
        const su2double dr = detInv * (f0 * a11 - f1 * a01);
        const su2double ds = detInv * (f1 * a00 - f0 * a10);

        parCoor[0] += dr;
        parCoor[1] += ds;

        /* Check for convergence. */
        if (fabs(dr) <= tolNewton && fabs(ds) <= tolNewton) converged = true;
        break;
      }

      case 3: {
        /* Three dimensional computation. Compute the values of the function
           and minus the Jacobian matrix. */
        su2double f0 = coor[0], f1 = coor[1], f2 = coor[2];
        su2double a00 = 0.0, a01 = 0.0, a02 = 0.0;
        su2double a10 = 0.0, a11 = 0.0, a12 = 0.0;
        su2double a20 = 0.0, a21 = 0.0, a22 = 0.0;

        for (unsigned short i = 0; i < nDOFs; ++i) {
          const su2double x = meshPoints[DOFs[i]].coor[0];
          const su2double y = meshPoints[DOFs[i]].coor[1];
          const su2double z = meshPoints[DOFs[i]].coor[2];

          f0 -= x * lagBasis[i];
          f1 -= y * lagBasis[i];
          f2 -= z * lagBasis[i];

          a00 += x * dLagBasis[0][i];
          a01 += x * dLagBasis[1][i];
          a02 += x * dLagBasis[2][i];
          a10 += y * dLagBasis[0][i];
          a11 += y * dLagBasis[1][i];
          a12 += y * dLagBasis[2][i];
          a20 += z * dLagBasis[0][i];
          a21 += z * dLagBasis[1][i];
          a22 += z * dLagBasis[2][i];
        }

        /* Compute the updates of the parametric values. As minus the
           Jacobian is computed, the updates should be added to parCoor. */
        const su2double detInv = 1.0 / (a00 * a11 * a22 - a00 * a12 * a21 - a01 * a10 * a22 + a01 * a12 * a20 +
                                        a02 * a10 * a21 - a02 * a11 * a20);
        const su2double dr = detInv * (a01 * a12 * f2 - a01 * a22 * f1 - a02 * a11 * f2 + a02 * a21 * f1 +
                                       a11 * a22 * f0 - a12 * a21 * f0);
        const su2double ds = -detInv * (a00 * a12 * f2 - a00 * a22 * f1 - a02 * a10 * f2 + a02 * a20 * f1 +
                                        a10 * a22 * f0 - a12 * a20 * f0);
        const su2double dt = detInv * (a00 * a11 * f2 - a00 * a21 * f1 - a01 * a10 * f2 + a01 * a20 * f1 +
                                       a10 * a21 * f0 - a11 * a20 * f0);
        parCoor[0] += dr;
        parCoor[1] += ds;
        parCoor[2] += dt;

        /* Check for convergence. */
        if (fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton) converged = true;
        break;
      }
    }

    /* Break the loop if the Newton algorithm converged. */
    if (converged) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if (itCount == maxIt) SU2_MPI::Error("Newton did not converge", CURRENT_FUNCTION);
}

void CMeshFEM_DG::InitStaticMeshMovement(const CConfig* config, const unsigned short Kind_Grid_Movement,
                                         const unsigned short iZone) {
  /*--- Get the reference values for the non-dimensionalization of the
        prescribed velocities. ---*/
  const su2double L_Ref = config->GetLength_Ref();
  const su2double Omega_Ref = config->GetOmega_Ref();
  const su2double Vel_Ref = config->GetVelocity_Ref();

  /*--- Make a distinction between the possibilities. ---*/
  switch (Kind_Grid_Movement) {
      /*-------------------------------------------------------------------------------------*/

    case ROTATING_FRAME: {
      /* Get the rotation rate and rotation center from config. */
      const su2double Center[] = {config->GetMotion_Origin(0), config->GetMotion_Origin(1),
                                  config->GetMotion_Origin(2)};
      const su2double Omega[] = {config->GetRotation_Rate(0) / Omega_Ref, config->GetRotation_Rate(1) / Omega_Ref,
                                 config->GetRotation_Rate(2) / Omega_Ref};

      /* Array used to store the distance to the rotation center. */
      su2double dist[] = {0.0, 0.0, 0.0};

      /* Loop over the owned volume elements. */
      for (unsigned long l = 0; l < nVolElemOwned; ++l) {
        /* Get the number of DOFs and integration points for this element. */
        const unsigned short ind = volElem[l].indStandardElement;
        const unsigned short nInt = standardElementsSol[ind].GetNIntegration();
        const unsigned short nDOFs = volElem[l].nDOFsSol;

        /* Determine the grid velocities in the integration points. */
        for (unsigned short i = 0; i < nInt; ++i) {
          const su2double* coor = volElem[l].coorIntegrationPoints.data() + i * nDim;
          su2double* gridVel = volElem[l].gridVelocities.data() + i * nDim;

          for (unsigned short iDim = 0; iDim < nDim; ++iDim) dist[iDim] = (coor[iDim] - Center[iDim]) / L_Ref;

          gridVel[0] = Omega[1] * dist[2] - Omega[2] * dist[1];
          gridVel[1] = Omega[2] * dist[0] - Omega[0] * dist[2];
          if (nDim == 3) gridVel[2] = Omega[0] * dist[1] - Omega[1] * dist[0];
        }

        /* Determine the grid velocities in the solution DOFs. */
        for (unsigned short i = 0; i < nDOFs; ++i) {
          const su2double* coor = volElem[l].coorSolDOFs.data() + i * nDim;
          su2double* gridVel = volElem[l].gridVelocitiesSolDOFs.data() + i * nDim;

          for (unsigned short iDim = 0; iDim < nDim; ++iDim) dist[iDim] = (coor[iDim] - Center[iDim]) / L_Ref;

          gridVel[0] = Omega[1] * dist[2] - Omega[2] * dist[1];
          gridVel[1] = Omega[2] * dist[0] - Omega[0] * dist[2];
          if (nDim == 3) gridVel[2] = Omega[0] * dist[1] - Omega[1] * dist[0];
        }
      }

      /* Loop over the internal matching faces. */
      for (unsigned long l = 0; l < matchingFaces.size(); ++l) {
        /* Get the number of integration points. */
        const unsigned short ind = matchingFaces[l].indStandardElement;
        const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();

        /* Determine the grid velocity in the integration points. */
        for (unsigned short i = 0; i < nInt; ++i) {
          const su2double* coor = matchingFaces[l].coorIntegrationPoints.data() + i * nDim;
          su2double* gridVel = matchingFaces[l].gridVelocities.data() + i * nDim;

          for (unsigned short iDim = 0; iDim < nDim; ++iDim) dist[iDim] = (coor[iDim] - Center[iDim]) / L_Ref;

          gridVel[0] = Omega[1] * dist[2] - Omega[2] * dist[1];
          gridVel[1] = Omega[2] * dist[0] - Omega[0] * dist[2];
          if (nDim == 3) gridVel[2] = Omega[0] * dist[1] - Omega[1] * dist[0];
        }
      }

      /* Loop over the physical boundaries. Exclude periodic boundaries. */
      for (unsigned short iMarker = 0; iMarker < boundaries.size(); ++iMarker) {
        if (!boundaries[iMarker].periodicBoundary) {
          /* Abbreviate the surface elements a bit easier. */
          vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;

          /* Check if this is a shroud boundary. */
          bool shroudBoundary = false;
          for (unsigned short iShroud = 0; iShroud < config->GetnMarker_Shroud(); ++iShroud) {
            if (boundaries[iMarker].markerTag == config->GetMarker_Shroud(iShroud)) {
              shroudBoundary = true;
              break;
            }
          }

          /* Loop over the boundary faces of this marker. */
          for (unsigned long l = 0; l < surfElem.size(); ++l) {
            /* Get the number of integration points. */
            const unsigned short ind = surfElem[l].indStandardElement;
            const unsigned short nInt = standardBoundaryFacesGrid[ind].GetNIntegration();

            /* For a shroud boundary the grid velocities must be set to zero. */
            if (shroudBoundary) {
              for (unsigned short i = 0; i < (nInt * nDim); ++i) surfElem[l].gridVelocities[i] = 0.0;
            } else {
              /* Normal boundary. Determine the grid velocity in the integration points. */
              for (unsigned short i = 0; i < nInt; ++i) {
                const su2double* coor = surfElem[l].coorIntegrationPoints.data() + i * nDim;
                su2double* gridVel = surfElem[l].gridVelocities.data() + i * nDim;

                for (unsigned short iDim = 0; iDim < nDim; ++iDim) dist[iDim] = (coor[iDim] - Center[iDim]) / L_Ref;

                gridVel[0] = Omega[1] * dist[2] - Omega[2] * dist[1];
                gridVel[1] = Omega[2] * dist[0] - Omega[0] * dist[2];
                if (nDim == 3) gridVel[2] = Omega[0] * dist[1] - Omega[1] * dist[0];
              }
            }
          }
        }
      }

      break;
    }

      /*-------------------------------------------------------------------------------------*/

    case STEADY_TRANSLATION: {
      /* Get the translation velocity from config. */
      const su2double vTrans[] = {config->GetTranslation_Rate(0) / Vel_Ref, config->GetTranslation_Rate(1) / Vel_Ref,
                                  config->GetTranslation_Rate(2) / Vel_Ref};

      /* Loop over the owned volume elements. */
      for (unsigned long l = 0; l < nVolElemOwned; ++l) {
        /* Get the number of DOFs and integration points for this element. */
        const unsigned short ind = volElem[l].indStandardElement;
        const unsigned short nInt = standardElementsSol[ind].GetNIntegration();
        const unsigned short nDOFs = volElem[l].nDOFsSol;

        /* Set the grid velocity in both the integration points and
           the sol DOFs. */
        for (unsigned short i = 0; i < nInt; ++i) {
          su2double* gridVel = volElem[l].gridVelocities.data() + i * nDim;
          for (unsigned short iDim = 0; iDim < nDim; ++iDim) gridVel[iDim] = vTrans[iDim];
        }

        for (unsigned short i = 0; i < nDOFs; ++i) {
          su2double* gridVel = volElem[l].gridVelocitiesSolDOFs.data() + i * nDim;
          for (unsigned short iDim = 0; iDim < nDim; ++iDim) gridVel[iDim] = vTrans[iDim];
        }
      }

      /* Loop over the internal matching faces. */
      for (unsigned long l = 0; l < matchingFaces.size(); ++l) {
        /* Get the number of integration points. */
        const unsigned short ind = matchingFaces[l].indStandardElement;
        const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();

        /* Set the grid velocity in the integration points. */
        for (unsigned short i = 0; i < nInt; ++i) {
          su2double* gridVel = matchingFaces[l].gridVelocities.data() + i * nDim;
          for (unsigned short iDim = 0; iDim < nDim; ++iDim) gridVel[iDim] = vTrans[iDim];
        }
      }

      /* Loop over the physical boundaries. Exclude periodic boundaries. */
      for (unsigned short iMarker = 0; iMarker < boundaries.size(); ++iMarker) {
        if (!boundaries[iMarker].periodicBoundary) {
          /* Abbreviate the surface elements a bit easier. */
          vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;

          /* Loop over the boundary faces of this marker. */
          for (unsigned long l = 0; l < surfElem.size(); ++l) {
            /* Get the number of integration points. */
            const unsigned short ind = surfElem[l].indStandardElement;
            const unsigned short nInt = standardBoundaryFacesGrid[ind].GetNIntegration();

            /* Set the grid velocity in the integration points. */
            for (unsigned short i = 0; i < nInt; ++i) {
              su2double* gridVel = surfElem[l].gridVelocities.data() + i * nDim;
              for (unsigned short iDim = 0; iDim < nDim; ++iDim) gridVel[iDim] = vTrans[iDim];
            }
          }
        }
      }

      break;
    }

      /*-------------------------------------------------------------------------------------*/

    default: /* Just to avoid a compiler warning. */
      break;
  }

  if (config->GetSurface_Movement(MOVING_WALL)) {
    /*--- Loop over the physical boundaries. Skip the periodic boundaries. ---*/
    for (unsigned short i = 0; i < boundaries.size(); ++i) {
      if (!boundaries[i].periodicBoundary) {
        /* Check if for this boundary a motion has been specified. */
        if (config->GetMarker_All_Moving(i) == YES) {
          /* Determine the prescribed translation velocity, rotation rate
               and rotation center. */
          const su2double Center[] = {config->GetMotion_Origin(0), config->GetMotion_Origin(1),
                                      config->GetMotion_Origin(2)};
          const su2double Omega[] = {config->GetRotation_Rate(0) / Omega_Ref, config->GetRotation_Rate(1) / Omega_Ref,
                                     config->GetRotation_Rate(2) / Omega_Ref};
          const su2double vTrans[] = {config->GetTranslation_Rate(0) / Vel_Ref,
                                      config->GetTranslation_Rate(1) / Vel_Ref,
                                      config->GetTranslation_Rate(2) / Vel_Ref};

          /* Easier storage of the surface elements and loop over them. */
          vector<CSurfaceElementFEM>& surfElem = boundaries[i].surfElem;

          for (unsigned long l = 0; l < surfElem.size(); ++l) {
            /* Determine the corresponding standard face element and get the
                 relevant information from it. Note that the standard element
                 of the solution must be taken and not of the grid. */
            const unsigned short ind = surfElem[l].indStandardElement;
            const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

            /* Loop over the number of integration points. */
            for (unsigned short j = 0; j < nInt; ++j) {
              /* Set the pointers for the coordinates and grid velocities
                   for this integration points. */
              const su2double* Coord = surfElem[l].coorIntegrationPoints.data() + j * nDim;
              su2double* gridVel = surfElem[l].gridVelocities.data() + j * nDim;

              /* Calculate non-dim. position from rotation center. */
              su2double r[] = {0.0, 0.0, 0.0};
              for (unsigned short iDim = 0; iDim < nDim; ++iDim) r[iDim] = (Coord[iDim] - Center[iDim]) / L_Ref;

              /* Cross Product of angular velocity and distance from center to
                   get the rotational velocity. Note that we are adding on the
                   velocity due to pure translation as well. Note that for the
                   2D case only Omega[2] can be non-zero. */
              su2double velGrid[] = {vTrans[0] + Omega[1] * r[2] - Omega[2] * r[1],
                                     vTrans[1] + Omega[2] * r[0] - Omega[0] * r[2],
                                     vTrans[2] + Omega[0] * r[1] - Omega[1] * r[0]};

              /* Store the grid velocities. */
              for (unsigned short iDim = 0; iDim < nDim; ++iDim) gridVel[iDim] = velGrid[iDim];
            }
          }
        }
      }
    }
  }
}

CDummyMeshFEM_DG::CDummyMeshFEM_DG(CConfig* config) : CMeshFEM_DG() {
  nZone = config->GetnZone();

  nPoint_P2PSend = new int[size]();
  nPoint_P2PRecv = new int[size]();

  nVertex = new unsigned long[config->GetnMarker_All()]();

  Tag_to_Marker = new string[config->GetnMarker_All()];

  for (unsigned short i = 0; i <= config->GetnLevels_TimeAccurateLTS(); i++) {
    nMatchingFacesWithHaloElem.push_back(0);
  }

  boundaries.resize(config->GetnMarker_All());

  nDim = CConfig::GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());
}
