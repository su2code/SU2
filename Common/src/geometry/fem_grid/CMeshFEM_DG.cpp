/*!
 * \file CMeshFEM_DG.cpp
 * \brief Implementations of the member functions of CMeshFEM_DG.
 * \author E. van der Weide
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

#include "../../../include/geometry/fem_grid/CMeshFEM_DG.hpp"
#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/toolboxes/fem/CFaceOfElement.hpp"
#include "../../../include/toolboxes/fem/CReorderElements.hpp"
#include "../../../include/toolboxes/fem/CSortFaces.hpp"
#include "../../../include/adt/CADTPointsOnlyClass.hpp"
#include "../../../include/fem/CFEMStandardHexVolumeSol.hpp"
#include "../../../include/fem/CFEMStandardPrismVolumeSol.hpp"
#include "../../../include/fem/CFEMStandardPyraVolumeSol.hpp"
#include "../../../include/fem/CFEMStandardQuadVolumeSol.hpp"
#include "../../../include/fem/CFEMStandardTetVolumeSol.hpp"
#include "../../../include/fem/CFEMStandardTriVolumeSol.hpp"
#include "../../../include/fem/CGemmFaceHex.hpp"
#include "../../../include/fem/CGemmFaceQuad.hpp"
#include "../../../include/fem/CGemmStandard.hpp"
#include "../../../include/fem/CFEMStandardLineAdjacentTriGrid.hpp"
#include "../../../include/fem/CFEMStandardLineAdjacentQuadGrid.hpp"
#include "../../../include/fem/CFEMStandardTriAdjacentTetGrid.hpp"
#include "../../../include/fem/CFEMStandardTriAdjacentPyraGrid.hpp"
#include "../../../include/fem/CFEMStandardTriAdjacentPrismGrid.hpp"
#include "../../../include/fem/CFEMStandardQuadAdjacentPyraGrid.hpp"
#include "../../../include/fem/CFEMStandardQuadAdjacentHexGrid.hpp"
#include "../../../include/fem/CFEMStandardQuadAdjacentPrismGrid.hpp"
#include "../../../include/fem/CFEMStandardLineAdjacentTriSol.hpp"
#include "../../../include/fem/CFEMStandardLineAdjacentQuadSol.hpp"
#include "../../../include/fem/CFEMStandardTriAdjacentTetSol.hpp"
#include "../../../include/fem/CFEMStandardTriAdjacentPyraSol.hpp"
#include "../../../include/fem/CFEMStandardTriAdjacentPrismSol.hpp"
#include "../../../include/fem/CFEMStandardQuadAdjacentPyraSol.hpp"
#include "../../../include/fem/CFEMStandardQuadAdjacentHexSol.hpp"
#include "../../../include/fem/CFEMStandardQuadAdjacentPrismSol.hpp"

/*---------------------------------------------------------------------*/
/*---          Public member functions of CMeshFEM_DG.              ---*/
/*---------------------------------------------------------------------*/

CMeshFEM_DG::~CMeshFEM_DG(void) {

  /*--- Release the memory of the standard elements. ---*/
  for(unsigned long i=0; i<standardVolumeElementsSolution.size(); ++i) {
    if( standardVolumeElementsSolution[i] ) delete standardVolumeElementsSolution[i];
    standardVolumeElementsSolution[i] = nullptr;
  }

  for(unsigned long i=0; i<standardSurfaceElementsSolution.size(); ++i) {
    if( standardSurfaceElementsSolution[i] ) delete standardSurfaceElementsSolution[i];
    standardSurfaceElementsSolution[i] = nullptr;
  }

  for(unsigned long i=0; i<standardInternalFaceGrid.size(); ++i) {
    if( standardInternalFaceGrid[i] ) delete standardInternalFaceGrid[i];
    standardInternalFaceGrid[i] = nullptr;
  }

  for(unsigned long i=0; i<standardInternalFaceSolution.size(); ++i) {
    if( standardInternalFaceSolution[i] ) delete standardInternalFaceSolution[i];
    standardInternalFaceSolution[i] = nullptr;
  }
}

CMeshFEM_DG::CMeshFEM_DG(CGeometry *geometry, CConfig *config)
  : CMeshFEM_Base(geometry, config) {

  /*--- Define the linear partitioning of the elements. ---*/
  CLinearPartitioner elemPartitioner(Global_nElem, 0);

  /*--- Determine a mapping from the global point ID to the local index
        of the points. ---*/
  map<unsigned long,unsigned long> globalPointIDToLocalInd;
  for(unsigned long i=0; i<geometry->GetnPoint(); ++i)
    globalPointIDToLocalInd[geometry->nodes->GetGlobalIndex(i)] = i;

  /*----------------------------------------------------------------------------*/
  /*--- Step 1: Determine the global offset for this rank in the numbering   ---*/
  /*---         of the solution DOFs of the elements. This value is needed   ---*/
  /*---         when writing the solution and much easier to determine on    ---*/
  /*---         the original linear partitioning than the final partitioning.---*/
  /*----------------------------------------------------------------------------*/

  /*--- Make a distinction between parallel and sequential mode. ---*/
#ifdef HAVE_MPI

  /*--- Determine the number of solution DOFs stored on the current MPI rank. ---*/
  long nSolDOFs = 0;
  for(unsigned long i=0; i<geometry->GetnElem(); ++i)
    nSolDOFs += geometry->elem[i]->GetNDOFsSol();

  /*--- Gather the number of solution DOFs of all ranks. ---*/
  vector<long> nSolDOFsAllRanks(size);
  SU2_MPI::Allgather(&nSolDOFs, 1, MPI_LONG, nSolDOFsAllRanks.data(), 1,
                     MPI_LONG, SU2_MPI::GetComm());

  /*--- Determine the offset for the current rank. ---*/
  long offsetSolDOFs = 0;
  for(int i=0; i<rank; ++i)
    offsetSolDOFs += nSolDOFsAllRanks[i];

  /*--- Release the memory of nSolDOFsAllRanks again. ---*/
  vector<long>().swap(nSolDOFsAllRanks);

#else
  /*--- Sequential mode. The offset of this rank is zero. ---*/
  long offsetSolDOFs = 0;
#endif

  /*----------------------------------------------------------------------------*/
  /*--- Step 2: Communicate the elements and the boundary elements to the    ---*/
  /*---         ranks where they will be stored during the computation.      ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Determine the ranks to which I have to send my elements. ---*/
  vector<int> sendToRank(size, 0);

  for(unsigned long i=0; i<geometry->GetnElem(); ++i) {
    sendToRank[geometry->elem[i]->GetColor()] = 1;
  }

  map<int,int> rankToIndCommBuf;
  for(int i=0; i<size; ++i) {
    if( sendToRank[i] ) {
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  /*--- Definition of the communication buffers, used to send the element data
        to the correct ranks. ---*/
  int nRankSend = (int)rankToIndCommBuf.size();
  vector<vector<short> >     shortSendBuf(nRankSend,  vector<short>(0));
  vector<vector<long>  >     longSendBuf(nRankSend,   vector<long>(0));
  vector<vector<su2double> > doubleSendBuf(nRankSend, vector<su2double>(0));

  /*--- The first element of longSendBuf will contain the number of elements, which
        are stored in the communication buffers. Initialize this value to 0. ---*/
  for(int i=0; i<nRankSend; ++i) longSendBuf[i].push_back(0);

  /*--- Determine the number of ranks, from which this rank will receive elements. ---*/
  int nRankRecv = nRankSend;

#ifdef HAVE_MPI
  vector<int> sizeRecv(size, 1);

  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- Loop over the local elements to fill the communication buffers with element data. ---*/
  for(unsigned long i=0; i<geometry->GetnElem(); ++i) {
    int ind = (int) geometry->elem[i]->GetColor();
    map<int,int>::const_iterator MI = rankToIndCommBuf.find(ind);
    ind = MI->second;

    ++longSendBuf[ind][0];   /* The number of elements in the buffers must be incremented. */

    shortSendBuf[ind].push_back(geometry->elem[i]->GetVTK_Type());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNPolyGrid());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNPolySol());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNDOFsGrid());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNDOFsSol());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetnFaces());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetTimeLevel());
    shortSendBuf[ind].push_back( (short) geometry->elem[i]->GetJacobianConsideredConstant());

    longSendBuf[ind].push_back(geometry->elem[i]->GetGlobalElemID());
    longSendBuf[ind].push_back(offsetSolDOFs);
    offsetSolDOFs += geometry->elem[i]->GetNDOFsSol();

    for(unsigned short j=0; j<geometry->elem[i]->GetNDOFsGrid(); ++j)
      longSendBuf[ind].push_back(geometry->elem[i]->GetNode(j));

    for(unsigned short j=0; j<geometry->elem[i]->GetnFaces(); ++j)
      longSendBuf[ind].push_back(geometry->elem[i]->GetNeighbor_Elements(j));

    for(unsigned short j=0; j<geometry->elem[i]->GetnFaces(); ++j) {
      shortSendBuf[ind].push_back(geometry->elem[i]->GetPeriodicIndex(j));
      shortSendBuf[ind].push_back( (short) geometry->elem[i]->GetJacobianConstantFace(j));
      shortSendBuf[ind].push_back( (short) geometry->elem[i]->GetOwnerFace(j));
    }

    doubleSendBuf[ind].push_back(geometry->elem[i]->GetLengthScale());
  }

  /*--- Determine for each rank to which I have to send elements the data of
        the corresponding nodes. ---*/
  for(int i=0; i<nRankSend; ++i) {

    /*--- Determine the vector with node IDs in the connectivity
          of the elements for this rank. ---*/
    vector<long> nodeIDs;

    unsigned long indL = 3;
    unsigned long indS = 3;
    for(long j=0; j<longSendBuf[i][0]; ++j) {
      short nDOFsGrid = shortSendBuf[i][indS], nFaces = shortSendBuf[i][indS+2];
      indS += 3*nFaces+8;

      for(short k=0; k<nDOFsGrid; ++k, ++indL)
        nodeIDs.push_back(longSendBuf[i][indL]);
      indL += nFaces+2;
    }

    /*--- Sort nodeIDs in increasing order and remove the double entities. ---*/
    std::sort(nodeIDs.begin(), nodeIDs.end());
    vector<long>::iterator lastNodeID = unique(nodeIDs.begin(), nodeIDs.end());
    nodeIDs.erase(lastNodeID, nodeIDs.end());

    /*--- Add the number of node IDs and the node IDs itself to longSendBuf[i]. ---*/
    longSendBuf[i].push_back(nodeIDs.size());
    longSendBuf[i].insert(longSendBuf[i].end(), nodeIDs.begin(), nodeIDs.end());

    /*--- Copy the coordinates to doubleSendBuf. ---*/
    for(unsigned long j=0; j<nodeIDs.size(); ++j) {
      map<unsigned long,unsigned long>::const_iterator LMI;
      LMI = globalPointIDToLocalInd.find(nodeIDs[j]);

      if(LMI == globalPointIDToLocalInd.end())
        SU2_MPI::Error("Entry not found in map", CURRENT_FUNCTION);

      unsigned long ind = LMI->second;
      for(unsigned short l=0; l<nDim; ++l)
        doubleSendBuf[i].push_back(geometry->nodes->GetCoord(ind, l));
    }
  }

  /*--- Loop over the boundaries to send the boundary data to the appropriate rank. ---*/
  nMarker = geometry->GetnMarker();
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /*--- Store the current indices in the longSendBuf, which are used to store the
          number of boundary elements sent to this rank. Initialize this value to 0. ---*/
    vector<long> indLongBuf(nRankSend);
    for(int i=0; i<nRankSend; ++i) {
      indLongBuf[i] = longSendBuf[i].size();
      longSendBuf[i].push_back(0);
    }

    /*--- Loop over the local boundary elements in geometry for this marker. ---*/
    for(unsigned long i=0; i<geometry->GetnElem_Bound(iMarker); ++i) {

      /*--- Determine the local ID of the corresponding domain element. ---*/
      unsigned long elemID = geometry->bound[iMarker][i]->GetDomainElement()
                           - elemPartitioner.GetFirstIndexOnRank(rank);

      if(elemID >= geometry->GetnElem())
        SU2_MPI::Error("Invalid element ID", CURRENT_FUNCTION);

      /*--- Determine to which rank this boundary element must be sent.
            That is the same as its corresponding domain element.
            Update the corresponding index in longSendBuf. ---*/
      int ind = (int) geometry->elem[elemID]->GetColor();
      map<int,int>::const_iterator MI = rankToIndCommBuf.find(ind);
      ind = MI->second;

      ++longSendBuf[ind][indLongBuf[ind]];

      /*--- Get the donor information for the wall function treatment. ---*/
      const unsigned short nDonors = geometry->bound[iMarker][i]->GetNDonorsWallFunctions();
      const unsigned long  *donors = geometry->bound[iMarker][i]->GetDonorsWallFunctions();

      /*--- Store the data for this boundary element in the communication buffers. ---*/
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetVTK_Type());
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNPolyGrid());
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNDOFsGrid());
      shortSendBuf[ind].push_back(nDonors);

      longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetDomainElement());
      longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetGlobalElemID());

      for(unsigned short j=0; j<geometry->bound[iMarker][i]->GetNDOFsGrid(); ++j)
        longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNode(j));

      for(unsigned short j=0; j<nDonors; ++j)
        longSendBuf[ind].push_back(donors[j]);
    }
  }

  /*--- Definition of the communication buffers, used to receive
        the element data from the other correct ranks. ---*/
  vector<vector<short> >     shortRecvBuf(nRankRecv,  vector<short>(0));
  vector<vector<long>  >     longRecvBuf(nRankRecv,   vector<long>(0));
  vector<vector<su2double> > doubleRecvBuf(nRankRecv, vector<su2double>(0));

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode. ---*/
  map<int,int>::const_iterator MI;

#ifdef HAVE_MPI

  /*--- Parallel mode. Send all the data using non-blocking sends. ---*/
  vector<SU2_MPI::Request> commReqs(3*nRankSend);
  MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {

    int dest = MI->first;
    SU2_MPI::Isend(shortSendBuf[i].data(), shortSendBuf[i].size(), MPI_SHORT,
                   dest, dest, SU2_MPI::GetComm(), &commReqs[3*i]);
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest+1, SU2_MPI::GetComm(), &commReqs[3*i+1]);
    SU2_MPI::Isend(doubleSendBuf[i].data(), doubleSendBuf[i].size(), MPI_DOUBLE,
                   dest, dest+2, SU2_MPI::GetComm(), &commReqs[3*i+2]);
  }

  /*--- Loop over the number of ranks from which I receive data. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /*--- Block until a message with shorts arrives from any processor.
          Determine the source and the size of the message. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_SHORT, &sizeMess);

    /*--- Allocate the memory for the short receive buffer and receive the message. ---*/
    shortRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(shortRecvBuf[i].data(), sizeMess, MPI_SHORT,
                  source, rank, SU2_MPI::GetComm(), &status);

    /*--- Block until the corresponding message with longs arrives, determine
          its size, allocate the memory and receive the message. ---*/
    SU2_MPI::Probe(source, rank+1, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);
    longRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(longRecvBuf[i].data(), sizeMess, MPI_LONG,
                  source, rank+1, SU2_MPI::GetComm(), &status);

    /*--- Idem for the message with doubles. ---*/
    SU2_MPI::Probe(source, rank+2, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_DOUBLE, &sizeMess);
    doubleRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(doubleRecvBuf[i].data(), sizeMess, MPI_DOUBLE,
                  source, rank+2, SU2_MPI::GetComm(), &status);
  }

  /*--- Complete the non-blocking sends. ---*/
  SU2_MPI::Waitall(3*nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Wild cards have been used in the communication,
        so synchronize the ranks to avoid problems. ---*/
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffers. ---*/
  shortRecvBuf[0]  = shortSendBuf[0];
  longRecvBuf[0]   = longSendBuf[0];
  doubleRecvBuf[0] = doubleSendBuf[0];

#endif

  /*--- Release the memory of the send buffers. To make sure that all
        the memory is deleted, the swap function is used. ---*/
  for(int i=0; i<nRankSend; ++i) {
    vector<short>().swap(shortSendBuf[i]);
    vector<long>().swap(longSendBuf[i]);
    vector<su2double>().swap(doubleSendBuf[i]);
  }

  /*--- Allocate the memory for the number of elements for every boundary
        marker and initialize them to zero. ---*/
  nElem_Bound = new unsigned long[nMarker];
  for(unsigned short i=0; i<nMarker; ++i)
    nElem_Bound[i] = 0;

  /*--- Determine the global element ID's of the elements stored on this rank.
        Sort them in increasing order, such that an easy search can be done.
        In the same loop determine the upper bound for the local nodes (without
        halos) and the number of boundary elements for every marker. ---*/
  nElem = nPoint = 0;
  for(int i=0; i<nRankRecv; ++i) nElem += longRecvBuf[i][0];

  vector<unsigned long> globalElemID;
  globalElemID.reserve(nElem);

  for(int i=0; i<nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for(long j=0; j<longRecvBuf[i][0]; ++j) {
      globalElemID.push_back(longRecvBuf[i][indL]);

      unsigned short nDOFsGrid = shortRecvBuf[i][indS+3];
      unsigned short nFaces    = shortRecvBuf[i][indS+5];
      indS += 3*nFaces + 8;
      indL += nDOFsGrid + nFaces + 2;
    }

    long nNodesThisRank = longRecvBuf[i][indL];
    nPoint += nNodesThisRank;
    indL   += nNodesThisRank+1;

    for(unsigned iMarker=0; iMarker<nMarker; ++iMarker) {
      long nBoundElemThisRank = longRecvBuf[i][indL]; ++indL;
      nElem_Bound[iMarker] += nBoundElemThisRank;

      for(long j=0; j<nBoundElemThisRank; ++j) {
        short nDOFsBoundElem      = shortRecvBuf[i][indS+2];
        short nDonorsWallFunction = shortRecvBuf[i][indS+3];
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
        see step 3 below. ---*/
  vector<CUnsignedLong2T>  haloElements;
  vector<CReorderElements> ownedElements;
  unsigned short maxTimeLevelLoc = 0;

  for(int i=0; i<nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for(long j=0; j<longRecvBuf[i][0]; ++j) {
      unsigned long  globalID    = longRecvBuf[i][indL];
      unsigned short VTK_Type    = shortRecvBuf[i][indS];
      unsigned short nPolySol    = shortRecvBuf[i][indS+2];
      unsigned short nDOFsGrid   = shortRecvBuf[i][indS+3];
      unsigned short nFaces      = shortRecvBuf[i][indS+5];
      unsigned short timeLevel   = shortRecvBuf[i][indS+6];
      bool           JacConstant = (bool) shortRecvBuf[i][indS+7];
      bool           commSol     = false;

      /*--- If the integration rule for constant and non-constant Jacobian
            elements is the same, simply set JacConstant to false, such that
            the elements with constant and non-constant Jacobians are
            considered the same. ---*/
      if( JacConstant ) {
        const unsigned short orderExactStraight = config->GetOrderExactIntegrationFEM(nPolySol, true);
        const unsigned short orderExactCurved   = config->GetOrderExactIntegrationFEM(nPolySol, false);
        if(orderExactStraight == orderExactCurved) JacConstant = false;
      }

      /*--- Update the local value of the maximum time level. ---*/
      maxTimeLevelLoc = max(maxTimeLevelLoc, timeLevel);

      /*--- Loop over the faces of this element to determine the halo elements
            and the information needed to reorder the owned elements. ---*/
      indS += 8;
      indL += nDOFsGrid + 2;
      for(unsigned short k=0; k<nFaces; ++k, indS+=3, ++indL) {

        if(longRecvBuf[i][indL] != -1) {  /* -1 indicates a boundary face. */

          /*--- Check if the neighbor of the face is also an owned element.
                Per definition an element for which a periodic transformation
                is carried out, is a halo element, even if the parent element
                is stored locally. ---*/
          bool neighborIsInternal = false;
          if(shortRecvBuf[i][indS] == -1) /* -1 indicates no periodic transformation. */
            neighborIsInternal = binary_search(globalElemID.begin(),
                                               globalElemID.end(),
                                               longRecvBuf[i][indL]);

          /*--- Check if this neighbor is not internal and if the element owns the face. ---*/
          if( !neighborIsInternal ) {
            if( shortRecvBuf[i][indS+2] ) {

              /*--- The face is owned by this element. As the neighboring element
                    is not owned, this implies that a halo element must be created. ---*/
              haloElements.push_back(CUnsignedLong2T(longRecvBuf[i][indL],
                                                     shortRecvBuf[i][indS]+1));  /* The +1, because haloElements */
            }                                                                    /* are unsigned longs.          */
            else {

              /*--- The face is not owned by this element and therefore it is owned
                    by the neighboring element on a different rank. Consequently the
                    solution of this element must be communicated. ---*/
              commSol = true;
            }
          }
        }
      }

      /*--- Store the required data for this element in ownedElements. ---*/
      ownedElements.push_back(CReorderElements(globalID, timeLevel, commSol,
                                               VTK_Type, nPolySol, JacConstant));
    }

    /*--- Skip the part with the node numbers. ---*/
    long nNodesThisRank = longRecvBuf[i][indL];
    indL += nNodesThisRank+1;

    /*--- Loop over the boundary markers. ---*/
    for(unsigned iMarker=0; iMarker<nMarker; ++iMarker) {
      long nBoundElemThisRank = longRecvBuf[i][indL]; ++indL;

      /*--- Loop over the boundary elements coming from this rank. ---*/
      for(long j=0; j<nBoundElemThisRank; ++j) {
        short nDOFsBoundElem       = shortRecvBuf[i][indS+2];
        short nDonorsWallFunctions = shortRecvBuf[i][indS+3];
        indS += 4;
        indL += nDOFsBoundElem + 2;

        /*--- Loop over the donors for the wall functions. ---*/
        for(short k=0; k<nDonorsWallFunctions; ++k, ++indL) {

          /*--- Check for an external donor. If external, store it in
                haloElements with no periodic transformation. ---*/
          if( !binary_search(globalElemID.begin(), globalElemID.end(),
                             longRecvBuf[i][indL]) )
            haloElements.push_back(CUnsignedLong2T(longRecvBuf[i][indL], 0));
        }
      }
    }
  }

  /*--- Sort the halo elements in increasing order and remove the double entities. ---*/
  sort(haloElements.begin(), haloElements.end());
  vector<CUnsignedLong2T>::iterator lastHalo = unique(haloElements.begin(), haloElements.end());
  haloElements.erase(lastHalo, haloElements.end());

  /*--- Determine the maximum global time level and possibly reset the number
        of time levels in config. ---*/
  unsigned short maxTimeLevelGlob = maxTimeLevelLoc;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&maxTimeLevelLoc, &maxTimeLevelGlob,
                     1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

  const unsigned short nTimeLevels = maxTimeLevelGlob+1;
  config->SetnLevels_TimeAccurateLTS(nTimeLevels);

  /*----------------------------------------------------------------------------*/
  /*--- Step 3: Find out on which rank the halo elements are stored.         ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Determine the number of owned elements and the total number of
        elements stored on this rank. ---*/
  nVolElemOwned = globalElemID.size();
  nVolElemTot   = nVolElemOwned + haloElements.size();

  /*--- Determine the map from the global element ID to the current storage
        sequence of ownedElements. ---*/
  map<unsigned long, unsigned long> mapGlobalElemIDToInd;
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = i;

  /*--- Determine to which ranks I have to send messages to find out the information
        of the halos stored on this rank. ---*/
  sendToRank.assign(size, 0);

  for(unsigned long i=0; i<haloElements.size(); ++i) {
    const unsigned long rankHalo = elemPartitioner.GetRankContainingIndex(haloElements[i].long0);
    sendToRank[rankHalo] = 1;
  }

  rankToIndCommBuf.clear();
  for(int i=0; i<size; ++i) {
    if( sendToRank[i] ) {
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  /*--- Resize the first index of the long send buffers for the communication of
        the halo data. ---*/
  nRankSend = (int)rankToIndCommBuf.size();
  longSendBuf.resize(nRankSend);

  /*--- Determine the number of ranks, from which this rank will receive elements. ---*/
  nRankRecv = nRankSend;

#ifdef HAVE_MPI
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- Loop over the local halo elements to fill the communication buffers. ---*/
  for(unsigned long i=0; i<haloElements.size(); ++i) {

    /*--- Determine the rank where this halo element was originally stored
          and convert this rank to the index in the send buffer. ---*/
    unsigned long ind = elemPartitioner.GetRankContainingIndex(haloElements[i].long0);
    MI = rankToIndCommBuf.find((int)ind);
    ind = MI->second;

    /*--- Store the global element ID and the periodic index in the long buffer.
          The subtraction of 1 is there to obtain the correct periodic index.
          In haloElements a +1 is added, because this variable is of unsigned long,
          which cannot handle negative numbers. */
    long perIndex = haloElements[i].long1 -1;

    longSendBuf[ind].push_back(haloElements[i].long0);
    longSendBuf[ind].push_back(perIndex);
  }

  /*--- Define a second set of long receive buffers, because the information from
        the first set is still needed later on. Also define the vector to
        store the ranks from which the messages came. ---*/
  vector<vector<long> > longSecondRecvBuf(nRankRecv, vector<long>(0));
  vector<int> sourceRank(nRankRecv);

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode. ---*/

#ifdef HAVE_MPI

  /*--- Parallel mode. Send all the data using non-blocking sends. ---*/
  commReqs.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {
    int dest = MI->first;
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest, SU2_MPI::GetComm(), &commReqs[i]);
  }

  /*--- Loop over the number of ranks from which I receive data. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /*--- Block until a message with longs arrives from any processor.
          Determine the source and the size of the message and receive it. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG,
                  sourceRank[i], rank, SU2_MPI::GetComm(), &status);
  }

  /*--- Complete the non-blocking sends. ---*/
  SU2_MPI::Waitall(nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for(int i=0; i<nRankRecv; ++i)
    longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /*--- Release the memory of the send buffers. To make sure that all the memory
        is deleted, the swap function is used. Afterwards resize the first index
        of the send buffers to nRankRecv, because this number of messages must
        be sent back to the sending ranks with halo information. ---*/
  for(int i=0; i<nRankSend; ++i)
    vector<long>().swap(longSendBuf[i]);

  longSendBuf.resize(nRankRecv);

#ifdef HAVE_MPI
  /*--- Resize the vector of the communication requests to the number of messages
        to be sent by this rank. Only in parallel node. ---*/
  commReqs.resize(nRankRecv);
#endif

  /*--- Loop over the receive buffers to fill and send the send buffers again. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /*--- Determine the number of elements present in longSecondRecvBuf[i] and
          reserve the memory for the send buffer. ---*/
    const long nElemBuf = longSecondRecvBuf[i].size()/2;
    longSendBuf[i].reserve(3*nElemBuf);

    /*--- Loop over the elements stored in the receive buffer. ---*/
    for(long j=0; j<nElemBuf; ++j) {

      /*--- Get the global element ID and periodic index from the receive buffer. ---*/
      const long globalID = longSecondRecvBuf[i][2*j];
      const long perInd   = longSecondRecvBuf[i][2*j+1];

      /*--- Determine the local index of the element in the original partitioning.
            Check if the index is valid. ---*/
      const long localID = globalID - elemPartitioner.GetFirstIndexOnRank(rank);
      if(localID < 0 || localID >= (long) elemPartitioner.GetSizeOnRank(rank) ) {
        ostringstream message;
        message << "Invalid local element ID " << localID 
                << " for original linear partitioning" << endl;
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /*--- Determine which rank owns this element and store everything in the
            send buffer. ---*/
      longSendBuf[i].push_back(globalID);
      longSendBuf[i].push_back(perInd);
      longSendBuf[i].push_back(geometry->elem[localID]->GetColor());
    }

     /*--- Release the memory of this receive buffer. ---*/
    vector<long>().swap(longSecondRecvBuf[i]);

    /*--- Send the send buffer back to the calling rank.
          Only in parallel mode of course. ---*/
#ifdef HAVE_MPI
    int dest = sourceRank[i];
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest+1, SU2_MPI::GetComm(), &commReqs[i]);
#endif
  }

  /*--- Resize the first index of the receive buffers to nRankSend, such that
        the requested halo information can be received. ---*/
  longSecondRecvBuf.resize(nRankSend);

  /*--- Receive the communication data from the correct ranks. Make a distinction
        between parallel and sequential mode. ---*/
#ifdef HAVE_MPI

  /*--- Parallel mode. Loop over the number of ranks from which I receive data
        in the return communication, i.e. nRankSend. ---*/
  for(int i=0; i<nRankSend; ++i) {

    /*--- Block until a message with longs arrives from any processor.
          Determine the source and the size of the message. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    /*--- Allocate the memory for the long receive buffer and receive the message. ---*/
    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG,
                  source, rank+1, SU2_MPI::GetComm(), &status);
  }

  /*--- Complete the non-blocking sends and synchronize the ranks, because
        wild cards have been used. ---*/
  SU2_MPI::Waitall(nRankRecv, commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for(int i=0; i<nRankRecv; ++i)
    longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /*--- Release the memory of the send buffers. To make sure that all
        the memory is deleted, the swap function is used. ---*/
  for(int i=0; i<nRankRecv; ++i)
    vector<long>().swap(longSendBuf[i]);

  /*--- Copy the data from the receive buffers into a class of CLong3T, such that
        it can be sorted in increasing order. Note that the rank of the element
        is stored first, followed by its global ID and last the periodic index. ---*/
  vector<CLong3T> haloData;
  for(int i=0; i<nRankSend; ++i) {
    const long nElemBuf = longSecondRecvBuf[i].size()/3;

    for(long j=0; j<nElemBuf; ++j) {
      const long j3 = 3*j;
      haloData.push_back(CLong3T(longSecondRecvBuf[i][j3+2], longSecondRecvBuf[i][j3],
                                 longSecondRecvBuf[i][j3+1]));
    }

    /*--- Release the memory of this receive buffer. ---*/
    vector<long>().swap(longSecondRecvBuf[i]);
  }

  /*--- Sort halo data in increasing order. ---*/
  sort(haloData.begin(), haloData.end());

  /*--- Determine the number of halo elements per rank in cumulative storage.
        The first element of this vector is nVolElemOwned, such that this vector
        contains the starting position in the vector volElem. ---*/
  vector<unsigned long> nHaloElemPerRank(size+1, 0);
  for(unsigned long i=0; i<haloData.size(); ++i)
    ++nHaloElemPerRank[haloData[i].long0+1];

  nHaloElemPerRank[0] = nVolElemOwned;
  for(int i=0; i<size; ++i)
    nHaloElemPerRank[i+1] += nHaloElemPerRank[i];

  if(nHaloElemPerRank[size] != nVolElemTot)
    SU2_MPI::Error("Inconsistency in total number of volume elements",
                   CURRENT_FUNCTION);

  /*--- Determine the number of ranks to which I have to send data in this cycle. ---*/
  sendToRank.assign(size, 0);
  rankToIndCommBuf.clear();
  for(int i=0; i<size; ++i) {
    if(nHaloElemPerRank[i+1] > nHaloElemPerRank[i]) {
      sendToRank[i] = 1;
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  nRankSend = (int)rankToIndCommBuf.size();

  /*--- Store the value of nRankSend for later use. ---*/
  const int nRankSendHaloInfo = nRankSend;

  /*--- Determine the number of ranks, from which this rank will receive elements. ---*/
  nRankRecv = nRankSend;

#ifdef HAVE_MPI
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- Copy the data to be sent to the send buffers. ---*/
  longSendBuf.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {
    int dest = MI->first;
    for(unsigned long j=nHaloElemPerRank[dest]; j<nHaloElemPerRank[dest+1]; ++j) {
      const unsigned long jj = j - nVolElemOwned;
      longSendBuf[i].push_back(haloData[jj].long1);
      longSendBuf[i].push_back(haloData[jj].long2);
    }
  }

  /*--- Resize the first index of the long receive buffer. ---*/
  longSecondRecvBuf.resize(nRankRecv);

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode. ---*/
#ifdef HAVE_MPI

  /*--- Parallel mode. Send all the data using non-blocking sends. ---*/
  commReqs.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {
    int dest = MI->first;
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest, SU2_MPI::GetComm(), &commReqs[i]);
  }

  /*--- Resize the vector to store the ranks from which the message came. ---*/
  sourceRank.resize(nRankRecv);

  /*--- Loop over the number of ranks from which I receive data. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /*--- Block until a message with longs arrives from any processor.
          Determine the source and the size of the message and receive it. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG,
                  sourceRank[i], rank, SU2_MPI::GetComm(), &status);
  }

  /*--- Complete the non-blocking sends and synchronize the ranks,
        because wild cards have been used. ---*/
  SU2_MPI::Waitall(nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for(int i=0; i<nRankSend; ++i)
    longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /*--- Release the memory of the send buffers. To make sure that all the memory
        is deleted, the swap function is used. ---*/
  for(int i=0; i<nRankSend; ++i)
    vector<long>().swap(longSendBuf[i]);

  /*--- Loop over the receive buffers to flag the locally owned elements for
        communication. Although the face information has already been used to
        do this when ownedElements are constructed, elements that are donors
        for the wall function treatment and are not direct neighbors may
        have been missed. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    const unsigned long nElemBuf = longSecondRecvBuf[i].size()/2;
    for(unsigned long j=0; j<nElemBuf; ++j) {
      const unsigned long elemID = longSecondRecvBuf[i][2*j];

      map<unsigned long, unsigned long>::iterator MMI = mapGlobalElemIDToInd.find(elemID);
      if(MMI == mapGlobalElemIDToInd.end())
        SU2_MPI::Error("Entry not found in mapGlobalElemIDToInd", CURRENT_FUNCTION);

      ownedElements[MMI->second].SetCommSolution(true);
    }
  }

  /*----------------------------------------------------------------------------*/
  /*--- Step 4: Determine the numbering of the owned elements. The following ---*/
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
  /*---           elements does not receive a contribution computed on a     ---*/
  /*---           different rank.                                            ---*/
  /*---         - A reverse Cuthill McKee renumbering takes place to obtain  ---*/
  /*---           better cache performance for the face residuals.           ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Allocate the memory for nVolElemOwnedPerTimeLevel and
        nVolElemInternalPerTimeLevel. The former uses cumulative storage format,
        while the latter uses regular storage format. Hence the +1 for
        nVolElemOwnedPerTimeLevel. ---*/
  nVolElemOwnedPerTimeLevel.resize(nTimeLevels+1, 0);
  nVolElemInternalPerTimeLevel.resize(nTimeLevels, 0);

  /*--- Loop over the owned elements and determine the number of elements
        per time level. ---*/
  for(vector<CReorderElements>::iterator OEI =ownedElements.begin();
                                         OEI!=ownedElements.end(); ++OEI) {
    ++nVolElemOwnedPerTimeLevel[OEI->GetTimeLevel()+1];
    if( !(OEI->GetCommSolution()) )
      ++nVolElemInternalPerTimeLevel[OEI->GetTimeLevel()];
  }

  /*--- Put nVolElemOwnedPerTimeLevel in cumulative storage format. ---*/
  for(unsigned short tLev=0; tLev<nTimeLevels; ++tLev)
    nVolElemOwnedPerTimeLevel[tLev+1] += nVolElemOwnedPerTimeLevel[tLev];

  /*--- Sort the elements of ownedElements in increasing order. ---*/
  sort(ownedElements.begin(), ownedElements.end());

  /*--- At the moment the owned elements are stored per time level,
        whether or not they are internal, per element type and global ID.
        The first two are essential for the code to work. The last two are
        not necessary and can be altered to decrease the bandwith of the
        Jacobian matrix of the discretization via a reverse Cuthill-McKee
        algorithm. An additional benefit would be an increase in cache
        performance, although this effect is most likely limited for
        high order elements. First determine the map from the global
        element ID to the current storage sequence of ownedElements. ---*/
  mapGlobalElemIDToInd.clear();
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = i;

  /*--- Define the starting positions for the different time levels for
        the internal elements and communication elements. ---*/
  vector<unsigned long> nIntElem(nTimeLevels), nCommElem(nTimeLevels);
  for(unsigned short tLev=0; tLev<nTimeLevels; ++tLev) {
    nIntElem[tLev]  = nVolElemOwnedPerTimeLevel[tLev];
    nCommElem[tLev] = nVolElemOwnedPerTimeLevel[tLev] + nVolElemInternalPerTimeLevel[tLev];
  }

  /*--- Create the graph of local elements. The halo elements are ignored. ---*/
  vector<vector<unsigned long> > neighElem(nVolElemOwned, vector<unsigned long>(0));

  nRankRecv = (int) longRecvBuf.size();
  for(int i=0; i<nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for(long j=0; j<longRecvBuf[i][0]; ++j) {
      unsigned long  globalID  = longRecvBuf[i][indL];
      unsigned short nDOFsGrid = shortRecvBuf[i][indS+3];
      unsigned short nFaces    = shortRecvBuf[i][indS+5];

      map<unsigned long, unsigned long>::iterator MMI = mapGlobalElemIDToInd.find(globalID);
      unsigned long ind = MMI->second;

      indS += 8;
      indL += nDOFsGrid + 2;
      for(unsigned short k=0; k<nFaces; ++k, indS+=3, ++indL) {
        if((longRecvBuf[i][indL] != -1) && (shortRecvBuf[i][indS] == -1)) { // Check for internal owned node.

          MMI = mapGlobalElemIDToInd.find(longRecvBuf[i][indL]);
          if(MMI != mapGlobalElemIDToInd.end()) neighElem[ind].push_back(MMI->second);
        }
      }
    }
  }

  /*--- Sort the neighbors of each element in increasing order. ---*/
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    sort(neighElem[i].begin(), neighElem[i].end());

  /*--- Define the vector, which contains the new numbering of the owned elements
        w.r.t. to the numbering currently stored in ownedElements. Note that signed
        longs are used for this purpose, because the vector is initialized with -1
        to indicate that no new number has been assigned yet. ---*/
  vector<long> oldElemToNewElem(nVolElemOwned, -1);

  /*--- While loop to carry out the renumbering. A while loop is present,
        because the local partition may not be contiguous. ---*/
  unsigned long nElemRenumbered = 0;
  while (nElemRenumbered < nVolElemOwned) {

    /*--- Determine the first element in the list that has not been renumbered. ---*/
    unsigned long indBeg;
    for(indBeg=0; indBeg<nVolElemOwned; ++indBeg)
      if(oldElemToNewElem[indBeg] == -1) break;

    /*--- Determine the time level the element indBeg and end index for the
          element range for renumbering. ---*/
    unsigned short timeLevel = ownedElements[indBeg].GetTimeLevel();

    unsigned long indEnd;
    if( ownedElements[indBeg].GetCommSolution() )
      indEnd = nVolElemOwnedPerTimeLevel[timeLevel+1];
    else
      indEnd = nVolElemOwnedPerTimeLevel[timeLevel]
             + nVolElemInternalPerTimeLevel[timeLevel];

    /*--- Determine the element in the range [indBeg,indEnd) with the least number
          of neighbors that has not been renumbered yet. This is the starting
          element for the current renumbering round. ---*/
    for(unsigned long i=(indBeg+1); i<indEnd; ++i) {
      if((oldElemToNewElem[i] == -1) &&
         (neighElem[i].size() < neighElem[indBeg].size())) indBeg = i;
    }

    /*--- Start of the Reverse Cuthil McKee renumbering. ---*/
    vector<unsigned long> frontElements(1, indBeg);
    while( frontElements.size() ) {

      /*--- Vector, which stores the front for the next round. ---*/
      vector<unsigned long> frontElementsNew;

      /*--- Loop over the elements of the current front. ---*/
      for(unsigned long i=0; i<frontElements.size(); ++i) {

        /*--- Carry out the renumbering for this element. ---*/
        const unsigned long iFront = frontElements[i];

        timeLevel = ownedElements[iFront].GetTimeLevel();

        if( ownedElements[iFront].GetCommSolution() )
          oldElemToNewElem[iFront] = nCommElem[timeLevel]++;
        else
          oldElemToNewElem[iFront] = nIntElem[timeLevel]++;

        /*--- Store the neighbors that have not been renumbered yet in the front
              for the next round. Set its index to -2 to indicate that the element
              is already on the new front. ---*/
        for(unsigned long j=0; j<neighElem[iFront].size(); ++j) {
          if(oldElemToNewElem[neighElem[iFront][j]] == -1) {
            frontElementsNew.push_back(neighElem[iFront][j]);
            oldElemToNewElem[neighElem[iFront][j]] = -2;
          }
        }
      }

      /*--- Update the counter nElemRenumbered. ---*/
      nElemRenumbered += frontElements.size();

      /*--- Sort frontElementsNew in increasing order. ---*/
      sort(frontElementsNew.begin(), frontElementsNew.end());

      /*--- Store the new front elements in frontElements for the next round. ---*/
      frontElements = frontElementsNew;
    }
  }

  if(nElemRenumbered != nVolElemOwned)
    SU2_MPI::Error("Something went wrong in the renumbering", CURRENT_FUNCTION);

  /*--- Determine the final mapping from the global element number to the local
        entry for the owned elements. First clear mapGlobalElemIDToInd before
        it can be used to store its correct content. ---*/
  mapGlobalElemIDToInd.clear();
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = oldElemToNewElem[i];

  /*----------------------------------------------------------------------------*/
  /*--- Step 5: Store the elements, nodes and boundary elements in the data  ---*/
  /*---         structures used by the FEM solver.                           ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Check in parallel mode for empty partitions. If present, print a warning.
        The solver is capable of handling empty partitions, but it may not be
        efficient. ---*/
#ifdef HAVE_MPI
  unsigned long thisPartitionEmpty = nVolElemOwned ? 0 : 1;
  unsigned long nEmptyPartitions = 0;

  SU2_MPI::Reduce(&thisPartitionEmpty, &nEmptyPartitions, 1,
                  MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());

  if(rank == MASTER_NODE && nEmptyPartitions) {
    cout << endl << "         WARNING" << endl;
    cout << "There are " << nEmptyPartitions << " empty partitions present." << endl;
    cout << "SU2 is able to handle this, but it may be inefficient." << endl << endl;
  }
#endif

  /*--- Allocate the memory for the volume elements and the boundaries
        and reserve the memory for the nodes and the surface elements
        of the boundaries. The reserve is done, because for these
        variables push_back is used. ---*/
  volElem.resize(nVolElemTot);
  meshPoints.reserve(nPoint);

  boundaries.resize(nMarker);
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    boundaries[iMarker].markerTag        = config->GetMarker_All_TagBound(iMarker);
    boundaries[iMarker].periodicBoundary = config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY;
    boundaries[iMarker].surfElem.reserve(nElem_Bound[iMarker]);
  }

  /*--- Copy the data from the communication buffers. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /*--- The data for the volume elements. Loop over these elements in the buffer. ---*/
    unsigned long indL = 1, indS = 0, indD = 0;
    for(long j=0; j<longRecvBuf[i][0]; ++j) {

      /*--- Determine the location in volElem where this data must be stored. ---*/
      unsigned long elemID = longRecvBuf[i][indL++];
      map<unsigned long, unsigned long>::iterator MMI = mapGlobalElemIDToInd.find(elemID);
      unsigned long ind = MMI->second;

      /*--- Store the data. ---*/
      volElem[ind].elemIsOwned        = true;
      volElem[ind].rankOriginal       = rank;
      volElem[ind].periodIndexToDonor = -1;

      volElem[ind].VTK_Type  = shortRecvBuf[i][indS++];
      volElem[ind].nPolyGrid = shortRecvBuf[i][indS++];
      volElem[ind].nPolySol  = shortRecvBuf[i][indS++];
      volElem[ind].nDOFsGrid = shortRecvBuf[i][indS++];
      volElem[ind].nDOFsSol  = shortRecvBuf[i][indS++];
      volElem[ind].nFaces    = shortRecvBuf[i][indS++];
      volElem[ind].timeLevel = shortRecvBuf[i][indS++];

      volElem[ind].JacIsConsideredConstant = (bool) shortRecvBuf[i][indS++];

      volElem[ind].elemIDGlobal        = elemID;
      volElem[ind].offsetDOFsSolGlobal = longRecvBuf[i][indL++];

      volElem[ind].nodeIDsGrid.resize(volElem[ind].nDOFsGrid);
      volElem[ind].JacFacesIsConsideredConstant.resize(volElem[ind].nFaces);
      volElem[ind].ElementOwnsFaces.resize(volElem[ind].nFaces);

      for(unsigned short k=0; k<volElem[ind].nDOFsGrid; ++k)
        volElem[ind].nodeIDsGrid[k] = longRecvBuf[i][indL++];

      for(unsigned short k=0; k<volElem[ind].nFaces; ++k) {
        long neighBorID = longRecvBuf[i][indL++];

        ++indS; // At this location the periodic index of the face is stored in
                // shortRecvBuf, which is not stored in volElem.
        volElem[ind].JacFacesIsConsideredConstant[k] = (bool) shortRecvBuf[i][indS++];
        volElem[ind].ElementOwnsFaces[k]             = (bool) shortRecvBuf[i][indS++];

        if(neighBorID == -1)
          volElem[ind].ElementOwnsFaces[k] = true;  // Boundary faces are always owned.
      }

      volElem[ind].lenScale = doubleRecvBuf[i][indD++];
    }

    /*--- The data for the nodes. Loop over these nodes in the buffer and store
          them in meshPoints. ---*/
    unsigned long nNodesThisRank = longRecvBuf[i][indL++];
    for(unsigned long j=0; j<nNodesThisRank; ++j) {
      CPointFEM thisPoint;
      thisPoint.globalID = longRecvBuf[i][indL++];
      thisPoint.periodIndexToDonor = -1;
      for(unsigned short k=0; k<nDim; ++k)
        thisPoint.coor[k] = doubleRecvBuf[i][indD++];

      meshPoints.push_back(thisPoint);
    }

    /*--- The data for the boundary markers. Loop over them. ---*/
    for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

      /*--- Loop over the number of boundary elements for this marker. ---*/
      unsigned long nElemThisRank = longRecvBuf[i][indL++];
      for(unsigned long j=0; j<nElemThisRank; ++j) {

        /*--- Copy the data to an object of the class CSurfaceElementFEM. ---*/
        CSurfaceElementFEM thisSurfElem;

        thisSurfElem.VTK_Type  = shortRecvBuf[i][indS++];
        thisSurfElem.nPolyGrid = shortRecvBuf[i][indS++];
        thisSurfElem.nDOFsGrid = shortRecvBuf[i][indS++];
        const short nDonors    = shortRecvBuf[i][indS++];

        thisSurfElem.volElemID         = longRecvBuf[i][indL++];
        thisSurfElem.boundElemIDGlobal = longRecvBuf[i][indL++];

        thisSurfElem.nodeIDsGrid.resize(thisSurfElem.nDOFsGrid);
        for(unsigned short k=0; k<thisSurfElem.nDOFsGrid; ++k)
          thisSurfElem.nodeIDsGrid[k] = longRecvBuf[i][indL++];

        indL += nDonors;

        /*--- Convert the global volume element ID to the local one.
              It is essential to do this before the sorting. ---*/
        map<unsigned long, unsigned long>::iterator MMI;
        MMI = mapGlobalElemIDToInd.find(thisSurfElem.volElemID);
        thisSurfElem.volElemID = MMI->second;

        /*--- Store the surface element in the data structure for this boundary. ---*/
        boundaries[iMarker].surfElem.push_back(thisSurfElem);
      }
    }
  }

  /*--- Sort meshPoints in increasing order and remove the double entities. ---*/
  sort(meshPoints.begin(), meshPoints.end());
  vector<CPointFEM>::iterator lastPoint = unique(meshPoints.begin(), meshPoints.end());
  meshPoints.erase(lastPoint, meshPoints.end());

  /*--- Clear the contents of the map globalPointIDToLocalInd and fill
        it with the information present in meshPoints. ---*/
  globalPointIDToLocalInd.clear();
  for(unsigned long i=0; i<meshPoints.size(); ++i)
    globalPointIDToLocalInd[meshPoints[i].globalID] = i;

  /*--- All the data from the receive buffers has been copied in the local
        data structures. Release the memory of the receive buffers. To make
        sure that all the memory is deleted, the swap function is used. ---*/
  for(int i=0; i<nRankRecv; ++i) {
    vector<short>().swap(shortRecvBuf[i]);
    vector<long>().swap(longRecvBuf[i]);
    vector<su2double>().swap(doubleRecvBuf[i]);
  }

  /*--- Sort the surface elements of the boundaries in increasing order. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker)
    sort(boundaries[iMarker].surfElem.begin(), boundaries[iMarker].surfElem.end());

  /*----------------------------------------------------------------------------*/
  /*--- Step 6: Obtain the information of the halo elements, which are       ---*/
  /*---         sorted per time level and afterwards per rank, where the     ---*/
  /*---         sequence is determined by the numbering on the sending rank. ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Resize the first index of the send buffers to nRankRecv, because this
        number of messages must be sent back to the sending ranks with halo
        information. ---*/
  nRankRecv = (int) longSecondRecvBuf.size();
  shortSendBuf.resize(nRankRecv);
  longSendBuf.resize(nRankRecv);
  doubleSendBuf.resize(nRankRecv);

#ifdef HAVE_MPI
  /*--- Resize the vector of the communication requests to the number of messages
        to be sent by this rank. Only in parallel node. ---*/
  commReqs.resize(3*nRankRecv);
#endif

  /*--- Loop over the receive buffers to fill the send buffers again. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /*--- Loop over the elements in this receive buffer to determine the local
          index on this rank. Note that also the periodic index must be stored,
          hence use an CUnsignedLong2T for this purpose. As -1 cannot be stored
          for an unsigned long a 1 is added to the periodic transformation. ---*/
    const unsigned long nElemBuf = longSecondRecvBuf[i].size()/2;
    vector<CUnsignedLong2T> elemBuf(nElemBuf);

    for(unsigned long j=0; j<nElemBuf; ++j) {
      const unsigned long j2 = 2*j;

      const unsigned long elemID = longSecondRecvBuf[i][j2];
      map<unsigned long, unsigned long>::iterator MMI = mapGlobalElemIDToInd.find(elemID);
      if(MMI == mapGlobalElemIDToInd.end())
        SU2_MPI::Error("Entry not found in mapGlobalElemIDToInd", CURRENT_FUNCTION);

      elemBuf[j].long0 = MMI->second;
      elemBuf[j].long1 = longSecondRecvBuf[i][j2+1] + 1;
    }

    /*--- Release the memory of the long receive buffer via the swap function
          and sort elemBuf in increasing order. ---*/
    vector<long>().swap(longSecondRecvBuf[i]);

    sort(elemBuf.begin(), elemBuf.end());

    /*--- Store the number of elements in the first element of the long send buffer. ---*/
    longSendBuf[i].push_back(nElemBuf);

    /*--- Vector with node IDs that must be returned to this calling rank.
          Note that also the periodic index must be stored, hence use an
         CUnsignedLong2T for this purpose. ---*/
    vector<CUnsignedLong2T> nodeIDs;

    /* Loop over the elements to fill the send buffers. */
    for(unsigned long j=0; j<nElemBuf; ++j) {

      /*--- Store the global element ID in the long buffer,
            the periodic index in the short send buffer and
            the length scale in the double send buffer. ---*/
      const unsigned long indV = elemBuf[j].long0;
      longSendBuf[i].push_back(volElem[indV].elemIDGlobal);

      const short perIndex = (short) elemBuf[j].long1 -1; // Note the -1, see above.
      shortSendBuf[i].push_back(perIndex);

      doubleSendBuf[i].push_back(volElem[indV].lenScale);

      /*--- Store the other relevant information of this element in the short
            and long communication buffers. Also store the node IDs and the
            periodic transformation in nodeIDs. ---*/
      shortSendBuf[i].push_back(volElem[indV].VTK_Type);
      shortSendBuf[i].push_back(volElem[indV].nPolyGrid);
      shortSendBuf[i].push_back(volElem[indV].nPolySol);
      shortSendBuf[i].push_back(volElem[indV].nDOFsGrid);
      shortSendBuf[i].push_back(volElem[indV].nDOFsSol);
      shortSendBuf[i].push_back(volElem[indV].nFaces);
      shortSendBuf[i].push_back(volElem[indV].timeLevel);

      for(unsigned short k=0; k<volElem[indV].nDOFsGrid; ++k) {
        longSendBuf[i].push_back(volElem[indV].nodeIDsGrid[k]);
        nodeIDs.push_back(CUnsignedLong2T(volElem[indV].nodeIDsGrid[k],
                                          elemBuf[j].long1));
      }

      for(unsigned short k=0; k<volElem[indV].nFaces; ++k)
        shortSendBuf[i].push_back((short) volElem[indV].JacFacesIsConsideredConstant[k]);
    }

    /*--- Sort nodeIDs in increasing order and remove the double entities. ---*/
    sort(nodeIDs.begin(), nodeIDs.end());
    vector<CUnsignedLong2T>::iterator lastNodeID = unique(nodeIDs.begin(), nodeIDs.end());
    nodeIDs.erase(lastNodeID, nodeIDs.end());

    /*--- Add the number of node IDs and the node IDs itself to longSendBuf[i]
          and the periodix index to shortSendBuf. Note again the -1 for the
          periodic index, because an unsigned long cannot represent -1, the
          value for the periodic index when no peridicity is present. ---*/
    longSendBuf[i].push_back(nodeIDs.size());
    for(unsigned long j=0; j<nodeIDs.size(); ++j) {
      longSendBuf[i].push_back(nodeIDs[j].long0);
      shortSendBuf[i].push_back( (short) nodeIDs[j].long1-1);
    }

    /*--- Copy the coordinates to doubleSendBuf. ---*/
    for(unsigned long j=0; j<nodeIDs.size(); ++j) {
      map<unsigned long,unsigned long>::const_iterator LMI;
      LMI = globalPointIDToLocalInd.find(nodeIDs[j].long0);

      if(LMI == globalPointIDToLocalInd.end())
        SU2_MPI::Error("Entry not found in map", CURRENT_FUNCTION);

      unsigned long ind = LMI->second;
      for(unsigned short l=0; l<nDim; ++l)
        doubleSendBuf[i].push_back(meshPoints[ind].coor[l]);
    }

    /*--- Send the communication buffers back to the calling rank.
          Only in parallel mode of course. ---*/
#ifdef HAVE_MPI
    int dest = sourceRank[i];
    SU2_MPI::Isend(shortSendBuf[i].data(), shortSendBuf[i].size(), MPI_SHORT,
                   dest, dest+1, SU2_MPI::GetComm(), &commReqs[3*i]);
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest+2, SU2_MPI::GetComm(), &commReqs[3*i+1]);
    SU2_MPI::Isend(doubleSendBuf[i].data(), doubleSendBuf[i].size(), MPI_DOUBLE,
                   dest, dest+3, SU2_MPI::GetComm(), &commReqs[3*i+2]);
#endif
  }

  /*--- Resize the first index of the receive buffers to nRankSendHaloInfo,
        such that the requested halo information can be received. ---*/
  nRankSend = nRankSendHaloInfo;
  shortRecvBuf.resize(nRankSend);
  longRecvBuf.resize(nRankSend);
  doubleRecvBuf.resize(nRankSend);

  /*--- Resize the vector to store the ranks from which the message came. ---*/
  sourceRank.resize(nRankSend);

  /*--- Receive the communication data from the correct ranks. Make a distinction
        between parallel and sequential mode. ---*/
#ifdef HAVE_MPI

  /*--- Parallel mode. Loop over the number of ranks from which I receive data
        in the return communication, i.e. nRankSend. ---*/
  for(int i=0; i<nRankSend; ++i) {

    /*--- Block until a message with shorts arrives from any processor.
          Determine the source and the size of the message. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_SHORT, &sizeMess);

    /*--- Allocate the memory for the short receive buffer and receive the message. ---*/
    shortRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(shortRecvBuf[i].data(), sizeMess, MPI_SHORT,
                  sourceRank[i], rank+1, SU2_MPI::GetComm(), &status);

    /*--- Block until the corresponding message with longs arrives, determine
          its size, allocate the memory and receive the message. ---*/
    SU2_MPI::Probe(sourceRank[i], rank+2, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);
    longRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(longRecvBuf[i].data(), sizeMess, MPI_LONG,
                  sourceRank[i], rank+2, SU2_MPI::GetComm(), &status);

    /*--- Idem for the message with doubles. ---*/
    SU2_MPI::Probe(sourceRank[i], rank+3, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_DOUBLE, &sizeMess);
    doubleRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(doubleRecvBuf[i].data(), sizeMess, MPI_DOUBLE,
                  sourceRank[i], rank+3, SU2_MPI::GetComm(), &status);
  }

  /*--- Complete the non-blocking sends and synchronize the ranks to
        avoid problems, because wild cards have been used. ---*/
  SU2_MPI::Waitall(3*nRankRecv, commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffers. Note that nRankSend is at most 1. ---*/
  for(int i=0; i<nRankSend; ++i) {
    sourceRank[i]    = i;
    shortRecvBuf[i]  = shortSendBuf[i];
    longRecvBuf[i]   = longSendBuf[i];
    doubleRecvBuf[i] = doubleSendBuf[i];
  }

#endif

  /*--- Release the memory of the send buffers. To make sure that all
        the memory is deleted, the swap function is used. ---*/
  for(int i=0; i<nRankRecv; ++i) {
    vector<short>().swap(shortSendBuf[i]);
    vector<long>().swap(longSendBuf[i]);
    vector<su2double>().swap(doubleSendBuf[i]);
  }

  /*----------------------------------------------------------------------------*/
  /*--- Step 7: Build the layer of halo elements from the information in the ---*/
  /*---         receive buffers shortRecvBuf, longRecvBuf and doubleRecvBuf. ---*/
  /*----------------------------------------------------------------------------*/

  /*--- The halo elements must be sorted first based on the time level, followed
        by the rank and finally the index in the receive buffers (which represents
        the sequence on the sending rank. This sorting can be accomplished by
        using a vector of CLong3T. The contents of this vector is build below. ---*/
  vector<CLong3T> haloElemInfo;
  haloElemInfo.reserve(nVolElemTot - nVolElemOwned);

  for(int i=0; i<nRankSend; ++i) {
    unsigned long indS = 0;

    for(long j=0; j<longRecvBuf[i][0]; ++j) {
      const unsigned short nFaces = shortRecvBuf[i][indS+6];
      haloElemInfo.push_back(CLong3T(shortRecvBuf[i][indS+7], sourceRank[i], j));
      indS += nFaces + 8;
    }
  }

  sort(haloElemInfo.begin(), haloElemInfo.end());

  /*--- Loop over the receive buffers to store the information of the
        halo elements and the halo points. ---*/
  vector<CPointFEM> haloPoints;
  for(int i=0; i<nRankSend; ++i) {

    /*--- Initialization of the indices in the communication buffers. ---*/
    unsigned long indL = 1, indS = 0, indD = 0;

    /*--- Loop over the halo elements received from this rank. ---*/
    for(long j=0; j<longRecvBuf[i][0]; ++j) {

      /*--- Create an object of CLong3T and search for its position
            in haloElemInfo to determine the index in volElem, where
            this element is stored. ---*/
      CLong3T thisElem(shortRecvBuf[i][indS+7], sourceRank[i], j);
      vector<CLong3T>::iterator low;
      low = lower_bound(haloElemInfo.begin(), haloElemInfo.end(), thisElem);

      unsigned long indV = low - haloElemInfo.begin();
      indV += nVolElemOwned;

      /*--- Retrieve the data from the communication buffers. ---*/
      volElem[indV].elemIDGlobal = longRecvBuf[i][indL++];
      volElem[indV].rankOriginal = sourceRank[i];

      volElem[indV].periodIndexToDonor = shortRecvBuf[i][indS++];
      volElem[indV].VTK_Type           = shortRecvBuf[i][indS++];
      volElem[indV].nPolyGrid          = shortRecvBuf[i][indS++];
      volElem[indV].nPolySol           = shortRecvBuf[i][indS++];
      volElem[indV].nDOFsGrid          = shortRecvBuf[i][indS++];
      volElem[indV].nDOFsSol           = shortRecvBuf[i][indS++];
      volElem[indV].nFaces             = shortRecvBuf[i][indS++];
      volElem[indV].timeLevel          = shortRecvBuf[i][indS++];

      volElem[indV].nodeIDsGrid.resize(volElem[indV].nDOFsGrid);
      for(unsigned short k=0; k<volElem[indV].nDOFsGrid; ++k)
        volElem[indV].nodeIDsGrid[k] = longRecvBuf[i][indL++];

      volElem[indV].JacFacesIsConsideredConstant.resize(volElem[indV].nFaces);
      for(unsigned short k=0; k<volElem[indV].nFaces; ++k)
        volElem[indV].JacFacesIsConsideredConstant[k] = (bool) shortRecvBuf[i][indS++];

      /*--- Give the member variables that are not obtained via communication their
            values. Some of these variables are not used for halo elements. ---*/
      volElem[indV].elemIsOwned             = false;
      volElem[indV].JacIsConsideredConstant = false;

      /*--- Halo elements do not own a face per definition. ---*/
      volElem[indV].ElementOwnsFaces.assign(volElem[indV].nFaces, false);

      /*--- Get the length scale from the double receive buffer. ---*/
      volElem[indV].lenScale = doubleRecvBuf[i][indD++];
    }

    /*--- Store the information of the points in haloPoints. ---*/
    const long nPointsThisRank = longRecvBuf[i][indL++];
    for(long j=0; j<nPointsThisRank; ++j) {
      CPointFEM thisPoint;
      thisPoint.globalID           = longRecvBuf[i][indL++];
      thisPoint.periodIndexToDonor = shortRecvBuf[i][indS++];
      for(unsigned short l=0; l<nDim; ++l)
        thisPoint.coor[l] = doubleRecvBuf[i][indD++];

      haloPoints.push_back(thisPoint);
    }

    /*--- The communication buffers from this rank are not needed anymore.
          Delete them using the swap function. ---*/
    vector<short>().swap(shortRecvBuf[i]);
    vector<long>().swap(longRecvBuf[i]);
    vector<su2double>().swap(doubleRecvBuf[i]);
  }

  /*--- Remove the duplicate entries from haloPoints. ---*/
  sort(haloPoints.begin(), haloPoints.end());
  lastPoint = unique(haloPoints.begin(), haloPoints.end());
  haloPoints.erase(lastPoint, haloPoints.end());

  /*--- Initialization of some variables to sort the halo points. ---*/
  Global_nPoint = geometry->GetGlobal_nPoint();
  unsigned long InvalidPointID = Global_nPoint + 10;
  short         InvalidPerInd  = SHRT_MAX;

  /*--- Search for the nonperiodic halo points in the local points to see
        if these points are already stored on this rank. If this is the
        case invalidate this halo and decrease the number of halo points.
        Afterwards remove the invalid halos from the vector. ---*/
  unsigned long nHaloPoints = haloPoints.size();
  for(unsigned long i=0; i<haloPoints.size(); ++i) {
    if(haloPoints[i].periodIndexToDonor != -1) break;  // Test for nonperiodic.

    if( binary_search(meshPoints.begin(), meshPoints.end(), haloPoints[i]) ) {
      haloPoints[i].globalID           = InvalidPointID;
      haloPoints[i].periodIndexToDonor = InvalidPerInd;
      --nHaloPoints;
    }
  }

  sort(haloPoints.begin(), haloPoints.end());
  haloPoints.resize(nHaloPoints);

  /*--- Increase the capacity of meshPoints, such that the halo points can be
        stored in there as well. Note that in case periodic points are present
        this is an upper bound. Add the non-periodic halo points to meshPoints. ---*/
  meshPoints.reserve(meshPoints.size() + nHaloPoints);

  for(unsigned long i=0; i<haloPoints.size(); ++i) {
    if(haloPoints[i].periodIndexToDonor != -1) break;  // Test for nonperiodic.

    meshPoints.push_back(haloPoints[i]);
  }

  /*--- Create a map from the global point ID and periodic index to the local
        index in the vector meshPoints. First store the points already present
        in meshPoints. ---*/
  map<CUnsignedLong2T, unsigned long> mapGlobalPointIDToInd;
  for(unsigned long i=0; i<meshPoints.size(); ++i) {
    CUnsignedLong2T globIndAndPer;
    globIndAndPer.long0 = meshPoints[i].globalID;
    globIndAndPer.long1 = meshPoints[i].periodIndexToDonor+1;  // Note the +1 again.

    mapGlobalPointIDToInd[globIndAndPer] = i;
  }

  /*--- Convert the global indices in the boundary connectivities to local ones.
        Note that the volume ID's already contain the local number. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    for(unsigned long i=0; i<boundaries[iMarker].surfElem.size(); ++i) {

      /*--- Convert the global node ID's to local values. Note that for these node
            ID's no periodic transformation can be present. ---*/
      for(unsigned short j=0; j<boundaries[iMarker].surfElem[i].nDOFsGrid; ++j) {
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
  for(unsigned long iLow=0; iLow<haloPoints.size(); ) {

    /*--- Determine the upper index for this periodic transformation. ---*/
    unsigned long iUpp;
    for(iUpp=iLow+1; iUpp<haloPoints.size(); ++iUpp)
      if(haloPoints[iUpp].periodIndexToDonor != haloPoints[iLow].periodIndexToDonor) break;

    /*--- Check for a true periodic index. ---*/
    short perIndex = haloPoints[iLow].periodIndexToDonor;
    if(perIndex != -1) {

      /*--- Easier storage of the surface elements. ---*/
      vector<CSurfaceElementFEM> &surfElem = boundaries[perIndex].surfElem;

      /*--- In the loop below the coordinates of the points of this local
            periodic boundary as well as a matching tolerance are determined.
            A vector of point ID's is also created, which is needed later on
            when it is checked whether or not a matching point is already
            stored in meshPoints. ---*/
      vector<long>          indInPoints(meshPoints.size(), -1);
      vector<unsigned long> IDsPoints;
      vector<su2double>     coordPoints;
      vector<su2double>     tolPoints;

      for(unsigned long j=0; j<surfElem.size(); ++j) {

        /*--- Determine the tolerance for equal points, which is a small value
              times the length scale of the adjacent volume element. ---*/
        const su2double tolElem = 1.e-2*volElem[surfElem[j].volElemID].lenScale;

        /*--- Loop over the nodes of this surface grid and update the points on
              this periodic boundary. ---*/
        for(unsigned short k=0; k<surfElem[j].nDOFsGrid; ++k) {
          unsigned long nn = surfElem[j].nodeIDsGrid[k];

          if(indInPoints[nn] == -1) {

            /*--- Point is not stored yet in pointsBoundary. Do so now. ---*/
            indInPoints[nn] = IDsPoints.size();
            IDsPoints.push_back(nn);
            tolPoints.push_back(tolElem);

            for(unsigned short l=0; l<nDim; ++l)
              coordPoints.push_back(meshPoints[nn].coor[l]);
          }
          else {

            /*--- Point is already stored. Update the tolerance. ---*/
            nn = indInPoints[nn];
            tolPoints[nn] = min(tolPoints[nn], tolElem);
          }
        }
      }

      /*--- Create a local ADT of the points on the periodic boundary. ---*/
      CADTPointsOnlyClass periodicADT(nDim, IDsPoints.size(), coordPoints.data(),
                                      IDsPoints.data(), false);

      /*--- Get the data for the periodic transformation to the donor. ---*/
      auto center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(perIndex));
      auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(perIndex));
      auto trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(perIndex));

      /*--- Compute the rotation matrix and translation vector for the
            transformation from the donor. This is the transpose of the
            transformation to the donor. ---*/

      /*--- Store (center-trans) as it is constant and will be added on. ---*/
      su2double translation[] = {center[0] - trans[0],
                                 center[1] - trans[1],
                                 center[2] - trans[2]};

      /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
      su2double theta = angles[0];
      su2double phi   = angles[1];
      su2double psi   = angles[2];

      su2double cosTheta = cos(theta), cosPhi = cos(phi), cosPsi = cos(psi);
      su2double sinTheta = sin(theta), sinPhi = sin(phi), sinPsi = sin(psi);

      /*--- Compute the rotation matrix. Note that the implicit
            ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
      su2double rotMatrix[3][3];
      rotMatrix[0][0] =  cosPhi*cosPsi;
      rotMatrix[0][1] =  cosPhi*sinPsi;
      rotMatrix[0][2] = -sinPhi;

      rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
      rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
      rotMatrix[1][2] = sinTheta*cosPhi;

      rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
      rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
      rotMatrix[2][2] = cosTheta*cosPhi;

      /*--- Loop over the halo points for this periodic transformation. ---*/
      for(unsigned long i=iLow; i<iUpp; ++i) {

        /*--- Apply the periodic transformation to the coordinates
              stored in this halo point. ---*/
        su2double dx =             haloPoints[i].coor[0] - center[0];
        su2double dy =             haloPoints[i].coor[1] - center[1];
        su2double dz = nDim == 3 ? haloPoints[i].coor[2] - center[2] : su2double(0.0);

        haloPoints[i].coor[0] = rotMatrix[0][0]*dx + rotMatrix[0][1]*dy
                              + rotMatrix[0][2]*dz + translation[0];
        haloPoints[i].coor[1] = rotMatrix[1][0]*dx + rotMatrix[1][1]*dy
                              + rotMatrix[1][2]*dz + translation[1];
        haloPoints[i].coor[2] = rotMatrix[2][0]*dx + rotMatrix[2][1]*dy
                              + rotMatrix[2][2]*dz + translation[2];

        /*--- Search for the nearest coordinate in the ADT. ---*/
        su2double dist;
        unsigned long pointID;
        int rankID;

        periodicADT.DetermineNearestNode(haloPoints[i].coor, dist,
                                         pointID, rankID);

        /*--- Check whether the distance is less equal to the tolerance for
              a matching point. ---*/
        const unsigned long nn = indInPoints[pointID];
        if(dist <= tolPoints[nn]) {

          /*--- The distance to the nearest point is less than the tolerance,
                hence this periodically transformed point is present on the
                boundary. Store it as such in the map mapGlobalPointIDToInd. ---*/
          CUnsignedLong2T globIndAndPer;
          globIndAndPer.long0 = haloPoints[i].globalID;
          globIndAndPer.long1 = haloPoints[i].periodIndexToDonor+1;  // Note the +1 again.

          mapGlobalPointIDToInd[globIndAndPer] = pointID;
        }
        else {

          /*--- The distance to the nearest point is larger than the tolerance,
                hence this periodically transformed point is not present yet on
                this rank. Store it in the mapping to the local points and
                create it in meshPoints. ---*/
          CUnsignedLong2T globIndAndPer;
          globIndAndPer.long0 = haloPoints[i].globalID;
          globIndAndPer.long1 = haloPoints[i].periodIndexToDonor+1;  // Note the +1 again.

          mapGlobalPointIDToInd[globIndAndPer] = meshPoints.size();
          meshPoints.push_back(haloPoints[i]);
        }
      }
    }

    /*--- Set iLow to iUpp for the next periodic transformation. ---*/
    iLow = iUpp;
  }

  /*--- Convert the global node numbering in the elements to a local numbering and
        determine the value of factTimeLevel. This is the number of local time steps
        of the element relative to the largest time step of an element in the mesh.
        This value can only differ from 1 when time accurate local time stepping is
        used. ---*/
  for(unsigned long i=0; i<nVolElemTot; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsGrid; ++j) {
      CUnsignedLong2T searchItem(volElem[i].nodeIDsGrid[j],
                                 volElem[i].periodIndexToDonor+1); // Again the +1.
      map<CUnsignedLong2T, unsigned long>::const_iterator LLMI;
      LLMI = mapGlobalPointIDToInd.find(searchItem);
      volElem[i].nodeIDsGrid[j] = LLMI->second;
    }

    const unsigned short diffTimeLevel = nTimeLevels - 1 - volElem[i].timeLevel;
    volElem[i].factTimeLevel = pow(2, diffTimeLevel);
  }

  /*--- Determine the number of halo elements per time level in cumulative
        storage format. ---*/
  nVolElemHaloPerTimeLevel.assign(nTimeLevels+1, 0);
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i)
    ++nVolElemHaloPerTimeLevel[volElem[i].timeLevel+1];

  nVolElemHaloPerTimeLevel[0] = nVolElemOwned;
  for(unsigned short i=0; i<nTimeLevels; ++i)
    nVolElemHaloPerTimeLevel[i+1] += nVolElemHaloPerTimeLevel[i];
}

void CMeshFEM_DG::CreateFaces(CConfig *config) {

  /*--- The master node writes a message. ---*/
  if(rank == MASTER_NODE) cout << "Creating face information." << endl;

  /*---------------------------------------------------------------------------*/
  /*--- Step 1: Determine the faces of the locally stored part of the grid. ---*/
  /*---------------------------------------------------------------------------*/

  /*--- Loop over the volume elements stored on this rank, including the halos. ---*/
  vector<CFaceOfElement> localFaces;

  for(unsigned long k=0; k<nVolElemTot; ++k) {

    /*--- Determine the corner points of all the faces of this element. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long  faceConn[6][4];

    volElem[k].GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Loop over the faces of this element, set the appropriate information,
          create a unique numbering and add the faces to localFaces. ---*/
    for(unsigned short i=0; i<nFaces; ++i) {
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for(unsigned short j=0; j<nPointsPerFace[i]; ++j)
        thisFace.cornerPoints[j] = faceConn[i][j];

      /*--- Copy the data from the volume element. ---*/
      thisFace.elemID0       =  k;
      thisFace.nPolyGrid0    =  volElem[k].nPolyGrid;
      thisFace.nPolySol0     =  volElem[k].nPolySol;
      thisFace.nDOFsElem0    =  volElem[k].nDOFsSol;
      thisFace.elemType0     =  volElem[k].VTK_Type;
      thisFace.faceID0       =  i;
      thisFace.faceIndicator = -2;   // Initialized to an invalid face.

      thisFace.JacFaceIsConsideredConstant = volElem[k].JacFacesIsConsideredConstant[i];
      thisFace.elem0IsOwner                = volElem[k].ElementOwnsFaces[i];

      /*--- Renumber the corner points of the face, but keep the orientation.
            Afterwards, add it to localFaces. ---*/
      thisFace.CreateUniqueNumberingWithOrientation();

      localFaces.push_back(thisFace);
    }
  }

  /*--- Sort the the local faces in increasing order. ---*/
  sort(localFaces.begin(), localFaces.end());

  /*--- Loop over the faces to merge the matching faces. ---*/
  for(unsigned long i=1; i<localFaces.size(); ++i) {

    /*--- Check for a matching face with the previous face in the vector.
          Note that the == operator only checks the node IDs. ---*/
    if(localFaces[i] == localFaces[i-1]) {

      /*--- Faces are matching. Check if it should be kept, i.e. if one
            of the elements owns the face. ---*/
      if(localFaces[i].elem0IsOwner || localFaces[i-1].elem0IsOwner) {

        /*--- Store the data for this matching face in faces[i-1]. ---*/
        if(localFaces[i].elemID0 < nVolElemTot) {
          localFaces[i-1].elemID0      = localFaces[i].elemID0;
          localFaces[i-1].nPolyGrid0   = localFaces[i].nPolyGrid0;
          localFaces[i-1].nPolySol0    = localFaces[i].nPolySol0;
          localFaces[i-1].nDOFsElem0   = localFaces[i].nDOFsElem0;
          localFaces[i-1].elemType0    = localFaces[i].elemType0;
          localFaces[i-1].faceID0      = localFaces[i].faceID0;
          localFaces[i-1].elem0IsOwner = localFaces[i].elem0IsOwner;
        }
        else {
          localFaces[i-1].elemID1    = localFaces[i].elemID1;
          localFaces[i-1].nPolyGrid1 = localFaces[i].nPolyGrid1;
          localFaces[i-1].nPolySol1  = localFaces[i].nPolySol1;
          localFaces[i-1].nDOFsElem1 = localFaces[i].nDOFsElem1;
          localFaces[i-1].elemType1  = localFaces[i].elemType1;
          localFaces[i-1].faceID1    = localFaces[i].faceID1;
        }

        /*--- Adapt the boolean to indicate whether or not the face has a constant
              Jacobian of the transformation, although in principle this info
              should be the same for both faces. --- */
        if(localFaces[i-1].JacFaceIsConsideredConstant !=
           localFaces[i].JacFaceIsConsideredConstant)
          localFaces[i-1].JacFaceIsConsideredConstant = false;

        /*--- Set this face indicator to -1 to indicate an internal face
              and set elem0IsOwner for localFaces[i] to false. ---*/
        localFaces[i-1].faceIndicator = -1;
        localFaces[i].elem0IsOwner    = false;
      }
    }
  }

  /*--- Loop over the boundary markers and its boundary elements to search for
        the corresponding faces in localFaces. These faces should be found.
        Note that periodic boundaries are skipped, because these are treated
        via the halo elements, which are already in place. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    if( !boundaries[iMarker].periodicBoundary ) {

      for(unsigned long k=0; k<boundaries[iMarker].surfElem.size(); ++k) {

        /*--- Determine the corner points of the face of this element. ---*/
        unsigned short nPointsPerFace;
        unsigned long  faceConn[4];

        boundaries[iMarker].surfElem[k].GetCornerPointsFace(nPointsPerFace, faceConn);

        /*--- Create an object of CFaceOfElement to carry out the search. ---*/
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace;
        for(unsigned short j=0; j<nPointsPerFace; ++j)
          thisFace.cornerPoints[j] = faceConn[j];
        thisFace.CreateUniqueNumberingWithOrientation();

        /*--- Search for thisFace in localFaces. It must be found. ---*/
        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        bool thisFaceFound = false;
        if(low != localFaces.end()) {
          if( !(thisFace < *low) ) thisFaceFound = true;
        }

        if( thisFaceFound ) {
          low->faceIndicator = iMarker;

          /*--- A few additional checks. ---*/
          bool side0IsBoundary = low->elemID0 < nVolElemTot;
          unsigned long elemID = side0IsBoundary ? low->elemID0    : low->elemID1;
          unsigned short nPoly = side0IsBoundary ? low->nPolyGrid0 : low->nPolyGrid1;

          if(elemID != boundaries[iMarker].surfElem[k].volElemID ||
             nPoly  != boundaries[iMarker].surfElem[k].nPolyGrid)
            SU2_MPI::Error(string("Element ID and/or polynomial degree do not match ") +
                           string("for this boundary element. This should not happen."),
                           CURRENT_FUNCTION);
        }
        else
          SU2_MPI::Error("Boundary face not found in localFaces. This should not happen.",
                         CURRENT_FUNCTION);
      }
    }
  }

  /*--- It is possible that owned non-matching faces are present in the list.
        These faces are indicated by an owned face and a faceIndicator of -2.
        To avoid that these faces are removed afterwards, set their
        faceIndicator to -1. ---*/
  for(unsigned long i=0; i<localFaces.size(); ++i) {
    if(localFaces[i].faceIndicator == -2 && localFaces[i].elem0IsOwner)
      localFaces[i].faceIndicator = -1;
  }

  /*--- Remove the invalid faces. This is accomplished by giving the face four
        points and global node ID's that are larger than the largest local point
        ID in the grid. In this way the sorting operator puts these faces at the
        end of the vector, see also the < operator of CFaceOfElement. ---*/
  unsigned long nFacesLoc = localFaces.size();
  for(unsigned long i=0; i<localFaces.size(); ++i) {
    if(localFaces[i].faceIndicator == -2) {
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

  /*--- Sort localFaces again, but now such that the boundary faces are numbered
        first, followed by the matching faces and at the end of localFaces the
        non-matching faces are stored. In order to carry out this sorting the
        functor CSortFaces is used for comparison. Within the categories
        the sorting depends on the time levels of the adjacent volume elements
        of the face. ---*/
  sort(localFaces.begin(), localFaces.end(),
       CSortFaces(nVolElemOwned, nVolElemTot, volElem.data()));

  /*--- Carry out a possible swap of side 0 and side 1 of the faces and renumber
        the nodes of the faces, such that the orientation of the faces matches
        the orientation of the face of the element on side 0. ---*/
  for(unsigned long i=0; i<localFaces.size(); ++i) {
    localFaces[i].SwapSidesIfNeeded(nVolElemOwned, nVolElemTot);
    localFaces[i].MatchOrientationElemSide0(volElem);
  }

  /*--- Determine the number of matching and non-matching internal faces.
        For the matching faces, determine these numbers per time level
        and also make a distinction between internal faces and faces that
        involve a halo element. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
  nMatchingFacesInternal.assign(nTimeLevels+1, 0);
  nMatchingFacesWithHaloElem.assign(nTimeLevels+1, 0);

  unsigned long nNonMatchingFaces = 0;
  for(unsigned long i=0; i<localFaces.size(); ++i) {
    if(localFaces[i].faceIndicator == -1) {
      const unsigned long e0 = localFaces[i].elemID0;
      const unsigned long e1 = localFaces[i].elemID1;

      if(e1 < nVolElemTot) {
        const unsigned short timeLevel = min(volElem[e0].timeLevel,
                                             volElem[e1].timeLevel);
        if(max(e0,e1) < nVolElemOwned) ++nMatchingFacesInternal[timeLevel+1];
        else                           ++nMatchingFacesWithHaloElem[timeLevel+1];
      }
      else ++nNonMatchingFaces;
    }
  }

  if( nNonMatchingFaces ) {
    ostringstream message;
    message << nNonMatchingFaces << " non-matching internal faces found.\n"
            << "This is not supported yet." << endl;
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /*--- Put nMatchingFacesInternal and nMatchingFacesWithHaloElem in
        cumulative storage format. ---*/
  for(unsigned short i=0; i<nTimeLevels; ++i)
    nMatchingFacesInternal[i+1] += nMatchingFacesInternal[i];

  nMatchingFacesWithHaloElem[0] = nMatchingFacesInternal[nTimeLevels];
  for(unsigned short i=0; i<nTimeLevels; ++i)
    nMatchingFacesWithHaloElem[i+1] += nMatchingFacesWithHaloElem[i];

  /*---------------------------------------------------------------------------*/
  /*--- Step 3: Create the local face based data structure for the internal ---*/
  /*---         faces. These are needed for the computation of the surface  ---*/
  /*---         integral in DG-FEM.                                         ---*/
  /*---------------------------------------------------------------------------*/

  /*--- Allocate the memory for the matching faces. ---*/
  matchingFaces.resize(nMatchingFacesWithHaloElem[nTimeLevels]);

  /*--- Initialize the boolean that an element is adjacent to a
        lower time level to false. ---*/
  for(unsigned long i=0; i<nVolElemTot; ++i)
    volElem[i].elemAdjLowerTimeLevel = false;

  /*--- Initialize the matching faces and flag the elements that share
        one or more faces with elements from a lower time level. ---*/
  unsigned long ii = 0;
  for(unsigned long i=0; i<localFaces.size(); ++i) {

    /*--- Check for a matching internal face. ---*/
    if(localFaces[i].faceIndicator == -1 && localFaces[i].elemID1 < nVolElemTot) {

      /*--- Abbreviate the adjacent element ID's and store them in matchingFaces.
            Update the counter ii afterwards. ---*/
      const unsigned long e0 = localFaces[i].elemID0;
      const unsigned long e1 = localFaces[i].elemID1;

      matchingFaces[ii].elemID0 = e0;
      matchingFaces[ii].elemID1 = e1;
      ++ii;

      /*--- Compare the time levels of the adjacent elements and set
            elemAdjLowerTimeLevel accordingly. ---*/
      if(volElem[e0].timeLevel > volElem[e1].timeLevel) volElem[e0].elemAdjLowerTimeLevel = true;
      if(volElem[e1].timeLevel > volElem[e0].timeLevel) volElem[e1].elemAdjLowerTimeLevel = true;
    }
  }

  /*--- Determine the list of elements per time level, which share one or more
        faces with elements of the lower time level. Make a distinction between
        owned and halo elements. ---*/
  ownedElemAdjLowTimeLevel.resize(nTimeLevels);
  haloElemAdjLowTimeLevel.resize(nTimeLevels);

  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    if( volElem[i].elemAdjLowerTimeLevel )
      ownedElemAdjLowTimeLevel[volElem[i].timeLevel].push_back(i);
  }

  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i) {
    if( volElem[i].elemAdjLowerTimeLevel )
      haloElemAdjLowTimeLevel[volElem[i].timeLevel].push_back(i);
  }

  /*---------------------------------------------------------------------------*/
  /*--- Step 4: Determine the number boundary surfaces per time level for   ---*/
  /*---         the physical boundaries.                                    ---*/
  /*---------------------------------------------------------------------------*/

  /*--- Loop over the boundary markers. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /*--- Initialize the number of surface elements per time level to zero. ---*/
    boundaries[iMarker].nSurfElem.assign(nTimeLevels+1, 0);

    /*--- The periodic boundaries are skipped, because these are not physical
          boundaries and are treated via the halo elements. These have already
          been created. ---*/
    if( !boundaries[iMarker].periodicBoundary ) {

      /*--- Easier storage of the surface elements for this boundary and
            loop over them. ---*/
      vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;

      for(unsigned long i=0; i<surfElem.size(); ++i) {

        /*--- Determine the time level of the adjacent element and increment
              the number of surface elements for this time level. The +1 is there,
              because this vector will be put in cumulative storage format
              afterwards. ---*/
        const unsigned short timeLevel = volElem[surfElem[i].volElemID].timeLevel;
        ++boundaries[iMarker].nSurfElem[timeLevel+1];
      }

      /*--- Put boundaries[iMarker].nSurfElem in cumulative storage. ---*/
      for(unsigned short i=0; i<nTimeLevels; ++i)
        boundaries[iMarker].nSurfElem[i+1] += boundaries[iMarker].nSurfElem[i];
    }
  } 

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Store for the volume elements the corresponding internal   ---*/
  /*---         faces and boundary faces.                                  ---*/
  /*--------------------------------------------------------------------------*/

  /* The internal faces. */
  for(unsigned long i=0; i<matchingFaces.size(); ++i) {
    volElem[matchingFaces[i].elemID0].internalFaceIDs.push_back(i);
    volElem[matchingFaces[i].elemID1].internalFaceIDs.push_back(i);
  }

  /* The boundary faces. */
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
    for(unsigned long i=0; i<surfElem.size(); ++i)
      volElem[surfElem[i].volElemID].boundaryFaceIDs.push_back(CUnsignedLong2T(iMarker,i));
  }

  /*---------------------------------------------------------------------------*/
  /*--- Step 6: Create the standard elements for the faces. It is called    ---*/
  /*---         from within this function, because the data of localFaces   ---*/
  /*---         is needed to construct these standard elements.             ---*/
  /*---------------------------------------------------------------------------*/

  if (rank == MASTER_NODE) cout << "Creating standard face elements." << endl;
  CreateStandardFaces(config, localFaces);
}

void CMeshFEM_DG::CreateStandardVolumeElements(CConfig *config) {

  /*--- Check if an incompressible solver is used. ---*/
  MAIN_SOLVER Kind_Solver = config->GetKind_Solver();
  const bool incompressible = (Kind_Solver == MAIN_SOLVER::FEM_INC_EULER) ||
                              (Kind_Solver == MAIN_SOLVER::FEM_INC_NAVIER_STOKES) ||
                              (Kind_Solver == MAIN_SOLVER::FEM_INC_RANS) ||
                              (Kind_Solver == MAIN_SOLVER::FEM_INC_LES);

  /*--- Define the number of variables per DOF to be solved, depending on
        the situation. If MKL is used this number must be known, such that
        a jitted gemm call can be constructed. Otherwise this is not necessary
        and this variable can change during runtime. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  const unsigned short nSolVar = incompressible ? nDim : nDim+2;
#else
  const unsigned short nSolVar = 1;
#endif

  /*--- Determine the number of variables for which memory must be
        allocated for the working variables. ---*/
  const unsigned short nAllocVar = incompressible ? nDim : nDim+2;

  /*--- Vector of four unsigned shorts per entity to determine the
        different element types in the locally stored volume elements. ---*/
  vector<CUnsignedShort4T> elemTypesGrid, elemTypesSol;
  elemTypesGrid.reserve(nVolElemTot);
  if( incompressible ) elemTypesSol.reserve(2*nVolElemTot);
  else                 elemTypesSol.reserve(nVolElemTot);

  /*--- Loop over the number of elements. ---*/
  for(unsigned long i=0; i<nVolElemTot; ++i) {

    /*--- Determine the polynomial order that must be integrated exactly. ---*/
    const unsigned short orderExact = config->GetOrderExactIntegrationFEM(volElem[i].nPolySol,
                                                                          volElem[i].JacIsConsideredConstant);

    /*--- Store the required information in elemTypesGrid and elemTypesSol. ---*/
    elemTypesGrid.push_back(CUnsignedShort4T(volElem[i].VTK_Type, volElem[i].nPolyGrid,
                                             volElem[i].nPolySol, orderExact));

    elemTypesSol.push_back(CUnsignedShort4T(volElem[i].VTK_Type, volElem[i].nPolySol, orderExact, nSolVar));
    if( incompressible )
      elemTypesSol.push_back(CUnsignedShort4T(volElem[i].VTK_Type, volElem[i].nPolySol-1, orderExact, 1));
  }

  /*--- Sort elemTypesGrid and elemTypesSol in increasing order
        and remove the multiple entities. ---*/
  vector<CUnsignedShort4T>::iterator lastEntry4T;
  sort(elemTypesGrid.begin(), elemTypesGrid.end());
  lastEntry4T = unique(elemTypesGrid.begin(), elemTypesGrid.end());
  elemTypesGrid.erase(lastEntry4T, elemTypesGrid.end());

  sort(elemTypesSol.begin(), elemTypesSol.end());
  lastEntry4T = unique(elemTypesSol.begin(), elemTypesSol.end());
  elemTypesSol.erase(lastEntry4T, elemTypesSol.end());

  /*--- Call the functions to actually create the standard elements. ---*/
  CreateStandardVolumeElementsGrid(elemTypesGrid, config->GetKind_FEM_GridDOFsLocation());
  CreateStandardVolumeElementsSolution(elemTypesSol, nAllocVar,
                                       config->GetKind_FEM_GridDOFsLocation());

  /*--- Loop again over the volume elements to set the pointers to the appropriate
        standard elements. ---*/
  for(unsigned long i=0; i<nVolElemTot; ++i) {

    const unsigned short orderExact = config->GetOrderExactIntegrationFEM(volElem[i].nPolySol,
                                                                          volElem[i].JacIsConsideredConstant);

    /*--- Create the CCUnsignedShort4T types to search the grid and solution elements. ---*/
    const CUnsignedShort4T gridType(volElem[i].VTK_Type, volElem[i].nPolyGrid,
                                    volElem[i].nPolySol, orderExact);
    const CUnsignedShort4T solType( volElem[i].VTK_Type, volElem[i].nPolySol,
                                    orderExact, nSolVar);
    const CUnsignedShort4T pType(   volElem[i].VTK_Type, volElem[i].nPolySol-1,
                                    orderExact, 1);

    /*--- Set the pointer for the standard element of the grid. ---*/
    unsigned long j;
    for(j=0; j<elemTypesGrid.size(); ++j)
      if(gridType == elemTypesGrid[j]) break;
    volElem[i].standardElemGrid = standardVolumeElementsGrid[j];

    /*--- Set the pointer for the standard element of the solution. ---*/
    for(j=0; j<elemTypesSol.size(); ++j)
      if(solType == elemTypesSol[j]) break;
    volElem[i].standardElemFlow = standardVolumeElementsSolution[j];

    /*--- Set the pointer for the standard element of the pressure, if needed. ---*/
    if( incompressible ) {
      for(j=0; j<elemTypesSol.size(); ++j)
        if(pType == elemTypesSol[j]) break;
      volElem[i].standardElemP = standardVolumeElementsSolution[j];
    }
  }
}

unsigned short CMeshFEM_DG::DetermineMaxNDOFs(void) {

  /*--- Loop over the standard volume elements for the solution and determine
        the maximum number of DOFs that are present. ---*/
  unsigned short nDOFsMax = 0;
  for(unsigned short l=0; l<standardVolumeElementsSolution.size(); ++l)
    nDOFsMax = max(nDOFsMax, standardVolumeElementsSolution[l]->GetNDOFsPad());

  /*--- Return the value of nDOFsMax. ---*/
  return nDOFsMax;
}

unsigned short CMeshFEM_DG::DetermineMaxNIntegration(void) {

  /*--- Loop over the standard volume elements for the solution and determine
        the maximum number of integration points that are present. ---*/
  unsigned short nIntMax = 0;
  for(unsigned short l=0; l<standardVolumeElementsSolution.size(); ++l)
    nIntMax = max(nIntMax, standardVolumeElementsSolution[l]->GetNIntegrationPad());

  /*--- Take the surface elements into account as well, although the number of
        integration points for the surfaces should be less than for the volume. ---*/
  for(unsigned short l=0; l<standardSurfaceElementsSolution.size(); ++l)
    nIntMax = max(nIntMax, standardSurfaceElementsSolution[l]->GetNIntegrationPad());

  /*--- Return the value of nIntMax. ---*/
  return nIntMax;
}

long CMeshFEM_DG::GetGlobal_to_Local_Point(unsigned long val_ipoint) const {

  /*--- Try to locate the give ID in the map Global_to_Local_Point.
        Return the local index if successfull and -1 otherwise. ---*/
  auto it = Global_to_Local_Point.find(val_ipoint);
  if(it != Global_to_Local_Point.cend()) return it->second;
  return -1;
}

void CMeshFEM_DG::InitStaticMeshMovement(const CConfig        *config,
                                         const unsigned short Kind_Grid_Movement,
                                         const unsigned short iZone) {

  /*---- Start of the OpenMP parallel region, if supported. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Get the reference values for the non-dimensionalization of the
          prescribed velocities. ---*/
    const su2double L_RefInv  = 1.0/config->GetLength_Ref();
    const su2double Omega_Ref = config->GetOmega_Ref();
    const su2double Vel_Ref   = config->GetVelocity_Ref();

    /*--- Make a distinction between the possibilities. ---*/
    switch( Kind_Grid_Movement ) {

      /*------------------------------------------------------------*/

      case ROTATING_FRAME: {

        /*--- Get the rotation rate and rotation center from config. ---*/
        const su2double Center[] = {config->GetMotion_Origin(0),
                                    config->GetMotion_Origin(1),
                                    config->GetMotion_Origin(2)};
        const su2double Omega[]  = {config->GetRotation_Rate(0)/Omega_Ref,
                                    config->GetRotation_Rate(1)/Omega_Ref,
                                    config->GetRotation_Rate(2)/Omega_Ref};

        /*--- Loop over the owned number of elements. ---*/
#ifdef HAVE_OMP
        const size_t omp_chunk_size = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
        SU2_OMP(for schedule(static,omp_chunk_size) SU2_NOWAIT)
        for(unsigned long l=0; l<nVolElemTot; ++l) {

          /*--- Determine the grid velocities in the integration points. ---*/
          const unsigned short nIntPad = volElem[l].standardElemGrid->GetNIntegrationPad();
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nIntPad; ++i) {

            su2double dist[3] = {0.0};
            dist[0] = (volElem[l].coorIntegrationPoints(i,0)-Center[0])*L_RefInv;
            dist[1] = (volElem[l].coorIntegrationPoints(i,1)-Center[1])*L_RefInv;
            if(nDim == 3) dist[2] = (volElem[l].coorIntegrationPoints(i,2)-Center[2])*L_RefInv;

            volElem[l].gridVelocitiesInt(i,0) = Omega[1]*dist[2] - Omega[2]*dist[1];
            volElem[l].gridVelocitiesInt(i,1) = Omega[2]*dist[0] - Omega[0]*dist[2];
            if(nDim == 3) volElem[l].gridVelocitiesInt(i,2) = Omega[0]*dist[1] - Omega[1]*dist[0];
          }

          /*--- Determine the grid velocities in the solution DOFs. ---*/
          const unsigned short nDOFsPad = volElem[l].standardElemGrid->GetNSolDOFsPad();
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i) {

            su2double dist[3] = {0.0};
            dist[0] = (volElem[l].coorSolDOFs(i,0)-Center[0])*L_RefInv;
            dist[1] = (volElem[l].coorSolDOFs(i,1)-Center[1])*L_RefInv;
            if(nDim == 3) dist[2] = (volElem[l].coorSolDOFs(i,2)-Center[2])*L_RefInv;

            volElem[l].gridVelocitiesSolDOFs(i,0) = Omega[1]*dist[2] - Omega[2]*dist[1];
            volElem[l].gridVelocitiesSolDOFs(i,1) = Omega[2]*dist[0] - Omega[0]*dist[2];
            if(nDim == 3) volElem[l].gridVelocitiesSolDOFs(i,2) = Omega[0]*dist[1] - Omega[1]*dist[0];
          }
        }

        /*--- Loop over the internally matching faces. ---*/
#ifdef HAVE_OMP
        const size_t omp_chunk_size_face = computeStaticChunkSize(matchingFaces.size(), omp_get_num_threads(), 64);
#endif
        SU2_OMP(for schedule(static,omp_chunk_size_face) SU2_NOWAIT)
        for(unsigned long l=0; l<matchingFaces.size(); ++l) {

          /*--- Determine the grid velocities in the integration points. ---*/
          const unsigned short nIntPad = matchingFaces[l].standardElemGrid->GetNIntegrationPad();
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nIntPad; ++i) {

            su2double dist[3] = {0.0};
            dist[0] = (matchingFaces[l].coorIntegrationPoints(i,0)-Center[0])*L_RefInv;
            dist[1] = (matchingFaces[l].coorIntegrationPoints(i,1)-Center[1])*L_RefInv;
            if(nDim == 3) dist[2] = (matchingFaces[l].coorIntegrationPoints(i,2)-Center[2])*L_RefInv;

            matchingFaces[l].gridVelocities(i,0) = Omega[1]*dist[2] - Omega[2]*dist[1];
            matchingFaces[l].gridVelocities(i,1) = Omega[2]*dist[0] - Omega[0]*dist[2];
            if(nDim == 3) matchingFaces[l].gridVelocities(i,2) = Omega[0]*dist[1] - Omega[1]*dist[0];
          }
        }

        /*--- Loop over the physical boundaries. Ignore the periodic boundaries. ---*/
        for(unsigned short iMarker=0; iMarker<boundaries.size(); ++iMarker) {
          if( !boundaries[iMarker].periodicBoundary ) {

            /*--- Easier storage of the surface elements and determine the chunk size
                  for the OpenMP parallelization, if supported. ---*/
            vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
#ifdef HAVE_OMP
            const size_t omp_chunk_size_surf = computeStaticChunkSize(surfElem.size(), omp_get_num_threads(), 64);
#endif
            /*--- Check if this is a shroud boundary. ---*/
            bool shroudBoundary = false;
            for(unsigned short iShroud=0; iShroud<config->GetnMarker_Shroud(); ++iShroud) {
              if(boundaries[iMarker].markerTag == config->GetMarker_Shroud(iShroud)) {
                shroudBoundary = true;
                break;
              }
            }

            /*--- For a shroud boundary the grid velocities must be set to zero. ---*/
            if( shroudBoundary ) {
              SU2_OMP(for schedule(static,omp_chunk_size_surf) SU2_NOWAIT)
              for(unsigned long l=0; l<surfElem.size(); ++l)
                surfElem[l].coorIntegrationPoints.setConstant(0.0);
            }
            else {

              /*--- Normal boundary. Loop over the number of boundary faces. ---*/
              SU2_OMP(for schedule(static,omp_chunk_size_surf) SU2_NOWAIT)
              for(unsigned long l=0; l<surfElem.size(); ++l) {

                /*--- Determine the grid velocities in the integration points. ---*/
                const unsigned short nIntPad = surfElem[l].standardElemGrid->GetNIntegrationPad();
                SU2_OMP_SIMD_IF_NOT_AD
                for(unsigned short i=0; i<nIntPad; ++i) {

                  su2double dist[3] = {0.0};
                  dist[0] = (surfElem[l].coorIntegrationPoints(i,0)-Center[0])*L_RefInv;
                  dist[1] = (surfElem[l].coorIntegrationPoints(i,1)-Center[1])*L_RefInv;
                  if(nDim == 3) dist[2] = (surfElem[l].coorIntegrationPoints(i,2)-Center[2])*L_RefInv;

                  surfElem[l].gridVelocities(i,0) = Omega[1]*dist[2] - Omega[2]*dist[1];
                  surfElem[l].gridVelocities(i,1) = Omega[2]*dist[0] - Omega[0]*dist[2];
                  if(nDim == 3) surfElem[l].gridVelocities(i,2) = Omega[0]*dist[1] - Omega[1]*dist[0];
                }
              }
            }
          }
        }

        break;
      }

      /*------------------------------------------------------------*/

      case STEADY_TRANSLATION: {

        /*--- Get the translation velocity from config. ---*/
        const su2double vTrans[] = {config->GetTranslation_Rate(0)/Vel_Ref,
                                    config->GetTranslation_Rate(1)/Vel_Ref,
                                    config->GetTranslation_Rate(2)/Vel_Ref};

      /*--- Loop over the owned number of elements. ---*/
#ifdef HAVE_OMP
        const size_t omp_chunk_size = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
        SU2_OMP(for schedule(static,omp_chunk_size) SU2_NOWAIT)
        for(unsigned long l=0; l<nVolElemTot; ++l) {

          /*--- Determine the grid velocities in the integration points. ---*/
          const unsigned short nIntPad = volElem[l].standardElemGrid->GetNIntegrationPad();
          for(unsigned short j=0; j<nDim; ++j ) {
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nIntPad; ++i)
              volElem[l].gridVelocitiesInt(i,j) = vTrans[j];
          }

          /*--- Determine the grid velocities in the solution DOFs. ---*/
          const unsigned short nDOFsPad = volElem[l].standardElemGrid->GetNSolDOFsPad();
          for(unsigned short j=0; j<nDim; ++j ) {
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              volElem[l].gridVelocitiesSolDOFs(i,j) = vTrans[j];
          }
        }

        /*--- Loop over the internally matching faces. ---*/
#ifdef HAVE_OMP
        const size_t omp_chunk_size_face = computeStaticChunkSize(matchingFaces.size(), omp_get_num_threads(), 64);
#endif
        SU2_OMP(for schedule(static,omp_chunk_size_face) SU2_NOWAIT)
        for(unsigned long l=0; l<matchingFaces.size(); ++l) {

          /*--- Determine the grid velocities in the integration points. ---*/
          const unsigned short nIntPad = matchingFaces[l].standardElemGrid->GetNIntegrationPad();
          for(unsigned short j=0; j<nDim; ++j ) {
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nIntPad; ++i)
              matchingFaces[l].gridVelocities(i,j) = vTrans[j];
          }
        }

        /*--- Loop over the physical boundaries. Ignore the periodic boundaries. ---*/
        for(unsigned short iMarker=0; iMarker<boundaries.size(); ++iMarker) {
          if( !boundaries[iMarker].periodicBoundary ) {

            /*--- Easier storage of the surface elements and loop over them. ---*/
            vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
#ifdef HAVE_OMP
            const size_t omp_chunk_size_surf = computeStaticChunkSize(surfElem.size(), omp_get_num_threads(), 64);
#endif
            SU2_OMP(for schedule(static,omp_chunk_size_surf) SU2_NOWAIT)
            for(unsigned long l=0; l<surfElem.size(); ++l) {

              /*--- Determine the grid velocities in the integration points. ---*/
              const unsigned short nIntPad = surfElem[l].standardElemGrid->GetNIntegrationPad();
              for(unsigned short j=0; j<nDim; ++j ) {
                SU2_OMP_SIMD_IF_NOT_AD
                for(unsigned short i=0; i<nIntPad; ++i)
                  surfElem[l].gridVelocities(i,j) = vTrans[j];
              }
            }
          }
        }

        break;
      }

      /*------------------------------------------------------------*/

      default:  /* Just to avoid a compiler warning. */
        break;
    }

    /*--- Check if moving walls are present. ---*/
    if( config->GetSurface_Movement(MOVING_WALL) ) {

      /*--- Loop over the physical boundaries. Skip the periodic boundaries. ---*/
      for(unsigned short iMarker=0; iMarker<boundaries.size(); ++iMarker) {
        if( !boundaries[iMarker].periodicBoundary ) {

          /*--- Check if for this boundary a motion has been specified. ---*/
          if(config->GetMarker_All_Moving(iMarker) == YES) {

            /*--- Determine the prescribed translation velocity, rotation rate
                  and rotation center. ---*/
            const su2double Center[] = {config->GetMotion_Origin(0),
                                        config->GetMotion_Origin(1),
                                        config->GetMotion_Origin(2)};
            const su2double Omega[]  = {config->GetRotation_Rate(0)/Omega_Ref,
                                        config->GetRotation_Rate(1)/Omega_Ref,
                                        config->GetRotation_Rate(2)/Omega_Ref};
            const su2double vTrans[] = {config->GetTranslation_Rate(0)/Vel_Ref,
                                        config->GetTranslation_Rate(1)/Vel_Ref,
                                        config->GetTranslation_Rate(2)/Vel_Ref};

            /*--- Easier storage of the surface elements and loop over them. ---*/
            vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;

#ifdef HAVE_OMP
            const size_t omp_chunk_size_surf = computeStaticChunkSize(surfElem.size(), omp_get_num_threads(), 64);
#endif
            SU2_OMP(for schedule(static,omp_chunk_size_surf) SU2_NOWAIT)
            for(unsigned long l=0; l<surfElem.size(); ++l) {

              /*--- Loop over the number of integration points. ---*/
              const unsigned short nIntPad = surfElem[l].standardElemGrid->GetNIntegrationPad();
              SU2_OMP_SIMD_IF_NOT_AD
              for(unsigned short i=0; i<nIntPad; ++i) {

                /*--- Calculate the non-dimensional distance from the
                      rotation center. ---*/
                su2double r[3] = {0.0};
                r[0] = (volElem[l].coorSolDOFs(i,0)-Center[0])*L_RefInv;
                r[1] = (volElem[l].coorSolDOFs(i,1)-Center[1])*L_RefInv;
                if(nDim == 3) r[2] = (volElem[l].coorSolDOFs(i,2)-Center[2])*L_RefInv;

                /*--- Cross Product of angular velocity and distance from center to
                      get the rotational velocity. Note that we are adding on the
                      velocity due to pure translation as well. Also note that for
                      that 2D case only Omega[2] can be non-zero. ---*/
                su2double velGrid[] = {vTrans[0] + Omega[1]*r[2] - Omega[2]*r[1],
                                       vTrans[1] + Omega[2]*r[0] - Omega[0]*r[2],
                                       vTrans[2] + Omega[0]*r[1] - Omega[1]*r[0]};

                /*--- Store the grid velocities. ---*/
                surfElem[l].gridVelocities(i,0) = velGrid[0];
                surfElem[l].gridVelocities(i,1) = velGrid[1];
                if(nDim == 3) surfElem[l].gridVelocities(i,2) = velGrid[2];
              }
            }
          }
        }
      }
    }

  }
  END_SU2_OMP_PARALLEL
}

void CMeshFEM_DG::MetricTermsSurfaceElements(CConfig *config) {

  /*--- The master node writes a message. ---*/
  if(rank == MASTER_NODE) cout << "Computing metric terms surface elements." << endl;

  /*---- Start of the OpenMP parallel region, if supported. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Loop over the internal matching faces. ---*/
#ifdef HAVE_OMP
    const size_t omp_chunk_size_face = computeStaticChunkSize(matchingFaces.size(), omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT(omp_chunk_size_face)
    for(unsigned long i=0; i<matchingFaces.size(); ++i) {

      /*--- Initialize the grid velocities. ---*/
      matchingFaces[i].InitGridVelocities(nDim);

      /*--- Compute the metric terms in the integration points. ---*/
      matchingFaces[i].MetricTermsIntegrationPoints(nDim, volElem);
    }
    END_SU2_OMP_FOR

    /*--- Loop over the physical boundaries. ---*/
    for(unsigned int l=0; l<boundaries.size(); ++l) {
      if( !boundaries[l].periodicBoundary ) {

        /*--- Loop over the boundary faces. ---*/
#ifdef HAVE_OMP
        const size_t omp_chunk_size_surf = computeStaticChunkSize(boundaries[l].surfElem.size(),
                                                                  omp_get_num_threads(), 64);
#endif
        SU2_OMP_FOR_STAT(omp_chunk_size_surf)
        for(unsigned long i=0; i<boundaries[l].surfElem.size(); ++i) {

          /*--- Initialize the grid velocities. ---*/
          boundaries[l].surfElem[i].InitGridVelocities(nDim);

          /*--- Compute the metric terms in the integration points. ---*/
          boundaries[l].surfElem[i].MetricTermsIntegrationPoints(nDim, volElem);
        }        
        END_SU2_OMP_FOR
      }
    }
  }
  END_SU2_OMP_PARALLEL
}

void CMeshFEM_DG::MetricTermsVolumeElements(CConfig *config) {

  /*--- The master node writes a message. ---*/
  if(rank == MASTER_NODE) cout << "Computing metric terms volume elements." << endl;

  /*--- Determine whether or not the LGL node distribution is used. ---*/
  const bool useLGL = config->GetKind_FEM_GridDOFsLocation() == LGL;

  /*--- Initialize the boolean DerMetricTerms to false to indicate that
        the derivatives of the metric terms are not needed. ---*/
  bool DerMetricTerms = false;

  /*--- Check if ADER-DG is used as explicit time integration scheme for
        the unsteady simulation. ---*/
  if((config->GetTime_Marching()           == TIME_MARCHING::TIME_STEPPING) &&
     (config->GetKind_TimeIntScheme_Flow() == ADER_DG)) {

    /*--- Check for a compressible Navier-Stokes simulation. ---*/
    MAIN_SOLVER solver = config->GetKind_Solver();
    if(solver == MAIN_SOLVER::FEM_NAVIER_STOKES || solver == MAIN_SOLVER::FEM_RANS || solver == MAIN_SOLVER::FEM_LES) {

      /*--- The derivatives of the metric terms are needed when a non-aliased
            predictor is used for ADER. ---*/
      if(config->GetKind_ADER_Predictor() == ADER_NON_ALIASED_PREDICTOR)
        DerMetricTerms = true;

      /*--- Determine the time coefficients in the iteration matrix of the ADER-DG
            predictor step. ---*/
      TimeCoefficientsPredictorADER_DG(config);
    }
  }

  /*--- Determine the chunk size for the OMP loop below. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif

  /*--- Initialize the number of elements with negative Jacobians. ---*/
  unsigned long nElemNegJac = 0;

  /*---- Start of the OpenMP parallel region, if supported. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Definition of the number of negative Jacobians for this thread.
          The reduction is handled manually to avoid complications
          with CODIPACK. ---*/
    unsigned long nNegJac = 0;

    /*--- Loop over all elements to set the coordinates. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for(unsigned long i=0; i<nVolElemTot; ++i)
      volElem[i].SetCoorGridDOFs(nDim, meshPoints);
    END_SU2_OMP_FOR

    /*--- Loop over the owned volume elements. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for(unsigned long i=0; i<nVolElemOwned; ++i) {

      /*--- Initialize the grid velocities. ---*/
      volElem[i].InitGridVelocities(nDim);

      /*--- Compute the metric terms in the integration points and
            in the solution DOFs. Update the number of bad elements
            if negative Jacobians are present. ---*/
      if( !(volElem[i].MetricTermsIntegrationPoints(useLGL, nDim)) ||
          !(volElem[i].MetricTermsSolDOFs(nDim)) ) ++nNegJac;

      /*--- Compute the derivatives of the metric terms in the
            integration points, if needed. ---*/
      if( DerMetricTerms )
        volElem[i].DerMetricTermsIntegrationPoints(nDim);
    }
    END_SU2_OMP_FOR

    /*--- Carry out the reduction over the threads. ---*/
    SU2_OMP_CRITICAL
    {
      nElemNegJac += nElemNegJac;
    }
    END_SU2_OMP_CRITICAL
  }
  END_SU2_OMP_PARALLEL

  /*--- Determine the global number of elements with negative Jacobians. ---*/
#ifdef HAVE_MPI
  unsigned long nElemNegJacLoc = nElemNegJac;

  SU2_MPI::Allreduce(&nElemNegJacLoc, &nElemNegJac, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- Terminate if there are elements with negative Jacobians. ---*/
  if(nElemNegJac > 0) {
    ostringstream message;
    message << "Found " <<  nElemNegJac << " elements with negative Jacobians.";
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }
}

void CMeshFEM_DG::SetGlobal_to_Local_Point(void) {

  /*--- The master node writes a message. ---*/
  if(rank == MASTER_NODE) cout << "Storing a mapping from global to local DOF index." << endl;

  /*--- Clear the current map and set the data. ---*/
  Global_to_Local_Point.clear();
  unsigned long ii = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j, ++ii) {
      Global_to_Local_Point[volElem[i].offsetDOFsSolGlobal+j] = ii;
    }
  }
}

void CMeshFEM_DG::SetSendReceive(const CConfig *config) {

  /*----------------------------------------------------------------------------*/
  /*--- Step 1: Determine the ranks from which this rank has to receive data ---*/
  /*---         during the actual communication of halo data, as well as the ---*/
  /*            data that must be communicated.                              ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Determine the ranks from which this rank will receive halo data. ---*/
  vector<int> recvFromRank(size, 0);
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i)
    recvFromRank[volElem[i].rankOriginal] = 1;

  map<int,int> rankToIndRecvBuf;
  for(int i=0; i<size; ++i) {
    if( recvFromRank[i] ) {
      int ind = (int)rankToIndRecvBuf.size();
      rankToIndRecvBuf[i] = ind;
    }
  }

  ranksRecv.resize(rankToIndRecvBuf.size());
  map<int,int>::const_iterator MI = rankToIndRecvBuf.begin();
  for(unsigned long i=0; i<rankToIndRecvBuf.size(); ++i, ++MI)
    ranksRecv[i] = MI->first;

  /*--- Define and determine the buffers to send the global indices of my halo
        elements to the appropriate ranks and the vectors which store the
        elements that I will receive from these ranks. ---*/
  vector<vector<unsigned long> > longBuf(rankToIndRecvBuf.size(), vector<unsigned long>(0));
  entitiesRecv.resize(rankToIndRecvBuf.size());

  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i) {
    MI = rankToIndRecvBuf.find(volElem[i].rankOriginal);
    longBuf[MI->second].push_back(volElem[i].elemIDGlobal);

    entitiesRecv[MI->second].push_back(i);
  }

  /*--- Determine the mapping from global element ID to local owned element ID. ---*/
  map<unsigned long,unsigned long> globalElemIDToLocalInd;
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    globalElemIDToLocalInd[volElem[i].elemIDGlobal] = i;

#ifdef HAVE_MPI

  /*--- Parallel mode. First determine the number of ranks to which this
        rank has to send halo data during the actual exchange. ---*/
  int nRankSend;
  vector<int> sizeReduce(size, 1);

  SU2_MPI::Reduce_scatter(recvFromRank.data(), &nRankSend, sizeReduce.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /*--- Resize ranksSend and the first index of entitiesSend to the number of
        ranks to which this rank has to send data. ---*/
  ranksSend.resize(nRankSend);
  entitiesSend.resize(nRankSend);

  /*--- Send all the data using non-blocking sends. ---*/
  vector<SU2_MPI::Request> commReqs(ranksRecv.size());

  for(unsigned long i=0; i<ranksRecv.size(); ++i) {
    int dest = ranksRecv[i];
    SU2_MPI::Isend(longBuf[i].data(), longBuf[i].size(), MPI_UNSIGNED_LONG,
                   dest, dest, SU2_MPI::GetComm(), &commReqs[i]);
  }

  /*--- Loop over the number of ranks from which I receive data about the
        global element ID's that I must send. ---*/
  for(int i=0; i<nRankSend; ++i) {

    /*--- Block until a message arrives and determine the source. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    ranksSend[i] = status.MPI_SOURCE;

    /*--- Determine the size of the message, allocate the memory for the
          receive buffer and receive the message. ---*/
    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    entitiesSend[i].resize(sizeMess);
    SU2_MPI::Recv(entitiesSend[i].data(), sizeMess, MPI_UNSIGNED_LONG,
                  ranksSend[i], rank, SU2_MPI::GetComm(), &status);

    /*--- Convert the global indices currently stored in entitiesSend[i]
          to local indices. ---*/
    for(int j=0; j<sizeMess; ++j) {
      map<unsigned long,unsigned long>::const_iterator LMI;
      LMI = globalElemIDToLocalInd.find(entitiesSend[i][j]);

      if(LMI == globalElemIDToLocalInd.end())
        SU2_MPI::Error("This should not happen", CURRENT_FUNCTION);

      entitiesSend[i][j] = LMI->second;
    }
  }

  /*--- Complete the non-blocking sends and synchronize the ranks, because
        wild cards have been used. ---*/
  SU2_MPI::Waitall(ranksRecv.size(), commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else
  /*--- Sequential mode. Resize ranksSend and the first index of entitiesSend to
        the number of ranks to which this rank has to send data. This number is
        only non-zero when periodic boundaries are present in the grid. ---*/
  ranksSend.resize(ranksRecv.size());
  entitiesSend.resize(ranksRecv.size());

  /*--- Convert the global element ID's of longBuf to local indices, which are
        stored in entitiesSend[0]. Note that an additional test for longBuf.size()
        is necessary to avoid problems. ---*/
  if( longBuf.size() ) {

    ranksSend[0] = MASTER_NODE;
    entitiesSend[0].resize(longBuf[0].size());

    for(unsigned long i=0; i<longBuf[0].size(); ++i) {
      map<unsigned long,unsigned long>::const_iterator LMI;
      LMI = globalElemIDToLocalInd.find(longBuf[0][i]);

      if(LMI == globalElemIDToLocalInd.end())
        SU2_MPI::Error("This should not happen", CURRENT_FUNCTION);

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
  map<short,unsigned short> mapRotationalPeriodicToInd;

  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {

      auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
      if(fabs(angles[0]) > 1.e-5 || fabs(angles[1]) > 1.e-5 || fabs(angles[2]) > 1.e-5) {

        unsigned short curSize = mapRotationalPeriodicToInd.size();
        mapRotationalPeriodicToInd[iMarker] = curSize;
      }
    }
  }

  /*--- Store the rotationally periodic indices in rotPerMarkers. ---*/
  rotPerMarkers.reserve(mapRotationalPeriodicToInd.size());
  for(map<short,unsigned short>::iterator SMI =mapRotationalPeriodicToInd.begin();
                                          SMI!=mapRotationalPeriodicToInd.end(); ++SMI)
    rotPerMarkers.push_back(SMI->first);

  /*--- Resize the first index of rotPerHalos to the correct size. ---*/
  rotPerHalos.resize(mapRotationalPeriodicToInd.size());

  /*--- Loop over the halo volume elements and store the indices of the
        rotationally periodic halo elements in rotPerHalos.     ---*/
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i) {
    if(volElem[i].periodIndexToDonor > -1) {
      map<short,unsigned short>::const_iterator SMI;
      SMI = mapRotationalPeriodicToInd.find(volElem[i].periodIndexToDonor);

      if(SMI != mapRotationalPeriodicToInd.end())
        rotPerHalos[SMI->second].push_back(i);
    }
  }
}

void CMeshFEM_DG::WallFunctionPreprocessing(CConfig *config) {

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Check whether wall functions are used at all.              ---*/
  /*--------------------------------------------------------------------------*/

  bool wallFunctions = false;
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX: {
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if(config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE)
          wallFunctions = true;
        break;
      }
      default:  /* Just to avoid a compiler warning. */
        break;
    }
  }

  /* If no wall functions are used, nothing needs to be done and a
     return can be made. */
  if( !wallFunctions ) return;

  /*--- The master node writes a message. ---*/
  if(rank == MASTER_NODE) cout << "Preprocessing for the wall functions. " << endl;

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Build the local ADT of the volume elements. The halo       ---*/
  /*---         elements are included, because these may also be donors.   ---*/
  /*---         Note that the ADT is built with the linear subelements.    ---*/
  /*---         This is done to avoid relatively many expensive Newton     ---*/
  /*---         solves for high order elements.                            ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Define the vectors, which store the mapping from the subelement to the
        parent element, subelement ID within the parent element, the element
        type and the connectivity of the subelements. ---*/
  vector<unsigned long>  parentElement;
  vector<unsigned short> subElementIDInParent;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned long>  elemConn;

  /*--- Loop over the locally stored volume elements (including halo elements)
        to create the connectivity of the subelements. ---*/
  for(unsigned long l=0; l<nVolElemTot; ++l) {

    /*--- Determine the necessary data from the corresponding standard element. ---*/
    unsigned short VTK_SubType[] = {volElem[l].standardElemGrid->GetVTK_SubType1(),
                                    volElem[l].standardElemGrid->GetVTK_SubType2()};
    unsigned short nSubElems[]       = {0, 0};
    unsigned short nDOFsPerSubElem[] = {0, 0};

    const unsigned short *connSubElems[] = {nullptr, nullptr};

    if(VTK_SubType[0] != NONE) {
      nSubElems[0]       = volElem[l].standardElemGrid->GetNSubElemsType1();
      nDOFsPerSubElem[0] = volElem[l].standardElemGrid->GetNDOFsPerSubElem(VTK_SubType[0]);
      connSubElems[0]    = volElem[l].standardElemGrid->GetSubConnType1();
    }

    if(VTK_SubType[1] != NONE) {
      nSubElems[1]       = volElem[l].standardElemGrid->GetNSubElemsType2();
      nDOFsPerSubElem[1] = volElem[l].standardElemGrid->GetNDOFsPerSubElem(VTK_SubType[1]);
      connSubElems[1]    = volElem[l].standardElemGrid->GetSubConnType2();
    }

    /*--- Abbreviate the grid DOFs of this element a bit easier. ---*/
    const vector<unsigned long> &nodeIDs = volElem[l].nodeIDsGrid;

    /*--- Loop over the number of subelements and store the required data. ---*/
    unsigned short jj = 0;
    for(unsigned short i=0; i<2; ++i) {
      unsigned short kk = 0;
      for(unsigned short j=0; j<nSubElems[i]; ++j, ++jj) {
        parentElement.push_back(l);
        subElementIDInParent.push_back(jj);
        VTK_TypeElem.push_back(VTK_SubType[i]);

        for(unsigned short k=0; k<nDOFsPerSubElem[i]; ++k, ++kk)
          elemConn.push_back(nodeIDs[connSubElems[i][kk]]);
      }
    }
  }

  /*--- Copy the coordinates in a vector that can be used by the ADT. ---*/
  vector<su2double> volCoor;
  volCoor.reserve(nDim*meshPoints.size());

  for(unsigned long l=0; l<meshPoints.size(); ++l) {
    for(unsigned short k=0; k<nDim; ++k)
      volCoor.push_back(meshPoints[l].coor[k]);
  }

  /*--- Build the local ADT. ---*/
  CADTElemClass localVolumeADT(nDim, volCoor, elemConn, VTK_TypeElem,
                               subElementIDInParent, parentElement, false);

  /*--- Release the memory of the vectors used to build the ADT. To make sure
        that all the memory is deleted, the swap function is used. ---*/
  vector<unsigned short>().swap(subElementIDInParent);
  vector<unsigned short>().swap(VTK_TypeElem);
  vector<unsigned long>().swap(parentElement);
  vector<unsigned long>().swap(elemConn);
  vector<su2double>().swap(volCoor);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3. Search for donor elements at the exchange locations in     ---*/
  /*---         the local elements.                                        ---*/
  /*--------------------------------------------------------------------------*/

  /*---- Start of the OpenMP parallel region, if supported. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Loop over the markers and select the ones for which a wall function
          treatment must be carried out. ---*/
    for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

      switch (config->GetMarker_All_KindBC(iMarker)) {
        case ISOTHERMAL:
        case HEAT_FLUX: {
          const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if(config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE){

            /*--- An LES wall model is used for this boundary marker. Determine
                  which wall model and allocate the memory for the member variable. ---*/
            SU2_OMP_SINGLE
            {
              switch (config->GetWallFunction_Treatment(Marker_Tag) ) {
                case WALL_FUNCTIONS::EQUILIBRIUM_MODEL: {
                  if(rank == MASTER_NODE)
                    cout << "Marker " << Marker_Tag << " uses an Equilibrium Wall Model." << endl;

                  boundaries[iMarker].wallModel = new CWallModel1DEQ(config, Marker_Tag);
                  break;
                }
                case WALL_FUNCTIONS::LOGARITHMIC_MODEL: {
                  if(rank == MASTER_NODE)
                    cout << "Marker " << Marker_Tag << " uses the Reichardt and Kader analytical laws for the Wall Model." << endl;

                  boundaries[iMarker].wallModel = new CWallModelLogLaw(config, Marker_Tag);
                  break;
                }
                default: {
                  SU2_MPI::Error("Wall function not present yet", CURRENT_FUNCTION);
                }
              }
            }
            END_SU2_OMP_SINGLE

            /*--- Retrieve the double information for this wall model. The height
                  of the exchange location is the first element of this array. ---*/
            const su2double *doubleInfo = config->GetWallFunction_DoubleInfo(Marker_Tag);

            /*--- Easier storage of the surface elements. ---*/
            vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;

            /*--- Determine the chunk size for the OMP loop below. ---*/
#ifdef HAVE_OMP
            const size_t omp_chunk_size = computeStaticChunkSize(surfElem.size(),
                                                                 omp_get_num_threads(), 64);
#endif
            /*--- Loop over the local boundary elements for this marker. ---*/
            SU2_OMP_FOR_DYN(omp_chunk_size)
            for(unsigned long l=0; l<surfElem.size(); ++l) {

              /*--- Easier storage of the number of integration points. ---*/
              const unsigned short nInt = surfElem[l].standardElemGrid->GetNIntegration();

              /*--- Allocate the memory for the memory to store the donors and
                    the parametric weights. The donor elements are stored in an
                    CUnsignedLong2T, such that they can be sorted. ---*/
              vector<CUnsignedLong2T> donorElements(nInt);
              ColMajorMatrix<su2double> parCoorInDonor(nInt,nDim);

              /*--- Loop over the integration points of the face. ---*/
              for(unsigned short i=0; i<nInt; ++i) {

                /*--- Determine the coordinates of the exchange point. ---*/
                su2double coorExchange[] = {0.0, 0.0, 0.0};  // To avoid a compiler warning.
                for(unsigned short iDim=0; iDim<nDim; ++iDim)
                  coorExchange[iDim] = surfElem[l].coorIntegrationPoints(i,iDim)
                                     - surfElem[l].metricNormalsFace(i,iDim)*doubleInfo[0];

                /*--- Search for the element, which contains the exchange location. ---*/
                unsigned short subElem;
                unsigned long  parElem;
                int            rank;
                su2double      parCoor[3], weightsInterpol[8];
                if( localVolumeADT.DetermineContainingElement(coorExchange, subElem,
                                                              parElem, rank, parCoor,
                                                              weightsInterpol) ) {

                  /*--- Subelement found that contains the exchange location. However,
                        what is needed is the location in the high order parent element.
                        Determine this. ---*/
                  HighOrderContainmentSearch(coorExchange, parElem, subElem,
                                             weightsInterpol, parCoor);

                  /*--- Store the info in donorElements and parCoorInDonor. ---*/
                  donorElements[i].long0 = parElem;
                  donorElements[i].long1 = i;
                  for(unsigned short iDim=0; iDim<nDim; ++iDim)
                    parCoorInDonor(i,iDim) = parCoor[iDim];
                }
                else {

                  /* No subelement found that contains the exchange location.
                     The partitioning is done such that this should not happen.
                     Print an error message and exit. */
                  SU2_MPI::Error(string("Exchange location not found in ADT. This should not happen"),
                                 CURRENT_FUNCTION);
                }
              }

              /*--- Sort donorElements in increasing order, such that the
                    integration points with the same donor are grouped together. ---*/
              std::sort(donorElements.begin(), donorElements.end());

              /*--- Store the donor information in the member variables of surfElem. ---*/
              surfElem[l].donorsWallFunction.push_back(donorElements[0].long0);
              surfElem[l].nIntPerWallFunctionDonor.push_back(0);
              surfElem[l].intPerWallFunctionDonor.resize(nInt);
              surfElem[l].intPerWallFunctionDonor[0] = donorElements[0].long1;

              for(unsigned short i=1; i<nInt; ++i) {
                if(donorElements[i].long0 != surfElem[l].donorsWallFunction.back()) {
                  surfElem[l].donorsWallFunction.push_back(donorElements[i].long0);
                  surfElem[l].nIntPerWallFunctionDonor.push_back(i);
                }
                surfElem[l].intPerWallFunctionDonor[i] = donorElements[i].long1;
              }

              surfElem[l].nIntPerWallFunctionDonor.push_back(nInt);

              /*--- Determine whether or not halo information is needed to apply
                    the boundary conditions for this surface. If needed, set
                    haloInfoNeededForBC of the boundary to true. As donorsWallFunction is
                    sorted in increasing order, it suffices to check the last element. ---*/
              if(surfElem[l].donorsWallFunction.back() >= nVolElemOwned)
                boundaries[iMarker].haloInfoNeededForBC = true;

              /*--- Allocate the memory of the first index of the interpolation
                    matrices for the donordata. ---*/
              surfElem[l].matWallFunctionDonor.resize(surfElem[l].donorsWallFunction.size());

              /*--- Loop over the different donors for the wall function data. ---*/
              for(unsigned long j=0; j<surfElem[l].donorsWallFunction.size(); ++j) {

                /*--- Easier storage of the donor and the vector for the
                      interpolation data for this donor. ---*/
                const unsigned long donor               = surfElem[l].donorsWallFunction[j];
                ColMajorMatrix<passivedouble> &matDonor = surfElem[l].matWallFunctionDonor[j];

                /*--- Determine the number of integration points for this donor. ---*/
                const unsigned short nIntThisDonor = surfElem[l].nIntPerWallFunctionDonor[j+1]
                                                   - surfElem[l].nIntPerWallFunctionDonor[j];

                /*--- Allocate the memory to store the parametric coordinates, which
                      must be handled by the current donor. Note that a passive double
                      must be used here. Set the appropriate values afterwards. ---*/
                vector<vector<passivedouble> > parCoor(nDim, vector<passivedouble>(nIntThisDonor));
                for(unsigned short i=surfElem[l].nIntPerWallFunctionDonor[j];
                                   i<surfElem[l].nIntPerWallFunctionDonor[j+1]; ++i) {
                  const unsigned short ii = surfElem[l].intPerWallFunctionDonor[i];
                  for(unsigned short k=0; k<nDim; ++k)
                    parCoor[k][i] = SU2_TYPE::GetValue(parCoorInDonor(ii,k));
                }

                /*--- Determine the values of the basis functions in these points.
                      Note that the standard element of the solution must be used
                      for this purpose, because the solution must be interpolated
                      in these exchange points. ---*/
                volElem[donor].standardElemFlow->BasisFunctionsInPoints(parCoor, matDonor);
              }
            }
            END_SU2_OMP_FOR
          }

          break;
        }
        default:  /* Just to avoid a compiler warning. */
          break;
      }
    }
  }
  END_SU2_OMP_PARALLEL
}

/*---------------------------------------------------------------------*/
/*---          Private member functions of CMeshFEM_DG.             ---*/
/*---------------------------------------------------------------------*/

void CMeshFEM_DG::CreateStandardFaces(CConfig                      *config,
                                      const vector<CFaceOfElement> &localFaces) {

  /*--- Determine whether or not the LGL node distribution is used. ---*/
  const bool useLGL = config->GetKind_FEM_GridDOFsLocation() == LGL;

  /*--- Check if an incompressible solver is used. ---*/
  MAIN_SOLVER Kind_Solver = config->GetKind_Solver();
  const bool incompressible = (Kind_Solver == MAIN_SOLVER::FEM_INC_EULER) ||
                              (Kind_Solver == MAIN_SOLVER::FEM_INC_NAVIER_STOKES) ||
                              (Kind_Solver == MAIN_SOLVER::FEM_INC_RANS) ||
                              (Kind_Solver == MAIN_SOLVER::FEM_INC_LES);

  /*--- Define the number of variables for each type of solve,
        depending on the situation. If MKL is used the actual number
        of variables per point must be known, such that a jitted
        gemm call can be constructed. Otherwise this is not necessary
        and this variable can change during runtime. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  const unsigned short nGridVar = nDim;
  const unsigned short nSolVar  = incompressible ? nDim : nDim+2;
#else
  const unsigned short nGridVar = 1;
  const unsigned short nSolVar  = 1;
#endif

  /*--- Determine the number of variables for which memory must be
        allocated for the working variables. ---*/
  const unsigned short nAllocVar = incompressible ? nDim : nDim+2;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Determine the types of surface standard elements needed    ---*/
  /*---         for the handling of the boundary conditions. This data is  ---*/
  /*---         stored in a vector of the class CUnsignedShort8T, which    ---*/
  /*---         can store 8 short integers as one entity. Find below the   ---*/
  /*---         explanation of the 8 parameters that define a surface.     ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Every standard face for a boundary surface is determined by 7
        parameters, namely
        1: VTK type of the face. Although this is not an independent parameter,
           it is taken into account for convenience.
        2: Number of integration points of the face.
        3: VTK type of the volume element.
        4: Polynomial degree of the volume element.
        5: Face ID of the volume element. This determines parameter 1.
        6: Orientation of the element w.r.t. the face. For a boundary surface
           this is always 0, but not for internally matching faces, for which
           standard faces are also used.
        7: The number of variables per point.
        8: Polynomial degree that must be integrated exactly.
        Define the variables to store this information for both the grid
        and solution. ---*/
  vector<CUnsignedShort8T> surfaceTypesGrid, surfaceTypesSol;

  /*--- Loop over all faces and select the boundary faces. ---*/
  for(unsigned long i=0; i<localFaces.size(); ++i) {
    if(localFaces[i].faceIndicator > -1) {

      /*--- Determine the polynomial degree that must be integrated exactly
            by the integration rule of the face. ---*/
      const unsigned short nPoly      = max(localFaces[i].nPolyGrid0, localFaces[i].nPolySol0);
      const bool           constJac   = localFaces[i].JacFaceIsConsideredConstant;
      const unsigned short orderExact = config->GetOrderExactIntegrationFEM(nPoly, constJac);

      /*--- Determine the type of the face according to the VTK convention. ---*/
      unsigned short VTK_Type;
      if(     localFaces[i].nCornerPoints == 2) VTK_Type = LINE;
      else if(localFaces[i].nCornerPoints == 3) VTK_Type = TRIANGLE;
      else                                      VTK_Type = QUADRILATERAL;

      /*--- Determine the number of integration points of the face.
            Note that for a quadrilateral the number of integration points
            in one direction is determined. ---*/
      const unsigned short nInt = CFEMStandardElementBase::GetNIntStatic(VTK_Type, orderExact);

      /*--- Set the number of grid and solution variables. If the volume element is a
            hexahedron or a quadrilateral, these numbers do not matter and are therefore
            set to 1. ---*/
      unsigned short mGridVar = nGridVar, mSolVar = nSolVar;
      if((localFaces[i].elemType0 == HEXAHEDRON) || (localFaces[i].elemType0 == QUADRILATERAL))
        mGridVar = mSolVar = 1;

      /*--- Store this face in surfaceTypesGrid and surfaceTypesSol. For an
            incompressible solver two solution types must be added, one for
            the regular velocities and one for the pressure. ---*/
      surfaceTypesGrid.push_back(CUnsignedShort8T(VTK_Type, nInt, localFaces[i].elemType0,
                                                  localFaces[i].nPolyGrid0, localFaces[i].faceID0,
                                                  0, mGridVar, orderExact));

      surfaceTypesSol.push_back(CUnsignedShort8T(VTK_Type, nInt, localFaces[i].elemType0,
                                                 localFaces[i].nPolySol0, localFaces[i].faceID0,
                                                 0, mSolVar, orderExact));
      if( incompressible )
        surfaceTypesSol.push_back(CUnsignedShort8T(VTK_Type, nInt, localFaces[i].elemType0,
                                                   localFaces[i].nPolySol0-1, localFaces[i].faceID0,
                                                   0, 1, orderExact));
    }
  }

  /*--- Copy the face types for the grid, sort them and remove the
        double entities. ---*/
  vector<CUnsignedShort8T> typesSurfaceGrid = surfaceTypesGrid;
  vector<CUnsignedShort8T>::iterator lastEntry8T;
  sort(typesSurfaceGrid.begin(), typesSurfaceGrid.end());
  lastEntry8T = unique(typesSurfaceGrid.begin(), typesSurfaceGrid.end());
  typesSurfaceGrid.erase(lastEntry8T, typesSurfaceGrid.end());

  /*--- Copy the face types for the solution, sort them and remove the
        double entities. ---*/
  vector<CUnsignedShort8T> typesSurfaceSol = surfaceTypesSol;
  sort(typesSurfaceSol.begin(), typesSurfaceSol.end());
  lastEntry8T = unique(typesSurfaceSol.begin(), typesSurfaceSol.end());
  typesSurfaceSol.erase(lastEntry8T, typesSurfaceSol.end());

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Determine the types of the standard elements for internal  ---*/
  /*---         matching faces. This data is stored in a vector of the     ---*/
  /*---         the class CUnsignedShort11T, which can store 11 short      ---*/
  /*---         integers as one entity. Find below the explanation of the  ---*/
  /*---         11 parameters that define an internal matching face.       ---*/
  /*--------------------------------------------------------------------------*/

  /*--- For the internally matching faces every situation gets its own standard
        element to avoid copying and indirect addressing during the residual
        computation. This means that quite a few standard internal faces may be
        present.  The standard element for internal matching faces is defined by
        11 parameters, which are
        1:  VTK type of the face. Although this is not an independent parameter,
            it is taken into account for convenience.
        2:  Number of integration points of the face.
        3:  VTK type of the volume element on side 0.
        4:  Polynomial degree of the volume element on side 0.
        5:  Face ID of the volume element on side 0. This determines parameter 1.
        6:  VTK type of the volume element on side 1.
        7:  Polynomial degree of the volume element on side 1.
        8:  Face ID of the volume element on side 1.
        9:  Orientation of the element on side 1 w.r.t. the face.
            Note that the orientation of the element on side 0 is not needed,
            because the face has been constructed to match this orientation,
            see the call to localFaces[i].MatchOrientationElemSide0.
        10: The number of variables per point.
        11: Polynomial degree that must be integrated exactly.
        Define the variables to store this information for both the grid
        and solution. ---*/
  vector<CUnsignedShort11T> faceTypesGrid, faceTypesSol;

  /*--- Loop over all faces and select the matching internal faces. ---*/
  for(unsigned long i=0; i<localFaces.size(); ++i) {
    if(localFaces[i].faceIndicator == -1 && localFaces[i].elemID1 < nVolElemTot) {

      /*--- Determine the polynomial degree that must be integrated exactly
            by the integration rule of the face. ---*/
      const unsigned short nPoly      = max(max(localFaces[i].nPolyGrid0, localFaces[i].nPolySol0),
                                            max(localFaces[i].nPolyGrid1, localFaces[i].nPolySol1));
      const bool           constJac   = localFaces[i].JacFaceIsConsideredConstant;
      const unsigned short orderExact = config->GetOrderExactIntegrationFEM(nPoly, constJac);

      /*--- Determine the type of the face according to the VTK convention. ---*/
      unsigned short VTK_Type;
      if(     localFaces[i].nCornerPoints == 2) VTK_Type = LINE;
      else if(localFaces[i].nCornerPoints == 3) VTK_Type = TRIANGLE;
      else                                      VTK_Type = QUADRILATERAL;

      /*--- Set the number of grid and solution variables. If both adjacent volume elements
            are hexahedra or quadrilaterals, these number do not matter and are therefore
            set to 1. ---*/
      unsigned short mGridVar = nGridVar, mSolVar = nSolVar;
      if(((localFaces[i].elemType0 == HEXAHEDRON)    && (localFaces[i].elemType1 == HEXAHEDRON)) ||
         ((localFaces[i].elemType0 == QUADRILATERAL) && (localFaces[i].elemType1 == QUADRILATERAL)))
        mGridVar = mSolVar = 1;

      /*--- Store the grid data for this face. ---*/
      CUnsignedShort11T thisFace;
      thisFace.short0  = VTK_Type;
      thisFace.short1  = CFEMStandardElementBase::GetNIntStatic(VTK_Type, orderExact);
      thisFace.short2  = localFaces[i].elemType0;
      thisFace.short3  = localFaces[i].nPolyGrid0;
      thisFace.short4  = localFaces[i].faceID0;
      thisFace.short5  = localFaces[i].elemType1;
      thisFace.short6  = localFaces[i].nPolyGrid1;
      thisFace.short7  = localFaces[i].faceID1;
      thisFace.short8  = localFaces[i].DetermineOrientationElemSide1(volElem);
      thisFace.short9  = mGridVar;
      thisFace.short10 = orderExact;

      /*--- Store this face in faceTypesGrid. ---*/
      faceTypesGrid.push_back(thisFace);

      /*--- Adapt this face for the solution and add it to faceTypesSol. ---*/
      thisFace.short3 = localFaces[i].nPolySol0;
      thisFace.short6 = localFaces[i].nPolySol1;
      thisFace.short9 = mSolVar;

      faceTypesSol.push_back(thisFace);

      /*--- In case of an incompressible solver, reduce the polynomial
            degree of the solution by 1, set the number of variables to 1
            and add it to faceTypesSol. ---*/
      if( incompressible ) {
        thisFace.short3 = localFaces[i].nPolySol0-1;
        thisFace.short6 = localFaces[i].nPolySol1-1;
        thisFace.short9 = 1;

        faceTypesSol.push_back(thisFace);
      }
    }
  }

  /*--- Copy the face types for the grid, sort them and remove the
        double entities. ---*/
  vector<CUnsignedShort11T> typesFaceGrid = faceTypesGrid;
  vector<CUnsignedShort11T>::iterator lastEntry11T;
  sort(typesFaceGrid.begin(), typesFaceGrid.end());
  lastEntry11T = unique(typesFaceGrid.begin(), typesFaceGrid.end());
  typesFaceGrid.erase(lastEntry11T, typesFaceGrid.end());

  /*--- Copy the face types for the solution, sort them and remove the
        double entities. ---*/
  vector<CUnsignedShort11T> typesFaceSol = faceTypesSol;
  sort(typesFaceSol.begin(), typesFaceSol.end());
  lastEntry11T = unique(typesFaceSol.begin(), typesFaceSol.end());
  typesFaceSol.erase(lastEntry11T, typesFaceSol.end());

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Every standard face for an internal matching face is a     ---*/
  /*---         combination of two standard surface faces, i.e. one for    ---*/
  /*---         each side of the internal matching face. The surface face  ---*/
  /*---         takes care of all the work and the internal matching face  ---*/
  /*---         is just the correct combination of two surface faces.      ---*/
  /*---         This information is added to the already stored standard   ---*/
  /*---         faces needed to handle the boundary conditions.            ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the standard internal matching faces for the grid and
        add the two surfaces to typesSurfaceGrid. Note that the element
        on side 0 has an orientation of 0 per definition. ---*/
  for(unsigned int i=0; i<typesFaceGrid.size(); ++i) {

    /*--- Correct the number of variables for side 0 if the adjacent element
          is a hexahedron or quadrilateral. Add the surface standard element
          to typesSurfaceGrid afterwards. ---*/
    unsigned short mGridVar = typesFaceGrid[i].short9;
    if((typesFaceGrid[i].short2 == HEXAHEDRON) || (typesFaceGrid[i].short2 == QUADRILATERAL))
      mGridVar = 1;

    typesSurfaceGrid.push_back(CUnsignedShort8T(typesFaceGrid[i].short0, typesFaceGrid[i].short1,
                                                typesFaceGrid[i].short2, typesFaceGrid[i].short3,
                                                typesFaceGrid[i].short4, 0,
                                                mGridVar,                typesFaceGrid[i].short10));

    /*--- Correct the number of variables for side 1 if the adjacent element
          is a hexahedron or quadrilateral. Add the surface standard element
          to typesSurfaceGrid afterwards. ---*/
    mGridVar = typesFaceGrid[i].short9;
    if((typesFaceGrid[i].short5 == HEXAHEDRON) || (typesFaceGrid[i].short5 == QUADRILATERAL))
      mGridVar = 1;

    typesSurfaceGrid.push_back(CUnsignedShort8T(typesFaceGrid[i].short0, typesFaceGrid[i].short1,
                                                typesFaceGrid[i].short5, typesFaceGrid[i].short6,
                                                typesFaceGrid[i].short7, typesFaceGrid[i].short8,
                                                mGridVar,                typesFaceGrid[i].short10));
  }

  /*--- Sort typesSurfaceGrid and remove the double entities. ---*/
  sort(typesSurfaceGrid.begin(), typesSurfaceGrid.end());
  lastEntry8T = unique(typesSurfaceGrid.begin(), typesSurfaceGrid.end());
  typesSurfaceGrid.erase(lastEntry8T, typesSurfaceGrid.end());

  /*--- Loop over the standard internal matching faces for the solution and
        add the two surfaces to typesSurfaceSol. Note that the element
        on side 0 has an orientation of 0 per definition. ---*/
  for(unsigned int i=0; i<typesFaceSol.size(); ++i) {

    /*--- Correct the number of variables for side 0 if the adjacent element
          is a hexahedron or quadrilateral. Add the surface standard element
          to typesSurfaceSol afterwards. ---*/
    unsigned short mSolVar = typesFaceSol[i].short9;
    if((typesFaceSol[i].short2 == HEXAHEDRON) || (typesFaceSol[i].short2 == QUADRILATERAL))
      mSolVar = 1;

    typesSurfaceSol.push_back(CUnsignedShort8T(typesFaceSol[i].short0, typesFaceSol[i].short1,
                                               typesFaceSol[i].short2, typesFaceSol[i].short3,
                                               typesFaceSol[i].short4, 0,
                                               mSolVar,                typesFaceSol[i].short10));

    /*--- Correct the number of variables for side 1 if the adjacent element
          is a hexahedron or quadrilateral. Add the surface standard element
          to typesSurfaceSol afterwards. ---*/
    mSolVar = typesFaceSol[i].short9;
    if((typesFaceSol[i].short5 == HEXAHEDRON) || (typesFaceSol[i].short5 == QUADRILATERAL))
      mSolVar = 1;
    typesSurfaceSol.push_back(CUnsignedShort8T(typesFaceSol[i].short0, typesFaceSol[i].short1,
                                               typesFaceSol[i].short5, typesFaceSol[i].short6,
                                               typesFaceSol[i].short7, typesFaceSol[i].short8,
                                               mSolVar,                typesFaceSol[i].short10));
  }

  /*--- Sort typesSurfaceSol and remove the double entities. ---*/
  sort(typesSurfaceSol.begin(), typesSurfaceSol.end());
  lastEntry8T = unique(typesSurfaceSol.begin(), typesSurfaceSol.end());
  typesSurfaceSol.erase(lastEntry8T, typesSurfaceSol.end());

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Determine the different GEMM calls present in the standard ---*/
  /*---         surface elements. This data is stored in a vector of the   ---*/
  /*---         the class CUnsignedShort4T, which can store 4 short        ---*/
  /*---         integers as one entity. Find below the explanation of the  ---*/
  /*---         4 parameters that define a GEMM call.                      ---*/
  /*--------------------------------------------------------------------------*/

  /*--- 1: VTK type of the adjacent volume element. Needed to check whether a
           tensor product can be used or a full gemm call. If this is not a
           tensor product element, this variable contains the information
           whether it is DOF2Int or Int2DOF. This is only relevant for solution
           GEMM's, because it is needed to compute the residual.
        2: M, the first matrix dimension of A and C in the gemm call.
           When a tensor product is used, these are the first
           dimensions of the tensor A.
        3: K, the first matrix dimension of B and second matrix dimension
           of A. When a tensor product is used these are the first dimensions
           of tensor B and the last of tensor A.
        4: N, the last dimension of B and C in the gemm call. This is the
           number of variables per point, if relevant. For a tensor product
           GEMM call it contains the information whether it is DOF2Int
           or Int2DOF. The latter is only relevant for the solution
           GEMM's, because it is needed to compute the residual.
        Define the variable to store this information. ---*/
  vector<CUnsignedShort4T> gemmTypes;

  /*--- Loop over the entries of typesSurfaceGrid to fill gemmTypes.
        Note that when no tensor product is used, the element type
        is not relevant and it contains the information that this is a
        gemm for DOFS_TO_INT. Also the number of integration points for
        a quadrilateral face must be corrected, because currently the
        number of integration points in 1D is stored.
        When a tensor product is used, the last variable of gemmTypes
        indicates this is a gemm call from DOFS_TO_INT. ---*/ 
  for(unsigned long i=0; i<typesSurfaceGrid.size(); ++i) {

    unsigned short VTK_Type_Elem = typesSurfaceGrid[i].short2;
    unsigned short nDOFs, short3;
    unsigned short nInt = typesSurfaceGrid[i].short1;
    if((VTK_Type_Elem == HEXAHEDRON) || (VTK_Type_Elem == QUADRILATERAL)) {
      nDOFs  = typesSurfaceGrid[i].short3+1;
      short3 = CGemmBase::DOFS_TO_INT;
    }
    else {
      if(typesSurfaceGrid[i].short0 == QUADRILATERAL) nInt *= nInt;
      nDOFs  = CFEMStandardElementBase::GetNDOFsStatic(VTK_Type_Elem,
                                                       typesSurfaceGrid[i].short3);
      short3 = typesSurfaceGrid[i].short6;
      VTK_Type_Elem = CGemmBase::DOFS_TO_INT;
    }

    gemmTypes.push_back(CUnsignedShort4T(VTK_Type_Elem, nInt, nDOFs, short3));
  }

  /*--- Loop over the entries of typesSurfaceSol to fill gemmTypes.
        Note that two gemm types are created, which correspond to
        the interpolation to the integration points from the DOFs
        of the element and vice versa to obtain the residual
        contribution from the surface integral. Also here, when no
        tensor product is used, the element type is not relevant and
        it contains the information whether it is DOFS_TO_INT or
        INT_TO_DOFS. Also the number of integration points for a
        quadrilateral must be corrected. When a tensor product is used,
        the last variable of gemmTypes indicates whether it is a gemm
        call from DOFS_TO_INT or from INT_TO_DOFS. ---*/ 
  for(unsigned long i=0; i<typesSurfaceSol.size(); ++i) {

    unsigned short short3_0 = typesSurfaceSol[i].short6;
    unsigned short short3_1 = typesSurfaceSol[i].short6;

    unsigned short VTK_Type_Elem_0 = typesSurfaceSol[i].short2;
    unsigned short VTK_Type_Elem_1 = typesSurfaceSol[i].short2;
    unsigned short nDOFs;
    unsigned short nInt = typesSurfaceSol[i].short1;
    if((VTK_Type_Elem_0 == HEXAHEDRON) || (VTK_Type_Elem_0 == QUADRILATERAL)) {
      nDOFs    = typesSurfaceSol[i].short3+1;
      short3_0 = CGemmBase::DOFS_TO_INT;
      short3_1 = CGemmBase::INT_TO_DOFS;
    }
    else {
      if(typesSurfaceSol[i].short0 == QUADRILATERAL) nInt *= nInt;
      nDOFs = CFEMStandardElementBase::GetNDOFsStatic(typesSurfaceSol[i].short2,
                                                      typesSurfaceSol[i].short3);
      VTK_Type_Elem_0 = CGemmBase::DOFS_TO_INT;
      VTK_Type_Elem_1 = CGemmBase::INT_TO_DOFS;
    }

    gemmTypes.push_back(CUnsignedShort4T(VTK_Type_Elem_0, nInt, nDOFs, short3_0));
    gemmTypes.push_back(CUnsignedShort4T(VTK_Type_Elem_1, nDOFs, nInt, short3_1));
  }

  /*--- Sort gemmTypes and remove the double entities. ---*/
  vector<CUnsignedShort4T>::iterator lastEntry4T;
  sort(gemmTypes.begin(), gemmTypes.end());
  lastEntry4T = unique(gemmTypes.begin(), gemmTypes.end());
  gemmTypes.erase(lastEntry4T, gemmTypes.end());

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Create the objects that carry out the gemm functionality   ---*/
  /*---         for the faces.                                             ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop to create the gemm types for the faces. ---*/
  gemmTypesFaces.resize(gemmTypes.size(), nullptr);
  for(unsigned long i=0; i<gemmTypes.size(); ++i) {

    /*--- Abbreviate the variables for readability. ---*/
    const unsigned short VTK_Type_Elem = gemmTypes[i].short0;
    const unsigned short M             = gemmTypes[i].short1;
    const unsigned short K             = gemmTypes[i].short2;
    const unsigned short N             = gemmTypes[i].short3;

    /*--- Determine the VTK type of the element and allocate the
          memory for the correct gemm type. Note that for a QUADRILATERAL and
          HEXAHEDRON N contains the information whether it is DOFS_TO_INT or
          INT_TO_DOFS, while for the default element type this information is
          stored in VTK_Type_Elem. ---*/
    switch( VTK_Type_Elem ) {
      case QUADRILATERAL:
        gemmTypesFaces[i] = new CGemmFaceQuad(M, N, K);
        break;
      case HEXAHEDRON:
        gemmTypesFaces[i] = new CGemmFaceHex(M, N, K);
        break;
      default:
        gemmTypesFaces[i] = new CGemmStandard(M, N, K, VTK_Type_Elem);
        break;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 6: Create the standard elements for the faces, both for the   ---*/
  /*---         grid and the solution.                                     ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop to create the standard surface elements for the grid. ---*/
  standardSurfaceElementsGrid.resize(typesSurfaceGrid.size(), nullptr);
  for(unsigned long i=0; i<typesSurfaceGrid.size(); ++i) {

    /*--- Abbreviate the variables for readability. ---*/
    unsigned short VTK_Type_Face = typesSurfaceGrid[i].short0;
    unsigned short nInt          = typesSurfaceGrid[i].short1;
    unsigned short VTK_Type_Elem = typesSurfaceGrid[i].short2;
    unsigned short nPoly         = typesSurfaceGrid[i].short3;
    unsigned short faceID_Elem   = typesSurfaceGrid[i].short4;
    unsigned short orientation   = typesSurfaceGrid[i].short5;
    unsigned short nVarPerPoint  = typesSurfaceGrid[i].short6;
    unsigned short orderExact    = typesSurfaceGrid[i].short7;

    /*--- Correct the VTK type of the element if no tensor product
          is used, because that is not relevant and correct the number
          of integration points for a quadrilateral. If a tensor product
          is used, nVarPerPoint is not used. Instead this entry
          indicates the tensor product is from DOFS_TO_INT for the
          grid. Also determine the number of DOFs. Again the number
          of DOFs in 1D for tensor product elements. ---*/
    unsigned short nDOFs, short3;
    if((VTK_Type_Elem == HEXAHEDRON) || (VTK_Type_Elem == QUADRILATERAL)) {
      nDOFs  = typesSurfaceGrid[i].short3+1;
      short3 = CGemmBase::DOFS_TO_INT;
    }
    else {
      if(VTK_Type_Face == QUADRILATERAL) nInt *= nInt;
      nDOFs  = CFEMStandardElementBase::GetNDOFsStatic(VTK_Type_Elem, nPoly);
      short3 = nVarPerPoint;
      VTK_Type_Elem = CGemmBase::DOFS_TO_INT;
    }

    /*--- Create the data for the gemm type used in this standard
          surface element. ---*/
    CUnsignedShort4T thisGemmType(VTK_Type_Elem, nInt, nDOFs, short3);

    /*--- Find this gemm type in gemmTypes. ---*/
    vector<CUnsignedShort4T>::const_iterator low;
    low = lower_bound(gemmTypes.begin(), gemmTypes.end(), thisGemmType);
    const unsigned long ind = low - gemmTypes.begin();

    /*--- Determine the type of the surface element. ---*/
    switch( VTK_Type_Face ) {
      case LINE: {

        /*--- 2D simulation. Determine the type of the adjacent volume
              element and allocate the appropriate standard element.
              Note that typesSurfaceGrid[i].short2 must be used. ---*/
        switch( typesSurfaceGrid[i].short2 ) {
          case TRIANGLE:
            standardSurfaceElementsGrid[i] = new CFEMStandardLineAdjacentTriGrid(nPoly, orderExact, 
                                                                                 faceID_Elem,
                                                                                 orientation, useLGL,
                                                                                 gemmTypesFaces[ind]);
            break;

          case QUADRILATERAL:
            standardSurfaceElementsGrid[i] = new CFEMStandardLineAdjacentQuadGrid(nPoly, orderExact,
                                                                                  faceID_Elem,
                                                                                  orientation, useLGL,
                                                                                  gemmTypesFaces[ind]);
            break;

          default:
            SU2_MPI::Error(string("Unknown adjacent volume element. This should not happen."),
                           CURRENT_FUNCTION);
        }

        break;
      }

      case TRIANGLE: {

        /*--- Triangular face. Determine the type of the adjacent volume
              element and allocate the appropriate standard element.
              Note that typesSurfaceGrid[i].short2 must be used. ---*/
        switch( typesSurfaceGrid[i].short2 ) {
          case TETRAHEDRON:
            standardSurfaceElementsGrid[i] = new CFEMStandardTriAdjacentTetGrid(nPoly, orderExact,
                                                                                faceID_Elem,
                                                                                orientation, useLGL,
                                                                                gemmTypesFaces[ind]);
            break;

          case PYRAMID:
            standardSurfaceElementsGrid[i] = new CFEMStandardTriAdjacentPyraGrid(nPoly, orderExact,
                                                                                 faceID_Elem,
                                                                                 orientation, useLGL,
                                                                                 gemmTypesFaces[ind]);
            break;

          case PRISM:
            standardSurfaceElementsGrid[i] = new CFEMStandardTriAdjacentPrismGrid(nPoly, orderExact,
                                                                                  faceID_Elem,
                                                                                  orientation, useLGL,
                                                                                  gemmTypesFaces[ind]);
            break;

          default:
            SU2_MPI::Error(string("Unknown adjacent volume element. This should not happen."),
                           CURRENT_FUNCTION);
        }

        break;
      }

      case QUADRILATERAL: {

        /*--- Quadrilateral face. Determine the type of the adjacent volume
              element and allocate the appropriate standard element.
              Note that typesSurfaceGrid[i].short2 must be used. ---*/
        switch( typesSurfaceGrid[i].short2 ) {
          case HEXAHEDRON:
            standardSurfaceElementsGrid[i] = new CFEMStandardQuadAdjacentHexGrid(nPoly, orderExact,
                                                                                 faceID_Elem,
                                                                                 orientation, useLGL,
                                                                                 gemmTypesFaces[ind]);
            break;

          case PYRAMID:
            standardSurfaceElementsGrid[i] = new CFEMStandardQuadAdjacentPyraGrid(nPoly, orderExact,
                                                                                  faceID_Elem,
                                                                                  orientation, useLGL,
                                                                                  gemmTypesFaces[ind]);
            break;

          case PRISM:
            standardSurfaceElementsGrid[i] = new CFEMStandardQuadAdjacentPrismGrid(nPoly, orderExact,
                                                                                   faceID_Elem,
                                                                                   orientation, useLGL,
                                                                                   gemmTypesFaces[ind]);
            break;

          default:
            SU2_MPI::Error(string("Unknown adjacent volume element. This should not happen."),
                           CURRENT_FUNCTION);
        }

        break;
      }

      default: 
        SU2_MPI::Error(string("Unknown surface element. This should not happen."),
                       CURRENT_FUNCTION);
    }
  }

  /*--- Loop to create the standard surface elements for the solution. ---*/
  standardSurfaceElementsSolution.resize(typesSurfaceSol.size(), nullptr);
  for(unsigned long i=0; i<typesSurfaceSol.size(); ++i) {

    /*--- Abbreviate the variables for readability. ---*/
    unsigned short VTK_Type_Face   = typesSurfaceSol[i].short0;
    unsigned short nInt            = typesSurfaceSol[i].short1;
    unsigned short VTK_Type_Elem_0 = typesSurfaceSol[i].short2;
    unsigned short VTK_Type_Elem_1 = typesSurfaceSol[i].short2;
    unsigned short nPoly           = typesSurfaceSol[i].short3;
    unsigned short faceID_Elem     = typesSurfaceSol[i].short4;
    unsigned short orientation     = typesSurfaceSol[i].short5;
    unsigned short nVarPerPoint    = typesSurfaceSol[i].short6;
    unsigned short orderExact      = typesSurfaceSol[i].short7;

    /*--- Correct the VTK type of the element if no tensor product
          is used, because here the information DOFS_TO_INT or
          INT_TO_DOFS is stored. Also correct the number of integration
          points for a quadrilateral. If a tensor product is used,
          nVarPerPoint is not used. Instead this entry
          indicates whether the tensor product is from DOFS_TO_INT
          or from INT_TO_DOFS. Also determine the number of DOFs.
          Again the number of DOFs in 1D for tensor product elements. ---*/
    unsigned short nDOFs, short3_0, short3_1;
    if((VTK_Type_Elem_0 == HEXAHEDRON) || (VTK_Type_Elem_0 == QUADRILATERAL)) {
      nDOFs    = typesSurfaceSol[i].short3+1;
      short3_0 = CGemmBase::DOFS_TO_INT;
      short3_1 = CGemmBase::INT_TO_DOFS;
    }
    else {
      if(VTK_Type_Face == QUADRILATERAL) nInt *= nInt;
      nDOFs = CFEMStandardElementBase::GetNDOFsStatic(VTK_Type_Elem_0, nPoly);
      VTK_Type_Elem_0 = CGemmBase::DOFS_TO_INT;
      VTK_Type_Elem_1 = CGemmBase::INT_TO_DOFS;
      short3_0 = short3_1 = nVarPerPoint;
    }

    /*--- Create the data for the gemm type used in this standard
          surface element. ---*/
    CUnsignedShort4T thisGemmType(VTK_Type_Elem_0, nInt, nDOFs, short3_0);

    /*--- Find this gemm type in gemmTypes. ---*/
    vector<CUnsignedShort4T>::const_iterator low;
    low = lower_bound(gemmTypes.begin(), gemmTypes.end(), thisGemmType);
    const unsigned long ind1 = low - gemmTypes.begin();

    /*--- Also search  for the gemm type with the number of integration points
          and number of DOFs swapped. ---*/
    swap(thisGemmType.short1, thisGemmType.short2);
    thisGemmType.short0 = VTK_Type_Elem_1;
    thisGemmType.short3 = short3_1;
    low = lower_bound(gemmTypes.begin(), gemmTypes.end(), thisGemmType);
    const unsigned long ind2 = low - gemmTypes.begin();

    /*--- Determine the type of the surface element. ---*/
    switch( VTK_Type_Face ) {
      case LINE: {

        /*--- 2D simulation. Determine the type of the adjacent volume
              element and allocate the appropriate standard element.
              Note that typesSurfaceSol[i].short2 must be used. ---*/
        switch( typesSurfaceSol[i].short2 ) {
          case TRIANGLE:
            standardSurfaceElementsSolution[i] = new CFEMStandardLineAdjacentTriSol(nPoly, orderExact,
                                                                                    faceID_Elem, orientation,
                                                                                    gemmTypesFaces[ind1],
                                                                                    gemmTypesFaces[ind2]);
            break;

          case QUADRILATERAL:
            standardSurfaceElementsSolution[i] = new CFEMStandardLineAdjacentQuadSol(nPoly, orderExact,
                                                                                     faceID_Elem, orientation,
                                                                                     gemmTypesFaces[ind1],
                                                                                     gemmTypesFaces[ind2]);
            break;

          default:
            SU2_MPI::Error(string("Unknown adjacent volume element. This should not happen."),
                           CURRENT_FUNCTION);
        }

        break;
      }

      case TRIANGLE: {

        /*--- Triangular face. Determine the type of the adjacent volume
              element and allocate the appropriate standard element.
              Note that typesSurfaceSol[i].short2 must be used. ---*/
        switch( typesSurfaceSol[i].short2 ) {
          case TETRAHEDRON:
            standardSurfaceElementsSolution[i] = new CFEMStandardTriAdjacentTetSol(nPoly, orderExact,
                                                                                   faceID_Elem, orientation,
                                                                                   gemmTypesFaces[ind1],
                                                                                   gemmTypesFaces[ind2]);
            break;

          case PYRAMID:
            standardSurfaceElementsSolution[i] = new CFEMStandardTriAdjacentPyraSol(nPoly, orderExact,
                                                                                    faceID_Elem, orientation,
                                                                                    gemmTypesFaces[ind1],
                                                                                    gemmTypesFaces[ind2]);
            break;

          case PRISM:
            standardSurfaceElementsSolution[i] = new CFEMStandardTriAdjacentPrismSol(nPoly, orderExact,
                                                                                     faceID_Elem, orientation,
                                                                                     gemmTypesFaces[ind1],
                                                                                     gemmTypesFaces[ind2]);
            break;

          default:
            SU2_MPI::Error(string("Unknown adjacent volume element. This should not happen."),
                           CURRENT_FUNCTION);
        }

        break;
      }

      case QUADRILATERAL: {

        /*--- Quadrilateral face. Determine the type of the adjacent volume
              element and allocate the appropriate standard element.
              Note that typesSurfaceSol[i].short2 must be used. ---*/
        switch( typesSurfaceSol[i].short2 ) {
          case HEXAHEDRON:
            standardSurfaceElementsSolution[i] = new CFEMStandardQuadAdjacentHexSol(nPoly, orderExact,
                                                                                    faceID_Elem, orientation,
                                                                                    gemmTypesFaces[ind1],
                                                                                    gemmTypesFaces[ind2]);
            break;

          case PYRAMID:
            standardSurfaceElementsSolution[i] = new CFEMStandardQuadAdjacentPyraSol(nPoly, orderExact,
                                                                                     faceID_Elem, orientation,
                                                                                     gemmTypesFaces[ind1],
                                                                                     gemmTypesFaces[ind2]);
            break;

          case PRISM:
            standardSurfaceElementsSolution[i] = new CFEMStandardQuadAdjacentPrismSol(nPoly, orderExact,
                                                                                      faceID_Elem, orientation,
                                                                                      gemmTypesFaces[ind1],
                                                                                      gemmTypesFaces[ind2]);
            break;

          default:
            SU2_MPI::Error(string("Unknown adjacent volume element. This should not happen."),
                           CURRENT_FUNCTION);
        }

        break;
      }

      default: 
        SU2_MPI::Error(string("Unknown surface element. This should not happen."),
                       CURRENT_FUNCTION);
    }

    /*--- Allocate the memory for the working variables. ---*/
    standardSurfaceElementsSolution[i]->AllocateWorkingVariables(nDim, nAllocVar, true);
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 7: Set the pointers to standard surface elements for the      ---*/
  /*---         faces belonging to physical boundaries.                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the physical boundaries. Exclude the periodic boundaries. ---*/
  unsigned long indGrid = 0, indSol = 0;
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    if( !boundaries[iMarker].periodicBoundary ) {

      /*--- Easier storage of the surfaces of this boundary and loop over them. ---*/
      vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
      for(unsigned long i=0; i<surfElem.size(); ++i) {

        /*--- Determine the index in the standard element for the grid that corresponds
              to this surface element and set the pointer accordingly. ---*/
        vector<CUnsignedShort8T>::const_iterator low;
        low = lower_bound(typesSurfaceGrid.begin(), typesSurfaceGrid.end(),
                          surfaceTypesGrid[indGrid++]);
        unsigned long ind = low - typesSurfaceGrid.begin();

        surfElem[i].standardElemGrid = standardSurfaceElementsGrid[ind];

        /*--- Determine the index in the standard element for the solution that corresponds
              to this surface element and set the pointer accordingly. ---*/
        low = lower_bound(typesSurfaceSol.begin(), typesSurfaceSol.end(),
                          surfaceTypesSol[indSol++]);
        ind = low - typesSurfaceSol.begin();

        surfElem[i].standardElemFlow = standardSurfaceElementsSolution[ind];

        /*--- When the incompressible equations are solved, also the standard element
              that corresponds to the pressure solution must be set. ---*/
        if( incompressible ) {
          low = lower_bound(typesSurfaceSol.begin(), typesSurfaceSol.end(),
                            surfaceTypesSol[indSol++]);
          ind = low - typesSurfaceSol.begin();

          surfElem[i].standardElemP = standardSurfaceElementsSolution[ind];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 8: Create the standard elements for internal matching faces.  ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop to create the standard elements for the internal matching
        faces for the grid. ---*/
  standardInternalFaceGrid.resize(typesFaceGrid.size(), nullptr);
  for(unsigned long i=0; i<typesFaceGrid.size(); ++i) {

    /*--- Create the data for the standard surface elements on side 0
          and side 1 of the internal matching face. ---*/
    CUnsignedShort8T surf0(typesFaceGrid[i].short0, typesFaceGrid[i].short1,
                           typesFaceGrid[i].short2, typesFaceGrid[i].short3,
                           typesFaceGrid[i].short4, 0,
                           typesFaceGrid[i].short9, typesFaceGrid[i].short10);
    CUnsignedShort8T surf1(typesFaceGrid[i].short0, typesFaceGrid[i].short1,
                           typesFaceGrid[i].short5, typesFaceGrid[i].short6,
                           typesFaceGrid[i].short7, typesFaceGrid[i].short8,
                           typesFaceGrid[i].short9, typesFaceGrid[i].short10);

    /*--- Correct the number of variables of the adjacent element is a
          quadrilateral or hexahedron, because this value does not matter. ---*/
    if((surf0.short2 == HEXAHEDRON) || (surf0.short2 == QUADRILATERAL)) surf0.short6 = 1;
    if((surf1.short2 == HEXAHEDRON) || (surf1.short2 == QUADRILATERAL)) surf1.short6 = 1;

    /*--- Determine the indices of surf0 and surf1 in the vector of standard
          elements for the surfaces, which are the same indices as in the
          vector typesSurfaceGrid. ---*/
    vector<CUnsignedShort8T>::const_iterator low;
    low = lower_bound(typesSurfaceGrid.begin(), typesSurfaceGrid.end(), surf0);
    const unsigned long ind0 = low - typesSurfaceGrid.begin();

    low = lower_bound(typesSurfaceGrid.begin(), typesSurfaceGrid.end(), surf1);
    const unsigned long ind1 = low - typesSurfaceGrid.begin();

    /*--- Create the standard element. ---*/
    standardInternalFaceGrid[i] = new CFEMStandardInternalFaceGrid(standardSurfaceElementsGrid[ind0],
                                                                   standardSurfaceElementsGrid[ind1]);
  }

  /*--- Loop to create the standard elements for the internal matching
        faces for the solution. ---*/
  standardInternalFaceSolution.resize(typesFaceSol.size(), nullptr);
  for(unsigned long i=0; i<typesFaceSol.size(); ++i) {

    /*--- Create the data for the standard surface elements on side 0
          and side 1 of the internal matching face. ---*/
    CUnsignedShort8T surf0(typesFaceSol[i].short0, typesFaceSol[i].short1,
                           typesFaceSol[i].short2, typesFaceSol[i].short3,
                           typesFaceSol[i].short4, 0,
                           typesFaceSol[i].short9, typesFaceSol[i].short10);
    CUnsignedShort8T surf1(typesFaceSol[i].short0, typesFaceSol[i].short1,
                           typesFaceSol[i].short5, typesFaceSol[i].short6,
                           typesFaceSol[i].short7, typesFaceSol[i].short8,
                           typesFaceSol[i].short9, typesFaceSol[i].short10);

    /*--- Correct the number of variables if the adjacent element is a
          quadrilateral or hexahedron, because this value does not matter. ---*/
    if((surf0.short2 == HEXAHEDRON) || (surf0.short2 == QUADRILATERAL)) surf0.short6 = 1;
    if((surf1.short2 == HEXAHEDRON) || (surf1.short2 == QUADRILATERAL)) surf1.short6 = 1;

    /*--- Determine the indices of surf0 and surf1 in the vector of standard
          elements for the surfaces, which are the same indices as in the
          vector typesSurfaceSol. ---*/
    vector<CUnsignedShort8T>::const_iterator low;
    low = lower_bound(typesSurfaceSol.begin(), typesSurfaceSol.end(), surf0);
    const unsigned long ind0 = low - typesSurfaceSol.begin();

    low = lower_bound(typesSurfaceSol.begin(), typesSurfaceSol.end(), surf1);
    const unsigned long ind1 = low - typesSurfaceSol.begin();

    /*--- Create the standard element. ---*/
    standardInternalFaceSolution[i] = new CFEMStandardInternalFaceSol(standardSurfaceElementsSolution[ind0],
                                                                      standardSurfaceElementsSolution[ind1]);
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 9: Set the pointers to standard internal faces for the        ---*/
  /*---         internal matching faces.                                   ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the internal matching faces. ---*/
  indSol = 0;
  for(unsigned long i=0; i<matchingFaces.size(); ++i) {

    /*--- Determine the index of faceTypesGrid in typesFaceGrid, which is also
          the index of the matching face in the standard matching faces.
          Set the pointer afterwards. ---*/
    vector<CUnsignedShort11T>::const_iterator low;
    low = lower_bound(typesFaceGrid.begin(), typesFaceGrid.end(), faceTypesGrid[i]);
    unsigned long ind = low - typesFaceGrid.begin();

    matchingFaces[i].standardElemGrid = standardInternalFaceGrid[ind];

    /*--- Determine the index of faceTypesSol in typesFaceSol, which is also
          the index of the matching face in the standard matching faces.
          Set the pointer afterwards. Note that it is possible that there is also
          a standard face for the pressure solution, hence the index indSol must
          be used instead of i. ---*/
    low = lower_bound(typesFaceSol.begin(), typesFaceSol.end(), faceTypesSol[indSol++]);
    ind = low - typesFaceSol.begin();

    matchingFaces[i].standardElemFlow = standardInternalFaceSolution[ind];

    /*--- When the incompressible equations are solved, also the standard element
          that corresponds to the pressure solution must be set. ---*/
    if( incompressible ) {
      low = lower_bound(typesFaceSol.begin(), typesFaceSol.end(), faceTypesSol[indSol++]);
      ind = low - typesFaceSol.begin();

      matchingFaces[i].standardElemP = standardInternalFaceSolution[ind];
    }
  }
}

void CMeshFEM_DG::CreateStandardVolumeElementsSolution(const vector<CUnsignedShort4T> &elemTypes,
                                                       const unsigned short           nAllocVar,
                                                       const unsigned short           locGridDOFs) {

  /*--- The master node writes a message. ---*/
  if(rank == MASTER_NODE) cout << "Creating standard volume elements." << endl;

  /*--- Allocate the memory for the pointers. ---*/
  standardVolumeElementsSolution.resize(elemTypes.size(), nullptr);

  /*--- Loop over the different element types for the grid. ---*/
  for(unsigned long i=0; i<elemTypes.size(); ++i) {

    /*--- Abbreviate the element type, polynomial degree, polynomial order that must be
          integrated exactly, and the number of solution variables for readability. ---*/
    const unsigned short VTK_Type   = elemTypes[i].short0;
    const unsigned short nPoly      = elemTypes[i].short1;
    const unsigned short orderExact = elemTypes[i].short2;
    const unsigned short nSolVar    = elemTypes[i].short3;

    /*--- Determine the element type and allocate the appropriate object. ---*/
    switch( VTK_Type ) {
      case TRIANGLE:
        standardVolumeElementsSolution[i] = new CFEMStandardTriVolumeSol(nPoly, orderExact, locGridDOFs,
                                                                         nSolVar);
        break;
      case QUADRILATERAL:
        standardVolumeElementsSolution[i] = new CFEMStandardQuadVolumeSol(nPoly, orderExact, locGridDOFs,
                                                                          nSolVar);
        break;
      case TETRAHEDRON:
        standardVolumeElementsSolution[i] = new CFEMStandardTetVolumeSol(nPoly, orderExact, locGridDOFs,
                                                                         nSolVar);
        break;
      case PYRAMID:
        standardVolumeElementsSolution[i] = new CFEMStandardPyraVolumeSol(nPoly, orderExact, locGridDOFs,
                                                                          nSolVar);
        break;
      case PRISM:
        standardVolumeElementsSolution[i] = new CFEMStandardPrismVolumeSol(nPoly, orderExact, locGridDOFs,
                                                                           nSolVar);
        break;
      case HEXAHEDRON:
        standardVolumeElementsSolution[i] = new CFEMStandardHexVolumeSol(nPoly, orderExact, locGridDOFs,
                                                                         nSolVar);
        break;
      default:  /*--- To avoid a compiler warning. ---*/
        SU2_MPI::Error(string("Unknown volume element. This should not happen"),
                       CURRENT_FUNCTION);
    }

    /*--- Allocate the memory for the working variables. ---*/
    standardVolumeElementsSolution[i]->AllocateWorkingVariables(nDim, nAllocVar, false);
  }
}

void CMeshFEM_DG::HighOrderContainmentSearch(const su2double      *coor,
                                             const unsigned long  parElem,
                                             const unsigned short subElem,
                                             const su2double      *weightsSubElem,
                                             su2double            *parCoor) {

  /*--- Definition of the maximum number of iterations in the Newton solver
        and the tolerance level. ---*/
  const unsigned short maxIt = 50;
  const su2double tolNewton  = 1.e-10;

  /*--- Create an initial guess for the parametric coordinates from
        interpolation in the linear sub-element of the parent element. ---*/
  volElem[parElem].standardElemGrid->InterpolCoorSubElem(subElem, weightsSubElem,
                                                         parCoor);

  /*--- Start of the Newton algorithm. Loop over the number of iterations. ---*/
  unsigned short itCount;
  for(itCount=0; itCount<maxIt; ++itCount) {

    /*--- Determine the coordinates and its gradients for the current
          values of the parametric coordinates. ---*/
    su2double x[3], dxdpar[3][3];
    volElem[parElem].standardElemGrid->EvalCoorAndGradCoor(volElem[parElem].coorGridDOFs,
                                                           parCoor, x, dxdpar);

    /*--- Compute the update, for which a distinction between
          2D and 3D must be made. ---*/
    bool converged = false;
    switch( nDim ) {
      case 2: {
        /*--- Two dimensional computation. Compute the values of the function
              and minus the Jacobian matrix. ---*/
        const su2double f0  = coor[0]-x[0], f1  = coor[1]-x[1];
        const su2double a00 = dxdpar[0][0], a01 = dxdpar[0][1],
                        a10 = dxdpar[1][0], a11 = dxdpar[1][1];

        /*--- Compute the updates of the parametric values. As minus the
              Jacobian is computed, the updates should be added to parCoor. ---*/
        const su2double detInv = 1.0/(a00*a11 - a01*a10);
        const su2double dr = detInv*(f0*a11 - f1*a01);
        const su2double ds = detInv*(f1*a00 - f0*a10);

        parCoor[0] += dr;
        parCoor[1] += ds;

        /*--- Check for convergence. ---*/
        if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton) converged = true;
        break;
      }

      case 3: {
        /*--- Three dimensional computation. Compute the values of the function
              and minus the Jacobian matrix. ---*/
        const su2double f0  = coor[0]-x[0], f1  = coor[1]-x[1], f2  = coor[2]-x[2];
        const su2double a00 = dxdpar[0][0], a01 = dxdpar[0][1], a02 = dxdpar[0][2],
                        a10 = dxdpar[1][0], a11 = dxdpar[1][1], a12 = dxdpar[1][2],
                        a20 = dxdpar[2][0], a21 = dxdpar[2][1], a22 = dxdpar[2][2];

        /*--- Compute the updates of the parametric values. As minus the
              Jacobian is computed, the updates should be added to parCoor. ---*/
        const su2double detInv = 1.0/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22
                               +      a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
        const su2double dr =  detInv*(a01*a12*f2 - a01*a22*f1 - a02*a11*f2
                           +          a02*a21*f1 + a11*a22*f0 - a12*a21*f0);
        const su2double ds = -detInv*(a00*a12*f2 - a00*a22*f1 - a02*a10*f2
                           +          a02*a20*f1 + a10*a22*f0 - a12*a20*f0);
        const su2double dt =  detInv*(a00*a11*f2 - a00*a21*f1 - a01*a10*f2
                           +          a01*a20*f1 + a10*a21*f0 - a11*a20*f0);
        parCoor[0] += dr;
        parCoor[1] += ds;
        parCoor[2] += dt;

        /* Check for convergence. */
        if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton)
          converged = true;
        break;
      }
    }

    /* Break the loop if the Newton algorithm converged. */
    if( converged ) break;
  }

  /*--- Terminate if the Newton algorithm did not converge. ---*/
  if(itCount == maxIt) SU2_MPI::Error(string("Newton did not converge"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::SetWallDistance(su2double val) {

  /*---- Start of the OpenMP parallel region, if supported. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Set the wall distance for the owned volume elements and
          the internal matching faces. ---*/
#ifdef HAVE_OMP
    const size_t omp_chunk_size_elem = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
    SU2_OMP(for schedule(static,omp_chunk_size_elem) SU2_NOWAIT)
    for(unsigned long l=0; l<nVolElemOwned; ++l)
      volElem[l].SetWallDistance(val);

#ifdef HAVE_OMP
    const size_t omp_chunk_size_face = computeStaticChunkSize(matchingFaces.size(), omp_get_num_threads(), 64);
#endif
    SU2_OMP(for schedule(static,omp_chunk_size_face) SU2_NOWAIT)
    for(unsigned long l=0; l<matchingFaces.size(); ++l)
      matchingFaces[l].SetWallDistance(val);

    /*--- Set the wall distance for the physical boundary markers. ---*/
    for(unsigned short iMarker=0; iMarker<boundaries.size(); ++iMarker) {
      if( !boundaries[iMarker].periodicBoundary ) {

        vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
#ifdef HAVE_OMP
        const size_t omp_chunk_size_surf = computeStaticChunkSize(surfElem.size(), omp_get_num_threads(), 64);
#endif
        SU2_OMP(for schedule(static,omp_chunk_size_surf) SU2_NOWAIT)
        for(unsigned long l=0; l<surfElem.size(); ++l)
          surfElem[l].SetWallDistance(val);
      }
    }
  }
  END_SU2_OMP_PARALLEL
}

void CMeshFEM_DG::SetWallDistance(CADTElemClass* WallADT,
                                  const CConfig *config,
                                  unsigned short iZone) {

  /*--- Return immediately if the tree is empty. ---*/
  if( WallADT->IsEmpty() ) return;

  /*---- Start of the OpenMP parallel region, if supported. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Determine the wall distances of the integration points and
          solution DOFs of the locally owned volume elements and of
          the integration points of the internal matching faces. ---*/
#ifdef HAVE_OMP
    const size_t omp_chunk_size_elem = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
    SU2_OMP(for schedule(static,omp_chunk_size_elem) SU2_NOWAIT)
    for(unsigned long l=0; l<nVolElemOwned; ++l)
      volElem[l].ComputeWallDistance(WallADT, nDim);

#ifdef HAVE_OMP
    const size_t omp_chunk_size_face = computeStaticChunkSize(matchingFaces.size(), omp_get_num_threads(), 64);
#endif
    SU2_OMP(for schedule(static,omp_chunk_size_face) SU2_NOWAIT)
    for(unsigned long l=0; l<matchingFaces.size(); ++l)
      matchingFaces[l].ComputeWallDistance(WallADT, nDim);

    /*--- Determine the wall distances of the integration points
          of the physical boundary markers. ---*/
    for(unsigned short iMarker=0; iMarker<boundaries.size(); ++iMarker) {
      if( !boundaries[iMarker].periodicBoundary ) {

        vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
#ifdef HAVE_OMP
        const size_t omp_chunk_size_surf = computeStaticChunkSize(surfElem.size(), omp_get_num_threads(), 64);
#endif
        SU2_OMP(for schedule(static,omp_chunk_size_surf) SU2_NOWAIT)
        for(unsigned long l=0; l<surfElem.size(); ++l)
          surfElem[l].ComputeWallDistance(WallADT, nDim);
      }
    }
  }
  END_SU2_OMP_PARALLEL
}

void CMeshFEM_DG::TimeCoefficientsPredictorADER_DG(CConfig *config) {

  /*--------------------------------------------------------------------------*/
  /*--- Determine the interpolation matrix from the time DOFs of ADER-DG   ---*/
  /*--- to the time integration points. Note that in many cases this       ---*/
  /*--- interpolation matrix is the identity matrix, i.e. the DOFs and the ---*/
  /*--- integration points coincide.                                       ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Determine the number of time DOFs in the predictor step of ADER as well
        as their location on the interval [-1..1]. ---*/
  const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();
  const passivedouble  *TimeDOFs = config->GetTimeDOFsADER_DG();

  /*--- Determine the number of time integration points of ADER as well as
        their location on the interval [-1..1] and their integration weights. ---*/
  unsigned short       nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();
  const passivedouble *TimeIntegrationPoints  = config->GetTimeIntegrationADER_DG();
  const passivedouble *TimeIntegrationWeights = config->GetWeightsIntegrationADER_DG();

  /*--- Store the time DOFs and the time integration points in a vector. ---*/
  vector<passivedouble> rTimeDOFs(nTimeDOFs);
  for(unsigned short i=0; i<nTimeDOFs; ++i)
    rTimeDOFs[i] = TimeDOFs[i];

  vector<passivedouble> rTimeIntPoints(nTimeIntegrationPoints);
  for(unsigned short i=0; i<nTimeIntegrationPoints; ++i)
    rTimeIntPoints[i] = TimeIntegrationPoints[i];

  /*--- Determine the Lagrangian basis functions in the time integration points. ---*/
  CFEMStandardLineBase timeElement;
  timeElement.LagBasisIntPointsLine(rTimeDOFs, rTimeIntPoints, false,
                                    timeInterpolDOFToIntegrationADER_DG);

  /*--------------------------------------------------------------------------*/
  /*--- Determine the coefficients that appear in the iteration matrix of  ---*/
  /*--- the prediction step of ADER-DG.                                    ---*/
  /*--------------------------------------------------------------------------*/
  
  /*--- Determine the Lagrangian interpolation functions for r == 1, i.e.
        the end of the time interval. ---*/
  vector<passivedouble> rEnd(1); rEnd[0] = 1.0;
  ColMajorMatrix<passivedouble> lEnd;
  timeElement.LagBasisIntPointsLine(rTimeDOFs, rEnd, false, lEnd);

  /*--- Determine the Lagrangian basis functions in the integration points. ---*/
  ColMajorMatrix<passivedouble> derLInt;
  timeElement.DerLagBasisIntPointsLine(rTimeDOFs, rTimeIntPoints, false, derLInt);

  /*--- Compute the time coefficients of the iteration matrix in
        the predictor step. ---*/
  timeCoefADER_DG.Initialize(nTimeDOFs);
  for(unsigned short j=0; j<nTimeDOFs; ++j) {
    for(unsigned short i=0; i<nTimeDOFs; ++i) {
      timeCoefADER_DG(i,j) = lEnd(0,i)*lEnd(0,j);
      for(unsigned short k=0; k<nTimeIntegrationPoints; ++k)
        timeCoefADER_DG(i,j) -= TimeIntegrationWeights[k]*derLInt(k,i)
                              * timeInterpolDOFToIntegrationADER_DG(k,j);
    }
  }

  /*--- Determine the inverse of timeCoefADER_DG, because this is needed in the
        iteration matrix of the predictor step of ADER-DG. ---*/
  timeCoefADER_DG.Invert();

  /*--------------------------------------------------------------------------*/
  /*--- Determine the interpolation matrix from the time DOFs of ADER-DG   ---*/
  /*--- of adjacent elements of a higher time level to the time            ---*/
  /*--- integration points.                                                ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Determine the location of the time integration points from the
        perspective of an adjacent element with twice the time step. Note that
        the time step of this adjacent element covers two fine time steps and
        hence the number of integration points double and the position is
        determined via the transformation xxi = 0.5(xi-1) and xi varies between
        -1 and 3. ---*/
  rTimeIntPoints.resize(2*nTimeIntegrationPoints);
  for(unsigned short i=0; i<nTimeIntegrationPoints; ++i)
    rTimeIntPoints[i] = 0.5*(TimeIntegrationPoints[i]-1.0);

  for(unsigned short i=0; i<nTimeIntegrationPoints; ++i) {
    const unsigned short ii = i + nTimeIntegrationPoints;
    rTimeIntPoints[ii] = 0.5*(TimeIntegrationPoints[i]+1.0);
  }

  /*--- Determine the Lagrangian basis functions in these time integration points. ---*/
  timeElement.LagBasisIntPointsLine(rTimeDOFs, rTimeIntPoints, false,
                                    timeInterpolAdjDOFToIntegrationADER_DG);
}
