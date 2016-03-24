/*!
 * \file fem_geometry_structure.cpp
 * \brief Main subroutines for creating the primal grid for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/fem_geometry_structure.hpp"

void CPointFEM::Copy(const CPointFEM &other) {
  globalID = other.globalID;
  coor[0]  = other.coor[0];
  coor[1]  = other.coor[1];
  coor[2]  = other.coor[2];
}

void CSurfaceElementFEM::Copy(const CSurfaceElementFEM &other) {
  VTK_Type           = other.VTK_Type;
  nPolyGrid          = other.nPolyGrid;
  nDOFsGrid          = other.nDOFsGrid;
  indStandardElement = other.indStandardElement;
  volElemID          = other.volElemID;
  boundElemIDGlobal  = other.boundElemIDGlobal;
  nodeIDsGrid        = other.nodeIDsGrid;
}

CMeshFEM::CMeshFEM(CGeometry *geometry, CConfig *config) {

  /*--- Determine the number of ranks and the current rank. ---*/
  int nRank = SINGLE_NODE;
  int rank  = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRank);
#endif

  /*--- Copy the number of dimensions. ---*/
  nDim = geometry->GetnDim();

  /*--- Determine a mapping from the global point ID to the local index
        of the points.            ---*/
  map<unsigned long,unsigned long> globalPointIDToLocalInd;
  for(unsigned i=0; i<geometry->local_node; ++i)
    globalPointIDToLocalInd[geometry->node[i]->GetGlobalIndex()] = i;

  /*----------------------------------------------------------------------------*/
  /*--- Step 1: Communicate the elements and the boundary elements to the    ---*/
  /*---         ranks where they will stored during the computation.         ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Determine the ranks to which I have to send my elements. ---*/
  vector<int> sendToRank(nRank, 0);

  for(unsigned long i=0; i<geometry->local_elem; i++) {
    sendToRank[geometry->elem[i]->GetColor()] = 1;
  }

  map<int,int> rankToIndCommBuf;
  for(int i=0; i<nRank; ++i) {
    if( sendToRank[i] ) {
      int ind = rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  /*--- Definition of the communication buffers, used to send the element data
        to the correct ranks.                ---*/
  int nRankSend = rankToIndCommBuf.size();
  vector<vector<short> >     shortSendBuf(nRankSend,  vector<short>(0));
  vector<vector<long>  >     longSendBuf(nRankSend,   vector<long>(0));
  vector<vector<su2double> > doubleSendBuf(nRankSend, vector<su2double>(0));

  /*--- The first element of longSendBuf will contain the number of elements, which
        are stored in the communication buffers. Initialize this value to 0. ---*/
  for(int i=0; i<nRankSend; ++i) longSendBuf[i].push_back(0);

  /*--- Determine the number of ranks, from which this rank will receive elements. ---*/
  int nRankRecv = 1;

#ifdef HAVE_MPI
  vector<int> sizeRecv(nRank, 1);

  MPI_Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(),
                     MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Loop over the local elements to fill the communication buffers with element data. ---*/
  for(unsigned long i=0; i<geometry->local_elem; i++) {
    int ind = geometry->elem[i]->GetColor();
    map<int,int>::const_iterator MI = rankToIndCommBuf.find(ind);
    ind = MI->second;

    ++longSendBuf[ind][0];   /* The number of elements in the buffers must be incremented. */

    shortSendBuf[ind].push_back(geometry->elem[i]->GetVTK_Type());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNPolyGrid());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNPolySol());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNDOFsGrid());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNDOFsSol());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetnFaces());
    shortSendBuf[ind].push_back( (short) geometry->elem[i]->GetJacobianConsideredConstant());

    longSendBuf[ind].push_back(geometry->elem[i]->GetGlobalElemID());
    longSendBuf[ind].push_back(geometry->elem[i]->GetGlobalOffsetDOFsSol());

    for(unsigned short j=0; j<geometry->elem[i]->GetNDOFsGrid(); ++j)
      longSendBuf[ind].push_back(geometry->elem[i]->GetNode(j));

    for(unsigned short j=0; j<geometry->elem[i]->GetnFaces(); ++j)
      longSendBuf[ind].push_back(geometry->elem[i]->GetNeighbor_Elements(j));

    for(unsigned short j=0; j<geometry->elem[i]->GetnFaces(); ++j) {
      shortSendBuf[ind].push_back(geometry->elem[i]->GetPeriodicIndex(j));
      shortSendBuf[ind].push_back( (short) geometry->elem[i]->GetJacobianConstantFace(j));
    }
  }

  /*--- Determine for each rank to which I have to send elements the data of
        the corresponding nodes.   ---*/
  for(int i=0; i<nRankSend; ++i) {

    /*--- Determine the vector with node IDs in the connectivity
          of the elements for this rank.   ---*/
    vector<long> nodeIDs;

    unsigned long indL = 3;
    unsigned long indS = 3;
    for(unsigned long j=0; j<longSendBuf[i][0]; ++j) {
      short nDOFsGrid = shortSendBuf[i][indS], nFaces = shortSendBuf[i][indS+2];
      indS += 2*nFaces+7;

      for(short k=0; k<nDOFsGrid; ++k, ++indL)
        nodeIDs.push_back(longSendBuf[i][indL]);
      indL += nFaces+2;
    }

    /*--- Sort nodeIDs in increasing order and remove the double entities. ---*/
    sort(nodeIDs.begin(), nodeIDs.end());
    vector<long>::iterator lastNodeID = unique(nodeIDs.begin(), nodeIDs.end());
    nodeIDs.erase(lastNodeID, nodeIDs.end());

    /*--- Add the number of node IDs and the node IDs itself to longSendBuf[i]. ---*/
    longSendBuf[i].push_back(nodeIDs.size());
    longSendBuf[i].insert(longSendBuf[i].end(), nodeIDs.begin(), nodeIDs.end());

    /*--- Copy the coordinates to doubleSendBuf. ---*/
    for(unsigned long j=0; j<nodeIDs.size(); ++j) {
      map<unsigned long,unsigned long>::const_iterator LMI;
      LMI = globalPointIDToLocalInd.find(nodeIDs[j]);

      if(LMI == globalPointIDToLocalInd.end()) {
        cout << "Entry not found in map in function CMeshFEM::CMeshFEM" << endl;
#ifndef HAVE_MPI
        exit(EXIT_FAILURE);
#else
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#endif
      }

      unsigned long ind = LMI->second;
      for(unsigned short l=0; l<nDim; ++l)
        doubleSendBuf[i].push_back(geometry->node[ind]->GetCoord(l));
    }
  }

  /*--- Loop over the boundaries to send the boundary data to the appropriate rank. ---*/
  nMarker = geometry->GetnMarker();
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /* Store the current indices in the longSendBuf, which are used to store the
       number of boundary elements sent to this rank. Initialize this value to 0. */
    vector<long> indLongBuf(nRankSend);
    for(int i=0; i<nRankSend; ++i) {
      indLongBuf[i] = longSendBuf[i].size();
      longSendBuf[i].push_back(0);
    }

    /* Loop over the local boundary elements in geometry for this marker. */
    for(unsigned long i=0; i<geometry->GetnElem_Bound(iMarker); ++i) {

      /* Determine the local ID of the corresponding domain element. */
      unsigned long elemID = geometry->bound[iMarker][i]->GetDomainElement()
                           - geometry->starting_node[rank];

      /* Determine to which rank this boundary element must be sent.
         That is the same as its corresponding domain element.
         Update the corresponding index in longSendBuf. */
      int ind = geometry->elem[elemID]->GetColor();
      map<int,int>::const_iterator MI = rankToIndCommBuf.find(ind);
      ind = MI->second;

      ++longSendBuf[ind][indLongBuf[ind]];

      /* Store the data for this boundary element in the communication buffers. */
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetVTK_Type());
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNPolyGrid());
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNDOFsGrid());

      longSendBuf[ind].push_back(elemID);
      longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetGlobalElemID());

      for(unsigned short j=0; j<geometry->bound[iMarker][i]->GetNDOFsGrid(); ++j)
        longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNode(j));
    }
  }

  /*--- Definition of the communication buffers, used to receive
        the element data from the other correct ranks.        ---*/
  vector<vector<short> >     shortRecvBuf(nRankRecv,  vector<short>(0));
  vector<vector<long>  >     longRecvBuf(nRankRecv,   vector<long>(0));
  vector<vector<su2double> > doubleRecvBuf(nRankRecv, vector<su2double>(0));

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/

#ifdef HAVE_MPI

  /*--- Parallel mode. Send all the data using non-blocking sends. ---*/
  vector<MPI_Request> commReqs(3*nRankSend);
  map<int,int>::const_iterator MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {

    int dest = MI->first;
    SU2_MPI::Isend(shortSendBuf[i].data(), shortSendBuf[i].size(), MPI_SHORT,
                   dest, dest, MPI_COMM_WORLD, &commReqs[3*i]);
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest+1, MPI_COMM_WORLD, &commReqs[3*i+1]);
    SU2_MPI::Isend(doubleSendBuf[i].data(), doubleSendBuf[i].size(), MPI_DOUBLE,
                   dest, dest+2, MPI_COMM_WORLD, &commReqs[3*i+2]);
  }

  /* Loop over the number of ranks from which I receive data. */
  for(int i=0; i<nRankRecv; ++i) {

    /* Block until a message with shorts arrives from any processor.
       Determine the source and the size of the message.   */
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    MPI_Get_count(&status, MPI_SHORT, &sizeMess);

    /* Allocate the memory for the short receive buffer and receive the message. */
    shortRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(shortRecvBuf[i].data(), sizeMess, MPI_SHORT,
                  source, rank, MPI_COMM_WORLD, &status);

    /* Block until the corresponding message with longs arrives, determine
       its size, allocate the memory and receive the message. */
    MPI_Probe(source, rank+1, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_LONG, &sizeMess);
    longRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(longRecvBuf[i].data(), sizeMess, MPI_LONG,
                  source, rank+1, MPI_COMM_WORLD, &status);

    /* Idem for the message with doubles. */
    MPI_Probe(source, rank+2, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &sizeMess);
    doubleRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(doubleRecvBuf[i].data(), sizeMess, MPI_DOUBLE,
                  source, rank+2, MPI_COMM_WORLD, &status);
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(3*nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  MPI_Barrier(MPI_COMM_WORLD);

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
        marker and initialize them to zero.     ---*/
  nElem_Bound = new unsigned long[nMarker];
  for(unsigned short i=0; i<nMarker; ++i)
    nElem_Bound[i] = 0;

  /*--- Determine the global element ID's of the elements stored on this rank.
        Sort them in increasing order, such that an easy search can be done.
        In the same loop determine the upper bound for the local nodes (without
        halos) and the number of boundary elements for every marker. ---*/
  local_elem = local_node = 0;
  for(int i=0; i<nRankRecv; ++i) local_elem += longRecvBuf[i][0];

  vector<unsigned long> globalElemID;
  globalElemID.reserve(local_elem);

  for(int i=0; i<nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for(unsigned long j=0; j<longRecvBuf[i][0]; ++j) {
      globalElemID.push_back(longRecvBuf[i][indL]);

      unsigned short nDOFsGrid = shortRecvBuf[i][indS+3];
      unsigned short nFaces    = shortRecvBuf[i][indS+5];
      indS += 2*nFaces + 7;
      indL += nDOFsGrid + nFaces + 2;
    }

    long nNodesThisRank = longRecvBuf[i][indL];
    local_node += nNodesThisRank;
    indL       += nNodesThisRank+1;

    for(unsigned iMarker=0; iMarker<nMarker; ++iMarker) {
      long nBoundElemThisRank = longRecvBuf[i][indL]; ++indL;
      nElem_Bound[iMarker] += nBoundElemThisRank;

      for(long j=0; j<nBoundElemThisRank; ++j) {
        short nDOFsBoundElem = shortRecvBuf[i][indS+2];
        indS += 3;
        indL += nDOFsBoundElem + 2;
      }
    }
  }

  sort(globalElemID.begin(), globalElemID.end());

  /*--- Determine the global element ID's of the halo elements. A vector of
        unsignedLong2T is used for this purpose, such that a possible periodic
        transformation can be taken into account. Neighbors with a periodic
        transformation will always become a halo element, even if the element
        is stored on this rank.             ---*/
  vector<unsignedLong2T> haloElements;

  for(int i=0; i<nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for(unsigned long j=0; j<longRecvBuf[i][0]; ++j) {
      unsigned short nDOFsGrid = shortRecvBuf[i][indS+3];
      unsigned short nFaces    = shortRecvBuf[i][indS+5];

      indS += 7;
      indL += nDOFsGrid + 2;
      for(unsigned short k=0; k<nFaces; ++k, indS+=2, ++indL) {
        if(longRecvBuf[i][indL] != -1) {
          bool neighborIsInternal = false;
          if(shortRecvBuf[i][indS] == -1) {
            neighborIsInternal = binary_search(globalElemID.begin(),
                                               globalElemID.end(),
                                               longRecvBuf[i][indL]);
          }

          if( !neighborIsInternal )
            haloElements.push_back(unsignedLong2T(longRecvBuf[i][indL],
                                                  shortRecvBuf[i][indS]+1));
        }
      }
    }
  }

  sort(haloElements.begin(), haloElements.end());
  vector<unsignedLong2T>::iterator lastHalo = unique(haloElements.begin(), haloElements.end());
  haloElements.erase(lastHalo, haloElements.end());

  /*----------------------------------------------------------------------------*/
  /*--- Step 2: Store the elements, nodes and boundary elements in the data  ---*/
  /*---         structures used by the FEM solver.                           ---*/
  /*----------------------------------------------------------------------------*/

  /* Determine the mapping from the global element number to the local entry.
     At the moment the sequence is based on the global element ID, but one can
     change this in order to have a better load balance for the threads. */
  map<unsigned long,  unsigned long> mapGlobalElemIDToInd;
  map<unsignedLong2T, unsigned long> mapGlobalHaloElemToInd;

  nVolElemOwned = globalElemID.size();
  nVolElemTot   = nVolElemOwned + haloElements.size();

  for(unsigned long i=0; i<nVolElemOwned; ++i)
    mapGlobalElemIDToInd[globalElemID[i]] = i;

  for(unsigned long i=0; i<haloElements.size(); ++i)
    mapGlobalHaloElemToInd[haloElements[i]] = nVolElemOwned + i;

  /*--- Allocate the memory for the volume elements, the nodes
        and the surface elements of the boundaries.    ---*/
  volElem.resize(nVolElemTot);
  meshPoints.reserve(local_node);

  boundaries.resize(nMarker);
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    boundaries[iMarker].markerTag = config->GetMarker_All_TagBound(iMarker);
    boundaries[iMarker].surfElem.reserve(nElem_Bound[iMarker]);
  }

  /*--- Copy the data from the communication buffers. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /* The data for the volume elements. Loop over these elements in the buffer. */
    unsigned long indL = 1, indS = 0, indD = 0;
    for(unsigned long j=0; j<longRecvBuf[i][0]; ++j) {

      /* Determine the location in volElem where this data must be stored. */
      unsigned long elemID = longRecvBuf[i][indL++];
      map<unsigned long, unsigned long>::iterator MI = mapGlobalElemIDToInd.find(elemID);
      unsigned long ind = MI->second;

      /* Store the data. */
      volElem[ind].elemIsOwned = true;
      volElem[ind].elemIsPeriodicHalo = false;

      volElem[ind].VTK_Type  = shortRecvBuf[i][indS++];
      volElem[ind].nPolyGrid = shortRecvBuf[i][indS++];
      volElem[ind].nPolySol  = shortRecvBuf[i][indS++];
      volElem[ind].nDOFsGrid = shortRecvBuf[i][indS++];
      volElem[ind].nDOFsSol  = shortRecvBuf[i][indS++];
      volElem[ind].nFaces    = shortRecvBuf[i][indS++];

      volElem[ind].JacIsConsideredConstant = (bool) shortRecvBuf[i][indS++];

      volElem[ind].elemIDGlobal        = elemID;
      volElem[ind].offsetDOFsSolGlobal = longRecvBuf[i][indL++];

      volElem[ind].nodeIDsGrid.resize(volElem[ind].nDOFsGrid);

      volElem[ind].neighborElemIDMatchingFaces.resize(volElem[ind].nFaces);
      volElem[ind].periodIndexNeighbors.resize(volElem[ind].nFaces);
      volElem[ind].JacFacesIsConsideredConstant.resize(volElem[ind].nFaces);

      for(unsigned short k=0; k<volElem[ind].nDOFsGrid; ++k)
        volElem[ind].nodeIDsGrid[k] = longRecvBuf[i][indL++];

      for(unsigned short k=0; k<volElem[ind].nFaces; ++k)
        volElem[ind].neighborElemIDMatchingFaces[k] = longRecvBuf[i][indL++];

      for(unsigned short k=0; k<volElem[ind].nFaces; ++k) {
        volElem[ind].periodIndexNeighbors[k] = shortRecvBuf[i][indS++];
        volElem[ind].JacFacesIsConsideredConstant[k] = (bool) shortRecvBuf[i][indS++];
      }
    }

    /* The data for the nodes. Loop over these nodes in the buffer and store
       them in meshPoints.             */
    unsigned long nNodesThisRank = longRecvBuf[i][indL++];
    for(unsigned long j=0; j<nNodesThisRank; ++j) {
      CPointFEM thisPoint;
      thisPoint.globalID = longRecvBuf[i][indL++];
      for(unsigned short k=0; k<nDim; ++k)
        thisPoint.coor[k] = doubleRecvBuf[i][indD++];

      meshPoints.push_back(thisPoint);
    }

    /* The data for the boundary markers. Loop over them. */
    for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

      unsigned long nElemThisRank = longRecvBuf[i][indL++];
      for(unsigned long j=0; j<nElemThisRank; ++j) {
        CSurfaceElementFEM thisSurfElem;

        thisSurfElem.VTK_Type  = shortRecvBuf[i][indS++];
        thisSurfElem.nPolyGrid = shortRecvBuf[i][indS++];
        thisSurfElem.nDOFsGrid = shortRecvBuf[i][indS++];

        thisSurfElem.volElemID         = longRecvBuf[i][indL++];
        thisSurfElem.boundElemIDGlobal = longRecvBuf[i][indL++];

        thisSurfElem.nodeIDsGrid.resize(thisSurfElem.nDOFsGrid);
        for(unsigned short k=0; k<thisSurfElem.nDOFsGrid; ++k)
          thisSurfElem.nodeIDsGrid[k] = longRecvBuf[i][indL++];

        boundaries[iMarker].surfElem.push_back(thisSurfElem);
      }
    }
  }

  /* Sort meshPoints in increasing order and remove the double entities. */
  sort(meshPoints.begin(), meshPoints.end());
  vector<CPointFEM>::iterator lastPoint = unique(meshPoints.begin(), meshPoints.end());
  meshPoints.erase(lastPoint, meshPoints.end());

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
  /*--- Step 3: Build the layer of halo elements.                            ---*/
  /*----------------------------------------------------------------------------*/

  MPI_Barrier(MPI_COMM_WORLD);
  cout << "CMeshFEM::CMeshFEM: Not implemented yet." << endl << flush;
  MPI_Barrier(MPI_COMM_WORLD);

#ifndef HAVE_MPI
  exit(EXIT_FAILURE);
#else
  MPI_Abort(MPI_COMM_WORLD,1);
  MPI_Finalize();
#endif

}

CMeshFEM_DG::CMeshFEM_DG(CGeometry *geometry, CConfig *config) : CMeshFEM(geometry, config) {

}
