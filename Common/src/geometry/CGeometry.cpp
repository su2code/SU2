/*!
 * \file CGeometry.cpp
 * \brief Implementation of the base geometry class.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/geometry/CGeometry.hpp"
#include "../../include/geometry/elements/CElement.hpp"

/*--- Cross product ---*/

#define CROSS(dest,v1,v2) \
(dest)[0] = (v1)[1]*(v2)[2] - (v1)[2]*(v2)[1];	\
(dest)[1] = (v1)[2]*(v2)[0] - (v1)[0]*(v2)[2];	\
(dest)[2] = (v1)[0]*(v2)[1] - (v1)[1]*(v2)[0];

/*--- Cross product ---*/

#define DOT(v1,v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1] + (v1)[2]*(v2)[2]);

/*--- a = b - c ---*/

#define SUB(dest,v1,v2) \
(dest)[0] = (v1)[0] - (v2)[0];	\
(dest)[1] = (v1)[1] - (v2)[1];	\
(dest)[2] = (v1)[2] - (v2)[2];

CGeometry::CGeometry(void) {

  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  nEdge      = 0;
  nPoint     = 0;
  nPointNode = 0;
  nElem      = 0;

  nElem_Bound         = NULL;
  Tag_to_Marker       = NULL;
  elem                = NULL;
  face                = NULL;
  bound               = NULL;
  node                = NULL;
  edge                = NULL;
  vertex              = NULL;
  nVertex             = NULL;
  newBound            = NULL;
  nNewElem_Bound      = NULL;
  Marker_All_SendRecv = NULL;

  XCoordList.clear();
  Xcoord_plane.clear();
  Ycoord_plane.clear();
  Zcoord_plane.clear();
  FaceArea_plane.clear();
  Plane_points.clear();

  /*--- Arrays for defining the linear partitioning ---*/

  beg_node = NULL;
  end_node = NULL;

  nPointLinear     = NULL;
  nPointCumulative = NULL;

  /*--- Containers for customized boundary conditions ---*/

  CustomBoundaryHeatFlux = NULL;      //Customized heat flux wall
  CustomBoundaryTemperature = NULL;   //Customized temperature wall

  /*--- MPI point-to-point data structures ---*/

  nP2PSend = 0;
  nP2PRecv = 0;

  countPerPoint = 0;

  bufD_P2PSend = NULL;
  bufD_P2PRecv = NULL;

  bufS_P2PSend = NULL;
  bufS_P2PRecv = NULL;

  req_P2PSend = NULL;
  req_P2PRecv = NULL;

  nPoint_P2PSend = NULL;
  nPoint_P2PRecv = NULL;

  Neighbors_P2PSend = NULL;
  Neighbors_P2PRecv = NULL;

  Local_Point_P2PSend = NULL;
  Local_Point_P2PRecv = NULL;

  /*--- MPI periodic data structures ---*/

  nPeriodicSend = 0;
  nPeriodicRecv = 0;

  countPerPeriodicPoint = 0;

  bufD_PeriodicSend = NULL;
  bufD_PeriodicRecv = NULL;

  bufS_PeriodicSend = NULL;
  bufS_PeriodicRecv = NULL;

  req_PeriodicSend = NULL;
  req_PeriodicRecv = NULL;

  nPoint_PeriodicSend = NULL;
  nPoint_PeriodicRecv = NULL;

  Neighbors_PeriodicSend = NULL;
  Neighbors_PeriodicRecv = NULL;

  Local_Point_PeriodicSend = NULL;
  Local_Point_PeriodicRecv = NULL;

  Local_Marker_PeriodicSend = NULL;
  Local_Marker_PeriodicRecv = NULL;

}

CGeometry::~CGeometry(void) {

  unsigned long iElem, iElem_Bound, iEdge, iFace, iPoint, iVertex;
  unsigned short iMarker;

  if (elem != NULL) {
    for (iElem = 0; iElem < nElem; iElem++)
      if (elem[iElem] != NULL) delete elem[iElem];
    delete[] elem;
  }

  if (bound != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (bound[iMarker][iElem_Bound] != NULL) delete bound[iMarker][iElem_Bound];
      }
      if (bound[iMarker] != NULL) delete [] bound[iMarker];
    }
    delete [] bound;
  }

  if (face != NULL) {
    for (iFace = 0; iFace < nFace; iFace ++)
      if (face[iFace] != NULL) delete face[iFace];
    delete[] face;
  }

  if (node != NULL) {
    for (iPoint = 0; iPoint < nPointNode; iPoint ++)
      if (node[iPoint] != NULL) delete node[iPoint];
    delete[] node;
  }


  if (edge != NULL) {
    for (iEdge = 0; iEdge < nEdge; iEdge ++)
      if (edge[iEdge] != NULL) delete edge[iEdge];
    delete[] edge;
  }

  if (vertex != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        if (vertex[iMarker][iVertex] != NULL) delete vertex[iMarker][iVertex];
      }
      if (vertex[iMarker] != NULL) delete [] vertex[iMarker];
    }
    delete [] vertex;
  }

  if (newBound != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (newBound[iMarker][iElem_Bound] != NULL) delete [] newBound[iMarker][iElem_Bound];
      }
      delete[] newBound[iMarker];
    }
    delete[] newBound;
  }

  if (nElem_Bound         != NULL) delete [] nElem_Bound;
  if (nVertex             != NULL) delete [] nVertex;
  if (nNewElem_Bound      != NULL) delete [] nNewElem_Bound;
  if (Marker_All_SendRecv != NULL) delete [] Marker_All_SendRecv;
  if (Tag_to_Marker       != NULL) delete [] Tag_to_Marker;

  if (beg_node != NULL) delete [] beg_node;
  if (end_node != NULL) delete [] end_node;
  if (nPointLinear      != NULL) delete [] nPointLinear;
  if (nPointCumulative  != NULL) delete [] nPointCumulative;

  if(CustomBoundaryHeatFlux != NULL){
    for(iMarker=0; iMarker < nMarker; iMarker++){
      if (CustomBoundaryHeatFlux[iMarker] != NULL) delete [] CustomBoundaryHeatFlux[iMarker];
    }
    delete [] CustomBoundaryHeatFlux;
  }

  if(CustomBoundaryTemperature != NULL){
    for(iMarker=0; iMarker < nMarker; iMarker++){
      if(CustomBoundaryTemperature[iMarker] != NULL) delete [] CustomBoundaryTemperature[iMarker];
    }
    delete [] CustomBoundaryTemperature;
  }

  /*--- Delete structures for MPI point-to-point communication. ---*/

  if (bufD_P2PRecv != NULL) delete [] bufD_P2PRecv;
  if (bufD_P2PSend != NULL) delete [] bufD_P2PSend;

  if (bufS_P2PRecv != NULL) delete [] bufS_P2PRecv;
  if (bufS_P2PSend != NULL) delete [] bufS_P2PSend;

  if (req_P2PSend != NULL) delete [] req_P2PSend;
  if (req_P2PRecv != NULL) delete [] req_P2PRecv;

  if (nPoint_P2PRecv != NULL) delete [] nPoint_P2PRecv;
  if (nPoint_P2PSend != NULL) delete [] nPoint_P2PSend;

  if (Neighbors_P2PSend != NULL) delete [] Neighbors_P2PSend;
  if (Neighbors_P2PRecv != NULL) delete [] Neighbors_P2PRecv;

  if (Local_Point_P2PSend != NULL) delete [] Local_Point_P2PSend;
  if (Local_Point_P2PRecv != NULL) delete [] Local_Point_P2PRecv;

  /*--- Delete structures for MPI periodic communication. ---*/

  if (bufD_PeriodicRecv != NULL) delete [] bufD_PeriodicRecv;
  if (bufD_PeriodicSend != NULL) delete [] bufD_PeriodicSend;

  if (bufS_PeriodicRecv != NULL) delete [] bufS_PeriodicRecv;
  if (bufS_PeriodicSend != NULL) delete [] bufS_PeriodicSend;

  if (req_PeriodicSend != NULL) delete [] req_PeriodicSend;
  if (req_PeriodicRecv != NULL) delete [] req_PeriodicRecv;

  if (nPoint_PeriodicRecv != NULL) delete [] nPoint_PeriodicRecv;
  if (nPoint_PeriodicSend != NULL) delete [] nPoint_PeriodicSend;

  if (Neighbors_PeriodicSend != NULL) delete [] Neighbors_PeriodicSend;
  if (Neighbors_PeriodicRecv != NULL) delete [] Neighbors_PeriodicRecv;

  if (Local_Point_PeriodicSend != NULL) delete [] Local_Point_PeriodicSend;
  if (Local_Point_PeriodicRecv != NULL) delete [] Local_Point_PeriodicRecv;

  if (Local_Marker_PeriodicSend != NULL) delete [] Local_Marker_PeriodicSend;
  if (Local_Marker_PeriodicRecv != NULL) delete [] Local_Marker_PeriodicRecv;

}

void CGeometry::PreprocessP2PComms(CGeometry *geometry,
                                   CConfig *config) {

  /*--- We start with the send and receive lists already available in
   the form of SEND_RECEIVE boundary markers. We will loop through
   these markers and establish the neighboring ranks and number of
   send/recv points per pair. We will store this information and set
   up persistent data structures so that we can reuse them throughout
   the calculation for any point-to-point communications. The goal
   is to break the non-blocking comms into InitiateComms() and
   CompleteComms() in separate routines so that we can overlap the
   communication and computation to hide the communication latency. ---*/

  /*--- Local variables. ---*/

  unsigned short iMarker;
  unsigned long  nVertexS, nVertexR, iVertex, MarkerS, MarkerR;

  int iRank, iSend, iRecv, count;

  /*--- Create some temporary structures for tracking sends/recvs. ---*/

  int *nPoint_Send_All = new int[size+1]; nPoint_Send_All[0] = 0;
  int *nPoint_Recv_All = new int[size+1]; nPoint_Recv_All[0] = 0;
  int *nPoint_Flag = new int[size];

  for (iRank = 0; iRank < size; iRank++) {
    nPoint_Send_All[iRank] = 0; nPoint_Recv_All[iRank] = 0; nPoint_Flag[iRank]= -1;
  }
  nPoint_Send_All[size] = 0; nPoint_Recv_All[size] = 0;

  /*--- Loop through all of our SEND_RECEIVE markers and track
   our sends with each rank. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      /*--- Get the destination rank and number of points to send. ---*/

      iRank    = config->GetMarker_All_SendRecv(iMarker)-1;
      nVertexS = geometry->nVertex[iMarker];

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nPoint_Flag[iRank] != (int)iMarker)) {
        nPoint_Flag[iRank]        = (int)iMarker;
        nPoint_Send_All[iRank+1] += nVertexS;
      }

    }
  }

  delete [] nPoint_Flag;

  /*--- Communicate the number of points to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nPoint_Send_All[1]), 1, MPI_INT,
                    &(nPoint_Recv_All[1]), 1, MPI_INT, MPI_COMM_WORLD);

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  nP2PSend = 0; nP2PRecv = 0;

  for (iRank = 0; iRank < size; iRank++) {
    if ((iRank != rank) && (nPoint_Send_All[iRank+1] > 0)) nP2PSend++;
    if ((iRank != rank) && (nPoint_Recv_All[iRank+1] > 0)) nP2PRecv++;

    nPoint_Send_All[iRank+1] += nPoint_Send_All[iRank];
    nPoint_Recv_All[iRank+1] += nPoint_Recv_All[iRank];
  }

  /*--- Allocate only as much memory as we need for the P2P neighbors. ---*/

  nPoint_P2PSend = new int[nP2PSend+1]; nPoint_P2PSend[0] = 0;
  nPoint_P2PRecv = new int[nP2PRecv+1]; nPoint_P2PRecv[0] = 0;

  Neighbors_P2PSend = new int[nP2PSend];
  Neighbors_P2PRecv = new int[nP2PRecv];

  iSend = 0; iRecv = 0;
  for (iRank = 0; iRank < size; iRank++) {

    if ((nPoint_Send_All[iRank+1] > nPoint_Send_All[iRank]) && (iRank != rank)) {
      Neighbors_P2PSend[iSend] = iRank;
      nPoint_P2PSend[iSend+1] = nPoint_Send_All[iRank+1];
      iSend++;
    }

    if ((nPoint_Recv_All[iRank+1] > nPoint_Recv_All[iRank]) && (iRank != rank)) {
      Neighbors_P2PRecv[iRecv] = iRank;
      nPoint_P2PRecv[iRecv+1] = nPoint_Recv_All[iRank+1];
      iRecv++;
    }

  }

  /*--- Create a reverse mapping of the message to the rank so that we
   can quickly access the correct data in the buffers when receiving
   messages dynamically. ---*/

  P2PSend2Neighbor.clear();
  for (iSend = 0; iSend < nP2PSend; iSend++)
    P2PSend2Neighbor[Neighbors_P2PSend[iSend]] = iSend;

  P2PRecv2Neighbor.clear();
  for (iRecv = 0; iRecv < nP2PRecv; iRecv++)
    P2PRecv2Neighbor[Neighbors_P2PRecv[iRecv]] = iRecv;

  delete [] nPoint_Send_All;
  delete [] nPoint_Recv_All;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  Local_Point_P2PSend = NULL;
  Local_Point_P2PSend = new unsigned long[nPoint_P2PSend[nP2PSend]];
  for (iSend = 0; iSend < nPoint_P2PSend[nP2PSend]; iSend++)
    Local_Point_P2PSend[iSend] = 0;

  Local_Point_P2PRecv = NULL;
  Local_Point_P2PRecv = new unsigned long[nPoint_P2PRecv[nP2PRecv]];
  for (iRecv = 0; iRecv < nPoint_P2PRecv[nP2PRecv]; iRecv++)
    Local_Point_P2PRecv[iRecv] = 0;

  /*--- We allocate the memory for communicating values in a later step
   once we know the maximum packet size that we need to communicate. This
   memory is deallocated and reallocated automatically in the case that
   the previously allocated memory is not sufficient. ---*/

  bufD_P2PSend = NULL;
  bufD_P2PRecv = NULL;

  bufS_P2PSend = NULL;
  bufS_P2PRecv = NULL;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nP2PSend > 0) {
    req_P2PSend   = new SU2_MPI::Request[nP2PSend];
  }
  if (nP2PRecv > 0) {
    req_P2PRecv   = new SU2_MPI::Request[nP2PRecv];
  }

  /*--- Build lists of local index values for send. ---*/

  count = 0;
  for (iSend = 0; iSend < nP2PSend; iSend++) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
          (config->GetMarker_All_SendRecv(iMarker) > 0)) {

        MarkerS  = iMarker;
        nVertexS = geometry->nVertex[MarkerS];
        iRank    = config->GetMarker_All_SendRecv(MarkerS)-1;

        if (iRank == Neighbors_P2PSend[iSend]) {
          for (iVertex = 0; iVertex < nVertexS; iVertex++) {
            Local_Point_P2PSend[count] = geometry->vertex[MarkerS][iVertex]->GetNode();
            count++;
          }
        }

      }
    }
  }

  /*--- Build lists of local index values for receive. ---*/

  count = 0;
  for (iRecv = 0; iRecv < nP2PRecv; iRecv++) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
          (config->GetMarker_All_SendRecv(iMarker) > 0)) {

        MarkerR  = iMarker+1;
        nVertexR = geometry->nVertex[MarkerR];
        iRank    = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

        if (iRank == Neighbors_P2PRecv[iRecv]) {
          for (iVertex = 0; iVertex < nVertexR; iVertex++) {
            Local_Point_P2PRecv[count] = geometry->vertex[MarkerR][iVertex]->GetNode();
            count++;
          }
        }

      }
    }
  }

  /*--- In the future, some additional data structures could be created
   here to separate the interior and boundary nodes in order to help
   further overlap computation and communication. ---*/

}

void CGeometry::AllocateP2PComms(unsigned short val_countPerPoint) {

  /*--- This routine is activated whenever we attempt to perform
   a point-to-point MPI communication with our neighbors but the
   memory buffer allocated is not large enough for the packet size.
   Therefore, we deallocate the previously allocated space and
   reallocate a large enough array. Note that after the first set
   communications, this routine will not need to be called again. ---*/

  int iSend, iRecv;

  /*--- Store the larger packet size to the class data. ---*/

  countPerPoint = val_countPerPoint;

  /*-- Deallocate and reallocate our su2double cummunication memory. ---*/

  if (bufD_P2PSend != NULL) delete [] bufD_P2PSend;

  bufD_P2PSend = new su2double[countPerPoint*nPoint_P2PSend[nP2PSend]];
  for (iSend = 0; iSend < countPerPoint*nPoint_P2PSend[nP2PSend]; iSend++)
    bufD_P2PSend[iSend] = 0.0;

  if (bufD_P2PRecv != NULL) delete [] bufD_P2PRecv;

  bufD_P2PRecv = new su2double[countPerPoint*nPoint_P2PRecv[nP2PRecv]];
  for (iRecv = 0; iRecv < countPerPoint*nPoint_P2PRecv[nP2PRecv]; iRecv++)
    bufD_P2PRecv[iRecv] = 0.0;

  if (bufS_P2PSend != NULL) delete [] bufS_P2PSend;

  bufS_P2PSend = new unsigned short[countPerPoint*nPoint_P2PSend[nP2PSend]];
  for (iSend = 0; iSend < countPerPoint*nPoint_P2PSend[nP2PSend]; iSend++)
    bufS_P2PSend[iSend] = 0;

  if (bufS_P2PRecv != NULL) delete [] bufS_P2PRecv;

  bufS_P2PRecv = new unsigned short[countPerPoint*nPoint_P2PRecv[nP2PRecv]];
  for (iRecv = 0; iRecv < countPerPoint*nPoint_P2PRecv[nP2PRecv]; iRecv++)
    bufS_P2PRecv[iRecv] = 0;

}

void CGeometry::PostP2PRecvs(CGeometry *geometry,
                             CConfig *config,
                             unsigned short commType,
                             bool val_reverse) {

  /*--- Local variables ---*/

  int iMessage, iRecv, offset, nPointP2P, count, source, tag;

  /*--- Launch the non-blocking recv's first. Note that we have stored
   the counts and sources, so we can launch these before we even load
   the data and send from the neighbor ranks. ---*/

  iMessage = 0;
  for (iRecv = 0; iRecv < nP2PRecv; iRecv++) {

    /*--- In some instances related to the adjoint solver, we need
     to reverse the direction of communications such that the normal
     send nodes become the recv nodes and vice-versa. ---*/

    if (val_reverse) {

      /*--- Compute our location in the buffer using the send data
       structure since we are reversing the comms. ---*/

      offset = countPerPoint*nPoint_P2PSend[iRecv];

      /*--- Take advantage of cumulative storage format to get the number
       of elems that we need to recv. Note again that we select the send
       points here as the recv points. ---*/

      nPointP2P = nPoint_P2PSend[iRecv+1] - nPoint_P2PSend[iRecv];

      /*--- Total count can include multiple pieces of data per element. ---*/

      count = countPerPoint*nPointP2P;

      /*--- Get the rank from which we receive the message. Note again
       that we use the send rank as the source instead of the recv rank. ---*/

      source = Neighbors_P2PSend[iRecv];
      tag    = source + 1;

      /*--- Post non-blocking recv for this proc. Note that we use the
       send buffer here too. This is important to make sure the arrays
       are the correct size. ---*/

      switch (commType) {
        case COMM_TYPE_DOUBLE:
          SU2_MPI::Irecv(&(bufD_P2PSend[offset]), count, MPI_DOUBLE,
                         source, tag, MPI_COMM_WORLD, &(req_P2PRecv[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_SHORT:
          SU2_MPI::Irecv(&(bufS_P2PSend[offset]), count, MPI_UNSIGNED_SHORT,
                         source, tag, MPI_COMM_WORLD, &(req_P2PRecv[iMessage]));
          break;
        default:
          SU2_MPI::Error("Unrecognized data type for point-to-point MPI comms.",
                         CURRENT_FUNCTION);
          break;
      }

    } else {

      /*--- Compute our location in the recv buffer. ---*/

      offset = countPerPoint*nPoint_P2PRecv[iRecv];

      /*--- Take advantage of cumulative storage format to get the number
       of elems that we need to recv. ---*/

      nPointP2P = nPoint_P2PRecv[iRecv+1] - nPoint_P2PRecv[iRecv];

      /*--- Total count can include multiple pieces of data per element. ---*/

      count = countPerPoint*nPointP2P;

      /*--- Get the rank from which we receive the message. ---*/

      source = Neighbors_P2PRecv[iRecv];
      tag    = source + 1;

      /*--- Post non-blocking recv for this proc. ---*/

      switch (commType) {
        case COMM_TYPE_DOUBLE:
          SU2_MPI::Irecv(&(bufD_P2PRecv[offset]), count, MPI_DOUBLE,
                         source, tag, MPI_COMM_WORLD, &(req_P2PRecv[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_SHORT:
          SU2_MPI::Irecv(&(bufS_P2PRecv[offset]), count, MPI_UNSIGNED_SHORT,
                         source, tag, MPI_COMM_WORLD, &(req_P2PRecv[iMessage]));
          break;
        default:
          SU2_MPI::Error("Unrecognized data type for point-to-point MPI comms.",
                         CURRENT_FUNCTION);
          break;
      }

    }

    /*--- Increment message counter. ---*/

    iMessage++;

  }

}

void CGeometry::PostP2PSends(CGeometry *geometry,
                             CConfig *config,
                             unsigned short commType,
                             int val_iSend,
                             bool val_reverse) {

  /*--- Local variables ---*/

  int iMessage, offset, nPointP2P, count, dest, tag;

  /*--- Post the non-blocking send as soon as the buffer is loaded. ---*/

  iMessage = val_iSend;

  /*--- In some instances related to the adjoint solver, we need
   to reverse the direction of communications such that the normal
   send nodes become the recv nodes and vice-versa. ---*/

  if (val_reverse) {

    /*--- Compute our location in the buffer using the recv data
     structure since we are reversing the comms. ---*/

    offset = countPerPoint*nPoint_P2PRecv[val_iSend];

    /*--- Take advantage of cumulative storage format to get the number
     of points that we need to send. Note again that we select the recv
     points here as the send points. ---*/

    nPointP2P = nPoint_P2PRecv[val_iSend+1] - nPoint_P2PRecv[val_iSend];

    /*--- Total count can include multiple pieces of data per element. ---*/

    count = countPerPoint*nPointP2P;

    /*--- Get the rank to which we send the message. Note again
     that we use the recv rank as the dest instead of the send rank. ---*/

    dest = Neighbors_P2PRecv[val_iSend];
    tag  = rank + 1;

    /*--- Post non-blocking send for this proc. Note that we use the
     send buffer here too. This is important to make sure the arrays
     are the correct size. ---*/

    switch (commType) {
      case COMM_TYPE_DOUBLE:
        SU2_MPI::Isend(&(bufD_P2PRecv[offset]), count, MPI_DOUBLE,
                       dest, tag, MPI_COMM_WORLD, &(req_P2PSend[iMessage]));
        break;
      case COMM_TYPE_UNSIGNED_SHORT:
        SU2_MPI::Isend(&(bufS_P2PRecv[offset]), count, MPI_UNSIGNED_SHORT,
                       dest, tag, MPI_COMM_WORLD, &(req_P2PSend[iMessage]));
        break;
      default:
        SU2_MPI::Error("Unrecognized data type for point-to-point MPI comms.",
                       CURRENT_FUNCTION);
        break;
    }

  } else {

    /*--- Compute our location in the send buffer. ---*/

    offset = countPerPoint*nPoint_P2PSend[val_iSend];

    /*--- Take advantage of cumulative storage format to get the number
     of points that we need to send. ---*/

    nPointP2P = nPoint_P2PSend[val_iSend+1] - nPoint_P2PSend[val_iSend];

    /*--- Total count can include multiple pieces of data per element. ---*/

    count = countPerPoint*nPointP2P;

    /*--- Get the rank to which we send the message. ---*/

    dest = Neighbors_P2PSend[val_iSend];
    tag  = rank + 1;

    /*--- Post non-blocking send for this proc. ---*/

    switch (commType) {
      case COMM_TYPE_DOUBLE:
        SU2_MPI::Isend(&(bufD_P2PSend[offset]), count, MPI_DOUBLE,
                       dest, tag, MPI_COMM_WORLD, &(req_P2PSend[iMessage]));
        break;
      case COMM_TYPE_UNSIGNED_SHORT:
        SU2_MPI::Isend(&(bufS_P2PSend[offset]), count, MPI_UNSIGNED_SHORT,
                       dest, tag, MPI_COMM_WORLD, &(req_P2PSend[iMessage]));
        break;
      default:
        SU2_MPI::Error("Unrecognized data type for point-to-point MPI comms.",
                       CURRENT_FUNCTION);
        break;
    }

  }

}

void CGeometry::InitiateComms(CGeometry *geometry,
                              CConfig *config,
                              unsigned short commType) {

  /*--- Local variables ---*/

  unsigned short iDim;
  unsigned short COUNT_PER_POINT = 0;
  unsigned short MPI_TYPE        = 0;

  unsigned long iPoint, msg_offset, buf_offset;

  int iMessage, iSend, nSend;

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  switch (commType) {
    case COORDINATES:
      COUNT_PER_POINT  = nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case GRID_VELOCITY:
      COUNT_PER_POINT  = nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case COORDINATES_OLD:
      if (config->GetTime_Marching() == DT_STEPPING_2ND)
        COUNT_PER_POINT  = nDim*2;
      else
        COUNT_PER_POINT  = nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case MAX_LENGTH:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case NEIGHBORS:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_UNSIGNED_SHORT;
      break;
    default:
      SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                     CURRENT_FUNCTION);
      break;
  }

  /*--- Check to make sure we have created a large enough buffer
   for these comms during preprocessing. This is only for the su2double
   buffer. It will be reallocated whenever we find a larger count
   per point. After the first cycle of comms, this should be inactive. ---*/

  if (COUNT_PER_POINT > geometry->countPerPoint) {
    geometry->AllocateP2PComms(COUNT_PER_POINT);
  }

  /*--- Set some local pointers to make access simpler. ---*/

  su2double *bufDSend      = geometry->bufD_P2PSend;
  unsigned short *bufSSend = geometry->bufS_P2PSend;

  su2double *vector = NULL;

  /*--- Load the specified quantity from the solver into the generic
   communication buffer in the geometry class. ---*/

  if (nP2PSend > 0) {

    /*--- Post all non-blocking recvs first before sends. ---*/

    geometry->PostP2PRecvs(geometry, config, MPI_TYPE, false);

    for (iMessage = 0; iMessage < nP2PSend; iMessage++) {

      /*--- Get the offset in the buffer for the start of this message. ---*/

      msg_offset = nPoint_P2PSend[iMessage];

      /*--- Total count can include multiple pieces of data per element. ---*/

      nSend = (nPoint_P2PSend[iMessage+1] - nPoint_P2PSend[iMessage]);

      for (iSend = 0; iSend < nSend; iSend++) {

        /*--- Get the local index for this communicated data. ---*/

        iPoint = geometry->Local_Point_P2PSend[msg_offset + iSend];

        /*--- Compute the offset in the recv buffer for this point. ---*/

        buf_offset = (msg_offset + iSend)*countPerPoint;

        switch (commType) {
          case COORDINATES:
            vector = node[iPoint]->GetCoord();
            for (iDim = 0; iDim < nDim; iDim++)
              bufDSend[buf_offset+iDim] = vector[iDim];
            break;
          case GRID_VELOCITY:
            vector = node[iPoint]->GetGridVel();
            for (iDim = 0; iDim < nDim; iDim++)
              bufDSend[buf_offset+iDim] = vector[iDim];
            break;
          case COORDINATES_OLD:
            vector = node[iPoint]->GetCoord_n();
            for (iDim = 0; iDim < nDim; iDim++) {
              bufDSend[buf_offset+iDim] = vector[iDim];
            }
            if (config->GetTime_Marching() == DT_STEPPING_2ND) {
              vector = node[iPoint]->GetCoord_n1();
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset+nDim+iDim] = vector[iDim];
              }
            }
            break;
          case MAX_LENGTH:
            bufDSend[buf_offset] = node[iPoint]->GetMaxLength();
            break;
          case NEIGHBORS:
            bufSSend[buf_offset] = geometry->node[iPoint]->GetnNeighbor();
            break;
          default:
            SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                           CURRENT_FUNCTION);
            break;
        }
      }

      /*--- Launch the point-to-point MPI send for this message. ---*/

      geometry->PostP2PSends(geometry, config, MPI_TYPE, iMessage, false);

    }
  }

}

void CGeometry::CompleteComms(CGeometry *geometry,
                              CConfig *config,
                              unsigned short commType) {

  /*--- Local variables ---*/

  unsigned short iDim;
  unsigned long iPoint, iRecv, nRecv, msg_offset, buf_offset;

  int ind, source, iMessage, jRecv;
  SU2_MPI::Status status;

  /*--- Set some local pointers to make access simpler. ---*/

  su2double *bufDRecv      = geometry->bufD_P2PRecv;
  unsigned short *bufSRecv = geometry->bufS_P2PRecv;

  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. Note that we
   recv and store the data in any order to take advantage of the
   non-blocking comms. ---*/

  if (nP2PRecv > 0) {

    for (iMessage = 0; iMessage < nP2PRecv; iMessage++) {

      /*--- For efficiency, recv the messages dynamically based on
       the order they arrive. ---*/

      SU2_MPI::Waitany(nP2PRecv, req_P2PRecv, &ind, &status);

      /*--- Once we have recv'd a message, get the source rank. ---*/

      source = status.MPI_SOURCE;

      /*--- We know the offsets based on the source rank. ---*/

      jRecv = P2PRecv2Neighbor[source];

      /*--- Get the offset in the buffer for the start of this message. ---*/

      msg_offset = nPoint_P2PRecv[jRecv];

      /*--- Get the number of packets to be received in this message. ---*/

      nRecv = nPoint_P2PRecv[jRecv+1] - nPoint_P2PRecv[jRecv];

      for (iRecv = 0; iRecv < nRecv; iRecv++) {

        /*--- Get the local index for this communicated data. ---*/

        iPoint = geometry->Local_Point_P2PRecv[msg_offset + iRecv];

        /*--- Compute the total offset in the recv buffer for this point. ---*/

        buf_offset = (msg_offset + iRecv)*countPerPoint;

        /*--- Store the data correctly depending on the quantity. ---*/

        switch (commType) {
          case COORDINATES:
            for (iDim = 0; iDim < nDim; iDim++)
              node[iPoint]->SetCoord(iDim, bufDRecv[buf_offset+iDim]);
            break;
          case GRID_VELOCITY:
            for (iDim = 0; iDim < nDim; iDim++)
              node[iPoint]->SetGridVel(iDim, bufDRecv[buf_offset+iDim]);
            break;
          case COORDINATES_OLD:
            node[iPoint]->SetCoord_n(&bufDRecv[buf_offset]);
            if (config->GetTime_Marching() == DT_STEPPING_2ND)
              node[iPoint]->SetCoord_n1(&bufDRecv[buf_offset+nDim]);
            break;
          case MAX_LENGTH:
            node[iPoint]->SetMaxLength(bufDRecv[buf_offset]);
            break;
          case NEIGHBORS:
            node[iPoint]->SetnNeighbor(bufSRecv[buf_offset]);
            break;
          default:
            SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                           CURRENT_FUNCTION);
            break;
        }
      }
    }

    /*--- Verify that all non-blocking point-to-point sends have finished.
     Note that this should be satisfied, as we have received all of the
     data in the loop above at this point. ---*/

#ifdef HAVE_MPI
    SU2_MPI::Waitall(nP2PSend, req_P2PSend, MPI_STATUS_IGNORE);
#endif

  }

}

void CGeometry::PreprocessPeriodicComms(CGeometry *geometry,
                                        CConfig *config) {

  /*--- We start with the send and receive lists already available in
   the form of stored periodic point-donor pairs. We will loop through
   these markers and establish the neighboring ranks and number of
   send/recv points per pair. We will store this information and set
   up persistent data structures so that we can reuse them throughout
   the calculation for any periodic boundary communications. The goal
   is to break the non-blocking comms into InitiatePeriodicComms() and
   CompletePeriodicComms() in separate routines so that we can overlap the
   communication and computation to hide the communication latency. ---*/

  /*--- Local variables. ---*/

  unsigned short iMarker;
  unsigned long iPoint, iVertex, iPeriodic;

  int iRank, iSend, iRecv, ii, jj;

  /*--- Create some temporary structures for tracking sends/recvs. ---*/

  int *nPoint_Send_All = new int[size+1]; nPoint_Send_All[0] = 0;
  int *nPoint_Recv_All = new int[size+1]; nPoint_Recv_All[0] = 0;
  int *nPoint_Flag     = new int[size];

  for (iRank = 0; iRank < size; iRank++) {
    nPoint_Send_All[iRank] = 0;
    nPoint_Recv_All[iRank] = 0;
    nPoint_Flag[iRank]= -1;
  }
  nPoint_Send_All[size] = 0; nPoint_Recv_All[size] = 0;

  /*--- Loop through all of our periodic markers and track
   our sends with each rank. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Get the current periodic point index. We only communicate
         the owned nodes on a rank, as the MPI comms will take care of
         the halos after completing the periodic comms. ---*/

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          /*--- Get the rank that holds the matching periodic point
           on the other marker in the periodic pair. ---*/

          iRank = (int)geometry->vertex[iMarker][iVertex]->GetDonorProcessor();

          /*--- If we have not visited this point last, increment our
           number of points that must be sent to a particular proc. ---*/

          if ((nPoint_Flag[iRank] != (int)iPoint)) {
            nPoint_Flag[iRank]    = (int)iPoint;
            nPoint_Send_All[iRank+1] += 1;
          }

        }
      }
    }
  }

  delete [] nPoint_Flag;

  /*--- Communicate the number of points to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many periodic points it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nPoint_Send_All[1]), 1, MPI_INT,
                    &(nPoint_Recv_All[1]), 1, MPI_INT, MPI_COMM_WORLD);

  /*--- Check how many messages we will be sending and receiving.
   Here we also put the counters into cumulative storage format to
   make the communications simpler. Note that we are allowing each
   rank to communicate to themselves in these counters, although
   it will not be done through MPI. ---*/

  nPeriodicSend = 0; nPeriodicRecv = 0;

  for (iRank = 0; iRank < size; iRank++) {
    if ((nPoint_Send_All[iRank+1] > 0)) nPeriodicSend++;
    if ((nPoint_Recv_All[iRank+1] > 0)) nPeriodicRecv++;

    nPoint_Send_All[iRank+1] += nPoint_Send_All[iRank];
    nPoint_Recv_All[iRank+1] += nPoint_Recv_All[iRank];
  }

  /*--- Allocate only as much memory as needed for the periodic neighbors. ---*/

  nPoint_PeriodicSend = new int[nPeriodicSend+1]; nPoint_PeriodicSend[0] = 0;
  nPoint_PeriodicRecv = new int[nPeriodicRecv+1]; nPoint_PeriodicRecv[0] = 0;

  Neighbors_PeriodicSend = new int[nPeriodicSend];
  Neighbors_PeriodicRecv = new int[nPeriodicRecv];

  iSend = 0; iRecv = 0;
  for (iRank = 0; iRank < size; iRank++) {
    if ((nPoint_Send_All[iRank+1] > nPoint_Send_All[iRank])) {
      Neighbors_PeriodicSend[iSend] = iRank;
      nPoint_PeriodicSend[iSend+1] = nPoint_Send_All[iRank+1];
      iSend++;
    }
    if ((nPoint_Recv_All[iRank+1] > nPoint_Recv_All[iRank])) {
      Neighbors_PeriodicRecv[iRecv] = iRank;
      nPoint_PeriodicRecv[iRecv+1] = nPoint_Recv_All[iRank+1];
      iRecv++;
    }
  }

  /*--- Create a reverse mapping of the message to the rank so that we
   can quickly access the correct data in the buffers when receiving
   messages dynamically later during the iterations. ---*/

  PeriodicSend2Neighbor.clear();
  for (iSend = 0; iSend < nPeriodicSend; iSend++)
    PeriodicSend2Neighbor[Neighbors_PeriodicSend[iSend]] = iSend;

  PeriodicRecv2Neighbor.clear();
  for (iRecv = 0; iRecv < nPeriodicRecv; iRecv++)
    PeriodicRecv2Neighbor[Neighbors_PeriodicRecv[iRecv]] = iRecv;

  delete [] nPoint_Send_All;
  delete [] nPoint_Recv_All;

  /*--- Allocate the memory to store the local index values for both
   the send and receive periodic points and periodic index. ---*/

  Local_Point_PeriodicSend = NULL;
  Local_Point_PeriodicSend = new unsigned long[nPoint_PeriodicSend[nPeriodicSend]];
  for (iSend = 0; iSend < nPoint_PeriodicSend[nPeriodicSend]; iSend++)
    Local_Point_PeriodicSend[iSend] = 0;

  Local_Marker_PeriodicSend = NULL;
  Local_Marker_PeriodicSend = new unsigned long[nPoint_PeriodicSend[nPeriodicSend]];
  for (iSend = 0; iSend < nPoint_PeriodicSend[nPeriodicSend]; iSend++)
    Local_Marker_PeriodicSend[iSend] = 0;

  Local_Point_PeriodicRecv = NULL;
  Local_Point_PeriodicRecv = new unsigned long[nPoint_PeriodicRecv[nPeriodicRecv]];
  for (iRecv = 0; iRecv < nPoint_PeriodicRecv[nPeriodicRecv]; iRecv++)
    Local_Point_PeriodicRecv[iRecv] = 0;

  Local_Marker_PeriodicRecv = NULL;
  Local_Marker_PeriodicRecv = new unsigned long[nPoint_PeriodicRecv[nPeriodicRecv]];
  for (iRecv = 0; iRecv < nPoint_PeriodicRecv[nPeriodicRecv]; iRecv++)
    Local_Marker_PeriodicRecv[iRecv] = 0;

  /*--- We allocate the buffers for communicating values in a later step
   once we know the maximum packet size that we need to communicate. This
   memory is deallocated and reallocated automatically in the case that
   the previously allocated memory is not sufficient. ---*/

  bufD_PeriodicSend = NULL;
  bufD_PeriodicRecv = NULL;

  bufS_PeriodicSend = NULL;
  bufS_PeriodicRecv = NULL;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nPeriodicSend > 0) {
    req_PeriodicSend   = new SU2_MPI::Request[nPeriodicSend];
  }
  if (nPeriodicRecv > 0) {
    req_PeriodicRecv   = new SU2_MPI::Request[nPeriodicRecv];
  }

  /*--- Allocate arrays for sending the periodic point index and marker
   index to the recv rank so that it can store the local values. Therefore,
   the recv rank can quickly loop through the buffers to unpack the data. ---*/

  unsigned short nPackets = 2;
  unsigned long *idSend = new unsigned long[nPoint_PeriodicSend[nPeriodicSend]*nPackets];
  for (iSend = 0; iSend < nPoint_PeriodicSend[nPeriodicSend]*nPackets; iSend++)
    idSend[iSend] = 0;

  /*--- Build the lists of local index and periodic marker index values. ---*/

  ii = 0; jj = 0;
  for (iSend = 0; iSend < nPeriodicSend; iSend++) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

          /*--- Get the current periodic point index. We only communicate
           the owned nodes on a rank, as the MPI comms will take care of
           the halos after completing the periodic comms. ---*/

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->node[iPoint]->GetDomain()) {

            /*--- Get the rank that holds the matching periodic point
             on the other marker in the periodic pair. ---*/

            iRank = (int)geometry->vertex[iMarker][iVertex]->GetDonorProcessor();

            /*--- If the rank for the current periodic point matches the
             rank of the current send message, then store the local point
             index on the matching periodic point and the periodic marker
             index to be communicated to the recv rank. ---*/

            if (iRank == Neighbors_PeriodicSend[iSend]) {
              Local_Point_PeriodicSend[ii]  = iPoint;
              Local_Marker_PeriodicSend[ii] = (unsigned long)iMarker;
              jj = ii*nPackets;
              idSend[jj] = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
              jj++;
              idSend[jj] = (unsigned long)iPeriodic;
              ii++;
            }

          }
        }
      }
    }
  }

  /*--- Allocate arrays for receiving the periodic point index and marker
   index to the recv rank so that it can store the local values. ---*/

  unsigned long *idRecv = new unsigned long[nPoint_PeriodicRecv[nPeriodicRecv]*nPackets];
  for (iRecv = 0; iRecv < nPoint_PeriodicRecv[nPeriodicRecv]*nPackets; iRecv++)
    idRecv[iRecv] = 0;

#ifdef HAVE_MPI

  int iMessage, offset, count, source, dest, tag;

  /*--- Launch the non-blocking recv's first. Note that we have stored
   the counts and sources, so we can launch these before we even load
   the data and send from the periodically matching ranks. ---*/

  iMessage = 0;
  for (iRecv = 0; iRecv < nPeriodicRecv; iRecv++) {

    /*--- Compute our location in the recv buffer. ---*/

    offset = nPackets*nPoint_PeriodicRecv[iRecv];

    /*--- Take advantage of cumulative storage format to get the number
     of elems that we need to recv. ---*/

    count = nPackets*(nPoint_PeriodicRecv[iRecv+1] - nPoint_PeriodicRecv[iRecv]);

    /*--- Get the rank from which we receive the message. ---*/

    source = Neighbors_PeriodicRecv[iRecv];
    tag    = source + 1;

    /*--- Post non-blocking recv for this proc. ---*/

    SU2_MPI::Irecv(&(static_cast<unsigned long*>(idRecv)[offset]),
                   count, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD,
                   &(req_PeriodicRecv[iMessage]));

    /*--- Increment message counter. ---*/

    iMessage++;

  }

  /*--- Post the non-blocking sends. ---*/

  iMessage = 0;
  for (iSend = 0; iSend < nPeriodicSend; iSend++) {

    /*--- Compute our location in the send buffer. ---*/

    offset = nPackets*nPoint_PeriodicSend[iSend];

    /*--- Take advantage of cumulative storage format to get the number
     of points that we need to send. ---*/

    count = nPackets*(nPoint_PeriodicSend[iSend+1] - nPoint_PeriodicSend[iSend]);

    /*--- Get the rank to which we send the message. ---*/

    dest = Neighbors_PeriodicSend[iSend];
    tag  = rank + 1;

    /*--- Post non-blocking send for this proc. ---*/

    SU2_MPI::Isend(&(static_cast<unsigned long*>(idSend)[offset]),
                   count, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD,
                   &(req_PeriodicSend[iMessage]));

    /*--- Increment message counter. ---*/

    iMessage++;

  }

  /*--- Wait for the non-blocking comms to complete. ---*/

  SU2_MPI::Waitall(nPeriodicSend, req_PeriodicSend, MPI_STATUS_IGNORE);
  SU2_MPI::Waitall(nPeriodicRecv, req_PeriodicRecv, MPI_STATUS_IGNORE);

#else

  /*--- Copy my own rank's data into the recv buffer directly in serial. ---*/

  int myStart, myFinal;
  for (int val_iSend = 0; val_iSend < nPeriodicSend; val_iSend++) {
    iRank   = geometry->PeriodicRecv2Neighbor[rank];
    iRecv   = geometry->nPoint_PeriodicRecv[iRank]*nPackets;
    myStart = nPoint_PeriodicSend[val_iSend]*nPackets;
    myFinal = nPoint_PeriodicSend[val_iSend+1]*nPackets;
    for (iSend = myStart; iSend < myFinal; iSend++) {
      idRecv[iRecv] = idSend[iSend];
      iRecv++;
    }
  }

#endif

  /*--- Store the local periodic point and marker index values in our
   data structures so we can quickly unpack data during the iterations. ---*/

  ii = 0;
  for (iRecv = 0; iRecv < nPoint_PeriodicRecv[nPeriodicRecv]; iRecv++) {
    Local_Point_PeriodicRecv[iRecv]  = idRecv[ii]; ii++;
    Local_Marker_PeriodicRecv[iRecv] = idRecv[ii]; ii++;
  }

  delete [] idSend;
  delete [] idRecv;

}

void CGeometry::AllocatePeriodicComms(unsigned short val_countPerPeriodicPoint) {

  /*--- This routine is activated whenever we attempt to perform
   a periodic MPI communication with our neighbors but the
   memory buffer allocated is not large enough for the packet size.
   Therefore, we deallocate the previously allocated arrays and
   reallocate a large enough array. Note that after the first set
   communications, this routine will not need to be called again. ---*/

  int iSend, iRecv, nSend, nRecv;

  /*--- Store the larger packet size to the class data. ---*/

  countPerPeriodicPoint = val_countPerPeriodicPoint;

  /*--- Store the total size of the send/recv arrays for clarity. ---*/

  nSend = countPerPeriodicPoint*nPoint_PeriodicSend[nPeriodicSend];
  nRecv = countPerPeriodicPoint*nPoint_PeriodicRecv[nPeriodicRecv];

  /*-- Deallocate and reallocate our cummunication memory. ---*/

  if (bufD_PeriodicSend != NULL) delete [] bufD_PeriodicSend;

  bufD_PeriodicSend = new su2double[nSend];
  for (iSend = 0; iSend < nSend; iSend++)
    bufD_PeriodicSend[iSend] = 0.0;

  if (bufD_PeriodicRecv != NULL) delete [] bufD_PeriodicRecv;

  bufD_PeriodicRecv = new su2double[nRecv];
  for (iRecv = 0; iRecv < nRecv; iRecv++)
    bufD_PeriodicRecv[iRecv] = 0.0;

  if (bufS_PeriodicSend != NULL) delete [] bufS_PeriodicSend;

  bufS_PeriodicSend = new unsigned short[nSend];
  for (iSend = 0; iSend < nSend; iSend++)
    bufS_PeriodicSend[iSend] = 0;

  if (bufS_PeriodicRecv != NULL) delete [] bufS_PeriodicRecv;

  bufS_PeriodicRecv = new unsigned short[nRecv];
  for (iRecv = 0; iRecv < nRecv; iRecv++)
    bufS_PeriodicRecv[iRecv] = 0;

}

void CGeometry::PostPeriodicRecvs(CGeometry *geometry,
                                  CConfig *config,
                                  unsigned short commType) {

  /*--- In parallel, communicate the data with non-blocking send/recv. ---*/

#ifdef HAVE_MPI

  /*--- Local variables ---*/

  int iMessage, iRecv, offset, nPointPeriodic, count, source, tag;

  /*--- Launch the non-blocking recv's first. Note that we have stored
   the counts and sources, so we can launch these before we even load
   the data and send from the neighbor ranks. ---*/

  iMessage = 0;
  for (iRecv = 0; iRecv < nPeriodicRecv; iRecv++) {

    /*--- Compute our location in the recv buffer. ---*/

    offset = countPerPeriodicPoint*nPoint_PeriodicRecv[iRecv];

    /*--- Take advantage of cumulative storage format to get the number
     of elems that we need to recv. ---*/

    nPointPeriodic = nPoint_PeriodicRecv[iRecv+1] - nPoint_PeriodicRecv[iRecv];

    /*--- Total count can include multiple pieces of data per element. ---*/

    count = countPerPeriodicPoint*nPointPeriodic;

    /*--- Get the rank from which we receive the message. ---*/

    source = Neighbors_PeriodicRecv[iRecv];
    tag    = source + 1;

    /*--- Post non-blocking recv for this proc. ---*/

    switch (commType) {
      case COMM_TYPE_DOUBLE:
        SU2_MPI::Irecv(&(static_cast<su2double*>(bufD_PeriodicRecv)[offset]),
                       count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD,
                       &(req_PeriodicRecv[iMessage]));
        break;
      case COMM_TYPE_UNSIGNED_SHORT:
        SU2_MPI::Irecv(&(static_cast<unsigned short*>(bufS_PeriodicRecv)[offset]),
                       count, MPI_UNSIGNED_SHORT, source, tag, MPI_COMM_WORLD,
                       &(req_PeriodicRecv[iMessage]));
        break;
      default:
        SU2_MPI::Error("Unrecognized data type for periodic MPI comms.",
                       CURRENT_FUNCTION);
        break;
    }

    /*--- Increment message counter. ---*/

    iMessage++;

  }

#endif

}

void CGeometry::PostPeriodicSends(CGeometry *geometry,
                                  CConfig *config,
                                  unsigned short commType,
                                  int val_iSend) {

  /*--- In parallel, communicate the data with non-blocking send/recv. ---*/

#ifdef HAVE_MPI

  /*--- Local variables ---*/

  int iMessage, offset, nPointPeriodic, count, dest, tag;

  /*--- Post the non-blocking send as soon as the buffer is loaded. ---*/

  iMessage = val_iSend;

  /*--- Compute our location in the send buffer. ---*/

  offset = countPerPeriodicPoint*nPoint_PeriodicSend[val_iSend];

  /*--- Take advantage of cumulative storage format to get the number
   of points that we need to send. ---*/

  nPointPeriodic = (nPoint_PeriodicSend[val_iSend+1] -
                    nPoint_PeriodicSend[val_iSend]);

  /*--- Total count can include multiple pieces of data per element. ---*/

  count = countPerPeriodicPoint*nPointPeriodic;

  /*--- Get the rank to which we send the message. ---*/

  dest = Neighbors_PeriodicSend[val_iSend];
  tag  = rank + 1;

  /*--- Post non-blocking send for this proc. ---*/

  switch (commType) {
    case COMM_TYPE_DOUBLE:
      SU2_MPI::Isend(&(static_cast<su2double*>(bufD_PeriodicSend)[offset]),
                     count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD,
                     &(req_PeriodicSend[iMessage]));
      break;
    case COMM_TYPE_UNSIGNED_SHORT:
      SU2_MPI::Isend(&(static_cast<unsigned short*>(bufS_PeriodicSend)[offset]),
                     count, MPI_UNSIGNED_SHORT, dest, tag, MPI_COMM_WORLD,
                     &(req_PeriodicSend[iMessage]));
      break;
    default:
      SU2_MPI::Error("Unrecognized data type for periodic MPI comms.",
                     CURRENT_FUNCTION);
      break;
  }

#else

  /*--- Copy my own rank's data into the recv buffer directly in serial. ---*/

  int iSend, myStart, myFinal, iRecv, iRank;
  iRank   = geometry->PeriodicRecv2Neighbor[rank];
  iRecv   = geometry->nPoint_PeriodicRecv[iRank]*countPerPeriodicPoint;
  myStart = nPoint_PeriodicSend[val_iSend]*countPerPeriodicPoint;
  myFinal = nPoint_PeriodicSend[val_iSend+1]*countPerPeriodicPoint;
  for (iSend = myStart; iSend < myFinal; iSend++) {
    switch (commType) {
      case COMM_TYPE_DOUBLE:
        bufD_PeriodicRecv[iRecv] =  bufD_PeriodicSend[iSend];
        break;
      case COMM_TYPE_UNSIGNED_SHORT:
        bufS_PeriodicRecv[iRecv] =  bufS_PeriodicSend[iSend];
        break;
      default:
        SU2_MPI::Error("Unrecognized data type for periodic MPI comms.",
                       CURRENT_FUNCTION);
        break;
    }
    iRecv++;
  }

#endif

}

su2double CGeometry::Point2Plane_Distance(su2double *Coord, su2double *iCoord, su2double *jCoord, su2double *kCoord) {
  su2double CrossProduct[3], iVector[3], jVector[3], distance, modulus;
  unsigned short iDim;

  for (iDim = 0; iDim < 3; iDim ++) {
    iVector[iDim] = jCoord[iDim] - iCoord[iDim];
    jVector[iDim] = kCoord[iDim] - iCoord[iDim];
  }

  CrossProduct[0] = iVector[1]*jVector[2] - iVector[2]*jVector[1];
  CrossProduct[1] = iVector[2]*jVector[0] - iVector[0]*jVector[2];
  CrossProduct[2] = iVector[0]*jVector[1] - iVector[1]*jVector[0];

  modulus = sqrt(CrossProduct[0]*CrossProduct[0]+CrossProduct[1]*CrossProduct[1]+CrossProduct[2]*CrossProduct[2]);

  distance = 0.0;
  for (iDim = 0; iDim < 3; iDim ++)
    distance += CrossProduct[iDim]*(Coord[iDim]-iCoord[iDim]);
  distance /= modulus;

  return distance;

}

long CGeometry::FindEdge(unsigned long first_point, unsigned long second_point) {
  unsigned long iPoint = 0;
  unsigned short iNode;
  for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
    iPoint = node[first_point]->GetPoint(iNode);
    if (iPoint == second_point) break;
  }

  if (iPoint == second_point) return node[first_point]->GetEdge(iNode);
  else {
    char buf[100];
    SPRINTF(buf, "Can't find the edge that connects %lu and %lu.", first_point, second_point);
    SU2_MPI::Error(buf, CURRENT_FUNCTION);
    return 0;
  }
}

bool CGeometry::CheckEdge(unsigned long first_point, unsigned long second_point) {
  unsigned long iPoint = 0;
  unsigned short iNode;
  for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
    iPoint = node[first_point]->GetPoint(iNode);
    if (iPoint == second_point) break;
  }

  if (iPoint == second_point) return true;
  else return false;

}

void CGeometry::SetEdges(void) {
  unsigned long iPoint, jPoint;
  long iEdge;
  unsigned short jNode, iNode;
  long TestEdge = 0;

  nEdge = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
      jPoint = node[iPoint]->GetPoint(iNode);
      for (jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
        if (node[jPoint]->GetPoint(jNode) == iPoint) {
          TestEdge = node[jPoint]->GetEdge(jNode);
          break;
        }
      if (TestEdge == -1) {
        node[iPoint]->SetEdge(nEdge, iNode);
        node[jPoint]->SetEdge(nEdge, jNode);
        nEdge++;
      }
    }

  edge = new CEdge*[nEdge];

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
      jPoint = node[iPoint]->GetPoint(iNode);
      iEdge = FindEdge(iPoint, jPoint);
      if (iPoint < jPoint) edge[iEdge] = new CEdge(iPoint, jPoint, nDim);
    }
}

void CGeometry::SetFaces(void) {
  //	unsigned long iPoint, jPoint, iFace;
  //	unsigned short jNode, iNode;
  //	long TestFace = 0;
  //
  //	nFace = 0;
  //	for (iPoint = 0; iPoint < nPoint; iPoint++)
  //		for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
  //			jPoint = node[iPoint]->GetPoint(iNode);
  //			for (jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
  //				if (node[jPoint]->GetPoint(jNode) == iPoint) {
  //					TestFace = node[jPoint]->GetFace(jNode);
  //					break;
  //				}
  //			if (TestFace == -1) {
  //				node[iPoint]->SetFace(nFace, iNode);
  //				node[jPoint]->SetFace(nFace, jNode);
  //				nFace++;
  //			}
  //		}
  //
  //	face = new CFace*[nFace];
  //
  //	for (iPoint = 0; iPoint < nPoint; iPoint++)
  //		for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
  //			jPoint = node[iPoint]->GetPoint(iNode);
  //			iFace = FindFace(iPoint, jPoint);
  //			if (iPoint < jPoint) face[iFace] = new CFace(iPoint, jPoint, nDim);
  //		}
}

void CGeometry::TestGeometry(void) {

  ofstream para_file;

  para_file.open("test_geometry.dat", ios::out);

  su2double *Normal = new su2double[nDim];

  for (unsigned long iEdge = 0; iEdge < nEdge; iEdge++) {
    para_file << "Edge index: " << iEdge << endl;
    para_file << "   Point index: " << edge[iEdge]->GetNode(0) << "\t" << edge[iEdge]->GetNode(1) << endl;
    edge[iEdge]->GetNormal(Normal);
    para_file << "      Face normal : ";
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      para_file << Normal[iDim] << "\t";
    para_file << endl;
  }

  para_file << endl;
  para_file << endl;
  para_file << endl;
  para_file << endl;

  for (unsigned short iMarker =0; iMarker < nMarker; iMarker++) {
    para_file << "Marker index: " << iMarker << endl;
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      para_file << "   Vertex index: " << iVertex << endl;
      para_file << "      Point index: " << vertex[iMarker][iVertex]->GetNode() << endl;
      para_file << "      Point coordinates : ";
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        para_file << node[vertex[iMarker][iVertex]->GetNode()]->GetCoord(iDim) << "\t";}
      para_file << endl;
      vertex[iMarker][iVertex]->GetNormal(Normal);
      para_file << "         Face normal : ";
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        para_file << Normal[iDim] << "\t";
      para_file << endl;
    }
  }

  delete [] Normal;

}

void CGeometry::SetSpline(vector<su2double> &x, vector<su2double> &y, unsigned long n, su2double yp1, su2double ypn, vector<su2double> &y2) {
  unsigned long i, k;
  su2double p, qn, sig, un, *u;

  u = new su2double [n];

  if (yp1 > 0.99e30)			// The lower boundary condition is set either to be "nat
    y2[0]=u[0]=0.0;			  // -ural"
  else {									// or else to have a specified first derivative.
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  for (i=2; i<=n-1; i++) {									//  This is the decomposition loop of the tridiagonal al-
    sig=(x[i-1]-x[i-2])/(x[i]-x[i-2]);		//	gorithm. y2 and u are used for tem-
    p=sig*y2[i-2]+2.0;										//	porary storage of the decomposed
    y2[i-1]=(sig-1.0)/p;										//	factors.

    su2double a1 = (y[i]-y[i-1])/(x[i]-x[i-1]); if (x[i] == x[i-1]) a1 = 1.0;
    su2double a2 = (y[i-1]-y[i-2])/(x[i-1]-x[i-2]); if (x[i-1] == x[i-2]) a2 = 1.0;
    u[i-1]= a1 - a2;
    u[i-1]=(6.0*u[i-1]/(x[i]-x[i-2])-sig*u[i-2])/p;

  }

  if (ypn > 0.99e30)						// The upper boundary condition is set either to be
    qn=un=0.0;									// "natural"
  else {												// or else to have a specified first derivative.
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-1; k>=1; k--)					// This is the backsubstitution loop of the tridiagonal
    y2[k-1]=y2[k-1]*y2[k]+u[k-1];	  // algorithm.

  delete[] u;

}

su2double CGeometry::GetSpline(vector<su2double>&xa, vector<su2double>&ya, vector<su2double>&y2a, unsigned long n, su2double x) {
  unsigned long klo, khi, k;
  su2double h, b, a, y;

  if (x < xa[0]) x = xa[0];       // Clip max and min values
  if (x > xa[n-1]) x = xa[n-1];

  klo = 1;										// We will find the right place in the table by means of
  khi = n;										// bisection. This is optimal if sequential calls to this
  while (khi-klo > 1) {			// routine are at random values of x. If sequential calls
    k = (khi+klo) >> 1;				// are in order, and closely spaced, one would do better
    if (xa[k-1] > x) khi = k;		// to store previous values of klo and khi and test if
    else klo=k;							// they remain appropriate on the next call.
  }								// klo and khi now bracket the input value of x
  h = xa[khi-1] - xa[klo-1];
  if (h == 0.0) h = EPS; // cout << "Bad xa input to routine splint" << endl;	// The xa?s must be distinct.
  a = (xa[khi-1]-x)/h;
  b = (x-xa[klo-1])/h;				// Cubic spline polynomial is now evaluated.
  y = a*ya[klo-1]+b*ya[khi-1]+((a*a*a-a)*y2a[klo-1]+(b*b*b-b)*y2a[khi-1])*(h*h)/6.0;

  return y;
}

bool CGeometry::SegmentIntersectsPlane(su2double *Segment_P0, su2double *Segment_P1, su2double Variable_P0, su2double Variable_P1,
                                                           su2double *Plane_P0, su2double *Plane_Normal, su2double *Intersection, su2double &Variable_Interp) {
  su2double u[3], v[3], Denominator, Numerator, Aux, ModU;
  su2double epsilon = 1E-6; // An epsilon is added to eliminate, as much as possible, the posibility of a line that intersects a point
  unsigned short iDim;

  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P1[iDim] - Segment_P0[iDim];
    v[iDim] = (Plane_P0[iDim]+epsilon) - Segment_P0[iDim];
  }

  ModU = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);

  Numerator = (Plane_Normal[0]+epsilon)*v[0] + (Plane_Normal[1]+epsilon)*v[1] + (Plane_Normal[2]+epsilon)*v[2];
  Denominator = (Plane_Normal[0]+epsilon)*u[0] + (Plane_Normal[1]+epsilon)*u[1] + (Plane_Normal[2]+epsilon)*u[2];

  if (fabs(Denominator) <= 0.0) return (false); // No intersection.

  Aux = Numerator / Denominator;

  if (Aux < 0.0 || Aux > 1.0) return (false); // No intersection.

  for (iDim = 0; iDim < 3; iDim++)
    Intersection[iDim] = Segment_P0[iDim] + Aux * u[iDim];


  /*--- Check that the intersection is in the segment ---*/

  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P0[iDim] - Intersection[iDim];
    v[iDim] = Segment_P1[iDim] - Intersection[iDim];
  }

  Variable_Interp = Variable_P0 + (Variable_P1 - Variable_P0)*sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])/ModU;

  Denominator = (Plane_Normal[0]+epsilon)*u[0] + (Plane_Normal[1]+epsilon)*u[1] + (Plane_Normal[2]+epsilon)*u[2];
  Numerator = (Plane_Normal[0]+epsilon)*v[0] + (Plane_Normal[1]+epsilon)*v[1] + (Plane_Normal[2]+epsilon)*v[2];

  Aux = Numerator * Denominator;

  if (Aux > 0.0) return (false); // Intersection outside the segment.

  return (true);

}

bool CGeometry::RayIntersectsTriangle(su2double orig[3], su2double dir[3],
                                      su2double vert0[3], su2double vert1[3], su2double vert2[3],
                                      su2double *intersect) {

  const passivedouble epsilon = 0.000001;
  su2double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  su2double det, inv_det, t, u, v;

  /*--- Find vectors for two edges sharing vert0 ---*/

  SUB(edge1, vert1, vert0);
  SUB(edge2, vert2, vert0);

  /*--- Begin calculating determinant - also used to calculate U parameter ---*/

  CROSS(pvec, dir, edge2);

  /*--- If determinant is near zero, ray lies in plane of triangle ---*/

  det = DOT(edge1, pvec);


  if (fabs(det) < epsilon) return(false);

  inv_det = 1.0 / det;

  /*--- Calculate distance from vert0 to ray origin ---*/

  SUB(tvec, orig, vert0);

  /*--- Calculate U parameter and test bounds ---*/

  u = inv_det * DOT(tvec, pvec);

  if (u < 0.0 || u > 1.0) return(false);

  /*--- prepare to test V parameter ---*/

  CROSS(qvec, tvec, edge1);

  /*--- Calculate V parameter and test bounds ---*/

  v = inv_det * DOT(dir, qvec);

  if (v < 0.0 || u + v > 1.0) return(false);

  /*--- Calculate t, ray intersects triangle ---*/

  t = inv_det * DOT(edge2, qvec);

  /*--- Compute the intersection point in cartesian coordinates ---*/

  intersect[0] = orig[0] + (t * dir[0]);
  intersect[1] = orig[1] + (t * dir[1]);
  intersect[2] = orig[2] + (t * dir[2]);

  return (true);

}

bool CGeometry::SegmentIntersectsLine(su2double point0[2], su2double point1[2], su2double vert0[2], su2double vert1[2]) {

  su2double det, diff0_A, diff0_B, diff1_A, diff1_B, intersect[2];

  diff0_A = point0[0] - point1[0];
  diff1_A = point0[1] - point1[1];

  diff0_B = vert0[0] - vert1[0];
  diff1_B = vert0[1] - vert1[1];

  det = (diff0_A)*(diff1_B) - (diff1_A)*(diff0_B);

  if (det == 0) return false;

  /*--- Compute point of intersection ---*/

  intersect[0] = ((point0[0]*point1[1] - point0[1]*point1[0])*diff0_B
                -(vert0[0]* vert1[1]  - vert0[1]* vert1[0])*diff0_A)/det;

  intersect[1] =  ((point0[0]*point1[1] - point0[1]*point1[0])*diff1_B
                  -(vert0[0]* vert1[1]  - vert0[1]* vert1[0])*diff1_A)/det;


  /*--- Check that the point is between the two surface points ---*/

  su2double dist0, dist1, length;

  dist0 = (intersect[0] - point0[0])*(intersect[0] - point0[0])
         +(intersect[1] - point0[1])*(intersect[1] - point0[1]);

  dist1 = (intersect[0] - point1[0])*(intersect[0] - point1[0])
         +(intersect[1] - point1[1])*(intersect[1] - point1[1]);

  length = diff0_A*diff0_A
          +diff1_A*diff1_A;

  if ( (dist0 > length) || (dist1 > length) ) {
    return false;
  }

  return true;
}

bool CGeometry::SegmentIntersectsTriangle(su2double point0[3], su2double point1[3],
                                          su2double vert0[3], su2double vert1[3], su2double vert2[3]) {

  su2double dir[3], intersect[3], u[3], v[3], edge1[3], edge2[3], Plane_Normal[3], Denominator, Numerator, Aux;

  SUB(dir, point1, point0);

  if (RayIntersectsTriangle(point0, dir, vert0, vert1, vert2, intersect)) {

    /*--- Check that the intersection is in the segment ---*/

    SUB(u, point0, intersect);
    SUB(v, point1, intersect);

    SUB(edge1, vert1, vert0);
    SUB(edge2, vert2, vert0);
    CROSS(Plane_Normal, edge1, edge2);

    Denominator = DOT(Plane_Normal, u);
    Numerator = DOT(Plane_Normal, v);

    Aux = Numerator * Denominator;

    /*--- Intersection outside the segment ---*/

    if (Aux > 0.0) return (false);

  }
  else {

    /*--- No intersection with the ray ---*/

    return (false);

  }

  /*--- Intersection inside the segment ---*/

  return (true);

}

void CGeometry::ComputeAirfoil_Section(su2double *Plane_P0, su2double *Plane_Normal,
                                       su2double MinXCoord, su2double MaxXCoord,
                                       su2double MinYCoord, su2double MaxYCoord,
                                       su2double MinZCoord, su2double MaxZCoord,
                                       su2double *FlowVariable,
                                       vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil,
                                       vector<su2double> &Zcoord_Airfoil, vector<su2double> &Variable_Airfoil,
                                       bool original_surface, CConfig *config) {

  AD_BEGIN_PASSIVE

  unsigned short iMarker, iNode, jNode, iDim, Index = 0;
  bool intersect;
  long Next_Edge = 0;
  unsigned long iPoint, jPoint, iElem, Trailing_Point, Airfoil_Point, iVertex, iEdge, PointIndex, jEdge;
  su2double Segment_P0[3] = {0.0, 0.0, 0.0}, Segment_P1[3] = {0.0, 0.0, 0.0}, Variable_P0 = 0.0, Variable_P1 = 0.0, Intersection[3] = {0.0, 0.0, 0.0}, Trailing_Coord,
  *VarCoord = NULL, Variable_Interp, v1[3] = {0.0, 0.0, 0.0}, v3[3] = {0.0, 0.0, 0.0}, CrossProduct = 1.0;
  bool Found_Edge;
  passivedouble Dist_Value;
  vector<su2double> Xcoord_Index0, Ycoord_Index0, Zcoord_Index0, Variable_Index0, Xcoord_Index1, Ycoord_Index1, Zcoord_Index1, Variable_Index1;
  vector<unsigned long> IGlobalID_Index0, JGlobalID_Index0, IGlobalID_Index1, JGlobalID_Index1, IGlobalID_Airfoil, JGlobalID_Airfoil;
  vector<unsigned short> Conection_Index0, Conection_Index1;
  vector<unsigned long> Duplicate;
  vector<unsigned long>::iterator it;
  su2double **Coord_Variation = NULL;
  vector<su2double> XcoordExtra, YcoordExtra, ZcoordExtra, VariableExtra;
  vector<unsigned long> IGlobalIDExtra, JGlobalIDExtra;
  vector<bool> AddExtra;
  unsigned long EdgeDonor;
  bool FoundEdge;

#ifdef HAVE_MPI
  unsigned long nLocalEdge, MaxLocalEdge, *Buffer_Send_nEdge, *Buffer_Receive_nEdge, nBuffer_Coord, nBuffer_Variable, nBuffer_GlobalID;
  int nProcessor, iProcessor;
  su2double *Buffer_Send_Coord, *Buffer_Receive_Coord;
  su2double *Buffer_Send_Variable, *Buffer_Receive_Variable;
  unsigned long *Buffer_Send_GlobalID, *Buffer_Receive_GlobalID;
#endif

  Xcoord_Airfoil.clear();
  Ycoord_Airfoil.clear();
  Zcoord_Airfoil.clear();
  Variable_Airfoil.clear();
  IGlobalID_Airfoil.clear();
  JGlobalID_Airfoil.clear();

  /*--- Set the right plane in 2D (note the change in Y-Z plane) ---*/

  if (nDim == 2) {
    Plane_P0[0] = 0.0;      Plane_P0[1] = 0.0;      Plane_P0[2] = 0.0;
    Plane_Normal[0] = 0.0;  Plane_Normal[1] = 1.0;  Plane_Normal[2] = 0.0;
  }

  /*--- Grid movement is stored using a vertices information,
   we should go from vertex to points ---*/

  if (original_surface == false) {

    Coord_Variation = new su2double *[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      Coord_Variation[iPoint] = new su2double [nDim];

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_GeoEval(iMarker) == YES) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          VarCoord = vertex[iMarker][iVertex]->GetVarCoord();
          iPoint = vertex[iMarker][iVertex]->GetNode();
          for (iDim = 0; iDim < nDim; iDim++)
            Coord_Variation[iPoint][iDim] = VarCoord[iDim];
        }
      }
    }

  }

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_GeoEval(iMarker) == YES) {

      for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {

        PointIndex=0;

        /*--- To decide if an element is going to be used or not should be done element based,
         The first step is to compute and average coordinate for the element ---*/

        su2double AveXCoord = 0.0;
        su2double AveYCoord = 0.0;
        su2double AveZCoord = 0.0;

        for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
          iPoint = bound[iMarker][iElem]->GetNode(iNode);
          AveXCoord += node[iPoint]->GetCoord(0);
          AveYCoord += node[iPoint]->GetCoord(1);
          if (nDim == 3) AveZCoord += node[iPoint]->GetCoord(2);
        }

        AveXCoord /= su2double(bound[iMarker][iElem]->GetnNodes());
        AveYCoord /= su2double(bound[iMarker][iElem]->GetnNodes());
        AveZCoord /= su2double(bound[iMarker][iElem]->GetnNodes());

        /*--- To only cut one part of the nacelle based on the cross product
         of the normal to the plane and a vector that connect the point
         with the center line ---*/

        CrossProduct = 1.0;

        if (config->GetGeo_Description() == NACELLE) {

          su2double Tilt_Angle = config->GetNacelleLocation(3)*PI_NUMBER/180;
          su2double Toe_Angle = config->GetNacelleLocation(4)*PI_NUMBER/180;

          /*--- Translate to the origin ---*/

          su2double XCoord_Trans = AveXCoord - config->GetNacelleLocation(0);
          su2double YCoord_Trans = AveYCoord - config->GetNacelleLocation(1);
          su2double ZCoord_Trans = AveZCoord - config->GetNacelleLocation(2);

          /*--- Apply tilt angle ---*/

          su2double XCoord_Trans_Tilt = XCoord_Trans*cos(Tilt_Angle) + ZCoord_Trans*sin(Tilt_Angle);
          su2double YCoord_Trans_Tilt = YCoord_Trans;
          su2double ZCoord_Trans_Tilt = ZCoord_Trans*cos(Tilt_Angle) - XCoord_Trans*sin(Tilt_Angle);

          /*--- Apply toe angle ---*/

          su2double YCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt*sin(Toe_Angle) + YCoord_Trans_Tilt*cos(Toe_Angle);
          su2double ZCoord_Trans_Tilt_Toe = ZCoord_Trans_Tilt;

          /*--- Undo plane rotation, we have already rotated the nacelle ---*/

          /*--- Undo tilt angle ---*/

          su2double XPlane_Normal_Tilt = Plane_Normal[0]*cos(-Tilt_Angle) + Plane_Normal[2]*sin(-Tilt_Angle);
          su2double YPlane_Normal_Tilt = Plane_Normal[1];
          su2double ZPlane_Normal_Tilt = Plane_Normal[2]*cos(-Tilt_Angle) - Plane_Normal[0]*sin(-Tilt_Angle);

          /*--- Undo toe angle ---*/

          su2double YPlane_Normal_Tilt_Toe = XPlane_Normal_Tilt*sin(-Toe_Angle) + YPlane_Normal_Tilt*cos(-Toe_Angle);
          su2double ZPlane_Normal_Tilt_Toe = ZPlane_Normal_Tilt;


          v1[1] = YCoord_Trans_Tilt_Toe - 0.0;
          v1[2] = ZCoord_Trans_Tilt_Toe - 0.0;
          v3[0] = v1[1]*ZPlane_Normal_Tilt_Toe-v1[2]*YPlane_Normal_Tilt_Toe;
          CrossProduct = v3[0] * 1.0;

        }

        for (unsigned short iFace = 0; iFace < bound[iMarker][iElem]->GetnFaces(); iFace++){
          iNode = bound[iMarker][iElem]->GetFaces(iFace,0);
          jNode = bound[iMarker][iElem]->GetFaces(iFace,1);
          iPoint = bound[iMarker][iElem]->GetNode(iNode);
          jPoint = bound[iMarker][iElem]->GetNode(jNode);

          if ((CrossProduct >= 0.0)
              && ((AveXCoord > MinXCoord) && (AveXCoord < MaxXCoord))
              && ((AveYCoord > MinYCoord) && (AveYCoord < MaxYCoord))
              && ((AveZCoord > MinZCoord) && (AveZCoord < MaxZCoord))) {

            Segment_P0[0] = 0.0;  Segment_P0[1] = 0.0;  Segment_P0[2] = 0.0;  Variable_P0 = 0.0;
            Segment_P1[0] = 0.0;  Segment_P1[1] = 0.0;  Segment_P1[2] = 0.0;  Variable_P1 = 0.0;


            for (iDim = 0; iDim < nDim; iDim++) {
              if (original_surface == true) {
                Segment_P0[iDim] = node[iPoint]->GetCoord(iDim);
                Segment_P1[iDim] = node[jPoint]->GetCoord(iDim);
              }
              else {
                Segment_P0[iDim] = node[iPoint]->GetCoord(iDim) + Coord_Variation[iPoint][iDim];
                Segment_P1[iDim] = node[jPoint]->GetCoord(iDim) + Coord_Variation[jPoint][iDim];
              }
            }

            if (FlowVariable != NULL) {
              Variable_P0 = FlowVariable[iPoint];
              Variable_P1 = FlowVariable[jPoint];
            }

            /*--- In 2D add the points directly (note the change between Y and Z coordinate) ---*/

            if (nDim == 2) {
              Xcoord_Index0.push_back(Segment_P0[0]);                     Xcoord_Index1.push_back(Segment_P1[0]);
              Ycoord_Index0.push_back(Segment_P0[2]);                     Ycoord_Index1.push_back(Segment_P1[2]);
              Zcoord_Index0.push_back(Segment_P0[1]);                     Zcoord_Index1.push_back(Segment_P1[1]);
              Variable_Index0.push_back(Variable_P0);                     Variable_Index1.push_back(Variable_P1);
              IGlobalID_Index0.push_back(node[iPoint]->GetGlobalIndex()); IGlobalID_Index1.push_back(node[jPoint]->GetGlobalIndex());
              JGlobalID_Index0.push_back(node[iPoint]->GetGlobalIndex()); JGlobalID_Index1.push_back(node[jPoint]->GetGlobalIndex());
              PointIndex++;
            }

            /*--- In 3D compute the intersection ---*/

            else if (nDim == 3) {
              intersect = SegmentIntersectsPlane(Segment_P0, Segment_P1, Variable_P0, Variable_P1, Plane_P0, Plane_Normal, Intersection, Variable_Interp);
              if (intersect == true) {
                if (PointIndex == 0) {
                  Xcoord_Index0.push_back(Intersection[0]);
                  Ycoord_Index0.push_back(Intersection[1]);
                  Zcoord_Index0.push_back(Intersection[2]);
                  Variable_Index0.push_back(Variable_Interp);
                  IGlobalID_Index0.push_back(node[iPoint]->GetGlobalIndex());
                  JGlobalID_Index0.push_back(node[jPoint]->GetGlobalIndex());
                }
                if (PointIndex == 1) {
                  Xcoord_Index1.push_back(Intersection[0]);
                  Ycoord_Index1.push_back(Intersection[1]);
                  Zcoord_Index1.push_back(Intersection[2]);
                  Variable_Index1.push_back(Variable_Interp);
                  IGlobalID_Index1.push_back(node[iPoint]->GetGlobalIndex());
                  JGlobalID_Index1.push_back(node[jPoint]->GetGlobalIndex());
                }
                PointIndex++;
              }
            }
          }
        }
      }
    }
  }

  if (original_surface == false) {
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      delete [] Coord_Variation[iPoint];
    delete [] Coord_Variation;
  }

#ifdef HAVE_MPI

  /*--- Copy the coordinates of all the points in the plane to the master node ---*/

  nLocalEdge = 0, MaxLocalEdge = 0;
  nProcessor = size;

  Buffer_Send_nEdge = new unsigned long [1];
  Buffer_Receive_nEdge = new unsigned long [nProcessor];

  nLocalEdge = Xcoord_Index0.size();

  Buffer_Send_nEdge[0] = nLocalEdge;

  SU2_MPI::Allreduce(&nLocalEdge, &MaxLocalEdge, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nEdge, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nEdge, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  Buffer_Send_Coord    = new su2double [MaxLocalEdge*6];
  Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalEdge*6];

  Buffer_Send_Variable    = new su2double [MaxLocalEdge*2];
  Buffer_Receive_Variable = new su2double [nProcessor*MaxLocalEdge*2];

  Buffer_Send_GlobalID    = new unsigned long [MaxLocalEdge*4];
  Buffer_Receive_GlobalID = new unsigned long [nProcessor*MaxLocalEdge*4];

  nBuffer_Coord    = MaxLocalEdge*6;
  nBuffer_Variable = MaxLocalEdge*2;
  nBuffer_GlobalID = MaxLocalEdge*4;

  for (iEdge = 0; iEdge < nLocalEdge; iEdge++) {
    Buffer_Send_Coord[iEdge*6 + 0] = Xcoord_Index0[iEdge];
    Buffer_Send_Coord[iEdge*6 + 1] = Ycoord_Index0[iEdge];
    Buffer_Send_Coord[iEdge*6 + 2] = Zcoord_Index0[iEdge];
    Buffer_Send_Coord[iEdge*6 + 3] = Xcoord_Index1[iEdge];
    Buffer_Send_Coord[iEdge*6 + 4] = Ycoord_Index1[iEdge];
    Buffer_Send_Coord[iEdge*6 + 5] = Zcoord_Index1[iEdge];

    Buffer_Send_Variable[iEdge*2 + 0] = Variable_Index0[iEdge];
    Buffer_Send_Variable[iEdge*2 + 1] = Variable_Index1[iEdge];

    Buffer_Send_GlobalID[iEdge*4 + 0] = IGlobalID_Index0[iEdge];
    Buffer_Send_GlobalID[iEdge*4 + 1] = JGlobalID_Index0[iEdge];
    Buffer_Send_GlobalID[iEdge*4 + 2] = IGlobalID_Index1[iEdge];
    Buffer_Send_GlobalID[iEdge*4 + 3] = JGlobalID_Index1[iEdge];
  }

  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_Variable, nBuffer_Variable, MPI_DOUBLE, Buffer_Receive_Variable, nBuffer_Variable, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_GlobalID, nBuffer_GlobalID, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalID, nBuffer_GlobalID, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /*--- Clean the vectors before adding the new vertices only to the master node ---*/

  Xcoord_Index0.clear();     Xcoord_Index1.clear();
  Ycoord_Index0.clear();     Ycoord_Index1.clear();
  Zcoord_Index0.clear();     Zcoord_Index1.clear();
  Variable_Index0.clear();   Variable_Index1.clear();
  IGlobalID_Index0.clear();  IGlobalID_Index1.clear();
  JGlobalID_Index0.clear();  JGlobalID_Index1.clear();

  /*--- Copy the boundary to the master node vectors ---*/

  if (rank == MASTER_NODE) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iEdge = 0; iEdge < Buffer_Receive_nEdge[iProcessor]; iEdge++) {
        Xcoord_Index0.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalEdge*6 + iEdge*6 + 0] );
        Ycoord_Index0.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalEdge*6 + iEdge*6 + 1] );
        Zcoord_Index0.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalEdge*6 + iEdge*6 + 2] );
        Xcoord_Index1.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalEdge*6 + iEdge*6 + 3] );
        Ycoord_Index1.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalEdge*6 + iEdge*6 + 4] );
        Zcoord_Index1.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalEdge*6 + iEdge*6 + 5] );

        Variable_Index0.push_back( Buffer_Receive_Variable[ iProcessor*MaxLocalEdge*2 + iEdge*2 + 0] );
        Variable_Index1.push_back( Buffer_Receive_Variable[ iProcessor*MaxLocalEdge*2 + iEdge*2 + 1] );

        IGlobalID_Index0.push_back( Buffer_Receive_GlobalID[ iProcessor*MaxLocalEdge*4 + iEdge*4 + 0] );
        JGlobalID_Index0.push_back( Buffer_Receive_GlobalID[ iProcessor*MaxLocalEdge*4 + iEdge*4 + 1] );
        IGlobalID_Index1.push_back( Buffer_Receive_GlobalID[ iProcessor*MaxLocalEdge*4 + iEdge*4 + 2] );
        JGlobalID_Index1.push_back( Buffer_Receive_GlobalID[ iProcessor*MaxLocalEdge*4 + iEdge*4 + 3] );

      }
    }
  }

  delete[] Buffer_Send_Coord;      delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_Variable;   delete[] Buffer_Receive_Variable;
  delete[] Buffer_Send_GlobalID;   delete[] Buffer_Receive_GlobalID;
  delete[] Buffer_Send_nEdge;    delete[] Buffer_Receive_nEdge;

#endif

  if ((rank == MASTER_NODE) && (Xcoord_Index0.size() != 0)) {

    /*--- Remove singular edges ---*/

    bool Remove;

    do { Remove = false;
      for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {

        if (((IGlobalID_Index0[iEdge] == IGlobalID_Index1[iEdge]) && (JGlobalID_Index0[iEdge] == JGlobalID_Index1[iEdge])) ||
            ((IGlobalID_Index0[iEdge] == JGlobalID_Index1[iEdge]) && (JGlobalID_Index0[iEdge] == IGlobalID_Index1[iEdge]))) {

          Xcoord_Index0.erase (Xcoord_Index0.begin() + iEdge);
          Ycoord_Index0.erase (Ycoord_Index0.begin() + iEdge);
          Zcoord_Index0.erase (Zcoord_Index0.begin() + iEdge);
          Variable_Index0.erase (Variable_Index0.begin() + iEdge);
          IGlobalID_Index0.erase (IGlobalID_Index0.begin() + iEdge);
          JGlobalID_Index0.erase (JGlobalID_Index0.begin() + iEdge);

          Xcoord_Index1.erase (Xcoord_Index1.begin() + iEdge);
          Ycoord_Index1.erase (Ycoord_Index1.begin() + iEdge);
          Zcoord_Index1.erase (Zcoord_Index1.begin() + iEdge);
          Variable_Index1.erase (Variable_Index1.begin() + iEdge);
          IGlobalID_Index1.erase (IGlobalID_Index1.begin() + iEdge);
          JGlobalID_Index1.erase (JGlobalID_Index1.begin() + iEdge);

          Remove = true; break;
        }
        if (Remove) break;
      }
    } while (Remove == true);

    /*--- Remove repeated edges computing distance, this could happend because the MPI ---*/

    do { Remove = false;
      for (iEdge = 0; iEdge < Xcoord_Index0.size()-1; iEdge++) {
        for (jEdge = iEdge+1; jEdge < Xcoord_Index0.size(); jEdge++) {

          /*--- Edges with the same orientation ---*/

          if ((((IGlobalID_Index0[iEdge] == IGlobalID_Index0[jEdge]) && (JGlobalID_Index0[iEdge] == JGlobalID_Index0[jEdge])) ||
               ((IGlobalID_Index0[iEdge] == JGlobalID_Index0[jEdge]) && (JGlobalID_Index0[iEdge] == IGlobalID_Index0[jEdge]))) &&
              (((IGlobalID_Index1[iEdge] == IGlobalID_Index1[jEdge]) && (JGlobalID_Index1[iEdge] == JGlobalID_Index1[jEdge])) ||
               ((IGlobalID_Index1[iEdge] == JGlobalID_Index1[jEdge]) && (JGlobalID_Index1[iEdge] == IGlobalID_Index1[jEdge])))) {

                Xcoord_Index0.erase (Xcoord_Index0.begin() + jEdge);
                Ycoord_Index0.erase (Ycoord_Index0.begin() + jEdge);
                Zcoord_Index0.erase (Zcoord_Index0.begin() + jEdge);
                Variable_Index0.erase (Variable_Index0.begin() + jEdge);
                IGlobalID_Index0.erase (IGlobalID_Index0.begin() + jEdge);
                JGlobalID_Index0.erase (JGlobalID_Index0.begin() + jEdge);

                Xcoord_Index1.erase (Xcoord_Index1.begin() + jEdge);
                Ycoord_Index1.erase (Ycoord_Index1.begin() + jEdge);
                Zcoord_Index1.erase (Zcoord_Index1.begin() + jEdge);
                Variable_Index1.erase (Variable_Index1.begin() + jEdge);
                IGlobalID_Index1.erase (IGlobalID_Index1.begin() + jEdge);
                JGlobalID_Index1.erase (JGlobalID_Index1.begin() + jEdge);

                Remove = true; break;

              }

          /*--- Edges with oposite orientation ---*/

          if ((((IGlobalID_Index0[iEdge] == IGlobalID_Index1[jEdge]) && (JGlobalID_Index0[iEdge] == JGlobalID_Index1[jEdge])) ||
               ((IGlobalID_Index0[iEdge] == JGlobalID_Index1[jEdge]) && (JGlobalID_Index0[iEdge] == IGlobalID_Index1[jEdge]))) &&
              (((IGlobalID_Index1[iEdge] == IGlobalID_Index0[jEdge]) && (JGlobalID_Index1[iEdge] == JGlobalID_Index0[jEdge])) ||
               ((IGlobalID_Index1[iEdge] == JGlobalID_Index0[jEdge]) && (JGlobalID_Index1[iEdge] == IGlobalID_Index0[jEdge])))) {

                Xcoord_Index0.erase (Xcoord_Index0.begin() + jEdge);
                Ycoord_Index0.erase (Ycoord_Index0.begin() + jEdge);
                Zcoord_Index0.erase (Zcoord_Index0.begin() + jEdge);
                Variable_Index0.erase (Variable_Index0.begin() + jEdge);
                IGlobalID_Index0.erase (IGlobalID_Index0.begin() + jEdge);
                JGlobalID_Index0.erase (JGlobalID_Index0.begin() + jEdge);

                Xcoord_Index1.erase (Xcoord_Index1.begin() + jEdge);
                Ycoord_Index1.erase (Ycoord_Index1.begin() + jEdge);
                Zcoord_Index1.erase (Zcoord_Index1.begin() + jEdge);
                Variable_Index1.erase (Variable_Index1.begin() + jEdge);
                IGlobalID_Index1.erase (IGlobalID_Index1.begin() + jEdge);
                JGlobalID_Index1.erase (JGlobalID_Index1.begin() + jEdge);

                Remove = true; break;
              }
          if (Remove) break;
        }
        if (Remove) break;
      }

    } while (Remove == true);

    if (Xcoord_Index0.size() != 1) {

      /*--- Rotate from the Y-Z plane to the X-Z plane to reuse the rest of subroutines  ---*/

      if (config->GetGeo_Description() == FUSELAGE) {
        su2double Angle = -0.5*PI_NUMBER;
        for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {
          su2double XCoord = Xcoord_Index0[iEdge]*cos(Angle) - Ycoord_Index0[iEdge]*sin(Angle);
          su2double YCoord = Ycoord_Index0[iEdge]*cos(Angle) + Xcoord_Index0[iEdge]*sin(Angle);
          su2double ZCoord = Zcoord_Index0[iEdge];
          Xcoord_Index0[iEdge] = XCoord; Ycoord_Index0[iEdge] = YCoord; Zcoord_Index0[iEdge] = ZCoord;
          XCoord = Xcoord_Index1[iEdge]*cos(Angle) - Ycoord_Index1[iEdge]*sin(Angle);
          YCoord = Ycoord_Index1[iEdge]*cos(Angle) + Xcoord_Index1[iEdge]*sin(Angle);
          ZCoord = Zcoord_Index1[iEdge];
          Xcoord_Index1[iEdge] = XCoord; Ycoord_Index1[iEdge] = YCoord; Zcoord_Index1[iEdge] = ZCoord;
        }
      }

      /*--- Rotate nacelle secction to a X-Z plane to reuse the rest of subroutines  ---*/


      if (config->GetGeo_Description() == NACELLE) {

        su2double Tilt_Angle = config->GetNacelleLocation(3)*PI_NUMBER/180;
        su2double Toe_Angle = config->GetNacelleLocation(4)*PI_NUMBER/180;
        su2double Theta_deg = atan2(Plane_Normal[1],-Plane_Normal[2])/PI_NUMBER*180 + 180;
        su2double Roll_Angle = 0.5*PI_NUMBER - Theta_deg*PI_NUMBER/180;

        su2double XCoord_Trans, YCoord_Trans, ZCoord_Trans, XCoord_Trans_Tilt, YCoord_Trans_Tilt, ZCoord_Trans_Tilt,
        XCoord_Trans_Tilt_Toe, YCoord_Trans_Tilt_Toe, ZCoord_Trans_Tilt_Toe, XCoord, YCoord, ZCoord;

        for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {

          /*--- First point of the edge ---*/

          /*--- Translate to the origin ---*/

          XCoord_Trans = Xcoord_Index0[iEdge] - config->GetNacelleLocation(0);
          YCoord_Trans = Ycoord_Index0[iEdge] - config->GetNacelleLocation(1);
          ZCoord_Trans = Zcoord_Index0[iEdge] - config->GetNacelleLocation(2);

          /*--- Apply tilt angle ---*/

          XCoord_Trans_Tilt = XCoord_Trans*cos(Tilt_Angle) + ZCoord_Trans*sin(Tilt_Angle);
          YCoord_Trans_Tilt = YCoord_Trans;
          ZCoord_Trans_Tilt = ZCoord_Trans*cos(Tilt_Angle) - XCoord_Trans*sin(Tilt_Angle);

          /*--- Apply toe angle ---*/

          XCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt*cos(Toe_Angle) - YCoord_Trans_Tilt*sin(Toe_Angle);
          YCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt*sin(Toe_Angle) + YCoord_Trans_Tilt*cos(Toe_Angle);
          ZCoord_Trans_Tilt_Toe = ZCoord_Trans_Tilt;

          /*--- Rotate to X-Z plane (roll) ---*/

          XCoord = XCoord_Trans_Tilt_Toe;
          YCoord = YCoord_Trans_Tilt_Toe*cos(Roll_Angle) - ZCoord_Trans_Tilt_Toe*sin(Roll_Angle);
          ZCoord = YCoord_Trans_Tilt_Toe*sin(Roll_Angle) + ZCoord_Trans_Tilt_Toe*cos(Roll_Angle);

          /*--- Update coordinates ---*/

          Xcoord_Index0[iEdge] = XCoord; Ycoord_Index0[iEdge] = YCoord; Zcoord_Index0[iEdge] = ZCoord;

          /*--- Second point of the edge ---*/

          /*--- Translate to the origin ---*/

          XCoord_Trans = Xcoord_Index1[iEdge] - config->GetNacelleLocation(0);
          YCoord_Trans = Ycoord_Index1[iEdge] - config->GetNacelleLocation(1);
          ZCoord_Trans = Zcoord_Index1[iEdge] - config->GetNacelleLocation(2);

          /*--- Apply tilt angle ---*/

          XCoord_Trans_Tilt = XCoord_Trans*cos(Tilt_Angle) + ZCoord_Trans*sin(Tilt_Angle);
          YCoord_Trans_Tilt = YCoord_Trans;
          ZCoord_Trans_Tilt = ZCoord_Trans*cos(Tilt_Angle) - XCoord_Trans*sin(Tilt_Angle);

          /*--- Apply toe angle ---*/

          XCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt*cos(Toe_Angle) - YCoord_Trans_Tilt*sin(Toe_Angle);
          YCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt*sin(Toe_Angle) + YCoord_Trans_Tilt*cos(Toe_Angle);
          ZCoord_Trans_Tilt_Toe = ZCoord_Trans_Tilt;

          /*--- Rotate to X-Z plane (roll) ---*/

          XCoord = XCoord_Trans_Tilt_Toe;
          YCoord = YCoord_Trans_Tilt_Toe*cos(Roll_Angle) - ZCoord_Trans_Tilt_Toe*sin(Roll_Angle);
          ZCoord = YCoord_Trans_Tilt_Toe*sin(Roll_Angle) + ZCoord_Trans_Tilt_Toe*cos(Roll_Angle);

          /*--- Update coordinates ---*/

          Xcoord_Index1[iEdge] = XCoord; Ycoord_Index1[iEdge] = YCoord; Zcoord_Index1[iEdge] = ZCoord;

        }
      }


      /*--- Identify the extreme of the curve and close it ---*/

      Conection_Index0.reserve(Xcoord_Index0.size()+1);
      Conection_Index1.reserve(Xcoord_Index0.size()+1);

      for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {
        Conection_Index0[iEdge] = 0;
        Conection_Index1[iEdge] = 0;
      }

      for (iEdge = 0; iEdge < Xcoord_Index0.size()-1; iEdge++) {
        for (jEdge = iEdge+1; jEdge < Xcoord_Index0.size(); jEdge++) {

          if (((IGlobalID_Index0[iEdge] == IGlobalID_Index0[jEdge]) && (JGlobalID_Index0[iEdge] == JGlobalID_Index0[jEdge])) ||
              ((IGlobalID_Index0[iEdge] == JGlobalID_Index0[jEdge]) && (JGlobalID_Index0[iEdge] == IGlobalID_Index0[jEdge])))
          { Conection_Index0[iEdge]++; Conection_Index0[jEdge]++; }

          if (((IGlobalID_Index0[iEdge] == IGlobalID_Index1[jEdge]) && (JGlobalID_Index0[iEdge] == JGlobalID_Index1[jEdge])) ||
              ((IGlobalID_Index0[iEdge] == JGlobalID_Index1[jEdge]) && (JGlobalID_Index0[iEdge] == IGlobalID_Index1[jEdge])))
          { Conection_Index0[iEdge]++; Conection_Index1[jEdge]++; }

          if (((IGlobalID_Index1[iEdge] == IGlobalID_Index0[jEdge]) && (JGlobalID_Index1[iEdge] == JGlobalID_Index0[jEdge])) ||
              ((IGlobalID_Index1[iEdge] == JGlobalID_Index0[jEdge]) && (JGlobalID_Index1[iEdge] == IGlobalID_Index0[jEdge])))
          { Conection_Index1[iEdge]++; Conection_Index0[jEdge]++; }

          if (((IGlobalID_Index1[iEdge] == IGlobalID_Index1[jEdge]) && (JGlobalID_Index1[iEdge] == JGlobalID_Index1[jEdge])) ||
              ((IGlobalID_Index1[iEdge] == JGlobalID_Index1[jEdge]) && (JGlobalID_Index1[iEdge] == IGlobalID_Index1[jEdge])))
          { Conection_Index1[iEdge]++; Conection_Index1[jEdge]++; }

        }
      }

      /*--- Connect extremes of the curves ---*/

      /*--- First: Identify the extremes of the curve in the extra vector  ---*/

      for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {
        if (Conection_Index0[iEdge] == 0) {
          XcoordExtra.push_back(Xcoord_Index0[iEdge]);
          YcoordExtra.push_back(Ycoord_Index0[iEdge]);
          ZcoordExtra.push_back(Zcoord_Index0[iEdge]);
          VariableExtra.push_back(Variable_Index0[iEdge]);
          IGlobalIDExtra.push_back(IGlobalID_Index0[iEdge]);
          JGlobalIDExtra.push_back(JGlobalID_Index0[iEdge]);
          AddExtra.push_back(true);
        }
        if (Conection_Index1[iEdge] == 0) {
          XcoordExtra.push_back(Xcoord_Index1[iEdge]);
          YcoordExtra.push_back(Ycoord_Index1[iEdge]);
          ZcoordExtra.push_back(Zcoord_Index1[iEdge]);
          VariableExtra.push_back(Variable_Index1[iEdge]);
          IGlobalIDExtra.push_back(IGlobalID_Index1[iEdge]);
          JGlobalIDExtra.push_back(JGlobalID_Index1[iEdge]);
          AddExtra.push_back(true);
        }
      }

      /*--- Second, if it is an open curve then find the closest point to an extreme to close it  ---*/

      if (XcoordExtra.size() > 1) {

        for (iEdge = 0; iEdge < XcoordExtra.size()-1; iEdge++) {

          su2double MinDist = 1E6; FoundEdge = false; EdgeDonor = 0;
          for (jEdge = iEdge+1; jEdge < XcoordExtra.size(); jEdge++) {
            Dist_Value = sqrt(pow(SU2_TYPE::GetValue(XcoordExtra[iEdge])-SU2_TYPE::GetValue(XcoordExtra[jEdge]), 2.0));
            if ((Dist_Value < MinDist) && (AddExtra[iEdge]) && (AddExtra[jEdge])) {
              EdgeDonor = jEdge; FoundEdge = true;
            }
          }

          if (FoundEdge) {

            /*--- Add first point of the new edge ---*/

            Xcoord_Index0.push_back (XcoordExtra[iEdge]);
            Ycoord_Index0.push_back (YcoordExtra[iEdge]);
            Zcoord_Index0.push_back (ZcoordExtra[iEdge]);
            Variable_Index0.push_back (VariableExtra[iEdge]);
            IGlobalID_Index0.push_back (IGlobalIDExtra[iEdge]);
            JGlobalID_Index0.push_back (JGlobalIDExtra[iEdge]);
            AddExtra[iEdge] = false;

            /*--- Add second (closest)  point of the new edge ---*/

            Xcoord_Index1.push_back (XcoordExtra[EdgeDonor]);
            Ycoord_Index1.push_back (YcoordExtra[EdgeDonor]);
            Zcoord_Index1.push_back (ZcoordExtra[EdgeDonor]);
            Variable_Index1.push_back (VariableExtra[EdgeDonor]);
            IGlobalID_Index1.push_back (IGlobalIDExtra[EdgeDonor]);
            JGlobalID_Index1.push_back (JGlobalIDExtra[EdgeDonor]);
            AddExtra[EdgeDonor] = false;

          }

        }

      }

      else if (XcoordExtra.size() == 1) {
        cout <<"There cutting system has failed, there is an incomplete curve (not used)." << endl;
      }

      /*--- Find and add the trailing edge to to the list
       and the contect the first point to the trailing edge ---*/

      Trailing_Point = 0; Trailing_Coord = Xcoord_Index0[0];
      for (iEdge = 1; iEdge < Xcoord_Index0.size(); iEdge++) {
        if (Xcoord_Index0[iEdge] > Trailing_Coord) {
          Trailing_Point = iEdge; Trailing_Coord = Xcoord_Index0[iEdge];
        }
      }

      Xcoord_Airfoil.push_back(Xcoord_Index0[Trailing_Point]);
      Ycoord_Airfoil.push_back(Ycoord_Index0[Trailing_Point]);
      Zcoord_Airfoil.push_back(Zcoord_Index0[Trailing_Point]);
      Variable_Airfoil.push_back(Variable_Index0[Trailing_Point]);
      IGlobalID_Airfoil.push_back(IGlobalID_Index0[Trailing_Point]);
      JGlobalID_Airfoil.push_back(JGlobalID_Index0[Trailing_Point]);

      Xcoord_Airfoil.push_back(Xcoord_Index1[Trailing_Point]);
      Ycoord_Airfoil.push_back(Ycoord_Index1[Trailing_Point]);
      Zcoord_Airfoil.push_back(Zcoord_Index1[Trailing_Point]);
      Variable_Airfoil.push_back(Variable_Index1[Trailing_Point]);
      IGlobalID_Airfoil.push_back(IGlobalID_Index1[Trailing_Point]);
      JGlobalID_Airfoil.push_back(JGlobalID_Index1[Trailing_Point]);

      Xcoord_Index0.erase (Xcoord_Index0.begin() + Trailing_Point);
      Ycoord_Index0.erase (Ycoord_Index0.begin() + Trailing_Point);
      Zcoord_Index0.erase (Zcoord_Index0.begin() + Trailing_Point);
      Variable_Index0.erase (Variable_Index0.begin() + Trailing_Point);
      IGlobalID_Index0.erase (IGlobalID_Index0.begin() + Trailing_Point);
      JGlobalID_Index0.erase (JGlobalID_Index0.begin() + Trailing_Point);

      Xcoord_Index1.erase (Xcoord_Index1.begin() + Trailing_Point);
      Ycoord_Index1.erase (Ycoord_Index1.begin() + Trailing_Point);
      Zcoord_Index1.erase (Zcoord_Index1.begin() + Trailing_Point);
      Variable_Index1.erase (Variable_Index1.begin() + Trailing_Point);
      IGlobalID_Index1.erase (IGlobalID_Index1.begin() + Trailing_Point);
      JGlobalID_Index1.erase (JGlobalID_Index1.begin() + Trailing_Point);


      /*--- Algorithm for adding the rest of the points ---*/

      do {

        /*--- Last added point in the list ---*/

        Airfoil_Point = Xcoord_Airfoil.size() - 1;

        /*--- Find the closest point  ---*/

        Found_Edge = false;

        for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {

          if (((IGlobalID_Index0[iEdge] == IGlobalID_Airfoil[Airfoil_Point]) && (JGlobalID_Index0[iEdge] == JGlobalID_Airfoil[Airfoil_Point])) ||
              ((IGlobalID_Index0[iEdge] == JGlobalID_Airfoil[Airfoil_Point]) && (JGlobalID_Index0[iEdge] == IGlobalID_Airfoil[Airfoil_Point]))) {
            Next_Edge = iEdge; Found_Edge = true; Index = 0; break;
          }

          if (((IGlobalID_Index1[iEdge] == IGlobalID_Airfoil[Airfoil_Point]) && (JGlobalID_Index1[iEdge] == JGlobalID_Airfoil[Airfoil_Point])) ||
              ((IGlobalID_Index1[iEdge] == JGlobalID_Airfoil[Airfoil_Point]) && (JGlobalID_Index1[iEdge] == IGlobalID_Airfoil[Airfoil_Point]))) {
            Next_Edge = iEdge; Found_Edge = true; Index = 1; break;
          }

        }

        /*--- Add and remove the next point to the list and the next point in the edge ---*/

        if (Found_Edge) {

          if (Index == 0) {
            Xcoord_Airfoil.push_back(Xcoord_Index1[Next_Edge]);
            Ycoord_Airfoil.push_back(Ycoord_Index1[Next_Edge]);
            Zcoord_Airfoil.push_back(Zcoord_Index1[Next_Edge]);
            Variable_Airfoil.push_back(Variable_Index1[Next_Edge]);
            IGlobalID_Airfoil.push_back(IGlobalID_Index1[Next_Edge]);
            JGlobalID_Airfoil.push_back(JGlobalID_Index1[Next_Edge]);
          }

          if (Index == 1) {
            Xcoord_Airfoil.push_back(Xcoord_Index0[Next_Edge]);
            Ycoord_Airfoil.push_back(Ycoord_Index0[Next_Edge]);
            Zcoord_Airfoil.push_back(Zcoord_Index0[Next_Edge]);
            Variable_Airfoil.push_back(Variable_Index0[Next_Edge]);
            IGlobalID_Airfoil.push_back(IGlobalID_Index0[Next_Edge]);
            JGlobalID_Airfoil.push_back(JGlobalID_Index0[Next_Edge]);
          }

          Xcoord_Index0.erase(Xcoord_Index0.begin() + Next_Edge);
          Ycoord_Index0.erase(Ycoord_Index0.begin() + Next_Edge);
          Zcoord_Index0.erase(Zcoord_Index0.begin() + Next_Edge);
          Variable_Index0.erase(Variable_Index0.begin() + Next_Edge);
          IGlobalID_Index0.erase(IGlobalID_Index0.begin() + Next_Edge);
          JGlobalID_Index0.erase(JGlobalID_Index0.begin() + Next_Edge);

          Xcoord_Index1.erase(Xcoord_Index1.begin() + Next_Edge);
          Ycoord_Index1.erase(Ycoord_Index1.begin() + Next_Edge);
          Zcoord_Index1.erase(Zcoord_Index1.begin() + Next_Edge);
          Variable_Index1.erase(Variable_Index1.begin() + Next_Edge);
          IGlobalID_Index1.erase(IGlobalID_Index1.begin() + Next_Edge);
          JGlobalID_Index1.erase(JGlobalID_Index1.begin() + Next_Edge);

        }
        else { break; }

      } while (Xcoord_Index0.size() != 0);

      /*--- Clean the vector before using them again for storing the upper or the lower side ---*/

      Xcoord_Index0.clear(); Ycoord_Index0.clear(); Zcoord_Index0.clear(); Variable_Index0.clear();  IGlobalID_Index0.clear();  JGlobalID_Index0.clear();
      Xcoord_Index1.clear(); Ycoord_Index1.clear(); Zcoord_Index1.clear(); Variable_Index1.clear();  IGlobalID_Index1.clear();  JGlobalID_Index1.clear();

    }

  }

  AD_END_PASSIVE

}

void CGeometry::RegisterCoordinates(CConfig *config) {
  unsigned short iDim;
  unsigned long iPoint;
  bool input = true;
  bool push_index = config->GetMultizone_Problem()? false : true;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      AD::RegisterInput(node[iPoint]->GetCoord()[iDim], push_index);
    }
    if(!push_index) {
      for (iDim = 0; iDim < nDim; iDim++) {
        node[iPoint]->SetIndex(input);
      }
    }
  }
}

void CGeometry::RegisterOutput_Coordinates(CConfig *config){
  unsigned short iDim;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    if(config->GetMultizone_Problem()) {
      for (iDim = 0; iDim < nDim; iDim++) {
        AD::RegisterOutput(node[iPoint]->GetCoord()[iDim]);
      }
    }
    else {
      for (iDim = 0; iDim < nDim; iDim++) {
        AD::RegisterOutput(node[iPoint]->GetCoord()[iDim]);
      }
    }
  }
}

void CGeometry::UpdateGeometry(CGeometry **geometry_container, CConfig *config) {

  unsigned short iMesh;

  geometry_container[MESH_0]->InitiateComms(geometry_container[MESH_0], config, COORDINATES);
  geometry_container[MESH_0]->CompleteComms(geometry_container[MESH_0], config, COORDINATES);
  if (config->GetGrid_Movement()){
    geometry_container[MESH_0]->InitiateComms(geometry_container[MESH_0], config, GRID_VELOCITY);
    geometry_container[MESH_0]->CompleteComms(geometry_container[MESH_0], config, GRID_VELOCITY);
  }

  geometry_container[MESH_0]->SetCoord_CG();
  geometry_container[MESH_0]->SetControlVolume(config, UPDATE);
  geometry_container[MESH_0]->SetBoundControlVolume(config, UPDATE);
  geometry_container[MESH_0]->SetMaxLength(config);

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    /*--- Update the control volume structures ---*/

    geometry_container[iMesh]->SetControlVolume(config,geometry_container[iMesh-1], UPDATE);
    geometry_container[iMesh]->SetBoundControlVolume(config,geometry_container[iMesh-1], UPDATE);
    geometry_container[iMesh]->SetCoord(geometry_container[iMesh-1]);

  }

  if (config->GetKind_Solver() == DISC_ADJ_RANS || config->GetKind_Solver() == DISC_ADJ_INC_RANS)
  geometry_container[MESH_0]->ComputeWall_Distance(config);

}

void CGeometry::SetCustomBoundary(CConfig *config) {

  unsigned short iMarker;
  unsigned long iVertex;
  string Marker_Tag;

  /*--- Initialize quantities for customized boundary conditions.
   * Custom values are initialized with the default values specified in the config (avoiding non physical values) ---*/
  CustomBoundaryTemperature = new su2double*[nMarker];
  CustomBoundaryHeatFlux = new su2double*[nMarker];

  for(iMarker=0; iMarker < nMarker; iMarker++){
    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    CustomBoundaryHeatFlux[iMarker] = NULL;
    CustomBoundaryTemperature[iMarker] = NULL;
    if(config->GetMarker_All_PyCustom(iMarker)){
      switch(config->GetMarker_All_KindBC(iMarker)){
        case HEAT_FLUX:
          CustomBoundaryHeatFlux[iMarker] = new su2double[nVertex[iMarker]];
          for(iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
            CustomBoundaryHeatFlux[iMarker][iVertex] = config->GetWall_HeatFlux(Marker_Tag);
          }
          break;
        case ISOTHERMAL:
          CustomBoundaryTemperature[iMarker] = new su2double[nVertex[iMarker]];
          for(iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
            CustomBoundaryTemperature[iMarker][iVertex] = config->GetIsothermal_Temperature(Marker_Tag);
          }
          break;
        case INLET_FLOW:
          // This case is handled in the solver class.
          break;
        default:
          cout << "WARNING: Marker " << Marker_Tag << " is not customizable. Using default behavior." << endl;
          break;
      }
    }
  }

}

void CGeometry::UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config){

  unsigned short iMGfine, iMGlevel, nMGlevel, iMarker;

  nMGlevel = config->GetnMGLevels();
  for (iMGlevel=1; iMGlevel <= nMGlevel; iMGlevel++){
    iMGfine = iMGlevel-1;
    for(iMarker = 0; iMarker< config->GetnMarker_All(); iMarker++){
      if(config->GetMarker_All_PyCustom(iMarker)){
        switch(config->GetMarker_All_KindBC(iMarker)){
          case HEAT_FLUX:
            geometry_container[iMGlevel]->SetMultiGridWallHeatFlux(geometry_container[iMGfine], iMarker);
            break;
          case ISOTHERMAL:
            geometry_container[iMGlevel]->SetMultiGridWallTemperature(geometry_container[iMGfine], iMarker);
            break;
          // Inlet flow handled in solver class.
          default: break;
        }
      }
    }
  }
}


void CGeometry::ComputeSurf_Straightness(CConfig *config,
                                         bool    print_on_screen) {

  bool RefUnitNormal_defined;
  unsigned short iDim,
                 iMarker,
                 iMarker_Global,
                 nMarker_Global = config->GetnMarker_CfgFile();
  unsigned long iVertex;
  constexpr passivedouble epsilon = 1.0e-6;
  su2double Area;
  string Local_TagBound,
         Global_TagBound;

  vector<su2double> Normal(nDim),
                    UnitNormal(nDim),
                    RefUnitNormal(nDim);

  /*--- Assume now that this boundary marker is straight. As soon as one
        AreaElement is found that is not aligend with a Reference then it is
        certain that the boundary marker is not straight and one can stop
        searching. Another possibility is that this process doesn't own
        any nodes of that boundary, in that case we also have to assume the
        boundary is straight.
        Any boundary type other than SYMMETRY_PLANE or EULER_WALL gets
        the value false (or see cases specified in the conditional below)
        which could be wrong. ---*/
  bound_is_straight.resize(nMarker);
  fill(bound_is_straight.begin(), bound_is_straight.end(), true);

  /*--- Loop over all local markers ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Local_TagBound = config->GetMarker_All_TagBound(iMarker);

    /*--- Marker has to be Symmetry or Euler. Additionally marker can't be a
          moving surface and Grid Movement Elasticity is forbidden as well. All
          other GridMovements are rigid. ---*/
    if ((config->GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE ||
         config->GetMarker_All_KindBC(iMarker) == EULER_WALL) &&
        config->GetMarker_Moving_Bool(Local_TagBound) == false &&
        config->GetKind_GridMovement() != ELASTICITY) {

      /*--- Loop over all global markers, and find the local-global pair via
            matching unique string tags. ---*/
      for (iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++) {

        Global_TagBound = config->GetMarker_CfgFile_TagBound(iMarker_Global);
        if (Local_TagBound == Global_TagBound) {

          RefUnitNormal_defined = false;
          iVertex = 0;

          while(bound_is_straight[iMarker] == true &&
                iVertex < nVertex[iMarker]) {

            vertex[iMarker][iVertex]->GetNormal(Normal.data());
            UnitNormal = Normal;

            /*--- Compute unit normal. ---*/
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);

            /*--- Negate for outward convention. ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              UnitNormal[iDim] /= -Area;

            /*--- Check if unit normal is within tolerance of the Reference unit normal.
                  Reference unit normal = first unit normal found. ---*/
            if(RefUnitNormal_defined) {
              for (iDim = 0; iDim < nDim; iDim++) {
                if( abs(RefUnitNormal[iDim] - UnitNormal[iDim]) > epsilon ) {
                  bound_is_straight[iMarker] = false;
                  break;
                }
              }
            } else {
              RefUnitNormal = UnitNormal; //deep copy of values
              RefUnitNormal_defined = true;
            }

          iVertex++;
          }//while iVertex
        }//if Local == Global
      }//for iMarker_Global
    } else {
      /*--- Enforce default value: false ---*/
      bound_is_straight[iMarker] = false;
    }//if sym or euler ...
  }//for iMarker

  /*--- Communicate results and print on screen. ---*/
  if(print_on_screen) {

    /*--- Additional vector which can later be MPI::Allreduce(d) to pring the results
          on screen as nMarker (local) can vary across ranks. Default 'true' as it can
          happen that a local rank does not contain an element of each surface marker.  ---*/
    vector<bool> bound_is_straight_Global(nMarker_Global, true);
    /*--- Match local with global tag bound and fill a Global Marker vector. ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Local_TagBound = config->GetMarker_All_TagBound(iMarker);
      for (iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++) {
        Global_TagBound = config->GetMarker_CfgFile_TagBound(iMarker_Global);

        if(Local_TagBound == Global_TagBound)
          bound_is_straight_Global[iMarker_Global] = bound_is_straight[iMarker];

      }//for iMarker_Global
    }//for iMarker

    vector<int> Buff_Send_isStraight(nMarker_Global),
                Buff_Recv_isStraight(nMarker_Global);

    /*--- Cast to int as std::vector<boolean> can be a special construct. MPI handling using <int>
          is more straight-forward. ---*/
    for (iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++)
      Buff_Send_isStraight[iMarker_Global] = static_cast<int> (bound_is_straight_Global[iMarker_Global]);

    /*--- Product of type <int>(bool) is equivalnt to a 'logical and' ---*/
    SU2_MPI::Allreduce(Buff_Send_isStraight.data(), Buff_Recv_isStraight.data(),
                       nMarker_Global, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

    /*--- Print results on screen. ---*/
    if(rank == MASTER_NODE) {
      for (iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++) {
        if (config->GetMarker_CfgFile_KindBC(config->GetMarker_CfgFile_TagBound(iMarker_Global)) == SYMMETRY_PLANE ||
          config->GetMarker_CfgFile_KindBC(config->GetMarker_CfgFile_TagBound(iMarker_Global)) == EULER_WALL) {

          cout << "Boundary marker " << config->GetMarker_CfgFile_TagBound(iMarker_Global) << " is";
          if(Buff_Recv_isStraight[iMarker_Global] == false) cout << " NOT";
          if(nDim == 2) cout << " a single straight." << endl;
          if(nDim == 3) cout << " a single plane." << endl;
        }//if sym or euler
      }//for iMarker_Global
    }//if rank==MASTER
  }//if print_on_scren

}


void CGeometry::ComputeSurf_Curvature(CConfig *config) {
  unsigned short iMarker, iNeigh_Point, iDim, iNode, iNeighbor_Nodes, Neighbor_Node;
  unsigned long Neighbor_Point, iVertex, iPoint, jPoint, iElem_Bound, iEdge, nLocalVertex, MaxLocalVertex , *Buffer_Send_nVertex, *Buffer_Receive_nVertex, TotalnPointDomain;
  int iProcessor, nProcessor;
  vector<unsigned long> Point_NeighborList, Elem_NeighborList, Point_Triangle, Point_Edge, Point_Critical;
  vector<unsigned long>::iterator it;
  su2double U[3] = {0.0,0.0,0.0}, V[3] = {0.0,0.0,0.0}, W[3] = {0.0,0.0,0.0}, Length_U, Length_V, Length_W, CosValue, Angle_Value, *K, *Angle_Defect, *Area_Vertex, *Angle_Alpha, *Angle_Beta, **NormalMeanK, MeanK, GaussK, MaxPrinK, cot_alpha, cot_beta, delta, X1, X2, X3, Y1, Y2, Y3, radius, *Buffer_Send_Coord, *Buffer_Receive_Coord, *Coord, Dist, MinDist, MaxK, MinK, SigmaK;
  bool *Check_Edge;

  bool fea = ((config->GetKind_Solver()==FEM_ELASTICITY) || (config->GetKind_Solver()==DISC_ADJ_FEM));

  /*--- Allocate surface curvature ---*/
  K = new su2double [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) K[iPoint] = 0.0;

  if (nDim == 2) {

    /*--- Loop over all the markers ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {

        /*--- Loop through all marker vertices again, this time also
         finding the neighbors of each node.---*/
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint  = vertex[iMarker][iVertex]->GetNode();

          if (node[iPoint]->GetDomain()) {
            /*--- Loop through neighbors. In 2-D, there should be 2 nodes on either
             side of this vertex that lie on the same surface. ---*/
            Point_Edge.clear();

            for (iNeigh_Point = 0; iNeigh_Point < node[iPoint]->GetnPoint(); iNeigh_Point++) {
              Neighbor_Point = node[iPoint]->GetPoint(iNeigh_Point);

              /*--- Check if this neighbor lies on the surface. If so,
               add to the list of neighbors. ---*/
              if (node[Neighbor_Point]->GetPhysicalBoundary()) {
                Point_Edge.push_back(Neighbor_Point);
              }

            }

            if (Point_Edge.size() == 2) {

              /*--- Compute the curvature using three points ---*/
              X1 = node[iPoint]->GetCoord(0);
              X2 = node[Point_Edge[0]]->GetCoord(0);
              X3 = node[Point_Edge[1]]->GetCoord(0);
              Y1 = node[iPoint]->GetCoord(1);
              Y2 = node[Point_Edge[0]]->GetCoord(1);
              Y3 = node[Point_Edge[1]]->GetCoord(1);

              radius = sqrt(((X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1))*
                            ((X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3))*
                            ((X3-X1)*(X3-X1) + (Y3-Y1)*(Y3-Y1)))/
              (2.0*fabs(X1*Y2+X2*Y3+X3*Y1-X1*Y3-X2*Y1-X3*Y2)+EPS);

              K[iPoint] = 1.0/radius;
              node[iPoint]->SetCurvature(K[iPoint]);
            }

          }

        }

      }

    }

  }

  else {

    Angle_Defect = new su2double [nPoint];
    Area_Vertex = new su2double [nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      Angle_Defect[iPoint] = 2*PI_NUMBER;
      Area_Vertex[iPoint] = 0.0;
    }

    Angle_Alpha = new su2double [nEdge];
    Angle_Beta = new su2double [nEdge];
    Check_Edge = new bool [nEdge];
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
      Angle_Alpha[iEdge] = 0.0;
      Angle_Beta[iEdge] = 0.0;
      Check_Edge[iEdge] = true;
    }

    NormalMeanK = new su2double *[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      NormalMeanK[iPoint] = new su2double [nDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        NormalMeanK[iPoint][iDim] = 0.0;
      }
    }

    /*--- Loop over all the markers ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {

        /*--- Loop over all the boundary elements ---*/
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {

          /*--- Only triangles ---*/
          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE) {

            /*--- Loop over all the nodes of the boundary element ---*/
            for (iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++) {

              iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);

              Point_Triangle.clear();

              for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                Neighbor_Point = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                Point_Triangle.push_back(Neighbor_Point);
              }

              iEdge = FindEdge(Point_Triangle[0], Point_Triangle[1]);

              for (iDim = 0; iDim < nDim; iDim++) {
                U[iDim] = node[Point_Triangle[0]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
                V[iDim] = node[Point_Triangle[1]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
              }

              W[0] = 0.5*(U[1]*V[2]-U[2]*V[1]); W[1] = -0.5*(U[0]*V[2]-U[2]*V[0]); W[2] = 0.5*(U[0]*V[1]-U[1]*V[0]);

              Length_U = 0.0; Length_V = 0.0; Length_W = 0.0; CosValue = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) { Length_U += U[iDim]*U[iDim]; Length_V += V[iDim]*V[iDim]; Length_W += W[iDim]*W[iDim]; }
              Length_U = sqrt(Length_U); Length_V = sqrt(Length_V); Length_W = sqrt(Length_W);
              for (iDim = 0; iDim < nDim; iDim++) { U[iDim] /= Length_U; V[iDim] /= Length_V; CosValue += U[iDim]*V[iDim]; }
              if (CosValue >= 1.0) CosValue = 1.0;
              if (CosValue <= -1.0) CosValue = -1.0;

              Angle_Value = acos(CosValue);
              Area_Vertex[iPoint] += Length_W;
              Angle_Defect[iPoint] -= Angle_Value;
              if (Angle_Alpha[iEdge] == 0.0) Angle_Alpha[iEdge] = Angle_Value;
              else Angle_Beta[iEdge] = Angle_Value;

            }
          }
        }
      }
    }

    /*--- Compute mean curvature ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE) {
            for (iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++) {
              iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);

              for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                jPoint = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);

                iEdge = FindEdge(iPoint, jPoint);

                if (Check_Edge[iEdge]) {

                  Check_Edge[iEdge] = false;

                  if (tan(Angle_Alpha[iEdge]) != 0.0) cot_alpha = 1.0/tan(Angle_Alpha[iEdge]); else cot_alpha = 0.0;
                  if (tan(Angle_Beta[iEdge]) != 0.0) cot_beta = 1.0/tan(Angle_Beta[iEdge]); else cot_beta = 0.0;

                  /*--- iPoint, and jPoint ---*/
                  for (iDim = 0; iDim < nDim; iDim++) {
                    if (Area_Vertex[iPoint] != 0.0) NormalMeanK[iPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[iPoint]->GetCoord(iDim) - node[jPoint]->GetCoord(iDim)) / Area_Vertex[iPoint];
                    if (Area_Vertex[jPoint] != 0.0) NormalMeanK[jPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[jPoint]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim)) / Area_Vertex[jPoint];
                  }
                }

              }
            }
          }
        }
      }
    }

    /*--- Compute Gauss, mean, max and min principal curvature,
     and set the list of critical points ---*/

    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint  = vertex[iMarker][iVertex]->GetNode();

          if (node[iPoint]->GetDomain()) {

            if (Area_Vertex[iPoint] != 0.0) GaussK = 3.0*Angle_Defect[iPoint]/Area_Vertex[iPoint];
            else GaussK = 0.0;

            MeanK = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              MeanK += NormalMeanK[iPoint][iDim]*NormalMeanK[iPoint][iDim];
            MeanK = sqrt(MeanK);

            delta = max((MeanK*MeanK - GaussK), 0.0);

            MaxPrinK = MeanK + sqrt(delta);

            /*--- Store the curvature value ---*/
            K[iPoint] = MaxPrinK;
            node[iPoint]->SetCurvature(K[iPoint]);
          }

        }
      }
    }

    delete [] Angle_Defect;
    delete [] Area_Vertex;
    delete [] Angle_Alpha;
    delete [] Angle_Beta;
    delete [] Check_Edge;

    for (iPoint = 0; iPoint < nPoint; iPoint++)
      delete [] NormalMeanK[iPoint];
    delete [] NormalMeanK;

  }

  /*--- Sharp edge detection is based in the statistical
   distribution of the curvature ---*/

  MaxK = K[0]; MinK = K[0]; MeanK = 0.0; TotalnPointDomain = 0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          MaxK = max(MaxK, fabs(K[iPoint]));
          MinK = min(MinK, fabs(K[iPoint]));
          MeanK += fabs(K[iPoint]);
          TotalnPointDomain++;
        }
      }
    }
  }

#ifdef HAVE_MPI
  su2double MyMeanK = MeanK; MeanK = 0.0;
  su2double MyMaxK = MaxK; MaxK = 0.0;
  unsigned long MynPointDomain = TotalnPointDomain; TotalnPointDomain = 0;
  SU2_MPI::Allreduce(&MyMeanK, &MeanK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyMaxK, &MaxK, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MynPointDomain, &TotalnPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Compute the mean ---*/
  MeanK /= su2double(TotalnPointDomain);

  /*--- Compute the standard deviation ---*/
  SigmaK = 0.0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          SigmaK += (fabs(K[iPoint]) - MeanK) * (fabs(K[iPoint]) - MeanK);
        }
      }
    }
  }

#ifdef HAVE_MPI
  su2double MySigmaK = SigmaK; SigmaK = 0.0;
  SU2_MPI::Allreduce(&MySigmaK, &SigmaK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  SigmaK = sqrt(SigmaK/su2double(TotalnPointDomain));

  if ((rank == MASTER_NODE) && (!fea))
    cout << "Max K: " << MaxK << ". Mean K: " << MeanK << ". Standard deviation K: " << SigmaK << "." << endl;

  Point_Critical.clear();

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          if (fabs(K[iPoint]) > MeanK + config->GetRefSharpEdges()*SigmaK) {
            Point_Critical.push_back(iPoint);
          }
        }
      }
    }
  }

  /*--- Variables and buffers needed for MPI ---*/

#ifdef HAVE_MPI
  SU2_MPI::Comm_size(MPI_COMM_WORLD, &nProcessor);
#else
  nProcessor = 1;
#endif

  Buffer_Send_nVertex    = new unsigned long [1];
  Buffer_Receive_nVertex = new unsigned long [nProcessor];

  /*--- Count the total number of critical edge nodes. ---*/

  nLocalVertex = Point_Critical.size();
  Buffer_Send_nVertex[0] = nLocalVertex;

  /*--- Communicate to all processors the total number of critical edge nodes. ---*/

#ifdef HAVE_MPI
  MaxLocalVertex = 0;
  SU2_MPI::Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  MaxLocalVertex = nLocalVertex;
  Buffer_Receive_nVertex[0] = nLocalVertex;
#endif


  /*--- Create and initialize to zero some buffers to hold the coordinates
   of the boundary nodes that are communicated from each partition (all-to-all). ---*/

  Buffer_Send_Coord     = new su2double [MaxLocalVertex*nDim];
  Buffer_Receive_Coord  = new su2double [nProcessor*MaxLocalVertex*nDim];

#ifdef HAVE_MPI
  unsigned long nBuffer               = MaxLocalVertex*nDim;
#endif

  for (iVertex = 0; iVertex < MaxLocalVertex; iVertex++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
    }
  }

  /*--- Retrieve and store the coordinates of the sharp edges boundary nodes on
   the local partition and broadcast them to all partitions. ---*/

  for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
    iPoint = Point_Critical[iVertex];
    for (iDim = 0; iDim < nDim; iDim++)
      Buffer_Send_Coord[iVertex*nDim+iDim] = node[iPoint]->GetCoord(iDim);
  }

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);
#else
  for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Buffer_Receive_Coord[iVertex*nDim+iDim] = Buffer_Send_Coord[iVertex*nDim+iDim];
    }
  }
#endif

  /*--- Loop over all interior mesh nodes on the local partition and compute
   the distances to each of the no-slip boundary nodes in the entire mesh.
   Store the minimum distance to the wall for each interior mesh node. ---*/

  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    Coord = node[iPoint]->GetCoord();

    MinDist = 1E20;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Dist = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Dist += (Coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex+iVertex)*nDim+iDim])*
          (Coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex+iVertex)*nDim+iDim]);
        }
        if (Dist!=0.0) Dist = sqrt(Dist);
        else Dist = 0.0;
        if (Dist < MinDist) MinDist = Dist;
      }
    }
    node[iPoint]->SetSharpEdge_Distance(MinDist);
  }

  /*--- Deallocate Max curvature ---*/
  delete[] K;

  /*--- Deallocate the buffers needed for the MPI communication. ---*/
  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_nVertex;
  delete[] Buffer_Receive_nVertex;

}

void CGeometry::FilterValuesAtElementCG(const vector<su2double> &filter_radius,
                                        const vector<pair<unsigned short,su2double> > &kernels,
                                        const unsigned short search_limit,
                                        const su2double *input_values,
                                        su2double *output_values) const
{
  /*--- Apply a filter to "input_values". The filter is an averaging process over the neighbourhood
  of each element, which is a circle in 2D and a sphere in 3D of radius "filter_radius".
  The filter is characterized by its kernel, i.e. how the weights are computed. Multiple kernels
  can be specified in which case they are applied sequentially (each one being applied to the
  output values of the previous filter. ---*/

  unsigned long iElem, iElem_global, limited_searches = 0;

  /*--- Initialize output values and check if we need to do any more work than that ---*/
  for (iElem=0; iElem<nElem; ++iElem)
    output_values[iElem] = input_values[iElem];

  if ( kernels.empty() ) return;

  /*--- FIRST: Gather the adjacency matrix, element centroids, volumes, and values on every
  processor, this is required because the filter reaches far into adjacent partitions. ---*/

  /*--- Adjacency matrix ---*/
  vector<unsigned long> neighbour_start;
  long *neighbour_idx = NULL;
  GetGlobalElementAdjacencyMatrix(neighbour_start,neighbour_idx);

  /*--- Element centroids and volumes ---*/
  su2double *cg_elem  = new su2double [Global_nElemDomain*nDim],
            *vol_elem = new su2double [Global_nElemDomain];

  /*--- Initialize ---*/
  for(iElem=0; iElem<Global_nElemDomain; ++iElem) {
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      cg_elem[nDim*iElem+iDim] = 0.0;
    vol_elem[iElem] = 0.0;
  }
  /*--- Populate ---*/
  for(iElem=0; iElem<nElem; ++iElem) {
    iElem_global = elem[iElem]->GetGlobalIndex();
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      cg_elem[nDim*iElem_global+iDim] = elem[iElem]->GetCG(iDim);
    vol_elem[iElem_global] = elem[iElem]->GetVolume();
  }
#ifdef HAVE_MPI
  /*--- Share with all processors ---*/
  {
    su2double *buffer = NULL, *tmp = NULL;

    buffer = new su2double [Global_nElemDomain*nDim];
    SU2_MPI::Allreduce(cg_elem,buffer,Global_nElemDomain*nDim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    tmp = cg_elem; cg_elem = buffer; delete [] tmp;

    buffer = new su2double [Global_nElemDomain];
    SU2_MPI::Allreduce(vol_elem,buffer,Global_nElemDomain,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    tmp = vol_elem; vol_elem = buffer; delete [] tmp;
  }

  /*--- Account for the duplication introduced by the halo elements and the
  reduction using MPI_SUM, which is required to maintain differentiabillity. ---*/
  vector<char> halo_detect(Global_nElemDomain);
  {
    vector<char> buffer(Global_nElemDomain,0);
    for(iElem=0; iElem<nElem; ++iElem) buffer[elem[iElem]->GetGlobalIndex()] = 1;
    MPI_Allreduce(buffer.data(),halo_detect.data(),Global_nElemDomain,MPI_CHAR,MPI_SUM,MPI_COMM_WORLD);
  }
  for(iElem=0; iElem<Global_nElemDomain; ++iElem) {
    su2double numRepeat = halo_detect[iElem];
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      cg_elem[nDim*iElem+iDim] /= numRepeat;
    vol_elem[iElem] /= numRepeat;
  }
#endif


  /*--- SECOND: Each processor performs the average for its elements. For each
  element we look for neighbours of neighbours of... until the distance to the
  closest newly found one is greater than the filter radius.  ---*/

  /*--- Inputs of a filter stage, like with CG and volumes, each processor needs to see everything ---*/
  su2double *work_values = new su2double [Global_nElemDomain];
  vector<bool> is_neighbor(Global_nElemDomain,false);

  for (unsigned long iKernel=0; iKernel<kernels.size(); ++iKernel)
  {
    unsigned short kernel_type = kernels[iKernel].first;
    su2double kernel_param = kernels[iKernel].second;
    su2double kernel_radius = filter_radius[iKernel];

    /*--- Synchronize work values ---*/
    /*--- Initialize ---*/
    for(iElem=0; iElem<Global_nElemDomain; ++iElem) work_values[iElem] = 0.0;
    /*--- Populate ---*/
    for(iElem=0; iElem<nElem; ++iElem)
      work_values[elem[iElem]->GetGlobalIndex()] = output_values[iElem];
#ifdef HAVE_MPI
    /*--- Share with all processors ---*/
    {
      su2double *buffer = new su2double [Global_nElemDomain], *tmp = NULL;
      SU2_MPI::Allreduce(work_values,buffer,Global_nElemDomain,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      tmp = work_values; work_values = buffer; delete [] tmp;
    }
    /*--- Account for duplication ---*/
    for(iElem=0; iElem<Global_nElemDomain; ++iElem) {
      su2double numRepeat = halo_detect[iElem];
      work_values[iElem] /= numRepeat;
    }
#endif

    /*--- Filter ---*/
    for (iElem=0; iElem<nElem; ++iElem)
    {
      /*--- Center of the search ---*/
      iElem_global = elem[iElem]->GetGlobalIndex();

      /*--- Find the neighbours of iElem ---*/
      vector<long> neighbours;
      limited_searches += !GetRadialNeighbourhood(iElem_global, SU2_TYPE::GetValue(kernel_radius),
                           search_limit, neighbour_start, neighbour_idx, cg_elem, neighbours, is_neighbor);

      /*--- Apply the kernel ---*/
      su2double weight = 0.0, numerator = 0.0, denominator = 0.0;

      switch ( kernel_type ) {
        /*--- distance-based kernels (weighted averages) ---*/
        case CONSTANT_WEIGHT_FILTER: case CONICAL_WEIGHT_FILTER: case GAUSSIAN_WEIGHT_FILTER:

          for (auto idx : neighbours)
          {
            su2double distance = 0.0;
            for (unsigned short iDim=0; iDim<nDim; ++iDim)
              distance += pow(cg_elem[nDim*iElem_global+iDim]-cg_elem[nDim*idx+iDim],2);
            distance = sqrt(distance);

            switch ( kernel_type ) {
              case CONSTANT_WEIGHT_FILTER: weight = 1.0; break;
              case CONICAL_WEIGHT_FILTER:  weight = kernel_radius-distance; break;
              case GAUSSIAN_WEIGHT_FILTER: weight = exp(-0.5*pow(distance/kernel_param,2)); break;
              default: break;
            }
            weight *= vol_elem[idx];
            numerator   += weight*work_values[idx];
            denominator += weight;
          }
          output_values[iElem] = numerator/denominator;
          break;

        /*--- morphology kernels (image processing) ---*/
        case DILATE_MORPH_FILTER: case ERODE_MORPH_FILTER:

          for (auto idx : neighbours)
          {
            switch ( kernel_type ) {
              case DILATE_MORPH_FILTER: numerator += exp(kernel_param*work_values[idx]); break;
              case ERODE_MORPH_FILTER:  numerator += exp(kernel_param*(1.0-work_values[idx])); break;
              default: break;
            }
            denominator += 1.0;
          }
          output_values[iElem] = log(numerator/denominator)/kernel_param;
          if ( kernel_type==ERODE_MORPH_FILTER ) output_values[iElem] = 1.0-output_values[iElem];
          break;

        default:
          SU2_MPI::Error("Unknown type of filter kernel",CURRENT_FUNCTION);
      }
    }
  }
  limited_searches /= kernels.size();
#ifdef HAVE_MPI
  unsigned long tmp = limited_searches;
  MPI_Reduce(&tmp,&limited_searches,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
#endif
  if (rank==MASTER_NODE && limited_searches>0)
    cout << "Warning: The filter radius was limited for " << limited_searches
         << " elements (" << limited_searches/(0.01*Global_nElemDomain) << "%).\n";

  delete [] neighbour_idx;
  delete [] cg_elem;
  delete [] vol_elem;
  delete [] work_values;
}

void CGeometry::GetGlobalElementAdjacencyMatrix(vector<unsigned long> &neighbour_start,
                                                long *&neighbour_idx) const
{
  if ( neighbour_idx != NULL )
    SU2_MPI::Error("neighbour_idx is expected to be NULL, stopping to avoid a potential memory leak",CURRENT_FUNCTION);

  unsigned long iElem, iElem_global;

  /*--- Determine how much space we need for the adjacency matrix by counting the
  neighbours of each element, i.e. its number of faces---*/
  unsigned short *nFaces_elem = new unsigned short [Global_nElemDomain];

  for(iElem=0; iElem<Global_nElemDomain; ++iElem) nFaces_elem[iElem] = 0;

  for(iElem=0; iElem<nElem; ++iElem) {
    iElem_global = elem[iElem]->GetGlobalIndex();
    nFaces_elem[iElem_global] = elem[iElem]->GetnFaces();
  }
#ifdef HAVE_MPI
  /*--- Share with all processors ---*/
  {
    unsigned short *buffer = new unsigned short [Global_nElemDomain], *tmp = NULL;
    MPI_Allreduce(nFaces_elem,buffer,Global_nElemDomain,MPI_UNSIGNED_SHORT,MPI_MAX,MPI_COMM_WORLD);
    /*--- swap pointers and delete old data to keep the same variable name after reduction ---*/
    tmp = nFaces_elem; nFaces_elem = buffer; delete [] tmp;
  }
#endif

  /*--- Vector with the addresses of the start of the neighbours of a given element.
  This is generated by a cumulative sum of the neighbour count. ---*/
  neighbour_start.resize(Global_nElemDomain+1);

  neighbour_start[0] = 0;
  for(iElem=0; iElem<Global_nElemDomain; ++iElem) {
    neighbour_start[iElem+1] = neighbour_start[iElem]+nFaces_elem[iElem];
  }
  delete [] nFaces_elem;

  /*--- Allocate ---*/
  unsigned long matrix_size = neighbour_start[Global_nElemDomain];
  neighbour_idx = new long [matrix_size];
  /*--- Initialize ---*/
  for(iElem=0; iElem<matrix_size; ++iElem) neighbour_idx[iElem] = -1;
  /*--- Populate ---*/
  for(iElem=0; iElem<nElem; ++iElem)
  {
    iElem_global = elem[iElem]->GetGlobalIndex();
    unsigned long start_pos = neighbour_start[iElem_global];

    for(unsigned short iFace=0; iFace<elem[iElem]->GetnFaces(); ++iFace)
    {
      long neighbour = elem[iElem]->GetNeighbor_Elements(iFace);

      if ( neighbour>=0 ) {
        neighbour_idx[start_pos+iFace] = elem[neighbour]->GetGlobalIndex();
      }
    }
  }
#ifdef HAVE_MPI
  /*--- Share with all processors ---*/
  {
    long *buffer = new long [matrix_size], *tmp = NULL;
    MPI_Allreduce(neighbour_idx,buffer,matrix_size,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
    tmp = neighbour_idx; neighbour_idx = buffer; delete [] tmp;
  }
#endif
}

bool CGeometry::GetRadialNeighbourhood(const unsigned long iElem_global,
                                       const passivedouble radius,
                                       size_t search_limit,
                                       const vector<unsigned long> &neighbour_start,
                                       const long *neighbour_idx,
                                       const su2double *cg_elem,
                                       vector<long> &neighbours,
                                       vector<bool> &is_neighbor) const
{
  /*--- Validate inputs if we are debugging. ---*/
  assert(neighbour_start.size() == Global_nElemDomain+1 &&
         neighbour_idx != nullptr && cg_elem != nullptr &&
         is_neighbor.size() == Global_nElemDomain && "invalid inputs");

  /*--- 0 search_limit means "unlimited" (it will probably
   stop once it gathers the entire domain, probably). ---*/
  if (!search_limit) search_limit = numeric_limits<size_t>::max();

  /*--- Center of the search ---*/
  neighbours.clear();
  neighbours.push_back(iElem_global);
  is_neighbor[iElem_global] = true;

  passivedouble X0[3] = {0.0, 0.0, 0.0};
  for (unsigned short iDim=0; iDim<nDim; ++iDim)
    X0[iDim] = SU2_TYPE::GetValue(cg_elem[nDim*iElem_global+iDim]);

  /*--- Loop stops when "neighbours" stops changing size, or degree reaches limit. ---*/
  bool finished = false;
  for (size_t degree=0, start=0; degree < search_limit && !finished; ++degree)
  {
    /*--- For each element of the last degree added consider its immediate
     neighbours, that are not already neighbours, as candidates. ---*/
    vector<long> candidates;

    for (auto it = neighbours.begin()+start; it!=neighbours.end(); ++it) {
      /*--- scan row of the adjacency matrix of element *it ---*/
      for (auto i = neighbour_start[*it]; i < neighbour_start[(*it)+1]; ++i) {
        auto idx = neighbour_idx[i];
        if (idx>=0) if (!is_neighbor[idx]) {
          candidates.push_back(idx);
          /*--- mark as neighbour for now to avoid duplicate candidates. ---*/
          is_neighbor[idx] = true;
        }
      }
    }
    /*--- update start position to fetch next degree candidates. ---*/
    start = neighbours.size();

    /*--- Add candidates within "radius" of X0, if none qualifies we are "finished". ---*/
    finished = true;
    for (auto idx : candidates)
    {
      /*--- passivedouble as we only need to compare "distance". ---*/
      passivedouble distance = 0.0;
      for (unsigned short iDim=0; iDim<nDim; ++iDim)
        distance += pow(X0[iDim]-SU2_TYPE::GetValue(cg_elem[nDim*idx+iDim]),2);

      if(distance < pow(radius,2)) {
        neighbours.push_back(idx);
        finished = false;
      }
      /*--- not a neighbour in the end. ---*/
      else is_neighbor[idx] = false;
    }
  }
  /*--- Restore the state of the working vector for next call. ---*/
  for(auto idx : neighbours) is_neighbor[idx] = false;

  return finished;
}

void CGeometry::SetElemVolume(CConfig *config)
{
  CElement *elements[4] = {NULL, NULL, NULL, NULL}, *element = NULL;

  /*--- Create a bank of elements to avoid instantiating inside loop ---*/
  if (nDim==2) {
    elements[0] = new CTRIA1();
    elements[1] = new CQUAD4();
  } else {
    elements[0] = new CTETRA1();
    elements[1] = new CPYRAM5();
    elements[2] = new CPRISM6();
    elements[3] = new CHEXA8();
  }

  /*--- Compute and store the volume of each "elem" ---*/
  for (unsigned long iElem=0; iElem<nElem; ++iElem)
  {
    /*--- Get the appropriate type of element ---*/
    switch (elem[iElem]->GetVTK_Type()) {
      case TRIANGLE:      element = elements[0]; break;
      case QUADRILATERAL: element = elements[1]; break;
      case TETRAHEDRON:   element = elements[0]; break;
      case PYRAMID:       element = elements[1]; break;
      case PRISM:         element = elements[2]; break;
      case HEXAHEDRON:    element = elements[3]; break;
      default:
        SU2_MPI::Error("Cannot compute the area/volume of a 1D element.",CURRENT_FUNCTION);
    }
    /*--- Set the nodal coordinates of the element ---*/
    for (unsigned short iNode=0; iNode<elem[iElem]->GetnNodes(); ++iNode) {
      unsigned long node_idx = elem[iElem]->GetNode(iNode);
      for (unsigned short iDim=0; iDim<nDim; ++iDim) {
        su2double coord = node[node_idx]->GetCoord(iDim);
        element->SetRef_Coord(iNode, iDim, coord);
      }
    }
    /*--- Compute ---*/
    if(nDim==2) elem[iElem]->SetVolume(element->ComputeArea());
    else        elem[iElem]->SetVolume(element->ComputeVolume());
  }

  delete elements[0];
  delete elements[1];
  if (nDim==3) {
    delete elements[2];
    delete elements[3];
  }
}

void CGeometry::UpdateBoundaries(CConfig *config){

  unsigned short iMarker;
  unsigned long iElem_Surface, iNode_Surface, Point_Surface;

  for (iMarker = 0; iMarker <config->GetnMarker_All(); iMarker++){
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {

        Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);

        node[Point_Surface]->SetPhysicalBoundary(false);
        node[Point_Surface]->SetSolidBoundary(false);

        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[Point_Surface]->SetPhysicalBoundary(true);

        if (config->GetSolid_Wall(iMarker))
          node[Point_Surface]->SetSolidBoundary(true);
      }
    }
  }

  /*--- Update the normal neighbors ---*/

  FindNormal_Neighbor(config);

  /*--- Compute wall distance ---- */

  ComputeWall_Distance(config);

}

void CGeometry::SetGeometryPlanes(CConfig *config) {

  bool loop_on;
  unsigned short iMarker = 0;
  su2double *Face_Normal = NULL, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
  unsigned long jVertex, iVertex, ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;

  /*--- Compute the total number of points on the near-field ---*/
  nVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      nVertex_Wall += nVertex[iMarker];


  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  Xcoord = new su2double[nVertex_Wall];
  Ycoord = new su2double[nVertex_Wall];
  if (nDim == 3) Zcoord = new su2double[nVertex_Wall];
  FaceArea = new su2double[nVertex_Wall];

  /*--- Copy the boundary information to an array ---*/
  iVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
        Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
        if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
        Face_Normal = vertex[iMarker][iVertex]->GetNormal();
        FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
        iVertex_Wall ++;
      }


  //vector<su2double> XCoordList;
  vector<su2double>::iterator IterXCoordList;

  for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
    XCoordList.push_back(Xcoord[iVertex]);

  sort( XCoordList.begin(), XCoordList.end());
  IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
  XCoordList.resize( IterXCoordList - XCoordList.begin() );

  /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
  Xcoord_plane.resize(XCoordList.size());
  Ycoord_plane.resize(XCoordList.size());
  if (nDim==3) Zcoord_plane.resize(XCoordList.size());
  FaceArea_plane.resize(XCoordList.size());
  Plane_points.resize(XCoordList.size());


  su2double dist_ratio;
  unsigned long iCoord;

  /*--- Distribute the values among the different PhiAngles ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (node[iPoint]->GetDomain()) {
      loop_on = true;
      for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
        dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
        if (dist_ratio >= 0 && dist_ratio <= 1.0) {
          if (dist_ratio <= 0.5) iCoord = ixCoord;
          else iCoord = ixCoord+1;
          Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
          Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
          if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
          FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
          Plane_points[iCoord].push_back(iPoint );
          loop_on = false;
        }
      }
    }
  }

  /*--- Order the arrays in ascending values of y ---*/
  /// TODO: Depending on the size of the arrays, this may not be a good way of sorting them.
  for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
    for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
      for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
        if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
          swap(Xcoord_plane[ixCoord][jVertex], Xcoord_plane[ixCoord][jVertex+1]);
          swap(Ycoord_plane[ixCoord][jVertex], Ycoord_plane[ixCoord][jVertex+1]);
          if (nDim==3) swap(Zcoord_plane[ixCoord][jVertex], Zcoord_plane[ixCoord][jVertex+1]);
          swap(Plane_points[ixCoord][jVertex], Plane_points[ixCoord][jVertex+1]);
          swap(FaceArea_plane[ixCoord][jVertex], FaceArea_plane[ixCoord][jVertex+1]);
        }

  /*--- Delete structures ---*/
  delete[] Xcoord; delete[] Ycoord;
  if (Zcoord != NULL) delete[] Zcoord;
  delete[] FaceArea;
}

void CGeometry::SetRotationalVelocity(CConfig *config, bool print) {

  unsigned long iPoint;
  unsigned short iDim;

  su2double RotVel[3] = {0.0,0.0,0.0}, Distance[3] = {0.0,0.0,0.0},
            Center[3] = {0.0,0.0,0.0}, Omega[3] = {0.0,0.0,0.0};

  /*--- Center of rotation & angular velocity vector from config ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    Center[iDim] = config->GetMotion_Origin(iDim);
    Omega[iDim]  = config->GetRotation_Rate(iDim)/config->GetOmega_Ref();
  }

  su2double L_Ref = config->GetLength_Ref();

  /*--- Print some information to the console ---*/

  if (rank == MASTER_NODE && print) {
    cout << " Rotational origin (x, y, z): ( " << Center[0] << ", " << Center[1];
    cout << ", " << Center[2] << " )\n";
    cout << " Angular velocity about x, y, z axes: ( " << Omega[0] << ", ";
    cout << Omega[1] << ", " << Omega[2] << " ) rad/s" << endl;
  }

  /*--- Loop over all nodes and set the rotational velocity ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Get the coordinates of the current node ---*/

    const su2double* Coord = node[iPoint]->GetCoord();

    /*--- Calculate the non-dim. distance from the rotation center ---*/

    for (iDim = 0; iDim < nDim; iDim++)
      Distance[iDim] = (Coord[iDim]-Center[iDim])/L_Ref;

    /*--- Calculate the angular velocity as omega X r ---*/

    RotVel[0] = Omega[1]*(Distance[2]) - Omega[2]*(Distance[1]);
    RotVel[1] = Omega[2]*(Distance[0]) - Omega[0]*(Distance[2]);
    RotVel[2] = Omega[0]*(Distance[1]) - Omega[1]*(Distance[0]);

    /*--- Store the grid velocity at this node ---*/

    node[iPoint]->SetGridVel(RotVel);

  }

}

void CGeometry::SetShroudVelocity(CConfig *config) {

  unsigned long iPoint, iVertex;
  unsigned short iMarker, iMarkerShroud;
  su2double RotVel[3] = {0.0,0.0,0.0};

  /*--- Loop over all vertex in the shroud marker and set the rotational velocity to 0.0 ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++){
    for(iMarkerShroud=0; iMarkerShroud < config->GetnMarker_Shroud(); iMarkerShroud++){
      if(config->GetMarker_Shroud(iMarkerShroud) == config->GetMarker_All_TagBound(iMarker)){
        for (iVertex = 0; iVertex  < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          node[iPoint]->SetGridVel(RotVel);
        }
      }
    }
  }
}

void CGeometry::SetTranslationalVelocity(CConfig *config, bool print) {

  su2double xDot[3] = {0.0,0.0,0.0};

  /*--- Get the translational velocity vector from config ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    xDot[iDim] = config->GetTranslation_Rate(iDim)/config->GetVelocity_Ref();

  /*--- Print some information to the console ---*/

  if (rank == MASTER_NODE && print) {
    cout << " Non-dim. translational velocity: ("
         << xDot[0] << ", " << xDot[1] << ", " << xDot[2] << ")." << endl;
  }

  /*--- Loop over all nodes and set the translational velocity ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetGridVel(xDot);

}

void CGeometry::SetGridVelocity(CConfig *config, unsigned long iter) {

  /*--- Get timestep and whether to use 1st or 2nd order backward finite differences ---*/

  bool FirstOrder = (config->GetTime_Marching() == DT_STEPPING_1ST);
  bool SecondOrder = (config->GetTime_Marching() == DT_STEPPING_2ND);

  su2double TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the velocity of each node in the volume mesh ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/

    const su2double *Coord_nM1 = node[iPoint]->GetCoord_n1();
    const su2double *Coord_n   = node[iPoint]->GetCoord_n();
    const su2double *Coord_nP1 = node[iPoint]->GetCoord();

    /*--- Compute and store mesh velocity with 1st or 2nd-order approximation ---*/

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {

      su2double GridVel = 0.0;

      if (FirstOrder)
        GridVel = (Coord_nP1[iDim] - Coord_n[iDim]) / TimeStep;

      if (SecondOrder)
        GridVel = (1.5*Coord_nP1[iDim] - 2.0*Coord_n[iDim] + 0.5*Coord_nM1[iDim]) / TimeStep;

      node[iPoint]->SetGridVel(iDim, GridVel);
    }
  }

}

const CCompressedSparsePatternUL& CGeometry::GetSparsePattern(ConnectivityType type, unsigned long fillLvl)
{
  bool fvm = (type == ConnectivityType::FiniteVolume);

  CCompressedSparsePatternUL* pattern = nullptr;

  if (fillLvl == 0)
    pattern = fvm? &finiteVolumeCSRFill0 : &finiteElementCSRFill0;
  else
    pattern = fvm? &finiteVolumeCSRFillN : &finiteElementCSRFillN;

  if (pattern->empty()) {
    *pattern = buildCSRPattern(*this, type, fillLvl);
    pattern->buildDiagPtr();
  }

  return *pattern;
}

const CEdgeToNonZeroMapUL& CGeometry::GetEdgeToSparsePatternMap(void)
{
  if (edgeToCSRMap.empty()) {
    if (finiteVolumeCSRFill0.empty()) {
      finiteVolumeCSRFill0 = buildCSRPattern(*this, ConnectivityType::FiniteVolume, 0ul);
    }
    edgeToCSRMap = mapEdgesToSparsePattern(*this, finiteVolumeCSRFill0);
  }
  return edgeToCSRMap;
}

const CCompressedSparsePatternUL& CGeometry::GetEdgeColoring(void)
{
  if (edgeColoring.empty()) {
    /*--- Create a temporary sparse pattern from the edges. ---*/
    /// TODO: Try to avoid temporary once grid information is made contiguous.
    su2vector<unsigned long> outerPtr(nEdge+1);
    su2vector<unsigned long> innerIdx(nEdge*2);

    for(unsigned long iEdge = 0; iEdge < nEdge; ++iEdge) {
      outerPtr(iEdge) = 2*iEdge;
      innerIdx(iEdge*2+0) = edge[iEdge]->GetNode(0);
      innerIdx(iEdge*2+1) = edge[iEdge]->GetNode(1);
    }
    outerPtr(nEdge) = 2*nEdge;

    CCompressedSparsePatternUL pattern(move(outerPtr), move(innerIdx));

    /*--- Color the edges. ---*/
    edgeColoring = colorSparsePattern(pattern, edgeColorGroupSize);

    if(edgeColoring.empty())
      SU2_MPI::Error("Edge coloring failed.", CURRENT_FUNCTION);
  }
  return edgeColoring;
}

const CCompressedSparsePatternUL& CGeometry::GetElementColoring(void)
{
  if (elemColoring.empty()) {
    /*--- Create a temporary sparse pattern from the elements. ---*/
    /// TODO: Try to avoid temporary once grid information is made contiguous.
    vector<unsigned long> outerPtr(nElem+1);
    vector<unsigned long> innerIdx; innerIdx.reserve(nElem);

    for(unsigned long iElem = 0; iElem < nElem; ++iElem) {
      outerPtr[iElem] = innerIdx.size();

      for(unsigned short iNode = 0; iNode < elem[iElem]->GetnNodes(); ++iNode) {
        innerIdx.push_back(elem[iElem]->GetNode(iNode));
      }
    }
    outerPtr[nElem] = innerIdx.size();

    CCompressedSparsePatternUL pattern(outerPtr, innerIdx);

    /*--- Color the elements. ---*/
    elemColoring = colorSparsePattern(pattern, elemColorGroupSize);

    if(elemColoring.empty())
      SU2_MPI::Error("Element coloring failed.", CURRENT_FUNCTION);
  }
  return elemColoring;
}
