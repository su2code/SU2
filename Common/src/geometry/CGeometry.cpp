/*!
 * \file CGeometry.cpp
 * \brief Implementation of the base geometry class.
 * \author F. Palacios, T. Economon
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

#include <unordered_set>

#include "../../include/geometry/CGeometry.hpp"
#include "../../include/geometry/elements/CElement.hpp"
#include "../../include/parallelization/omp_structure.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"
#include "../../include/toolboxes/ndflattener.hpp"

CGeometry::CGeometry() : size(SU2_MPI::GetSize()), rank(SU2_MPI::GetRank()) {}

CGeometry::~CGeometry() {
  unsigned long iElem, iElem_Bound, iVertex;
  unsigned short iMarker;

  if (elem != nullptr) {
    for (iElem = 0; iElem < nElem; iElem++) delete elem[iElem];
    delete[] elem;
  }

  if (bound != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        delete bound[iMarker][iElem_Bound];
      }
      delete[] bound[iMarker];
    }
    delete[] bound;
  }

  delete nodes;

  delete edges;

  if (vertex != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        delete vertex[iMarker][iVertex];
      }
      delete[] vertex[iMarker];
    }
    delete[] vertex;
  }

  delete[] nElem_Bound;
  delete[] nVertex;
  delete[] Marker_All_SendRecv;
  delete[] Tag_to_Marker;

  delete[] beg_node;
  delete[] end_node;
  delete[] nPointLinear;
  delete[] nPointCumulative;

  if (CustomBoundaryHeatFlux != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete[] CustomBoundaryHeatFlux[iMarker];
    }
    delete[] CustomBoundaryHeatFlux;
  }

  if (CustomBoundaryTemperature != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete[] CustomBoundaryTemperature[iMarker];
    }
    delete[] CustomBoundaryTemperature;
  }

  /*--- Delete structures for MPI point-to-point communication. ---*/

  delete[] bufD_P2PRecv;
  delete[] bufD_P2PSend;

  delete[] bufS_P2PRecv;
  delete[] bufS_P2PSend;

  delete[] req_P2PSend;
  delete[] req_P2PRecv;

  delete[] nPoint_P2PRecv;
  delete[] nPoint_P2PSend;

  delete[] Neighbors_P2PSend;
  delete[] Neighbors_P2PRecv;

  delete[] Local_Point_P2PSend;
  delete[] Local_Point_P2PRecv;

  /*--- Delete structures for MPI periodic communication. ---*/

  delete[] bufD_PeriodicRecv;
  delete[] bufD_PeriodicSend;

  delete[] bufS_PeriodicRecv;
  delete[] bufS_PeriodicSend;

  delete[] req_PeriodicSend;
  delete[] req_PeriodicRecv;

  delete[] nPoint_PeriodicRecv;
  delete[] nPoint_PeriodicSend;

  delete[] Neighbors_PeriodicSend;
  delete[] Neighbors_PeriodicRecv;

  delete[] Local_Point_PeriodicSend;
  delete[] Local_Point_PeriodicRecv;

  delete[] Local_Marker_PeriodicSend;
  delete[] Local_Marker_PeriodicRecv;
}

void CGeometry::PreprocessP2PComms(CGeometry* geometry, CConfig* config) {
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
  unsigned long nVertexS, nVertexR, iVertex, MarkerS, MarkerR;

  int iRank, iSend, iRecv, count;

  /*--- Create some temporary structures for tracking sends/recvs. ---*/

  int* nPoint_Send_All = new int[size + 1];
  nPoint_Send_All[0] = 0;
  int* nPoint_Recv_All = new int[size + 1];
  nPoint_Recv_All[0] = 0;
  int* nPoint_Flag = new int[size];

  for (iRank = 0; iRank < size; iRank++) {
    nPoint_Send_All[iRank] = 0;
    nPoint_Recv_All[iRank] = 0;
    nPoint_Flag[iRank] = -1;
  }
  nPoint_Send_All[size] = 0;
  nPoint_Recv_All[size] = 0;

  /*--- Loop through all of our SEND_RECEIVE markers and track
   our sends with each rank. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      /*--- Get the destination rank and number of points to send. ---*/

      iRank = config->GetMarker_All_SendRecv(iMarker) - 1;
      nVertexS = geometry->nVertex[iMarker];

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nPoint_Flag[iRank] != (int)iMarker)) {
        nPoint_Flag[iRank] = (int)iMarker;
        nPoint_Send_All[iRank + 1] += nVertexS;
      }
    }
  }

  delete[] nPoint_Flag;

  /*--- Communicate the number of points to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nPoint_Send_All[1]), 1, MPI_INT, &(nPoint_Recv_All[1]), 1, MPI_INT, SU2_MPI::GetComm());

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  nP2PSend = 0;
  nP2PRecv = 0;

  for (iRank = 0; iRank < size; iRank++) {
    if ((iRank != rank) && (nPoint_Send_All[iRank + 1] > 0)) nP2PSend++;
    if ((iRank != rank) && (nPoint_Recv_All[iRank + 1] > 0)) nP2PRecv++;

    nPoint_Send_All[iRank + 1] += nPoint_Send_All[iRank];
    nPoint_Recv_All[iRank + 1] += nPoint_Recv_All[iRank];
  }

  /*--- Allocate only as much memory as we need for the P2P neighbors. ---*/

  nPoint_P2PSend = new int[nP2PSend + 1];
  nPoint_P2PSend[0] = 0;
  nPoint_P2PRecv = new int[nP2PRecv + 1];
  nPoint_P2PRecv[0] = 0;

  Neighbors_P2PSend = new int[nP2PSend];
  Neighbors_P2PRecv = new int[nP2PRecv];

  iSend = 0;
  iRecv = 0;
  for (iRank = 0; iRank < size; iRank++) {
    if ((nPoint_Send_All[iRank + 1] > nPoint_Send_All[iRank]) && (iRank != rank)) {
      Neighbors_P2PSend[iSend] = iRank;
      nPoint_P2PSend[iSend + 1] = nPoint_Send_All[iRank + 1];
      iSend++;
    }

    if ((nPoint_Recv_All[iRank + 1] > nPoint_Recv_All[iRank]) && (iRank != rank)) {
      Neighbors_P2PRecv[iRecv] = iRank;
      nPoint_P2PRecv[iRecv + 1] = nPoint_Recv_All[iRank + 1];
      iRecv++;
    }
  }

  /*--- Create a reverse mapping of the message to the rank so that we
   can quickly access the correct data in the buffers when receiving
   messages dynamically. ---*/

  P2PSend2Neighbor.clear();
  for (iSend = 0; iSend < nP2PSend; iSend++) P2PSend2Neighbor[Neighbors_P2PSend[iSend]] = iSend;

  P2PRecv2Neighbor.clear();
  for (iRecv = 0; iRecv < nP2PRecv; iRecv++) P2PRecv2Neighbor[Neighbors_P2PRecv[iRecv]] = iRecv;

  delete[] nPoint_Send_All;
  delete[] nPoint_Recv_All;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  Local_Point_P2PSend = nullptr;
  Local_Point_P2PSend = new unsigned long[nPoint_P2PSend[nP2PSend]];
  for (iSend = 0; iSend < nPoint_P2PSend[nP2PSend]; iSend++) Local_Point_P2PSend[iSend] = 0;

  Local_Point_P2PRecv = nullptr;
  Local_Point_P2PRecv = new unsigned long[nPoint_P2PRecv[nP2PRecv]];
  for (iRecv = 0; iRecv < nPoint_P2PRecv[nP2PRecv]; iRecv++) Local_Point_P2PRecv[iRecv] = 0;

  /*--- We allocate the memory for communicating values in a later step
   once we know the maximum packet size that we need to communicate. This
   memory is deallocated and reallocated automatically in the case that
   the previously allocated memory is not sufficient. ---*/

  bufD_P2PSend = nullptr;
  bufD_P2PRecv = nullptr;

  bufS_P2PSend = nullptr;
  bufS_P2PRecv = nullptr;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nP2PSend > 0) {
    req_P2PSend = new SU2_MPI::Request[nP2PSend];
  }
  if (nP2PRecv > 0) {
    req_P2PRecv = new SU2_MPI::Request[nP2PRecv];
  }

  /*--- Build lists of local index values for send. ---*/

  count = 0;
  for (iSend = 0; iSend < nP2PSend; iSend++) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)) {
        MarkerS = iMarker;
        nVertexS = geometry->nVertex[MarkerS];
        iRank = config->GetMarker_All_SendRecv(MarkerS) - 1;

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
      if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)) {
        MarkerR = iMarker + 1;
        nVertexR = geometry->nVertex[MarkerR];
        iRank = abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;

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

void CGeometry::AllocateP2PComms(unsigned short countPerPoint) {
  /*--- This routine is activated whenever we attempt to perform
   a point-to-point MPI communication with our neighbors but the
   memory buffer allocated is not large enough for the packet size.
   Therefore, we deallocate the previously allocated space and
   reallocate a large enough array. Note that after the first set
   communications, this routine will not need to be called again. ---*/

  if (countPerPoint <= maxCountPerPoint) return;

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Store the larger packet size to the class data. ---*/

    maxCountPerPoint = countPerPoint;

    /*-- Deallocate and reallocate our su2double cummunication memory. ---*/

    delete[] bufD_P2PSend;
    bufD_P2PSend = new su2double[maxCountPerPoint * nPoint_P2PSend[nP2PSend]]();

    delete[] bufD_P2PRecv;
    bufD_P2PRecv = new su2double[maxCountPerPoint * nPoint_P2PRecv[nP2PRecv]]();

    delete[] bufS_P2PSend;
    bufS_P2PSend = new unsigned short[maxCountPerPoint * nPoint_P2PSend[nP2PSend]]();

    delete[] bufS_P2PRecv;
    bufS_P2PRecv = new unsigned short[maxCountPerPoint * nPoint_P2PRecv[nP2PRecv]]();
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CGeometry::PostP2PRecvs(CGeometry* geometry, const CConfig* config, unsigned short commType,
                             unsigned short countPerPoint, bool val_reverse) const {
  /*--- Launch the non-blocking recv's first. Note that we have stored
   the counts and sources, so we can launch these before we even load
   the data and send from the neighbor ranks. ---*/

  SU2_OMP_MASTER
  for (int iRecv = 0; iRecv < nP2PRecv; iRecv++) {
    const auto iMessage = iRecv;

    /*--- In some instances related to the adjoint solver, we need
     to reverse the direction of communications such that the normal
     send nodes become the recv nodes and vice-versa. ---*/

    if (val_reverse) {
      /*--- Compute our location in the buffer using the send data
       structure since we are reversing the comms. ---*/

      auto offset = countPerPoint * nPoint_P2PSend[iRecv];

      /*--- Take advantage of cumulative storage format to get the number
       of elems that we need to recv. Note again that we select the send
       points here as the recv points. ---*/

      auto nPointP2P = nPoint_P2PSend[iRecv + 1] - nPoint_P2PSend[iRecv];

      /*--- Total count can include multiple pieces of data per element. ---*/

      auto count = countPerPoint * nPointP2P;

      /*--- Get the rank from which we receive the message. Note again
       that we use the send rank as the source instead of the recv rank. ---*/

      auto source = Neighbors_P2PSend[iRecv];
      auto tag = source + 1;

      /*--- Post non-blocking recv for this proc. Note that we use the
       send buffer here too. This is important to make sure the arrays
       are the correct size. ---*/

      switch (commType) {
        case COMM_TYPE_DOUBLE:
          SU2_MPI::Irecv(&(bufD_P2PSend[offset]), count, MPI_DOUBLE, source, tag, SU2_MPI::GetComm(),
                         &(req_P2PRecv[iRecv]));
          break;
        case COMM_TYPE_UNSIGNED_SHORT:
          SU2_MPI::Irecv(&(bufS_P2PSend[offset]), count, MPI_UNSIGNED_SHORT, source, tag, SU2_MPI::GetComm(),
                         &(req_P2PRecv[iRecv]));
          break;
        default:
          SU2_MPI::Error("Unrecognized data type for point-to-point MPI comms.", CURRENT_FUNCTION);
          break;
      }

    } else {
      /*--- Compute our location in the recv buffer. ---*/

      auto offset = countPerPoint * nPoint_P2PRecv[iRecv];

      /*--- Take advantage of cumulative storage format to get the number
       of elems that we need to recv. ---*/

      auto nPointP2P = nPoint_P2PRecv[iRecv + 1] - nPoint_P2PRecv[iRecv];

      /*--- Total count can include multiple pieces of data per element. ---*/

      auto count = countPerPoint * nPointP2P;

      /*--- Get the rank from which we receive the message. ---*/

      auto source = Neighbors_P2PRecv[iRecv];
      auto tag = source + 1;

      /*--- Post non-blocking recv for this proc. ---*/

      switch (commType) {
        case COMM_TYPE_DOUBLE:
          SU2_MPI::Irecv(&(bufD_P2PRecv[offset]), count, MPI_DOUBLE, source, tag, SU2_MPI::GetComm(),
                         &(req_P2PRecv[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_SHORT:
          SU2_MPI::Irecv(&(bufS_P2PRecv[offset]), count, MPI_UNSIGNED_SHORT, source, tag, SU2_MPI::GetComm(),
                         &(req_P2PRecv[iMessage]));
          break;
        default:
          SU2_MPI::Error("Unrecognized data type for point-to-point MPI comms.", CURRENT_FUNCTION);
          break;
      }
    }
  }
  END_SU2_OMP_MASTER
}

void CGeometry::PostP2PSends(CGeometry* geometry, const CConfig* config, unsigned short commType,
                             unsigned short countPerPoint, int val_iSend, bool val_reverse) const {
  /*--- Post the non-blocking send as soon as the buffer is loaded. ---*/

  /*--- In some instances related to the adjoint solver, we need
   to reverse the direction of communications such that the normal
   send nodes become the recv nodes and vice-versa. ---*/

  SU2_OMP_MASTER
  if (val_reverse) {
    /*--- Compute our location in the buffer using the recv data
     structure since we are reversing the comms. ---*/

    auto offset = countPerPoint * nPoint_P2PRecv[val_iSend];

    /*--- Take advantage of cumulative storage format to get the number
     of points that we need to send. Note again that we select the recv
     points here as the send points. ---*/

    auto nPointP2P = nPoint_P2PRecv[val_iSend + 1] - nPoint_P2PRecv[val_iSend];

    /*--- Total count can include multiple pieces of data per element. ---*/

    auto count = countPerPoint * nPointP2P;

    /*--- Get the rank to which we send the message. Note again
     that we use the recv rank as the dest instead of the send rank. ---*/

    auto dest = Neighbors_P2PRecv[val_iSend];
    auto tag = rank + 1;

    /*--- Post non-blocking send for this proc. Note that we use the
     send buffer here too. This is important to make sure the arrays
     are the correct size. ---*/

    switch (commType) {
      case COMM_TYPE_DOUBLE:
        SU2_MPI::Isend(&(bufD_P2PRecv[offset]), count, MPI_DOUBLE, dest, tag, SU2_MPI::GetComm(),
                       &(req_P2PSend[val_iSend]));
        break;
      case COMM_TYPE_UNSIGNED_SHORT:
        SU2_MPI::Isend(&(bufS_P2PRecv[offset]), count, MPI_UNSIGNED_SHORT, dest, tag, SU2_MPI::GetComm(),
                       &(req_P2PSend[val_iSend]));
        break;
      default:
        SU2_MPI::Error("Unrecognized data type for point-to-point MPI comms.", CURRENT_FUNCTION);
        break;
    }

  } else {
    /*--- Compute our location in the send buffer. ---*/

    auto offset = countPerPoint * nPoint_P2PSend[val_iSend];

    /*--- Take advantage of cumulative storage format to get the number
     of points that we need to send. ---*/

    auto nPointP2P = nPoint_P2PSend[val_iSend + 1] - nPoint_P2PSend[val_iSend];

    /*--- Total count can include multiple pieces of data per element. ---*/

    auto count = countPerPoint * nPointP2P;

    /*--- Get the rank to which we send the message. ---*/

    auto dest = Neighbors_P2PSend[val_iSend];
    auto tag = rank + 1;

    /*--- Post non-blocking send for this proc. ---*/

    switch (commType) {
      case COMM_TYPE_DOUBLE:
        SU2_MPI::Isend(&(bufD_P2PSend[offset]), count, MPI_DOUBLE, dest, tag, SU2_MPI::GetComm(),
                       &(req_P2PSend[val_iSend]));
        break;
      case COMM_TYPE_UNSIGNED_SHORT:
        SU2_MPI::Isend(&(bufS_P2PSend[offset]), count, MPI_UNSIGNED_SHORT, dest, tag, SU2_MPI::GetComm(),
                       &(req_P2PSend[val_iSend]));
        break;
      default:
        SU2_MPI::Error("Unrecognized data type for point-to-point MPI comms.", CURRENT_FUNCTION);
        break;
    }
  }
  END_SU2_OMP_MASTER
}

void CGeometry::GetCommCountAndType(const CConfig* config, unsigned short commType, unsigned short& COUNT_PER_POINT,
                                    unsigned short& MPI_TYPE) const {
  switch (commType) {
    case COORDINATES:
      COUNT_PER_POINT = nDim;
      MPI_TYPE = COMM_TYPE_DOUBLE;
      break;
    case GRID_VELOCITY:
      COUNT_PER_POINT = nDim;
      MPI_TYPE = COMM_TYPE_DOUBLE;
      break;
    case COORDINATES_OLD:
      if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
        COUNT_PER_POINT = nDim * 2;
      else
        COUNT_PER_POINT = nDim;
      MPI_TYPE = COMM_TYPE_DOUBLE;
      break;
    case MAX_LENGTH:
      COUNT_PER_POINT = 1;
      MPI_TYPE = COMM_TYPE_DOUBLE;
      break;
    case NEIGHBORS:
      COUNT_PER_POINT = 1;
      MPI_TYPE = COMM_TYPE_UNSIGNED_SHORT;
      break;
    default:
      SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.", CURRENT_FUNCTION);
      break;
  }
}

void CGeometry::InitiateComms(CGeometry* geometry, const CConfig* config, unsigned short commType) const {
  if (nP2PSend == 0) return;

  /*--- Local variables ---*/

  unsigned short iDim;
  unsigned short COUNT_PER_POINT = 0;
  unsigned short MPI_TYPE = 0;

  unsigned long iPoint, msg_offset, buf_offset;

  int iMessage, iSend, nSend;

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  GetCommCountAndType(config, commType, COUNT_PER_POINT, MPI_TYPE);

  /*--- Check to make sure we have created a large enough buffer
   for these comms during preprocessing. This is only for the su2double
   buffer. It will be reallocated whenever we find a larger count
   per point. After the first cycle of comms, this should be inactive. ---*/

  geometry->AllocateP2PComms(COUNT_PER_POINT);

  /*--- Set some local pointers to make access simpler. ---*/

  su2double* bufDSend = geometry->bufD_P2PSend;
  unsigned short* bufSSend = geometry->bufS_P2PSend;

  su2double* vector = nullptr;

  /*--- Load the specified quantity from the solver into the generic
   communication buffer in the geometry class. ---*/

  /*--- Post all non-blocking recvs first before sends. ---*/

  geometry->PostP2PRecvs(geometry, config, MPI_TYPE, COUNT_PER_POINT, false);

  for (iMessage = 0; iMessage < nP2PSend; iMessage++) {
    /*--- Get the offset in the buffer for the start of this message. ---*/

    msg_offset = nPoint_P2PSend[iMessage];

    /*--- Total count can include multiple pieces of data per element. ---*/

    nSend = (nPoint_P2PSend[iMessage + 1] - nPoint_P2PSend[iMessage]);

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (iSend = 0; iSend < nSend; iSend++) {
      /*--- Get the local index for this communicated data. ---*/

      iPoint = geometry->Local_Point_P2PSend[msg_offset + iSend];

      /*--- Compute the offset in the recv buffer for this point. ---*/

      buf_offset = (msg_offset + iSend) * COUNT_PER_POINT;

      switch (commType) {
        case COORDINATES:
          vector = nodes->GetCoord(iPoint);
          for (iDim = 0; iDim < nDim; iDim++) bufDSend[buf_offset + iDim] = vector[iDim];
          break;
        case GRID_VELOCITY:
          vector = nodes->GetGridVel(iPoint);
          for (iDim = 0; iDim < nDim; iDim++) bufDSend[buf_offset + iDim] = vector[iDim];
          break;
        case COORDINATES_OLD:
          vector = nodes->GetCoord_n(iPoint);
          for (iDim = 0; iDim < nDim; iDim++) {
            bufDSend[buf_offset + iDim] = vector[iDim];
          }
          if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) {
            vector = nodes->GetCoord_n1(iPoint);
            for (iDim = 0; iDim < nDim; iDim++) {
              bufDSend[buf_offset + nDim + iDim] = vector[iDim];
            }
          }
          break;
        case MAX_LENGTH:
          bufDSend[buf_offset] = nodes->GetMaxLength(iPoint);
          break;
        case NEIGHBORS:
          bufSSend[buf_offset] = geometry->nodes->GetnNeighbor(iPoint);
          break;
        default:
          SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.", CURRENT_FUNCTION);
          break;
      }
    }
    END_SU2_OMP_FOR

    /*--- Launch the point-to-point MPI send for this message. ---*/

    geometry->PostP2PSends(geometry, config, MPI_TYPE, COUNT_PER_POINT, iMessage, false);
  }
}

void CGeometry::CompleteComms(CGeometry* geometry, const CConfig* config, unsigned short commType) {
  if (nP2PRecv == 0) return;

  /*--- Local variables ---*/

  unsigned short iDim, COUNT_PER_POINT = 0, MPI_TYPE = 0;
  unsigned long iPoint, iRecv, nRecv, msg_offset, buf_offset;

  int ind, source, iMessage, jRecv;

  /*--- Global status so all threads can see the result of Waitany. ---*/
  static SU2_MPI::Status status;

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  GetCommCountAndType(config, commType, COUNT_PER_POINT, MPI_TYPE);

  /*--- Set some local pointers to make access simpler. ---*/

  const su2double* bufDRecv = geometry->bufD_P2PRecv;
  const unsigned short* bufSRecv = geometry->bufS_P2PRecv;

  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. Note that we
   recv and store the data in any order to take advantage of the
   non-blocking comms. ---*/

  for (iMessage = 0; iMessage < nP2PRecv; iMessage++) {
    /*--- For efficiency, recv the messages dynamically based on
     the order they arrive. ---*/

    SU2_OMP_SAFE_GLOBAL_ACCESS(SU2_MPI::Waitany(nP2PRecv, req_P2PRecv, &ind, &status);)

    /*--- Once we have recv'd a message, get the source rank. ---*/

    source = status.MPI_SOURCE;

    /*--- We know the offsets based on the source rank. ---*/

    jRecv = P2PRecv2Neighbor[source];

    /*--- Get the offset in the buffer for the start of this message. ---*/

    msg_offset = nPoint_P2PRecv[jRecv];

    /*--- Get the number of packets to be received in this message. ---*/

    nRecv = nPoint_P2PRecv[jRecv + 1] - nPoint_P2PRecv[jRecv];

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (iRecv = 0; iRecv < nRecv; iRecv++) {
      /*--- Get the local index for this communicated data. ---*/

      iPoint = geometry->Local_Point_P2PRecv[msg_offset + iRecv];

      /*--- Compute the total offset in the recv buffer for this point. ---*/

      buf_offset = (msg_offset + iRecv) * COUNT_PER_POINT;

      /*--- Store the data correctly depending on the quantity. ---*/

      switch (commType) {
        case COORDINATES:
          for (iDim = 0; iDim < nDim; iDim++) nodes->SetCoord(iPoint, iDim, bufDRecv[buf_offset + iDim]);
          break;
        case GRID_VELOCITY:
          for (iDim = 0; iDim < nDim; iDim++) nodes->SetGridVel(iPoint, iDim, bufDRecv[buf_offset + iDim]);
          break;
        case COORDINATES_OLD:
          nodes->SetCoord_n(iPoint, &bufDRecv[buf_offset]);
          if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
            nodes->SetCoord_n1(iPoint, &bufDRecv[buf_offset + nDim]);
          break;
        case MAX_LENGTH:
          nodes->SetMaxLength(iPoint, bufDRecv[buf_offset]);
          break;
        case NEIGHBORS:
          nodes->SetnNeighbor(iPoint, bufSRecv[buf_offset]);
          break;
        default:
          SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.", CURRENT_FUNCTION);
          break;
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- Verify that all non-blocking point-to-point sends have finished.
   Note that this should be satisfied, as we have received all of the
   data in the loop above at this point. ---*/

#ifdef HAVE_MPI
  SU2_OMP_SAFE_GLOBAL_ACCESS(SU2_MPI::Waitall(nP2PSend, req_P2PSend, MPI_STATUS_IGNORE);)
#endif
}

void CGeometry::PreprocessPeriodicComms(CGeometry* geometry, CConfig* config) {
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

  int* nPoint_Send_All = new int[size + 1];
  nPoint_Send_All[0] = 0;
  int* nPoint_Recv_All = new int[size + 1];
  nPoint_Recv_All[0] = 0;
  int* nPoint_Flag = new int[size];

  for (iRank = 0; iRank < size; iRank++) {
    nPoint_Send_All[iRank] = 0;
    nPoint_Recv_All[iRank] = 0;
    nPoint_Flag[iRank] = -1;
  }
  nPoint_Send_All[size] = 0;
  nPoint_Recv_All[size] = 0;

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

        if (geometry->nodes->GetDomain(iPoint)) {
          /*--- Get the rank that holds the matching periodic point
           on the other marker in the periodic pair. ---*/

          iRank = (int)geometry->vertex[iMarker][iVertex]->GetDonorProcessor();

          /*--- If we have not visited this point last, increment our
           number of points that must be sent to a particular proc. ---*/

          if ((nPoint_Flag[iRank] != (int)iPoint)) {
            nPoint_Flag[iRank] = (int)iPoint;
            nPoint_Send_All[iRank + 1] += 1;
          }
        }
      }
    }
  }

  delete[] nPoint_Flag;

  /*--- Communicate the number of points to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many periodic points it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nPoint_Send_All[1]), 1, MPI_INT, &(nPoint_Recv_All[1]), 1, MPI_INT, SU2_MPI::GetComm());

  /*--- Check how many messages we will be sending and receiving.
   Here we also put the counters into cumulative storage format to
   make the communications simpler. Note that we are allowing each
   rank to communicate to themselves in these counters, although
   it will not be done through MPI. ---*/

  nPeriodicSend = 0;
  nPeriodicRecv = 0;

  for (iRank = 0; iRank < size; iRank++) {
    if ((nPoint_Send_All[iRank + 1] > 0)) nPeriodicSend++;
    if ((nPoint_Recv_All[iRank + 1] > 0)) nPeriodicRecv++;

    nPoint_Send_All[iRank + 1] += nPoint_Send_All[iRank];
    nPoint_Recv_All[iRank + 1] += nPoint_Recv_All[iRank];
  }

  /*--- Allocate only as much memory as needed for the periodic neighbors. ---*/

  nPoint_PeriodicSend = new int[nPeriodicSend + 1];
  nPoint_PeriodicSend[0] = 0;
  nPoint_PeriodicRecv = new int[nPeriodicRecv + 1];
  nPoint_PeriodicRecv[0] = 0;

  Neighbors_PeriodicSend = new int[nPeriodicSend];
  Neighbors_PeriodicRecv = new int[nPeriodicRecv];

  iSend = 0;
  iRecv = 0;
  for (iRank = 0; iRank < size; iRank++) {
    if ((nPoint_Send_All[iRank + 1] > nPoint_Send_All[iRank])) {
      Neighbors_PeriodicSend[iSend] = iRank;
      nPoint_PeriodicSend[iSend + 1] = nPoint_Send_All[iRank + 1];
      iSend++;
    }
    if ((nPoint_Recv_All[iRank + 1] > nPoint_Recv_All[iRank])) {
      Neighbors_PeriodicRecv[iRecv] = iRank;
      nPoint_PeriodicRecv[iRecv + 1] = nPoint_Recv_All[iRank + 1];
      iRecv++;
    }
  }

  /*--- Create a reverse mapping of the message to the rank so that we
   can quickly access the correct data in the buffers when receiving
   messages dynamically later during the iterations. ---*/

  PeriodicSend2Neighbor.clear();
  for (iSend = 0; iSend < nPeriodicSend; iSend++) PeriodicSend2Neighbor[Neighbors_PeriodicSend[iSend]] = iSend;

  PeriodicRecv2Neighbor.clear();
  for (iRecv = 0; iRecv < nPeriodicRecv; iRecv++) PeriodicRecv2Neighbor[Neighbors_PeriodicRecv[iRecv]] = iRecv;

  delete[] nPoint_Send_All;
  delete[] nPoint_Recv_All;

  /*--- Allocate the memory to store the local index values for both
   the send and receive periodic points and periodic index. ---*/

  Local_Point_PeriodicSend = nullptr;
  Local_Point_PeriodicSend = new unsigned long[nPoint_PeriodicSend[nPeriodicSend]];
  for (iSend = 0; iSend < nPoint_PeriodicSend[nPeriodicSend]; iSend++) Local_Point_PeriodicSend[iSend] = 0;

  Local_Marker_PeriodicSend = nullptr;
  Local_Marker_PeriodicSend = new unsigned long[nPoint_PeriodicSend[nPeriodicSend]];
  for (iSend = 0; iSend < nPoint_PeriodicSend[nPeriodicSend]; iSend++) Local_Marker_PeriodicSend[iSend] = 0;

  Local_Point_PeriodicRecv = nullptr;
  Local_Point_PeriodicRecv = new unsigned long[nPoint_PeriodicRecv[nPeriodicRecv]];
  for (iRecv = 0; iRecv < nPoint_PeriodicRecv[nPeriodicRecv]; iRecv++) Local_Point_PeriodicRecv[iRecv] = 0;

  Local_Marker_PeriodicRecv = nullptr;
  Local_Marker_PeriodicRecv = new unsigned long[nPoint_PeriodicRecv[nPeriodicRecv]];
  for (iRecv = 0; iRecv < nPoint_PeriodicRecv[nPeriodicRecv]; iRecv++) Local_Marker_PeriodicRecv[iRecv] = 0;

  /*--- We allocate the buffers for communicating values in a later step
   once we know the maximum packet size that we need to communicate. This
   memory is deallocated and reallocated automatically in the case that
   the previously allocated memory is not sufficient. ---*/

  bufD_PeriodicSend = nullptr;
  bufD_PeriodicRecv = nullptr;

  bufS_PeriodicSend = nullptr;
  bufS_PeriodicRecv = nullptr;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nPeriodicSend > 0) {
    req_PeriodicSend = new SU2_MPI::Request[nPeriodicSend];
  }
  if (nPeriodicRecv > 0) {
    req_PeriodicRecv = new SU2_MPI::Request[nPeriodicRecv];
  }

  /*--- Allocate arrays for sending the periodic point index and marker
   index to the recv rank so that it can store the local values. Therefore,
   the recv rank can quickly loop through the buffers to unpack the data. ---*/

  unsigned short nPackets = 2;
  auto* idSend = new unsigned long[nPoint_PeriodicSend[nPeriodicSend] * nPackets];
  for (iSend = 0; iSend < nPoint_PeriodicSend[nPeriodicSend] * nPackets; iSend++) idSend[iSend] = 0;

  /*--- Build the lists of local index and periodic marker index values. ---*/

  ii = 0;
  jj = 0;
  for (iSend = 0; iSend < nPeriodicSend; iSend++) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          /*--- Get the current periodic point index. We only communicate
           the owned nodes on a rank, as the MPI comms will take care of
           the halos after completing the periodic comms. ---*/

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->nodes->GetDomain(iPoint)) {
            /*--- Get the rank that holds the matching periodic point
             on the other marker in the periodic pair. ---*/

            iRank = (int)geometry->vertex[iMarker][iVertex]->GetDonorProcessor();

            /*--- If the rank for the current periodic point matches the
             rank of the current send message, then store the local point
             index on the matching periodic point and the periodic marker
             index to be communicated to the recv rank. ---*/

            if (iRank == Neighbors_PeriodicSend[iSend]) {
              Local_Point_PeriodicSend[ii] = iPoint;
              Local_Marker_PeriodicSend[ii] = (unsigned long)iMarker;
              jj = ii * nPackets;
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

  auto* idRecv = new unsigned long[nPoint_PeriodicRecv[nPeriodicRecv] * nPackets];
  for (iRecv = 0; iRecv < nPoint_PeriodicRecv[nPeriodicRecv] * nPackets; iRecv++) idRecv[iRecv] = 0;

#ifdef HAVE_MPI

  int iMessage, offset, count, source, dest, tag;

  /*--- Launch the non-blocking recv's first. Note that we have stored
   the counts and sources, so we can launch these before we even load
   the data and send from the periodically matching ranks. ---*/

  iMessage = 0;
  for (iRecv = 0; iRecv < nPeriodicRecv; iRecv++) {
    /*--- Compute our location in the recv buffer. ---*/

    offset = nPackets * nPoint_PeriodicRecv[iRecv];

    /*--- Take advantage of cumulative storage format to get the number
     of elems that we need to recv. ---*/

    count = nPackets * (nPoint_PeriodicRecv[iRecv + 1] - nPoint_PeriodicRecv[iRecv]);

    /*--- Get the rank from which we receive the message. ---*/

    source = Neighbors_PeriodicRecv[iRecv];
    tag = source + 1;

    /*--- Post non-blocking recv for this proc. ---*/

    SU2_MPI::Irecv(&(static_cast<unsigned long*>(idRecv)[offset]), count, MPI_UNSIGNED_LONG, source, tag,
                   SU2_MPI::GetComm(), &(req_PeriodicRecv[iMessage]));

    /*--- Increment message counter. ---*/

    iMessage++;
  }

  /*--- Post the non-blocking sends. ---*/

  iMessage = 0;
  for (iSend = 0; iSend < nPeriodicSend; iSend++) {
    /*--- Compute our location in the send buffer. ---*/

    offset = nPackets * nPoint_PeriodicSend[iSend];

    /*--- Take advantage of cumulative storage format to get the number
     of points that we need to send. ---*/

    count = nPackets * (nPoint_PeriodicSend[iSend + 1] - nPoint_PeriodicSend[iSend]);

    /*--- Get the rank to which we send the message. ---*/

    dest = Neighbors_PeriodicSend[iSend];
    tag = rank + 1;

    /*--- Post non-blocking send for this proc. ---*/

    SU2_MPI::Isend(&(static_cast<unsigned long*>(idSend)[offset]), count, MPI_UNSIGNED_LONG, dest, tag,
                   SU2_MPI::GetComm(), &(req_PeriodicSend[iMessage]));

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
    iRank = geometry->PeriodicRecv2Neighbor[rank];
    iRecv = geometry->nPoint_PeriodicRecv[iRank] * nPackets;
    myStart = nPoint_PeriodicSend[val_iSend] * nPackets;
    myFinal = nPoint_PeriodicSend[val_iSend + 1] * nPackets;
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
    Local_Point_PeriodicRecv[iRecv] = idRecv[ii];
    ii++;
    Local_Marker_PeriodicRecv[iRecv] = idRecv[ii];
    ii++;
  }

  delete[] idSend;
  delete[] idRecv;
}

void CGeometry::AllocatePeriodicComms(unsigned short countPerPeriodicPoint) {
  /*--- This routine is activated whenever we attempt to perform
   a periodic MPI communication with our neighbors but the
   memory buffer allocated is not large enough for the packet size.
   Therefore, we deallocate the previously allocated arrays and
   reallocate a large enough array. Note that after the first set
   communications, this routine will not need to be called again. ---*/

  if (countPerPeriodicPoint <= maxCountPerPeriodicPoint) return;

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Store the larger packet size to the class data. ---*/

    maxCountPerPeriodicPoint = countPerPeriodicPoint;

    /*--- Store the total size of the send/recv arrays for clarity. ---*/

    auto nSend = countPerPeriodicPoint * nPoint_PeriodicSend[nPeriodicSend];
    auto nRecv = countPerPeriodicPoint * nPoint_PeriodicRecv[nPeriodicRecv];

    /*-- Deallocate and reallocate our cummunication memory. ---*/

    delete[] bufD_PeriodicSend;
    bufD_PeriodicSend = new su2double[nSend]();

    delete[] bufD_PeriodicRecv;
    bufD_PeriodicRecv = new su2double[nRecv]();

    delete[] bufS_PeriodicSend;
    bufS_PeriodicSend = new unsigned short[nSend]();

    delete[] bufS_PeriodicRecv;
    bufS_PeriodicRecv = new unsigned short[nRecv]();
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CGeometry::PostPeriodicRecvs(CGeometry* geometry, const CConfig* config, unsigned short commType,
                                  unsigned short countPerPeriodicPoint) {
  /*--- In parallel, communicate the data with non-blocking send/recv. ---*/

#ifdef HAVE_MPI

  /*--- Launch the non-blocking recv's first. Note that we have stored
   the counts and sources, so we can launch these before we even load
   the data and send from the neighbor ranks. ---*/

  SU2_OMP_MASTER
  for (int iRecv = 0; iRecv < nPeriodicRecv; iRecv++) {
    /*--- Compute our location in the recv buffer. ---*/

    auto offset = countPerPeriodicPoint * nPoint_PeriodicRecv[iRecv];

    /*--- Take advantage of cumulative storage format to get the number
     of elems that we need to recv. ---*/

    auto nPointPeriodic = nPoint_PeriodicRecv[iRecv + 1] - nPoint_PeriodicRecv[iRecv];

    /*--- Total count can include multiple pieces of data per element. ---*/

    auto count = countPerPeriodicPoint * nPointPeriodic;

    /*--- Get the rank from which we receive the message. ---*/

    auto source = Neighbors_PeriodicRecv[iRecv];
    auto tag = source + 1;

    /*--- Post non-blocking recv for this proc. ---*/

    switch (commType) {
      case COMM_TYPE_DOUBLE:
        SU2_MPI::Irecv(&(static_cast<su2double*>(bufD_PeriodicRecv)[offset]), count, MPI_DOUBLE, source, tag,
                       SU2_MPI::GetComm(), &(req_PeriodicRecv[iRecv]));
        break;
      case COMM_TYPE_UNSIGNED_SHORT:
        SU2_MPI::Irecv(&(static_cast<unsigned short*>(bufS_PeriodicRecv)[offset]), count, MPI_UNSIGNED_SHORT, source,
                       tag, SU2_MPI::GetComm(), &(req_PeriodicRecv[iRecv]));
        break;
      default:
        SU2_MPI::Error("Unrecognized data type for periodic MPI comms.", CURRENT_FUNCTION);
        break;
    }
  }
  END_SU2_OMP_MASTER

#endif
}

void CGeometry::PostPeriodicSends(CGeometry* geometry, const CConfig* config, unsigned short commType,
                                  unsigned short countPerPeriodicPoint, int val_iSend) const {
  /*--- In parallel, communicate the data with non-blocking send/recv. ---*/

#ifdef HAVE_MPI
  SU2_OMP_MASTER {
    /*--- Post the non-blocking send as soon as the buffer is loaded. ---*/

    /*--- Compute our location in the send buffer. ---*/

    auto offset = countPerPeriodicPoint * nPoint_PeriodicSend[val_iSend];

    /*--- Take advantage of cumulative storage format to get the number
     of points that we need to send. ---*/

    auto nPointPeriodic = (nPoint_PeriodicSend[val_iSend + 1] - nPoint_PeriodicSend[val_iSend]);

    /*--- Total count can include multiple pieces of data per element. ---*/

    auto count = countPerPeriodicPoint * nPointPeriodic;

    /*--- Get the rank to which we send the message. ---*/

    auto dest = Neighbors_PeriodicSend[val_iSend];
    auto tag = rank + 1;

    /*--- Post non-blocking send for this proc. ---*/

    switch (commType) {
      case COMM_TYPE_DOUBLE:
        SU2_MPI::Isend(&(static_cast<su2double*>(bufD_PeriodicSend)[offset]), count, MPI_DOUBLE, dest, tag,
                       SU2_MPI::GetComm(), &(req_PeriodicSend[val_iSend]));
        break;
      case COMM_TYPE_UNSIGNED_SHORT:
        SU2_MPI::Isend(&(static_cast<unsigned short*>(bufS_PeriodicSend)[offset]), count, MPI_UNSIGNED_SHORT, dest, tag,
                       SU2_MPI::GetComm(), &(req_PeriodicSend[val_iSend]));
        break;
      default:
        SU2_MPI::Error("Unrecognized data type for periodic MPI comms.", CURRENT_FUNCTION);
        break;
    }
  }
  END_SU2_OMP_MASTER
#else

  /*--- Copy my own rank's data into the recv buffer directly in serial. ---*/

  int myStart, myFinal, iRecv, iRank;
  iRank = geometry->PeriodicRecv2Neighbor[rank];
  iRecv = geometry->nPoint_PeriodicRecv[iRank] * countPerPeriodicPoint;
  myStart = nPoint_PeriodicSend[val_iSend] * countPerPeriodicPoint;
  myFinal = nPoint_PeriodicSend[val_iSend + 1] * countPerPeriodicPoint;

  switch (commType) {
    case COMM_TYPE_DOUBLE:
      parallelCopy(myFinal - myStart, &bufD_PeriodicSend[myStart], &bufD_PeriodicRecv[iRecv]);
      break;
    case COMM_TYPE_UNSIGNED_SHORT:
      parallelCopy(myFinal - myStart, &bufS_PeriodicSend[myStart], &bufS_PeriodicRecv[iRecv]);
      break;
    default:
      SU2_MPI::Error("Unrecognized data type for periodic MPI comms.", CURRENT_FUNCTION);
      break;
  }

#endif
}

void CGeometry::SetEdges() {
  nEdge = 0;
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    for (auto iNode = 0u; iNode < nodes->GetnPoint(iPoint); iNode++) {
      auto jPoint = nodes->GetPoint(iPoint, iNode);
      for (auto jNode = 0u; jNode < nodes->GetnPoint(jPoint); jNode++) {
        if (nodes->GetPoint(jPoint, jNode) == iPoint) {
          auto TestEdge = nodes->GetEdge(jPoint, jNode);
          if (TestEdge == -1) {
            nodes->SetEdge(iPoint, nEdge, iNode);
            nodes->SetEdge(jPoint, nEdge, jNode);
            nEdge++;
          }
          break;
        }
      }
    }
  }

  edges = new CEdge(nEdge, nDim);

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    for (auto jPoint : nodes->GetPoints(iPoint)) {
      if (iPoint < jPoint) {
        auto iEdge = FindEdge(iPoint, jPoint);
        edges->SetNodes(iEdge, iPoint, jPoint);
      }
    }
  }
  edges->SetPaddingNodes();
}

void CGeometry::SetFaces() {
  //  unsigned long iPoint, jPoint, iFace;
  //  unsigned short jNode, iNode;
  //  long TestFace = 0;
  //
  //  nFace = 0;
  //  for (iPoint = 0; iPoint < nPoint; iPoint++)
  //    for (iNode = 0; iNode < nodes->GetnPoint(iPoint); iNode++) {
  //      jPoint = nodes->GetPoint(iPoint, iNode);
  //      for (jNode = 0; jNode < nodes->GetnPoint(jPoint); jNode++)
  //        if (nodes->GetPoint(jPoint, jNode) == iPoint) {
  //          TestFace = nodes->GetFace(jPoint, jNode);
  //          break;
  //        }
  //      if (TestFace == -1) {
  //        nodes->SetFace(iPoint, nFace, iNode);
  //        nodes->SetFace(jPoint, nFace, jNode);
  //        nFace++;
  //      }
  //    }
  //
  //  face = new CFace*[nFace];
  //
  //  for (iPoint = 0; iPoint < nPoint; iPoint++)
  //    for (iNode = 0; iNode < nodes->GetnPoint(iPoint); iNode++) {
  //      jPoint = nodes->GetPoint(iPoint, iNode);
  //      iFace = FindFace(iPoint, jPoint);
  //      if (iPoint < jPoint) face[iFace] = new CFace(iPoint, jPoint, nDim);
  //    }
}

void CGeometry::TestGeometry() const {
  ofstream para_file;

  para_file.open("test_geometry.dat", ios::out);

  auto* Normal = new su2double[nDim];

  for (unsigned long iEdge = 0; iEdge < nEdge; iEdge++) {
    para_file << "Edge index: " << iEdge << endl;
    para_file << "   Point index: " << edges->GetNode(iEdge, 0) << "\t" << edges->GetNode(iEdge, 1) << endl;
    edges->GetNormal(iEdge, Normal);
    para_file << "      Face normal : ";
    for (unsigned short iDim = 0; iDim < nDim; iDim++) para_file << Normal[iDim] << "\t";
    para_file << endl;
  }

  para_file << endl;
  para_file << endl;
  para_file << endl;
  para_file << endl;

  for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++) {
    para_file << "Marker index: " << iMarker << endl;
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      para_file << "   Vertex index: " << iVertex << endl;
      para_file << "      Point index: " << vertex[iMarker][iVertex]->GetNode() << endl;
      para_file << "      Point coordinates : ";
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        para_file << nodes->GetCoord(vertex[iMarker][iVertex]->GetNode(), iDim) << "\t";
      }
      para_file << endl;
      vertex[iMarker][iVertex]->GetNormal(Normal);
      para_file << "         Face normal : ";
      for (unsigned short iDim = 0; iDim < nDim; iDim++) para_file << Normal[iDim] << "\t";
      para_file << endl;
    }
  }

  delete[] Normal;
}

bool CGeometry::SegmentIntersectsPlane(const su2double* Segment_P0, const su2double* Segment_P1, su2double Variable_P0,
                                       su2double Variable_P1, const su2double* Plane_P0, const su2double* Plane_Normal,
                                       su2double* Intersection, su2double& Variable_Interp) {
  su2double u[3], v[3], Denominator, Numerator, Aux, ModU;
  su2double epsilon =
      1E-6;  // An epsilon is added to eliminate, as much as possible, the posibility of a line that intersects a point
  unsigned short iDim;

  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P1[iDim] - Segment_P0[iDim];
    v[iDim] = (Plane_P0[iDim] + epsilon) - Segment_P0[iDim];
  }

  ModU = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

  Numerator =
      (Plane_Normal[0] + epsilon) * v[0] + (Plane_Normal[1] + epsilon) * v[1] + (Plane_Normal[2] + epsilon) * v[2];
  Denominator =
      (Plane_Normal[0] + epsilon) * u[0] + (Plane_Normal[1] + epsilon) * u[1] + (Plane_Normal[2] + epsilon) * u[2];

  if (fabs(Denominator) <= 0.0) return (false);  // No intersection.

  Aux = Numerator / Denominator;

  if (Aux < 0.0 || Aux > 1.0) return (false);  // No intersection.

  for (iDim = 0; iDim < 3; iDim++) Intersection[iDim] = Segment_P0[iDim] + Aux * u[iDim];

  /*--- Check that the intersection is in the segment ---*/

  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P0[iDim] - Intersection[iDim];
    v[iDim] = Segment_P1[iDim] - Intersection[iDim];
  }

  Variable_Interp = Variable_P0 + (Variable_P1 - Variable_P0) * sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) / ModU;

  Denominator =
      (Plane_Normal[0] + epsilon) * u[0] + (Plane_Normal[1] + epsilon) * u[1] + (Plane_Normal[2] + epsilon) * u[2];
  Numerator =
      (Plane_Normal[0] + epsilon) * v[0] + (Plane_Normal[1] + epsilon) * v[1] + (Plane_Normal[2] + epsilon) * v[2];

  Aux = Numerator * Denominator;

  if (Aux > 0.0) return (false);  // Intersection outside the segment.

  return (true);
}

bool CGeometry::RayIntersectsTriangle(const su2double orig[3], const su2double dir[3], const su2double vert0[3],
                                      const su2double vert1[3], const su2double vert2[3], su2double* intersect) {
  const passivedouble epsilon = 0.000001;
  su2double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  su2double det, inv_det, t, u, v;

  /*--- Find vectors for two edges sharing vert0 ---*/

  GeometryToolbox::Distance(3, vert1, vert0, edge1);
  GeometryToolbox::Distance(3, vert2, vert0, edge2);

  /*--- Begin calculating determinant - also used to calculate U parameter ---*/

  GeometryToolbox::CrossProduct(dir, edge2, pvec);

  /*--- If determinant is near zero, ray lies in plane of triangle ---*/

  det = GeometryToolbox::DotProduct(3, edge1, pvec);

  if (fabs(det) < epsilon) return (false);

  inv_det = 1.0 / det;

  /*--- Calculate distance from vert0 to ray origin ---*/

  GeometryToolbox::Distance(3, orig, vert0, tvec);

  /*--- Calculate U parameter and test bounds ---*/

  u = inv_det * GeometryToolbox::DotProduct(3, tvec, pvec);

  if (u < 0.0 || u > 1.0) return (false);

  /*--- prepare to test V parameter ---*/

  GeometryToolbox::CrossProduct(tvec, edge1, qvec);

  /*--- Calculate V parameter and test bounds ---*/

  v = inv_det * GeometryToolbox::DotProduct(3, dir, qvec);

  if (v < 0.0 || u + v > 1.0) return (false);

  /*--- Calculate t, ray intersects triangle ---*/

  t = inv_det * GeometryToolbox::DotProduct(3, edge2, qvec);

  /*--- Compute the intersection point in cartesian coordinates ---*/

  intersect[0] = orig[0] + (t * dir[0]);
  intersect[1] = orig[1] + (t * dir[1]);
  intersect[2] = orig[2] + (t * dir[2]);

  return (true);
}

bool CGeometry::SegmentIntersectsLine(const su2double point0[2], const su2double point1[2], const su2double vert0[2],
                                      const su2double vert1[2]) {
  su2double det, diff0_A, diff0_B, diff1_A, diff1_B, intersect[2];

  diff0_A = point0[0] - point1[0];
  diff1_A = point0[1] - point1[1];

  diff0_B = vert0[0] - vert1[0];
  diff1_B = vert0[1] - vert1[1];

  det = (diff0_A) * (diff1_B) - (diff1_A) * (diff0_B);

  if (det == 0) return false;

  /*--- Compute point of intersection ---*/

  intersect[0] = ((point0[0] * point1[1] - point0[1] * point1[0]) * diff0_B -
                  (vert0[0] * vert1[1] - vert0[1] * vert1[0]) * diff0_A) /
                 det;

  intersect[1] = ((point0[0] * point1[1] - point0[1] * point1[0]) * diff1_B -
                  (vert0[0] * vert1[1] - vert0[1] * vert1[0]) * diff1_A) /
                 det;

  /*--- Check that the point is between the two surface points ---*/

  su2double dist0, dist1, length;

  dist0 =
      (intersect[0] - point0[0]) * (intersect[0] - point0[0]) + (intersect[1] - point0[1]) * (intersect[1] - point0[1]);

  dist1 =
      (intersect[0] - point1[0]) * (intersect[0] - point1[0]) + (intersect[1] - point1[1]) * (intersect[1] - point1[1]);

  length = diff0_A * diff0_A + diff1_A * diff1_A;

  return (dist0 <= length) && (dist1 <= length);
}

bool CGeometry::SegmentIntersectsTriangle(su2double point0[3], const su2double point1[3], su2double vert0[3],
                                          su2double vert1[3], su2double vert2[3]) {
  su2double dir[3], intersect[3], u[3], v[3], edge1[3], edge2[3], Plane_Normal[3], Denominator, Numerator, Aux;

  GeometryToolbox::Distance(3, point1, point0, dir);

  if (RayIntersectsTriangle(point0, dir, vert0, vert1, vert2, intersect)) {
    /*--- Check that the intersection is in the segment ---*/

    GeometryToolbox::Distance(3, point0, intersect, u);
    GeometryToolbox::Distance(3, point1, intersect, v);

    GeometryToolbox::Distance(3, vert1, vert0, edge1);
    GeometryToolbox::Distance(3, vert2, vert0, edge2);
    GeometryToolbox::CrossProduct(edge1, edge2, Plane_Normal);

    Denominator = GeometryToolbox::DotProduct(3, Plane_Normal, u);
    Numerator = GeometryToolbox::DotProduct(3, Plane_Normal, v);

    Aux = Numerator * Denominator;

    /*--- Intersection outside the segment ---*/

    if (Aux > 0.0) return (false);

  } else {
    /*--- No intersection with the ray ---*/

    return (false);
  }

  /*--- Intersection inside the segment ---*/

  return (true);
}

void CGeometry::ComputeAirfoil_Section(su2double* Plane_P0, su2double* Plane_Normal, su2double MinXCoord,
                                       su2double MaxXCoord, su2double MinYCoord, su2double MaxYCoord,
                                       su2double MinZCoord, su2double MaxZCoord, const su2double* FlowVariable,
                                       vector<su2double>& Xcoord_Airfoil, vector<su2double>& Ycoord_Airfoil,
                                       vector<su2double>& Zcoord_Airfoil, vector<su2double>& Variable_Airfoil,
                                       bool original_surface, CConfig* config) {
  const bool wasActive = AD::BeginPassive();

  unsigned short iMarker, iNode, jNode, iDim, Index = 0;
  bool intersect;
  long Next_Edge = 0;
  unsigned long iPoint, jPoint, iElem, Trailing_Point, Airfoil_Point, iVertex, iEdge, PointIndex, jEdge;
  su2double Segment_P0[3] = {0.0, 0.0, 0.0}, Segment_P1[3] = {0.0, 0.0, 0.0}, Variable_P0 = 0.0, Variable_P1 = 0.0,
            Intersection[3] = {0.0, 0.0, 0.0}, Trailing_Coord, *VarCoord = nullptr, Variable_Interp,
            v1[3] = {0.0, 0.0, 0.0}, v3[3] = {0.0, 0.0, 0.0}, CrossProduct = 1.0;
  bool Found_Edge;
  passivedouble Dist_Value;
  vector<su2double> Xcoord_Index0, Ycoord_Index0, Zcoord_Index0, Variable_Index0, Xcoord_Index1, Ycoord_Index1,
      Zcoord_Index1, Variable_Index1;
  vector<unsigned long> IGlobalID_Index0, JGlobalID_Index0, IGlobalID_Index1, JGlobalID_Index1, IGlobalID_Airfoil,
      JGlobalID_Airfoil;
  vector<unsigned short> Conection_Index0, Conection_Index1;
  vector<unsigned long> Duplicate;
  su2double** Coord_Variation = nullptr;
  vector<su2double> XcoordExtra, YcoordExtra, ZcoordExtra, VariableExtra;
  vector<unsigned long> IGlobalIDExtra, JGlobalIDExtra;
  vector<bool> AddExtra;
  unsigned long EdgeDonor;
  bool FoundEdge;

#ifdef HAVE_MPI
  unsigned long nLocalEdge, MaxLocalEdge, *Buffer_Send_nEdge, *Buffer_Receive_nEdge, nBuffer_Coord, nBuffer_Variable,
      nBuffer_GlobalID;
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
    Plane_P0[0] = 0.0;
    Plane_P0[1] = 0.0;
    Plane_P0[2] = 0.0;
    Plane_Normal[0] = 0.0;
    Plane_Normal[1] = 1.0;
    Plane_Normal[2] = 0.0;
  }

  /*--- Grid movement is stored using a vertices information,
   we should go from vertex to points ---*/

  if (!original_surface) {
    Coord_Variation = new su2double*[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) Coord_Variation[iPoint] = new su2double[nDim];

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_GeoEval(iMarker) == YES) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          VarCoord = vertex[iMarker][iVertex]->GetVarCoord();
          iPoint = vertex[iMarker][iVertex]->GetNode();
          for (iDim = 0; iDim < nDim; iDim++) Coord_Variation[iPoint][iDim] = VarCoord[iDim];
        }
      }
    }
  }

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_GeoEval(iMarker) == YES) {
      for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
        PointIndex = 0;

        /*--- To decide if an element is going to be used or not should be done element based,
         The first step is to compute and average coordinate for the element ---*/

        su2double AveXCoord = 0.0;
        su2double AveYCoord = 0.0;
        su2double AveZCoord = 0.0;

        for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
          iPoint = bound[iMarker][iElem]->GetNode(iNode);
          AveXCoord += nodes->GetCoord(iPoint, 0);
          AveYCoord += nodes->GetCoord(iPoint, 1);
          if (nDim == 3) AveZCoord += nodes->GetCoord(iPoint, 2);
        }

        AveXCoord /= su2double(bound[iMarker][iElem]->GetnNodes());
        AveYCoord /= su2double(bound[iMarker][iElem]->GetnNodes());
        AveZCoord /= su2double(bound[iMarker][iElem]->GetnNodes());

        /*--- To only cut one part of the nacelle based on the cross product
         of the normal to the plane and a vector that connect the point
         with the center line ---*/

        CrossProduct = 1.0;

        if (config->GetGeo_Description() == NACELLE) {
          su2double Tilt_Angle = config->GetNacelleLocation(3) * PI_NUMBER / 180;
          su2double Toe_Angle = config->GetNacelleLocation(4) * PI_NUMBER / 180;

          /*--- Translate to the origin ---*/

          su2double XCoord_Trans = AveXCoord - config->GetNacelleLocation(0);
          su2double YCoord_Trans = AveYCoord - config->GetNacelleLocation(1);
          su2double ZCoord_Trans = AveZCoord - config->GetNacelleLocation(2);

          /*--- Apply tilt angle ---*/

          su2double XCoord_Trans_Tilt = XCoord_Trans * cos(Tilt_Angle) + ZCoord_Trans * sin(Tilt_Angle);
          su2double YCoord_Trans_Tilt = YCoord_Trans;
          su2double ZCoord_Trans_Tilt = ZCoord_Trans * cos(Tilt_Angle) - XCoord_Trans * sin(Tilt_Angle);

          /*--- Apply toe angle ---*/

          su2double YCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt * sin(Toe_Angle) + YCoord_Trans_Tilt * cos(Toe_Angle);
          su2double ZCoord_Trans_Tilt_Toe = ZCoord_Trans_Tilt;

          /*--- Undo plane rotation, we have already rotated the nacelle ---*/

          /*--- Undo tilt angle ---*/

          su2double XPlane_Normal_Tilt = Plane_Normal[0] * cos(-Tilt_Angle) + Plane_Normal[2] * sin(-Tilt_Angle);
          su2double YPlane_Normal_Tilt = Plane_Normal[1];
          su2double ZPlane_Normal_Tilt = Plane_Normal[2] * cos(-Tilt_Angle) - Plane_Normal[0] * sin(-Tilt_Angle);

          /*--- Undo toe angle ---*/

          su2double YPlane_Normal_Tilt_Toe =
              XPlane_Normal_Tilt * sin(-Toe_Angle) + YPlane_Normal_Tilt * cos(-Toe_Angle);
          su2double ZPlane_Normal_Tilt_Toe = ZPlane_Normal_Tilt;

          v1[1] = YCoord_Trans_Tilt_Toe - 0.0;
          v1[2] = ZCoord_Trans_Tilt_Toe - 0.0;
          v3[0] = v1[1] * ZPlane_Normal_Tilt_Toe - v1[2] * YPlane_Normal_Tilt_Toe;
          CrossProduct = v3[0] * 1.0;
        }

        for (unsigned short iFace = 0; iFace < bound[iMarker][iElem]->GetnFaces(); iFace++) {
          iNode = bound[iMarker][iElem]->GetFaces(iFace, 0);
          jNode = bound[iMarker][iElem]->GetFaces(iFace, 1);
          iPoint = bound[iMarker][iElem]->GetNode(iNode);
          jPoint = bound[iMarker][iElem]->GetNode(jNode);

          if ((CrossProduct >= 0.0) && ((AveXCoord > MinXCoord) && (AveXCoord < MaxXCoord)) &&
              ((AveYCoord > MinYCoord) && (AveYCoord < MaxYCoord)) &&
              ((AveZCoord > MinZCoord) && (AveZCoord < MaxZCoord))) {
            Segment_P0[0] = 0.0;
            Segment_P0[1] = 0.0;
            Segment_P0[2] = 0.0;
            Variable_P0 = 0.0;
            Segment_P1[0] = 0.0;
            Segment_P1[1] = 0.0;
            Segment_P1[2] = 0.0;
            Variable_P1 = 0.0;

            for (iDim = 0; iDim < nDim; iDim++) {
              if (original_surface) {
                Segment_P0[iDim] = nodes->GetCoord(iPoint, iDim);
                Segment_P1[iDim] = nodes->GetCoord(jPoint, iDim);
              } else {
                Segment_P0[iDim] = nodes->GetCoord(iPoint, iDim) + Coord_Variation[iPoint][iDim];
                Segment_P1[iDim] = nodes->GetCoord(jPoint, iDim) + Coord_Variation[jPoint][iDim];
              }
            }

            if (FlowVariable != nullptr) {
              Variable_P0 = FlowVariable[iPoint];
              Variable_P1 = FlowVariable[jPoint];
            }

            /*--- In 2D add the points directly (note the change between Y and Z coordinate) ---*/

            if (nDim == 2) {
              Xcoord_Index0.push_back(Segment_P0[0]);
              Xcoord_Index1.push_back(Segment_P1[0]);
              Ycoord_Index0.push_back(Segment_P0[2]);
              Ycoord_Index1.push_back(Segment_P1[2]);
              Zcoord_Index0.push_back(Segment_P0[1]);
              Zcoord_Index1.push_back(Segment_P1[1]);
              Variable_Index0.push_back(Variable_P0);
              Variable_Index1.push_back(Variable_P1);
              IGlobalID_Index0.push_back(nodes->GetGlobalIndex(iPoint));
              IGlobalID_Index1.push_back(nodes->GetGlobalIndex(jPoint));
              JGlobalID_Index0.push_back(nodes->GetGlobalIndex(iPoint));
              JGlobalID_Index1.push_back(nodes->GetGlobalIndex(jPoint));
              PointIndex++;
            }

            /*--- In 3D compute the intersection ---*/

            else if (nDim == 3) {
              intersect = SegmentIntersectsPlane(Segment_P0, Segment_P1, Variable_P0, Variable_P1, Plane_P0,
                                                 Plane_Normal, Intersection, Variable_Interp);
              if (intersect) {
                if (PointIndex == 0) {
                  Xcoord_Index0.push_back(Intersection[0]);
                  Ycoord_Index0.push_back(Intersection[1]);
                  Zcoord_Index0.push_back(Intersection[2]);
                  Variable_Index0.push_back(Variable_Interp);
                  IGlobalID_Index0.push_back(nodes->GetGlobalIndex(iPoint));
                  JGlobalID_Index0.push_back(nodes->GetGlobalIndex(jPoint));
                }
                if (PointIndex == 1) {
                  Xcoord_Index1.push_back(Intersection[0]);
                  Ycoord_Index1.push_back(Intersection[1]);
                  Zcoord_Index1.push_back(Intersection[2]);
                  Variable_Index1.push_back(Variable_Interp);
                  IGlobalID_Index1.push_back(nodes->GetGlobalIndex(iPoint));
                  JGlobalID_Index1.push_back(nodes->GetGlobalIndex(jPoint));
                }
                PointIndex++;
              }
            }
          }
        }
      }
    }
  }

  if (!original_surface) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) delete[] Coord_Variation[iPoint];
    delete[] Coord_Variation;
  }

#ifdef HAVE_MPI

  /*--- Copy the coordinates of all the points in the plane to the master node ---*/

  nLocalEdge = 0, MaxLocalEdge = 0;
  nProcessor = size;

  Buffer_Send_nEdge = new unsigned long[1];
  Buffer_Receive_nEdge = new unsigned long[nProcessor];

  nLocalEdge = Xcoord_Index0.size();

  Buffer_Send_nEdge[0] = nLocalEdge;

  SU2_MPI::Allreduce(&nLocalEdge, &MaxLocalEdge, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Allgather(Buffer_Send_nEdge, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nEdge, 1, MPI_UNSIGNED_LONG,
                     SU2_MPI::GetComm());

  Buffer_Send_Coord = new su2double[MaxLocalEdge * 6];
  Buffer_Receive_Coord = new su2double[nProcessor * MaxLocalEdge * 6];

  Buffer_Send_Variable = new su2double[MaxLocalEdge * 2];
  Buffer_Receive_Variable = new su2double[nProcessor * MaxLocalEdge * 2];

  Buffer_Send_GlobalID = new unsigned long[MaxLocalEdge * 4];
  Buffer_Receive_GlobalID = new unsigned long[nProcessor * MaxLocalEdge * 4];

  nBuffer_Coord = MaxLocalEdge * 6;
  nBuffer_Variable = MaxLocalEdge * 2;
  nBuffer_GlobalID = MaxLocalEdge * 4;

  for (iEdge = 0; iEdge < nLocalEdge; iEdge++) {
    Buffer_Send_Coord[iEdge * 6 + 0] = Xcoord_Index0[iEdge];
    Buffer_Send_Coord[iEdge * 6 + 1] = Ycoord_Index0[iEdge];
    Buffer_Send_Coord[iEdge * 6 + 2] = Zcoord_Index0[iEdge];
    Buffer_Send_Coord[iEdge * 6 + 3] = Xcoord_Index1[iEdge];
    Buffer_Send_Coord[iEdge * 6 + 4] = Ycoord_Index1[iEdge];
    Buffer_Send_Coord[iEdge * 6 + 5] = Zcoord_Index1[iEdge];

    Buffer_Send_Variable[iEdge * 2 + 0] = Variable_Index0[iEdge];
    Buffer_Send_Variable[iEdge * 2 + 1] = Variable_Index1[iEdge];

    Buffer_Send_GlobalID[iEdge * 4 + 0] = IGlobalID_Index0[iEdge];
    Buffer_Send_GlobalID[iEdge * 4 + 1] = JGlobalID_Index0[iEdge];
    Buffer_Send_GlobalID[iEdge * 4 + 2] = IGlobalID_Index1[iEdge];
    Buffer_Send_GlobalID[iEdge * 4 + 3] = JGlobalID_Index1[iEdge];
  }

  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE,
                     SU2_MPI::GetComm());
  SU2_MPI::Allgather(Buffer_Send_Variable, nBuffer_Variable, MPI_DOUBLE, Buffer_Receive_Variable, nBuffer_Variable,
                     MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(Buffer_Send_GlobalID, nBuffer_GlobalID, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalID,
                     nBuffer_GlobalID, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  /*--- Clean the vectors before adding the new vertices only to the master node ---*/

  Xcoord_Index0.clear();
  Xcoord_Index1.clear();
  Ycoord_Index0.clear();
  Ycoord_Index1.clear();
  Zcoord_Index0.clear();
  Zcoord_Index1.clear();
  Variable_Index0.clear();
  Variable_Index1.clear();
  IGlobalID_Index0.clear();
  IGlobalID_Index1.clear();
  JGlobalID_Index0.clear();
  JGlobalID_Index1.clear();

  /*--- Copy the boundary to the master node vectors ---*/

  if (rank == MASTER_NODE) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iEdge = 0; iEdge < Buffer_Receive_nEdge[iProcessor]; iEdge++) {
        Xcoord_Index0.push_back(Buffer_Receive_Coord[iProcessor * MaxLocalEdge * 6 + iEdge * 6 + 0]);
        Ycoord_Index0.push_back(Buffer_Receive_Coord[iProcessor * MaxLocalEdge * 6 + iEdge * 6 + 1]);
        Zcoord_Index0.push_back(Buffer_Receive_Coord[iProcessor * MaxLocalEdge * 6 + iEdge * 6 + 2]);
        Xcoord_Index1.push_back(Buffer_Receive_Coord[iProcessor * MaxLocalEdge * 6 + iEdge * 6 + 3]);
        Ycoord_Index1.push_back(Buffer_Receive_Coord[iProcessor * MaxLocalEdge * 6 + iEdge * 6 + 4]);
        Zcoord_Index1.push_back(Buffer_Receive_Coord[iProcessor * MaxLocalEdge * 6 + iEdge * 6 + 5]);

        Variable_Index0.push_back(Buffer_Receive_Variable[iProcessor * MaxLocalEdge * 2 + iEdge * 2 + 0]);
        Variable_Index1.push_back(Buffer_Receive_Variable[iProcessor * MaxLocalEdge * 2 + iEdge * 2 + 1]);

        IGlobalID_Index0.push_back(Buffer_Receive_GlobalID[iProcessor * MaxLocalEdge * 4 + iEdge * 4 + 0]);
        JGlobalID_Index0.push_back(Buffer_Receive_GlobalID[iProcessor * MaxLocalEdge * 4 + iEdge * 4 + 1]);
        IGlobalID_Index1.push_back(Buffer_Receive_GlobalID[iProcessor * MaxLocalEdge * 4 + iEdge * 4 + 2]);
        JGlobalID_Index1.push_back(Buffer_Receive_GlobalID[iProcessor * MaxLocalEdge * 4 + iEdge * 4 + 3]);
      }
    }
  }

  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_Variable;
  delete[] Buffer_Receive_Variable;
  delete[] Buffer_Send_GlobalID;
  delete[] Buffer_Receive_GlobalID;
  delete[] Buffer_Send_nEdge;
  delete[] Buffer_Receive_nEdge;

#endif

  if ((rank == MASTER_NODE) && (!Xcoord_Index0.empty())) {
    /*--- Remove singular edges ---*/

    bool Remove;

    do {
      Remove = false;
      for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {
        if (((IGlobalID_Index0[iEdge] == IGlobalID_Index1[iEdge]) &&
             (JGlobalID_Index0[iEdge] == JGlobalID_Index1[iEdge])) ||
            ((IGlobalID_Index0[iEdge] == JGlobalID_Index1[iEdge]) &&
             (JGlobalID_Index0[iEdge] == IGlobalID_Index1[iEdge]))) {
          Xcoord_Index0.erase(Xcoord_Index0.begin() + iEdge);
          Ycoord_Index0.erase(Ycoord_Index0.begin() + iEdge);
          Zcoord_Index0.erase(Zcoord_Index0.begin() + iEdge);
          Variable_Index0.erase(Variable_Index0.begin() + iEdge);
          IGlobalID_Index0.erase(IGlobalID_Index0.begin() + iEdge);
          JGlobalID_Index0.erase(JGlobalID_Index0.begin() + iEdge);

          Xcoord_Index1.erase(Xcoord_Index1.begin() + iEdge);
          Ycoord_Index1.erase(Ycoord_Index1.begin() + iEdge);
          Zcoord_Index1.erase(Zcoord_Index1.begin() + iEdge);
          Variable_Index1.erase(Variable_Index1.begin() + iEdge);
          IGlobalID_Index1.erase(IGlobalID_Index1.begin() + iEdge);
          JGlobalID_Index1.erase(JGlobalID_Index1.begin() + iEdge);

          Remove = true;
          break;
        }
        if (Remove) break;
      }
    } while (Remove);

    /*--- Remove repeated edges computing distance, this could happend because the MPI ---*/

    do {
      Remove = false;
      for (iEdge = 0; iEdge < Xcoord_Index0.size() - 1; iEdge++) {
        for (jEdge = iEdge + 1; jEdge < Xcoord_Index0.size(); jEdge++) {
          /*--- Edges with the same orientation ---*/

          if ((((IGlobalID_Index0[iEdge] == IGlobalID_Index0[jEdge]) &&
                (JGlobalID_Index0[iEdge] == JGlobalID_Index0[jEdge])) ||
               ((IGlobalID_Index0[iEdge] == JGlobalID_Index0[jEdge]) &&
                (JGlobalID_Index0[iEdge] == IGlobalID_Index0[jEdge]))) &&
              (((IGlobalID_Index1[iEdge] == IGlobalID_Index1[jEdge]) &&
                (JGlobalID_Index1[iEdge] == JGlobalID_Index1[jEdge])) ||
               ((IGlobalID_Index1[iEdge] == JGlobalID_Index1[jEdge]) &&
                (JGlobalID_Index1[iEdge] == IGlobalID_Index1[jEdge])))) {
            Xcoord_Index0.erase(Xcoord_Index0.begin() + jEdge);
            Ycoord_Index0.erase(Ycoord_Index0.begin() + jEdge);
            Zcoord_Index0.erase(Zcoord_Index0.begin() + jEdge);
            Variable_Index0.erase(Variable_Index0.begin() + jEdge);
            IGlobalID_Index0.erase(IGlobalID_Index0.begin() + jEdge);
            JGlobalID_Index0.erase(JGlobalID_Index0.begin() + jEdge);

            Xcoord_Index1.erase(Xcoord_Index1.begin() + jEdge);
            Ycoord_Index1.erase(Ycoord_Index1.begin() + jEdge);
            Zcoord_Index1.erase(Zcoord_Index1.begin() + jEdge);
            Variable_Index1.erase(Variable_Index1.begin() + jEdge);
            IGlobalID_Index1.erase(IGlobalID_Index1.begin() + jEdge);
            JGlobalID_Index1.erase(JGlobalID_Index1.begin() + jEdge);

            Remove = true;
            break;
          }

          /*--- Edges with oposite orientation ---*/

          if ((((IGlobalID_Index0[iEdge] == IGlobalID_Index1[jEdge]) &&
                (JGlobalID_Index0[iEdge] == JGlobalID_Index1[jEdge])) ||
               ((IGlobalID_Index0[iEdge] == JGlobalID_Index1[jEdge]) &&
                (JGlobalID_Index0[iEdge] == IGlobalID_Index1[jEdge]))) &&
              (((IGlobalID_Index1[iEdge] == IGlobalID_Index0[jEdge]) &&
                (JGlobalID_Index1[iEdge] == JGlobalID_Index0[jEdge])) ||
               ((IGlobalID_Index1[iEdge] == JGlobalID_Index0[jEdge]) &&
                (JGlobalID_Index1[iEdge] == IGlobalID_Index0[jEdge])))) {
            Xcoord_Index0.erase(Xcoord_Index0.begin() + jEdge);
            Ycoord_Index0.erase(Ycoord_Index0.begin() + jEdge);
            Zcoord_Index0.erase(Zcoord_Index0.begin() + jEdge);
            Variable_Index0.erase(Variable_Index0.begin() + jEdge);
            IGlobalID_Index0.erase(IGlobalID_Index0.begin() + jEdge);
            JGlobalID_Index0.erase(JGlobalID_Index0.begin() + jEdge);

            Xcoord_Index1.erase(Xcoord_Index1.begin() + jEdge);
            Ycoord_Index1.erase(Ycoord_Index1.begin() + jEdge);
            Zcoord_Index1.erase(Zcoord_Index1.begin() + jEdge);
            Variable_Index1.erase(Variable_Index1.begin() + jEdge);
            IGlobalID_Index1.erase(IGlobalID_Index1.begin() + jEdge);
            JGlobalID_Index1.erase(JGlobalID_Index1.begin() + jEdge);

            Remove = true;
            break;
          }
          if (Remove) break;
        }
        if (Remove) break;
      }

    } while (Remove);

    if (Xcoord_Index0.size() != 1) {
      /*--- Rotate from the Y-Z plane to the X-Z plane to reuse the rest of subroutines  ---*/

      if (config->GetGeo_Description() == FUSELAGE) {
        su2double Angle = -0.5 * PI_NUMBER;
        for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {
          su2double XCoord = Xcoord_Index0[iEdge] * cos(Angle) - Ycoord_Index0[iEdge] * sin(Angle);
          su2double YCoord = Ycoord_Index0[iEdge] * cos(Angle) + Xcoord_Index0[iEdge] * sin(Angle);
          su2double ZCoord = Zcoord_Index0[iEdge];
          Xcoord_Index0[iEdge] = XCoord;
          Ycoord_Index0[iEdge] = YCoord;
          Zcoord_Index0[iEdge] = ZCoord;
          XCoord = Xcoord_Index1[iEdge] * cos(Angle) - Ycoord_Index1[iEdge] * sin(Angle);
          YCoord = Ycoord_Index1[iEdge] * cos(Angle) + Xcoord_Index1[iEdge] * sin(Angle);
          ZCoord = Zcoord_Index1[iEdge];
          Xcoord_Index1[iEdge] = XCoord;
          Ycoord_Index1[iEdge] = YCoord;
          Zcoord_Index1[iEdge] = ZCoord;
        }
      }

      /*--- Rotate nacelle secction to a X-Z plane to reuse the rest of subroutines  ---*/

      if (config->GetGeo_Description() == NACELLE) {
        su2double Tilt_Angle = config->GetNacelleLocation(3) * PI_NUMBER / 180;
        su2double Toe_Angle = config->GetNacelleLocation(4) * PI_NUMBER / 180;
        su2double Theta_deg = atan2(Plane_Normal[1], -Plane_Normal[2]) / PI_NUMBER * 180 + 180;
        su2double Roll_Angle = 0.5 * PI_NUMBER - Theta_deg * PI_NUMBER / 180;

        su2double XCoord_Trans, YCoord_Trans, ZCoord_Trans, XCoord_Trans_Tilt, YCoord_Trans_Tilt, ZCoord_Trans_Tilt,
            XCoord_Trans_Tilt_Toe, YCoord_Trans_Tilt_Toe, ZCoord_Trans_Tilt_Toe, XCoord, YCoord, ZCoord;

        for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {
          /*--- First point of the edge ---*/

          /*--- Translate to the origin ---*/

          XCoord_Trans = Xcoord_Index0[iEdge] - config->GetNacelleLocation(0);
          YCoord_Trans = Ycoord_Index0[iEdge] - config->GetNacelleLocation(1);
          ZCoord_Trans = Zcoord_Index0[iEdge] - config->GetNacelleLocation(2);

          /*--- Apply tilt angle ---*/

          XCoord_Trans_Tilt = XCoord_Trans * cos(Tilt_Angle) + ZCoord_Trans * sin(Tilt_Angle);
          YCoord_Trans_Tilt = YCoord_Trans;
          ZCoord_Trans_Tilt = ZCoord_Trans * cos(Tilt_Angle) - XCoord_Trans * sin(Tilt_Angle);

          /*--- Apply toe angle ---*/

          XCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt * cos(Toe_Angle) - YCoord_Trans_Tilt * sin(Toe_Angle);
          YCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt * sin(Toe_Angle) + YCoord_Trans_Tilt * cos(Toe_Angle);
          ZCoord_Trans_Tilt_Toe = ZCoord_Trans_Tilt;

          /*--- Rotate to X-Z plane (roll) ---*/

          XCoord = XCoord_Trans_Tilt_Toe;
          YCoord = YCoord_Trans_Tilt_Toe * cos(Roll_Angle) - ZCoord_Trans_Tilt_Toe * sin(Roll_Angle);
          ZCoord = YCoord_Trans_Tilt_Toe * sin(Roll_Angle) + ZCoord_Trans_Tilt_Toe * cos(Roll_Angle);

          /*--- Update coordinates ---*/

          Xcoord_Index0[iEdge] = XCoord;
          Ycoord_Index0[iEdge] = YCoord;
          Zcoord_Index0[iEdge] = ZCoord;

          /*--- Second point of the edge ---*/

          /*--- Translate to the origin ---*/

          XCoord_Trans = Xcoord_Index1[iEdge] - config->GetNacelleLocation(0);
          YCoord_Trans = Ycoord_Index1[iEdge] - config->GetNacelleLocation(1);
          ZCoord_Trans = Zcoord_Index1[iEdge] - config->GetNacelleLocation(2);

          /*--- Apply tilt angle ---*/

          XCoord_Trans_Tilt = XCoord_Trans * cos(Tilt_Angle) + ZCoord_Trans * sin(Tilt_Angle);
          YCoord_Trans_Tilt = YCoord_Trans;
          ZCoord_Trans_Tilt = ZCoord_Trans * cos(Tilt_Angle) - XCoord_Trans * sin(Tilt_Angle);

          /*--- Apply toe angle ---*/

          XCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt * cos(Toe_Angle) - YCoord_Trans_Tilt * sin(Toe_Angle);
          YCoord_Trans_Tilt_Toe = XCoord_Trans_Tilt * sin(Toe_Angle) + YCoord_Trans_Tilt * cos(Toe_Angle);
          ZCoord_Trans_Tilt_Toe = ZCoord_Trans_Tilt;

          /*--- Rotate to X-Z plane (roll) ---*/

          XCoord = XCoord_Trans_Tilt_Toe;
          YCoord = YCoord_Trans_Tilt_Toe * cos(Roll_Angle) - ZCoord_Trans_Tilt_Toe * sin(Roll_Angle);
          ZCoord = YCoord_Trans_Tilt_Toe * sin(Roll_Angle) + ZCoord_Trans_Tilt_Toe * cos(Roll_Angle);

          /*--- Update coordinates ---*/

          Xcoord_Index1[iEdge] = XCoord;
          Ycoord_Index1[iEdge] = YCoord;
          Zcoord_Index1[iEdge] = ZCoord;
        }
      }

      /*--- Identify the extreme of the curve and close it ---*/

      Conection_Index0.reserve(Xcoord_Index0.size() + 1);
      Conection_Index1.reserve(Xcoord_Index0.size() + 1);

      for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {
        Conection_Index0[iEdge] = 0;
        Conection_Index1[iEdge] = 0;
      }

      for (iEdge = 0; iEdge < Xcoord_Index0.size() - 1; iEdge++) {
        for (jEdge = iEdge + 1; jEdge < Xcoord_Index0.size(); jEdge++) {
          if (((IGlobalID_Index0[iEdge] == IGlobalID_Index0[jEdge]) &&
               (JGlobalID_Index0[iEdge] == JGlobalID_Index0[jEdge])) ||
              ((IGlobalID_Index0[iEdge] == JGlobalID_Index0[jEdge]) &&
               (JGlobalID_Index0[iEdge] == IGlobalID_Index0[jEdge]))) {
            Conection_Index0[iEdge]++;
            Conection_Index0[jEdge]++;
          }

          if (((IGlobalID_Index0[iEdge] == IGlobalID_Index1[jEdge]) &&
               (JGlobalID_Index0[iEdge] == JGlobalID_Index1[jEdge])) ||
              ((IGlobalID_Index0[iEdge] == JGlobalID_Index1[jEdge]) &&
               (JGlobalID_Index0[iEdge] == IGlobalID_Index1[jEdge]))) {
            Conection_Index0[iEdge]++;
            Conection_Index1[jEdge]++;
          }

          if (((IGlobalID_Index1[iEdge] == IGlobalID_Index0[jEdge]) &&
               (JGlobalID_Index1[iEdge] == JGlobalID_Index0[jEdge])) ||
              ((IGlobalID_Index1[iEdge] == JGlobalID_Index0[jEdge]) &&
               (JGlobalID_Index1[iEdge] == IGlobalID_Index0[jEdge]))) {
            Conection_Index1[iEdge]++;
            Conection_Index0[jEdge]++;
          }

          if (((IGlobalID_Index1[iEdge] == IGlobalID_Index1[jEdge]) &&
               (JGlobalID_Index1[iEdge] == JGlobalID_Index1[jEdge])) ||
              ((IGlobalID_Index1[iEdge] == JGlobalID_Index1[jEdge]) &&
               (JGlobalID_Index1[iEdge] == IGlobalID_Index1[jEdge]))) {
            Conection_Index1[iEdge]++;
            Conection_Index1[jEdge]++;
          }
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
        for (iEdge = 0; iEdge < XcoordExtra.size() - 1; iEdge++) {
          su2double MinDist = 1E6;
          FoundEdge = false;
          EdgeDonor = 0;
          for (jEdge = iEdge + 1; jEdge < XcoordExtra.size(); jEdge++) {
            Dist_Value =
                sqrt(pow(SU2_TYPE::GetValue(XcoordExtra[iEdge]) - SU2_TYPE::GetValue(XcoordExtra[jEdge]), 2.0));
            if ((Dist_Value < MinDist) && (AddExtra[iEdge]) && (AddExtra[jEdge])) {
              EdgeDonor = jEdge;
              FoundEdge = true;
            }
          }

          if (FoundEdge) {
            /*--- Add first point of the new edge ---*/

            Xcoord_Index0.push_back(XcoordExtra[iEdge]);
            Ycoord_Index0.push_back(YcoordExtra[iEdge]);
            Zcoord_Index0.push_back(ZcoordExtra[iEdge]);
            Variable_Index0.push_back(VariableExtra[iEdge]);
            IGlobalID_Index0.push_back(IGlobalIDExtra[iEdge]);
            JGlobalID_Index0.push_back(JGlobalIDExtra[iEdge]);
            AddExtra[iEdge] = false;

            /*--- Add second (closest)  point of the new edge ---*/

            Xcoord_Index1.push_back(XcoordExtra[EdgeDonor]);
            Ycoord_Index1.push_back(YcoordExtra[EdgeDonor]);
            Zcoord_Index1.push_back(ZcoordExtra[EdgeDonor]);
            Variable_Index1.push_back(VariableExtra[EdgeDonor]);
            IGlobalID_Index1.push_back(IGlobalIDExtra[EdgeDonor]);
            JGlobalID_Index1.push_back(JGlobalIDExtra[EdgeDonor]);
            AddExtra[EdgeDonor] = false;
          }
        }

      }

      else if (XcoordExtra.size() == 1) {
        cout << "There cutting system has failed, there is an incomplete curve (not used)." << endl;
      }

      /*--- Find and add the trailing edge to to the list
       and the contect the first point to the trailing edge ---*/

      Trailing_Point = 0;
      Trailing_Coord = Xcoord_Index0[0];
      for (iEdge = 1; iEdge < Xcoord_Index0.size(); iEdge++) {
        if (Xcoord_Index0[iEdge] > Trailing_Coord) {
          Trailing_Point = iEdge;
          Trailing_Coord = Xcoord_Index0[iEdge];
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

      Xcoord_Index0.erase(Xcoord_Index0.begin() + Trailing_Point);
      Ycoord_Index0.erase(Ycoord_Index0.begin() + Trailing_Point);
      Zcoord_Index0.erase(Zcoord_Index0.begin() + Trailing_Point);
      Variable_Index0.erase(Variable_Index0.begin() + Trailing_Point);
      IGlobalID_Index0.erase(IGlobalID_Index0.begin() + Trailing_Point);
      JGlobalID_Index0.erase(JGlobalID_Index0.begin() + Trailing_Point);

      Xcoord_Index1.erase(Xcoord_Index1.begin() + Trailing_Point);
      Ycoord_Index1.erase(Ycoord_Index1.begin() + Trailing_Point);
      Zcoord_Index1.erase(Zcoord_Index1.begin() + Trailing_Point);
      Variable_Index1.erase(Variable_Index1.begin() + Trailing_Point);
      IGlobalID_Index1.erase(IGlobalID_Index1.begin() + Trailing_Point);
      JGlobalID_Index1.erase(JGlobalID_Index1.begin() + Trailing_Point);

      /*--- Algorithm for adding the rest of the points ---*/

      do {
        /*--- Last added point in the list ---*/

        Airfoil_Point = Xcoord_Airfoil.size() - 1;

        /*--- Find the closest point  ---*/

        Found_Edge = false;

        for (iEdge = 0; iEdge < Xcoord_Index0.size(); iEdge++) {
          if (((IGlobalID_Index0[iEdge] == IGlobalID_Airfoil[Airfoil_Point]) &&
               (JGlobalID_Index0[iEdge] == JGlobalID_Airfoil[Airfoil_Point])) ||
              ((IGlobalID_Index0[iEdge] == JGlobalID_Airfoil[Airfoil_Point]) &&
               (JGlobalID_Index0[iEdge] == IGlobalID_Airfoil[Airfoil_Point]))) {
            Next_Edge = iEdge;
            Found_Edge = true;
            Index = 0;
            break;
          }

          if (((IGlobalID_Index1[iEdge] == IGlobalID_Airfoil[Airfoil_Point]) &&
               (JGlobalID_Index1[iEdge] == JGlobalID_Airfoil[Airfoil_Point])) ||
              ((IGlobalID_Index1[iEdge] == JGlobalID_Airfoil[Airfoil_Point]) &&
               (JGlobalID_Index1[iEdge] == IGlobalID_Airfoil[Airfoil_Point]))) {
            Next_Edge = iEdge;
            Found_Edge = true;
            Index = 1;
            break;
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

        } else {
          break;
        }

      } while (!Xcoord_Index0.empty());

      /*--- Clean the vector before using them again for storing the upper or the lower side ---*/

      Xcoord_Index0.clear();
      Ycoord_Index0.clear();
      Zcoord_Index0.clear();
      Variable_Index0.clear();
      IGlobalID_Index0.clear();
      JGlobalID_Index0.clear();
      Xcoord_Index1.clear();
      Ycoord_Index1.clear();
      Zcoord_Index1.clear();
      Variable_Index1.clear();
      IGlobalID_Index1.clear();
      JGlobalID_Index1.clear();
    }
  }

  AD::EndPassive(wasActive);
}

void CGeometry::RegisterCoordinates() const {
  const bool input = true;

  SU2_OMP_FOR_STAT(roundUpDiv(nPoint, omp_get_num_threads()))
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    nodes->RegisterCoordinates(iPoint, input);
  }
  END_SU2_OMP_FOR
}

void CGeometry::UpdateGeometry(CGeometry** geometry_container, CConfig* config) {
  geometry_container[MESH_0]->InitiateComms(geometry_container[MESH_0], config, COORDINATES);
  geometry_container[MESH_0]->CompleteComms(geometry_container[MESH_0], config, COORDINATES);

  geometry_container[MESH_0]->SetControlVolume(config, UPDATE);
  geometry_container[MESH_0]->SetBoundControlVolume(config, UPDATE);
  geometry_container[MESH_0]->SetMaxLength(config);

  for (unsigned short iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    /*--- Update the control volume structures ---*/

    geometry_container[iMesh]->SetControlVolume(geometry_container[iMesh - 1], UPDATE);
    geometry_container[iMesh]->SetBoundControlVolume(geometry_container[iMesh - 1], UPDATE);
    geometry_container[iMesh]->SetCoord(geometry_container[iMesh - 1]);
  }

  /*--- Compute the global surface areas for all markers. ---*/
  geometry_container[MESH_0]->ComputeSurfaceAreaCfgFile(config);
}

void CGeometry::SetCustomBoundary(CConfig* config) {
  unsigned short iMarker;
  unsigned long iVertex;
  string Marker_Tag;

  /*--- Initialize quantities for customized boundary conditions.
   * Custom values are initialized with the default values specified in the config (avoiding non physical values) ---*/
  CustomBoundaryTemperature = new su2double*[nMarker];
  CustomBoundaryHeatFlux = new su2double*[nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    CustomBoundaryHeatFlux[iMarker] = nullptr;
    CustomBoundaryTemperature[iMarker] = nullptr;
    if (config->GetMarker_All_PyCustom(iMarker)) {
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case HEAT_FLUX:
          CustomBoundaryHeatFlux[iMarker] = new su2double[nVertex[iMarker]];
          for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
            CustomBoundaryHeatFlux[iMarker][iVertex] = config->GetWall_HeatFlux(Marker_Tag);
          }
          break;
        case ISOTHERMAL:
          CustomBoundaryTemperature[iMarker] = new su2double[nVertex[iMarker]];
          for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
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

void CGeometry::UpdateCustomBoundaryConditions(CGeometry** geometry_container, CConfig* config) {
  unsigned short iMGfine, iMGlevel, nMGlevel, iMarker;

  nMGlevel = config->GetnMGLevels();
  for (iMGlevel = 1; iMGlevel <= nMGlevel; iMGlevel++) {
    iMGfine = iMGlevel - 1;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_PyCustom(iMarker)) {
        switch (config->GetMarker_All_KindBC(iMarker)) {
          case HEAT_FLUX:
            geometry_container[iMGlevel]->SetMultiGridWallHeatFlux(geometry_container[iMGfine], iMarker);
            break;
          case ISOTHERMAL:
            geometry_container[iMGlevel]->SetMultiGridWallTemperature(geometry_container[iMGfine], iMarker);
            break;
          // Inlet flow handled in solver class.
          default:
            break;
        }
      }
    }
  }
}

void CGeometry::ComputeSurfaceAreaCfgFile(const CConfig* config){
    SU2_OMP_MASTER{const auto nMarker_Global = config->GetnMarker_CfgFile();
SurfaceAreaCfgFile.resize(nMarker_Global);
vector<su2double> LocalSurfaceArea(nMarker_Global, 0.0);

/*--- Loop over all local markers ---*/
for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++) {
  const auto Local_TagBound = config->GetMarker_All_TagBound(iMarker);

  /*--- Loop over all global markers, and find the local-global pair via
        matching unique string tags. ---*/
  for (unsigned short iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++) {
    const auto Global_TagBound = config->GetMarker_CfgFile_TagBound(iMarker_Global);
    if (Local_TagBound == Global_TagBound) {
      for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) {
        const auto iPoint = vertex[iMarker][iVertex]->GetNode();

        if (!nodes->GetDomain(iPoint)) continue;

        const auto AreaNormal = vertex[iMarker][iVertex]->GetNormal();
        const auto Area = GeometryToolbox::Norm(nDim, AreaNormal);

        LocalSurfaceArea[iMarker_Global] += Area;
      }  // for iVertex
    }    // if Local == Global
  }      // for iMarker_Global
}  // for iMarker

SU2_MPI::Allreduce(LocalSurfaceArea.data(), SurfaceAreaCfgFile.data(), SurfaceAreaCfgFile.size(), MPI_DOUBLE, MPI_SUM,
                   SU2_MPI::GetComm());
}
END_SU2_OMP_MASTER
}

su2double CGeometry::GetSurfaceArea(const CConfig* config, unsigned short val_marker) const {
  /*---Find the precomputed marker surface area by local-global string-matching. ---*/
  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  for (unsigned short iMarker_Global = 0; iMarker_Global < config->GetnMarker_CfgFile(); iMarker_Global++) {
    const auto Global_TagBound = config->GetMarker_CfgFile_TagBound(iMarker_Global);

    if (Marker_Tag == Global_TagBound) return SurfaceAreaCfgFile[iMarker_Global];
  }

  SU2_MPI::Error("Unable to match local-marker with cfg-marker for Surface Area.", CURRENT_FUNCTION);
  return 0.0;
}

void CGeometry::ComputeSurf_Straightness(CConfig* config, bool print_on_screen) {
  bool RefUnitNormal_defined;
  unsigned short iDim, iMarker, iMarker_Global, nMarker_Global = config->GetnMarker_CfgFile();
  unsigned long iVertex;
  constexpr passivedouble epsilon = 1.0e-6;
  su2double Area;
  string Local_TagBound, Global_TagBound;

  vector<su2double> Normal(nDim), UnitNormal(nDim), RefUnitNormal(nDim);

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
        !config->GetMarker_Moving_Bool(Local_TagBound) && !config->GetMarker_Deform_Mesh_Bool(Local_TagBound)) {
      /*--- Loop over all global markers, and find the local-global pair via
            matching unique string tags. ---*/
      for (iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++) {
        Global_TagBound = config->GetMarker_CfgFile_TagBound(iMarker_Global);
        if (Local_TagBound == Global_TagBound) {
          RefUnitNormal_defined = false;
          iVertex = 0;

          while (bound_is_straight[iMarker] && iVertex < nVertex[iMarker]) {
            vertex[iMarker][iVertex]->GetNormal(Normal.data());
            UnitNormal = Normal;

            /*--- Compute unit normal. ---*/
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim] * Normal[iDim];
            Area = sqrt(Area);

            /*--- Negate for outward convention. ---*/
            for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] /= -Area;

            /*--- Check if unit normal is within tolerance of the Reference unit normal.
                  Reference unit normal = first unit normal found. ---*/
            if (RefUnitNormal_defined) {
              for (iDim = 0; iDim < nDim; iDim++) {
                if (abs(RefUnitNormal[iDim] - UnitNormal[iDim]) > epsilon) {
                  bound_is_straight[iMarker] = false;
                  break;
                }
              }
            } else {
              RefUnitNormal = UnitNormal;  // deep copy of values
              RefUnitNormal_defined = true;
            }

            iVertex++;
          }  // while iVertex
        }    // if Local == Global
      }      // for iMarker_Global
    } else {
      /*--- Enforce default value: false ---*/
      bound_is_straight[iMarker] = false;
    }  // if sym or euler ...
  }    // for iMarker

  /*--- Communicate results and print on screen. ---*/
  if (print_on_screen) {
    /*--- Additional vector which can later be MPI::Allreduce(d) to pring the results
          on screen as nMarker (local) can vary across ranks. Default 'true' as it can
          happen that a local rank does not contain an element of each surface marker.  ---*/
    vector<bool> bound_is_straight_Global(nMarker_Global, true);
    /*--- Match local with global tag bound and fill a Global Marker vector. ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Local_TagBound = config->GetMarker_All_TagBound(iMarker);
      for (iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++) {
        Global_TagBound = config->GetMarker_CfgFile_TagBound(iMarker_Global);

        if (Local_TagBound == Global_TagBound) bound_is_straight_Global[iMarker_Global] = bound_is_straight[iMarker];

      }  // for iMarker_Global
    }    // for iMarker

    vector<int> Buff_Send_isStraight(nMarker_Global), Buff_Recv_isStraight(nMarker_Global);

    /*--- Cast to int as std::vector<boolean> can be a special construct. MPI handling using <int>
          is more straight-forward. ---*/
    for (iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++)
      Buff_Send_isStraight[iMarker_Global] = static_cast<int>(bound_is_straight_Global[iMarker_Global]);

    /*--- Product of type <int>(bool) is equivalnt to a 'logical and' ---*/
    SU2_MPI::Allreduce(Buff_Send_isStraight.data(), Buff_Recv_isStraight.data(), nMarker_Global, MPI_INT, MPI_PROD,
                       SU2_MPI::GetComm());

    /*--- Print results on screen. ---*/
    if (rank == MASTER_NODE) {
      for (iMarker_Global = 0; iMarker_Global < nMarker_Global; iMarker_Global++) {
        if (config->GetMarker_CfgFile_KindBC(config->GetMarker_CfgFile_TagBound(iMarker_Global)) == SYMMETRY_PLANE ||
            config->GetMarker_CfgFile_KindBC(config->GetMarker_CfgFile_TagBound(iMarker_Global)) == EULER_WALL) {
          cout << "Boundary marker " << config->GetMarker_CfgFile_TagBound(iMarker_Global) << " is";
          if (!static_cast<bool>(Buff_Recv_isStraight[iMarker_Global])) cout << " NOT";
          if (nDim == 2) cout << " a single straight." << endl;
          if (nDim == 3) cout << " a single plane." << endl;
        }  // if sym or euler
      }    // for iMarker_Global
    }      // if rank==MASTER
  }        // if print_on_scren
}

void CGeometry::ComputeSurf_Curvature(CConfig* config) {
  unsigned short iMarker, iNeigh_Point, iDim, iNode, iNeighbor_Nodes, Neighbor_Node;
  unsigned long Neighbor_Point, iVertex, iPoint, jPoint, iElem_Bound, iEdge, nLocalVertex, MaxLocalVertex,
      *Buffer_Send_nVertex, *Buffer_Receive_nVertex, TotalnPointDomain;
  vector<unsigned long> Point_NeighborList, Elem_NeighborList, Point_Triangle, Point_Edge, Point_Critical;
  su2double U[3] = {0.0}, V[3] = {0.0}, W[3] = {0.0}, Length_U, Length_V, Length_W, CosValue, Angle_Value, *K,
            *Angle_Defect, *Area_Vertex, *Angle_Alpha, *Angle_Beta, **NormalMeanK, MeanK, GaussK, MaxPrinK, cot_alpha,
            cot_beta, delta, X1, X2, X3, Y1, Y2, Y3, radius, *Buffer_Send_Coord, *Buffer_Receive_Coord, *Coord, Dist,
            MinDist, MaxK, MinK, SigmaK;
  bool* Check_Edge;

  /*--- Allocate surface curvature ---*/
  K = new su2double[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) K[iPoint] = 0.0;

  if (nDim == 2) {
    /*--- Loop over all the markers ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        /*--- Loop through all marker vertices again, this time also
         finding the neighbors of each node.---*/
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();

          if (nodes->GetDomain(iPoint)) {
            /*--- Loop through neighbors. In 2-D, there should be 2 nodes on either
             side of this vertex that lie on the same surface. ---*/
            Point_Edge.clear();

            for (iNeigh_Point = 0; iNeigh_Point < nodes->GetnPoint(iPoint); iNeigh_Point++) {
              Neighbor_Point = nodes->GetPoint(iPoint, iNeigh_Point);

              /*--- Check if this neighbor lies on the surface. If so,
               add to the list of neighbors. ---*/
              if (nodes->GetPhysicalBoundary(Neighbor_Point)) {
                Point_Edge.push_back(Neighbor_Point);
              }
            }

            if (Point_Edge.size() == 2) {
              /*--- Compute the curvature using three points ---*/
              X1 = nodes->GetCoord(iPoint, 0);
              X2 = nodes->GetCoord(Point_Edge[0], 0);
              X3 = nodes->GetCoord(Point_Edge[1], 0);
              Y1 = nodes->GetCoord(iPoint, 1);
              Y2 = nodes->GetCoord(Point_Edge[0], 1);
              Y3 = nodes->GetCoord(Point_Edge[1], 1);

              radius = sqrt(((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1)) *
                            ((X2 - X3) * (X2 - X3) + (Y2 - Y3) * (Y2 - Y3)) *
                            ((X3 - X1) * (X3 - X1) + (Y3 - Y1) * (Y3 - Y1))) /
                       (2.0 * fabs(X1 * Y2 + X2 * Y3 + X3 * Y1 - X1 * Y3 - X2 * Y1 - X3 * Y2) + EPS);

              K[iPoint] = 1.0 / radius;
              nodes->SetCurvature(iPoint, K[iPoint]);
            }
          }
        }
      }
    }

  }

  else {
    Angle_Defect = new su2double[nPoint];
    Area_Vertex = new su2double[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      Angle_Defect[iPoint] = 2 * PI_NUMBER;
      Area_Vertex[iPoint] = 0.0;
    }

    Angle_Alpha = new su2double[nEdge];
    Angle_Beta = new su2double[nEdge];
    Check_Edge = new bool[nEdge];
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
      Angle_Alpha[iEdge] = 0.0;
      Angle_Beta[iEdge] = 0.0;
      Check_Edge[iEdge] = true;
    }

    NormalMeanK = new su2double*[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      NormalMeanK[iPoint] = new su2double[nDim];
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

              for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode);
                   iNeighbor_Nodes++) {
                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                Neighbor_Point = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                Point_Triangle.push_back(Neighbor_Point);
              }

              iEdge = FindEdge(Point_Triangle[0], Point_Triangle[1]);

              for (iDim = 0; iDim < nDim; iDim++) {
                U[iDim] = nodes->GetCoord(Point_Triangle[0], iDim) - nodes->GetCoord(iPoint, iDim);
                V[iDim] = nodes->GetCoord(Point_Triangle[1], iDim) - nodes->GetCoord(iPoint, iDim);
              }

              W[0] = 0.5 * (U[1] * V[2] - U[2] * V[1]);
              W[1] = -0.5 * (U[0] * V[2] - U[2] * V[0]);
              W[2] = 0.5 * (U[0] * V[1] - U[1] * V[0]);

              Length_U = 0.0;
              Length_V = 0.0;
              Length_W = 0.0;
              CosValue = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                Length_U += U[iDim] * U[iDim];
                Length_V += V[iDim] * V[iDim];
                Length_W += W[iDim] * W[iDim];
              }
              Length_U = sqrt(Length_U);
              Length_V = sqrt(Length_V);
              Length_W = sqrt(Length_W);
              for (iDim = 0; iDim < nDim; iDim++) {
                U[iDim] /= Length_U;
                V[iDim] /= Length_V;
                CosValue += U[iDim] * V[iDim];
              }
              if (CosValue >= 1.0) CosValue = 1.0;
              if (CosValue <= -1.0) CosValue = -1.0;

              Angle_Value = acos(CosValue);
              Area_Vertex[iPoint] += Length_W;
              Angle_Defect[iPoint] -= Angle_Value;
              if (Angle_Alpha[iEdge] == 0.0)
                Angle_Alpha[iEdge] = Angle_Value;
              else
                Angle_Beta[iEdge] = Angle_Value;
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

              for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode);
                   iNeighbor_Nodes++) {
                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                jPoint = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);

                iEdge = FindEdge(iPoint, jPoint);

                if (Check_Edge[iEdge]) {
                  Check_Edge[iEdge] = false;

                  if (tan(Angle_Alpha[iEdge]) != 0.0)
                    cot_alpha = 1.0 / tan(Angle_Alpha[iEdge]);
                  else
                    cot_alpha = 0.0;
                  if (tan(Angle_Beta[iEdge]) != 0.0)
                    cot_beta = 1.0 / tan(Angle_Beta[iEdge]);
                  else
                    cot_beta = 0.0;

                  /*--- iPoint, and jPoint ---*/
                  for (iDim = 0; iDim < nDim; iDim++) {
                    if (Area_Vertex[iPoint] != 0.0)
                      NormalMeanK[iPoint][iDim] += 3.0 * (cot_alpha + cot_beta) *
                                                   (nodes->GetCoord(iPoint, iDim) - nodes->GetCoord(jPoint, iDim)) /
                                                   Area_Vertex[iPoint];
                    if (Area_Vertex[jPoint] != 0.0)
                      NormalMeanK[jPoint][iDim] += 3.0 * (cot_alpha + cot_beta) *
                                                   (nodes->GetCoord(jPoint, iDim) - nodes->GetCoord(iPoint, iDim)) /
                                                   Area_Vertex[jPoint];
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
          iPoint = vertex[iMarker][iVertex]->GetNode();

          if (nodes->GetDomain(iPoint)) {
            if (Area_Vertex[iPoint] != 0.0)
              GaussK = 3.0 * Angle_Defect[iPoint] / Area_Vertex[iPoint];
            else
              GaussK = 0.0;

            MeanK = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) MeanK += NormalMeanK[iPoint][iDim] * NormalMeanK[iPoint][iDim];
            MeanK = sqrt(MeanK);

            delta = max((MeanK * MeanK - GaussK), 0.0);

            MaxPrinK = MeanK + sqrt(delta);

            /*--- Store the curvature value ---*/
            K[iPoint] = MaxPrinK;
            nodes->SetCurvature(iPoint, K[iPoint]);
          }
        }
      }
    }

    delete[] Angle_Defect;
    delete[] Area_Vertex;
    delete[] Angle_Alpha;
    delete[] Angle_Beta;
    delete[] Check_Edge;

    for (iPoint = 0; iPoint < nPoint; iPoint++) delete[] NormalMeanK[iPoint];
    delete[] NormalMeanK;
  }

  /*--- Sharp edge detection is based in the statistical
   distribution of the curvature ---*/

  MaxK = K[0];
  MinK = K[0];
  MeanK = 0.0;
  TotalnPointDomain = 0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (nodes->GetDomain(iPoint)) {
          MaxK = max(MaxK, fabs(K[iPoint]));
          MinK = min(MinK, fabs(K[iPoint]));
          MeanK += fabs(K[iPoint]);
          TotalnPointDomain++;
        }
      }
    }
  }

  su2double MyMeanK = MeanK;
  MeanK = 0.0;
  su2double MyMaxK = MaxK;
  MaxK = 0.0;
  unsigned long MynPointDomain = TotalnPointDomain;
  TotalnPointDomain = 0;
  SU2_MPI::Allreduce(&MyMeanK, &MeanK, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MyMaxK, &MaxK, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MynPointDomain, &TotalnPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  /*--- Compute the mean ---*/
  MeanK /= su2double(TotalnPointDomain);

  /*--- Compute the standard deviation ---*/
  SigmaK = 0.0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (nodes->GetDomain(iPoint)) {
          SigmaK += (fabs(K[iPoint]) - MeanK) * (fabs(K[iPoint]) - MeanK);
        }
      }
    }
  }

  su2double MySigmaK = SigmaK;
  SigmaK = 0.0;
  SU2_MPI::Allreduce(&MySigmaK, &SigmaK, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  SigmaK = sqrt(SigmaK / su2double(TotalnPointDomain));

  if (rank == MASTER_NODE)
    cout << "Max K: " << MaxK << ". Mean K: " << MeanK << ". Standard deviation K: " << SigmaK << "." << endl;

  Point_Critical.clear();

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (nodes->GetDomain(iPoint)) {
          if (fabs(K[iPoint]) > MeanK + config->GetRefSharpEdges() * SigmaK) {
            Point_Critical.push_back(iPoint);
          }
        }
      }
    }
  }

  /*--- Variables and buffers needed for MPI ---*/

  Buffer_Send_nVertex = new unsigned long[1];
  Buffer_Receive_nVertex = new unsigned long[size];

  /*--- Count the total number of critical edge nodes. ---*/

  nLocalVertex = Point_Critical.size();
  Buffer_Send_nVertex[0] = nLocalVertex;

  /*--- Communicate to all processors the total number of critical edge nodes. ---*/

  MaxLocalVertex = 0;
  SU2_MPI::Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG,
                     SU2_MPI::GetComm());

  /*--- Create and initialize to zero some buffers to hold the coordinates
   of the boundary nodes that are communicated from each partition (all-to-all). ---*/

  const unsigned long nBuffer = MaxLocalVertex * nDim;
  Buffer_Send_Coord = new su2double[nBuffer]();
  Buffer_Receive_Coord = new su2double[size * nBuffer];

  /*--- Retrieve and store the coordinates of the sharp edges boundary nodes on
   the local partition and broadcast them to all partitions. ---*/

  for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
    iPoint = Point_Critical[iVertex];
    for (iDim = 0; iDim < nDim; iDim++) Buffer_Send_Coord[iVertex * nDim + iDim] = nodes->GetCoord(iPoint, iDim);
  }

  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE,
                     SU2_MPI::GetComm());

  /*--- Loop over all interior mesh nodes on the local partition and compute
   the distances to each of the no-slip boundary nodes in the entire mesh.
   Store the minimum distance to the wall for each interior mesh node. ---*/

  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    Coord = nodes->GetCoord(iPoint);

    MinDist = 1E20;
    for (int iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Dist = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Dist += (Coord[iDim] - Buffer_Receive_Coord[(iProcessor * MaxLocalVertex + iVertex) * nDim + iDim]) *
                  (Coord[iDim] - Buffer_Receive_Coord[(iProcessor * MaxLocalVertex + iVertex) * nDim + iDim]);
        }
        if (Dist != 0.0)
          Dist = sqrt(Dist);
        else
          Dist = 0.0;
        if (Dist < MinDist) MinDist = Dist;
      }
    }
    nodes->SetSharpEdge_Distance(iPoint, MinDist);
  }

  /*--- Deallocate Max curvature ---*/
  delete[] K;

  /*--- Deallocate the buffers needed for the MPI communication. ---*/
  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_nVertex;
  delete[] Buffer_Receive_nVertex;
}

void CGeometry::FilterValuesAtElementCG(const vector<su2double>& filter_radius,
                                        const vector<pair<ENUM_FILTER_KERNEL, su2double>>& kernels,
                                        const unsigned short search_limit, su2double* values) const {
  /*--- Apply a filter to "input_values". The filter is an averaging process over the neighbourhood
  of each element, which is a circle in 2D and a sphere in 3D of radius "filter_radius".
  The filter is characterized by its kernel, i.e. how the weights are computed. Multiple kernels
  can be specified in which case they are applied sequentially (each one being applied to the
  output values of the previous filter. ---*/

  /*--- Check if we need to do any work. ---*/
  if (kernels.empty()) return;

  /*--- FIRST: Gather the adjacency matrix, element centroids, volumes, and values on every
  processor, this is required because the filter reaches far into adjacent partitions. ---*/

  /*--- Adjacency matrix ---*/
  vector<unsigned long> neighbour_start;
  long* neighbour_idx = nullptr;
  GetGlobalElementAdjacencyMatrix(neighbour_start, neighbour_idx);

  /*--- Element centroids and volumes. ---*/
  auto *cg_elem = new su2double[Global_nElemDomain * nDim], *vol_elem = new su2double[Global_nElemDomain];
#ifdef HAVE_MPI
  /*--- Number of subdomain each point is part of. ---*/
  vector<char> halo_detect(Global_nElemDomain);
#endif

  /*--- Inputs of a filter stage, like with CG and volumes, each processor needs to see everything. ---*/
  auto* work_values = new su2double[Global_nElemDomain];

  /*--- When gathering the neighborhood of each element we use a vector of booleans to indicate
  whether an element is already added to the list of neighbors (one vector per thread). ---*/
  vector<vector<bool>> is_neighbor(omp_get_max_threads());

  /*--- Begin OpenMP parallel section, count total number of searches for which
  the recursion limit is reached and the full neighborhood is not considered. ---*/
  unsigned long limited_searches = 0;

  SU2_OMP_PARALLEL_(reduction(+ : limited_searches)) {
    /*--- Initialize ---*/
    SU2_OMP_FOR_STAT(256)
    for (auto iElem = 0ul; iElem < Global_nElemDomain; ++iElem) {
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) cg_elem[nDim * iElem + iDim] = 0.0;
      vol_elem[iElem] = 0.0;
    }
    END_SU2_OMP_FOR

    /*--- Populate ---*/
    SU2_OMP_FOR_STAT(256)
    for (auto iElem = 0ul; iElem < nElem; ++iElem) {
      auto iElem_global = elem[iElem]->GetGlobalIndex();
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) cg_elem[nDim * iElem_global + iDim] = elem[iElem]->GetCG(iDim);
      vol_elem[iElem_global] = elem[iElem]->GetVolume();
    }
    END_SU2_OMP_FOR

#ifdef HAVE_MPI
    /*--- Account for the duplication introduced by the halo elements and the
    reduction using MPI_SUM, which is required to maintain differentiabillity. ---*/
    SU2_OMP_FOR_STAT(256)
    for (auto iElem = 0ul; iElem < Global_nElemDomain; ++iElem) halo_detect[iElem] = 0;
    END_SU2_OMP_FOR

    SU2_OMP_FOR_STAT(256)
    for (auto iElem = 0ul; iElem < nElem; ++iElem) halo_detect[elem[iElem]->GetGlobalIndex()] = 1;
    END_SU2_OMP_FOR

    /*--- Share with all processors ---*/
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      auto* dbl_buffer = new su2double[Global_nElemDomain * nDim];
      SU2_MPI::Allreduce(cg_elem, dbl_buffer, Global_nElemDomain * nDim, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      swap(dbl_buffer, cg_elem);
      delete[] dbl_buffer;

      dbl_buffer = new su2double[Global_nElemDomain];
      SU2_MPI::Allreduce(vol_elem, dbl_buffer, Global_nElemDomain, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      swap(dbl_buffer, vol_elem);
      delete[] dbl_buffer;

      vector<char> char_buffer(Global_nElemDomain);
      MPI_Allreduce(halo_detect.data(), char_buffer.data(), Global_nElemDomain, MPI_CHAR, MPI_SUM, SU2_MPI::GetComm());
      halo_detect.swap(char_buffer);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    SU2_OMP_FOR_STAT(256)
    for (auto iElem = 0ul; iElem < Global_nElemDomain; ++iElem) {
      su2double numRepeat = halo_detect[iElem];
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) cg_elem[nDim * iElem + iDim] /= numRepeat;
      vol_elem[iElem] /= numRepeat;
    }
    END_SU2_OMP_FOR
#endif

    /*--- SECOND: Each processor performs the average for its elements. For each
    element we look for neighbours of neighbours of... until the distance to the
    closest newly found one is greater than the filter radius.  ---*/

    is_neighbor[omp_get_thread_num()].resize(Global_nElemDomain, false);

    for (unsigned long iKernel = 0; iKernel < kernels.size(); ++iKernel) {
      auto kernel_type = kernels[iKernel].first;
      su2double kernel_param = kernels[iKernel].second;
      su2double kernel_radius = filter_radius[iKernel];

      /*--- Synchronize work values ---*/
      /*--- Initialize ---*/
      SU2_OMP_FOR_STAT(256)
      for (auto iElem = 0ul; iElem < Global_nElemDomain; ++iElem) work_values[iElem] = 0.0;
      END_SU2_OMP_FOR

      /*--- Populate ---*/
      SU2_OMP_FOR_STAT(256)
      for (auto iElem = 0ul; iElem < nElem; ++iElem) work_values[elem[iElem]->GetGlobalIndex()] = values[iElem];
      END_SU2_OMP_FOR

#ifdef HAVE_MPI
      /*--- Share with all processors ---*/
      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
        auto* buffer = new su2double[Global_nElemDomain];
        SU2_MPI::Allreduce(work_values, buffer, Global_nElemDomain, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
        swap(buffer, work_values);
        delete[] buffer;
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS

      /*--- Account for duplication ---*/
      SU2_OMP_FOR_STAT(256)
      for (auto iElem = 0ul; iElem < Global_nElemDomain; ++iElem) {
        su2double numRepeat = halo_detect[iElem];
        work_values[iElem] /= numRepeat;
      }
      END_SU2_OMP_FOR
#endif

      /*--- Filter ---*/
      SU2_OMP_FOR_DYN(128)
      for (auto iElem = 0ul; iElem < nElem; ++iElem) {
        int thread = omp_get_thread_num();

        /*--- Center of the search ---*/
        auto iElem_global = elem[iElem]->GetGlobalIndex();

        /*--- Find the neighbours of iElem ---*/
        vector<long> neighbours;
        limited_searches +=
            !GetRadialNeighbourhood(iElem_global, SU2_TYPE::GetValue(kernel_radius), search_limit, neighbour_start,
                                    neighbour_idx, cg_elem, neighbours, is_neighbor[thread]);
        /*--- Apply the kernel ---*/
        su2double weight = 0.0, numerator = 0.0, denominator = 0.0;

        switch (kernel_type) {
          /*--- distance-based kernels (weighted averages) ---*/
          case ENUM_FILTER_KERNEL::CONSTANT_WEIGHT:
          case ENUM_FILTER_KERNEL::CONICAL_WEIGHT:
          case ENUM_FILTER_KERNEL::GAUSSIAN_WEIGHT:

            for (auto idx : neighbours) {
              su2double distance = 0.0;
              for (unsigned short iDim = 0; iDim < nDim; ++iDim)
                distance += pow(cg_elem[nDim * iElem_global + iDim] - cg_elem[nDim * idx + iDim], 2);
              distance = sqrt(distance);

              switch (kernel_type) {
                case ENUM_FILTER_KERNEL::CONSTANT_WEIGHT:
                  weight = 1.0;
                  break;
                case ENUM_FILTER_KERNEL::CONICAL_WEIGHT:
                  weight = kernel_radius - distance;
                  break;
                case ENUM_FILTER_KERNEL::GAUSSIAN_WEIGHT:
                  weight = exp(-0.5 * pow(distance / kernel_param, 2));
                  break;
                default:
                  break;
              }
              weight *= vol_elem[idx];
              numerator += weight * work_values[idx];
              denominator += weight;
            }
            values[iElem] = numerator / denominator;
            break;

          /*--- morphology kernels (image processing) ---*/
          case ENUM_FILTER_KERNEL::DILATE_MORPH:
          case ENUM_FILTER_KERNEL::ERODE_MORPH:

            for (auto idx : neighbours) {
              switch (kernel_type) {
                case ENUM_FILTER_KERNEL::DILATE_MORPH:
                  numerator += exp(kernel_param * work_values[idx]);
                  break;
                case ENUM_FILTER_KERNEL::ERODE_MORPH:
                  numerator += exp(kernel_param * (1.0 - work_values[idx]));
                  break;
                default:
                  break;
              }
              denominator += 1.0;
            }
            values[iElem] = log(numerator / denominator) / kernel_param;
            if (kernel_type == ENUM_FILTER_KERNEL::ERODE_MORPH) values[iElem] = 1.0 - values[iElem];
            break;

          default:
            SU2_MPI::Error("Unknown type of filter kernel", CURRENT_FUNCTION);
        }
      }
      END_SU2_OMP_FOR
    }
  }
  END_SU2_OMP_PARALLEL

  limited_searches /= kernels.size();

  unsigned long tmp = limited_searches;
  SU2_MPI::Reduce(&tmp, &limited_searches, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());

  if (rank == MASTER_NODE && limited_searches > 0)
    cout << "Warning: The filter radius was limited for " << limited_searches << " elements ("
         << limited_searches / (0.01 * Global_nElemDomain) << "%).\n";

  delete[] neighbour_idx;
  delete[] cg_elem;
  delete[] vol_elem;
  delete[] work_values;
}

void CGeometry::GetGlobalElementAdjacencyMatrix(vector<unsigned long>& neighbour_start, long*& neighbour_idx) const {
  if (neighbour_idx != nullptr)
    SU2_MPI::Error("neighbour_idx is expected to be NULL, stopping to avoid a potential memory leak", CURRENT_FUNCTION);

  /*--- Determine how much space we need for the adjacency matrix by counting the
  neighbours of each element, i.e. its number of faces---*/
  auto* nFaces_elem = new unsigned short[Global_nElemDomain];

  SU2_OMP_PARALLEL {
    SU2_OMP_FOR_STAT(256)
    for (auto iElem = 0ul; iElem < Global_nElemDomain; ++iElem) nFaces_elem[iElem] = 0;
    END_SU2_OMP_FOR

    SU2_OMP_FOR_STAT(256)
    for (auto iElem = 0ul; iElem < nElem; ++iElem) {
      auto iElem_global = elem[iElem]->GetGlobalIndex();
      nFaces_elem[iElem_global] = elem[iElem]->GetnFaces();
    }
    END_SU2_OMP_FOR
  }
  END_SU2_OMP_PARALLEL
#ifdef HAVE_MPI
  /*--- Share with all processors ---*/
  {
    auto* buffer = new unsigned short[Global_nElemDomain];
    MPI_Allreduce(nFaces_elem, buffer, Global_nElemDomain, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
    /*--- swap pointers and delete old data to keep the same variable name after reduction ---*/
    swap(buffer, nFaces_elem);
    delete[] buffer;
  }
#endif

  /*--- Vector with the addresses of the start of the neighbours of a given element.
  This is generated by a cumulative sum of the neighbour count. ---*/
  neighbour_start.resize(Global_nElemDomain + 1);

  neighbour_start[0] = 0;
  for (auto iElem = 0ul; iElem < Global_nElemDomain; ++iElem) {
    neighbour_start[iElem + 1] = neighbour_start[iElem] + nFaces_elem[iElem];
  }
  delete[] nFaces_elem;

  /*--- Allocate ---*/
  unsigned long matrix_size = neighbour_start[Global_nElemDomain];
  neighbour_idx = new long[matrix_size];

  SU2_OMP_PARALLEL {
    /*--- Initialize ---*/
    SU2_OMP_FOR_STAT(256)
    for (auto iElem = 0ul; iElem < matrix_size; ++iElem) neighbour_idx[iElem] = -1;
    END_SU2_OMP_FOR

    /*--- Populate ---*/
    SU2_OMP_FOR_STAT(128)
    for (auto iElem = 0ul; iElem < nElem; ++iElem) {
      auto iElem_global = elem[iElem]->GetGlobalIndex();
      auto start_pos = neighbour_start[iElem_global];

      for (unsigned short iFace = 0; iFace < elem[iElem]->GetnFaces(); ++iFace) {
        long neighbour = elem[iElem]->GetNeighbor_Elements(iFace);

        if (neighbour >= 0) {
          neighbour_idx[start_pos + iFace] = elem[neighbour]->GetGlobalIndex();
        }
      }
    }
    END_SU2_OMP_FOR
  }
  END_SU2_OMP_PARALLEL
#ifdef HAVE_MPI
  /*--- Share with all processors ---*/
  {
    long* buffer = new long[matrix_size];
    MPI_Allreduce(neighbour_idx, buffer, matrix_size, MPI_LONG, MPI_MAX, SU2_MPI::GetComm());
    swap(buffer, neighbour_idx);
    delete[] buffer;
  }
#endif
}

bool CGeometry::GetRadialNeighbourhood(const unsigned long iElem_global, const passivedouble radius,
                                       size_t search_limit, const vector<unsigned long>& neighbour_start,
                                       const long* neighbour_idx, const su2double* cg_elem, vector<long>& neighbours,
                                       vector<bool>& is_neighbor) const {
  /*--- Validate inputs if we are debugging. ---*/
  assert(neighbour_start.size() == Global_nElemDomain + 1 && neighbour_idx != nullptr && cg_elem != nullptr &&
         is_neighbor.size() == Global_nElemDomain && "invalid inputs");

  /*--- 0 search_limit means "unlimited" (it will probably
   stop once it gathers the entire domain, probably). ---*/
  if (!search_limit) search_limit = numeric_limits<size_t>::max();

  /*--- Center of the search ---*/
  neighbours.clear();
  neighbours.push_back(iElem_global);
  is_neighbor[iElem_global] = true;

  passivedouble X0[3] = {0.0, 0.0, 0.0};
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) X0[iDim] = SU2_TYPE::GetValue(cg_elem[nDim * iElem_global + iDim]);

  /*--- Loop stops when "neighbours" stops changing size, or degree reaches limit. ---*/
  bool finished = false;
  for (size_t degree = 0, start = 0; degree < search_limit && !finished; ++degree) {
    /*--- For each element of the last degree added consider its immediate
     neighbours, that are not already neighbours, as candidates. ---*/
    vector<long> candidates;

    for (auto it = neighbours.begin() + start; it != neighbours.end(); ++it) {
      /*--- scan row of the adjacency matrix of element *it ---*/
      for (auto i = neighbour_start[*it]; i < neighbour_start[(*it) + 1]; ++i) {
        auto idx = neighbour_idx[i];
        if (idx >= 0)
          if (!is_neighbor[idx]) {
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
    for (auto idx : candidates) {
      /*--- passivedouble as we only need to compare "distance". ---*/
      passivedouble distance = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; ++iDim)
        distance += pow(X0[iDim] - SU2_TYPE::GetValue(cg_elem[nDim * idx + iDim]), 2);

      if (distance < pow(radius, 2)) {
        neighbours.push_back(idx);
        finished = false;
      }
      /*--- not a neighbour in the end. ---*/
      else
        is_neighbor[idx] = false;
    }
  }
  /*--- Restore the state of the working vector for next call. ---*/
  for (auto idx : neighbours) is_neighbor[idx] = false;

  return finished;
}

void CGeometry::SetElemVolume() {
  SU2_OMP_PARALLEL {
    /*--- Create a bank of elements to avoid instantiating inside loop. ---*/
    CElement* elements[4] = {nullptr, nullptr, nullptr, nullptr};

    if (nDim == 2) {
      elements[0] = new CTRIA1();
      elements[1] = new CQUAD4();
    } else {
      elements[0] = new CTETRA1();
      elements[1] = new CPYRAM5();
      elements[2] = new CPRISM6();
      elements[3] = new CHEXA8();
    }

    /*--- Compute and store the volume of each "elem". ---*/
    SU2_OMP_FOR_DYN(128)
    for (unsigned long iElem = 0; iElem < nElem; ++iElem) {
      /*--- Get the appropriate type of element. ---*/
      CElement* element = nullptr;
      switch (elem[iElem]->GetVTK_Type()) {
        case TRIANGLE:
          element = elements[0];
          break;
        case QUADRILATERAL:
          element = elements[1];
          break;
        case TETRAHEDRON:
          element = elements[0];
          break;
        case PYRAMID:
          element = elements[1];
          break;
        case PRISM:
          element = elements[2];
          break;
        case HEXAHEDRON:
          element = elements[3];
          break;
        default:
          SU2_MPI::Error("Cannot compute the area/volume of a 1D element.", CURRENT_FUNCTION);
      }
      /*--- Set the nodal coordinates of the element. ---*/
      for (unsigned short iNode = 0; iNode < elem[iElem]->GetnNodes(); ++iNode) {
        unsigned long node_idx = elem[iElem]->GetNode(iNode);
        for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
          su2double coord = nodes->GetCoord(node_idx, iDim);
          element->SetRef_Coord(iNode, iDim, coord);
        }
      }
      /*--- Compute ---*/
      if (nDim == 2)
        elem[iElem]->SetVolume(element->ComputeArea());
      else
        elem[iElem]->SetVolume(element->ComputeVolume());
    }
    END_SU2_OMP_FOR

    delete elements[0];
    delete elements[1];
    if (nDim == 3) {
      delete elements[2];
      delete elements[3];
    }
  }
  END_SU2_OMP_PARALLEL
}

void CGeometry::SetRotationalVelocity(const CConfig* config, bool print) {
  unsigned long iPoint;
  unsigned short iDim;

  su2double GridVel[3] = {0.0, 0.0, 0.0}, Distance[3] = {0.0, 0.0, 0.0}, Center[3] = {0.0, 0.0, 0.0},
            Omega[3] = {0.0, 0.0, 0.0}, xDot[3] = {0.0, 0.0, 0.0};

  /*--- Center of rotation & angular velocity vector from config ---*/

  for (iDim = 0; iDim < 3; iDim++) {
    Center[iDim] = config->GetMotion_Origin(iDim);
    Omega[iDim] = config->GetRotation_Rate(iDim) / config->GetOmega_Ref();
    xDot[iDim] = config->GetTranslation_Rate(iDim) / config->GetVelocity_Ref();
  }

  su2double L_Ref = config->GetLength_Ref();

  /*--- Print some information to the console ---*/

  if (rank == MASTER_NODE && print) {
    cout << " Rotational origin (x, y, z): ( " << Center[0] << ", " << Center[1];
    cout << ", " << Center[2] << " )\n";
    cout << " Angular velocity about x, y, z axes: ( " << Omega[0] << ", ";
    cout << Omega[1] << ", " << Omega[2] << " ) rad/s" << endl;
    cout << " Translational velocity in x, y, z direction: (" << xDot[0] << ", " << xDot[1] << ", " << xDot[2] << ")."
         << endl;
  }

  /*--- Loop over all nodes and set the rotational and translational velocity ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    /*--- Get the coordinates of the current node ---*/

    const su2double* Coord = nodes->GetCoord(iPoint);

    /*--- Calculate the non-dim. distance from the rotation center ---*/

    for (iDim = 0; iDim < nDim; iDim++) Distance[iDim] = (Coord[iDim] - Center[iDim]) / L_Ref;

    /*--- Calculate the angular velocity as omega X r and add translational velocity ---*/

    GridVel[0] = Omega[1] * (Distance[2]) - Omega[2] * (Distance[1]) + xDot[0];
    GridVel[1] = Omega[2] * (Distance[0]) - Omega[0] * (Distance[2]) + xDot[1];
    GridVel[2] = Omega[0] * (Distance[1]) - Omega[1] * (Distance[0]) + xDot[2];

    /*--- Store the grid velocity at this node ---*/

    nodes->SetGridVel(iPoint, GridVel);
  }
}

void CGeometry::SetShroudVelocity(const CConfig* config) {
  unsigned long iPoint, iVertex;
  unsigned short iMarker, iMarkerShroud;
  su2double RotVel[3] = {0.0, 0.0, 0.0};

  /*--- Loop over all vertex in the shroud marker and set the rotational velocity to 0.0 ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iMarkerShroud = 0; iMarkerShroud < config->GetnMarker_Shroud(); iMarkerShroud++) {
      if (config->GetMarker_Shroud(iMarkerShroud) == config->GetMarker_All_TagBound(iMarker)) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          nodes->SetGridVel(iPoint, RotVel);
        }
      }
    }
  }
}

void CGeometry::SetTranslationalVelocity(const CConfig* config, bool print) {
  su2double xDot[3] = {0.0, 0.0, 0.0};

  /*--- Get the translational velocity vector from config ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    xDot[iDim] = config->GetTranslation_Rate(iDim) / config->GetVelocity_Ref();

  /*--- Print some information to the console ---*/

  if (rank == MASTER_NODE && print) {
    cout << " Non-dim. translational velocity: (" << xDot[0] << ", " << xDot[1] << ", " << xDot[2] << ")." << endl;
  }

  /*--- Loop over all nodes and set the translational velocity ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) nodes->SetGridVel(iPoint, xDot);
}

void CGeometry::SetWallVelocity(const CConfig* config, bool print) {
  const su2double L_Ref = config->GetLength_Ref();
  const su2double Omega_Ref = config->GetOmega_Ref();
  const su2double Vel_Ref = config->GetVelocity_Ref();

  /*--- Store grid velocity for each node on the moving surface(s).
   Sum and store the x, y, & z velocities due to translation and rotation. ---*/

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Moving(iMarker) != YES) continue;

    /*--- Identify iMarker from the list of those under MARKER_MOVING
     Get prescribed wall speed from config for this marker. ---*/

    const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    const auto jMarker = config->GetMarker_Moving(Marker_Tag);

    su2double xDot[MAXNDIM], Center[MAXNDIM], Omega[MAXNDIM];

    for (auto iDim = 0u; iDim < MAXNDIM; iDim++) {
      Center[iDim] = config->GetMarkerMotion_Origin(jMarker, iDim);
      Omega[iDim] = config->GetMarkerRotationRate(jMarker, iDim) / Omega_Ref;
      xDot[iDim] = config->GetMarkerTranslationRate(jMarker, iDim) / Vel_Ref;
    }

    if (rank == MASTER_NODE && print) {
      cout << " Storing grid velocity for marker: ";
      cout << Marker_Tag << ".\n";
      cout << " Translational velocity: (" << xDot[0] * config->GetVelocity_Ref() << ", "
           << xDot[1] * config->GetVelocity_Ref();
      cout << ", " << xDot[2] * config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == SI)
        cout << ") m/s.\n";
      else
        cout << ") ft/s.\n";
      cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
      cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
      cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
    }

    for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) {
      const auto iPoint = vertex[iMarker][iVertex]->GetNode();

      /*--- Calculate non-dim. position from rotation center ---*/
      su2double r[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++) r[iDim] = (nodes->GetCoord(iPoint, iDim) - Center[iDim]) / L_Ref;

      /*--- Cross Product of angular velocity and distance from center to
       get the rotational velocity. Note that we are adding on the velocity
       due to pure translation as well. ---*/

      su2double GridVel[MAXNDIM];
      GeometryToolbox::CrossProduct(Omega, r, GridVel);

      for (auto iDim = 0u; iDim < nDim; iDim++) nodes->SetGridVel(iPoint, iDim, xDot[iDim] + GridVel[iDim]);
    }
  }
}

void CGeometry::SetGridVelocity(const CConfig* config) {
  /*--- Get timestep and whether to use 1st or 2nd order backward finite differences ---*/

  bool FirstOrder = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
  bool SecondOrder = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);

  su2double TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the velocity of each node in the volume mesh ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    /*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/

    const su2double* Coord_nM1 = nodes->GetCoord_n1(iPoint);
    const su2double* Coord_n = nodes->GetCoord_n(iPoint);
    const su2double* Coord_nP1 = nodes->GetCoord(iPoint);

    /*--- Compute and store mesh velocity with 1st or 2nd-order approximation ---*/

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      su2double GridVel = 0.0;

      if (FirstOrder) GridVel = (Coord_nP1[iDim] - Coord_n[iDim]) / TimeStep;

      if (SecondOrder) GridVel = (1.5 * Coord_nP1[iDim] - 2.0 * Coord_n[iDim] + 0.5 * Coord_nM1[iDim]) / TimeStep;

      nodes->SetGridVel(iPoint, iDim, GridVel);
    }
  }
}

const CCompressedSparsePatternUL& CGeometry::GetSparsePattern(ConnectivityType type, unsigned long fillLvl) {
  bool fvm = (type == ConnectivityType::FiniteVolume);

  CCompressedSparsePatternUL* pattern = nullptr;

  if (fillLvl == 0)
    pattern = fvm ? &finiteVolumeCSRFill0 : &finiteElementCSRFill0;
  else
    pattern = fvm ? &finiteVolumeCSRFillN : &finiteElementCSRFillN;

  if (pattern->empty()) {
    *pattern = buildCSRPattern(*this, type, fillLvl);
    pattern->buildDiagPtr();
  }

  return *pattern;
}

const CEdgeToNonZeroMapUL& CGeometry::GetEdgeToSparsePatternMap() {
  if (edgeToCSRMap.empty()) {
    if (finiteVolumeCSRFill0.empty()) {
      finiteVolumeCSRFill0 = buildCSRPattern(*this, ConnectivityType::FiniteVolume, 0ul);
    }
    edgeToCSRMap = mapEdgesToSparsePattern(*this, finiteVolumeCSRFill0);
  }
  return edgeToCSRMap;
}

const su2vector<unsigned long>& CGeometry::GetTransposeSparsePatternMap(ConnectivityType type) {
  /*--- Yes the const cast is weird but it is still better than repeating code. ---*/
  auto& pattern = const_cast<CCompressedSparsePatternUL&>(GetSparsePattern(type));
  pattern.buildTransposePtr();
  return pattern.transposePtr();
}

const CCompressedSparsePatternUL& CGeometry::GetEdgeColoring(su2double* efficiency) {
  /*--- Check for dry run mode with dummy geometry. ---*/
  if (nEdge == 0) return edgeColoring;

  /*--- Build if required. ---*/
  if (edgeColoring.empty()) {
    /*--- When not using threading use the natural coloring to reduce overhead. ---*/
    if (omp_get_max_threads() == 1) {
      SetNaturalEdgeColoring();
      if (efficiency != nullptr) *efficiency = 1.0;  // by definition
      return edgeColoring;
    }

    /*--- Create a temporary sparse pattern from the edges. ---*/
    su2vector<unsigned long> outerPtr(nEdge + 1);
    su2vector<unsigned long> innerIdx(nEdge * 2);

    for (unsigned long iEdge = 0; iEdge < nEdge; ++iEdge) {
      outerPtr(iEdge) = 2 * iEdge;
      innerIdx(iEdge * 2 + 0) = edges->GetNode(iEdge, 0);
      innerIdx(iEdge * 2 + 1) = edges->GetNode(iEdge, 1);
    }
    outerPtr(nEdge) = 2 * nEdge;

    CCompressedSparsePatternUL pattern(move(outerPtr), move(innerIdx));

    /*--- Color the edges. ---*/
    constexpr bool balanceColors = true;
    edgeColoring = colorSparsePattern(pattern, edgeColorGroupSize, balanceColors);

    /*--- If the coloring fails use the natural coloring. This is a
     *    "soft" failure as this "bad" coloring should be detected
     *    downstream and a fallback strategy put in place. ---*/
    if (edgeColoring.empty()) SetNaturalEdgeColoring();
  }

  if (efficiency != nullptr) {
    *efficiency = coloringEfficiency(edgeColoring, omp_get_max_threads(), edgeColorGroupSize);
  }
  return edgeColoring;
}

void CGeometry::SetNaturalEdgeColoring() {
  if (nEdge == 0) return;
  edgeColoring = createNaturalColoring(nEdge);
  /*--- In parallel, set the group size to nEdge to protect client code. ---*/
  if (omp_get_max_threads() > 1) edgeColorGroupSize = nEdge;
}

const CCompressedSparsePatternUL& CGeometry::GetElementColoring(su2double* efficiency) {
  /*--- Check for dry run mode with dummy geometry. ---*/
  if (nElem == 0) return elemColoring;

  /*--- Build if required. ---*/
  if (elemColoring.empty()) {
    /*--- When not using threading use the natural coloring. ---*/
    if (omp_get_max_threads() == 1) {
      SetNaturalElementColoring();
      if (efficiency != nullptr) *efficiency = 1.0;  // by definition
      return elemColoring;
    }

    /*--- Create a temporary sparse pattern from the elements. ---*/
    vector<unsigned long> outerPtr(nElem + 1);
    vector<unsigned long> innerIdx;
    innerIdx.reserve(nElem);

    for (unsigned long iElem = 0; iElem < nElem; ++iElem) {
      outerPtr[iElem] = innerIdx.size();

      for (unsigned short iNode = 0; iNode < elem[iElem]->GetnNodes(); ++iNode) {
        innerIdx.push_back(elem[iElem]->GetNode(iNode));
      }
    }
    outerPtr[nElem] = innerIdx.size();

    CCompressedSparsePatternUL pattern(outerPtr, innerIdx);

    /*--- Color the elements. ---*/
    constexpr bool balanceColors = true;
    elemColoring = colorSparsePattern(pattern, elemColorGroupSize, balanceColors);

    /*--- Same as for the edge coloring. ---*/
    if (elemColoring.empty()) SetNaturalElementColoring();
  }

  if (efficiency != nullptr) {
    *efficiency = coloringEfficiency(elemColoring, omp_get_max_threads(), elemColorGroupSize);
  }
  return elemColoring;
}

void CGeometry::SetNaturalElementColoring() {
  if (nElem == 0) return;
  elemColoring = createNaturalColoring(nElem);
  /*--- In parallel, set the group size to nElem to protect client code. ---*/
  if (omp_get_max_threads() > 1) elemColorGroupSize = nElem;
}

void CGeometry::ColorMGLevels(unsigned short nMGLevels, const CGeometry* const* geometry) {
  using tColor = uint8_t;
  constexpr auto nColor = numeric_limits<tColor>::max();

  if (nMGLevels) CoarseGridColor_.resize(nPoint, nMGLevels) = 0;

  for (auto iMesh = nMGLevels; iMesh >= 1; --iMesh) {
    /*--- Color the coarse points. ---*/
    vector<tColor> color;
    const auto& adjacency = geometry[iMesh]->nodes->GetPoints();
    if (colorSparsePattern<tColor, nColor>(adjacency, 1, false, &color).empty()) continue;

    /*--- Propagate colors to fine mesh. ---*/
    for (auto step = 0u; step < iMesh; ++step) {
      auto coarseMesh = geometry[iMesh - 1 - step];
      if (step)
        for (auto iPoint = 0ul; iPoint < coarseMesh->GetnPoint(); ++iPoint)
          CoarseGridColor_(iPoint, step) = CoarseGridColor_(coarseMesh->nodes->GetParent_CV(iPoint), step - 1);
      else
        for (auto iPoint = 0ul; iPoint < coarseMesh->GetnPoint(); ++iPoint)
          CoarseGridColor_(iPoint, step) = color[coarseMesh->nodes->GetParent_CV(iPoint)];
    }
  }
}

const CGeometry::CLineletInfo& CGeometry::GetLineletInfo(const CConfig* config) const {
  auto& li = lineletInfo;
  if (!li.linelets.empty() || nPoint == 0) return li;

  li.lineletIdx.resize(nPoint, CLineletInfo::NO_LINELET);

  /*--- Estimate number of linelets. ---*/

  unsigned long nLinelet = 0;
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetSolid_Wall(iMarker) || config->GetMarker_All_KindBC(iMarker) == DISPLACEMENT_BOUNDARY) {
      nLinelet += nVertex[iMarker];
    }
  }

  /*--- If the domain contains well defined Linelets ---*/

  unsigned long maxNPoints = 0, sumNPoints = 0;

  if (nLinelet != 0) {
    /*--- Define the basic linelets, starting from each vertex, preventing duplication of points. ---*/

    li.linelets.resize(nLinelet);
    nLinelet = 0;
    for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetSolid_Wall(iMarker) || config->GetMarker_All_KindBC(iMarker) == DISPLACEMENT_BOUNDARY) {
        for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) {
          const auto iPoint = vertex[iMarker][iVertex]->GetNode();
          if (li.lineletIdx[iPoint] == CLineletInfo::NO_LINELET && nodes->GetDomain(iPoint)) {
            li.linelets[nLinelet].push_back(iPoint);
            li.lineletIdx[iPoint] = nLinelet;
            ++nLinelet;
          }
        }
      }
    }
    li.linelets.resize(nLinelet);

    /*--- Create the linelet structure. ---*/

    nLinelet = 0;
    for (auto& linelet : li.linelets) {
      while (linelet.size() < CLineletInfo::MAX_LINELET_POINTS) {
        const auto iPoint = linelet.back();

        /*--- Compute the value of the max and min weights to detect if this region is isotropic. ---*/

        su2double max_weight = 0.0, min_weight = std::numeric_limits<su2double>::max();
        for (auto iNode = 0u; iNode < nodes->GetnPoint(iPoint); iNode++) {
          const auto jPoint = nodes->GetPoint(iPoint, iNode);
          const auto iEdge = nodes->GetEdge(iPoint, iNode);
          const auto* normal = edges->GetNormal(iEdge);
          const su2double area = GeometryToolbox::Norm(nDim, normal);
          const su2double volume_iPoint = nodes->GetVolume(iPoint);
          const su2double volume_jPoint = nodes->GetVolume(jPoint);
          const su2double weight = 0.5 * area * (1.0 / volume_iPoint + 1.0 / volume_jPoint);
          max_weight = max(max_weight, weight);
          min_weight = min(min_weight, weight);
        }

        /*--- Isotropic, stop this linelet. ---*/
        if (min_weight / max_weight > CLineletInfo::ALPHA_ISOTROPIC()) break;

        /*--- Otherwise, add the closest valid neighbor. ---*/

        su2double min_dist2 = std::numeric_limits<su2double>::max();
        auto next_Point = iPoint;
        const auto* iCoord = nodes->GetCoord(iPoint);

        for (const auto jPoint : nodes->GetPoints(iPoint)) {
          if (li.lineletIdx[jPoint] == CLineletInfo::NO_LINELET && nodes->GetDomain(jPoint)) {
            const auto* jCoord = nodes->GetCoord(jPoint);
            const su2double d2 = GeometryToolbox::SquaredDistance(nDim, iCoord, jCoord);
            su2double cosTheta = 1;
            if (linelet.size() > 1) {
              const auto* kCoord = nodes->GetCoord(linelet[linelet.size() - 2]);
              su2double dij[3] = {0.0}, dki[3] = {0.0};
              GeometryToolbox::Distance(nDim, iCoord, kCoord, dki);
              GeometryToolbox::Distance(nDim, jCoord, iCoord, dij);
              cosTheta = GeometryToolbox::DotProduct(3, dki, dij) / sqrt(d2 * GeometryToolbox::SquaredNorm(nDim, dki));
            }
            if (d2 < min_dist2 && cosTheta > 0.7071) {
              next_Point = jPoint;
              min_dist2 = d2;
            }
          }
        }

        /*--- Did not find a suitable point. ---*/
        if (next_Point == iPoint) break;

        linelet.push_back(next_Point);
        li.lineletIdx[next_Point] = nLinelet;
      }
      ++nLinelet;

      maxNPoints = max<unsigned long>(maxNPoints, linelet.size());
      sumNPoints += linelet.size();
    }
  }

  /*--- Average linelet size over all ranks. ---*/

  unsigned long globalNPoints, globalNLineLets;
  SU2_MPI::Allreduce(&sumNPoints, &globalNPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&nLinelet, &globalNLineLets, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  if (rank == MASTER_NODE) {
    std::cout << "Computed linelet structure, "
              << static_cast<unsigned long>(passivedouble(globalNPoints) / globalNLineLets)
              << " points in each line (average)." << std::endl;
  }

  /*--- Color the linelets for OpenMP parallelization and visualization. ---*/

  /*--- Adjacency between lines, computed from point neighbors. ---*/
  std::vector<std::vector<unsigned long>> adjacency(nLinelet);
  for (auto iLine = 0ul; iLine < nLinelet; ++iLine) {
    std::unordered_set<unsigned long> neighbors;
    for (const auto iPoint : li.linelets[iLine]) {
      neighbors.insert(iPoint);
      for (const auto jPoint : nodes->GetPoints(iPoint)) {
        neighbors.insert(jPoint);
      }
    }
    adjacency[iLine].reserve(neighbors.size());
    for (const auto iPoint : neighbors) {
      adjacency[iLine].push_back(iPoint);
    }
  }

  const auto coloring =
      colorSparsePattern<uint8_t, std::numeric_limits<uint8_t>::max()>(CCompressedSparsePatternUL(adjacency), 1, true);
  const auto nColors = coloring.getOuterSize();

  /*--- Sort linelets by color. ---*/
  std::vector<std::vector<unsigned long>> sortedLinelets;
  sortedLinelets.reserve(nLinelet);
  li.colorOffsets.reserve(nColors + 1);
  li.colorOffsets.push_back(0);

  for (auto iColor = 0ul; iColor < nColors; ++iColor) {
    for (const auto iLine : coloring.getInnerIter(iColor)) {
      /*--- Store the new linelet index for its points. ---*/
      for (const auto iPoint : li.linelets[iLine]) {
        li.lineletIdx[iPoint] = sortedLinelets.size();
      }
      sortedLinelets.push_back(std::move(li.linelets[iLine]));
    }
    li.colorOffsets.push_back(sortedLinelets.size());
  }
  li.linelets = std::move(sortedLinelets);

  /*--- For visualization, offset colors to avoid coloring across ranks. ---*/
  std::vector<unsigned long> allNColors(size);
  SU2_MPI::Allgather(&nColors, 1, MPI_UNSIGNED_LONG, allNColors.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
  unsigned long offset = 0;
  for (int i = 0; i < rank; ++i) offset += allNColors[i];

  /*--- Finally, transfer colors to points, using "0" as "no linelet". ---*/
  li.lineletColor.resize(nPoint, 0);
  for (auto iColor = 0ul; iColor < nColors; ++iColor) {
    for (auto iLine = li.colorOffsets[iColor]; iLine < li.colorOffsets[iColor + 1]; ++iLine) {
      for (const auto iPoint : li.linelets[iLine]) {
        li.lineletColor[iPoint] = 1 + offset + iColor;
      }
    }
  }

  return li;
}

void CGeometry::ComputeWallDistance(const CConfig* const* config_container, CGeometry**** geometry_container) {
  int nZone = config_container[ZONE_0]->GetnZone();
  bool allEmpty = true;
  vector<bool> wallDistanceNeeded(nZone, false);

  for (int iInst = 0; iInst < config_container[ZONE_0]->GetnTimeInstances(); iInst++) {
    for (int iZone = 0; iZone < nZone; iZone++) {
      /*--- Check if a zone needs the wall distance and store a boolean ---*/

      MAIN_SOLVER kindSolver = config_container[iZone]->GetKind_Solver();
      if (kindSolver == MAIN_SOLVER::RANS || kindSolver == MAIN_SOLVER::INC_RANS ||
          kindSolver == MAIN_SOLVER::DISC_ADJ_RANS || kindSolver == MAIN_SOLVER::DISC_ADJ_INC_RANS ||
          kindSolver == MAIN_SOLVER::FEM_LES || kindSolver == MAIN_SOLVER::FEM_RANS) {
        wallDistanceNeeded[iZone] = true;
      }

      /*--- Set the wall distances in all zones to the numerical limit.
       * This is necessary, because before a computed distance is set, it will be checked
       * whether the new distance is smaller than the currently stored one. ---*/
      CGeometry* geometry = geometry_container[iZone][iInst][MESH_0];
      if (wallDistanceNeeded[iZone]) geometry->SetWallDistance(numeric_limits<su2double>::max());
    }

    /*--- Loop over all zones and compute the ADT based on the viscous walls in that zone ---*/
    for (int iZone = 0; iZone < nZone; iZone++) {
      unique_ptr<CADTElemClass> WallADT =
          geometry_container[iZone][iInst][MESH_0]->ComputeViscousWallADT(config_container[iZone]);
      if (WallADT && !WallADT->IsEmpty()) {
        allEmpty = false;
        /*--- Inner loop over all zones to update the wall distances.
         * It might happen that there is a closer viscous wall in zone iZone for points in zone jZone. ---*/
        for (int jZone = 0; jZone < nZone; jZone++) {
          if (wallDistanceNeeded[jZone])
            geometry_container[jZone][iInst][MESH_0]->SetWallDistance(WallADT.get(), config_container[jZone], iZone);
        }
      }
    }

    /*--- If there are no viscous walls in the entire domain, set distances to zero ---*/
    if (allEmpty) {
      for (int iZone = 0; iZone < nZone; iZone++) {
        CGeometry* geometry = geometry_container[iZone][iInst][MESH_0];
        geometry->SetWallDistance(0.0);
      }
    }
    /*--- Otherwise, set wall roughnesses. ---*/
    if (!allEmpty) {
      /*--- Store all wall roughnesses in a common data structure. ---*/
      // [iZone][iMarker] -> roughness, for this rank
      auto roughness_f = make_pair(nZone, [config_container, geometry_container, iInst](unsigned long iZone) {
        const CConfig* config = config_container[iZone];
        const auto nMarker = geometry_container[iZone][iInst][MESH_0]->GetnMarker();

        return make_pair(nMarker, [config](unsigned long iMarker) {
          return config->GetWallRoughnessProperties(config->GetMarker_All_TagBound(iMarker)).second;
        });
      });
      NdFlattener<2> roughness_local(roughness_f);
      // [rank][iZone][iMarker] -> roughness
      NdFlattener<3> roughness_global(Nd_MPI_Environment(), roughness_local);
      // use it to update roughnesses
      for (int jZone = 0; jZone < nZone; jZone++) {
        if (wallDistanceNeeded[jZone] && config_container[jZone]->GetnRoughWall() > 0) {
          geometry_container[jZone][iInst][MESH_0]->nodes->SetWallRoughness(roughness_global);
        }
      }
    }
  }
}
