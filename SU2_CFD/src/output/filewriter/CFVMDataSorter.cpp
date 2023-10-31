/*!
 * \file CFVMDataSorter.cpp
 * \brief Datasorter class for FVM solvers.
 * \author T. Albring
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

#include "../../../include/output/filewriter/CFVMDataSorter.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include <numeric>

CFVMDataSorter::CFVMDataSorter(CConfig *config, CGeometry *geometry, const vector<string> &valFieldNames) :
  CParallelDataSorter(config, valFieldNames){

  nDim = geometry->GetnDim();

  std::vector<unsigned long> globalID;

  nGlobalPointBeforeSort = geometry->GetGlobal_nPointDomain();
  nLocalPointsBeforeSort  = geometry->GetnPointDomain();

  Local_Halo.resize(geometry->GetnPoint());

  for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){

    /*--- Store the global IDs ---*/

    globalID.push_back(geometry->nodes->GetGlobalIndex(iPoint));

    Local_Halo[iPoint] = !geometry->nodes->GetDomain(iPoint);
  }


  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. ---*/

  SetHaloPoints(geometry, config);

  /*--- Create the linear partitioner --- */

  linearPartitioner.Initialize(nGlobalPointBeforeSort, 0);

  /*--- Prepare the send buffers ---*/

  PrepareSendBuffers(globalID);

}

void CFVMDataSorter::SetHaloPoints(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, iVertex;
  unsigned short iMarker;
  int SendRecv, RecvFrom;
  bool notHalo;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- We need to keep one copy of overlapping halo cells. ---*/

        notHalo = ((SendRecv < 0) && (rank > RecvFrom));

        /*--- If we found either of these types of nodes, flag them to be kept. ---*/

        if (notHalo) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }
}

void CFVMDataSorter::SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) {

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/

  nElemPerType.fill(0);

  SortVolumetricConnectivity(config, geometry, TRIANGLE,      val_sort);
  SortVolumetricConnectivity(config, geometry, QUADRILATERAL, val_sort);
  SortVolumetricConnectivity(config, geometry, TETRAHEDRON,   val_sort);
  SortVolumetricConnectivity(config, geometry, HEXAHEDRON,    val_sort);
  SortVolumetricConnectivity(config, geometry, PRISM,         val_sort);
  SortVolumetricConnectivity(config, geometry, PYRAMID,       val_sort);

  SetTotalElements();

  connectivitySorted = true;

}

void CFVMDataSorter::SortVolumetricConnectivity(CConfig *config,
                                         CGeometry *geometry,
                                         unsigned short Elem_Type,
                                         bool val_sort) {

  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT = 0;
  unsigned long iPoint, jPoint;
  unsigned long nElem_Total = 0, Global_Index;

  int *Conn_Elem  = nullptr;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  switch (Elem_Type) {
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- We start with the connectivity distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first nPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  int *nElem_Flag = new int[size]();

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Cum[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Cum[size] = 0;

  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++ ) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {

        /*--- Get the index of the current point. ---*/

        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->nodes->GetGlobalIndex(iPoint);

        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/

        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->nodes->GetGlobalIndex(jPoint);
          if (newID < Global_Index) Global_Index = newID;
        }

        /*--- Search for the processor that owns this point. If we are
         sorting the elements, we use the linear partitioning to find
         the rank, otherwise, we simply have the current rank load its
         own elements into the connectivity data structure. ---*/

        if (val_sort) {
          iProcessor = linearPartitioner.GetRankContainingIndex(Global_Index);
        } else {
          iProcessor = rank;
        }


        /*--- If we have not visited this element yet, increment our
         number of elements that must be sent to a particular proc. ---*/

        if ((nElem_Flag[iProcessor] != ii)) {
          nElem_Flag[iProcessor] = ii;
          nElem_Send[iProcessor+1]++;
        }

      }
    }
  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Cum[1]), 1, MPI_INT, SU2_MPI::GetComm());

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Cum[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Cum[ii+1] += nElem_Cum[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/

  auto connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]]();

  /*--- Allocate arrays for storing halo flags. ---*/

  auto haloSend = new unsigned short[nElem_Send[size]]();

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  auto *index = new unsigned long[size]();
  for (int ii=0; ii < size; ii++) index[ii] = NODES_PER_ELEMENT*nElem_Send[ii];

  auto *haloIndex = new unsigned long[size]();
  for (int ii=0; ii < size; ii++) haloIndex[ii] = nElem_Send[ii];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {

        /*--- Get the index of the current point. ---*/

        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->nodes->GetGlobalIndex(iPoint);

        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/

        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->nodes->GetGlobalIndex(jPoint);
          if (newID < Global_Index) Global_Index = newID;
        }

        /*--- Search for the processor that owns this point. If we are
         sorting the elements, we use the linear partitioning to find
         the rank, otherwise, we simply have the current rank load its
         own elements into the connectivity data structure. ---*/

        if (val_sort) {
          iProcessor = linearPartitioner.GetRankContainingIndex(Global_Index);
        } else {
          iProcessor = rank;
        }

        /*--- Load connectivity into the buffer for sending ---*/

        if (nElem_Flag[iProcessor] != ii) {

          nElem_Flag[iProcessor] = ii;
          unsigned long nn = index[iProcessor];
          unsigned long mm = haloIndex[iProcessor];

          /*--- Load the connectivity values. ---*/

          for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
            iPoint = geometry->elem[ii]->GetNode(kk);
            connSend[nn] = geometry->nodes->GetGlobalIndex(iPoint); nn++;

            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. Note that just checking
             whether the point is a halo point is not enough, since we want to keep
             elements on one side of the send receive boundary. ---*/

            if (Local_Halo[iPoint]) haloSend[mm] = true;
          }

          /*--- Increment the index by the message length ---*/

          index[iProcessor]    += NODES_PER_ELEMENT;
          haloIndex[iProcessor]++;

        }
      }
    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] haloIndex;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  auto connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Cum[size]]();

  auto haloRecv = new unsigned short[nElem_Cum[size]]();

#ifdef HAVE_MPI

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  /*--- Launch the non-blocking recv's for the connectivity. ---*/

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Cum[ii+1] > nElem_Cum[ii])) {
      int ll     = NODES_PER_ELEMENT*nElem_Cum[ii];
      int kk     = nElem_Cum[ii+1] - nElem_Cum[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     SU2_MPI::GetComm(), &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = NODES_PER_ELEMENT*nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     SU2_MPI::GetComm(), &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the halo flags. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Cum[ii+1] > nElem_Cum[ii])) {
      int ll     = nElem_Cum[ii];
      int kk     = nElem_Cum[ii+1] - nElem_Cum[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(haloRecv[ll]), count, MPI_UNSIGNED_SHORT, source, tag,
                     SU2_MPI::GetComm(), &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the halo flags. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(haloSend[ll]), count, MPI_UNSIGNED_SHORT, dest, tag,
                     SU2_MPI::GetComm(), &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = NODES_PER_ELEMENT*nElem_Cum[rank];
  int ll = NODES_PER_ELEMENT*nElem_Send[rank];
  int kk = NODES_PER_ELEMENT*nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];

  mm = nElem_Cum[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) haloRecv[mm] = haloSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. Note that we add 1 here
   to the connectivity for vizualization packages. First, allocate
   appropriate amount of memory for this section. ---*/

  if (nElem_Cum[size] > 0) Conn_Elem = new int[NODES_PER_ELEMENT*nElem_Cum[size]]();
  int count = 0; nElem_Total = 0;
  for (int ii = 0; ii < nElem_Cum[size]; ii++) {
    if (!haloRecv[ii]) {
      nElem_Total++;
      for (int jj = 0; jj < NODES_PER_ELEMENT; jj++) {
        Conn_Elem[count] = (int)connRecv[ii*NODES_PER_ELEMENT+jj] + 1;
        count++;
      }
    }
  }

  nElemPerType[TypeMap.at(Elem_Type)] = nElem_Total;

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  switch (Elem_Type) {
    case TRIANGLE:
      delete [] Conn_Tria_Par;
      Conn_Tria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      delete [] Conn_Quad_Par;
      Conn_Quad_Par = Conn_Elem;
      break;
    case TETRAHEDRON:
      delete [] Conn_Tetr_Par;
      Conn_Tetr_Par = Conn_Elem;
      break;
    case HEXAHEDRON:
      delete [] Conn_Hexa_Par;
      Conn_Hexa_Par = Conn_Elem;
      break;
    case PRISM:
      delete [] Conn_Pris_Par;
      Conn_Pris_Par = Conn_Elem;
      break;
    case PYRAMID:
      delete [] Conn_Pyra_Par;
      Conn_Pyra_Par = Conn_Elem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }

  /*--- Free temporary memory from communications ---*/

  delete [] connSend;
  delete [] connRecv;
  delete [] haloSend;
  delete [] haloRecv;
  delete [] nElem_Flag;

}
