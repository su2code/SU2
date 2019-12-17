/*!
 * \file CParallelDataSorter.cpp
 * \brief Datasorter base class.
 * \author T. Albring
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

#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include <cassert>

CParallelDataSorter::CParallelDataSorter(CConfig *config, unsigned short nFields){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  GlobalField_Counter = nFields;

  nParallel_Hexa = 0;
  nParallel_Line = 0;
  nParallel_Quad = 0;
  nParallel_Tetr = 0;
  nParallel_Pris = 0;
  nParallel_Pyra = 0;
  nParallel_Tria = 0;

  Conn_Line_Par = NULL;
  Conn_Hexa_Par = NULL;
  Conn_Pris_Par = NULL;
  Conn_Quad_Par = NULL;
  Conn_Tetr_Par = NULL;
  Conn_Tria_Par = NULL;
  Conn_Pyra_Par = NULL;

  nPoint_Send  = NULL;
  nPoint_Recv  = NULL;
  Index        = NULL;
  connSend     = NULL;
  dataBuffer   = NULL;
  passiveDoubleBuffer = NULL;
  doubleBuffer = NULL;
  idSend       = NULL;
  nSends = 0;
  nRecvs = 0;

  nLocalPoint_Sort  = 0;
  nGlobalPoint_Sort = 0;

  nPoint_Send = new int[size+1]();
  nPoint_Recv = new int[size+1]();

  linearPartitioner = NULL;

}

CParallelDataSorter::~CParallelDataSorter(){

  if (nPoint_Send != NULL) delete [] nPoint_Send;
  if (nPoint_Recv != NULL) delete [] nPoint_Recv;

  /*--- Deallocate memory for connectivity data on each processor. ---*/

  if (nParallel_Line > 0 && Conn_Line_Par != NULL) delete [] Conn_Line_Par;
  if (nParallel_Tria > 0 && Conn_Tria_Par != NULL) delete [] Conn_Tria_Par;
  if (nParallel_Quad > 0 && Conn_Quad_Par != NULL) delete [] Conn_Quad_Par;
  if (nParallel_Tetr > 0 && Conn_Tetr_Par != NULL) delete [] Conn_Tetr_Par;
  if (nParallel_Hexa > 0 && Conn_Hexa_Par != NULL) delete [] Conn_Hexa_Par;
  if (nParallel_Pris > 0 && Conn_Pris_Par != NULL) delete [] Conn_Pris_Par;
  if (nParallel_Pyra > 0 && Conn_Pyra_Par != NULL) delete [] Conn_Pyra_Par;

  if (connSend != NULL) delete [] connSend;

  if (dataBuffer != NULL) delete [] dataBuffer;
}


unsigned long CParallelDataSorter::GetnElem(GEO_TYPE type){

  switch (type) {
    case LINE:
      return nParallel_Line;
      break;
    case TRIANGLE:
      return nParallel_Tria;
      break;
    case QUADRILATERAL:
      return nParallel_Quad;
      break;
    case TETRAHEDRON:
      return nParallel_Tetr;
      break;
    case HEXAHEDRON:
      return nParallel_Hexa;
      break;
    case PRISM:
      return nParallel_Pris;
      break;
    case PYRAMID:
      return nParallel_Pyra;
      break;
    default:
      break;
  }

  SU2_MPI::Error("GEO_TYPE not found", CURRENT_FUNCTION);

  return 0;
}

void CParallelDataSorter::SortOutputData() {

  int VARS_PER_POINT = GlobalField_Counter;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/


  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]]();

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the global IDs. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = VARS_PER_POINT*nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = VARS_PER_POINT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(doubleBuffer[ll]), count, MPI_DOUBLE, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = VARS_PER_POINT*nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = VARS_PER_POINT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_DOUBLE, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = VARS_PER_POINT*nPoint_Recv[rank];
  int ll = VARS_PER_POINT*nPoint_Send[rank];
  int kk = VARS_PER_POINT*nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) doubleBuffer[mm] = connSend[nn];

  mm = nPoint_Recv[rank];
  ll = nPoint_Send[rank];
  kk = nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

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

  /*--- Note, passiveDoubleBuffer and doubleBuffer point to the same address.
   * This is the reason why we have to do the following copy/reordering in two steps. ---*/

  /*--- Step 1: Extract the underlying double value --- */

  if (!std::is_same<su2double, passivedouble>::value){
    for (int jj = 0; jj < VARS_PER_POINT*nPoint_Recv[size]; jj++){
      const passivedouble tmpVal = SU2_TYPE::GetValue(doubleBuffer[jj]);
      passiveDoubleBuffer[jj] = tmpVal;
      /*--- For some AD datatypes a call of the destructor is
       *  necessary to properly delete the AD type ---*/
      doubleBuffer[jj].~su2double();
    }
  }

  /*--- Step 2: Reorder the data in the buffer --- */

  passivedouble *tmpBuffer = new passivedouble[nPoint_Recv[size]];
  for (int jj = 0; jj < VARS_PER_POINT; jj++){
    for (int ii = 0; ii < nPoint_Recv[size]; ii++){
      tmpBuffer[idRecv[ii]] = passiveDoubleBuffer[ii*VARS_PER_POINT+jj];
    }
    for (int ii = 0; ii < nPoint_Recv[size]; ii++){
      passiveDoubleBuffer[ii*VARS_PER_POINT+jj] = tmpBuffer[ii];
    }
  }

  delete [] tmpBuffer;

  /*--- Store the total number of local points my rank has for
   the current section after completing the communications. ---*/

  nParallel_Poin = nPoint_Recv[size];

  /*--- Reduce the total number of points we will write in the output files. ---*/

#ifndef HAVE_MPI
  nGlobal_Poin_Par = nParallel_Poin;
#else
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Free temporary memory from communications ---*/

  delete [] idRecv;
}

void CParallelDataSorter::PrepareSendBuffers(std::vector<unsigned long>& globalID){

  unsigned long iPoint;
  unsigned short iProcessor;

  int VARS_PER_POINT = GlobalField_Counter;

  /*--- We start with the grid nodes distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many nodes we must send to each other rank in order to
   have all nodes sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first ~ nGlobalPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++ ) {

    iProcessor = linearPartitioner->GetRankContainingIndex(globalID[iPoint]);

    /*--- If we have not visited this node yet, increment our
       number of elements that must be sent to a particular proc. ---*/

    nPoint_Send[iProcessor+1]++;
  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif

  /*--- Prepare to send coordinates. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  nSends = 0; nRecvs = 0;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nPoint_Recv[ii+1] > 0)) nRecvs++;

    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/

  connSend = NULL;
  connSend = new su2double[VARS_PER_POINT*nPoint_Send[size]]();

  /*--- Allocate the data buffer to hold the sorted data. We have to make it large enough
   * to hold passivedoubles and su2doubles ---*/
  unsigned short maxSize = max(sizeof(passivedouble), sizeof(su2double));
  dataBuffer = new char[VARS_PER_POINT*nPoint_Recv[size]*maxSize];

  /*--- doubleBuffer and passiveDouble buffer use the same memory allocated above using the dataBuffer. ---*/

  doubleBuffer = reinterpret_cast<su2double*>(dataBuffer);
  passiveDoubleBuffer = reinterpret_cast<passivedouble*>(dataBuffer);

  /*--- Allocate arrays for sending the global ID. ---*/

  idSend = new unsigned long[nPoint_Send[size]]();

  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size]();
  for (int ii=0; ii < size; ii++) index[ii] = VARS_PER_POINT*nPoint_Send[ii];

  unsigned long *idIndex = new unsigned long[size]();
  for (int ii=0; ii < size; ii++) idIndex[ii] = nPoint_Send[ii];

  Index = new unsigned long[nLocalPoint_Sort]();

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++) {

    iProcessor = linearPartitioner->GetRankContainingIndex(globalID[iPoint]);

    /*--- Load the global ID (minus offset) for sorting the
         points once they all reach the correct processor. ---*/

    unsigned long nn = idIndex[iProcessor];
    idSend[nn] = globalID[iPoint] - linearPartitioner->GetFirstIndexOnRank(iProcessor);

    /*--- Store the index this point has in the send buffer ---*/

    Index[iPoint] = index[iProcessor];

    /*--- Increment the index by the message length ---*/

    index[iProcessor]  += VARS_PER_POINT;
    idIndex[iProcessor]++;


  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] idIndex;
}

unsigned long CParallelDataSorter::GetElem_Connectivity(GEO_TYPE type, unsigned long iElem, unsigned long iNode) {

  switch (type) {
    case LINE:
      return Conn_Line_Par[iElem*N_POINTS_LINE + iNode];
      break;
    case TRIANGLE:
      return Conn_Tria_Par[iElem*N_POINTS_TRIANGLE + iNode];
      break;
    case QUADRILATERAL:
      return Conn_Quad_Par[iElem*N_POINTS_QUADRILATERAL + iNode];
      break;
    case TETRAHEDRON:
      return Conn_Tetr_Par[iElem*N_POINTS_TETRAHEDRON + iNode];
      break;
    case HEXAHEDRON:
      return Conn_Hexa_Par[iElem*N_POINTS_HEXAHEDRON + iNode];
      break;
    case PRISM:
      return Conn_Pris_Par[iElem*N_POINTS_PRISM + iNode];
      break;
    case PYRAMID:
      return Conn_Pyra_Par[iElem*N_POINTS_PYRAMID + iNode];
      break;
    default:
      break;
  }

  SU2_MPI::Error("GEO_TYPE not found", CURRENT_FUNCTION);

  return 0;
}

