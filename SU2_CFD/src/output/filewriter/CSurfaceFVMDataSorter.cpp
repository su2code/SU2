/*!
 * \file CSurfaceFVMDataSorter.cpp
 * \brief Datasorter for FVM surfaces.
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

#include "../../../include/output/filewriter/CSurfaceFVMDataSorter.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"


CSurfaceFVMDataSorter::CSurfaceFVMDataSorter(CConfig *config, CGeometry *geometry, unsigned short nFields, CFVMDataSorter* volume_sorter) : CParallelDataSorter(config, nFields){

  this->volume_sorter = volume_sorter;

  connectivity_sorted = false;

  nGlobalPoint_Sort = geometry->GetGlobal_nPointDomain();
  nLocalPoint_Sort  = geometry->GetnPointDomain();

  /*--- Create the linear partitioner --- */

  linearPartitioner = new CLinearPartitioner(nGlobalPoint_Sort, 0);

}

CSurfaceFVMDataSorter::~CSurfaceFVMDataSorter(){

  if (linearPartitioner != NULL) delete linearPartitioner;
  delete [] passiveDoubleBuffer;

}

void CSurfaceFVMDataSorter::SortOutputData() {

  unsigned long iProcessor;
  unsigned long iPoint, iElem;
  unsigned long Global_Index;

  int VARS_PER_POINT = GlobalField_Counter;
  int *Local_Halo = NULL;
  int iNode, count;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif


  /*--- Prepare to check and communicate the nodes that each proc has
   locally from the surface connectivity. ---*/

  int *nElem_Send = new int[size+1](); nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1](); nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size]();

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  /*--- Loop through our local line elements and check where each
   of the grid nodes resides based on global index. ---*/

  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_Line_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Reset out flags and then loop through our local triangle surface
   elements performing the same check for where each grid node resides. ---*/

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_Tria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_Tria_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Reset out flags and then loop through our local quad surface
   elements performing the same check for where each grid node resides. ---*/

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_Quad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_Quad_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many nodes it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  /*--- Allocate arrays for sending the global ID. ---*/

  unsigned long *idSend = new unsigned long[nElem_Send[size]]();

  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *idIndex = new unsigned long[size]();
  for (int ii=0; ii < size; ii++) idIndex[ii] = nElem_Send[ii];

  /*--- Now loop back through the local connectivities for the surface
   elements and load up the global IDs for sending to their home proc. ---*/

  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_Line_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- Load global ID into the buffer for sending ---*/

      if (nElem_Flag[iProcessor] != iNode) {

        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];

        /*--- Load the connectivity values. ---*/

        idSend[nn] = Global_Index; nn++;

        /*--- Increment the index by the message length ---*/

        idIndex[iProcessor]++;

      }

    }
  }

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_Tria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_Tria_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- Load global ID into the buffer for sending ---*/

      if (nElem_Flag[iProcessor] != iNode) {

        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];

        /*--- Load the connectivity values. ---*/

        idSend[nn] = Global_Index; nn++;

        /*--- Increment the index by the message length ---*/

        idIndex[iProcessor]++;

      }

    }
  }

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_Quad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_Quad_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- Load global ID into the buffer for sending ---*/

      if (nElem_Flag[iProcessor] != iNode) {

        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];

        /*--- Load the connectivity values. ---*/

        idSend[nn] = Global_Index; nn++;

        /*--- Increment the index by the message length ---*/

        idIndex[iProcessor]++;

      }

    }
  }

  /*--- Allocate the memory that we need for receiving the global IDs
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *idRecv = new unsigned long[nElem_Recv[size]]();

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[nSends];
  recv_req = new SU2_MPI::Request[nRecvs];

  /*--- Launch the non-blocking recv's for the global IDs. ---*/

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = nElem_Recv[rank];
  int ll = nElem_Send[rank];
  int kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Each proc now knows which is its local grid nodes from     ---*/
  /*---         the entire volume solution are part of the surface. We     ---*/
  /*---         now apply a mask to extract just those points on the       ---*/
  /*---         surface. We also need to perform a renumbering so that     ---*/
  /*---         the surface data (nodes and connectivity) have their       ---*/
  /*---         own global numbering. This is important for writing        ---*/
  /*---         output files in a later routine.                           ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Create a local data structure that acts as a mask to extract the
   set of points within the local set that are on the surface. ---*/

  int *surfPoint = new int[volume_sorter->GetnPoints()];
  for (iPoint = 0; iPoint < volume_sorter->GetnPoints(); iPoint++) surfPoint[iPoint] = -1;

  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    surfPoint[(int)idRecv[ii] - volume_sorter->GetNodeBegin(rank)] = (int)idRecv[ii];
  }

  /*--- First, add up the number of surface points I have on my rank. ---*/

  nParallel_Poin = 0;
  Renumber2Global.clear();

  for (iPoint = 0; iPoint < volume_sorter->GetnPoints(); iPoint++) {
    if (surfPoint[iPoint] != -1) {

      /*--- Save the global index values for CSV output. ---*/

      Renumber2Global[nParallel_Poin] = surfPoint[iPoint];

      /*--- Increment total number of surface points found locally. ---*/

      nParallel_Poin++;

    }
  }

  /*--- Communicate this number of local surface points to all other
   processors so that it can be used to create offsets for the new
   global numbering for the surface points. ---*/

  int *nPoint_Send = new int[size+1](); nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1](); nPoint_Recv[0] = 0;

  for (int ii=1; ii < size+1; ii++) nPoint_Send[ii]= (int)nParallel_Poin;

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif

  /*--- Go to cumulative storage format to compute the offsets. ---*/

  for (int ii = 0; ii < size; ii++) {
    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }

  /*--- Now that we know the number of local surface points that we have,
   we can allocate the new data structure to hold these points alone. Here,
   we also copy the data for those points from our volume data structure. ---*/

  if (passiveDoubleBuffer == nullptr){
    passiveDoubleBuffer = new passivedouble[nParallel_Poin*VARS_PER_POINT];
  }
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    count = 0;
    for (int ii = 0; ii < (int)volume_sorter->GetnPoints(); ii++) {
      if (surfPoint[ii] !=-1) {
        passiveDoubleBuffer[count*VARS_PER_POINT + jj] = volume_sorter->GetData(jj,ii);
        count++;
      }
    }
  }
  /*--- Reduce the total number of surf points we have. This will be
   needed for writing the surface solution files later. ---*/

#ifndef HAVE_MPI
  nGlobal_Poin_Par = nParallel_Poin;
#else
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Now that we know every proc's global offset for the number of
   surface points, we can create the new global numbering. Here, we
   create a new mapping using two arrays, which will need to be
   communicated. We use our mask again here.  ---*/

  unsigned long *globalP = new unsigned long[nParallel_Poin]();
  unsigned long *renumbP = new unsigned long[nParallel_Poin]();

  count = 0;
  for (iPoint = 0; iPoint < volume_sorter->GetnPoints(); iPoint++) {
    if (surfPoint[iPoint] != -1) {
      globalP[count] = surfPoint[iPoint];
      renumbP[count] = count + nPoint_Recv[rank];
      count++;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Communicate the arrays with the new global surface point   ---*/
  /*---         numbering to the procs that hold the connectivity for      ---*/
  /*---         each element. This will be done in two phases. First,      ---*/
  /*---         we send the arrays around to the other procs based on      ---*/
  /*---         the linear partitioning for the elems. This gets us        ---*/
  /*---         most of the way there, however, due to the type of         ---*/
  /*---         linear partitioning for the elements, there may exist      ---*/
  /*---         elements that have nodes outside of the linear part.       ---*/
  /*---         bounds. This is because the elems are distributed based    ---*/
  /*---         on the node with the smallest global ID.                   ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Reset our flags and counters ---*/

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  /*--- Loop through my local surface nodes, find which proc the global
   value lives on, then communicate the global ID and remumbered value. ---*/

  for (int ii = 0; ii < (int)nParallel_Poin; ii++) {

    Global_Index = globalP[ii];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if ((nElem_Flag[iProcessor] != ii)) {
      nElem_Flag[iProcessor] = ii;
      nElem_Send[iProcessor+1]++;
    }

  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  nSends = 0; nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  /*--- Allocate memory to hold the globals that we are
   sending. ---*/

  unsigned long *globalSend = NULL;
  globalSend = new unsigned long[nElem_Send[size]]();

  /*--- Allocate memory to hold the renumbering that we are
   sending. ---*/

  unsigned long *renumbSend = NULL;
  renumbSend = new unsigned long[nElem_Send[size]]();

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size]();
  for (int ii=0; ii < size; ii++) index[ii] = nElem_Send[ii];

  /*--- Loop back through and load up the buffers for the global IDs
   and their new renumbering values. ---*/

  for (int ii = 0; ii < (int)nParallel_Poin; ii++) {

    Global_Index = globalP[ii];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

    if (nElem_Flag[iProcessor] != ii) {

      nElem_Flag[iProcessor] = ii;
      unsigned long nn = index[iProcessor];

      globalSend[nn] = Global_Index;
      renumbSend[nn] = renumbP[ii];

      /*--- Increment the index by the message length ---*/

      index[iProcessor]++;

    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;

  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *globalRecv = NULL;
  globalRecv = new unsigned long[nElem_Recv[size]]();

  unsigned long *renumbRecv = NULL;
  renumbRecv = new unsigned long[nElem_Recv[size]]();

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  /*--- Launch the non-blocking recv's for the global ID. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(globalRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global ID. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(globalSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking recv's for the renumbered ID. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(renumbRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the renumbered ID. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(renumbSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }

#endif

  /*--- Load our own procs data into the buffers directly. ---*/

  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) globalRecv[mm] = globalSend[nn];

  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) renumbRecv[mm] = renumbSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*-- Now update my local connectivitiy for the surface with the new
   numbering. Create a new mapping for global -> renumber for nodes. Note
   the adding of 1 back in here for the eventual viz. purposes. ---*/

  map<unsigned long,unsigned long> Global2Renumber;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    Global2Renumber[globalRecv[ii]] = renumbRecv[ii] + 1;
  }


  /*--- The final step is one last pass over all elements to check
   for points outside of the linear partitions of the elements. Again,
   note that elems were distributed based on their smallest global ID,
   so some nodes of the elem may have global IDs lying outside of the
   linear partitioning. We need to recover the mapping for these
   outliers. We loop over all local surface elements to find these. ---*/

  vector<unsigned long>::iterator it;
  vector<unsigned long> outliers;

  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {

      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_Line_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- Store the global ID if it is outside our own linear partition. ---*/

      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }

    }
  }

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_Tria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {

      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_Tria_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- Store the global ID if it is outside our own linear partition. ---*/

      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }

    }
  }

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_Quad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {

      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_Quad_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

      /*--- Store the global ID if it is outside our own linear partition. ---*/

      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }

    }
  }

  /*--- Create a unique list of global IDs that fall outside of our procs
   linear partition. ---*/

  sort(outliers.begin(), outliers.end());
  it = unique(outliers.begin(), outliers.end());
  outliers.resize(it - outliers.begin());

  /*--- Now loop over the outliers and communicate to those procs that
   hold the new numbering for our outlier points. We need to ask for the
   new numbering from these procs. ---*/

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (int ii = 0; ii < (int)outliers.size(); ii++) {

    Global_Index = outliers[ii];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if ((nElem_Flag[iProcessor] != ii)) {
      nElem_Flag[iProcessor] = ii;
      nElem_Send[iProcessor+1]++;
    }

  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  nSends = 0; nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  delete [] idSend;
  idSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    idSend[ii] = 0;

  /*--- Reset our index variable for reuse. ---*/

  for (int ii=0; ii < size; ii++) idIndex[ii] = nElem_Send[ii];

  /*--- Loop over the outliers again and load up the global IDs. ---*/

  for (int ii = 0; ii < (int)outliers.size(); ii++) {

    Global_Index = outliers[ii];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if ((nElem_Flag[iProcessor] != ii)) {

      nElem_Flag[iProcessor] = ii;
      unsigned long nn = idIndex[iProcessor];

      /*--- Load the global ID values. ---*/

      idSend[nn] = Global_Index; nn++;

      /*--- Increment the index by the message length ---*/

      idIndex[iProcessor]++;

    }
  }

  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  delete [] idRecv;
  idRecv = new unsigned long[nElem_Recv[size]]();

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[nSends];
  recv_req = new SU2_MPI::Request[nRecvs];

  /*--- Launch the non-blocking recv's for the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- The procs holding the outlier grid nodes now have the global IDs
   that they need to have their renumbering shared. ---*/

  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
      if (idRecv[ii] == globalP[iPoint]) {
        idRecv[ii] = renumbP[iPoint];
      }
    }
  }

  /*--- Now simply reverse the last communication to give the renumbered IDs
   back to the owner of the outlier points. Note everything is flipped. ---*/

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[nRecvs];
  recv_req = new SU2_MPI::Request[nSends];

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking recv's for the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  mm = nElem_Send[rank];
  ll = nElem_Recv[rank];
  kk = nElem_Recv[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idSend[mm] = idRecv[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Add the renumbering for the outliers to the map from before carrying
   the global -> renumber transformation. Note that by construction,
   nElem_Send[ii] == outliers.size(). We also add in the 1 for viz. here. ---*/

  for (int ii = 0; ii < nElem_Send[size]; ii++) {
    Global2Renumber[outliers[ii]] = idSend[ii] + 1;
  }

  /*--- We can now overwrite the local connectivity for our surface elems
   using our completed map with the new global renumbering. Whew!! Note
   the -1 when accessing the conn from the map. ---*/

  for (iElem = 0; iElem < nParallel_Line; iElem++) {
    iNode = (int)iElem*N_POINTS_LINE;
    Conn_Line_Par[iNode+0] = (int)Global2Renumber[Conn_Line_Par[iNode+0]-1];
    Conn_Line_Par[iNode+1] = (int)Global2Renumber[Conn_Line_Par[iNode+1]-1];
  }

  for (iElem = 0; iElem < nParallel_Tria; iElem++) {
    iNode = (int)iElem*N_POINTS_TRIANGLE;
    Conn_Tria_Par[iNode+0] = (int)Global2Renumber[Conn_Tria_Par[iNode+0]-1];
    Conn_Tria_Par[iNode+1] = (int)Global2Renumber[Conn_Tria_Par[iNode+1]-1];
    Conn_Tria_Par[iNode+2] = (int)Global2Renumber[Conn_Tria_Par[iNode+2]-1];
  }

  for (iElem = 0; iElem < nParallel_Quad; iElem++) {
    iNode = (int)iElem*N_POINTS_QUADRILATERAL;
    Conn_Quad_Par[iNode+0] = (int)Global2Renumber[Conn_Quad_Par[iNode+0]-1];
    Conn_Quad_Par[iNode+1] = (int)Global2Renumber[Conn_Quad_Par[iNode+1]-1];
    Conn_Quad_Par[iNode+2] = (int)Global2Renumber[Conn_Quad_Par[iNode+2]-1];
    Conn_Quad_Par[iNode+3] = (int)Global2Renumber[Conn_Quad_Par[iNode+3]-1];
  }

  /*--- Free temporary memory ---*/

  delete [] idIndex;
  delete [] surfPoint;
  delete [] globalP;
  delete [] renumbP;

  delete [] idSend;
  delete [] idRecv;
  delete [] globalSend;
  delete [] globalRecv;
  delete [] renumbSend;
  delete [] renumbRecv;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Local_Halo;
  delete [] nPoint_Send;
  delete [] nPoint_Recv;

}

void CSurfaceFVMDataSorter::SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) {

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/

  /*--- Sort volumetric grid connectivity. ---*/

  SortSurfaceConnectivity(config, geometry, LINE         );
  SortSurfaceConnectivity(config, geometry, TRIANGLE     );
  SortSurfaceConnectivity(config, geometry, QUADRILATERAL);


  unsigned long nTotal_Surf_Elem = nParallel_Line + nParallel_Tria + nParallel_Quad;
#ifndef HAVE_MPI
  nGlobal_Elem_Par   = nTotal_Surf_Elem;
#else
  SU2_MPI::Allreduce(&nTotal_Surf_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  connectivity_sorted = true;

}

void CSurfaceFVMDataSorter::SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT;
  unsigned long iPoint, jPoint;
  unsigned long nElem_Total = 0, Global_Index;

  unsigned long iMarker;

  int *Conn_Elem  = NULL;

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
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      NODES_PER_ELEMENT = 0;
      break;
  }

  /*--- We start with the connectivity distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first nPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {

      for (int ii = 0; ii < (int)geometry->GetnElem_Bound(iMarker); ii++) {

        if (geometry->bound[iMarker][ii]->GetVTK_Type() == Elem_Type) {
          for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {

            /*--- Get the index of the current point. ---*/

            iPoint = geometry->bound[iMarker][ii]->GetNode(jj);
            Global_Index = geometry->node[iPoint]->GetGlobalIndex();

            /*--- Search for the lowest global index in this element. We
             send the element to the processor owning the range that includes
             the lowest global index value. ---*/

            for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
              jPoint = geometry->bound[iMarker][ii]->GetNode(kk);
              unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
              if (newID < Global_Index) Global_Index = newID;
            }

            /*--- Search for the processor that owns this point ---*/

            iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

            /*--- If we have not visited this element yet, increment our
             number of elements that must be sent to a particular proc. ---*/

            if ((nElem_Flag[iProcessor] != ii)) {
              nElem_Flag[iProcessor] = ii;
              nElem_Send[iProcessor+1]++;
            }

          }
        }
      }
    }
  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/

  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Send[size]; ii++)
    connSend[ii] = 0;

  /*--- Allocate arrays for storing halo flags. ---*/

  unsigned short *haloSend = new unsigned short[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    haloSend[ii] = false;

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = NODES_PER_ELEMENT*nElem_Send[ii];

  unsigned long *haloIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) haloIndex[ii] = nElem_Send[ii];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {

      for (int ii = 0; ii < (int)geometry->GetnElem_Bound(iMarker); ii++) {

        if (geometry->bound[iMarker][ii]->GetVTK_Type() == Elem_Type) {
          for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {

            /*--- Get the index of the current point. ---*/

            iPoint = geometry->bound[iMarker][ii]->GetNode(jj);
            Global_Index = geometry->node[iPoint]->GetGlobalIndex();

            /*--- Search for the lowest global index in this element. We
             send the element to the processor owning the range that includes
             the lowest global index value. ---*/

            for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
              jPoint = geometry->bound[iMarker][ii]->GetNode(kk);
              unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
              if (newID < Global_Index) Global_Index = newID;
            }

            /*--- Search for the processor that owns this point ---*/

            iProcessor = linearPartitioner->GetRankContainingIndex(Global_Index);

            /*--- Load connectivity into the buffer for sending ---*/

            if (nElem_Flag[iProcessor] != ii) {

              nElem_Flag[iProcessor] = ii;
              unsigned long nn = index[iProcessor];
              unsigned long mm = haloIndex[iProcessor];

              /*--- Load the connectivity values. ---*/

              for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
                iPoint = geometry->bound[iMarker][ii]->GetNode(kk);
                connSend[nn] = geometry->node[iPoint]->GetGlobalIndex(); nn++;

                /*--- Check if this is a halo node. If so, flag this element
                 as a halo cell. We will use this later to sort and remove
                 any duplicates from the connectivity list. ---*/

                if (volume_sorter->GetHalo(iPoint)) haloSend[mm] = true;

              }

              /*--- Increment the index by the message length ---*/

              index[iProcessor]    += NODES_PER_ELEMENT;
              haloIndex[iProcessor]++;

            }
          }
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

  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Recv[size]; ii++)
    connRecv[ii] = 0;

  unsigned short *haloRecv = new unsigned short[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    haloRecv[ii] = false;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  /*--- Launch the non-blocking recv's for the connectivity. ---*/

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = NODES_PER_ELEMENT*nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
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
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the halo flags. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(haloRecv[ll]), count, MPI_UNSIGNED_SHORT, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
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
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = NODES_PER_ELEMENT*nElem_Recv[rank];
  int ll = NODES_PER_ELEMENT*nElem_Send[rank];
  int kk = NODES_PER_ELEMENT*nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];

  mm = nElem_Recv[rank];
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

  if (nElem_Recv[size] > 0) Conn_Elem = new int[NODES_PER_ELEMENT*nElem_Recv[size]];
  int count = 0; nElem_Total = 0;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    if (!haloRecv[ii]) {
      nElem_Total++;
      for (int jj = 0; jj < NODES_PER_ELEMENT; jj++) {
        Conn_Elem[count] = (int)connRecv[ii*NODES_PER_ELEMENT+jj] + 1;
        count++;
      }
    }
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  switch (Elem_Type) {
    case LINE:
      nParallel_Line = nElem_Total;
      if (Conn_Line_Par != NULL) delete [] Conn_Line_Par;
      if (nParallel_Line > 0) Conn_Line_Par = Conn_Elem;
      break;
    case TRIANGLE:
      nParallel_Tria = nElem_Total;
      if (Conn_Tria_Par != NULL) delete [] Conn_Tria_Par;
      if (nParallel_Tria > 0) Conn_Tria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      nParallel_Quad = nElem_Total;
      if (Conn_Quad_Par != NULL) delete [] Conn_Quad_Par;
      if (nParallel_Quad > 0) Conn_Quad_Par = Conn_Elem;
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
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;

}
