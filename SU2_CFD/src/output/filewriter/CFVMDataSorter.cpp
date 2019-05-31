#include "../../../include/output/filewriter/CFVMDataSorter.hpp"


CFVMDataSorter::CFVMDataSorter(CConfig *config, CGeometry *geometry, unsigned short nFields, std::vector<std::vector<su2double> >& Local_Data) : CParallelDataSorter(config, nFields){
 
  this->Local_Data = &Local_Data;
  
  unsigned long iPoint;
  
  /*--- Reset point sorting counters ---*/

  nGlobalPoint_Sort = 0;
  nLocalPoint_Sort  = 0;

  /*--- Search all send/recv boundaries on this partition for any periodic
     nodes that were part of the original domain. We want to recover these
     for visualization purposes. ---*/
  
  unsigned long iVertex;
  bool isPeriodic;
  
  Local_Halo_Sort = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo_Sort[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      
      /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo_Sort[iPoint] = false;
      }
    }
  }
  
  /*--- Sum total number of nodes that belong to the domain ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo_Sort[iPoint] == false)
      nLocalPoint_Sort++;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nGlobalPoint_Sort = nLocalPoint_Sort;
#endif
  
  
  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  beg_node = new unsigned long[size];
  end_node = new unsigned long[size];

  nPoint_Lin = new unsigned long[size];
  nPoint_Cum = new unsigned long[size+1];

  unsigned long total_points = 0;
  for (int ii = 0; ii < size; ii++) {
    nPoint_Lin[ii] = nGlobalPoint_Sort/size;
    total_points  += nPoint_Lin[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long remainder = nGlobalPoint_Sort - total_points;
  for (unsigned long ii = 0; ii < remainder; ii++) {
    nPoint_Lin[ii]++;
  }

  /*--- Store the local number of nodes on each proc in the linear
   partitioning, the beginning/end index, and the linear partitioning
   within an array in cumulative storage format. ---*/

  beg_node[0] = 0;
  end_node[0] = beg_node[0] + nPoint_Lin[0];
  nPoint_Cum[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    beg_node[ii]   = end_node[ii-1];
    end_node[ii]   = beg_node[ii] + nPoint_Lin[ii];
    nPoint_Cum[ii] = nPoint_Cum[ii-1] + nPoint_Lin[ii-1];
  }
  nPoint_Cum[size] = nGlobalPoint_Sort;
  
}

CFVMDataSorter::~CFVMDataSorter(){
  
  delete [] beg_node;
  delete [] end_node;
  delete [] nPoint_Cum;
  delete [] nPoint_Lin;
  
  delete [] Local_Halo_Sort;
  
}



void CFVMDataSorter::SortOutputData(CConfig *config, CGeometry *geometry) {
  
  unsigned long iProcessor;
  unsigned long iPoint, Global_Index, nTotalPoint;
  
  int VARS_PER_POINT = GlobalField_Counter;
  
#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint_Sort, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint_Sort;
#endif
  
  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;
  
  /*--- We start with the grid nodes distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many nodes we must send to each other rank in order to
   have all nodes sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first ~ nGlobalPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/
  
  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nPoint_Send[ii] = 0;
    nPoint_Recv[ii] = 0;
    nPoint_Flag[ii]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++ ) {
    
    /*--- We only write interior points and recovered periodic points. ---*/
    
    if (!Local_Halo_Sort[iPoint]) {
      
      /*--- Get the global index of the current point. ---*/
      
      Global_Index = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear[iProcessor])
        while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this node yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        nPoint_Flag[iProcessor] = (int)iPoint;
        nPoint_Send[iProcessor+1]++;
      }
      
    }
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
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nPoint_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nPoint_Recv[ii+1] > 0)) nRecvs++;
    
    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }
  
  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/
  
  su2double *connSend = NULL;
  connSend = new su2double[VARS_PER_POINT*nPoint_Send[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Send[size]; ii++)
    connSend[ii] = 0;
  
  /*--- Allocate arrays for sending the global ID. ---*/
  
  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (int ii = 0; ii < nPoint_Send[size]; ii++)
    idSend[ii] = 0;
  
  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = VARS_PER_POINT*nPoint_Send[ii];
  
  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nPoint_Send[ii];
  
  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- We only write interior points and recovered periodic points. ---*/
    
    if (!Local_Halo_Sort[iPoint]) {
      
      /*--- Get the index of the current point. ---*/
      
      Global_Index = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Search for the processor that owns this point. ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear[iProcessor])
        while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
      
      /*--- Load node coordinates into the buffer for sending. ---*/
      
      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        
        nPoint_Flag[iProcessor] = (int)iPoint;
        unsigned long nn = index[iProcessor];
        
        /*--- Load the data values. ---*/
        
        for (unsigned short kk = 0; kk < VARS_PER_POINT; kk++) {
          connSend[nn] = (*Local_Data)[iPoint][kk]; nn++;
        }
        
        /*--- Load the global ID (minus offset) for sorting the
         points once they all reach the correct processor. ---*/
        
        nn = idIndex[iProcessor];
        idSend[nn] = Global_Index - starting_node[iProcessor];
        
        /*--- Increment the index by the message length ---*/
        
        index[iProcessor]  += VARS_PER_POINT;
        idIndex[iProcessor]++;
        
      }
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  delete [] idIndex;
  
  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  su2double *connRecv = NULL;
  connRecv = new su2double[VARS_PER_POINT*nPoint_Recv[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Recv[size]; ii++)
    connRecv[ii] = 0;
  
  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (int ii = 0; ii < nPoint_Recv[size]; ii++)
    idRecv[ii] = 0;
  
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
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_DOUBLE, source, tag,
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
  
  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];
  
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
  
  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate the
   appropriate amount of memory for this section. ---*/
  
  Parallel_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Data[jj] = new su2double[nPoint_Recv[size]];
    for (int ii = 0; ii < nPoint_Recv[size]; ii++) {
      Parallel_Data[jj][idRecv[ii]] = connRecv[ii*VARS_PER_POINT+jj];
    }
  }
  
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
  
  delete [] connSend;
  delete [] connRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;

  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;
  
}


void CFVMDataSorter::SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) {

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/
  
  /*--- Sort volumetric grid connectivity. ---*/

  if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
    cout <<"Sorting volumetric grid connectivity." << endl;
  
  SortVolumetricConnectivity(config, geometry, TRIANGLE,      val_sort);
  SortVolumetricConnectivity(config, geometry, QUADRILATERAL, val_sort);
  SortVolumetricConnectivity(config, geometry, TETRAHEDRON,   val_sort);
  SortVolumetricConnectivity(config, geometry, HEXAHEDRON,    val_sort);
  SortVolumetricConnectivity(config, geometry, PRISM,         val_sort);
  SortVolumetricConnectivity(config, geometry, PYRAMID,       val_sort);  
  
  
  /*--- Reduce the total number of cells we will be writing in the output files. ---*/
  
  unsigned long nTotal_Elem = nParallel_Tria + nParallel_Quad + nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;
#ifndef HAVE_MPI
  nGlobal_Elem_Par = nTotal_Elem;
#else
  SU2_MPI::Allreduce(&nTotal_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
}

void CFVMDataSorter::SortVolumetricConnectivity(CConfig *config,
                                         CGeometry *geometry,
                                         unsigned short Elem_Type,
                                         bool val_sort) {
  
  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT = 0;
  unsigned long iPoint, jPoint, kPoint, nLocalPoint, nTotalPoint;
  unsigned long nElem_Total = 0, Global_Index;
  
  unsigned long iVertex, iMarker;
  int SendRecv, RecvFrom;
  
  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;
  
  int *Local_Halo = NULL;
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
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  
  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found. ---*/
  
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/
        
        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }
        
        notPeriodic = (isPeriodic && (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          
          kPoint = (iProcessor+1)*maxAddedPeriodic;
          
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
        
      }
    }
  }
  
  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/
  
  nLocalPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;
  
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
  
  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++ ) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
        
        /*--- Get the index of the current point. ---*/
        
        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/
        
        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
          if (newID < Global_Index) Global_Index = newID;
        }
        
        /*--- Search for the processor that owns this point. If we are
         sorting the elements, we use the linear partitioning to find
         the rank, otherwise, we simply have the current rank load its
         own elements into the connectivity data structure. ---*/
        
        if (val_sort) {
          iProcessor = Global_Index/npoint_procs[0];
          if (iProcessor >= (unsigned long)size)
            iProcessor = (unsigned long)size-1;
          if (Global_Index >= nPoint_Linear[iProcessor])
            while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
          else
            while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
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
  
  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
        
        /*--- Get the index of the current point. ---*/
        
        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/
        
        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
          if (newID < Global_Index) Global_Index = newID;
        }
       
        /*--- Search for the processor that owns this point. If we are
         sorting the elements, we use the linear partitioning to find
         the rank, otherwise, we simply have the current rank load its
         own elements into the connectivity data structure. ---*/
        
        if (val_sort) {
          iProcessor = Global_Index/npoint_procs[0];
          if (iProcessor >= (unsigned long)size)
            iProcessor = (unsigned long)size-1;
          if (Global_Index >= nPoint_Linear[iProcessor])
            while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
          else
            while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
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
            connSend[nn] = geometry->node[iPoint]->GetGlobalIndex(); nn++;
            
            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. ---*/
            
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
    case TRIANGLE:
      nParallel_Tria = nElem_Total;
      if (nParallel_Tria > 0) Conn_Tria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      nParallel_Quad = nElem_Total;
      if (nParallel_Quad > 0) Conn_Quad_Par = Conn_Elem;
      break;
    case TETRAHEDRON:
      nParallel_Tetr = nElem_Total;
      if (nParallel_Tetr > 0) Conn_Tetr_Par = Conn_Elem;
      break;
    case HEXAHEDRON:
      nParallel_Hexa = nElem_Total;
      if (nParallel_Hexa > 0) Conn_Hexa_Par = Conn_Elem;
      break;
    case PRISM:
      nParallel_Pris = nElem_Total;
      if (nParallel_Pris > 0) Conn_Pris_Par = Conn_Elem;
      break;
    case PYRAMID:
      nParallel_Pyra = nElem_Total;
      if (nParallel_Pyra > 0) Conn_Pyra_Par = Conn_Elem;
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
  delete [] Local_Halo;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic; 
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;

}
