#include "../../../include/output/filewriter/CFEMDataSorter.hpp"
#include "../../../Common/include/fem_geometry_structure.hpp"

CFEMDataSorter::CFEMDataSorter(CConfig *config, CGeometry *geometry, unsigned short nFields) : CParallelDataSorter(config, nFields){
 
  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  
  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();
  
  /*--- Create the map from the global DOF ID to the local index. ---*/

  vector<unsigned long> globalID;
  
  /*--- Update the solution by looping over the owned volume elements. ---*/
  
  for(unsigned long l=0; l<nVolElemOwned; ++l) {
    
    /* Count up the number of local points we have for allocating storage. */
    
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);
      
      nLocalPoint_Sort++;
    }
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nGlobalPoint_Sort = nLocalPoint_Sort;
#endif

  /*--- Create a linear partition --- */
  
  linearPartitioner = new CLinearPartitioner(nGlobalPoint_Sort, 0);
  
  /*--- Prepare the send buffers ---*/
  
  PrepareSendBuffers(globalID);
  
}

CFEMDataSorter::~CFEMDataSorter(){

  if (connSend != NULL)    delete [] connSend;
  if (Index != NULL)       delete [] Index;
  if (idSend != NULL)      delete [] idSend;
  if (linearPartitioner != NULL) delete linearPartitioner;
  
}


void CFEMDataSorter::SortOutputData() {

  /* For convenience, set the total number of variables stored at each DOF. */

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

  su2double *connRecv = NULL;
  connRecv = new su2double[VARS_PER_POINT*nPoint_Recv[size]]();
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Recv[size]; ii++)
    connRecv[ii] = 0;

  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]]();
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
  
  delete [] connSend;
  connSend = NULL;

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate the
   appropriate amount of memory for this section. ---*/

  Parallel_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Data[jj] = new su2double[nPoint_Recv[size]]();
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

  delete [] connRecv;
  delete [] idRecv;
  
}


void CFEMDataSorter::SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) {

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/
  
  /*--- Sort volumetric grid connectivity. ---*/

  if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
    cout <<"Sorting volumetric grid connectivity." << endl;
  
  SortVolumetricConnectivity(config, geometry, TRIANGLE     );
  SortVolumetricConnectivity(config, geometry, QUADRILATERAL);
  SortVolumetricConnectivity(config, geometry, TETRAHEDRON  );
  SortVolumetricConnectivity(config, geometry, HEXAHEDRON   );
  SortVolumetricConnectivity(config, geometry, PRISM        );
  SortVolumetricConnectivity(config, geometry, PYRAMID      );
  
  
  /*--- Reduce the total number of cells we will be writing in the output files. ---*/
  
  unsigned long nTotal_Elem = nParallel_Tria + nParallel_Quad + nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;
#ifndef HAVE_MPI
  nGlobal_Elem_Par = nTotal_Elem;
#else
  SU2_MPI::Allreduce(&nTotal_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  connectivity_sorted = true;
  
}

void CFEMDataSorter::SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  /* Determine the number of nodes for this element type. */
  unsigned short NODES_PER_ELEMENT = 0;
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

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  const CFEMStandardElement *standardElementsSol = DGGeometry->GetStandardElementsSol();

  /*--- Determine the number of sub-elements on this rank. ---*/
  unsigned long nSubElem_Local = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Determine the necessary data from the corresponding standard elem. */
    const unsigned short ind       = volElem[i].indStandardElement;
    const unsigned short VTK_Type1 = standardElementsSol[ind].GetVTK_Type1();
    const unsigned short VTK_Type2 = standardElementsSol[ind].GetVTK_Type2();

     /* Only store the linear sub elements if they are of
        the current type that we are storing. */
     if(Elem_Type == VTK_Type1) nSubElem_Local += standardElementsSol[ind].GetNSubElemsType1();
     if(Elem_Type == VTK_Type2) nSubElem_Local += standardElementsSol[ind].GetNSubElemsType2();
  }

  /* Allocate the memory to store the connectivity if the size is
     larger than zero. */
  int *Conn_SubElem = NULL;
  if(nSubElem_Local > 0) Conn_SubElem = new int[nSubElem_Local*NODES_PER_ELEMENT];

  /*--- Loop again over the local volume elements and store the global
        connectivities of the sub-elements. Note one is added to the
        index value, because visualization softwares typically use
        1-based indexing. ---*/
  unsigned long kNode = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Determine the necessary data from the corresponding standard elem. */
    const unsigned short ind       = volElem[i].indStandardElement;
    const unsigned short VTK_Type1 = standardElementsSol[ind].GetVTK_Type1();
    const unsigned short VTK_Type2 = standardElementsSol[ind].GetVTK_Type2();

    /* Check if the first sub-element is of the required type. */
    if(Elem_Type == VTK_Type1) {

      /* Get the number of sub-elements and the local connectivity of
         the sub-elements. */
      const unsigned short nSubElems     = standardElementsSol[ind].GetNSubElemsType1();
      const unsigned short *connSubElems = standardElementsSol[ind].GetSubConnType1();

      /* Store the global connectivities. */
      const unsigned short kk = NODES_PER_ELEMENT*nSubElems;
      for(unsigned short k=0; k<kk; ++k, ++kNode)
        Conn_SubElem[kNode] = connSubElems[k] + volElem[i].offsetDOFsSolGlobal + 1;
    }

    /* Check if the second sub-element is of the required type. */
    if(Elem_Type == VTK_Type2) {

      /* Get the number of sub-elements and the local connectivity of
         the sub-elements. */
      const unsigned short nSubElems     = standardElementsSol[ind].GetNSubElemsType2();
      const unsigned short *connSubElems = standardElementsSol[ind].GetSubConnType2();

      /* Store the global connectivities. */
      const unsigned short kk = NODES_PER_ELEMENT*nSubElems;
      for(unsigned short k=0; k<kk; ++k, ++kNode)
        Conn_SubElem[kNode] = connSubElems[k] + volElem[i].offsetDOFsSolGlobal + 1;
    }
  }

  /*--- Store the particular global element count in the class data,
        and set the class data pointer to the connectivity array. ---*/
  switch (Elem_Type) {
    case TRIANGLE:
      nParallel_Tria = nSubElem_Local;
      Conn_Tria_Par = Conn_SubElem;
      break;
    case QUADRILATERAL:
      nParallel_Quad = nSubElem_Local;
      Conn_Quad_Par = Conn_SubElem;
      break;
    case TETRAHEDRON:
      nParallel_Tetr = nSubElem_Local;
      Conn_Tetr_Par = Conn_SubElem;
      break;
    case HEXAHEDRON:
      nParallel_Hexa = nSubElem_Local;
      Conn_Hexa_Par = Conn_SubElem;
      break;
    case PRISM:
      nParallel_Pris = nSubElem_Local;
      Conn_Pris_Par = Conn_SubElem;
      break;
    case PYRAMID:
      nParallel_Pyra = nSubElem_Local;
      Conn_Pyra_Par = Conn_SubElem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }
}
