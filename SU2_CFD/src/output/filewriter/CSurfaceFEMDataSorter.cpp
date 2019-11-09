#include "../../../include/output/filewriter/CSurfaceFEMDataSorter.hpp"
#include "../../../../Common/include/fem_geometry_structure.hpp"


CSurfaceFEMDataSorter::CSurfaceFEMDataSorter(CConfig *config, CGeometry *geometry, unsigned short nFields, CFEMDataSorter* volume_sorter) : CParallelDataSorter(config, nFields){
 
  this->volume_sorter = volume_sorter;
  
  connectivity_sorted = false;
  
  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  
  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();
  
  /*--- Update the solution by looping over the owned volume elements. ---*/
  
  for(unsigned long l=0; l<nVolElemOwned; ++l) {
    
    /* Count up the number of local points we have for allocating storage. */
    
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      nLocalPoint_Sort++;
    }
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nGlobalPoint_Sort = nLocalPoint_Sort;
#endif
  
  /*--- Create the linear partitioner --- */
  
  linearPartitioner = new CLinearPartitioner(nGlobalPoint_Sort, 0);

}

CSurfaceFEMDataSorter::~CSurfaceFEMDataSorter(){

  if (linearPartitioner != NULL) delete linearPartitioner;
  delete [] passiveDoubleBuffer;
  
}



void CSurfaceFEMDataSorter::SortOutputData() {
  
  if (!connectivity_sorted){
    SU2_MPI::Error("Connectivity must be sorted before sorting output data", CURRENT_FUNCTION);
  }
  
  const int VARS_PER_POINT = GlobalField_Counter;

  /*---------------------------------------------------*/
  /*--- Step 1: Determine the global DOF ID's of the   */
  /*---         locally stored surface connectivities. */
  /*---------------------------------------------------*/

  /* Loop over the surface connectivities and store global
     DOF ID's in a vector. Subtract 1, because the stored
     connectivities are 1 based. */
  globalSurfaceDOFIDs.clear();
  globalSurfaceDOFIDs.reserve(nParallel_Line*N_POINTS_LINE +
                              nParallel_Tria*N_POINTS_TRIANGLE +
                              nParallel_Quad*N_POINTS_QUADRILATERAL);

  for(unsigned long i=0; i<(nParallel_Line*N_POINTS_LINE); ++i) {
    const unsigned long globalID = Conn_Line_Par[i]-1;
    globalSurfaceDOFIDs.push_back(globalID);
  }

  for(unsigned long i=0; i<(nParallel_Tria*N_POINTS_TRIANGLE); ++i) {
    const unsigned long globalID = Conn_Tria_Par[i]-1;
    globalSurfaceDOFIDs.push_back(globalID);
  }

  for(unsigned long i=0; i<(nParallel_Quad*N_POINTS_QUADRILATERAL); ++i) {
    const unsigned long globalID = Conn_Quad_Par[i]-1;
    globalSurfaceDOFIDs.push_back(globalID);
  }

  /* Sort globalSurfaceDOFIDs in increasing order and remove the
     multiple entries. */
  sort(globalSurfaceDOFIDs.begin(), globalSurfaceDOFIDs.end());
  vector<unsigned long>::iterator lastEntry;
  lastEntry = unique(globalSurfaceDOFIDs.begin(), globalSurfaceDOFIDs.end());
  globalSurfaceDOFIDs.erase(lastEntry, globalSurfaceDOFIDs.end());

  /*-------------------------------------------------------------------*/
  /*--- Step 2: Communicate the information of globalSurfaceDOFIDs  ---*/
  /*---         to the ranks, which actually store this information ---*/
  /*---         in Parallel_Data.                                   ---*/
  /*-------------------------------------------------------------------*/

  /* Allocate the memory for the first index of the communication buffers. */
  vector<vector<unsigned long> > sendBuf(size, vector<unsigned long>(0));

  /* Loop over the entries of globalSurfaceDOFIDs and fill the
     communication buffers accordingly. */
  for(unsigned long i=0; i<globalSurfaceDOFIDs.size(); ++i) {

    /* Search for the processor that owns this point. */
    unsigned long iProcessor = linearPartitioner->GetRankContainingIndex(globalSurfaceDOFIDs[i]);

    /* Store the global ID in the send buffer for iProcessor. */
    sendBuf[iProcessor].push_back(globalSurfaceDOFIDs[i]);
  }

  /* Determine the number of DOFs to be sent to each processor. */
  int nRankSend = 0;
  vector<unsigned long> nDOFSend(size);
  for(int i=0; i<size; ++i) {
    nDOFSend[i] = sendBuf[i].size();
    if(nDOFSend[i] && (i != rank)) ++nRankSend;
  }

  /* Communicate nDOFSend using Alltoall. */
  vector<unsigned long> nDOFRecv(size);

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(nDOFSend.data(), 1, MPI_UNSIGNED_LONG,
                    nDOFRecv.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  nDOFRecv[rank] = nDOFSend[rank];
#endif

  /* Determine the number of messages this rank will receive. */
  int nRankRecv = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFRecv[i] && (i != rank)) ++nRankRecv;
  }

  /* Define the receive buffers. */
  vector<vector<unsigned long> > recvBuf(size, vector<unsigned long>(0));

#ifdef HAVE_MPI
  /* Launch the non-blocking sends. Do not send anything to myself. */
  vector<SU2_MPI::Request> sendReq(nRankSend);
  nRankSend = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFSend[i] && (i != rank)) {
      SU2_MPI::Isend(sendBuf[i].data(), nDOFSend[i], MPI_UNSIGNED_LONG,
                     i, rank, MPI_COMM_WORLD, &sendReq[nRankSend]);
      ++nRankSend;
    }
  }

  /* Launch the non-blocking receives. Do not receive anything from myself. */
  vector<SU2_MPI::Request> recvReq(nRankRecv);
  nRankRecv = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFRecv[i] && (i != rank)) {
      recvBuf[i].resize(nDOFRecv[i]);
      SU2_MPI::Irecv(recvBuf[i].data(), nDOFRecv[i], MPI_UNSIGNED_LONG,
                     i, i, MPI_COMM_WORLD, &recvReq[nRankRecv]);
      ++nRankRecv;
    }
  }

#endif

  /* Copy my own data. */
  recvBuf[rank] = sendBuf[rank];

#ifdef HAVE_MPI
  /* Complete the non-blocking communication. */
  SU2_MPI::Waitall(nRankSend, sendReq.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Waitall(nRankRecv, recvReq.data(), MPI_STATUSES_IGNORE);
#endif

  /*-------------------------------------------------------------------*/
  /*--- Step 2: Determine the number of locally stored surface DOFs ---*/
  /*---         and store the data in Parallel_Surf_Data.           ---*/
  /*-------------------------------------------------------------------*/

  /* Store all received surface ID in one vector and sort it. */
  globalSurfaceDOFIDs.resize(0);
  for(int i=0; i<size; ++i) {
    if( nDOFRecv[i] )
      globalSurfaceDOFIDs.insert(globalSurfaceDOFIDs.end(),
                                  recvBuf[i].begin(), recvBuf[i].end());
  }

  sort(globalSurfaceDOFIDs.begin(), globalSurfaceDOFIDs.end());
  lastEntry = unique(globalSurfaceDOFIDs.begin(), globalSurfaceDOFIDs.end());
  globalSurfaceDOFIDs.erase(lastEntry, globalSurfaceDOFIDs.end());

  /* Allocate the memory for Parallel_Surf_Data. */
  nParallel_Poin = globalSurfaceDOFIDs.size();
  
  if (passiveDoubleBuffer == nullptr){
    passiveDoubleBuffer = new passivedouble[nParallel_Poin*VARS_PER_POINT];
  }
  
  /* Determine the local index of the global surface DOFs and
     copy the data into Parallel_Surf_Data. */
  for(unsigned long i=0; i<nParallel_Poin; ++i) {
    const unsigned long ii = globalSurfaceDOFIDs[i] - linearPartitioner->GetCumulativeSizeBeforeRank(rank);

    for(int jj=0; jj<VARS_PER_POINT; jj++)
      passiveDoubleBuffer[i*VARS_PER_POINT+jj] = volume_sorter->GetData(jj,ii);
  }

  /*--- Reduce the total number of surf points we have. This will be
        needed for writing the surface solution files later. ---*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nGlobal_Poin_Par = nParallel_Poin;
#endif

  /*-------------------------------------------------------------------*/
  /*--- Step 3: Modify the surface connectivities, such that only   ---*/
  /*---         the data of surface DOFs needs to be written.       ---*/
  /*-------------------------------------------------------------------*/

  /* Determine the offset for my surface DOFs. */
  unsigned long offsetSurfaceDOFs = 0;
#ifdef HAVE_MPI
  vector<unsigned long> nSurfaceDOFsRanks(size);

  SU2_MPI::Allgather(&nParallel_Poin, 1, MPI_UNSIGNED_LONG,
                     nSurfaceDOFsRanks.data(), 1, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);

  for(int i=0; i<rank; ++i) offsetSurfaceDOFs += nSurfaceDOFsRanks[i];
#endif

  /* Create the map from the global volume numbering to the new global
     surface numbering. */
  map<unsigned long, unsigned long> mapGlobalVol2Surf;
  for(unsigned long i=0; i<nParallel_Poin; ++i)
    mapGlobalVol2Surf[globalSurfaceDOFIDs[i]] = offsetSurfaceDOFs + i;
   

  /* Fill the receive buffers with the modified global surface DOF numbers,
     such that this information can be returned to the original ranks. */
  for(int i=0; i<size; ++i) {
    for(unsigned long j=0; j<nDOFRecv[i]; ++j)
      recvBuf[i][j] = mapGlobalVol2Surf.find(recvBuf[i][j])->second;
  }

  /* Copy the original send buffers, because that information is
     still needed. */
  vector<vector<unsigned long> > originalSendBuf = sendBuf;

#ifdef HAVE_MPI
  /* Launch the non-blocking sends for the reverse communication, where the
     receive buffers must be used for sending. Do not send anything to myself. */
  nRankRecv = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFRecv[i] && (i != rank)) {
      SU2_MPI::Isend(recvBuf[i].data(), nDOFRecv[i], MPI_UNSIGNED_LONG,
                     i, rank+1, MPI_COMM_WORLD, &recvReq[nRankRecv]);
      ++nRankRecv;
    }
  }

  /* Launch the non-blocking receives for the reverse communication, where the
     send buffers must be used for receiving.  Do not receive anything from myself. */
  nRankSend = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFSend[i] && (i != rank)) {
      SU2_MPI::Irecv(sendBuf[i].data(), nDOFSend[i], MPI_UNSIGNED_LONG,
                     i, i+1, MPI_COMM_WORLD, &sendReq[nRankSend]);
      ++nRankSend;
    }
  }

#endif

  /* Copy my own data. */
  sendBuf[rank] = recvBuf[rank];

#ifdef HAVE_MPI
  /* Complete the non-blocking communication. */
  SU2_MPI::Waitall(nRankRecv, recvReq.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Waitall(nRankSend, sendReq.data(), MPI_STATUSES_IGNORE);
#endif

  /* Clear the map mapGlobalVol2Surf and fill it with the data
     needed to modify the surface connectivity. Note that 1 is added,
     because the visualization softwares typically use 1 based indexing. */
  mapGlobalVol2Surf.clear();

  for(int i=0; i<size; ++i) {
    for(unsigned long j=0; j<nDOFSend[i]; ++j)
      mapGlobalVol2Surf[originalSendBuf[i][j]+1] = sendBuf[i][j]+1;
  }

  /* Modify the locally stored surface connectivities. */
  for(unsigned long i=0; i<(nParallel_Line*N_POINTS_LINE); ++i)
    Conn_Line_Par[i] = mapGlobalVol2Surf.find(Conn_Line_Par[i])->second;

  for(unsigned long i=0; i<(nParallel_Tria*N_POINTS_TRIANGLE); ++i)
    Conn_Tria_Par[i] = mapGlobalVol2Surf.find(Conn_Tria_Par[i])->second;

  for(unsigned long i=0; i<(nParallel_Quad*N_POINTS_QUADRILATERAL); ++i)
    Conn_Quad_Par[i] = mapGlobalVol2Surf.find(Conn_Quad_Par[i])->second;
}

void CSurfaceFEMDataSorter::SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) {
  
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



void CSurfaceFEMDataSorter::SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  /* Determine the number of nodes for this element type. */
    unsigned short NODES_PER_ELEMENT = 0;
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
    }
    
    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
          geometrical information for the FEM DG solver. ---*/
    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  
    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();
  
    const CBoundaryFEM *boundaries = DGGeometry->GetBoundaries();
    const CFEMStandardBoundaryFace *standardBoundaryFacesSol = DGGeometry->GetStandardBoundaryFacesSol();
  
    /*--- Create the map from the global DOF ID to the local index.
          Note one is added to the index value, because visualization
          softwares typically use 1-based indexing. ---*/
    vector<unsigned long> globalID;
    for(unsigned long l=0; l<nVolElemOwned; ++l) {
      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
        const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j + 1;
        globalID.push_back(globalIndex);
      }
    }
  
    /*--- Determine the number of sub-elements on this rank by looping
          over the surface elements of the boundary markers that must
          be plotted. ---*/
    unsigned long nSubElem_Local = 0;
    for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); ++iMarker) {
      if( !boundaries[iMarker].periodicBoundary ) {
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          const vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
          for(unsigned long i=0; i<surfElem.size(); ++i) {
            const unsigned short ind      = surfElem[i].indStandardElement;
            const unsigned short VTK_Type = standardBoundaryFacesSol[ind].GetVTK_Type();
            if(Elem_Type == VTK_Type) nSubElem_Local += standardBoundaryFacesSol[ind].GetNSubFaces();
          }
        }
      }
    }
  
    /* Allocate the memory to store the connectivity if the size is
       larger than zero. */
    int *Conn_SubElem = NULL;
    if(nSubElem_Local > 0) Conn_SubElem = new int[nSubElem_Local*NODES_PER_ELEMENT]();
  
    /*--- Repeat the loop over the surface elements of the boundary markers
          that must be plotted, but now store the connectivity. ---*/
    unsigned long kNode = 0;
    for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); ++iMarker) {
      if( !boundaries[iMarker].periodicBoundary ) {
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          const vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
  
          /* Loop over the surface elements of this boundary marker. */
          for(unsigned long i=0; i<surfElem.size(); ++i) {
  
            /* Check if this is the element type to be stored. */
            const unsigned short ind      = surfElem[i].indStandardElement;
            const unsigned short VTK_Type = standardBoundaryFacesSol[ind].GetVTK_Type();
            if(Elem_Type == VTK_Type) {
  
              /* Get the number of sub-elements and the local connectivity of
                 the sub-elements. */
              const unsigned short nSubFaces     = standardBoundaryFacesSol[ind].GetNSubFaces();
              const unsigned short *connSubFaces = standardBoundaryFacesSol[ind].GetSubFaceConn();
  
              /* Store the global connectivities. */
              const unsigned short kk = NODES_PER_ELEMENT*nSubFaces;
              for(unsigned short k=0; k<kk; ++k, ++kNode)
                Conn_SubElem[kNode] = globalID[surfElem[i].DOFsSolFace[connSubFaces[k]]];
            }
          }
        }
      }
    }
  
    /*--- Store the particular global element count in the class data,
          and set the class data pointer to the connectivity array. ---*/
    switch (Elem_Type) {
      case LINE:
        nParallel_Line = nSubElem_Local;
        if (Conn_Line_Par != NULL) delete [] Conn_Line_Par;
        Conn_Line_Par = Conn_SubElem;
        break;
      case TRIANGLE:
        nParallel_Tria = nSubElem_Local;
        if (Conn_Tria_Par != NULL) delete [] Conn_Tria_Par;        
        Conn_Tria_Par = Conn_SubElem;
        break;
      case QUADRILATERAL:
        nParallel_Quad = nSubElem_Local;
        if (Conn_Quad_Par != NULL) delete [] Conn_Quad_Par;       
        Conn_Quad_Par = Conn_SubElem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
        break;
    }
    
}
