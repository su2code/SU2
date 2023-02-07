/*!
 * \file CCGNSMeshReaderFVM.cpp
 * \brief Class that reads a single zone of a CGNS mesh file from disk into
 *        linear partitions across all ranks.
 * \author T. Economon
 * \version 7.5.1 "Blackbird"
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CCGNSMeshReaderFVM.hpp"

CCGNSMeshReaderFVM::CCGNSMeshReaderFVM(CConfig        *val_config,
                                       unsigned short val_iZone,
                                       unsigned short val_nZone)
: CCGNSMeshReaderBase(val_config, val_iZone, val_nZone) {

#ifdef HAVE_CGNS
  OpenCGNSFile(config->GetMesh_FileName());

  /*--- Read the basic information about the database and zone(s). ---*/
  ReadCGNSDatabaseMetadata();
  ReadCGNSZoneMetadata();

  /*--- Read the point coordinates into linear partitions. ---*/
  ReadCGNSPointCoordinates();

  /*--- Loop over all sections to access the grid connectivity. We
   treat the interior and boundary elements with separate routines.
   If we have found that this is a boundary section (we assume
   that internal cells and boundary cells do not exist in the same
   section together), the master node reads the boundary section.
   Otherwise, all ranks read and communicate the interior sections. ---*/
  ReadCGNSSectionMetadata();
  numberOfMarkers = 0;
  for (int s = 0; s < nSections; s++) {
    if (isInterior[s]) {
      ReadCGNSVolumeSection(s);
    } else {
      numberOfMarkers++;
      ReadCGNSSurfaceSection(s);
    }
  }

  /*--- We have extracted all CGNS data. Close the CGNS file. ---*/
  if (cg_close(cgnsFileID)) cg_error_exit();

  /*--- Put our CGNS data into the class data for the mesh reader. ---*/
  ReformatCGNSVolumeConnectivity();
  ReformatCGNSSurfaceConnectivity();
#else
  SU2_MPI::Error(string(" SU2 built without CGNS support. \n") +
                 string(" To use CGNS, build SU2 accordingly."),
                 CURRENT_FUNCTION);
#endif
}

#ifdef HAVE_CGNS
void CCGNSMeshReaderFVM::ReadCGNSVolumeSection(int val_section) {

  /*--- In this routine, each rank will read a chunk of the element
   connectivity for a single specified section of the CGNS mesh file.
   All operations are executed in parallel here and the reading of the
   section proceeds based on a linear partitioning of the elements
   across all ranks in the calculation. We will use partial reads of
   the CGNS section from the CGNS API to accomplish this. Once each
   rank has a linear chunk of the mesh, we will redistribute the
   connectivity to match the linear partitioning of the grid points,
   not the elements, since the points control the overall partitioning. ---*/

  int nbndry, parent_flag, npe, iProcessor;
  unsigned long iElem = 0, iPoint = 0, iNode = 0, jNode = 0;
  cgsize_t startE, endE;
  ElementType_t elemType;
  char sectionName[CGNS_STRING_SIZE];

  /*--- Read the connectivity details for this section.
   Store the total number of elements in this section
   to be used later for memory allocation. ---*/

  if (cg_section_read(cgnsFileID, cgnsBase, cgnsZone, val_section+1, sectionName,
                      &elemType, &startE, &endE, &nbndry, &parent_flag))
    cg_error_exit();

  /*--- Compute element linear partitioning ---*/

  unsigned long element_count = (endE-startE+1);
  CLinearPartitioner elementPartitioner(element_count,startE,true);

  /*--- Store the number of elements that this rank is responsible for
   in the current section. ---*/

  nElems[val_section] = elementPartitioner.GetSizeOnRank(rank);

  /*--- Allocate some memory for the handling the connectivity
   and auxiliary data that we need to communicate. ---*/

  vector<cgsize_t> elemTypes(nElems[val_section],   0);
  vector<cgsize_t> nPoinPerElem(nElems[val_section],0);
  vector<cgsize_t> elemGlobalID(nElems[val_section],0);

  /*--- Determine the size of the vector needed to read the connectivity
   data from the CGNS file. Only call the CGNS API if we have a non-zero
   number of elements on this rank. ---*/

  cgsize_t sizeNeeded = 0, sizeOffset = 0;
  if (nElems[val_section] > 0) {
    if (cg_ElementPartialSize(cgnsFileID, cgnsBase, cgnsZone, val_section+1,
                              (cgsize_t)elementPartitioner.GetFirstIndexOnRank(rank),
                              (cgsize_t)elementPartitioner.GetLastIndexOnRank(rank),
                              &sizeNeeded) != CG_OK)
    cg_error_exit();
  }

  /*--- Allocate the memory for the connectivity, the offset if needed
   and read the data. ---*/

  vector<cgsize_t> connElemCGNS(sizeNeeded,0);
  if (elemType == MIXED || elemType == NFACE_n || elemType == NGON_n) {
    sizeOffset = nElems[val_section]+1;
  }
  vector<cgsize_t> connOffsetCGNS(sizeOffset,0);

  /*--- Retrieve the connectivity information and store. Note that
   we are only accessing our rank's piece of the data here in the
   partial read function in the CGNS API. Only call the CGNS API
   if we have a non-zero number of elements on this rank. ---*/

  if (nElems[val_section] > 0) {
    if (elemType == MIXED || elemType == NFACE_n || elemType == NGON_n) {
      if (cg_poly_elements_partial_read(cgnsFileID, cgnsBase, cgnsZone, val_section+1,
                                        (cgsize_t)elementPartitioner.GetFirstIndexOnRank(rank),
                                        (cgsize_t)elementPartitioner.GetLastIndexOnRank(rank),
                                        connElemCGNS.data(),
                                        connOffsetCGNS.data(), NULL) != CG_OK)
      cg_error_exit();
    } else {
      if (cg_elements_partial_read(cgnsFileID, cgnsBase, cgnsZone, val_section+1,
                                   (cgsize_t)elementPartitioner.GetFirstIndexOnRank(rank),
                                   (cgsize_t)elementPartitioner.GetLastIndexOnRank(rank),
                                   connElemCGNS.data(), NULL) != CG_OK)
      cg_error_exit();
    }
  }

  /*--- Print some information to the console. ---*/

  if (rank == MASTER_NODE) {
    cout << "Loading volume section " << string(sectionName);
    cout <<  " from file." << endl;
  }

  /*--- Find the number of nodes required to represent
   this type of element. ---*/

  if (cg_npe(elemType, &npe)) cg_error_exit();

  /*--- Check whether the sections contains a mixture of multiple
   element types, which will require special handling to get the
   element type one-by-one when reading. ---*/

  bool isMixed = false;
  if (elemType == MIXED) {
    isMixed = true;
  }

  /*--- Loop through all of the elements in this section to get more
   information and to decide whether it has interior elements. ---*/

  unsigned long counterCGNS = 0;
  for (iElem = 0; iElem < nElems[val_section]; iElem++) {

    ElementType_t iElemType = elemType;

    /*--- If we have a mixed element section, we need to check the elem
     type one-by-one. We also must manually advance the counter to the
     next element's position in the buffer. ---*/

    if (isMixed) {
      iElemType = ElementType_t(connElemCGNS[counterCGNS]);
      npe       = connOffsetCGNS[iElem+1]-connOffsetCGNS[iElem]-1;
      counterCGNS++;
      for (int jj = 0; jj < npe; jj++) counterCGNS++;
    }

    /*--- Store the number of points per element for the current elem. ---*/

    nPoinPerElem[iElem] = npe;

    /*--- Store the global ID for this element. Note the -1 to move
     from CGNS convention to SU2 convention. We also subtract off
     an additional offset in case we have found boundary sections
     prior to this one, in order to keep the internal element global
     IDs indexed starting from zero. ---*/

    elemGlobalID[iElem] = (elementPartitioner.GetFirstIndexOnRank(rank) +
                           iElem - elemOffset[val_section]);

    /* Get the VTK type for this element. */

    int vtk_type;
    string elem_name = GetCGNSElementType(iElemType, vtk_type);
    elemTypes[iElem] = vtk_type;

  }

  /*--- Force free the memory for the conn offset from the CGNS file. ---*/

  vector<cgsize_t>().swap(connOffsetCGNS);

  /*--- These are internal elems. Allocate memory on each proc. ---*/

  vector<cgsize_t> connElemTemp(nElems[val_section]*SU2_CONN_SIZE,0);

  /*--- Copy the connectivity into the larger array with a standard
   format per element: [globalID vtkType n0 n1 n2 n3 n4 n5 n6 n7 n8]. ---*/

  counterCGNS = 0;
  for (iElem = 0; iElem < nElems[val_section]; iElem++) {

    /*--- Store the conn in chunks of SU2_CONN_SIZE for simplicity. ---*/

    unsigned long nn = iElem*SU2_CONN_SIZE;

    /*--- First, store the global element ID and the VTK type. ---*/

    connElemTemp[nn] = elemGlobalID[iElem]; nn++;
    connElemTemp[nn] = elemTypes[iElem];    nn++;

    /*--- Store the connectivity values. Note we subtract one from
     the CGNS 1-based convention. We may also need to remove the first
     entry if this is a mixed element section. ---*/

    if (isMixed) counterCGNS++;
    for (iNode = 0; iNode < (unsigned long)nPoinPerElem[iElem]; iNode++) {
      connElemTemp[nn] = connElemCGNS[counterCGNS + iNode] - 1; nn++;
    }
    counterCGNS += nPoinPerElem[iElem];

  }

  /*--- Force free the memory for the conn from the CGNS file. ---*/

  vector<cgsize_t>().swap(connElemCGNS);

  /*--- We now have the connectivity stored in linearly partitioned
   chunks. We need to loop through and decide how many elements we
   must send to each rank in order to have all elements that
   surround a particular "owned" node on each rank (i.e., elements
   will appear on multiple ranks). First, initialize a counter
   and flag. ---*/

  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    nElem_Send[iProcessor] = 0;
    nElem_Recv[iProcessor] = 0;
    nElem_Flag[iProcessor]= -1;
  }
  nElem_Send[size] = 0;
  nElem_Recv[size] = 0;

  /*--- Create a partitioner object to find the owning rank of points. ---*/

  CLinearPartitioner pointPartitioner(numberOfGlobalPoints,0);

  for (iElem = 0; iElem < nElems[val_section]; iElem++) {
    for (iNode = 0; iNode < (unsigned long)nPoinPerElem[iElem]; iNode++) {

      /*--- Get the index of the current point. ---*/

      iPoint = connElemTemp[iElem*SU2_CONN_SIZE + SU2_CONN_SKIP + iNode];

      /*--- Search for the processor that owns this point. ---*/

      iProcessor = pointPartitioner.GetRankContainingIndex(iPoint);

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nElem_Flag[iProcessor] != (int)iElem)) {
        nElem_Flag[iProcessor] = iElem;
        nElem_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, SU2_MPI::GetComm());

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  unsigned long nSends = 0, nRecvs = 0;
  for (iProcessor = 0; iProcessor < size; iProcessor++)
    nElem_Flag[iProcessor] = -1;

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if ((iProcessor != rank) && (nElem_Send[iProcessor+1] > 0)) nSends++;
    if ((iProcessor != rank) && (nElem_Recv[iProcessor+1] > 0)) nRecvs++;

    nElem_Send[iProcessor+1] += nElem_Send[iProcessor];
    nElem_Recv[iProcessor+1] += nElem_Recv[iProcessor];
  }

  /*--- Allocate memory to hold the connectivity that we are
   sending. Note that we are also sending the global ID and the
   VTK element type in the first and second positions, respectively.
   We have assumed a constant message size of a hex element (8 nodes)
   + 2 extra values for the ID and VTK. ---*/

  unsigned long *connSend = NULL, iSend = 0;
  unsigned long sendSize = (unsigned long)SU2_CONN_SIZE*nElem_Send[size];
  connSend = new unsigned long[sendSize];
  for (iSend = 0; iSend < sendSize; iSend++)
    connSend[iSend] = 0;

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  vector<unsigned long> index(size);
  for (iProcessor = 0; iProcessor < size; iProcessor++)
    index[iProcessor] = SU2_CONN_SIZE*nElem_Send[iProcessor];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iElem = 0; iElem < (unsigned long)nElems[val_section]; iElem++) {
    for (iNode = 0; iNode < (unsigned long)nPoinPerElem[iElem]; iNode++) {

      /*--- Get the index of the current point. ---*/

      iPoint = connElemTemp[iElem*SU2_CONN_SIZE + SU2_CONN_SKIP + iNode];

      /*--- Search for the processor that owns this point ---*/

      iProcessor = pointPartitioner.GetRankContainingIndex(iPoint);

      /*--- Load connectivity into the buffer for sending ---*/

      if (nElem_Flag[iProcessor] != (int)iElem) {

        nElem_Flag[iProcessor] = iElem;
        unsigned long nn = index[iProcessor];

        /*--- Load the VTK type first into the conn array,
         then the connectivity vals, and last, the global ID. ---*/

        for (jNode = 0; jNode < SU2_CONN_SIZE; jNode++) {
          connSend[nn] = connElemTemp[iElem*SU2_CONN_SIZE + jNode]; nn++;
        }

        /*--- Increment the index by the message length ---*/

        index[iProcessor] += SU2_CONN_SIZE;

      }

    }
  }

  /*--- Force free memory after loading up the send buffer. ---*/

  vector<cgsize_t>().swap(connElemTemp);
  vector<cgsize_t>().swap(elemTypes);
  vector<cgsize_t>().swap(nPoinPerElem);
  vector<cgsize_t>().swap(elemGlobalID);

  vector<unsigned long>().swap(index);

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *connRecv = NULL, iRecv = 0;
  unsigned long recvSize = (unsigned long)SU2_CONN_SIZE*nElem_Recv[size];
  connRecv = new unsigned long[recvSize];
  for (iRecv = 0; iRecv < recvSize; iRecv++)
    connRecv[iRecv] = 0;

  /*--- Allocate memory for the MPI requests if we will communicate. ---*/

  SU2_MPI::Request *connSendReq = NULL;
  SU2_MPI::Request *connRecvReq = NULL;

  if (nSends > 0) {
    connSendReq = new SU2_MPI::Request[nSends];
  }
  if (nRecvs > 0) {
    connRecvReq = new SU2_MPI::Request[nRecvs];
  }

  /*--- Launch the non-blocking sends and receives. ---*/

  InitiateCommsAll(connSend, nElem_Send, connSendReq,
                   connRecv, nElem_Recv, connRecvReq,
                   SU2_CONN_SIZE, COMM_TYPE_UNSIGNED_LONG);

  /*--- Copy the current rank's data into the recv buffer directly. ---*/

  iRecv = SU2_CONN_SIZE*nElem_Recv[rank];
  unsigned long myStart = SU2_CONN_SIZE*nElem_Send[rank];
  unsigned long myFinal = SU2_CONN_SIZE*nElem_Send[rank+1];
  for (iSend = myStart; iSend < myFinal; iSend++) {
    connRecv[iRecv] = connSend[iSend];
    iRecv++;
  }

  /*--- Complete the non-blocking communications. ---*/

  CompleteCommsAll(nSends, connSendReq, nRecvs, connRecvReq);

  /*--- Store the connectivity for this rank in the proper data
   structure. First, allocate the appropriate amount of memory
   for this section, then write the recv'd values. ---*/

  if (nElem_Recv[size] > 0) {
    connElems[val_section].resize(nElem_Recv[size]*SU2_CONN_SIZE,0);
    unsigned long count = 0;
    for (iElem = 0; iElem < (unsigned long)nElem_Recv[size]; iElem++) {
      for (iNode = 0; iNode < SU2_CONN_SIZE; iNode++) {
        unsigned long nn = iElem*SU2_CONN_SIZE+iNode;
        connElems[val_section][count] = (cgsize_t)connRecv[nn];
        count++;
      }
    }

    /*--- Store the total number of elements the current rank
     now has for the current section after completing the comms. ---*/

    nElems[val_section] = nElem_Recv[size];

  } else {

    /*--- The current rank did not recv any elements from this
     section. Set the count to zero and nullify the data structure. ---*/

    nElems[val_section]    = 0;
    connElems[val_section].resize(0);

  }

  /*--- Free temporary memory from communications ---*/

  if (connSendReq != NULL) delete [] connSendReq;
  if (connRecvReq != NULL) delete [] connRecvReq;

  delete [] connSend;
  delete [] connRecv;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;

}

void CCGNSMeshReaderFVM::ReadCGNSSurfaceSection(int val_section) {

  /*--- In this routine, we access a CGNS surface section and have the
   master rank load all of the surface conn. This can help avoid issues
   where there are fewer elements than ranks on a surface. This is later
   linearly partitioned. A limitation of this approach is that there
   could be a memory bottleneck for extremely large grids, but this has
   not been reached yet in practice. In that event, we can also read these
   surface sections with a linear partitioning. ---*/

  int nbndry, parent_flag, npe;
  unsigned long iElem = 0, iNode = 0;
  cgsize_t startE, endE, ElementDataSize;
  ElementType_t elemType;
  char sectionName[CGNS_STRING_SIZE];

  if (rank == MASTER_NODE) {

    /*--- Allocate some memory for the handling the connectivity
     and auxiliary data that we are need to communicate. ---*/

    vector<cgsize_t> connElemCGNS(nElems[val_section]*SU2_CONN_SIZE,0);
    vector<unsigned short> elemTypes(nElems[val_section],0);
    vector<unsigned short> nPoinPerElem(nElems[val_section],0);
    vector<unsigned long>  elemGlobalID(nElems[val_section],0);

    /*--- Read the section info again ---*/

    if (cg_section_read(cgnsFileID, cgnsBase, cgnsZone, val_section+1,
                        sectionName, &elemType, &startE, &endE, &nbndry,
                        &parent_flag))
      cg_error_exit();

    /*--- Print some information to the console. ---*/

    cout << "Loading surface section " << string(sectionName);
    cout <<  " from file." << endl;

    /*--- Store the number of elems (all on the master). ---*/

    nElems[val_section] = (endE-startE+1);

    /*--- Read and store the total amount of data that will be
     listed when reading this section. ---*/

    if (cg_ElementDataSize(cgnsFileID, cgnsBase, cgnsZone, val_section+1,
                           &ElementDataSize))
      cg_error_exit();

    /*--- Find the number of nodes required to represent
     this type of element. ---*/

    if (cg_npe(elemType, &npe)) cg_error_exit();

    /*--- Check whether the sections contains a mixture of multiple
     element types, which will require special handling to get the
     element type one-by-one when reading. ---*/

    bool isMixed = false;
    if (elemType == MIXED) {
      isMixed = true;
    }

    /*--- Allocate memory for accessing the connectivity and to
     store it in the proper data structure for post-processing. ---*/

    vector<cgsize_t> connElemTemp(ElementDataSize,0);

    /*--- Retrieve the connectivity information and store. ---*/

    if (elemType == MIXED || elemType == NGON_n || elemType == NFACE_n) {
      vector<cgsize_t> connOffsetTemp(nElems[val_section]+1, 0);
      if (cg_poly_elements_partial_read(cgnsFileID, cgnsBase, cgnsZone,
                                        val_section+1, startE, endE,
                                        connElemTemp.data(),
                                        connOffsetTemp.data(), NULL) != CG_OK)
        cg_error_exit();
    } else {
      if (cg_elements_read(cgnsFileID, cgnsBase, cgnsZone, val_section+1,
                           connElemTemp.data(), NULL))
        cg_error_exit();
    }

    /*--- Allocate the memory for the data structure used to carry
     the connectivity for this section. ---*/

    connElems[val_section].resize(nElems[val_section]*SU2_CONN_SIZE,0);

    unsigned long counterCGNS = 0;
    for (iElem = 0; iElem < nElems[val_section]; iElem++) {

      ElementType_t iElemType = elemType;

      /*--- If we have a mixed element section, we need to check the elem
       type one-by-one. We also must manually advance the counter. ---*/

      if (isMixed) {
        iElemType = ElementType_t(connElemTemp[counterCGNS]);
        counterCGNS++;
      }

      /*--- Get the VTK type for this element. ---*/

      int vtk_type;
      string elem_name = GetCGNSElementType(iElemType, vtk_type);

      /*--- Get the number of nodes per element. ---*/

      cg_npe(iElemType, &npe);

      /*--- Load the surface element connectivity into the SU2 data
       structure with format: [globalID VTK n1 n2 n3 n4 n5 n6 n7 n8].
       We do not need a global ID for the surface elements, so we
       simply set that to zero to maintain the same data structure
       format as the interior elements. Note that we subtract 1 to
       move from the CGNS 1-based indexing to SU2's zero-based. ---*/

      connElems[val_section][iElem*SU2_CONN_SIZE+0] = 0;
      connElems[val_section][iElem*SU2_CONN_SIZE+1] = vtk_type;
      for (iNode = 0; iNode < (unsigned long)npe; iNode++) {
        unsigned long nn = iElem*SU2_CONN_SIZE+SU2_CONN_SKIP+iNode;
        connElems[val_section][nn] = connElemTemp[counterCGNS] - 1;
        counterCGNS++;
      }

    }

  } else {

    /*--- We are not the master, so we resize to zero for safety. ---*/

    nElems[val_section] = 0;
    connElems[val_section].resize(0);

  }

}

void CCGNSMeshReaderFVM::ReformatCGNSVolumeConnectivity() {

  /*--- Loop to store total number of elements we have locally.
   This number includes repeats across ranks due to redistribution
   according to the linear partitioning of the grid nodes. ---*/

  numberOfLocalElements = 0;
  for (int s = 0; s < nSections; s++)
    if (isInterior[s]) numberOfLocalElements += nElems[s];

  /* Put our CGNS data into the class data structures for the mesh reader */

  localVolumeElementConnectivity.resize(numberOfLocalElements*SU2_CONN_SIZE);
  unsigned long count = 0;
  for (int s = 0; s < nSections; s++) {
    if (isInterior[s]) {
      for (unsigned long iElem = 0; iElem < nElems[s]; iElem++) {
        for (unsigned long iNode = 0; iNode < SU2_CONN_SIZE; iNode++) {
          unsigned long nn = iElem*SU2_CONN_SIZE+iNode;
          localVolumeElementConnectivity[count] = (unsigned long)connElems[s][nn];
          count++;
        }
      }
      vector<cgsize_t>().swap(connElems[s]);
    }
  }

}

void CCGNSMeshReaderFVM::ReformatCGNSSurfaceConnectivity() {

  /*--- Prepare the class data for the marker names and connectivity. ---*/

  markerNames.resize(numberOfMarkers);
  surfaceElementConnectivity.resize(numberOfMarkers);

  int markerCount  = 0;
  int elementCount = 0;
  for (int s = 0; s < nSections; s++) {
    if (!isInterior[s]) {

      /*--- Store the tag for this marker. Remove any whitespaces from
       the marker names found in the CGNS file to avoid any issues. ---*/

      string Marker_Tag = string(sectionNames[s].data());
      Marker_Tag.erase(remove(Marker_Tag.begin(), Marker_Tag.end(),' '),
                       Marker_Tag.end());
      markerNames[markerCount] = Marker_Tag;

      /*--- The master node alone stores the connectivity. ---*/

      if (rank == MASTER_NODE) {
        surfaceElementConnectivity[markerCount].resize(nElems[s]*SU2_CONN_SIZE);
        elementCount = 0;
        for (unsigned long iElem = 0; iElem < nElems[s]; iElem++) {
          for (unsigned long iNode = 0; iNode < SU2_CONN_SIZE; iNode++) {
            unsigned long nn = iElem*SU2_CONN_SIZE+iNode;
            surfaceElementConnectivity[markerCount][elementCount] = (unsigned long)connElems[s][nn];
            elementCount++;
          }
        }
        vector<cgsize_t>().swap(connElems[s]);
      }
      markerCount++;
    }
  }

}
#endif

void CCGNSMeshReaderFVM::InitiateCommsAll(void *bufSend,
                                          const int *nElemSend,
                                          SU2_MPI::Request *sendReq,
                                          void *bufRecv,
                                          const int *nElemRecv,
                                          SU2_MPI::Request *recvReq,
                                          unsigned short countPerElem,
                                          unsigned short commType) {

  /*--- Local variables ---*/

  int iMessage, iProc, offset, nElem, count, source, dest, tag;

  /*--- Launch the non-blocking recv's first. ---*/

  iMessage = 0;
  for (iProc = 0; iProc < size; iProc++) {

    /*--- Post recv's only if another proc is sending us data. We do
     not communicate with ourselves or post recv's for zero length
     messages to keep overhead down. ---*/

    if ((nElemRecv[iProc+1] > nElemRecv[iProc]) && (iProc != rank)) {

      /*--- Compute our location in the recv buffer. ---*/

      offset = countPerElem*nElemRecv[iProc];

      /*--- Take advantage of cumulative storage format to get the number
       of elems that we need to recv. ---*/

      nElem = nElemRecv[iProc+1] - nElemRecv[iProc];

      /*--- Total count can include multiple pieces of data per element. ---*/

      count = countPerElem*nElem;

      /*--- Post non-blocking recv for this proc. ---*/

      source = iProc; tag = iProc + 1;

      switch (commType) {
        case COMM_TYPE_DOUBLE:
          SU2_MPI::Irecv(&(static_cast<su2double*>(bufRecv)[offset]),
                         count, MPI_DOUBLE, source, tag, SU2_MPI::GetComm(),
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_LONG:
          SU2_MPI::Irecv(&(static_cast<unsigned long*>(bufRecv)[offset]),
                         count, MPI_UNSIGNED_LONG, source, tag, SU2_MPI::GetComm(),
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_LONG:
          SU2_MPI::Irecv(&(static_cast<long*>(bufRecv)[offset]),
                         count, MPI_LONG, source, tag, SU2_MPI::GetComm(),
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_SHORT:
          SU2_MPI::Irecv(&(static_cast<unsigned short*>(bufRecv)[offset]),
                         count, MPI_UNSIGNED_SHORT, source, tag, SU2_MPI::GetComm(),
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_CHAR:
          SU2_MPI::Irecv(&(static_cast<char*>(bufRecv)[offset]),
                         count, MPI_CHAR, source, tag, SU2_MPI::GetComm(),
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_SHORT:
          SU2_MPI::Irecv(&(static_cast<short*>(bufRecv)[offset]),
                         count, MPI_SHORT, source, tag, SU2_MPI::GetComm(),
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_INT:
          SU2_MPI::Irecv(&(static_cast<int*>(bufRecv)[offset]),
                         count, MPI_INT, source, tag, SU2_MPI::GetComm(),
                         &(recvReq[iMessage]));
          break;
        default:
          break;
      }

      /*--- Increment message counter. ---*/

      iMessage++;

    }
  }

  /*--- Launch the non-blocking sends next. ---*/

  iMessage = 0;
  for (iProc = 0; iProc < size; iProc++) {

    /*--- Post sends only if we are sending another proc data. We do
     not communicate with ourselves or post sends for zero length
     messages to keep overhead down. ---*/

    if ((nElemSend[iProc+1] > nElemSend[iProc]) && (iProc != rank)) {

      /*--- Compute our location in the send buffer. ---*/

      offset = countPerElem*nElemSend[iProc];

      /*--- Take advantage of cumulative storage format to get the number
       of elems that we need to send. ---*/

      nElem = nElemSend[iProc+1] - nElemSend[iProc];

      /*--- Total count can include multiple pieces of data per element. ---*/

      count = countPerElem*nElem;

      /*--- Post non-blocking send for this proc. ---*/

      dest = iProc; tag = rank + 1;

      switch (commType) {
        case COMM_TYPE_DOUBLE:
          SU2_MPI::Isend(&(static_cast<su2double*>(bufSend)[offset]),
                         count, MPI_DOUBLE, dest, tag, SU2_MPI::GetComm(),
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_LONG:
          SU2_MPI::Isend(&(static_cast<unsigned long*>(bufSend)[offset]),
                         count, MPI_UNSIGNED_LONG, dest, tag, SU2_MPI::GetComm(),
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_LONG:
          SU2_MPI::Isend(&(static_cast<long*>(bufSend)[offset]),
                         count, MPI_LONG, dest, tag, SU2_MPI::GetComm(),
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_SHORT:
          SU2_MPI::Isend(&(static_cast<unsigned short*>(bufSend)[offset]),
                         count, MPI_UNSIGNED_SHORT, dest, tag, SU2_MPI::GetComm(),
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_CHAR:
          SU2_MPI::Isend(&(static_cast<char*>(bufSend)[offset]),
                         count, MPI_CHAR, dest, tag, SU2_MPI::GetComm(),
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_SHORT:
          SU2_MPI::Isend(&(static_cast<short*>(bufSend)[offset]),
                         count, MPI_SHORT, dest, tag, SU2_MPI::GetComm(),
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_INT:
          SU2_MPI::Isend(&(static_cast<int*>(bufSend)[offset]),
                         count, MPI_INT, dest, tag, SU2_MPI::GetComm(),
                         &(sendReq[iMessage]));
          break;
        default:
          break;
      }

      /*--- Increment message counter. ---*/

      iMessage++;

    }
  }

}

void CCGNSMeshReaderFVM::CompleteCommsAll(int nSends,
                                          SU2_MPI::Request *sendReq,
                                          int nRecvs,
                                          SU2_MPI::Request *recvReq) {

  /*--- Local variables ---*/

  int ind, iSend, iRecv;
  SU2_MPI::Status status;

  /*--- Wait for the non-blocking sends to complete. ---*/

  for (iSend = 0; iSend < nSends; iSend++)
    SU2_MPI::Waitany(nSends, sendReq, &ind, &status);

  /*--- Wait for the non-blocking recvs to complete. ---*/

  for (iRecv = 0; iRecv < nRecvs; iRecv++)
    SU2_MPI::Waitany(nRecvs, recvReq, &ind, &status);

}
