/*!
 * \file CCGNSMeshReaderFEM.cpp
 * \brief Class that reads a single zone of a CGNS mesh file from disk into
 *        linear partitions across all ranks.
 * \author T. Economon
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../include/geometry/meshreader/CCGNSMeshReaderFEM.hpp"
#include "../../../include/geometry/meshreader/CCGNSElementType.hpp"

CCGNSMeshReaderFEM::CCGNSMeshReaderFEM(const CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone)
    : CCGNSMeshReaderBase(val_config, val_iZone, val_nZone) {
#ifdef HAVE_CGNS
  OpenCGNSFile(config->GetMesh_FileName());

  /*--- Read the basic information about the database and zone(s). ---*/
  ReadCGNSDatabaseMetadata();
  ReadCGNSZoneMetadata();

  /*--- Read the basic information about the sections. ---*/
  ReadCGNSSectionMetadata();

  /*--- Read the volume connectivity and distribute it
        linearly over the MPI ranks. ---*/
  ReadCGNSVolumeElementConnectivity();

  /*--- Read the coordinates of the points and communicate the ones that
        are needed on this MPI rank. ---*/
  ReadCGNSPointCoordinates();
  CommPointCoordinates();

  /*--- Read the surface connectivity and store the surface elements whose
        corresponding volume element is stored on this MPI rank. ---*/
  ReadCGNSSurfaceElementConnectivity();

  /*--- We have extracted all CGNS data. Close the CGNS file. ---*/
  if (cg_close(cgnsFileID)) cg_error_exit();

#else
  SU2_MPI::Error(string(" SU2 built without CGNS support. \n") + string(" To use CGNS, build SU2 accordingly."),
                 CURRENT_FUNCTION);
#endif
}

CCGNSMeshReaderFEM::~CCGNSMeshReaderFEM() = default;

#ifdef HAVE_CGNS
void CCGNSMeshReaderFEM::ReadCGNSConnectivityRangeSection(const int val_section, const unsigned long val_firstIndex,
                                                          const unsigned long val_lastIndex, unsigned long& elemCount,
                                                          unsigned long& localElemCount,
                                                          vector<unsigned long>& localConn) {
  /*--- Read the connectivity details for this section. ---*/
  int nbndry, parent_flag;
  cgsize_t startE, endE;
  ElementType_t elemType;
  char sectionName[CGNS_STRING_SIZE];

  if (cg_section_read(cgnsFileID, cgnsBase, cgnsZone, val_section + 1, sectionName, &elemType, &startE, &endE, &nbndry,
                      &parent_flag))
    cg_error_exit();

  /*--- Determine the number of elements in this section and update
        the element counters accordingly. ---*/
  const unsigned long nElemSection = (endE - startE + 1);
  const unsigned long elemCountOld = elemCount;
  elemCount += nElemSection;

  /*--- Check for overlap with the element range this rank is responsible for. ---*/
  const unsigned long indBegOverlap = max(elemCountOld, val_firstIndex);
  const unsigned long indEndOverlap = min(elemCount, val_lastIndex);

  if (indEndOverlap > indBegOverlap) {
    /*--- This rank must read element data from this connectivity section.
          Determine the offset relative to the start of this section and
          the number of elements to be read by this rank. ---*/
    const unsigned long offsetRank = indBegOverlap - elemCountOld;
    const unsigned long nElemRank = indEndOverlap - indBegOverlap;
    nElems[val_section] = nElemRank;

    /*--- Determine the index range to be read for this rank. ---*/
    const cgsize_t iBeg = startE + offsetRank;
    const cgsize_t iEnd = iBeg + nElemRank - 1;

    /*--- Determine the size of the vector needed to read
          the connectivity data from the CGNS file. ---*/
    cgsize_t sizeNeeded;
    if (cg_ElementPartialSize(cgnsFileID, cgnsBase, cgnsZone, val_section + 1, iBeg, iEnd, &sizeNeeded) != CG_OK)
      cg_error_exit();

    /*--- Allocate the memory for the connectivity and read the data. ---*/
    vector<cgsize_t> connCGNSVec(sizeNeeded);
    if (elemType == MIXED) {
      vector<cgsize_t> connCGNSOffsetVec(iEnd - iBeg + 2);
      if (cg_poly_elements_partial_read(cgnsFileID, cgnsBase, cgnsZone, val_section + 1, iBeg, iEnd, connCGNSVec.data(),
                                        connCGNSOffsetVec.data(), NULL) != CG_OK)
        cg_error_exit();

    } else {
      if (cg_elements_partial_read(cgnsFileID, cgnsBase, cgnsZone, val_section + 1, iBeg, iEnd, connCGNSVec.data(),
                                   NULL) != CG_OK)
        cg_error_exit();
    }

    /*--- Define the variables needed for the conversion of the CGNS
          format to the SU2 internal format. ---*/
    CCGNSElementType CGNSElem;
    std::vector<unsigned long> connSU2;
    ElementType_t typeElem = elemType;
    const cgsize_t* connCGNS = connCGNSVec.data();

    /*--- Loop over the elements just read. ---*/
    for (unsigned long i = 0; i < nElemRank; ++i, ++localElemCount) {
      /*--- Determine the element type for this element if this is a mixed
            connectivity and set the pointer to the actual connectivity data. ---*/
      if (elemType == MIXED) {
        typeElem = (ElementType_t)connCGNS[0];
        ++connCGNS;
      }

      /*--- Convert the CGNS connectivity to SU2 connectivity. ---*/
      const unsigned long globalID = val_firstIndex + localElemCount;
      CGNSElem.CGNSToSU2(typeElem, globalID, connCGNS, connSU2);

      /*--- Update the pointer for the connectivity for the next
            element and store the SU2 connectivity in localConn. ---*/
      connCGNS += connSU2[3];

      localConn.insert(localConn.end(), connSU2.begin(), connSU2.end());
    }
  }
}

void CCGNSMeshReaderFEM::ReadCGNSVolumeElementConnectivity(void) {
  /*--- Write a message that the volume elements are loaded. ---*/
  if (rank == MASTER_NODE) {
    if (size > SINGLE_NODE)
      cout << "Loading volume elements into linear partitions." << endl;
    else
      cout << "Loading volume elements." << endl;
  }

  /*--- Get a partitioner to help with linear partitioning. ---*/
  CLinearPartitioner elemPartitioner(numberOfGlobalElements, 0);

  /*--- Determine the index of the first and last element to be stored
        on this rank and the number of local elements. ---*/
  const unsigned long firstIndex = elemPartitioner.GetFirstIndexOnRank(rank);
  const unsigned long lastIndex = elemPartitioner.GetLastIndexOnRank(rank);
  numberOfLocalElements = elemPartitioner.GetSizeOnRank(rank);

  /*--- Loop over the section and check for a section with volume elements. ---*/
  unsigned long elemCount = 0, localElemCount = 0;
  for (int s = 0; s < nSections; ++s) {
    if (isInterior[s]) {
      /*--- Read the connectivity of this section and store the
            data in localVolumeElementConnectivity. ---*/
      ReadCGNSConnectivityRangeSection(s, firstIndex, lastIndex, elemCount, localElemCount,
                                       localVolumeElementConnectivity);
    }
  }
}

void CCGNSMeshReaderFEM::ReadCGNSSurfaceElementConnectivity(void) {
  /*--- Write a message that the surface elements are loaded. ---*/
  if (rank == MASTER_NODE) cout << "Loading surface elements." << endl;

  /*--- Determine the vector to hold the faces of the local elements. ---*/
  vector<CFaceOfElement> localFaces;
  DetermineFacesVolumeElements(localFaces);

  /*--- Determine the number of markers. ---*/
  numberOfMarkers = 0;
  for (int s = 0; s < nSections; s++)
    if (!isInterior[s]) ++numberOfMarkers;

  /*--- Allocate the memory for the number of markers and local surface elements.
        Also allocate the first index of surfaceElementConnectivity. ---*/
  markerNames.resize(numberOfMarkers);
  numberOfLocalSurfaceElements.resize(numberOfMarkers);
  surfaceElementConnectivity.resize(numberOfMarkers);

  /*--- Loop over the number of sections and check for a surface. ---*/
  int markerCount = 0;
  for (int s = 0; s < nSections; ++s) {
    if (!isInterior[s]) {
      /*--- Create the marker name. ---*/
      string Marker_Tag = string(sectionNames[s].data());
      Marker_Tag.erase(remove(Marker_Tag.begin(), Marker_Tag.end(), ' '), Marker_Tag.end());
      markerNames[markerCount] = Marker_Tag;

      /*--- Call the function ReadCGNSSurfaceSection to carry out
            the actual reading and storing of the required faces. ---*/
      ReadCGNSSurfaceSection(s, localFaces, numberOfLocalSurfaceElements[markerCount],
                             surfaceElementConnectivity[markerCount]);

      /*--- Update the marker counter. ---*/
      ++markerCount;
    }
  }
}

void CCGNSMeshReaderFEM::ReadCGNSSurfaceSection(const int val_section, const vector<CFaceOfElement>& localFaces,
                                                unsigned long& nSurfElem, vector<unsigned long>& surfConn) {
  /*--- Initialize nSurfElem to zero. ---*/
  nSurfElem = 0;

  /*--- Read the connectivity details for this section and determine the number
        of elements present in this section. ---*/
  int nbndry, parent_flag;
  cgsize_t startE, endE;
  ElementType_t elemType;
  char sectionName[CGNS_STRING_SIZE];

  if (cg_section_read(cgnsFileID, cgnsBase, cgnsZone, val_section + 1, sectionName, &elemType, &startE, &endE, &nbndry,
                      &parent_flag))
    cg_error_exit();

  const unsigned long nElemSection = (endE - startE + 1);

  /*--- Determine the number of chunks used for the reading of the surface
        elements. This is done to avoid a memory bottleneck for extremely
        big cases. For reasonably sized grids this connectivity can be
        read in a single call. ---*/
  unsigned long nChunks = nElemSection / localFaces.size();
  if (nChunks * localFaces.size() != nElemSection) ++nChunks;
  const unsigned long nElemChunk = nElemSection / nChunks;

  /*--- Loop over the number of chunks. ---*/
  for (unsigned long iChunk = 0; iChunk < nChunks; ++iChunk) {
    /*--- Determine the start and end index for this chunk. ---*/
    const cgsize_t iBeg = startE + iChunk * nElemChunk;
    const cgsize_t iEnd = (iChunk == (nChunks - 1)) ? endE : iBeg + nElemChunk - 1;

    /*--- Determine the size of the vector needed to read
          the connectivity data from the CGNS file. ---*/
    cgsize_t sizeNeeded;
    if (cg_ElementPartialSize(cgnsFileID, cgnsBase, cgnsZone, val_section + 1, iBeg, iEnd, &sizeNeeded) != CG_OK)
      cg_error_exit();

    /*--- Allocate the memory for the connectivity and read the data. ---*/
    vector<cgsize_t> connCGNSVec(sizeNeeded);
    if (elemType == MIXED) {
      vector<cgsize_t> connCGNSOffsetVec(iEnd - iBeg + 2);
      if (cg_poly_elements_partial_read(cgnsFileID, cgnsBase, cgnsZone, val_section + 1, iBeg, iEnd, connCGNSVec.data(),
                                        connCGNSOffsetVec.data(), NULL) != CG_OK)
        cg_error_exit();

    } else {
      if (cg_elements_partial_read(cgnsFileID, cgnsBase, cgnsZone, val_section + 1, iBeg, iEnd, connCGNSVec.data(),
                                   NULL) != CG_OK)
        cg_error_exit();
    }

    /*--- Define the variables needed for the conversion of the CGNS
          format to the SU2 internal format. ---*/
    CCGNSElementType CGNSElem;
    std::vector<unsigned long> connSU2;
    ElementType_t typeElem = elemType;
    const cgsize_t* connCGNS = connCGNSVec.data();

    /*--- Loop over the elements just read. ---*/
    for (cgsize_t i = iBeg; i <= iEnd; ++i) {
      /*--- Determine the element type for this element if this is a mixed
            connectivity and set the pointer to the actual connectivity data. ---*/
      if (elemType == MIXED) {
        typeElem = (ElementType_t)connCGNS[0];
        ++connCGNS;
      }

      /*--- Convert the CGNS connectivity to SU2 connectivity and update
            the pointer for the next surface element. ---*/
      const unsigned long globalID = i - 1;
      CGNSElem.CGNSToSU2(typeElem, globalID, connCGNS, connSU2);
      connCGNS += connSU2[3];

      /*--- Easier storage of the VTK type, polynomial degree
            and number of DOFs of the surface element. ---*/
      const unsigned short VTK_Type = static_cast<unsigned short>(connSU2[0]);
      const unsigned short nPolyGrid = static_cast<unsigned short>(connSU2[1]);
      const unsigned short nDOFsGrid = static_cast<unsigned short>(connSU2[3]);

      /*--- Make a distinction between the possible element surface types and
            determine the corner points in local numbering of the element. ---*/
      const unsigned short nDOFEdgeGrid = nPolyGrid + 1;

      CFaceOfElement thisFace;
      thisFace.cornerPoints[0] = 0;
      thisFace.cornerPoints[1] = nPolyGrid;

      switch (VTK_Type) {
        case LINE:
          thisFace.nCornerPoints = 2;
          break;

        case TRIANGLE:
          thisFace.nCornerPoints = 3;
          thisFace.cornerPoints[2] = nDOFsGrid - 1;
          break;

        case QUADRILATERAL:
          thisFace.nCornerPoints = 4;
          thisFace.cornerPoints[2] = static_cast<unsigned long>(nPolyGrid) * nDOFEdgeGrid;
          thisFace.cornerPoints[3] = nDOFsGrid - 1;
          break;

        default:
          ostringstream message;
          message << "Unsupported FEM boundary element value, " << typeElem << ", in surface section " << sectionName;
          SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /*--- Convert the local numbering of thisFace to global numbering
            and create a unique numbering of corner points. ---*/
      for (unsigned short j = 0; j < thisFace.nCornerPoints; ++j)
        thisFace.cornerPoints[j] = connSU2[thisFace.cornerPoints[j] + 5];
      thisFace.CreateUniqueNumbering();

      /*--- Check if this boundary face must be stored on this rank. ---*/
      vector<CFaceOfElement>::const_iterator low;
      low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);
      if (low != localFaces.end()) {
        if (!(thisFace < *low)) {
          /*--- Update the number of local surface elements. ---*/
          ++nSurfElem;

          /*--- Store the meta data in the first 5 locations of connSU2. ---*/
          connSU2[0] = VTK_Type;
          connSU2[1] = nPolyGrid;
          connSU2[2] = nDOFsGrid;
          connSU2[3] = globalID;      // Global surface elem ID.
          connSU2[4] = low->elemID0;  // Global volume elem ID.

          /*--- Store the connectivity data in surfConn. ---*/
          surfConn.insert(surfConn.end(), connSU2.begin(), connSU2.end());
        }
      }
    }
  }
}
#endif

void CCGNSMeshReaderFEM::CommPointCoordinates(void) {
  /*--- Determine the linear partitioning of the points. ---*/
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints, 0);

  /*--- Loop over the local elements to determine the global
        point IDs to be stored on this rank. --*/
  unsigned long ind = 0;
  for (unsigned long i = 0; i < numberOfLocalElements; ++i) {
    /*--- Store the number of grid DOFs for this element and
          skip the meta data for this element (5 entries). ---*/
    const unsigned long nDOFsGrid = localVolumeElementConnectivity[ind + 3];
    ind += 5;

    /*--- Copy the connectivity to globalPointIDs. ---*/
    unsigned long* conn = localVolumeElementConnectivity.data() + ind;
    ind += nDOFsGrid;
    globalPointIDs.insert(globalPointIDs.end(), conn, conn + nDOFsGrid);
  }

  /*--- Sort the globalPointIDs and remove the duplicate entries. ---*/
  sort(globalPointIDs.begin(), globalPointIDs.end());
  vector<unsigned long>::iterator lastNode;
  lastNode = unique(globalPointIDs.begin(), globalPointIDs.end());
  globalPointIDs.erase(lastNode, globalPointIDs.end());

  /*--- Determine the number of locally stored points. ---*/
  numberOfLocalPoints = globalPointIDs.size();

  /*--- This rest of this function only needs to done when MPI is used. ---*/
#ifdef HAVE_MPI

  /*--- Store the global ID's of the points in such a way that they can
        be sent to the rank that actually stores the coordinates.. ---*/
  vector<vector<unsigned long> > pointBuf(size, vector<unsigned long>(0));

  for (unsigned long i = 0; i < globalPointIDs.size(); ++i)
    pointBuf[pointPartitioner.GetRankContainingIndex(globalPointIDs[i])].push_back(globalPointIDs[i]);

  /*--- Determine the total number of ranks to which this rank will send
        a message and also determine the number of ranks from which this
        rank will receive a message. Furthermore, determine the starting
        indices where data from the different ranks should be stored in
        localPointCoordinates. ---*/
  int nRankSend = 0;
  vector<int> sendToRank(size, 0);
  vector<unsigned long> startingIndRanksInPoint(size + 1);
  startingIndRanksInPoint[0] = 0;

  for (int i = 0; i < size; ++i) {
    startingIndRanksInPoint[i + 1] = startingIndRanksInPoint[i] + pointBuf[i].size();

    if (pointBuf[i].size()) {
      ++nRankSend;
      sendToRank[i] = 1;
    }
  }

  int nRankRecv;
  vector<int> sizeRecv(size, 1);
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /*--- Send out the messages with the global point numbers. Use nonblocking
        sends to avoid deadlock. ---*/
  vector<SU2_MPI::Request> sendReqs(nRankSend);
  nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (pointBuf[i].size()) {
      SU2_MPI::Isend(pointBuf[i].data(), pointBuf[i].size(), MPI_UNSIGNED_LONG, i, i, SU2_MPI::GetComm(),
                     &sendReqs[nRankSend]);
      ++nRankSend;
    }
  }

  /*--- Define the communication buffer for the coordinates and the vector
        for the return communication requests. */
  vector<SU2_MPI::Request> returnReqs(nRankRecv);
  vector<vector<su2double> > coorReturnBuf(nRankRecv, vector<su2double>(0));

  /*--- Get the first index of the points as well as the number of points
        currently stored on this rank. ---*/
  const unsigned long firstIndex = pointPartitioner.GetFirstIndexOnRank(rank);
  const unsigned long nPointsRead = pointPartitioner.GetSizeOnRank(rank);

  /*--- Loop over the number of ranks from which this rank receives global
        point numbers that should be stored on this rank. ---*/
  for (int i = 0; i < nRankRecv; ++i) {
    /* Block until a message arrives. Determine the source and size
       of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /*--- Allocate the memory for a buffer to receive this message and also
          for the buffer to return to coordinates. ---*/
    vector<unsigned long> pointRecvBuf(sizeMess);
    coorReturnBuf[i].resize(static_cast<size_t>(dimension) * sizeMess);

    /* Receive the message using a blocking receive. */
    SU2_MPI::Recv(pointRecvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, source, rank, SU2_MPI::GetComm(), &status);

    /*--- Loop over the nodes just received and fill the return communication
          buffer with the coordinates of the requested nodes. ---*/
    for (int j = 0; j < sizeMess; ++j) {
      const int jj = dimension * j;
      const long kk = pointRecvBuf[j] - firstIndex;
      if (kk < 0 || kk >= static_cast<long>(nPointsRead))
        SU2_MPI::Error("Invalid point requested. This should not happen.", CURRENT_FUNCTION);

      for (int k = 0; k < dimension; ++k) coorReturnBuf[i][jj + k] = localPointCoordinates[k][kk];
    }

    /* Send the buffer just filled back to the requesting rank.
       Use a non-blocking send to avoid deadlock. */
    SU2_MPI::Isend(coorReturnBuf[i].data(), coorReturnBuf[i].size(), MPI_DOUBLE, source, source + 1, SU2_MPI::GetComm(),
                   &returnReqs[i]);
  }

  /*--- Resize the second indices of localPointCoordinates. ---*/
  for (int k = 0; k < dimension; ++k) localPointCoordinates[k].resize(numberOfLocalPoints);

  /*--- Loop over the ranks from which this rank has requested coordinates. ---*/
  for (int i = 0; i < nRankSend; ++i) {
    /* Block until a message arrives. Determine the source of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank + 1, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    /* Allocate the memory for the coordinate receive buffer. */
    vector<su2double> coorRecvBuf(dimension * pointBuf[source].size());

    /* Receive the message using a blocking receive. */
    SU2_MPI::Recv(coorRecvBuf.data(), coorRecvBuf.size(), MPI_DOUBLE, source, rank + 1, SU2_MPI::GetComm(), &status);

    /*--- Loop over the points just received. ---*/
    for (unsigned long j = 0; j < pointBuf[source].size(); ++j) {
      /*--- Copy the data into the correct location of localPointCoordinates. ---*/
      const unsigned long jj = dimension * j;
      const unsigned long kk = startingIndRanksInPoint[source] + j;

      for (int k = 0; k < dimension; ++k) localPointCoordinates[k][kk] = SU2_TYPE::GetValue(coorRecvBuf[jj + k]);
    }
  }

  /* Complete the non-blocking sends of both rounds. */
  SU2_MPI::Waitall(sendReqs.size(), sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Waitall(returnReqs.size(), returnReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Wild cards have been used in the communication,
        so synchronize the ranks to avoid problems. ---*/
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif
}
