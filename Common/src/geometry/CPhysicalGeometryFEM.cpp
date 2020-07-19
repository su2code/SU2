/*!
 * \file CPhysicalGeometry.cpp
 * \brief Implementation of the FEM part physical geometry class.
 * \author F. Palacios, T. Economon
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../include/toolboxes/fem/CMatchingFace.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridFEM.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"
#include "../../include/geometry/CPhysicalGeometry.hpp"

void CPhysicalGeometry::LoadLinearlyPartitionedPointsFEM(CConfig *config, CMeshReader *mesh) {

  /*--- Get the partitioned coordinates and their global IDs from the mesh object. ---*/
  const auto &gridCoords     = mesh->GetLocalPointCoordinates();
  const auto &globalPointIDs = mesh->GetGlobalPointIDs();

  /*--- Initialize point counts and the grid node data structure. ---*/
  nPointNode = nPoint;
  nodes = new CPoint(nPoint, nDim);

  /*--- Loop over the points and set the coordinates and global index. ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iDim = 0; iDim < nDim; ++iDim)
      nodes->SetCoord(iPoint, iDim, gridCoords[iDim][iPoint]);
    nodes->SetGlobalIndex(iPoint, globalPointIDs[iPoint]);
  }
}

void CPhysicalGeometry::LoadLinearlyPartitionedVolumeElementsFEM(CConfig *config, CMeshReader *mesh) {

  /*--- Reset the global to local element mapping. ---*/
  Global_to_Local_Elem.clear();

  /*--- Get the volume connectivity from the mesh object. ---*/
  const auto &dataElems = mesh->GetLocalVolumeElementConnectivity();

  /*--- Allocate space for the interior elements in our SU2 data
        structure. Note that we only instantiate our rank's local set. ---*/
  elem = new CPrimalGrid*[nElem] ();

  /*--- Loop over all of the internal, local volumetric elements. ---*/
  unsigned long ind = 0;
  for (unsigned long jElem=0; jElem<nElem; ++jElem) {

    /*--- Create a FEM element from the data dataElems. ---*/
    const auto *dataElem = dataElems.data() + ind;
    elem[jElem] = new CPrimalGridFEM(dataElem);

    /*--- Store the global to local mapping in Global_to_Local_Elem. ---*/
    Global_to_Local_Elem[dataElem[4]] = jElem;

    /*--- Update ind for the next element. ---*/
    ind += dataElem[3] + 5;
  }
}

void CPhysicalGeometry::LoadLinearlyPartitionedSurfaceElementsFEM(CConfig *config, CMeshReader *mesh) {

  /*--- Store the number of markers and print to the screen. ---*/
  nMarker = mesh->GetNumberOfMarkers();
  config->SetnMarker_All(nMarker);
  if (rank == MASTER_NODE)
    cout << nMarker << " surface markers." << endl;

  /*--- Create the data structure for boundary elements. ---*/
  bound         = new CPrimalGrid**[nMarker];
  nElem_Bound   = new unsigned long [nMarker];
  Tag_to_Marker = new string [config->GetnMarker_Max()];

  /*--- Retrieve the name of the surface markers as well as
        the number of surface elements for every marker. ---*/
  const auto &sectionNames       = mesh->GetMarkerNames();
  const auto &nSurfElemPerMarker = mesh->GetNumberOfSurfaceElementsAllMarkers();

  /*--- Loop over all sections that we extracted from the grid file
        that were identified as boundary element sections so that we can
        store those elements into our SU2 data structures. ---*/
  for (int iMarker = 0; iMarker < nMarker; ++iMarker) {

    /*--- Get the string name and set the number of surface elements
          for this marker. ---*/
    string Marker_Tag    = sectionNames[iMarker];
    nElem_Bound[iMarker] = nSurfElemPerMarker[iMarker];

    /*--- Allocate the memory of the pointers for the surface
          elements for this marker. ---*/
    bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];

    /*--- Retrieve the boundary element data for this marker. ---*/
    const auto &dataElems = mesh->GetSurfaceElementConnectivityForMarker(iMarker);

    /*--- Loop over the number of boundary elements for this marker. ---*/
    unsigned long ind = 0;
    for (unsigned long jElem=0; jElem<nElem_Bound[iMarker]; ++jElem) {

      /*--- Create a boundary FEM element from the data dataElems. ---*/
      const auto *dataElem = dataElems.data() + ind;
      bound[iMarker][jElem] = new CPrimalGridBoundFEM(dataElem);

      /*--- Update ind for the next element. ---*/
      ind += dataElem[2] + 5;
    }

    /*--- Update config file lists in order to store the boundary
          information for this marker in the correct place. ---*/
    Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
    config->SetMarker_All_TagBound(iMarker, Marker_Tag);
    config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
    config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
    config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
    config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
    config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
    config->SetMarker_All_Analyze(iMarker, config->GetMarker_CfgFile_Analyze(Marker_Tag));
    config->SetMarker_All_ZoneInterface(iMarker, config->GetMarker_CfgFile_ZoneInterface(Marker_Tag));
    config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
    config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
    config->SetMarker_All_Deform_Mesh(iMarker, config->GetMarker_CfgFile_Deform_Mesh(Marker_Tag));
    config->SetMarker_All_Fluid_Load(iMarker, config->GetMarker_CfgFile_Fluid_Load(Marker_Tag));
    config->SetMarker_All_PyCustom(iMarker, config->GetMarker_CfgFile_PyCustom(Marker_Tag));
    config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
    config->SetMarker_All_SendRecv(iMarker, NONE);
    config->SetMarker_All_Turbomachinery(iMarker, config->GetMarker_CfgFile_Turbomachinery(Marker_Tag));
    config->SetMarker_All_TurbomachineryFlag(iMarker, config->GetMarker_CfgFile_TurbomachineryFlag(Marker_Tag));
    config->SetMarker_All_MixingPlaneInterface(iMarker, config->GetMarker_CfgFile_MixingPlaneInterface(Marker_Tag));
  }
}

void CPhysicalGeometry::SetColorFEMGrid_Parallel(CConfig *config) {

  /*--- Initialize the color of the elements. ---*/
  for(unsigned long i=0; i<nElem; ++i)
    elem[i]->SetColor(0);

  /*--- Determine the matching faces of the local elements. ---*/
  vector<CFaceOfElement> localFaces;
  DetermineMatchingFacesFEMGrid(config, localFaces);

  /*--- In the function above the periodic boundaries are not found.
        A different treatment must be used in order to find these. ---*/
  DeterminePeriodicFacesFEMGrid(config, localFaces);

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CPhysicalGeometry::DetermineMatchingFacesFEMGrid(const CConfig          *config,
                                                      vector<CFaceOfElement> &localFaces) {

  /*--- Loop over all elements to determine the faces of the elements. ---*/
  for(unsigned long k=0; k<nElem; ++k) {

    /*--- Get the global IDs of the corner points of all the faces of this element. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long  faceConn[6][4];

    elem[k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Loop over the faces of the element and add them to localFaces. ---*/
    for(unsigned short i=0; i<nFaces; ++i) {

      /*--- Set the information for this face. ---*/
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for(unsigned short j=0; j<nPointsPerFace[i]; ++j)
        thisFace.cornerPoints[j] = faceConn[i][j];
      thisFace.elemID0    = elem[k]->GetGlobalElemID();
      thisFace.nPolySol0  = elem[k]->GetNPolySol();
      thisFace.nDOFsElem0 = elem[k]->GetNDOFsSol();
      thisFace.elemType0  = elem[k]->GetVTK_Type();

      /*--- Create a unique number for the face and add it to localFaces. ---*/
      thisFace.CreateUniqueNumbering();
      localFaces.push_back(thisFace);
    }
  }

  /*--- Sort localFaces in increasing order. */
  sort(localFaces.begin(), localFaces.end());

  /*--- Loop over the boundary markers and its faces in order to flag the
        physical boundary faces in localFaces. Note that use is made of the
        overloaded function GetCornerPointsAllFaces, which explains the
        dimensions of the variables used in the function call. Also note that
        the periodic boundaries are excluded, because they are not physical.  ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {
      for(unsigned long k=0; k<nElem_Bound[iMarker]; ++k) {

        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long  faceConn[6][4];
        bound[iMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for(unsigned short j=0; j<nPointsPerFace[0]; ++j)
          thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        /*--- Store the corresponding element ID for the boundary element and
              invalidate the face by setting the element ID to an invalid value. ---*/
        bound[iMarker][k]->SetDomainElement(low->elemID0);
        low->elemID0 = Global_nElem + 10;
      }
    }
  }

  /*--- Loop over the faces and check for double entries. If a double entry is
        found, the elemID from the second entry is copied to the first entry,
        the polynomial degree is updated, and the second entry is invalidated. ---*/
  unsigned long nFacesLoc = localFaces.size();
  for(unsigned long i=1; i<nFacesLoc; ++i) {
    if(localFaces[i] == localFaces[i-1]) {
      localFaces[i-1].elemID1    = localFaces[i].elemID0;
      localFaces[i-1].nPolySol1  = localFaces[i].nPolySol0;
      localFaces[i-1].nDOFsElem1 = localFaces[i].nDOFsElem0;
      localFaces[i-1].elemType1  = localFaces[i].elemType0;
      localFaces[i].elemID0      = Global_nElem + 10;
    }
  }

  /*--- Remove the invalidated faces. This is accomplished by giving the
        face four points a global node ID that is larger than the largest
        point ID in the grid. In this way the sorting operator puts
        these faces at the end of the vector, see also the < operator
        of CFaceOfElement.                                         ---*/
  unsigned long nFacesLocOr = nFacesLoc;
  for(unsigned long i=0; i<nFacesLocOr; ++i) {
    if(localFaces[i].elemID0 > Global_nElem) {
      localFaces[i].nCornerPoints = 4;
      localFaces[i].cornerPoints[0] = Global_nPoint;
      localFaces[i].cornerPoints[1] = Global_nPoint;
      localFaces[i].cornerPoints[2] = Global_nPoint;
      localFaces[i].cornerPoints[3] = Global_nPoint;
      --nFacesLoc;
    }
  }

  sort(localFaces.begin(), localFaces.end());
  localFaces.resize(nFacesLoc);

  /*--- The remaining part of this function is only needed if MPI is used. ---*/ 
#ifdef HAVE_MPI

  /*--- Determine and store the faces, for which only one neighboring element
        was found. For these faces the other neighbor might be stored on a
        different rank (unless there are periodic or non-matching interfaces). ---*/  
  vector<CFaceOfElement> localFacesComm;
  for(unsigned long i=0; i<nFacesLoc; ++i)
    if(localFaces[i].elemID1 > Global_nElem) localFacesComm.push_back(localFaces[i]);

  /*--- Determine the maximum global point ID that occurs in localFacesComm
        of all ranks. Note that only the first point is taken into account,
        because this point determines the rank where the face is sent to. ---*/
  unsigned long nFacesLocComm = localFacesComm.size();
  unsigned long maxPointIDLoc = 0;
  for(unsigned long i=0; i<nFacesLocComm; ++i)
    maxPointIDLoc = max(maxPointIDLoc, localFacesComm[i].cornerPoints[0]);

  unsigned long maxPointID;
  SU2_MPI::Allreduce(&maxPointIDLoc, &maxPointID, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  ++maxPointID;

  /*--- Create a linear partition for the points that occur in
        the faces of localFacesComm. ---*/
  CLinearPartitioner facePartitioner(maxPointID,0);  

  /*--- Define the send buffers for the faces. ---*/
  vector<vector<unsigned long> > sendBufFace(size, vector<unsigned long>(0));

  /*--- Loop over the local faces to be communicated and fill
        the entries of sendBufFace. ---*/
  for(unsigned long i=0; i<nFacesLocComm; ++i) {
    const unsigned long rankFace = facePartitioner.GetRankContainingIndex(localFacesComm[i].cornerPoints[0]);

    sendBufFace[rankFace].push_back(localFacesComm[i].nCornerPoints);
    sendBufFace[rankFace].push_back(localFacesComm[i].cornerPoints[0]);
    sendBufFace[rankFace].push_back(localFacesComm[i].cornerPoints[1]);
    sendBufFace[rankFace].push_back(localFacesComm[i].cornerPoints[2]);
    sendBufFace[rankFace].push_back(localFacesComm[i].cornerPoints[3]);
    sendBufFace[rankFace].push_back(localFacesComm[i].elemID0);
    sendBufFace[rankFace].push_back(localFacesComm[i].nPolySol0);
    sendBufFace[rankFace].push_back(localFacesComm[i].nDOFsElem0);
    sendBufFace[rankFace].push_back(localFacesComm[i].elemType0);
  }

  /*--- Determine the number of ranks from which I receive a message. */
  vector<unsigned long> counter(size, 0);
  unsigned long nMessSend = 0;
  for(int i=0; i<size; ++i) {
    if( sendBufFace[i].size() ) {counter[i] = 1; ++nMessSend;}
  }

  vector<int> sizeRecv(size, 1);

  unsigned long nMessRecv;
  SU2_MPI::Reduce_scatter(counter.data(), &nMessRecv, sizeRecv.data(),
                          MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  /*--- Send the data using nonblocking sends. ---*/
  vector<SU2_MPI::Request> commReqs(max(nMessSend,nMessRecv));

  nMessSend = 0;
  for(int i=0; i<size; ++i) {
    if( sendBufFace[i].size() ) {
      SU2_MPI::Isend(sendBufFace[i].data(), sendBufFace[i].size(),
                     MPI_UNSIGNED_LONG, i, i, MPI_COMM_WORLD,
                     &commReqs[nMessSend]);
      ++nMessSend;
    }
  }

  /*--- Loop over the number of ranks from which faces are received.
        Receive the messages and store the faces in facesRecv. ---*/
  vector<vector<CFaceOfElement> > facesRecv(nMessRecv, vector<CFaceOfElement>(0));
  vector<int> rankRecv(nMessRecv);

  for(unsigned long i=0; i<nMessRecv; ++i) {
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
    rankRecv[i] = status.MPI_SOURCE;
    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  rankRecv[i], rank, MPI_COMM_WORLD, &status);

    const int nFacesRecv = sizeMess/9;
    facesRecv[i].resize(nFacesRecv);

    unsigned long ii = 0;
    for(int j=0; j<nFacesRecv; ++j, ii+=9) {
      facesRecv[i][j].nCornerPoints   = recvBuf[ii];
      facesRecv[i][j].cornerPoints[0] = recvBuf[ii+1];
      facesRecv[i][j].cornerPoints[1] = recvBuf[ii+2];
      facesRecv[i][j].cornerPoints[2] = recvBuf[ii+3];
      facesRecv[i][j].cornerPoints[3] = recvBuf[ii+4];
      facesRecv[i][j].elemID0         = recvBuf[ii+5];
      facesRecv[i][j].nPolySol0       = recvBuf[ii+6];
      facesRecv[i][j].nDOFsElem0      = recvBuf[ii+7];
      facesRecv[i][j].elemType0       = recvBuf[ii+8];
    }
  }

  /*--- The exact contents of facesRecv is needed when communicating back
        the data to the original ranks, so the data is copied to
        localFacesComm, which is not needed anymore. ---*/
  localFacesComm.clear();
  for(unsigned long i=0; i<nMessRecv; ++i)
    localFacesComm.insert(localFacesComm.end(), facesRecv[i].begin(), facesRecv[i].end());

  /*--- Sort localFacesComm in increasing order and determine the
        double entries. ---*/
  sort(localFacesComm.begin(), localFacesComm.end());

  nFacesLocComm = localFacesComm.size();
  for(unsigned long i=1; i<localFacesComm.size(); ++i) {
    if(localFacesComm[i] == localFacesComm[i-1]) {
      localFacesComm[i-1].elemID1    = localFacesComm[i].elemID0;
      localFacesComm[i-1].nPolySol1  = localFacesComm[i].nPolySol0;
      localFacesComm[i-1].nDOFsElem1 = localFacesComm[i].nDOFsElem0;
      localFacesComm[i-1].elemType1  = localFacesComm[i].elemType0;

      localFacesComm[i].nCornerPoints = 4;
      localFacesComm[i].cornerPoints[0] = Global_nPoint;
      localFacesComm[i].cornerPoints[1] = Global_nPoint;
      localFacesComm[i].cornerPoints[2] = Global_nPoint;
      localFacesComm[i].cornerPoints[3] = Global_nPoint;
      --nFacesLocComm;
    }
  }

  /*--- Remove the invalidated faces. ---*/
  sort(localFacesComm.begin(), localFacesComm.end());
  localFacesComm.resize(nFacesLocComm);

  /*--- Complete the first round of non-blocking sends. ---*/
  SU2_MPI::Waitall(nMessSend, commReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Send the data back to the requesting ranks. ---*/
  for(unsigned long i=0; i<nMessRecv; ++i) {
    sendBufFace[i].resize(9*facesRecv[i].size());

    unsigned long ii = 0;
    for(unsigned long j=0; j<facesRecv[i].size(); ++j, ii+=9) {

      /*--- Copy the first 5 elements of the received face in the
            send buffer. ---*/
      sendBufFace[i][ii]   = facesRecv[i][j].nCornerPoints;
      sendBufFace[i][ii+1] = facesRecv[i][j].cornerPoints[0];
      sendBufFace[i][ii+2] = facesRecv[i][j].cornerPoints[1];
      sendBufFace[i][ii+3] = facesRecv[i][j].cornerPoints[2];
      sendBufFace[i][ii+4] = facesRecv[i][j].cornerPoints[3];

      /*--- Search for this face in localFacesComm. ---*/
      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFacesComm.begin(), localFacesComm.end(),
                        facesRecv[i][j]);

      /*--- Fill the remainder of the send buffer, depending on
            whether the received face corresponds to the first
            or second element. ---*/
      if(facesRecv[i][j].elemID0 == low->elemID0) {
        sendBufFace[i][ii+5] = low->elemID1;
        sendBufFace[i][ii+6] = low->nPolySol1;
        sendBufFace[i][ii+7] = low->nDOFsElem1;
        sendBufFace[i][ii+8] = low->elemType1;
      }
      else {
        sendBufFace[i][ii+5] = low->elemID0;
        sendBufFace[i][ii+6] = low->nPolySol0;
        sendBufFace[i][ii+7] = low->nDOFsElem0;
        sendBufFace[i][ii+8] = low->elemType0;
      }
    }

    /*--- Send the data using a non-blocking communication. ---*/
    SU2_MPI::Isend(sendBufFace[i].data(), sendBufFace[i].size(), MPI_UNSIGNED_LONG,
                   rankRecv[i], rankRecv[i]+1, MPI_COMM_WORLD, &commReqs[i]);
  }

  /*--- Loop over the ranks to which I originally sent my face data.
        The return data contains information about the neighboring element. ---*/
  for(unsigned long i=0; i<nMessSend; ++i) {

    /*--- Wait until a message arrives and determine the source and
          size of the message. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, MPI_COMM_WORLD, &status);
    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /*--- Allocate the memory for the receive buffer and receive the data. ---*/
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  status.MPI_SOURCE, rank+1, MPI_COMM_WORLD, &status);

    /*--- Loop over the number of received faces. ---*/
    sizeMess /= 9;
    unsigned long jj = 0;
    for(int j=0; j<sizeMess; ++j, jj+=9) {

      /*--- Create an object of type CFaceOfElement to search in the
            locally stored faces. ---*/
      CFaceOfElement thisFace;
      thisFace.nCornerPoints   = recvBuf[jj];
      thisFace.cornerPoints[0] = recvBuf[jj+1];
      thisFace.cornerPoints[1] = recvBuf[jj+2];
      thisFace.cornerPoints[2] = recvBuf[jj+3];
      thisFace.cornerPoints[3] = recvBuf[jj+4];

      /*-- Search in localFaces for this face and set the information
           of the neighboring element. ---*/
      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);
      low->elemID1    = recvBuf[jj+5];
      low->nPolySol1  = recvBuf[jj+6];
      low->nDOFsElem1 = recvBuf[jj+7];
      low->elemType1  = recvBuf[jj+8];
    }
  }

  SU2_MPI::Waitall(nMessRecv, commReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Wild cards have been used in the communication, so
        synchronize the ranks to avoid problems.          ---*/
  SU2_MPI::Barrier(MPI_COMM_WORLD);

#endif
}

void CPhysicalGeometry::DeterminePeriodicFacesFEMGrid(const CConfig          *config,
                                                      vector<CFaceOfElement> &localFaces) {

  /*--- Determine a mapping from the global point ID to the local index
        of the points.            ---*/
  map<unsigned long,unsigned long> globalPointIDToLocalInd;
  for(unsigned i=0; i<nPoint; ++i) {
    globalPointIDToLocalInd[nodes->GetGlobalIndex(i)] = i;
  }

  /*--- Loop over the number of markers present in the grid and check for a periodic one. ---*/
  for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {

      /*--- Determine the donor marker and the transformation from the
            current marker to the donor marker.    ---*/
      const unsigned short jMarker = config->GetMarker_Periodic_Donor(config->GetMarker_All_TagBound(iMarker));

      auto center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
      auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
      auto trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));

      /*--- Store (center+trans) as it is constant and will be added on. ---*/
      const su2double translation[] = {center[0] + trans[0],
                                       center[1] + trans[1],
                                       center[2] + trans[2]};

      /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
      const su2double theta = angles[0];
      const su2double phi   = angles[1];
      const su2double psi   = angles[2];

      const su2double cosTheta = cos(theta), cosPhi = cos(phi), cosPsi = cos(psi);
      const su2double sinTheta = sin(theta), sinPhi = sin(phi), sinPsi = sin(psi);

      /*--- Compute the rotation matrix. Note that the implicit
            ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
      su2double rotMatrix[3][3];
      rotMatrix[0][0] =  cosPhi*cosPsi;
      rotMatrix[1][0] =  cosPhi*sinPsi;
      rotMatrix[2][0] = -sinPhi;

      rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
      rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
      rotMatrix[2][1] = sinTheta*cosPhi;

      rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
      rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
      rotMatrix[2][2] = cosTheta*cosPhi;

      /*--- Define the vector to store the faces of the donor. Initialize its
            size to the number of local donor faces.    ---*/
      vector<CMatchingFace> facesDonor(nElem_Bound[jMarker]);

      /*------------------------------------------------------------------*/
      /*--- Step 1: Store the information of the local faces of the donor
                    marker in the variables defined above.             ---*/
      /*------------------------------------------------------------------*/

      /*--- Loop over the local elements of the donor marker. ---*/
      for(unsigned long k=0; k<nElem_Bound[jMarker]; ++k) {

        /*--- Get the connectivity of this face. The reason for the used
              function arguments for GetCornerPointsAllFaces, is that
              GetCornerPointsAllFaces is an overloaded function.   ---*/
        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long  faceConn[6][4];
        bound[jMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        /*--- Search for this face in localFaces. It must be present. ---*/
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for(unsigned short j=0; j<nPointsPerFace[0]; ++j)
          thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        /*--- Store the relevant data in facesDonor. ---*/
        facesDonor[k].nDim          = nDim;
        facesDonor[k].nCornerPoints = nPointsPerFace[0];
        facesDonor[k].elemID        = low->elemID0;
        facesDonor[k].nPoly         = low->nPolySol0;
        facesDonor[k].nDOFsElem     = low->nDOFsElem0;
        facesDonor[k].elemType      = low->elemType0;

        for(unsigned short j=0; j<nPointsPerFace[0]; ++j) {
          map<unsigned long,unsigned long>::const_iterator MI;
          MI = globalPointIDToLocalInd.find(faceConn[0][j]);
          unsigned long ind = MI->second;

          for(unsigned l=0; l<nDim; ++l)
            facesDonor[k].cornerCoor[j][l] = nodes->GetCoord(ind, l);
        }

        /*--- Create the tolerance for this face and sort the coordinates. ---*/
        facesDonor[k].SortFaceCoordinates();
      }

      /*------------------------------------------------------------------*/
      /*--- Step 2: In parallel mode the data of the donor marker is
                    gathered on all ranks.                             ---*/
      /*------------------------------------------------------------------*/

#ifdef HAVE_MPI

      /*--- Check if this is indeed a parallel simulation. ---*/
      if(size > 1) {

        /*--- Allocate the memory for the size arrays in Allgatherv. ---*/
        vector<int> recvCounts(size), displs(size);

        /*--- Create the values of recvCounts for the gather of the facesDonor. ---*/
        int sizeLocal = facesDonor.size();

        SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1,
                           MPI_INT, MPI_COMM_WORLD);

        /*--- Create the data for the vector displs from the known values of
              recvCounts. Also determine the total size of the data.   ---*/
        displs[0] = 0;
        for(int i=1; i<size; ++i) displs[i] = displs[i-1] + recvCounts[i-1];

        int sizeGlobal = displs.back() + recvCounts.back();

        /*--- SU2_MPI does not support the communication of derived data types,
              at least not at the moment when this was coded. Therefore the data
              from facesDonor is copied into buffers of primitive types, which
              are communicated. ---*/
        vector<unsigned short> shortLocBuf(5*sizeLocal);
        vector<unsigned long>  longLocBuf(sizeLocal);
        vector<su2double>      doubleLocBuf(13*sizeLocal);

        unsigned long cS=0, cL=0, cD=0;
        for(vector<CMatchingFace>::const_iterator fIt =facesDonor.begin();
                                                  fIt!=facesDonor.end(); ++fIt) {
          shortLocBuf[cS++] = fIt->nCornerPoints;
          shortLocBuf[cS++] = fIt->nDim;
          shortLocBuf[cS++] = fIt->nPoly;
          shortLocBuf[cS++] = fIt->nDOFsElem;
          shortLocBuf[cS++] = fIt->elemType;

          longLocBuf[cL++] = fIt->elemID;

          doubleLocBuf[cD++] = fIt->cornerCoor[0][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[0][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[0][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[1][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[1][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[1][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[2][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[2][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[2][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[3][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[3][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[3][2];

          doubleLocBuf[cD++] = fIt->tolForMatching;
        }

        /*--- Gather the faces from all ranks to all ranks. Use Allgatherv
              for this purpose.                  ---*/
        vector<unsigned short> shortGlobBuf(5*sizeGlobal);
        vector<unsigned long>  longGlobBuf(sizeGlobal);
        vector<su2double>      doubleGlobBuf(13*sizeGlobal);

        SU2_MPI::Allgatherv(longLocBuf.data(), longLocBuf.size(), MPI_UNSIGNED_LONG,
                            longGlobBuf.data(), recvCounts.data(), displs.data(),
                            MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

        for(int i=0; i<size; ++i) {
          recvCounts[i] *= 5; displs[i] *= 5;
        }

        SU2_MPI::Allgatherv(shortLocBuf.data(), shortLocBuf.size(), MPI_UNSIGNED_SHORT,
                            shortGlobBuf.data(), recvCounts.data(), displs.data(),
                            MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);

        for(int i=0; i<size; ++i) {
          recvCounts[i] /=  5; displs[i] /=  5;
          recvCounts[i] *= 13; displs[i] *= 13;
        }

        SU2_MPI::Allgatherv(doubleLocBuf.data(), doubleLocBuf.size(), MPI_DOUBLE,
                            doubleGlobBuf.data(), recvCounts.data(), displs.data(),
                            MPI_DOUBLE, MPI_COMM_WORLD);

        /*--- Copy the data back into facesDonor, which will contain the
              global information after the copies. ---*/
        facesDonor.resize(sizeGlobal);
        cS = cL = cD = 0;

        for(vector<CMatchingFace>::iterator fIt =facesDonor.begin();
                                            fIt!=facesDonor.end(); ++fIt) {
          fIt->nCornerPoints = shortGlobBuf[cS++];
          fIt->nDim          = shortGlobBuf[cS++];
          fIt->nPoly         = shortGlobBuf[cS++];
          fIt->nDOFsElem     = shortGlobBuf[cS++];
          fIt->elemType      = shortGlobBuf[cS++];

          fIt->elemID = longGlobBuf[cL++];

          fIt->cornerCoor[0][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[0][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[0][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[1][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[1][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[1][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[2][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[2][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[2][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[3][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[3][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[3][2] = doubleGlobBuf[cD++];

          fIt->tolForMatching = doubleGlobBuf[cD++];
        }
      }
#endif

      /*--- Sort facesDonor in increasing order. ---*/
      sort(facesDonor.begin(), facesDonor.end());

      /*------------------------------------------------------------------*/
      /*--- Step 3: Find for the current marker the required data in the
                    variables for the donor marker.                    ---*/
      /*------------------------------------------------------------------*/

      /*--- Loop over the local faces of this boundary marker. ---*/
      for(unsigned long k=0; k<nElem_Bound[iMarker]; ++k) {

        /*--- Get the connectivity of this face. The reason for the used
              function arguments for GetCornerPointsAllFaces, is that
              GetCornerPointsAllFaces is an overloaded function.   ---*/
        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long  faceConn[6][4];
        bound[iMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        /*--- Search for this face in localFaces. It must be present. ---*/
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for(unsigned short j=0; j<nPointsPerFace[0]; ++j)
          thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        /*--- Indicate that this face is a periodic face. This is accomplished
              by setting periodicIndex to iMarker + 1.       ---*/
        low->periodicIndex = iMarker + 1;

        /*--- Store the data for this boundary element also in thisMatchingFace,
              such that a search can be carried out in donorFaces. Note that the
              periodic transformation must be applied to the coordinates. ---*/
        CMatchingFace thisMatchingFace;
        thisMatchingFace.nDim          = nDim;
        thisMatchingFace.nCornerPoints = nPointsPerFace[0];
        thisMatchingFace.elemID        = low->elemID0;
        thisMatchingFace.nPoly         = low->nPolySol0;
        thisMatchingFace.nDOFsElem     = low->nDOFsElem0;
        thisMatchingFace.elemType      = low->elemType0;

        for(unsigned short j=0; j<nPointsPerFace[0]; ++j) {
          map<unsigned long,unsigned long>::const_iterator MI;
          MI = globalPointIDToLocalInd.find(faceConn[0][j]);
          unsigned long ind = MI->second;
          const su2double *coor = nodes->GetCoord(ind);

          const su2double dx =             coor[0] - center[0];
          const su2double dy =             coor[1] - center[1];
          const su2double dz = nDim == 3 ? coor[2] - center[2] : su2double(0.0);

          thisMatchingFace.cornerCoor[j][0] = rotMatrix[0][0]*dx + rotMatrix[0][1]*dy
                                            + rotMatrix[0][2]*dz + translation[0];
          thisMatchingFace.cornerCoor[j][1] = rotMatrix[1][0]*dx + rotMatrix[1][1]*dy
                                            + rotMatrix[1][2]*dz + translation[1];
          thisMatchingFace.cornerCoor[j][2] = rotMatrix[2][0]*dx + rotMatrix[2][1]*dy
                                            + rotMatrix[2][2]*dz + translation[2];
        }

        /*--- Create the tolerance for this face and sort the coordinates. ---*/
        thisMatchingFace.SortFaceCoordinates();

        /*--- Check if thisMatchingFace is present in facesDonor. If so, set the
              missing information in the face corresponding to the iterator low. ---*/
        vector<CMatchingFace>::const_iterator donorLow;
        donorLow = lower_bound(facesDonor.begin(), facesDonor.end(), thisMatchingFace);

        if(donorLow != facesDonor.end()) {
          if( !(thisMatchingFace < *donorLow) ) {
            low->elemID1            = donorLow->elemID;
            low->nPolySol1          = donorLow->nPoly;
            low->nDOFsElem1         = donorLow->nDOFsElem;
            low->elemType1          = donorLow->elemType;
            low->periodicIndexDonor = jMarker + 1;
          }
        }
      }
    }
  } 
}
