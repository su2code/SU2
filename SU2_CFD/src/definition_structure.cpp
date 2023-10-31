/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD
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


#include "../include/definition_structure.hpp"


void Partition_Analysis(CGeometry *geometry, CConfig *config) {

  /*--- This routine does a quick and dirty output of the total
   vertices, ghost vertices, total elements, ghost elements, etc.,
   so that we can analyze the partition quality. ---*/

  unsigned short nMarker = config->GetnMarker_All();
  unsigned short iMarker, iNodes, MarkerS, MarkerR;
  unsigned long iElem, iPoint, nVertexS, nVertexR;
  unsigned long nNeighbors = 0, nSendTotal = 0, nRecvTotal = 0;
  unsigned long nPointTotal=0, nPointGhost=0, nElemTotal=0;
  unsigned long nElemHalo=0, nEdge=0, nElemBound=0;
  int iRank;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;

#ifdef HAVE_MPI
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);
#endif

  nPointTotal = geometry->GetnPoint();
  nPointGhost = geometry->GetnPoint() - geometry->GetnPointDomain();
  nElemTotal  = geometry->GetnElem();
  nEdge       = geometry->GetnEdge();

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    nElemBound  += geometry->GetnElem_Bound(iMarker);
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      MarkerS = iMarker;  MarkerR = iMarker+1;
      nVertexS = geometry->nVertex[MarkerS];
      nVertexR = geometry->nVertex[MarkerR];
      nNeighbors++;
      nSendTotal += nVertexS;
      nRecvTotal += nVertexR;
    }
  }

  bool *isHalo = new bool[geometry->GetnElem()];
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    isHalo[iElem] = false;
    for (iNodes = 0; iNodes < geometry->elem[iElem]->GetnNodes(); iNodes++) {
      iPoint = geometry->elem[iElem]->GetNode(iNodes);
      if (!geometry->nodes->GetDomain(iPoint)) isHalo[iElem] = true;
    }
  }

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if (isHalo[iElem]) nElemHalo++;
  }

  unsigned long *row_ptr = nullptr, nnz;
  unsigned short *nNeigh = nullptr;
  vector<unsigned long> vneighs;

  /*--- Don't delete *row_ptr, *col_ind because they are
   asigned to the Jacobian structure. ---*/

  /*--- Compute the number of neighbors ---*/

  nNeigh = new unsigned short [geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    // +1 -> to include diagonal element
    nNeigh[iPoint] = (geometry->nodes->GetnPoint(iPoint)+1);
  }

  /*--- Create row_ptr structure, using the number of neighbors ---*/

  row_ptr = new unsigned long [geometry->GetnPoint()+1];
  row_ptr[0] = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    row_ptr[iPoint+1] = row_ptr[iPoint] + nNeigh[iPoint];
  nnz = row_ptr[geometry->GetnPoint()];

  delete [] row_ptr;
  delete [] nNeigh;

  /*--- Now put this info into a CSV file for processing ---*/

  ofstream Profile_File;
  Profile_File.precision(15);

  if (rank == MASTER_NODE) {
    /*--- Prepare and open the file ---*/
    Profile_File.open("partitioning.csv");
    /*--- Create the CSV header ---*/
    Profile_File << R"("Rank", "nNeighbors", "nPointTotal", "nEdge", "nPointGhost", "nSendTotal", "nRecvTotal", "nElemTotal", "nElemBoundary", "nElemHalo", "nnz")" << endl;
    Profile_File.close();
  }
  SU2_MPI::Barrier(SU2_MPI::GetComm());

  /*--- Loop through the map and write the results to the file ---*/

  for (iRank = 0; iRank < size; iRank++) {
    if (rank == iRank) {
      Profile_File.open("partitioning.csv", ios::out | ios::app);
      Profile_File << rank << ", " << nNeighbors << ", " << nPointTotal << ", " << nEdge << "," << nPointGhost << ", " << nSendTotal << ", " << nRecvTotal << ", " << nElemTotal << "," << nElemBound << ", " << nElemHalo << ", " << nnz << endl;
      Profile_File.close();
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  delete [] isHalo;

}

void Partition_Analysis_FEM(CGeometry *geometry, CConfig *config) {

  /*--- This routine does a quick and dirty output of the total
   vertices, ghost vertices, total elements, ghost elements, etc.,
   so that we can analyze the partition quality. ---*/

  unsigned long nNeighSend = 0, nNeighRecv     = 0;
  unsigned long nElemOwned = 0, nElemSendTotal = 0, nElemRecvTotal = 0;
  unsigned long nDOFOwned  = 0, nDOFSendTotal  = 0, nDOFRecvTotal  = 0;

  int iRank;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;

#ifdef HAVE_MPI
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);
#endif

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  auto *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  /*--- Determine the number of owned elements and DOFs. ---*/
  nElemOwned = nVolElemOwned;
  for(unsigned long l=0; l<nVolElemOwned; ++l) {
    nDOFOwned += volElem[l].nDOFsSol;
  }

  /*--- Get the communication information from DG_Geometry. Note that for a
   FEM DG discretization the communication entities of FEMGeometry contain
   the volume elements. ---*/
  const vector<int>                    &ranksSend    = DGGeometry->GetRanksSend();
  const vector<int>                    &ranksRecv    = DGGeometry->GetRanksRecv();
  const vector<vector<unsigned long> > &elementsSend = DGGeometry->GetEntitiesSend();
  const vector<vector<unsigned long> > &elementsRecv = DGGeometry->GetEntitiesRecv();

  nNeighSend = ranksSend.size();
  nNeighRecv = ranksRecv.size();

  /*--- Determine the total number of elements and DOFS to be send. ---*/
  for(unsigned long i=0; i<ranksSend.size(); ++i) {

    const auto nElemSend = (unsigned int)elementsSend[i].size();

    nElemSendTotal += nElemSend;

    for(unsigned int j=0; j<nElemSend; ++j) {
      const unsigned long jj = elementsSend[i][j];
      nDOFSendTotal += volElem[jj].nDOFsSol;
    }
  }

  /*--- Determine the total number of elements and DOFS to be received. ---*/
  for(unsigned long i=0; i<ranksRecv.size(); ++i) {

    const auto nElemRecv = (unsigned int)elementsRecv[i].size();

    nElemRecvTotal += nElemRecv;

    for(unsigned int j=0; j<nElemRecv; ++j) {
      const unsigned long jj = elementsRecv[i][j];
      nDOFRecvTotal += volElem[jj].nDOFsSol;
    }

  }

  /*--- Now put this info into a CSV file for processing ---*/

  ofstream Profile_File;
  Profile_File.precision(15);

  if (rank == MASTER_NODE) {
    /*--- Prepare and open the file ---*/
    Profile_File.open("partitioning.csv");
    /*--- Create the CSV header ---*/
    Profile_File << R"("Rank", "nNeighSend",  "nNeighRecv", "nElemOwned", "nElemSendTotal", "nElemRecvTotal", "nDOFOwned", "nDOFSendTotal", "nDOFRecvTotal")" << endl;
    Profile_File.close();
  }
  SU2_MPI::Barrier(SU2_MPI::GetComm());

  /*--- Loop through the map and write the results to the file ---*/

  for (iRank = 0; iRank < size; iRank++) {
    if (rank == iRank) {
      Profile_File.open("partitioning.csv", ios::out | ios::app);
      Profile_File << rank << ", " << nNeighSend << ", " << nNeighRecv << ", " << nElemOwned << ", "
                   << nElemSendTotal << ", " << nElemRecvTotal << ", " << nDOFOwned << ", "
                   << nDOFSendTotal << ", " << nDOFRecvTotal << endl;
      Profile_File.close();
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

}
