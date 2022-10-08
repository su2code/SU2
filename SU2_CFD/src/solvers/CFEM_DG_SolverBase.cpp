/*!
 * \file CFEM_DG_SolverBase.cpp
 * \brief Main subroutines for the base class of the finite element flow solvers.
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 7.1.1 "Blackbird"
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


#include "../../include/solvers/CFEM_DG_SolverBase.hpp"
#include "../../../Common/include/toolboxes/CGraphColoringStructure.hpp"

CFEM_DG_SolverBase::CFEM_DG_SolverBase(void)
  : CSolver() {
}

CFEM_DG_SolverBase::CFEM_DG_SolverBase(CGeometry      *geometry,
                                       CConfig        *config,
                                       unsigned short iMesh)
  : CSolver() {

  /*--- Store the multigrid level and retrieve the number of dimensions
        and the number of markers. ---*/
  MGLevel = iMesh;
  nDim    = geometry->GetnDim();
  nMarker = config->GetnMarker_All();

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  if( !DGGeometry) SU2_MPI::Error(string("Dynamic cast failed"), CURRENT_FUNCTION);

  nVolElemTot   = DGGeometry->GetNVolElemTot();
  nVolElemOwned = DGGeometry->GetNVolElemOwned();
  volElem       = DGGeometry->GetVolElem();

  nVolElemOwnedPerTimeLevel    = DGGeometry->GetNVolElemOwnedPerTimeLevel();
  nVolElemInternalPerTimeLevel = DGGeometry->GetNVolElemInternalPerTimeLevel();
  nVolElemHaloPerTimeLevel     = DGGeometry->GetNVolElemHaloPerTimeLevel();

  ownedElemAdjLowTimeLevel = DGGeometry->GetOwnedElemAdjLowTimeLevel();
  haloElemAdjLowTimeLevel  = DGGeometry->GetHaloElemAdjLowTimeLevel();

  nMeshPoints = DGGeometry->GetNMeshPoints();
  meshPoints  = DGGeometry->GetMeshPoints();

  nMatchingInternalFacesWithHaloElem = DGGeometry->GetNMatchingFacesWithHaloElem();
  nMatchingInternalFacesLocalElem    = DGGeometry->GetNMatchingFacesInternal();
  matchingInternalFaces              = DGGeometry->GetMatchingFaces();

  boundaries = DGGeometry->GetBoundaries();

  timeCoefADER_DG                        = DGGeometry->GetTimeCoefADER_DG();
  timeInterpolDOFToIntegrationADER_DG    = DGGeometry->GetTimeInterpolDOFToIntegrationADER_DG();
  timeInterpolAdjDOFToIntegrationADER_DG = DGGeometry->GetTimeInterpolAdjDOFToIntegrationADER_DG();

  /*--- Determine the number of owned and total number of DOFs stored
        on this rank. Note that for this the number of DOFs of the
        standard element of the flow is taken. For incompressible flows
        the number of DOFs is different for the pressure, but that is
        taken into account. ---*/
  nDOFsLocOwned = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    nDOFsLocOwned += volElem[i].standardElemFlow->GetNDOFs();

  nDOFsLocTot = nDOFsLocOwned;
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i)
    nDOFsLocTot += volElem[i].standardElemFlow->GetNDOFs();

  /*--- Determine the global number of DOFs in the simulation. ---*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nDOFsLocOwned, &nDOFsGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#else
  nDOFsGlobal = nDOFsLocOwned;
#endif

  /*--- Store the number of DOFs in the geometry class in case of restart. ---*/
  geometry->SetnPointDomain(nDOFsLocOwned);
  geometry->SetGlobal_nPointDomain(nDOFsGlobal);

  /*--- Determine the maximum number of integration points and DOFs
        used on this rank. ---*/
  const unsigned int nIntegrationMax = DGGeometry->DetermineMaxNIntegration();
  const unsigned int nDOFsMax        = DGGeometry->DetermineMaxNDOFs();

  /*--- Allocate the memory for the aerodynamic coefficients. ---*/
  InvCoeff.allocate(nMarker);
  SurfaceInvCoeff.allocate(nMarker);
  ViscCoeff.allocate(nMarker);
  SurfaceViscCoeff.allocate(nMarker);
  SurfaceCoeff.allocate(nMarker);

  /*--- Initialize the total coefficients. ---*/
  AllBoundInvCoeff.setZero();
  AllBoundViscCoeff.setZero();
  TotalCoeff.setZero();

  /*--- Determine for all elements the adjacent internal matching faces,
        i.e. the faces that contribute to the residual of the elements. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
  adjMatchingFacesOfVolElem.resize(nVolElemTot);
  for(unsigned long i=0; i<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++i) {
    adjMatchingFacesOfVolElem[matchingInternalFaces[i].elemID0].push_back(i);
    adjMatchingFacesOfVolElem[matchingInternalFaces[i].elemID1].push_back(i);
  }

  /*--- Determine for all elements the adjacent physical surfaces, i.e.
        on physical boundaries, that contribute to the residual of
        the elements. Note that for the boundary surfaces the halo
        elements do not need to be considered. ---*/
  adjSurfacesOfVolElem.resize(nVolElemOwned);
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    adjSurfacesOfVolElem[i].resize(nMarker);

  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    if( !(boundaries[iMarker].periodicBoundary) ) {
      vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
      for(unsigned long i=0; i<surfElem.size(); ++i)
        adjSurfacesOfVolElem[surfElem[i].volElemID][iMarker].push_back(i);
    }
  }

  /*--- Check if the exact Jacobian of the spatial discretization must be
        determined. If so, the color of each DOF must be determined, which
        is converted to the DOFs for each color. Note that for incompressible
        this is a slight overhead, because the pressure has less DOFs than
        the velocity. ---*/
  if( config->GetJacobian_Spatial_Discretization_Only() ) {

    /*--- Write a message that the graph coloring is performed. ---*/
    if(rank == MASTER_NODE)
      cout << "Creating the vertex colors of the graph. " << std::endl;

    /*--- Determine the graph of the stencil of every local DOF. ---*/
    DetermineGraphDOFs(DGGeometry, config);

    /*--- Carry out the vertex coloring of the graph. ---*/
    CGraphColoringStructure graphColoring;
    vector<int> colorLocalDOFs;
    graphColoring.GraphVertexColoring(config, nDOFsPerRank, nonZeroEntriesJacobian,
                                      nGlobalColors, colorLocalDOFs);

    /*--- Write a message that the all volume DOFs have been colored. ---*/
    if(rank == MASTER_NODE)
      cout << "There are " << nGlobalColors
           << " colors present in the graph for " << nDOFsPerRank.back()
           << " DOFs." << std::endl;

    /*--- Determine the meta data needed for the computation of the Jacobian. ---*/
    MetaDataJacobianComputation(DGGeometry, colorLocalDOFs);
  }

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "DG FLOW";
}

CFEM_DG_SolverBase::~CFEM_DG_SolverBase(void) {

  for(auto& model : FluidModel) delete model;
}

void CFEM_DG_SolverBase::DetermineGraphDOFs(const CMeshFEM_DG *DGGeometry,
                                            CConfig           *config) {

  /*-------------------------------------------------------------------*/
  /* Step 1: Determine the number of owned DOFs per rank in cumulative */
  /*         storage format.                                           */
  /*********************************************************************/

  /*--- Gather the number of DOFs from all ranks. */
  nDOFsPerRank.resize(size+1);
  nDOFsPerRank[0] = 0;

#ifdef HAVE_MPI
  SU2_MPI::Allgather(&nDOFsLocOwned, 1, MPI_UNSIGNED_LONG, &nDOFsPerRank[1], 1,
                     MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
#else
  nDOFsPerRank[1] = nDOFsLocOwned;
#endif

  /*--- Put nDOFsPerRank in cumulative storage format. ---*/
  for(int i=0; i<size; ++i)
    nDOFsPerRank[i+1] += nDOFsPerRank[i];

  /*-------------------------------------------------------------------*/
  /* Step 2: Determine the global numbering of the DOFs, especially    */
  /*         of the externals.                                         */
  /*********************************************************************/

  /*--- Determine the local offset of the DOFs of the elements. ---*/
  vector<unsigned long> offsetDOFsElemLoc(nVolElemTot+1);
  offsetDOFsElemLoc[0] = 0;
  for(unsigned long i=0; i<nVolElemTot; ++i)
    offsetDOFsElemLoc[i+1] = offsetDOFsElemLoc[i]
                           + volElem[i].standardElemFlow->GetNDOFs();

  /*--- Determine the global offset of the DOF numbering for the
        locally owned elements. Allocate it for all elements,
        because this information is communicated. ---*/
  vector<unsigned long> offsetDOFsElem(nVolElemTot+1);
  for(unsigned long i=0; i<=nVolElemOwned; ++i)
    offsetDOFsElem[i] = offsetDOFsElemLoc[i] + nDOFsPerRank[rank];

  /*--- Get the communication information from DG_Geometry. Note that for a
        FEM DG discretization the communication entities of DGGeometry contain
        the volume elements. ---*/
  const vector<int>                    &ranksSend    = DGGeometry->GetRanksSend();
  const vector<int>                    &ranksRecv    = DGGeometry->GetRanksRecv();
  const vector<vector<unsigned long> > &elementsSend = DGGeometry->GetEntitiesSend();
  const vector<vector<unsigned long> > &elementsRecv = DGGeometry->GetEntitiesRecv();

#ifdef HAVE_MPI

  /*--- Parallel implementation. Allocate the memory for the send buffers and
        loop over the number of ranks to which I have to send data. Note that
        self communication is not excluded here. ---*/
  vector<vector<unsigned long> > sendBuf;
  sendBuf.resize(ranksSend.size());

  vector<SU2_MPI::Request> sendReqs(ranksSend.size());

  for(unsigned i=0; i<ranksSend.size(); ++i) {

    /*--- Fill the send buffer for this rank. ---*/
    for(unsigned long j=0; j<elementsSend[i].size(); ++j)
      sendBuf[i].push_back(offsetDOFsElem[elementsSend[i][j]]);

    /*--- Send the data using non-blocking sends to avoid deadlock. ---*/
    int dest = ranksSend[i];
    SU2_MPI::Isend(sendBuf[i].data(), sendBuf[i].size(), MPI_UNSIGNED_LONG,
                   dest, dest, SU2_MPI::GetComm(), &sendReqs[i]);
  }

  /*--- Create a map of the receive rank to the index in ranksRecv. ---*/
  map<int,int> rankToIndRecvBuf;
  for(int i=0; i<(int) ranksRecv.size(); ++i)
    rankToIndRecvBuf[ranksRecv[i]] = i;

  /*--- Loop over the number of ranks from which I receive data. ---*/
  for(unsigned long i=0; i<ranksRecv.size(); ++i) {

    /*--- Block until a message arrives and determine the source
          and size of the message. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /*--- Allocate the memory for the receive buffer and receive the data. ---*/
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank, SU2_MPI::GetComm(), &status);

    /*--- Determine the actual index of this rank in ranksRecv. ---*/
    map<int,int>::const_iterator MI = rankToIndRecvBuf.find(source);
    source = MI->second;

    /*--- Loop over the data received and set the value of
          offsetDOFsElem accordingly. ---*/
    for(int j=0; j<sizeMess; ++j)
      offsetDOFsElem[elementsRecv[source][j]] = recvBuf[j];
  }

  /*--- Complete the non-blocking sends. ---*/
  SU2_MPI::Waitall(ranksSend.size(), sendReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Wild cards have been used in the communication,
        so synchronize the ranks to avoid problems. ---*/
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential implementation. Halo's may be present due to periodic
        boundary conditions. A loop over the receiving ranks is carried out
        to cover both situations, i.e. halo's present and halo's not present. ---*/
  for(unsigned long i=0; i<ranksRecv.size(); ++i)
    for(unsigned long j=0; j<elementsRecv[i].size(); ++j)
      offsetDOFsElem[elementsRecv[i][j]] = offsetDOFsElem[elementsSend[i][j]];

#endif

  /*-------------------------------------------------------------------*/
  /* Step 3: Determine the entries of the graph for the locally stored */
  /*         data. Note that this data must also be determined for the */
  /*         halo DOFs, because the faces between partitions are       */
  /*         stored only once.                                         */
  /*********************************************************************/

  /*--- Allocate the memory for the first index of nonZeroEntriesJacobian.
        During the construction of this data, also the halo DOFs are needed,
        hence it is allocate for all DOFs stored on this rank. ---*/
  nonZeroEntriesJacobian.resize(nDOFsLocTot);

  /*--- The residual of the DOFs of a volume element depends on all the
        DOFs of that volume element. Set these dependencies. Only needed
        for the owned elements. ---*/
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    const unsigned short nDOFs = volElem[i].standardElemFlow->GetNDOFs();
    for(unsigned short j=0; j<nDOFs; ++j) {
      const unsigned long jj = offsetDOFsElemLoc[i] + j;
      for(unsigned short k=0; k<nDOFs; ++k)
        nonZeroEntriesJacobian[jj].push_back(offsetDOFsElem[i] + k);
    }
  }

  /*--- Loop over the internal matching faces to set the dependencies of the
        DOFs w.r.t. the DOFs of the neighbors. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
  for(unsigned long i=0; i<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++i) {

    /*--- Easier storage of the elements on both sides and the
          number of DOFs on both elements. ---*/
    const unsigned long elem0 = matchingInternalFaces[i].elemID0;
    const unsigned long elem1 = matchingInternalFaces[i].elemID1;

    const unsigned short nDOFs0 = volElem[elem0].standardElemFlow->GetNDOFs();
    const unsigned short nDOFs1 = volElem[elem1].standardElemFlow->GetNDOFs();

    /*--- Double loop over the DOFs of both elements and set
          the dependencies accordingly. ---*/
    for(unsigned short j=0; j<nDOFs0; ++j) {
      const unsigned long jj = offsetDOFsElemLoc[elem0] + j;
      for(unsigned short k=0; k<nDOFs1; ++k)
        nonZeroEntriesJacobian[jj].push_back(offsetDOFsElem[elem1]+k);
    }

    for(unsigned short j=0; j<nDOFs1; ++j) {
      const unsigned long jj = offsetDOFsElemLoc[elem1] + j;
      for(unsigned short k=0; k<nDOFs0; ++k)
        nonZeroEntriesJacobian[jj].push_back(offsetDOFsElem[elem0]+k);
    }
  }

  /*-------------------------------------------------------------------*/
  /* Step 4: Carry out the communication for the halo DOFs, such that  */
  /*         all the graph information is obtained for the owned DOFs. */
  /*********************************************************************/

#ifdef HAVE_MPI

  /*--- Parallel implementation. Allocate the memory for the inverse send buffers
        and loop over the number of ranks to which I have to send data for the
        inverse communication pattern. Note that self communication is not
        excluded here, because efficiency is not really important for this function. ---*/
  vector<vector<unsigned long> > invSendBuf;
  invSendBuf.resize(ranksRecv.size());

  vector<SU2_MPI::Request> invSendReqs(ranksRecv.size());

  for(unsigned long i=0; i<ranksRecv.size(); ++i) {

    /*--- Fill the inverse send buffer for this rank. ---*/
    for(unsigned long j=0; j<elementsRecv[i].size(); ++j) {
      const unsigned long jj = elementsRecv[i][j];
      const unsigned short nDOFs = volElem[jj].standardElemFlow->GetNDOFs();
      for(unsigned short k=0; k<nDOFs; ++k) {
        const unsigned long kk = offsetDOFsElemLoc[jj] + k;
        invSendBuf[i].push_back(nonZeroEntriesJacobian[kk].size());
        invSendBuf[i].insert(invSendBuf[i].end(),
                             nonZeroEntriesJacobian[kk].begin(),
                             nonZeroEntriesJacobian[kk].end());
      }
    }

    /*--- Send the data using non-blocking sends to avoid deadlock. ---*/
    int dest = ranksRecv[i];
    SU2_MPI::Isend(invSendBuf[i].data(), invSendBuf[i].size(), MPI_UNSIGNED_LONG,
                   dest, dest+1, SU2_MPI::GetComm(), &invSendReqs[i]);
  }

  /*--- Create a map of the inverse receive (i.e. the original send) rank
        to the index in ranksSend. ---*/
  map<int,int> rankToIndSendBuf;
  for(int i=0; i<(int) ranksSend.size(); ++i)
    rankToIndSendBuf[ranksSend[i]] = i;

  /*--- Loop over the number of ranks from which I receive data
        in the inverse communication. ---*/
  for(unsigned long i=0; i<ranksSend.size(); ++i) {

    /* Block until a message arrives and determine the source and size
       of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /*--- Allocate the memory for the receive buffer, receive the message
          and determine the actual index of this rank in ranksSend. ---*/
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank+1, SU2_MPI::GetComm(), &status);

    map<int,int>::const_iterator MI = rankToIndSendBuf.find(source);
    source = MI->second;

    /*--- Loop over the elements for which the data was just received. ---*/
    int ii = 0;
    for(unsigned long j=0; j<elementsSend[source].size(); ++j) {
      const unsigned long jj = elementsSend[source][j];
      const unsigned short nDOFs = volElem[jj].standardElemFlow->GetNDOFs();
      for(unsigned short k=0; k<nDOFs; ++k) {
        const unsigned long kk       = offsetDOFsElemLoc[jj] + k;
        const unsigned long nEntries = recvBuf[ii++];

        nonZeroEntriesJacobian[kk].insert(nonZeroEntriesJacobian[kk].end(),
                                          recvBuf.data() + ii,
                                          recvBuf.data() + ii + nEntries);
        ii += nEntries;
      }
    }
  }

  /*--- Complete the non-blocking sends. ---*/
  SU2_MPI::Waitall(ranksRecv.size(), invSendReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Wild cards have been used in the communication,
        so synchronize the ranks to avoid problems. ---*/
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else
  /*--- Sequential implementation. Just add the data of the halo DOFs
        to the corresponding owned DOFs. Note that only when periodic
        boundary conditions are present, something will happen. ---*/
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    for(unsigned long j=0; j<elementsSend[i].size(); ++j) {
      const unsigned long elemR = elementsRecv[i][j];
      const unsigned long elemS = elementsSend[i][j];

      const unsigned short nDOFs = volElem[elemR].standardElemFlow->GetNDOFs();
      for(unsigned short k=0; k<nDOFs; ++k) {
        const unsigned long kR = offsetDOFsElemLoc[elemR] + k;
        const unsigned long kS = offsetDOFsElemLoc[elemS] + k;

        nonZeroEntriesJacobian[kS].insert(nonZeroEntriesJacobian[kS].end(),
                                          nonZeroEntriesJacobian[kR].begin(),
                                          nonZeroEntriesJacobian[kR].end());
      }
    }
  }

#endif

  /*-------------------------------------------------------------------*/
  /* Step 5: Sort the entries of nonZeroEntriesJacobian in increasing  */
  /*         order and remove possible double entries (caused by       */
  /*         periodic boundaries under certain circumstances). Delete  */
  /*         the memory of the halo data, because it is not needed     */
  /*         anymore.                                                  */
  /*********************************************************************/

  /*--- Loop over the locally owned DOFs and remove the possible double
        entries. This is accomplished by first sorting the data followed
        by the removal of duplicate entries using the unique function. ---*/
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
    sort(nonZeroEntriesJacobian[i].begin(), nonZeroEntriesJacobian[i].end());
    vector<unsigned long>::iterator lastEntry;
    lastEntry = unique(nonZeroEntriesJacobian[i].begin(), nonZeroEntriesJacobian[i].end());
    nonZeroEntriesJacobian[i].erase(lastEntry, nonZeroEntriesJacobian[i].end());
  }

  /*--- Delete the memory of the halo data. To make sure that the memory
        is physically deleted, use the swap function. ---*/
  for(unsigned long i=nDOFsLocOwned; i<nDOFsLocTot; ++i)
    vector<unsigned long>().swap(nonZeroEntriesJacobian[i]);
}

void CFEM_DG_SolverBase::MetaDataJacobianComputation(const CMeshFEM_DG *DGGeometry,
                                                     const vector<int> &colorLocalDOFs) {

  /*--------------------------------------------------------------------------*/
  /*--- Part 1: Convert the coloring information of the DOFs to the        ---*/
  /*---         reverse information, namely the DOFs for each color.       ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the locally owned DOFs for each color. */
  localDOFsPerColor.resize(nGlobalColors);

  for(unsigned long i=0; i<nDOFsLocOwned; ++i)
    localDOFsPerColor[colorLocalDOFs[i]].push_back(i);

  /*--------------------------------------------------------------------------*/
  /*--- Part 2: Create the mapping from the global matrix index to the     ---*/
  /*---         color. In parallel this is a bit more work than one might  ---*/
  /*---         think, because the parallel data structures do not store   ---*/
  /*---         all the matrix entities of the local DOFs. Hence, this     ---*/
  /*---         must be reconstructed via communication.                   ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Get the send and receive ranks for DGGeometry. ---*/
  const vector<int> &ranksSend = DGGeometry->GetRanksSend();
  const vector<int> &ranksRecv = DGGeometry->GetRanksRecv();

  /*--- Store the ranks this rank communicates with in a map. These are both
        sending and receiving ranks. Make sure to exclude my own rank. ---*/
  map<int,int> rankCommToInd;
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    if(ranksSend[i] != rank) {
      const int ind = (int)rankCommToInd.size();
      rankCommToInd[ranksSend[i]] = ind;
    }
  }

  for(unsigned long i=0; i<ranksRecv.size(); ++i) {
    map<int,int>::const_iterator MI = rankCommToInd.find(ranksRecv[i]);
    if((MI == rankCommToInd.end()) &&(ranksRecv[i] != rank)) {
      const int ind = (int)rankCommToInd.size();
      rankCommToInd[ranksRecv[i]] = ind;
    }
  }

  /*--- Allocate the memory for the first index of the send buffers. ---*/
  vector<vector<unsigned long> > sendBuf(rankCommToInd.size(), vector<unsigned long>(0));

  /*--- Define the map from the global matrix index to the color and the set to
        keep track whether an index has already been treated. ---*/
  map<unsigned long, int> mapMatrixIndToColor;
  set<unsigned long> setTreatedIndices;

  /*--- Loop over the owned DOFs and its matrix entries to create the contents of
        either mapMatrixIndToColor (owned DOFs) or sendBuf (unowned DOFs). ---*/
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
    for(unsigned long j=0; j<nonZeroEntriesJacobian[i].size(); ++j) {
      const unsigned long jj = nonZeroEntriesJacobian[i][j];

      /*--- Check if jj has not been stored yet. ---*/
      pair<set<unsigned long>::iterator,bool> attempt = setTreatedIndices.insert(jj);
      if( attempt.second ) {

        /*--- Determine the rank where this DOF is stored. ---*/
        vector<unsigned long>::iterator low;
        low = lower_bound(nDOFsPerRank.begin(), nDOFsPerRank.end(), jj);
        int rankDOF = (int)(low - nDOFsPerRank.begin());
        if(*low > jj) --rankDOF;

        /*--- If rankDOF is the current rank, create the entry in
              mapMatrixIndToColor. Otherwise store jj in the appropriate send
              buffer. ---*/
        if(rankDOF == rank) {
          mapMatrixIndToColor[jj] = colorLocalDOFs[jj-nDOFsPerRank[rank]];
        }
        else {
          map<int,int>::const_iterator MI = rankCommToInd.find(rankDOF);
          sendBuf[MI->second].push_back(jj);
        }
      }
    }
  }

#ifdef HAVE_MPI

  /*--- Loop over the ranks to which this rank has to send data.
        Use non-blocking sends to avoid deadlock. ---*/
  vector<SU2_MPI::Request> sendReqs(rankCommToInd.size());

  map<int,int>::const_iterator MI = rankCommToInd.begin();
  for(unsigned long i=0; i<rankCommToInd.size(); ++i, ++MI) {
    const int dest = MI->first;
    const int ind  = MI->second;

    SU2_MPI::Isend(sendBuf[ind].data(), sendBuf[ind].size(), MPI_UNSIGNED_LONG,
                   dest, dest+2, SU2_MPI::GetComm(), &sendReqs[i]);
  }

  /*--- Loop over the ranks from which I receive data to be processed. The number
        of ranks is equal to the number of ranks to which I just sent data. ---*/
  vector<SU2_MPI::Request> sendReturnReqs(rankCommToInd.size());
  vector<vector<int> > sendReturnBuf(rankCommToInd.size(), vector<int>(0));
  for(unsigned long i=0; i<rankCommToInd.size(); ++i) {

    /*--- Block until a message arrives and determine the source and size
          of the message. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+2, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /*--- Allocate the memory for the receive buffer as well as for the
          return send buffer. Receive the message afterwards. ---*/
    vector<unsigned long> recvBuf(sizeMess);
    sendReturnBuf[i].resize(sizeMess);

    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank+2, SU2_MPI::GetComm(), &status);

    /*--- Loop over the data just received and fill the return send buffer
          with the color of the DOFs. ---*/
    for(int j=0; j<sizeMess; ++j) {
      const unsigned long jj = recvBuf[j] - nDOFsPerRank[rank];
      if(jj >= nDOFsLocOwned)
        SU2_MPI::Error(string("This DOF should be owned, but it is not. Should not happen."),
                       CURRENT_FUNCTION);
      sendReturnBuf[i][j] = colorLocalDOFs[jj];
    }

    /*--- Send the return buffer back to the calling rank. Again use non-blocking
          sends to avoid deadlock. ---*/
    SU2_MPI::Isend(sendReturnBuf[i].data(), sendReturnBuf[i].size(), MPI_INT,
                   source, source+3, SU2_MPI::GetComm(), &sendReturnReqs[i]);
  }

  /*--- Complete the first round of non-blocking sends. ---*/
  SU2_MPI::Waitall(sendReqs.size(), sendReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Loop over the ranks from which I receive the return information. ---*/
  for(unsigned long i=0; i<rankCommToInd.size(); ++i) {

    /*--- Block until a message arrives and determine the source of the message
          and its index in the original send buffers. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+3, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    MI = rankCommToInd.find(source);
    const int ind = MI->second;

    /*--- Allocate the memory for the receive buffer and receive the message using
          a blocking receive. ---*/
    vector<int> recvBuf(sendBuf[ind].size());
    SU2_MPI::Recv(recvBuf.data(), recvBuf.size(), MPI_INT,
                  source, rank+3, SU2_MPI::GetComm(), &status);

    /*--- Loop over the data just received and add them to the map
          mapMatrixIndToColor. ---*/
    for(unsigned long j=0; j<sendBuf[ind].size(); ++j)
      mapMatrixIndToColor[sendBuf[ind][j]] = recvBuf[j];
  }

  /*--- Complete the second round of non-blocking sends. ---*/
  SU2_MPI::Waitall(sendReturnReqs.size(), sendReturnReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Wild cards have been used in the communication,
        so synchronize the ranks to avoid problems. ---*/
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Part 3: For each locally owned DOF it is determined which matrix   ---*/
  /*---         entry is computed for which color. This information is     ---*/
  /*---         stored in the vector of vectors colorToIndEntriesJacobian. ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the ownded DOFs. ---*/
  colorToIndEntriesJacobian.resize(nDOFsLocOwned);
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {

    /*--- Initialize the elements of colorToIndEntriesJacobian to -1 to indicate
          that the color does not have a matrix entity . ---*/
    colorToIndEntriesJacobian[i].assign(nGlobalColors, -1);

    /*--- Loop over the matrix entries of this DOF. ---*/
    for(int j=0; j<(int) nonZeroEntriesJacobian[i].size(); ++j) {

      /*--- Search for the entry in mapMatrixIndToColor and determine the
            color of this entry. ---*/
      map<unsigned long, int>::const_iterator MMI;
      MMI = mapMatrixIndToColor.find(nonZeroEntriesJacobian[i][j]);
      if(MMI == mapMatrixIndToColor.end())
        SU2_MPI::Error(string("Matrix entry not found in mapMatrixIndToColor. Should not happen."),
                       CURRENT_FUNCTION);
      const int color = MMI->second;

      /*--- Set the correct information in colorToIndEntriesJacobian. ---*/
      colorToIndEntriesJacobian[i][color] = j;
    }
  }
}

void CFEM_DG_SolverBase::Prepare_MPI_Communication(const CMeshFEM_DG *DGGeometry,
                                                   CConfig           *config) {

  /*--- Get the communication information from DG_Geometry. Note that for a
        FEM DG discretization the communication entities of FEMGeometry contain
        the volume elements. ---*/
  const vector<int>                    &ranksSend    = DGGeometry->GetRanksSend();
  const vector<int>                    &ranksRecv    = DGGeometry->GetRanksRecv();
  const vector<vector<unsigned long> > &elementsSend = DGGeometry->GetEntitiesSend();
  const vector<vector<unsigned long> > &elementsRecv = DGGeometry->GetEntitiesRecv();                                                    

  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Create the triple vector, which contains elements to be    ---*/
  /*---         sent and received per time level. Note that when time      ---*/
  /*---         accurate local time stepping is not used, this is just a   ---*/
  /*---         copy of elementsSend and elementsRecv.                     ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Allocate the first two indices of the triple vectors. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

  vector<vector<vector<unsigned long> > > elemSendPerTimeLevel, elemRecvPerTimeLevel;
  elemSendPerTimeLevel.resize(nTimeLevels);
  elemRecvPerTimeLevel.resize(nTimeLevels);

  for(unsigned short i=0; i<nTimeLevels; ++i) {
    elemSendPerTimeLevel[i].resize(elementsSend.size());
    elemRecvPerTimeLevel[i].resize(elementsRecv.size());
  }

  /*--- Loop over the send data and set the corresponding data
        in elemSendPerTimeLevel. ---*/
  for(unsigned long i=0; i<elementsSend.size(); ++i) {
    for(unsigned long j=0; j<elementsSend[i].size(); ++j) {
      const unsigned long ii = elementsSend[i][j];
      elemSendPerTimeLevel[volElem[ii].timeLevel][i].push_back(ii);
    }
  }

  /*--- Loop over the receive data and set the corresponding data
        in recvSendPerTimeLevel. ---*/
  for(unsigned long i=0; i<elementsRecv.size(); ++i) {
    for(unsigned long j=0; j<elementsRecv[i].size(); ++j) {
      const unsigned long ii = elementsRecv[i][j];
      elemRecvPerTimeLevel[volElem[ii].timeLevel][i].push_back(ii);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Find out whether or not self communication is present and  ---*/
  /*---         set the corresponding data for elementsSendSelfComm and    ---*/
  /*---         elementsRecvSelfComm for all time levels. This can only be ---*/
  /*---         the case when periodic boundaries are present in the grid. ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Allocate the first index of elementsSendSelfComm and elementsRecvSelfComm. ---*/
  elementsSendSelfComm.resize(nTimeLevels);
  elementsRecvSelfComm.resize(nTimeLevels);

  /*--- Loop over the ranks to which data is sent and copy the elements for self
        communication for all time levels. ---*/
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    if(ranksSend[i] == rank) {
      for(unsigned short j=0; j<nTimeLevels; ++j)
        elementsSendSelfComm[j] = elemSendPerTimeLevel[j][i];
    }
  }

  /*--- Loop over the ranks from which data is received and copy the elements
        for self communication for all time levels. ---*/
  for(unsigned long i=0; i<ranksRecv.size(); ++i) {
    if(ranksRecv[i] == rank) {
      for(unsigned short j=0; j<nTimeLevels; ++j)
        elementsRecvSelfComm[j] = elemRecvPerTimeLevel[j][i];
    }
  }

#ifdef HAVE_MPI

  /*--------------------------------------------------------------------------*/
  /*--- Step 3. Determine the MPI communication patterns for               ---*/
  /*---         all time levels.                                           ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Determine the number of time DOFs that must be communicated.
        For non-ADER schemes this is simply set to 1. ---*/
  unsigned short nTimeDOFs = 1;
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG)
    nTimeDOFs = config->GetnTimeDOFsADER_DG();

  /*--- Determine the number of items per DOF that must be communicated. ---*/
  const unsigned short nItemsPerDOF = nTimeDOFs*nVar;

  /*--- Allocate the memory for the first index of the vectors that
        determine the MPI communication patterns. ---*/
  commRequests.resize(nTimeLevels);
  elementsRecvMPIComm.resize(nTimeLevels);
  elementsSendMPIComm.resize(nTimeLevels);
  ranksRecvMPI.resize(nTimeLevels);
  ranksSendMPI.resize(nTimeLevels);
  commRecvBuf.resize(nTimeLevels);
  commSendBuf.resize(nTimeLevels);

  /*--- Loop over the time levels. ---*/
  for(unsigned short level=0; level<nTimeLevels; ++level) {

    /*--- Determine the number of ranks from which data will be received
          for this time level. Self communication is excluded. ---*/
    int nRankRecv = 0;
    for(unsigned long i=0; i<ranksRecv.size(); ++i) {
      if((ranksRecv[i] != rank) && elemRecvPerTimeLevel[level][i].size()) ++nRankRecv;
    }

    /*--- Determine the number of ranks to which data will be send
          for this time level. Self communication is excluded. ---*/
    int nRankSend = 0;
    for(unsigned long i=0; i<ranksSend.size(); ++i) {
      if((ranksSend[i] != rank) && elemSendPerTimeLevel[level][i].size()) ++nRankSend;
    }

    /*--- Allocate the memory for the second index of the vectors that
          determine the MPI communication. ---*/
    commRequests[level].resize(nRankRecv+nRankSend);
    elementsRecvMPIComm[level].resize(nRankRecv);
    elementsSendMPIComm[level].resize(nRankSend);
    ranksRecvMPI[level].resize(nRankRecv);
    ranksSendMPI[level].resize(nRankSend);
    commRecvBuf[level].resize(nRankRecv);
    commSendBuf[level].resize(nRankSend);

    /*--- Determine the receive information. ---*/
    nRankRecv = 0;
    for(unsigned long i=0; i<ranksRecv.size(); ++i) {
      if((ranksRecv[i] != rank) && elemRecvPerTimeLevel[level][i].size()) {

        /*--- Copy the elements to be received and the rank from where they come. ---*/
        elementsRecvMPIComm[level][nRankRecv] = elemRecvPerTimeLevel[level][i];
        ranksRecvMPI[level][nRankRecv] = ranksRecv[i];

        /*--- Determine the size of the receive buffer and allocate the memory. ---*/
        unsigned long sizeBuf = 0;
        for(unsigned long j=0; j<elementsRecvMPIComm[level][nRankRecv].size(); ++j) {
          const unsigned long jj = elementsRecvMPIComm[level][nRankRecv][j];
          sizeBuf += volElem[jj].standardElemFlow->GetNDOFs();
        }

        sizeBuf *= nItemsPerDOF;
        commRecvBuf[level][nRankRecv].resize(sizeBuf);

        /*--- Update nRankRecv. ---*/
        ++nRankRecv;
      }
    }

    /*--- Determine the send information. ---*/
    nRankSend = 0;
    for(unsigned long i=0; i<ranksSend.size(); ++i) {
      if((ranksSend[i] != rank) && elemSendPerTimeLevel[level][i].size()) {

        /*--- Copy the elements to be sent and the rank to where they are sent. ---*/
        elementsSendMPIComm[level][nRankSend] = elemSendPerTimeLevel[level][i];
        ranksSendMPI[level][nRankSend] = ranksSend[i];

        /*--- Determine the size of the send buffer and allocate the memory. ---*/
        unsigned long sizeBuf = 0;
        for(unsigned long j=0; j<elementsSendMPIComm[level][nRankSend].size(); ++j) {
          const unsigned long jj = elementsSendMPIComm[level][nRankSend][j];
          sizeBuf += volElem[jj].standardElemFlow->GetNDOFs();
        }

        sizeBuf *= nItemsPerDOF;
        commSendBuf[level][nRankSend].resize(sizeBuf);

        /*--- Update nRankSend. ---*/
        ++nRankSend;
      }
    }
  }

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 4. Store the information for the rotational periodic          ---*/
  /*---         corrections for the halo elements.                         ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Get the data for the rotational periodic halos from FEMGeometry. ---*/
  const vector<unsigned short> &markersRotPer = DGGeometry->GetRotPerMarkers();
  const vector<vector<unsigned long> > &rotPerHalos = DGGeometry->GetRotPerHalos();

  /*--- Allocate the memory for the first and second index of
        halosRotationalPeriodicity. Also allocate the memory to store the
        rotation matrices of the periodic transformations. ---*/
  const unsigned short nRotPerMarkers = markersRotPer.size();
  halosRotationalPeriodicity.resize(nTimeLevels);
  for(unsigned short i=0; i<nTimeLevels; ++i)
    halosRotationalPeriodicity[i].resize(nRotPerMarkers);

  rotationMatricesPeriodicity.resize(9*nRotPerMarkers);

  /*--- Loop over the rotational periodic transformations. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<nRotPerMarkers; ++i) {

   /*--- Get the rotation angles from config for this marker. ---*/
   const unsigned short pInd = markersRotPer[i];
   auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(pInd));

   /*--- Determine the rotation matrix from the donor to the halo elements.
         This is the transpose of the rotation matrix from the halo to the
         donor elements, which is stored in periodic angles of the marker. ---*/
    const su2double theta = angles[0], phi = angles[1], psi = angles[2];

    const su2double cosTheta = cos(theta), cosPhi = cos(phi), cosPsi = cos(psi);
    const su2double sinTheta = sin(theta), sinPhi = sin(phi), sinPsi = sin(psi);

    rotationMatricesPeriodicity[ii++] =  cosPhi*cosPsi;
    rotationMatricesPeriodicity[ii++] =  cosPhi*sinPsi;
    rotationMatricesPeriodicity[ii++] = -sinPhi;

    rotationMatricesPeriodicity[ii++] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
    rotationMatricesPeriodicity[ii++] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
    rotationMatricesPeriodicity[ii++] = sinTheta*cosPhi;

    rotationMatricesPeriodicity[ii++] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
    rotationMatricesPeriodicity[ii++] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
    rotationMatricesPeriodicity[ii++] = cosTheta*cosPhi;

    /*--- Loop over the elements of this periodic transformation and store them
          in the appropriate location of halosRotationalPeriodicity. ---*/
    for(unsigned long j=0; j<rotPerHalos[i].size(); ++j) {
      const unsigned long jj = rotPerHalos[i][j];
      halosRotationalPeriodicity[volElem[jj].timeLevel][i].push_back(jj);
    }
  }
}

void CFEM_DG_SolverBase::DetermineCurrentPInPSequencing(CConfig *config) {

  if (config->GetUsePSequencing_DG() ) {

    const unsigned long nIterPSequencing = config->GetnIterPSequencing_DG();
    if (nIterPSequencing > 0) {
      const unsigned long currentIter = config->GetInnerIter();
      unsigned long current_p = currentIter/nIterPSequencing;
      if (current_p > 100) current_p = 100;
      currentPInPSequencing = static_cast<unsigned short>(current_p);
    }
  }
}
