/*!
 * \file CFEM_DG_SolverBase.cpp
 * \brief Main subroutines for the base class of the finite element flow solvers.
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 7.1.0 "Blackbird"
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
  SU2_MPI::Allreduce(&nDOFsLocOwned, &nDOFsGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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

  /*--- Determine the size of the work array (per variable).
        Assume a viscous simulation. ---*/
  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = max(sizeFluxes, nDOFsMax);

  unsigned int sizeGradSolInt = nIntegrationMax*nDim*nDOFsMax;

  sizeWorkArray = sizeFluxes + sizeGradSolInt + 4*nIntegrationMax;

  /*--- Check if the size suffices when ADER is used. ---*/
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {

    const unsigned int nTimeDOFs   = config->GetnTimeDOFsADER_DG();
    unsigned int sizePredictorADER = 4*nDOFsMax*nTimeDOFs + nDOFsMax;

    if(config->GetKind_ADER_Predictor() == ADER_ALIASED_PREDICTOR)
      sizePredictorADER += nDim*nDOFsMax + nDim*nDim*nIntegrationMax;
    else
      sizePredictorADER += (nDim+1)*max(nIntegrationMax, nDOFsMax)
                         + nDim*nDim*max(nIntegrationMax,nDOFsMax);

    sizeWorkArray = max(sizeWorkArray, sizePredictorADER);
  }

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

  /*--- Check if the symmetrizing terms are present. ---*/
  symmetrizingTermsPresent = false;
  if(config->GetViscous() && (fabs(config->GetTheta_Interior_Penalty_DGFEM()) > 1.e-8))
    symmetrizingTermsPresent = true;

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
  }

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

CFEM_DG_SolverBase::~CFEM_DG_SolverBase(void) {

  delete FluidModel;
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
                     MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
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

  /*--- Determine the global offset of the DOF numbering for the
        locally owned elements. Allocate it for all elements,
        because this information is communicated. ---*/
  vector<unsigned long> offsetDOFsElem(nVolElemTot+1);
  offsetDOFsElem[0] = nDOFsPerRank[rank];
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    offsetDOFsElem[i+1] = offsetDOFsElem[i] + volElem[i].standardElemFlow->GetNDOFs();

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
                   dest, dest, MPI_COMM_WORLD, &sendReqs[i]);
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
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /*--- Allocate the memory for the receive buffer and receive the data. ---*/
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank, MPI_COMM_WORLD, &status);

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
  SU2_MPI::Barrier(MPI_COMM_WORLD);

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
}
