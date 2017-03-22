/*!
 * \file solution_direct_mean_fem.cpp
 * \brief Main subroutines for solving finite element flow problems (Euler, Navier-Stokes, etc.).
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/solver_structure.hpp"

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(void) : CSolver() {

  /*--- Basic array initialization ---*/

  FluidModel = NULL;

  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;

  /*--- Surface-based array initialization ---*/
  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  Cauchy_Serie = NULL;

  /*--- Initialization of the boolean symmetrizingTermsPresent. ---*/
  symmetrizingTermsPresent = true;

  /*--- Initialize boolean for deallocating MPI comm. structure. ---*/
  mpiCommsPresent = false;

}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CConfig *config, unsigned short val_nDim, unsigned short iMesh) : CSolver() {

  /*--- Dummy solver constructor that calls the SetNondim. routine in
   order to load the flow non-dim. information into the config class.
   This is needed to complete a partitioning for time-accurate local
   time stepping that depends on the flow state. ---*/

  nDim = val_nDim;
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  SetNondimensionalization(config, iMesh, false);

  /*--- Basic array initialization ---*/

  FluidModel = NULL;

  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;

  /*--- Surface-based array initialization ---*/
  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  Cauchy_Serie = NULL;

  /*--- Initialization of the boolean symmetrizingTermsPresent. ---*/
  symmetrizingTermsPresent = true;

  /*--- Initialize boolean for deallocating MPI comm. structure. ---*/
  mpiCommsPresent = false;

}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

  /*--- Determine the restart information. ---*/
  const bool restart = (config->GetRestart() || config->GetRestart_Flow());
  string filename = config->GetSolution_FlowFileName();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Array initialization ---*/
  FluidModel = NULL;

  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL; CEff_Inv = NULL;
  CMx_Inv = NULL;   CMy_Inv = NULL;   CMz_Inv = NULL;
  CFx_Inv = NULL;   CFy_Inv = NULL;   CFz_Inv = NULL;

  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL;   Surface_CFy_Inv = NULL;   Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL;   Surface_CMy_Inv = NULL;   Surface_CMz_Inv = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL;   Surface_CFy = NULL;   Surface_CFz = NULL;
  Surface_CMx = NULL;   Surface_CMy = NULL;   Surface_CMz = NULL;

  Cauchy_Serie = NULL;

  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure. ---*/
  nDim    = geometry->GetnDim();
  nMarker = config->GetnMarker_All();

  const bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);

  if( compressible ) nVar = nDim + 2;
  else               nVar = nDim + 1;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  nVolElemTot   = DGGeometry->GetNVolElemTot();
  nVolElemOwned = DGGeometry->GetNVolElemOwned();
  volElem       = DGGeometry->GetVolElem();

  nVolElemOwnedPerTimeLevel    = DGGeometry->GetNVolElemOwnedPerTimeLevel();
  nVolElemInternalPerTimeLevel = DGGeometry->GetNVolElemInternalPerTimeLevel();

  nMeshPoints = DGGeometry->GetNMeshPoints();
  meshPoints  = DGGeometry->GetMeshPoints();

  nMatchingInternalFacesWithHaloElem = DGGeometry->GetNMatchingFacesWithHaloElem();
  nMatchingInternalFacesLocalElem    = DGGeometry->GetNMatchingFacesInternal();
  matchingInternalFaces              = DGGeometry->GetMatchingFaces();

  boundaries = DGGeometry->GetBoundaries();

  nStandardBoundaryFacesSol = DGGeometry->GetNStandardBoundaryFacesSol();
  nStandardElementsSol      = DGGeometry->GetNStandardElementsSol();
  nStandardMatchingFacesSol = DGGeometry->GetNStandardMatchingFacesSol();

  standardBoundaryFacesSol = DGGeometry->GetStandardBoundaryFacesSol();
  standardElementsSol      = DGGeometry->GetStandardElementsSol();
  standardMatchingFacesSol = DGGeometry->GetStandardMatchingFacesSol();

  LagrangianBeginTimeIntervalADER_DG  = DGGeometry->GetLagrangianBeginTimeIntervalADER_DG();
  timeInterpolDOFToIntegrationADER_DG = DGGeometry->GetTimeInterpolDOFToIntegrationADER_DG();

  /*--- Determine the maximum number of integration points used. Usually this
        is for the volume integral, but to avoid problems the faces are also
        taken into account.  ---*/
  nIntegrationMax = 0;
  for(unsigned short i=0; i<nStandardBoundaryFacesSol; ++i) {
    const unsigned short nInt = standardBoundaryFacesSol[i].GetNIntegration();
    nIntegrationMax = max(nIntegrationMax, nInt);
  }

  for(unsigned short i=0; i<nStandardElementsSol; ++i) {
    const unsigned short nInt = standardElementsSol[i].GetNIntegration();
    nIntegrationMax = max(nIntegrationMax, nInt);
  }

  for(unsigned short i=0; i<nStandardMatchingFacesSol; ++i) {
    const unsigned short nInt = standardMatchingFacesSol[i].GetNIntegration();
    nIntegrationMax = max(nIntegrationMax, nInt);
  }

  /*--- Determine the maximum number of DOFs used. This is for the volume elements.
        Note that also the element adjacent to side 1 of the matching faces must
        be taken into account, because this could be an external element. */
  nDOFsMax = 0;
  for(unsigned short i=0; i<nStandardElementsSol; ++i) {
    const unsigned short nDOFs = standardElementsSol[i].GetNDOFs();
    nDOFsMax = max(nDOFsMax, nDOFs);
  }

  for(unsigned short i=0; i<nStandardMatchingFacesSol; ++i) {
    const unsigned short nDOFs = standardMatchingFacesSol[i].GetNDOFsElemSide1();
    nDOFsMax = max(nDOFsMax, nDOFs);
  }

  /*--- Make sure that nIntegrationMax and nDOFsMax are such that the allocations
        of temporary memory is 64 byte aligned.    ---*/
  if( nIntegrationMax%8 ) nIntegrationMax += 8 - nIntegrationMax%8;
  if( nDOFsMax%8 )        nDOFsMax        += 8 - nDOFsMax%8;

  /*--- Determine the size of the vector VecTmpMemory and allocate its memory. ---*/
  unsigned int sizeVecTmp;
  if( config->GetViscous() ) {

    /* Viscous simulation. */
    unsigned int sizeFluxes = nIntegrationMax*nDim;
    sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

    const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

    sizeVecTmp = 2*nIntegrationMax*(2 + nVar) + sizeFluxes + sizeGradSolInt
               + max(nIntegrationMax,nDOFsMax)*nVar;
  }
  else {

    /* Inviscid simulation. */
    unsigned int sizeVol = nVar*nIntegrationMax*(nDim+1);
    unsigned int sizeSur = nVar*(2*nIntegrationMax + max(nIntegrationMax,nDOFsMax));

    sizeVecTmp = max(sizeVol, sizeSur);
  }

  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {

    /*--- ADER-DG scheme. Determine the size needed for the predictor step
          and make sure that sizeVecTmp is big enough. This size depends
          whether an aliased or a non-aliased predictor step is used. Note
          that the size estimates are for viscous computations. ---*/
    const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();

    unsigned int sizePredictorADER = 4*nVar*nDOFsMax*nTimeDOFs
                                   +   nVar*nDOFsMax;

    if(config->GetKind_ADER_Predictor() == ADER_ALIASED_PREDICTOR)
      sizePredictorADER += nDim*nVar*nDOFsMax + nDim*nDim*nVar*nIntegrationMax;
    else
      sizePredictorADER += (nDim+1)*nVar*max(nIntegrationMax, nDOFsMax)
                         + nDim*nDim*nVar*max(nIntegrationMax,nDOFsMax);

    sizeVecTmp = max(sizeVecTmp, sizePredictorADER);
  }

  VecTmpMemory.resize(sizeVecTmp);

  /*--- Perform the non-dimensionalization for the flow equations using the
        specified reference values. ---*/
  SetNondimensionalization(config, iMesh, true);

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual_RMS = new su2double[nVar];     for(unsigned short iVar=0; iVar<nVar; ++iVar) Residual_RMS[iVar] = 0.0;
  Residual_Max = new su2double[nVar];     for(unsigned short iVar=0; iVar<nVar; ++iVar) Residual_Max[iVar] = 0.0;
  Point_Max    = new unsigned long[nVar]; for(unsigned short iVar=0; iVar<nVar; ++iVar) Point_Max[iVar]    = 0;

  Point_Max_Coord = new su2double*[nVar];
  for (unsigned short iVar=0; iVar<nVar; ++iVar) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for(unsigned short iDim=0; iDim<nDim; ++iDim) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Non-dimensional coefficients ---*/
  CD_Inv   = new su2double[nMarker];
  CL_Inv   = new su2double[nMarker];
  CSF_Inv  = new su2double[nMarker];
  CMx_Inv  = new su2double[nMarker];
  CMy_Inv  = new su2double[nMarker];
  CMz_Inv  = new su2double[nMarker];
  CEff_Inv = new su2double[nMarker];
  CFx_Inv  = new su2double[nMarker];
  CFy_Inv  = new su2double[nMarker];
  CFz_Inv  = new su2double[nMarker];

  Surface_CL_Inv   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Inv   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CL       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz      = new su2double[config->GetnMarker_Monitoring()];

  /*--- Init total coefficients ---*/
  Total_CD   = 0.0; Total_CL  = 0.0; Total_CSF = 0.0;
  Total_CMx  = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
  Total_CFx  = 0.0; Total_CFy = 0.0; Total_CFz = 0.0;
  Total_CEff = 0.0;

  /*--- Read farfield conditions ---*/
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Mach_Inf        = config->GetMach();

  /*--- Set the conservative variables of the free-stream. ---*/
  ConsVarFreeStream.resize(nVar);
  if( compressible ) {
    ConsVarFreeStream[0] = Density_Inf;
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      ConsVarFreeStream[iDim+1] = Density_Inf*Velocity_Inf[iDim];
    ConsVarFreeStream[nVar-1] = Density_Inf*Energy_Inf;
  }
  else {
    ConsVarFreeStream[0] = Pressure_Inf;
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      ConsVarFreeStream[iDim+1] = Density_Inf*Velocity_Inf[iDim];
  }

  /*--- Determine the total number of DOFs stored on this rank and allocate the memory
        to store the conservative variables. ---*/
  nDOFsLocOwned = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) nDOFsLocOwned += volElem[i].nDOFsSol;

  nDOFsLocTot = nDOFsLocOwned;
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i) nDOFsLocTot += volElem[i].nDOFsSol;

  VecSolDOFs.resize(nVar*nDOFsLocTot);
  VecSolDOFsOld.resize(nVar*nDOFsLocOwned);

  /*--- Check for the ADER-DG time integration scheme and allocate the memory
        for the additional vectors. ---*/
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {

    const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();

    VecTotResDOFsADER.resize(nVar*nDOFsLocOwned);
    VecSolDOFsPredictorADER.resize(nTimeDOFs*nVar*nDOFsLocOwned);
  }
  else {

    /*--- Runge Kutta type of time integration schemes. Allocate the memory to
          possibly store the new solution. ---*/
    if(config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT)
      VecSolDOFsNew.resize(nVar*nDOFsLocOwned);
  }

  /*--- Determine the global number of DOFs. ---*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nDOFsLocOwned, &nDOFsGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nDOFsGlobal = nDOFsLocOwned;
#endif

  /*--- Allocate the memory to store the time steps, residuals, etc. ---*/
  VecDeltaTime.resize(nVolElemOwned);
  VecResDOFs.resize(nVar*nDOFsLocTot);
  nEntriesResFaces.assign(nDOFsLocTot+1, 0);
  startLocResFacesMarkers.resize(nMarker);

  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

  startLocResInternalFacesLocalElem.assign(nTimeLevels+1, 0);
  startLocResInternalFacesWithHaloElem.assign(nTimeLevels+1, 0);

  /*--- Determine the size of the vector to store residuals that come from the
        integral over the faces and determine the number of entries in this
        vector for the local DOFs. ---*/
  symmetrizingTermsPresent = false;
  if(config->GetViscous() && (fabs(config->GetTheta_Interior_Penalty_DGFEM()) > 1.e-8))
    symmetrizingTermsPresent = true;

  /*--- First the internal matching faces. ---*/
  unsigned long sizeVecResFaces = 0;
  for(unsigned long i=0; i<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++i) {

    /* The terms that only contribute to the DOFs located on the face. */
    const unsigned short ind = matchingInternalFaces[i].indStandardElement;
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    sizeVecResFaces += nDOFsFace0;
    for(unsigned short j=0; j<nDOFsFace0; ++j)
      ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolFaceSide0[j]+1];

    sizeVecResFaces += nDOFsFace1;
    for(unsigned short j=0; j<nDOFsFace1; ++j)
      ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolFaceSide1[j]+1];

    /* The symmetrizing terms, if present, contribute to all
       the DOFs of the adjacent elements. */
    if( symmetrizingTermsPresent ) {
      const unsigned short nDOFsElem0 = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
      const unsigned short nDOFsElem1 = standardMatchingFacesSol[ind].GetNDOFsElemSide1();

      sizeVecResFaces += nDOFsElem0;
      for(unsigned short j=0; j<nDOFsElem0; ++j)
        ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolElementSide0[j]+1];

      sizeVecResFaces += nDOFsElem1;
      for(unsigned short j=0; j<nDOFsElem1; ++j)
        ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolElementSide1[j]+1];
    }

    /* Determine the time level of this face and store the position of the
       residual in the appropriate entry. */
    const unsigned long  elem0     = matchingInternalFaces[i].elemID0;
    const unsigned long  elem1     = matchingInternalFaces[i].elemID1;
    const unsigned short timeLevel = min(volElem[elem0].timeLevel,
                                         volElem[elem1].timeLevel);

    if(i < nMatchingInternalFacesWithHaloElem[0] )
      startLocResInternalFacesLocalElem[timeLevel+1] = sizeVecResFaces;
    else
      startLocResInternalFacesWithHaloElem[timeLevel+1] = sizeVecResFaces;
  }

  /* Set the uninitialized values of startLocResInternalFacesLocalElem. */
  for(unsigned short i=1; i<=nTimeLevels; ++i) {
    if(startLocResInternalFacesLocalElem[i] == 0)
      startLocResInternalFacesLocalElem[i] = startLocResInternalFacesLocalElem[i-1];
  }

  /* Set the uninitialized values of startLocResInternalFacesWithHaloElem. */
  startLocResInternalFacesWithHaloElem[0] = startLocResInternalFacesLocalElem[nTimeLevels];

  for(unsigned short i=1; i<=nTimeLevels; ++i) {
    if(startLocResInternalFacesWithHaloElem[i] == 0)
      startLocResInternalFacesWithHaloElem[i] = startLocResInternalFacesWithHaloElem[i-1];
  }

  /* The physical boundary faces. Exclude the periodic boundaries,
     because these are not physical boundaries. */
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    startLocResFacesMarkers[iMarker].assign(nTimeLevels+1, 0);
    startLocResFacesMarkers[iMarker][0] = sizeVecResFaces;

    if( !(boundaries[iMarker].periodicBoundary) ) {

      /* Easier storage of the variables for this boundary. */
      const unsigned long      nSurfElem = boundaries[iMarker].surfElem.size();
      const CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

      /*--- Loop over the surface elements and update the required data. ---*/
      for(unsigned long i=0; i<nSurfElem; ++i) {
        const unsigned short ind       = surfElem[i].indStandardElement;
        const unsigned short nDOFsFace = standardBoundaryFacesSol[ind].GetNDOFsFace();

        /* The terms that only contribute to the DOFs located on the face. */
        sizeVecResFaces += nDOFsFace;
        for(unsigned short j=0; j<nDOFsFace; ++j)
          ++nEntriesResFaces[surfElem[i].DOFsSolFace[j]+1];

        /* The symmetrizing terms, if present, contribute to all
           the DOFs of the adjacent elements. */
        if( symmetrizingTermsPresent ) {
          const unsigned short nDOFsElem = standardBoundaryFacesSol[ind].GetNDOFsElem();

          sizeVecResFaces += nDOFsElem;
          for(unsigned short j=0; j<nDOFsElem; ++j)
            ++nEntriesResFaces[surfElem[i].DOFsSolElement[j]+1];
        }

        /* Determine the time level of the adjacent element and store the position of the
           residual in the appropriate entry. */
        const unsigned short timeLevel = volElem[surfElem[i].volElemID].timeLevel;
        startLocResFacesMarkers[iMarker][timeLevel+1] = sizeVecResFaces;
      }
    }

    /* Set the unitialized values of startLocResFacesMarkers[iMarker]. */
    for(unsigned short i=1; i<=nTimeLevels; ++i) {
      if(startLocResFacesMarkers[iMarker][i] == 0)
        startLocResFacesMarkers[iMarker][i] = startLocResFacesMarkers[iMarker][i-1];
    }
  }

  /*--- Put nEntriesResFaces in cumulative storage format and allocate the
        memory for entriesResFaces and VecResFaces. ---*/
  for(unsigned long i=0; i<nDOFsLocTot; ++i)
    nEntriesResFaces[i+1] += nEntriesResFaces[i];

  entriesResFaces.resize(nEntriesResFaces[nDOFsLocTot]);
  VecResFaces.resize(nVar*sizeVecResFaces);

  /*--- Repeat the loops over the internal and boundary faces, but now store
        the enties in entriesResFaces. A counter variable is needed to keep
        track of the appropriate location in entriesResFaces. ---*/
  vector<unsigned long> counterEntries = nEntriesResFaces;

  /* First the loop over the internal matching faces. */
  sizeVecResFaces = 0;
  for(unsigned long i=0; i<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++i) {

    /* The terms that only contribute to the DOFs located on the face. */
    const unsigned short ind = matchingInternalFaces[i].indStandardElement;
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    for(unsigned short j=0; j<nDOFsFace0; ++j) {
      unsigned long jj    = counterEntries[matchingInternalFaces[i].DOFsSolFaceSide0[j]]++;
      entriesResFaces[jj] = sizeVecResFaces++;
    }

    for(unsigned short j=0; j<nDOFsFace1; ++j) {
      unsigned long jj    = counterEntries[matchingInternalFaces[i].DOFsSolFaceSide1[j]]++;
      entriesResFaces[jj] = sizeVecResFaces++;
    }

    /* The symmetrizing terms, if present, contribute to all
       the DOFs of the adjacent elements. */
    if( symmetrizingTermsPresent ) {
      const unsigned short nDOFsElem0 = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
      const unsigned short nDOFsElem1 = standardMatchingFacesSol[ind].GetNDOFsElemSide1();

      for(unsigned short j=0; j<nDOFsElem0; ++j) {
        unsigned long jj    = counterEntries[matchingInternalFaces[i].DOFsSolElementSide0[j]]++;
        entriesResFaces[jj] = sizeVecResFaces++;
      }

      for(unsigned short j=0; j<nDOFsElem1; ++j) {
        unsigned long jj    = counterEntries[matchingInternalFaces[i].DOFsSolElementSide1[j]]++;
        entriesResFaces[jj] = sizeVecResFaces++;
      }
    }
  }

  /* And the physical boundary faces. Exclude the periodic boundaries,
     because these are not physical boundaries. */
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    if( !(boundaries[iMarker].periodicBoundary) ) {

      /* Easier storage of the variables for this boundary. */
      const unsigned long      nSurfElem = boundaries[iMarker].surfElem.size();
      const CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

      /*--- Loop over the surface elements to set entriesResFaces. ---*/
      for(unsigned long i=0; i<nSurfElem; ++i) {

        /* The terms that only contribute to the DOFs located on the face. */
        const unsigned short ind       = surfElem[i].indStandardElement;
        const unsigned short nDOFsFace = standardBoundaryFacesSol[ind].GetNDOFsFace();

        for(unsigned short j=0; j<nDOFsFace; ++j) {
          unsigned long jj    = counterEntries[surfElem[i].DOFsSolFace[j]]++;
          entriesResFaces[jj] = sizeVecResFaces++;
        }

        /* The symmetrizing terms, if present, contribute to all
           the DOFs of the adjacent elements. */
        if( symmetrizingTermsPresent ) {
          const unsigned short nDOFsElem = standardBoundaryFacesSol[ind].GetNDOFsElem();

          for(unsigned short j=0; j<nDOFsElem; ++j) {
            unsigned long jj    = counterEntries[surfElem[i].DOFsSolElement[j]]++;
            entriesResFaces[jj] = sizeVecResFaces++;
          }
        }
      }
    }
  }

  /*--- Check for a restart and set up the variables at each node
        appropriately. Coarse multigrid levels will be intitially set to
        the farfield values bc the solver will immediately interpolate
        the solution from the finest mesh to the coarser levels. ---*/
  if (!restart || (iMesh != MESH_0)) {

    /*--- Start the solution from the free-stream state ---*/

    unsigned long ii = 0;
    for(unsigned long i=0; i<nDOFsLocTot; ++i) {
      for(unsigned short j=0; j<nVar; ++j, ++ii) {
        VecSolDOFs[ii] = ConsVarFreeStream[j];
      }
    }

  } else {

    /*--- Open the restart file, throw an error if this fails. ---*/
    ifstream restart_file;
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file " << filename.data() << "!!!"<< endl;
      exit(EXIT_FAILURE);
    }

    /*--- Create the map from the global DOF ID to the local index. ---*/
    map<unsigned long, unsigned long> mapGlobal2Local;

    unsigned long ii = 0;
    for(unsigned long i=0; i<nVolElemOwned; ++i) {
      for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j, ++ii) {
        mapGlobal2Local[volElem[i].offsetDOFsSolGlobal+j] = ii;
      }
    }

    /*--- The first line is the header ---*/
    string text_line;
    getline (restart_file, text_line);

    /*--- Read all lines in the restart file ---*/
    unsigned long iDOF_Global = 0, nDOF_Read = 0;
    while (getline (restart_file, text_line)) {
      istringstream point_line(text_line);

      /*--- Check if this DOF must be stored on this rank. ---*/
      map<unsigned long, unsigned long>::const_iterator MI;
      MI = mapGlobal2Local.find(iDOF_Global);
      if(MI != mapGlobal2Local.end()) {

        /*--- This DOF must be stored on this rank. Retrieve the local index
              and read the data from file. ---*/
        const unsigned long iDOF_Local = nVar*MI->second;

        unsigned long index;
        point_line >> index;
        for(unsigned short i=0; i<nDim; ++i) {
          su2double dull_val;
          point_line >> dull_val;
        }

        for(unsigned short i=0; i<nVar; ++i) {
          point_line >> VecSolDOFs[iDOF_Local+i];
        }

        /*--- Update the local counter nDOF_Read. ---*/
        ++nDOF_Read;
      }

      /*--- Update the counter iDOF_Global. ---*/
      ++iDOF_Global;
    }

    /*--- Detect a wrong solution file ---*/
    unsigned short rbuf_NotMatching = 0;
    if(nDOF_Read < nDOFsLocOwned) rbuf_NotMatching = 1;

#ifdef HAVE_MPI
    unsigned short sbuf_NotMatching = rbuf_NotMatching;
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, MPI_COMM_WORLD);
#endif

    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
        cout << "It could be empty lines at the end of the file." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }

    /*--- Close the restart file ---*/
    restart_file.close();
  }

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  unsigned long nBadDOFs = 0;

  if( compressible ) {

    for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
      const unsigned long ii = nVar*i;
      su2double DensityInv = 1.0/VecSolDOFs[ii];

      su2double Velocity2 = 0.0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double vel = VecSolDOFs[ii+iDim]*DensityInv;
        Velocity2 += vel*vel;
      }

      su2double StaticEnergy = VecSolDOFs[ii+nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(VecSolDOFs[ii], StaticEnergy);
      su2double Pressure = FluidModel->GetPressure();
      su2double Temperature = FluidModel->GetTemperature();

      /*--- Use the values at the infinity if the state is not physical. ---*/
      if((Pressure < 0.0) || (VecSolDOFs[ii] < 0.0) || (Temperature < 0.0)) {
        for(unsigned short j=0; j<nVar; ++j) {
          VecSolDOFs[ii+j] = ConsVarFreeStream[j];
        }

        ++nBadDOFs;
      }
    }
  }

  /*--- Warning message about non-physical points ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    unsigned long nBadDOFsLoc = nBadDOFs;
    SU2_MPI::Reduce(&nBadDOFsLoc, &nBadDOFs, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#endif

    if((rank == MASTER_NODE) && (nBadDOFs != 0))
      cout << "Warning. The initial solution contains "<< nBadDOFs << " DOFs that are not physical." << endl;
  }

  /*--- Set up the persistent communication for the conservative variables and
        the reverse communication for the residuals of the halo elements. ---*/
  Prepare_MPI_Communication(DGGeometry, config);

  /*--- Initialize boolean for deallocating MPI comm. structure. ---*/
  mpiCommsPresent = true;

  /*--- Perform the MPI communication of the solution. ---*/
  Initiate_MPI_Communication();
  Complete_MPI_Communication();
}

CFEM_DG_EulerSolver::~CFEM_DG_EulerSolver(void) {

  if(FluidModel != NULL) delete FluidModel;

  /*--- Array deallocation ---*/
  if (CD_Inv != NULL)           delete [] CD_Inv;
  if (CL_Inv != NULL)           delete [] CL_Inv;
  if (CSF_Inv != NULL)          delete [] CSF_Inv;
  if (CMx_Inv != NULL)          delete [] CMx_Inv;
  if (CMy_Inv != NULL)          delete [] CMy_Inv;
  if (CMz_Inv != NULL)          delete [] CMz_Inv;
  if (CFx_Inv != NULL)          delete [] CFx_Inv;
  if (CFy_Inv != NULL)          delete [] CFy_Inv;
  if (CFz_Inv != NULL)          delete [] CFz_Inv;
  if (Surface_CL_Inv != NULL)   delete [] Surface_CL_Inv;
  if (Surface_CD_Inv != NULL)   delete [] Surface_CD_Inv;
  if (Surface_CSF_Inv != NULL)  delete [] Surface_CSF_Inv;
  if (Surface_CEff_Inv != NULL) delete [] Surface_CEff_Inv;
  if (Surface_CFx_Inv != NULL)  delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv != NULL)  delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv != NULL)  delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;
  if (Surface_CL != NULL)       delete [] Surface_CL;
  if (Surface_CD != NULL)       delete [] Surface_CD;
  if (Surface_CSF != NULL)      delete [] Surface_CSF;
  if (Surface_CEff != NULL)     delete [] Surface_CEff;
  if (Surface_CFx != NULL)      delete [] Surface_CFx;
  if (Surface_CFy != NULL)      delete [] Surface_CFy;
  if (Surface_CFz != NULL)      delete [] Surface_CFz;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)         delete [] CEff_Inv;

  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;

#ifdef HAVE_MPI

  if (mpiCommsPresent) {

    /*--- Release the memory of the persistent communication and the derived
          data types. ---*/
    for(int i=0; i<nCommRequests; ++i) MPI_Request_free(&commRequests[i]);
    for(int i=0; i<nCommRequests; ++i) MPI_Type_free(&commTypes[i]);

    for(int i=0; i<nCommRequests; ++i) MPI_Request_free(&reverseCommRequests[i]);
    for(unsigned int i=0; i<reverseCommTypes.size(); ++i)
       MPI_Type_free(&reverseCommTypes[i]);
  }

#endif
}

void CFEM_DG_EulerSolver::SetNondimensionalization(CConfig        *config,
                                                   unsigned short iMesh,
                                                   const bool     writeOutput) {

  su2double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0,
  Energy_FreeStream = 0.0, ModVel_FreeStreamND = 0.0, Velocity_Reynolds = 0.0,
  Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
  Temperature_Ref = 0.0, Time_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Energy_Ref= 0.0,
  Froude = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;

  unsigned short iDim;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Local variables ---*/

  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double Reynolds         = config->GetReynolds();
  bool unsteady           = (config->GetUnsteady_Simulation() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool turbulent          = (config->GetKind_Solver() == FEM_RANS) || (config->GetKind_Solver() == FEM_LES);
  bool tkeNeeded          = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool standard_air       = (config->GetKind_FluidModel() == STANDARD_AIR);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);

      FluidModel = new CIdealGas(1.4, config->GetGas_Constant());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case IDEAL_GAS:

      FluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case VW_GAS:

      FluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                       config->GetPressure_Critical(), config->GetTemperature_Critical());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case PR_GAS:

      FluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                     config->GetTemperature_Critical(), config->GetAcentric_Factor());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

  }

  Mach2Vel_FreeStream = FluidModel->GetSoundSpeed();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }

  /*--- Compute the modulus of the free stream velocity ---*/

  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Viscous initialization ---*/

  if (viscous) {

    /*--- Reynolds based initialization ---*/

    if (reynolds_init) {

      /*--- First, check if there is mesh motion. If yes, use the Mach
       number relative to the body to initialize the flow. ---*/

      if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
      else Velocity_Reynolds = ModVel_FreeStream;

      /*--- Change of measurement system, hard coded value working only with STANDAR AIR model ---*/

      if (standard_air) {
        if (config->GetSystemMeasurements() == SI) {
          config->SetMu_RefND(1.716E-5);
          config->SetMu_SND(110.4);
          config->SetMu_Temperature_RefND(273.15);
        }
        if (config->GetSystemMeasurements() == US) {
          config->SetMu_RefND(3.62E-7);
          config->SetMu_SND(198.72);
          config->SetMu_Temperature_RefND(518.7);
        }
      }

      /*--- For viscous flows, pressure will be computed from a density
       that is found from the Reynolds number. The viscosity is computed
       from the dimensional version of Sutherland's law ---*/

      FluidModel->SetLaminarViscosityModel(config);

      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);

      Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
      config->SetDensity_FreeStream(Density_FreeStream);
      FluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
      Pressure_FreeStream = FluidModel->GetPressure();
      config->SetPressure_FreeStream(Pressure_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Thermodynamics quantities based initialization ---*/

    else {

      FluidModel->SetLaminarViscosityModel(config);
      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Turbulence kinetic energy ---*/

    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }
  else {

    /*--- For inviscid flow, energy is calculated from the specified
     FreeStream quantities using the proper gas law. ---*/

    Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

  }

  /*-- Compute the freestream energy. ---*/

  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute non dimensional quantities. By definition,
   Lref is one because we have converted the grid to meters. ---*/

  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
    Density_Ref       = 1.0;
    Temperature_Ref   = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
    Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;                        config->SetForce_Ref(Force_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);         config->SetFroude(Froude);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/

  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);

  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);


  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);

  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);

  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/(Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStream(Omega_FreeStream);

  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/(Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);

  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/

  /*--- Delete the original (dimensional) FluidModel object before replacing. ---*/

  delete FluidModel;

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:
      FluidModel = new CIdealGas(1.4, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case IDEAL_GAS:
      FluidModel = new CIdealGas(Gamma, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case VW_GAS:
      FluidModel = new CVanDerWaalsGas(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                       config->GetTemperature_Critical()/config->GetTemperature_Ref());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case PR_GAS:
      FluidModel = new CPengRobinson(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                     config->GetTemperature_Critical()/config->GetTemperature_Ref(), config->GetAcentric_Factor());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

  }

  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

  if (viscous) {

    /*--- Constant viscosity model ---*/
    config->SetMu_ConstantND(config->GetMu_ConstantND()/Viscosity_Ref);

    /*--- Sutherland's model ---*/

    config->SetMu_RefND(config->GetMu_RefND()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_SND()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_RefND()/config->GetTemperature_Ref());

    /* constant thermal conductivity model */
    config->SetKt_ConstantND(config->GetKt_ConstantND()/Conductivity_Ref);

    FluidModel->SetLaminarViscosityModel(config);
    FluidModel->SetThermalConductivityModel(config);

  }

  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is required and if this is the master node and first domain ---*/

  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && writeOutput) {

    cout.precision(6);

    if (viscous) {
      cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
      cout << "based on the free-stream temperature and a density computed" << endl;
      cout << "from the Reynolds number." << endl;
    } else {
      cout << "Inviscid flow: Computing density based on free-stream" << endl;
      cout << "temperature and pressure using the ideal gas law." << endl;
    }

    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;

    cout <<"-- Input conditions:"<< endl;

    switch (config->GetKind_FluidModel()) {

      case STANDARD_AIR:
        cout << "Fluid Model: STANDARD_AIR "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant();
        if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;

      case IDEAL_GAS:
        cout << "Fluid Model: IDEAL_GAS "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;

      case VW_GAS:
        cout << "Fluid Model: Van der Waals "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;

      case PR_GAS:
        cout << "Fluid Model: Peng-Robinson "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;

    }
    if (viscous) {
      switch (config->GetKind_ViscosityModel()) {

        case CONSTANT_VISCOSITY:
          cout << "Viscosity Model: CONSTANT_VISCOSITY  "<< endl;
          cout << "Laminar Viscosity: " << config->GetMu_ConstantND()*Viscosity_Ref;
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          break;

        case SUTHERLAND:
          cout << "Viscosity Model: SUTHERLAND "<< endl;
          cout << "Ref. Laminar Viscosity: " << config->GetMu_RefND()*Viscosity_Ref;
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Ref. Temperature: " << config->GetMu_Temperature_RefND()*config->GetTemperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Sutherland Constant: "<< config->GetMu_SND()*config->GetTemperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          cout << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND()<< endl;
          cout << "Sutherland constant (non-dim): "<< config->GetMu_SND()<< endl;
          break;

      }
      switch (config->GetKind_ConductivityModel()) {

        case CONSTANT_PRANDTL:
          cout << "Conductivity Model: CONSTANT_PRANDTL  "<< endl;
          cout << "Prandtl: " << config->GetPrandtl_Lam()<< endl;
          break;

        case CONSTANT_CONDUCTIVITY:
          cout << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< endl;
          cout << "Molecular Conductivity: " << config->GetKt_ConstantND()*Conductivity_Ref<< " W/m^2.K." << endl;
          cout << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< endl;
          break;

      }
    }

    cout << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Free-stream total pressure: " << config->GetPressure_FreeStream() * pow( 1.0+Mach*Mach*0.5*(Gamma-1.0), Gamma/(Gamma-1.0) );
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Free-stream temperature: " << config->GetTemperature_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;

    cout << "Free-stream density: " << config->GetDensity_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;

    if (nDim == 2) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
    }
    if (config->GetSystemMeasurements() == SI) cout << " m/s. ";
    else if (config->GetSystemMeasurements() == US) cout << " ft/s. ";

    cout << "Magnitude: "	<< config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;

    cout << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;

    if (viscous) {
      cout << "Free-stream viscosity: " << config->GetViscosity_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
        cout << "Free-stream specific dissipation: " << config->GetOmega_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " 1/s." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " 1/s." << endl;
      }
    }

    if (unsteady) { cout << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << endl; }

    /*--- Print out reference values. ---*/

    cout <<"-- Reference values:"<< endl;

    cout << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;

    cout << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Reference temperature: " << config->GetTemperature_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;

    cout << "Reference density: " << config->GetDensity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;

    cout << "Reference velocity: " << config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;

    cout << "Reference energy per unit mass: " << config->GetEnergy_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;

    if (viscous) {
      cout << "Reference viscosity: " << config->GetViscosity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      cout << "Reference conductivity: " << config->GetConductivity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " W/m^2.K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf/ft.s.R." << endl;
    }


    if (unsteady) cout << "Reference time: " << config->GetTime_Ref() <<" s." << endl;

    /*--- Print out resulting non-dim values here. ---*/

    cout << "-- Resulting non-dimensional state:" << endl;
    cout << "Mach number (non-dim): " << config->GetMach() << endl;
    if (viscous) {
      cout << "Reynolds number (non-dim): " << config->GetReynolds() <<". Re length: " << config->GetLength_Reynolds();
      if (config->GetSystemMeasurements() == SI) cout << " m." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " ft." << endl;
    }

    cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
    cout << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << endl;

    cout << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << endl;

    cout << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << endl;

    if (nDim == 2) {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    cout << "Magnitude: "	 << config->GetModVel_FreeStreamND() << endl;

    cout << "Free-stream total energy per unit mass (non-dim): " << config->GetEnergy_FreeStreamND() << endl;

    if (viscous) {
      cout << "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << endl;
        cout << "Free-stream specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << endl;
      }
    }

    if (unsteady) {
      cout << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << endl;
      cout << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << endl;
    }

    cout << endl;

  }
}

void CFEM_DG_EulerSolver::Prepare_MPI_Communication(const CMeshFEM *FEMGeometry,
                                                    CConfig        *config) {

  /*--- Get the communication information from DG_Geometry. Note that for a
        FEM DG discretization the communication entities of FEMGeometry contain
        the volume elements. ---*/
  const vector<int>                    &ranksSend    = FEMGeometry->GetRanksSend();
  const vector<int>                    &ranksRecv    = FEMGeometry->GetRanksRecv();
  const vector<vector<unsigned long> > &elementsSend = FEMGeometry->GetEntitiesSend();
  const vector<vector<unsigned long> > &elementsRecv = FEMGeometry->GetEntitiesRecv();

  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Find out whether or not self communication is present.     ---*/
  /*---         This can only be the case when periodic boundaries are     ---*/
  /*---         present in the grid.                                       ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the rank inside MPI_COMM_WORLD. */
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* Store the send information of the self communication in
     elementsSendSelfComm, if present. */
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    if(ranksSend[i] == rank)
      elementsSendSelfComm = elementsSend[i];
  }

  /* Store the receive information of the self communication in
     elementsRecvSelfComm, if present. */
  for(unsigned long i=0; i<ranksRecv.size(); ++i) {
    if(ranksRecv[i] == rank)
      elementsRecvSelfComm = elementsRecv[i];
  }

#ifdef HAVE_MPI

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Set up the persistent MPI communication for the halo data. ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Determine the number of communication requests that take place
        for the exchange of the halo data. These requests are both send and
        receive requests. Allocate the necessary memory for the communication
        of the halo data. ---*/
  nCommRequests = ranksSend.size() + ranksRecv.size();
  if( elementsRecvSelfComm.size() ) --nCommRequests;
  if( elementsSendSelfComm.size() ) --nCommRequests;

  commRequests.resize(nCommRequests);
  commTypes.resize(nCommRequests);

  /*--- Loop over the ranks to which this rank has to send halo data and create
        the send requests. Exclude self communication. ---*/
  unsigned int nn = 0;
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    if(ranksSend[i] != rank) {

      /*--- Determine the derived data type for sending the data. ---*/
      const unsigned int nElemSend = elementsSend[i].size();
      vector<int> blockLen(nElemSend), displ(nElemSend);

      for(unsigned int j=0; j<nElemSend; ++j) {
        const unsigned long jj = elementsSend[i][j];
        blockLen[j] = nVar*volElem[jj].nDOFsSol;
        displ[j]    = nVar*volElem[jj].offsetDOFsSolLocal;
      }

      MPI_Type_indexed(nElemSend, blockLen.data(), displ.data(),
                       MPI_DOUBLE, &commTypes[nn]);
      MPI_Type_commit(&commTypes[nn]);

      /* Create the communication request for this send operation and
         update the counter nn for the next request. */
      MPI_Send_init(VecSolDOFs.data(), 1, commTypes[nn], ranksSend[i],
                    ranksSend[i], MPI_COMM_WORLD, &commRequests[nn]);
      ++nn;
    }
  }

  /*--- Loop over the ranks from which this rank has to receive data and create
        the receive requests. Exclude self communication. ---*/
  for(unsigned long i=0; i<ranksRecv.size(); ++i) {
    if(ranksRecv[i] != rank) {

      /*--- Determine the derived data type for receiving the data. ---*/
      const unsigned int nElemRecv = elementsRecv[i].size();
      vector<int> blockLen(nElemRecv), displ(nElemRecv);

      for(unsigned int j=0; j<nElemRecv; ++j) {
        const unsigned long jj = elementsRecv[i][j];
        blockLen[j] = nVar*volElem[jj].nDOFsSol;
        displ[j]    = nVar*volElem[jj].offsetDOFsSolLocal;
      }

      MPI_Type_indexed(nElemRecv, blockLen.data(), displ.data(),
                       MPI_DOUBLE, &commTypes[nn]);
      MPI_Type_commit(&commTypes[nn]);

      /* Create the communication request for this receive operation and
         update the counter nn for the next request. */
      MPI_Recv_init(VecSolDOFs.data(), 1, commTypes[nn], ranksRecv[i],
                    rank, MPI_COMM_WORLD, &commRequests[nn]);
      ++nn;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3. Set up the persistent reverse MPI communication for the    ---*/
  /*---         residuals.                                                 ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the number of derived data types necessary in the reverse
     communication. This is the number of ranks to which this rank has to send
     halo data, which corresponds to the original receive pattern. Exclude
     self communication. */
  nn = ranksRecv.size();
  if( elementsRecvSelfComm.size() ) --nn;

  /* Allocate the memory for the reverse communication requests. The sending and
     receiving are reversed, but the total number of requests is the same. Also
     allocate the memory for the derived data types for the reverse communcation
     which are only needed for the sending of the halo data. */
  reverseCommRequests.resize(nCommRequests);
  reverseCommTypes.resize(nn);

  /*--- Loop over the ranks to which this rank has to send residual data and
        create the send requests. Exclude self communication. ---*/
  nn = 0;
  for(unsigned long i=0; i<ranksRecv.size(); ++i) {
    if(ranksRecv[i] != rank) {

      /*--- Determine the derived data type for sending the data. ---*/
      const unsigned int nElemRecv = elementsRecv[i].size();
      vector<int> blockLen(nElemRecv), displ(nElemRecv);

      for(unsigned int j=0; j<nElemRecv; ++j) {
        const unsigned long jj = elementsRecv[i][j];
        blockLen[j] = nVar*volElem[jj].nDOFsSol;
        displ[j]    = nVar*volElem[jj].offsetDOFsSolLocal;
      }

      MPI_Type_indexed(nElemRecv, blockLen.data(), displ.data(),
                       MPI_DOUBLE, &reverseCommTypes[nn]);
      MPI_Type_commit(&reverseCommTypes[nn]);

      /* Create the communication request for this send operation and
         update the counter nn for the next request. */
      MPI_Send_init(VecResDOFs.data(), 1, reverseCommTypes[nn], ranksRecv[i],
                    ranksRecv[i], MPI_COMM_WORLD, &reverseCommRequests[nn]);
      ++nn;
    }
  }

  /*--- Determine the number of receive buffers needed to receive the externally
        computed part of the residuals. Allocate its first index as well as
        the first index of reverseElementsRecv, which stores the corresponding
        elements where the residual must be stored. */
  unsigned int mm = ranksSend.size();
  if( elementsSendSelfComm.size() ) --mm;

  reverseCommRecvBuf.resize(mm);
  reverseElementsRecv.resize(mm);

  /*--- Loop over the ranks from which I receive residual data and
        create the receive requests. Exclude self communication. ---*/
  mm = 0;
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    if(ranksSend[i] != rank) {

      /* Copy the data of elementsSend into reverseElementsRecv. */
      reverseElementsRecv[mm] = elementsSend[i];

      /* Determine the size of the receive buffer and allocate the memory. */
      int sizeBuf = 0;
      for(unsigned long j=0; j<elementsSend[i].size(); ++j) {
        const unsigned long jj = elementsSend[i][j];
        sizeBuf += volElem[jj].nDOFsSol;
      }

      sizeBuf *= nVar;
      reverseCommRecvBuf[mm].resize(sizeBuf);

      /* Create the communication request for this receive operation and
         update the counters mm and nn for the next request. */
      MPI_Recv_init(reverseCommRecvBuf[mm].data(), sizeBuf, MPI_DOUBLE,
                    ranksSend[i], rank, MPI_COMM_WORLD, &reverseCommRequests[nn]);
      ++mm;
      ++nn;
    }
  }

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 4. Store the information for the rotational periodic          ---*/
  /*---         corrections for the halo elements.                         ---*/
  /*--------------------------------------------------------------------------*/

  /* Get the data for the rotational periodic halos from FEMGeometry. */
  vector<unsigned short> markersRotationalPeriodicity = FEMGeometry->GetRotPerMarkers();

  halosRotationalPeriodicity = FEMGeometry->GetRotPerHalos();

  /* Determine the number of rotational periodic transformations and allocate
     the memory to store the corresponding rotation matrices. */
  const unsigned short nRotPerMarkers = markersRotationalPeriodicity.size();
  rotationMatricesPeriodicity.resize(9*nRotPerMarkers);

  /* Loop over the rotational periodic transformations and store the rotation
     matrices in rotationMatricesPeriodicity. */
  unsigned int ii = 0;
  for(unsigned short i=0; i<nRotPerMarkers; ++i) {

   /* Get the rotation angles from config for this marker. */
   const unsigned short pInd = markersRotationalPeriodicity[i];
   su2double *angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(pInd));

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
  }
}

void CFEM_DG_EulerSolver::Initiate_MPI_Communication(void) {

  /*--- Start the MPI communication, if needed. ---*/

#ifdef HAVE_MPI
  if( nCommRequests )
    MPI_Startall(nCommRequests, commRequests.data());
#endif
}

void CFEM_DG_EulerSolver::Complete_MPI_Communication(void) {

  /*-----------------------------------------------------------------------*/
  /*---               Carry out the self communication.                 ---*/
  /*-----------------------------------------------------------------------*/

  /* Easier storage of the number of variables times the size of su2double. */
  const unsigned long sizeEntity = nVar*sizeof(su2double);

  /* Loop over the number of elements involved and copy the data of the DOFs. */
  const unsigned long nSelfElements = elementsSendSelfComm.size();
  for(unsigned long i=0; i<nSelfElements; ++i) {
    const unsigned long indS   = nVar*volElem[elementsSendSelfComm[i]].offsetDOFsSolLocal;
    const unsigned long indR   = nVar*volElem[elementsRecvSelfComm[i]].offsetDOFsSolLocal;
    const unsigned long nBytes = volElem[elementsSendSelfComm[i]].nDOFsSol * sizeEntity;

    memcpy(VecSolDOFs.data()+indR, VecSolDOFs.data()+indS, nBytes);
  }

  /*-----------------------------------------------------------------------*/
  /*---         Complete the MPI communication, if needed.              ---*/
  /*-----------------------------------------------------------------------*/

#ifdef HAVE_MPI
  if( nCommRequests )
    SU2_MPI::Waitall(nCommRequests, commRequests.data(), MPI_STATUSES_IGNORE);
#endif

  /*------------------------------------------------------------------------*/
  /*--- Correct the vector quantities in the rotational periodic halo's. ---*/
  /*------------------------------------------------------------------------*/

  /*--- Loop over the markers for which a rotational periodic
        correction must be applied to the momentum variables. ---*/
  unsigned int ii = 0;
  const unsigned short nRotPerMarkers = halosRotationalPeriodicity.size();
  for(unsigned short k=0; k<nRotPerMarkers; ++k) {

    /* Easier storage of the rotational matrix. */
    su2double rotMatrix[3][3];

    rotMatrix[0][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[0][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[0][2] = rotationMatricesPeriodicity[ii++];

    rotMatrix[1][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][2] = rotationMatricesPeriodicity[ii++];

    rotMatrix[2][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][2] = rotationMatricesPeriodicity[ii++];

    /* Determine the number of elements for this transformation and
       loop over them. */
    const unsigned long nHaloElem = halosRotationalPeriodicity[k].size();
    for(unsigned long j=0; j<nHaloElem; ++j) {

      /* Easier storage of the halo index and loop over its DOFs. */
      const unsigned long ind = halosRotationalPeriodicity[k][j];
      for(unsigned short i=0; i<volElem[ind].nDOFsSol; ++i) {

        /* Determine the pointer in VecSolDOFs where the solution of this DOF
           is stored and copy the momentum variables. Note that a rotational
           correction can only take place for a 3D simulation. */
        su2double *sol = VecSolDOFs.data() + nVar*(volElem[ind].offsetDOFsSolLocal + i);

        const su2double ru = sol[1], rv = sol[2], rw = sol[3];

        /* Correct the momentum variables. */
        sol[1] = rotMatrix[0][0]*ru + rotMatrix[0][1]*rv + rotMatrix[0][2]*rw;
        sol[2] = rotMatrix[1][0]*ru + rotMatrix[1][1]*rv + rotMatrix[1][2]*rw;
        sol[3] = rotMatrix[2][0]*ru + rotMatrix[2][1]*rv + rotMatrix[2][2]*rw;
      }
    }
  }
}

void CFEM_DG_EulerSolver::Initiate_MPI_ReverseCommunication(void) {

  /*------------------------------------------------------------------------*/
  /*---            Accumulate the residuals of the halo DOFs.            ---*/
  /*------------------------------------------------------------------------*/

  /*--- Initialize the residuals of the halo elements to zero. ---*/
  for(unsigned long i=nVar*nDOFsLocOwned; i<nVar*nDOFsLocTot; ++i)
    VecResDOFs[i] = 0.0;

  /*--- Loop over the external DOFs to accumulate the local data. ---*/
  for(unsigned long i=nDOFsLocOwned; i<nDOFsLocTot; ++i) {

    su2double *resDOF = VecResDOFs.data() + nVar*i;

    for(unsigned long j=nEntriesResFaces[i]; j<nEntriesResFaces[i+1]; ++j) {
      const su2double *resFace = VecResFaces.data() + nVar*entriesResFaces[j];
      for(unsigned short k=0; k<nVar; ++k)
        resDOF[k] += resFace[k];
    }
  }

  /*------------------------------------------------------------------------*/
  /*--- Correct the vector residuals in the rotational periodic halo's.  ---*/
  /*------------------------------------------------------------------------*/

  /*--- Loop over the markers for which a rotational periodic
        correction must be applied to the momentum variables. ---*/
  unsigned int ii = 0;
  const unsigned short nRotPerMarkers = halosRotationalPeriodicity.size();
  for(unsigned short k=0; k<nRotPerMarkers; ++k) {

    /* Easier storage of the transpose of the rotational matrix. */
    su2double rotMatrix[3][3];

    rotMatrix[0][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][0] = rotationMatricesPeriodicity[ii++];

    rotMatrix[0][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][1] = rotationMatricesPeriodicity[ii++];

    rotMatrix[0][2] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][2] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][2] = rotationMatricesPeriodicity[ii++];

    /* Determine the number of elements for this transformation and
       loop over them. */
    const unsigned long nHaloElem = halosRotationalPeriodicity[k].size();
    for(unsigned long j=0; j<nHaloElem; ++j) {

      /* Easier storage of the halo index and loop over its DOFs. */
      const unsigned long ind = halosRotationalPeriodicity[k][j];
      for(unsigned short i=0; i<volElem[ind].nDOFsSol; ++i) {

        /* Determine the pointer in VecResDOFs where the residual of this DOF
           is stored and copy the momentum residuals. Note that a rotational
           correction can only take place for a 3D simulation. */
        su2double *res = VecResDOFs.data() + nVar*(volElem[ind].offsetDOFsSolLocal + i);

        const su2double ru = res[1], rv = res[2], rw = res[3];

        /* Correct the momentum variables. */
        res[1] = rotMatrix[0][0]*ru + rotMatrix[0][1]*rv + rotMatrix[0][2]*rw;
        res[2] = rotMatrix[1][0]*ru + rotMatrix[1][1]*rv + rotMatrix[1][2]*rw;
        res[3] = rotMatrix[2][0]*ru + rotMatrix[2][1]*rv + rotMatrix[2][2]*rw;
      }
    }
  }

  /*--- Start the MPI communication, if needed. ---*/

#ifdef HAVE_MPI
  if( nCommRequests )
    MPI_Startall(nCommRequests, reverseCommRequests.data());
#endif

  /*-----------------------------------------------------------------------*/
  /*---               Carry out the self communication.                 ---*/
  /*-----------------------------------------------------------------------*/

  /*--- Loop over the number of elements involved and update the residuals
        of the owned DOFs. ---*/
  const unsigned long nSelfElements = elementsSendSelfComm.size();
  for(unsigned long i=0; i<nSelfElements; ++i) {

    su2double *resOwned = VecResDOFs.data()
                        + nVar*volElem[elementsSendSelfComm[i]].offsetDOFsSolLocal;
    su2double *resHalo  = VecResDOFs.data()
                        + nVar*volElem[elementsRecvSelfComm[i]].offsetDOFsSolLocal;

    const unsigned short nRes = nVar*volElem[elementsSendSelfComm[i]].nDOFsSol;
    for(unsigned short j=0; j<nRes; ++j)
      resOwned[j] += resHalo[j];
  }
}

void CFEM_DG_EulerSolver::Complete_MPI_ReverseCommunication(void) {

#ifdef HAVE_MPI

  /*-----------------------------------------------------------------------*/
  /*---         Complete the MPI communication, if needed.              ---*/
  /*-----------------------------------------------------------------------*/

  if( nCommRequests )
    SU2_MPI::Waitall(nCommRequests, reverseCommRequests.data(), MPI_STATUSES_IGNORE);

  /*---------------------------------------------------------------------------*/
  /*--- Update the residuals of the owned DOFs with the data just received. ---*/
  /*---------------------------------------------------------------------------*/

  /*--- Loop over the ranks from which this rank received residual data. As
        the reverse communication pattern is used, ranksSend must be used.
        Exclude self communication. ---*/
  for(unsigned long i=0; i<reverseElementsRecv.size(); ++i) {

    /* Loop over the elements that must be updated. */
    unsigned long nn = 0;
    for(unsigned long j=0; j<reverseElementsRecv[i].size(); ++j) {
      const unsigned long jj = reverseElementsRecv[i][j];

      /* Easier storage of the starting residual of this element, loop over the
         DOFs of this element and update the residuals of the DOFs. */
      su2double *res = VecResDOFs.data() + nVar*volElem[jj].offsetDOFsSolLocal;
      const unsigned short nRes = nVar*volElem[jj].nDOFsSol;
      for(unsigned short k=0; k<nRes; ++k, ++nn)
        res[k] += reverseCommRecvBuf[i][nn];
    }
  }

#endif
}

void CFEM_DG_EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

  /* Initialize solutionSet. If a solution is set below, this boolean must be
     set to true, such that the solution is communicated to the halo elements
     at the end of this function. */
  bool solutionSet = false;

#ifdef INVISCID_VORTEX

  solutionSet = true;
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* Write a message that the solution is initialized for the inviscid vortex
     test case. */
  if(rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Solution is initialized for the inviscid vortex test case!!!" << endl;
    cout << endl << flush;
  }

  /* The initial conditions are set to the solution of the inviscid vortex,
     which is an exact solution of the Euler equations. The initialization
     below is valid for both 2D and 3D. For the 3D case the z-direction is
     assumed to be the direction in which the solution does not change.
     First set the parameters, which define this test case. */

  const su2double MachVortex  =  0.5;     // Mach number of the undisturbed flow.
  const su2double x0Vortex    = -0.5;     // Initial x-coordinate of the vortex center.
  const su2double y0Vortex    =  0.0;     // Initial y-coordinate of the vortex center.
  const su2double RVortex     =  0.1;     // Radius of the vortex.
  const su2double epsVortex   =  1.0;     // Strength of the vortex.
  const su2double thetaVortex =  0.0;     // Advection angle (in degrees) of the vortex.

  /* Compute the free stream velocities in x- and y-direction. */
  const su2double VelInf = MachVortex*sqrt(Gamma);
  const su2double uInf   = VelInf*cos(thetaVortex*PI_NUMBER/180.0);
  const su2double vInf   = VelInf*sin(thetaVortex*PI_NUMBER/180.0);

  /* Useful coefficients in which Gamma is present. */
  const su2double ovGm1    = 1.0/Gamma_Minus_One;
  const su2double gamOvGm1 = Gamma*ovGm1;

  /* Loop over the owned elements. */
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Loop over the DOFs of this element. */
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {

      // Set the pointer to the solution of this DOF and to the
      // coordinates of its corresponding node ID of the grid.
      su2double *solDOF = VecSolDOFs.data() + nVar*(volElem[i].offsetDOFsSolLocal + j);

      const unsigned long ind = volElem[i].nodeIDsGrid[j];
      const su2double *coor   = meshPoints[ind].coor;

      /* Compute the coordinates relative to the center of the vortex. */
      const su2double dx = coor[0] - x0Vortex;
      const su2double dy = coor[1] - y0Vortex;

      /* Compute the components of the velocity. */
      su2double f  = 1.0 - (dx*dx + dy*dy)/(RVortex*RVortex);
      su2double t1 = epsVortex*dy*exp(0.5*f)/(2.0*PI_NUMBER*RVortex);
      su2double u  = uInf - VelInf*t1;

      t1          = epsVortex*dx*exp(0.5*f)/(2.0*PI_NUMBER*RVortex);
      su2double v = vInf + VelInf*t1;

      /* Compute the density and the pressure. */
      t1 = 1.0 - epsVortex*epsVortex*Gamma_Minus_One
         *       MachVortex*MachVortex*exp(f)/(8.0*PI_NUMBER*PI_NUMBER);

      su2double rho = pow(t1,ovGm1);
      su2double p   = pow(t1,gamOvGm1);

      /* Compute the conservative variables. Note that both 2D and 3D
         cases are treated correctly. */
      solDOF[0]      = rho;
      solDOF[1]      = rho*u;
      solDOF[2]      = rho*v;
      solDOF[3]      = 0.0;
      solDOF[nVar-1] = p*ovGm1 + 0.5*rho*(u*u + v*v);
    }
  }

#elif RINGLEB

  /* The initial conditions are set to the exact solution of the Ringleb flow.
     The reason for doing so, is that the Ringleb flow is an isolated solution
     of the Euler equations. If the initialization is too far off from the
     final solution, shocks develop, which may destabilize the solution and it
     is impossible to obtain a converged solution. */

  solutionSet = true;
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* Write a message that the solution is initialized for the Ringleb test case. */
  if(rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Solution is initialized for the Ringleb test case!!!" << endl;
    cout << endl << flush;
  }

  /* Loop over the owned elements. */
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Loop over the DOFs of this element. */
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {

      /* Set the pointer to the solution of this DOF and to the
         coordinates of its corresponding node ID of the grid. */
      su2double *solDOF = VecSolDOFs.data() + nVar*(volElem[i].offsetDOFsSolLocal + j);

      const unsigned long ind = volElem[i].nodeIDsGrid[j];
      const su2double *coor   = meshPoints[ind].coor;

      /* Compute the conservative flow variables of the Ringleb solution for the
         given coordinates. Note that it is possible to run this case in both 2D
         and 3D, where the z-direction is assumed to be the inactive direction. */
      RinglebSolution(coor, solDOF);
    }
  }

#elif TAYLOR_GREEN

  solutionSet = true;
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* Write a message that the solution is initialized for the Taylor-Green vortex
     test case. */
  if(rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Solution is initialized for the Taylor-Green vortex test case!!!" << endl;
    cout << endl << flush;
  }

  /* The initial conditions are set for the Taylor-Green vortex case, which
   is a DNS case that features vortex breakdown into turbulence. These
   particular settings are for the typical Re = 1600 case (M = 0.08) with
   an initial temperature of 300 K. Note that this condition works in both
   2D and 3D. */

  const su2double tgvLength   = 1.0;     // Taylor-Green length scale.
  const su2double tgvVelocity = 1.0;     // Taylor-Green velocity.
  const su2double tgvDensity  = 1.0;     // Taylor-Green density.
  const su2double tgvPressure = 100.0;   // Taylor-Green pressure.

  /* Useful coefficient in which Gamma is present. */
  const su2double ovGm1    = 1.0/Gamma_Minus_One;

  /* Loop over the owned elements. */
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Loop over the DOFs of this element. */
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {

      // Set the pointer to the solution of this DOF and to the
      // coordinates of its corresponding node ID of the grid.
      su2double *solDOF = VecSolDOFs.data() + nVar*(volElem[i].offsetDOFsSolLocal + j);

      const unsigned long ind = volElem[i].nodeIDsGrid[j];
      const su2double *coor   = meshPoints[ind].coor;

      su2double coorZ = 0.0;
      if (nDim == 3) coorZ = coor[2];

      /* Compute the primitive variables. */
      su2double rho = tgvDensity;
      su2double u   =  tgvVelocity * (sin(coor[0]/tgvLength)*
                                      cos(coor[1]/tgvLength)*
                                      cos(coorZ  /tgvLength));
      su2double v   = -tgvVelocity * (cos(coor[0]/tgvLength)*
                                      sin(coor[1]/tgvLength)*
                                      cos(coorZ  /tgvLength));
      su2double factorA = cos(2.0*coorZ/tgvLength) + 2.0;
      su2double factorB = cos(2.0*coor[0]/tgvLength) + cos(2.0*coor[1]/tgvLength);
      su2double p   = tgvPressure+tgvDensity*(pow(tgvVelocity,2.0)/16.0)*factorA*factorB;

      /* Compute the conservative variables. Note that both 2D and 3D
       cases are treated correctly. */
      solDOF[0]      = rho;
      solDOF[1]      = rho*u;
      solDOF[2]      = rho*v;
      solDOF[3]      = 0.0;
      solDOF[nVar-1] = p*ovGm1 + 0.5*rho*(u*u + v*v);
    }
  }

#endif

  /*--- If the solution was set in this function, perform the MPI
        communication of the solution. ---*/
  if( solutionSet ) {
    Initiate_MPI_Communication();
    Complete_MPI_Communication();
  }
}

void CFEM_DG_EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long ErrorCounter = 0;

  /*-----------------------------------------------------------------------------*/
  /*--- Check for non-physical points. Only needed for a compressible solver. ---*/
  /*-----------------------------------------------------------------------------*/

  if(config->GetKind_Regime() == COMPRESSIBLE) {

    /*--- Loop over the owned DOFs and check for non-physical points. ---*/
    for(unsigned long i=0; i<nDOFsLocOwned; ++i) {

      su2double *solDOF = VecSolDOFs.data() + nVar*i;

      const su2double DensityInv = 1.0/solDOF[0];
      su2double Velocity2 = 0.0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double vel = solDOF[iDim]*DensityInv;
        Velocity2 += vel*vel;
      }

      su2double StaticEnergy = solDOF[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
      su2double Pressure = FluidModel->GetPressure();
      su2double Temperature = FluidModel->GetTemperature();

      if((Pressure < 0.0) || (solDOF[0] < 0.0) || (Temperature < 0.0)) {

        /* Reset the state to the free-stream state. This usually does not work. */
        for(unsigned short j=0; j<nVar; ++j) solDOF[j] = ConsVarFreeStream[j];

        /* Update the error counter. */
        ++ErrorCounter;
      }
    }

    /*--- Collect the number of non-physical points for this iteration. ---*/
    if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
      unsigned long MyErrorCounter = ErrorCounter;
      SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
      if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
    }
  }
}

void CFEM_DG_EulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iMesh) { }

void CFEM_DG_EulerSolver::Set_OldSolution(CGeometry *geometry) {

  memcpy(VecSolDOFsOld.data(), VecSolDOFs.data(), VecSolDOFsOld.size()*sizeof(su2double));
}

void CFEM_DG_EulerSolver::Set_NewSolution(CGeometry *geometry) {

  memcpy(VecSolDOFsNew.data(), VecSolDOFs.data(), VecSolDOFsNew.size()*sizeof(su2double));
}

void CFEM_DG_EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /* Check whether or not a time stepping scheme is used. */
  const bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  /* Initialize the minimum and maximum time step. */
  Min_Delta_Time = 1.e25; Max_Delta_Time = 0.0;

  /* Easier storage of the CFL number. Note that if we are using explicit
   time stepping, the regular CFL condition has been overwritten with the
   unsteady CFL condition in the config post-processing (if non-zero). */

  const su2double CFL = config->GetCFL(iMesh);

  /*--- Explicit time stepping with imposed time step (eventually will
   allow for local time stepping with this value imposed as the time
   for syncing the cells). If the unsteady CFL is set to zero (default),
   it uses the defined unsteady time step, otherwise it computes the time
   step based on the provided unsteady CFL. Note that the regular CFL
   option in the config is always ignored with time stepping. ---*/

  if (time_stepping && (config->GetUnst_CFL() == 0.0)) {

    /*--- Loop over the owned volume elements and set the fixed dt. ---*/
    for(unsigned long i=0; i<nVolElemOwned; ++i)
      VecDeltaTime[i] = config->GetDelta_UnstTimeND();

  } else {

    /*--- Check for a compressible solver. ---*/
    if(config->GetKind_Regime() == COMPRESSIBLE) {

      /*--- Loop over the owned volume elements. ---*/
      for(unsigned long i=0; i<nVolElemOwned; ++i) {

        /*--- Loop over the DOFs of this element and determine
              the maximum wave speed. ---*/
        su2double charVel2Max = 0.0;
        for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
          const su2double *solDOF = VecSolDOFs.data() + nVar*(volElem[i].offsetDOFsSolLocal + j);

          /* Compute the velocities. */
          su2double velAbs[3];
          const su2double DensityInv = 1.0/solDOF[0];
          su2double Velocity2 = 0.0;
          for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
            const su2double vel = solDOF[iDim]*DensityInv;
            velAbs[iDim-1] = fabs(vel);
            Velocity2 += vel*vel;
          }

          /*--- Compute the maximum value of the wave speed. This is a rather
           conservative estimate. ---*/
          const su2double StaticEnergy = solDOF[nDim+1]*DensityInv - 0.5*Velocity2;
          FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
          const su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
          const su2double SoundSpeed  = sqrt(fabs(SoundSpeed2));

          su2double charVel2 = 0.0;
          for(unsigned short iDim=0; iDim<nDim; ++iDim) {
            const su2double rad = velAbs[iDim] + SoundSpeed;
            charVel2 += rad*rad;
          }

          charVel2Max = max(charVel2Max, charVel2);
        }

        /*--- Compute the time step for the element and update the minimum and
              maximum value. Take the factor for time accurate local time
              stepping into account for the minimum and maximum. Note that in
              the length scale the polynomial degree must be taken into account
              for the high order element. ---*/
        const unsigned short ind = volElem[i].indStandardElement;
        unsigned short nPoly = standardElementsSol[ind].GetNPoly();
        if(nPoly == 0) nPoly = 1;

        const su2double lenScaleInv = nPoly/volElem[i].lenScale;
        const su2double dtInv       = lenScaleInv*sqrt(charVel2Max);

        VecDeltaTime[i] = CFL/dtInv;

        Min_Delta_Time = min(Min_Delta_Time, volElem[i].factTimeLevel*VecDeltaTime[i]);
        Max_Delta_Time = max(Max_Delta_Time, volElem[i].factTimeLevel*VecDeltaTime[i]);
      }
    }
    else {

      /*--- Incompressible solver. ---*/

      int rank = MASTER_NODE;
#ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      if(rank == MASTER_NODE) {
        cout << "In function CFEM_DG_EulerSolver::SetTime_Step" << endl;
        cout << "Incompressible solver not implemented yet" << endl;
      }

#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#else
      exit(EXIT_FAILURE);
#endif

    }

    /*--- Compute the max and the min dt (in parallel). Note that we only
     do this for steady calculations if the high verbosity is set, but we
     always perform the reduction for unsteady calculations where the CFL
     limit is used to set the global time step. ---*/
    if ((config->GetConsole_Output_Verb() == VERB_HIGH) || time_stepping) {
#ifdef HAVE_MPI
      su2double rbuf_time = Min_Delta_Time;
      SU2_MPI::Allreduce(&rbuf_time, &Min_Delta_Time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      rbuf_time = Max_Delta_Time;
      SU2_MPI::Allreduce(&rbuf_time, &Max_Delta_Time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    }

    /*--- For explicit time stepping with an unsteady CFL imposed, use the
          minimum delta time of the entire mesh. As Min_Delta_Time is scaled to
          the time step of the largest time level, a correction must be used
          for the time level when time accurate local time stepping is used. ---*/
    if (time_stepping) {
      for(unsigned long i=0; i<nVolElemOwned; ++i)
        VecDeltaTime[i] = Min_Delta_Time/volElem[i].factTimeLevel;
    }
  }
}

void CFEM_DG_EulerSolver::CheckTimeSynchronization(CConfig         *config,
                                                   const su2double TimeSync,
                                                   su2double       &timeEvolved,
                                                   bool            &syncTimeReached) {

  /* Check if this is the first time this check is carried out
     and determine the new time evolved. */
  const bool firstTime = timeEvolved == 0.0;
  timeEvolved         += Min_Delta_Time;

  /*--- Check for a (too) small a value for the synchronization time and
        print a warning if this happens. ---*/
  if(firstTime && timeEvolved >= 1.5*TimeSync) {

    int rank = MASTER_NODE;
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if(rank == MASTER_NODE) {
      cout << endl << "              WARNING" << endl;
      cout << "The specified synchronization time is " << timeEvolved/TimeSync
           << " times smaller than the time step for stability" << endl;
      cout << "This is inefficient!!!!!" << endl << endl;
    }
  }

  /*--- If the current value of timeEvolved is larger or equal than the
        synchronization time, syncTimeReached is set to true and a correction
        to the time step is carried out. The factor for the time accurate local
        time stepping must be taken into account. If the synchronization time
        has not been reached yet, syncTimeReached is set to false. ---*/
  if(timeEvolved >= TimeSync) {
    syncTimeReached = true;
    const su2double newDeltaTime = Min_Delta_Time + (TimeSync - timeEvolved);

    for(unsigned long i=0; i<nVolElemOwned; ++i)
      VecDeltaTime[i] = newDeltaTime/volElem[i].factTimeLevel;
  }
  else syncTimeReached = false;
}

void CFEM_DG_EulerSolver::ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                                    CNumerics **numerics, CConfig *config,
                                                    unsigned short iMesh, unsigned short RunTime_EqSystem) {

  su2double tick = 0.0;

  /*--- Preprocessing ---*/
  config->Tick(&tick);
  Preprocessing(geometry, solver_container, config, iMesh, 0, RunTime_EqSystem, false);
  config->Tock(tick,"Preprocessing",2);

  /*--- Set the old solution, carry out the predictor step and set the number
        of time integration points and the corresponding weights. ---*/
  Set_OldSolution(geometry);

  config->Tick(&tick);
  ADER_DG_PredictorStep(config, 0);
  config->Tock(tick,"ADER_DG_PredictorStep",2);

  const unsigned short nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();
  const su2double     *timeIntegrationWeights = config->GetWeightsIntegrationADER_DG();

  /* Loop over the number of integration points in time. */
  for(unsigned short iTime=0; iTime<nTimeIntegrationPoints; iTime++) {

    /* The predictor solution must be interpolated in time to the correct
       time location of the current integration point. */
    config->Tick(&tick);
    ADER_DG_TimeInterpolatePredictorSol(config, iTime);
    config->Tock(tick,"ADER_DG_TimeInterpolatePredictorSol",3);

    /* Compute the artificial viscosity for shock capturing in DG. */
    config->Tick(&tick);
    Shock_Capturing_DG(geometry, solver_container, numerics[CONV_TERM], config, iMesh, 0);
    config->Tock(tick,"Shock_Capturing",3);

    /*--- Compute the volume portion of the residual. ---*/
    config->Tick(&tick);
    Volume_Residual(geometry, solver_container, numerics[CONV_TERM], config, iMesh, 0);
    config->Tock(tick,"Volume_Residual",3);

    /*--- Compute source term residuals ---*/
    config->Tick(&tick);
    Source_Residual(geometry, solver_container, numerics[SOURCE_FIRST_TERM], numerics[SOURCE_SECOND_TERM], config, iMesh);
    config->Tock(tick,"Source_Residual",3);

    /*--- Boundary conditions ---*/
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case EULER_WALL:
          config->Tick(&tick);
          BC_Euler_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
          config->Tock(tick,"BC_Euler_Wall",3);
          break;
        case FAR_FIELD:
          config->Tick(&tick);
          BC_Far_Field(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
          config->Tock(tick,"BC_Far_Field",3);
          break;
        case SYMMETRY_PLANE:
          config->Tick(&tick);
          BC_Sym_Plane(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
          config->Tock(tick,"BC_Sym_Plane",3);
          break;
        case INLET_FLOW:
          config->Tick(&tick);
          BC_Inlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
          config->Tock(tick,"BC_Inlet",3);
          break;
        case OUTLET_FLOW:
          config->Tick(&tick);
          BC_Outlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
          config->Tock(tick,"BC_Outlet",3);
          break;
        case ISOTHERMAL:
          config->Tick(&tick);
          BC_Isothermal_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
          config->Tock(tick,"BC_Isothermal_Wall",3);
          break;
        case HEAT_FLUX:
          config->Tick(&tick);
          BC_HeatFlux_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
          config->Tock(tick,"BC_HeatFlux_Wall",3);
          break;
        case CUSTOM_BOUNDARY:
          config->Tick(&tick);
          BC_Custom(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
          config->Tock(tick,"BC_Custom",3);
          break;
        case PERIODIC_BOUNDARY:  // Nothing to be done for a periodic boundary.
          break;
        default:
          cout << "BC not implemented." << endl;
#ifndef HAVE_MPI
          exit(EXIT_FAILURE);
#else
          MPI_Abort(MPI_COMM_WORLD,1);
          MPI_Finalize();
#endif
      }
    }

    /*--- Compute surface portion of the residual. ---*/
    config->Tick(&tick);
    Surface_Residual(geometry, solver_container, numerics[CONV_TERM], config, iMesh, 0);
    config->Tock(tick,"Surface_Residual",3);

    /*--- Accumulate the space time residual. ---*/
    AccumulateSpaceTimeResidualADER(iTime, timeIntegrationWeights[iTime]);
  }

  /*--- Multiply the residual by the (lumped) mass matrix, to obtain the final value. ---*/
  config->Tick(&tick);
  MultiplyResidualByInverseMassMatrix(config, true);
  config->Tock(tick,"MultiplyResidualByInverseMassMatrix",3);

  /*--- Perform the time integration ---*/
  config->Tick(&tick);
  ADER_DG_Iteration(geometry, solver_container, config, 0);
  config->Tock(tick,"ADER_DG_Iteration",3);

  /*--- Postprocessing ---*/
  config->Tick(&tick);
  Postprocessing(geometry, solver_container, config, iMesh);
  config->Tock(tick,"Postprocessing",2);
}

void CFEM_DG_EulerSolver::ADER_DG_PredictorStep(CConfig *config, unsigned short iStep) {

  /*----------------------------------------------------------------------*/
  /*---        Get the data of the ADER time integration scheme.       ---*/
  /*----------------------------------------------------------------------*/

  const unsigned short nTimeDOFs              = config->GetnTimeDOFsADER_DG();
  const unsigned short nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();
  const su2double     *timeIntegrationWeights = config->GetWeightsIntegrationADER_DG();
  const bool          useAliasedPredictor     = config->GetKind_ADER_Predictor() == ADER_ALIASED_PREDICTOR;

  /*----------------------------------------------------------------------*/
  /*--- Determine the reference values for the conservative variables. ---*/
  /*----------------------------------------------------------------------*/

  /* Initialization to zero. */
  su2double URef[] = {0.0, 0.0, 0.0, 0.0, 0.0};

  /* Loop over the owned DOFs to determine the maximum values of the
     conservative variables on this rank. */
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + i*nVar;

    for(unsigned short j=0; j<nVar; ++j) {
      const su2double solAbs = fabs(solDOF[j]);
      if(solAbs > URef[j]) URef[j] = solAbs;
    }
  }

  /*--- Determine the maximum values on all ranks in case of a
        parallel computation. */
#ifdef HAVE_MPI
  su2double URefLoc[5];
  for(unsigned short i=0; i<nVar; ++i) URefLoc[i] = URef[i];

  SU2_MPI::Allreduce(URefLoc, URef, nVar, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  /*--- Determine the maximum scale of the momentum variables and adapt the
        corresponding values of URef accordingly. ---*/
  su2double momRef = URef[1];
  for(unsigned short i=2; i<=nDim; ++i) momRef  = max(momRef, URef[i]);
  for(unsigned short i=1; i<=nDim; ++i) URef[i] = momRef;

  /*--------------------------------------------------------------------------*/
  /*--- Loop over the owned elements to compute the predictor solution.    ---*/
  /*--- For the predictor solution only an integration over the element is ---*/
  /*--- performed to obtain the weak formulation, i.e. no integration by   ---*/
  /*--- parts. As a consequence there is no surface term and also no       ---*/
  /*--- coupling with neighboring elements. This is also the reason why    ---*/
  /*--- the ADER scheme is only conditionally stable.                      ---*/
  /*--------------------------------------------------------------------------*/

  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Easier storage of the number of spatial DOFs. */
    const unsigned short nDOFs = volElem[l].nDOFsSol;

    /* Set the pointers for the working variables. */
    su2double *resInt = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    su2double *solPred = VecTmpMemory.data();
    su2double *solOld  = solPred + nVar*nDOFs*nTimeDOFs;
    su2double *resSol  = solOld  + nVar*nDOFs*nTimeDOFs;
    su2double *resTot  = resSol  + nVar*nDOFs*nTimeDOFs;
    su2double *solInt  = resTot  + nVar*nDOFs*nTimeDOFs;
    su2double *work    = solInt  + nVar*nDOFs;

    /* Initialize the predictor solution to the current solution. */
    const su2double    *solCur = VecSolDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;
    const unsigned long nBytes = nVar*nDOFs*sizeof(su2double);

    for(unsigned short j=0; j<nTimeDOFs; ++j) {
      su2double *solPredTimeInd = solPred + j*nVar*nDOFs;
      memcpy(solPredTimeInd, solCur, nBytes);
    }

    /*-------------------------------------------------------------------------*/
    /*--- Compute the contribution of the current solution to the residual. ---*/
    /*--- As this a constant contribution, it is stored in resSol.          ---*/
    /*-------------------------------------------------------------------------*/

    /* First compute the product of the mass matrix and the current solution. */
    DenseMatrixProduct(nDOFs, nVar, nDOFs, volElem[l].massMatrix.data(),
                       solCur, resSol);

    /* Loop over the number of time DOFs and multiply this product by the  */
    /* value of the Lagrangian interpolation function. Note that this loop */
    /* starts at 1, such that the initial value is not touched yet.        */
    for(unsigned short j=1; j<nTimeDOFs; ++j) {

      su2double *resSolTimeInd = resSol + j*nVar*nDOFs;
      for(unsigned short i=0; i<(nVar*nDOFs); ++i)
        resSolTimeInd[i] = resSol[i]*LagrangianBeginTimeIntervalADER_DG[j];
    }

    /* Multiply the values of resSol for the first time DOF with the */
    /* value of its corresponding Lagrangian interpolation function. */
    for(unsigned short i=0; i<(nVar*nDOFs); ++i)
      resSol[i] *= LagrangianBeginTimeIntervalADER_DG[0];

    /*-------------------------------------------------------------------------*/
    /*--- Iterative algorithm to compute the predictor solution for all     ---*/
    /*--- the time DOFs of this element simultaneously.                     ---*/
    /*-------------------------------------------------------------------------*/

    for(;;) {

      /* Initialize the total residual to resSol and store the current solution
         in solOld, such that the updates can be determined after the new
         solution has been computed. */
      memcpy(resTot, resSol,  nTimeDOFs*nBytes);
      memcpy(solOld, solPred, nTimeDOFs*nBytes);

      /* Loop over the number of integration points in time to compute the
         space time integral of the spatial derivatives of the fluxes. */
      for(unsigned short intPoint=0; intPoint<nTimeIntegrationPoints; ++intPoint) {

        /*--------------------------------------------------------------------*/
        /*--- Interpolate the predictor solution to the current time       ---*/
        /*--- integration point. It is likely that this integration        ---*/
        /*--- coincides with the location of one of the time DOFs. When    ---*/
        /*--- this is the case the interpolation boils down to a copy.     ---*/
        /*--------------------------------------------------------------------*/

        /* Store the interpolation data for this time integration point. */
        const su2double *DOFToThisTimeInt = timeInterpolDOFToIntegrationADER_DG
                                          + intPoint*nTimeDOFs;

        /* Initialize the interpolated solution to zero. */
        for(unsigned short i=0; i<(nVar*nDOFs); ++i) solInt[i] = 0.0;

        /* Carry out the actual interpolation. */
        for(unsigned short j=0; j<nTimeDOFs; ++j) {

          const su2double *solPredTimeInd = solPred + j*nVar*nDOFs;
          for(unsigned short i=0; i<(nVar*nDOFs); ++i)
            solInt[i] += DOFToThisTimeInt[j]*solPredTimeInd[i];
        }

        /*--------------------------------------------------------------------*/
        /*--- Compute the spatial residual of the predictor step for the   ---*/
        /*--- current time integration point. A distinction is             ---*/
        /*--- made between an aliased and a non-aliased evaluation of the  ---*/
        /*--- predictor residual.                                          ---*/
        /*--------------------------------------------------------------------*/

        if( useAliasedPredictor )
          ADER_DG_AliasedPredictorResidual(config, &volElem[l], solInt, resInt, work);
        else
          ADER_DG_NonAliasedPredictorResidual(config, &volElem[l], solInt, resInt, work);

        /*--------------------------------------------------------------------*/
        /*--- Update the total residual with the residual of the current   ---*/
        /*--- integration point in time. Note the minus sign, because the  ---*/
        /*--- residual is put on the RHS of the equation.                  ---*/
        /*--------------------------------------------------------------------*/

        /* Loop over all the time DOFs. */
        for(unsigned short j=0; j<nTimeDOFs; ++j) {

          /* Determine the multiplication factor for this time DOF. The factor
             0.5 is present, because the length of the interval [-1..1] is 2. */
          const su2double w = 0.5*timeIntegrationWeights[intPoint]
                            * VecDeltaTime[l]*DOFToThisTimeInt[j];

          /* Update the residual of this time DOF. */
          su2double *res = resTot + j*nVar*nDOFs;
          for(unsigned short i=0; i<(nVar*nDOFs); ++i)
            res[i] -= w*resInt[i];
        }
      }

      /* Solve for the new values of solPred, which are obtained by
         carrying out the matrix product iterMat X resTot. */
      DenseMatrixProduct(nDOFs*nTimeDOFs, nVar, nDOFs*nTimeDOFs,
                         volElem[l].ADERIterationMatrix.data(),
                         resTot, solPred);

      /* Compute the L2 norm of the updates. */
      su2double L2[5];
      for(unsigned short i=0; i<nVar; ++i) L2[i] = 0.0;

      for(unsigned short j=0; j<(nDOFs*nTimeDOFs); ++j) {
        for(unsigned short i=0; i<nVar; ++i) {
          const unsigned short ind  = j*nVar + i;
          const su2double      diff = solPred[ind] - solOld[ind];
          L2[i] += diff*diff;
        }
      }

      for(unsigned short i=0; i<nVar; ++i) L2[i] = sqrt(L2[i]/(nDOFs*nTimeDOFs));

      /* Check for convergence. */
      bool converged = true;
      for(unsigned short i=0; i<nVar; ++i) {
        if(L2[i] > 1.e-6*URef[i]) converged = false;
      }

      if( converged ) break;
    }

    /* Store the predictor solution in the correct location of
       VecSolDOFsPredictorADER. */
    for(unsigned short j=0; j<nTimeDOFs; ++j) {
      su2double *solCur      = VecSolDOFsPredictorADER.data()
                             + nVar*(j*nDOFsLocOwned + volElem[l].offsetDOFsSolLocal);
      su2double *solPredTime = solPred + j*nVar*nDOFs;

      memcpy(solCur, solPredTime, nVar*nDOFs*sizeof(su2double));
    }
  }
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual(CConfig           *config,
                                                           CVolumeElementFEM *elem,
                                                           const su2double   *sol,
                                                           su2double         *res,
                                                           su2double         *work) {

  /*--- Get the necessary information from the standard element. ---*/
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *matDerBasisInt         = matBasisInt + nDOFs*nInt;
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  /* Set the pointers for fluxesDOF, gradFluxesInt and divFlux. Note that the
     same array can be used for the first and last array, because divFlux is
     needed after fluxesDOF. */
  su2double *fluxesDOF     = work;
  su2double *gradFluxesInt = fluxesDOF + nDOFs*nVar*nDim;
  su2double *divFlux       = work;

  /* Determine the offset between the r-derivatives and s-derivatives, which is
     also the offset between s- and t-derivatives, of the fluxes. */
  const unsigned short offDeriv = nVar*nInt*nDim;

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /*-- Compute the Cartesian fluxes in the DOFs. ---*/
  for(unsigned short i=0; i<nDOFs; ++i) {

    /* Set the pointers to the location where the solution of this DOF is
       stored and the location where the Cartesian fluxes are stored. */
    const su2double *solDOF = sol       + i*nVar;
    su2double       *fluxes = fluxesDOF + i*nVar*nDim;

    /*--- Compute the velocities and pressure in this DOF. ---*/
    const su2double DensityInv = 1.0/solDOF[0];
    su2double vel[3], Velocity2 = 0.0;
    for(unsigned short j=0; j<nDim; ++j) {
      vel[j]     = solDOF[j+1]*DensityInv;
      Velocity2 += vel[j]*vel[j];
    }

    const su2double TotalEnergy  = solDOF[nDim+1]*DensityInv;
    const su2double StaticEnergy = TotalEnergy - 0.5*Velocity2;

    FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
    const su2double Pressure = FluidModel->GetPressure();

    /* Loop over the number of dimensions for the number of fluxes. */
    unsigned short ll = 0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim) {

      /* Mass flux for the current direction. */
      fluxes[ll++] = solDOF[iDim+1];

      /* Momentum fluxes for the current direction. */
      for(unsigned short jDim=0; jDim<nDim; ++jDim)
        fluxes[ll+jDim] = solDOF[iDim+1]*vel[jDim];
      fluxes[ll+iDim]  += Pressure;

      ll += nDim;

      /* Energy flux for the current direction. */
      fluxes[ll++] = (solDOF[nDim+1] + Pressure)*vel[iDim];
    }
  }

  /* Compute the derivatives of the Cartesian fluxes w.r.t. the parametric
     coordinates in the integration points. */
  DenseMatrixProduct(nInt*nDim, nVar*nDim, nDOFs, matDerBasisInt, fluxesDOF,
                     gradFluxesInt);

  /*--- Loop over the integration points to compute the divergence of the inviscid
        fluxes in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the data of the derivatives of the
       fluxes of this integration point starts. */
    const su2double *gradFluxes = gradFluxesInt + nVar*nDim*i;

    /* Easier storage of the metric terms in this integration point. The +1
       is present, because the first element of the metric terms is the
       Jacobian in the integration point. Also set the point where the divergence
       of the flux terms is stored for the current integration point. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint + 1;
    su2double       *divFluxInt  = divFlux + nVar*i;

    /* Initialize the divergence to zero. */
    for(unsigned short j=0; j<nVar; ++j) divFluxInt[j] = 0.0;

    /* Loop over the nDim parametric coordinates. */
    unsigned short ll = 0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim) {

      /* Set the pointer to derivatives of this parametric coordinate. */
      const su2double *derFluxes = gradFluxes + iDim*offDeriv;

      /* Loop over the number of Cartesian dimensions. */
      for(unsigned short jDim=0; jDim<nDim; ++jDim, ++ll) {

        /* Set the pointer to the derivatives of the Cartesian flux in
           the jDim direction. */
        const su2double *derFlux = derFluxes + jDim*nVar;

        /* Update the divergence for all equations. */
        for(unsigned short j=0; j<nVar; ++j)
          divFluxInt[j] += derFlux[j]*metricTerms[ll];
      }
    }

    /* Multiply the divergence with the integration weight. */
    for(unsigned short j=0; j<nVar; ++j) divFluxInt[j] *= weights[i];
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual(CConfig           *config,
                                                              CVolumeElementFEM *elem,
                                                              const su2double   *sol,
                                                              su2double         *res,
                                                              su2double         *work) {

  /* Set the pointers for solAndGradInt and divFlux to work. The same array
     can be used for both help arrays. */
  su2double *solAndGradInt = work;
  su2double *divFlux       = work;

  /*--- Get the necessary information from the standard element. ---*/
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  /* Determine the offset between the solution variables and the r-derivatives,
     which is also the offset between the r- and s-derivatives and the offset
     between s- and t-derivatives. */
  const unsigned short offDeriv = nVar*nInt;

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /* Compute the solution and the derivatives w.r.t. the parametric coordinates
     in the integration points. */
  DenseMatrixProduct(nInt*(nDim+1), nVar, nDOFs, matBasisInt, sol, solAndGradInt);

  /*--- Loop over the integration points to compute the divergence of the inviscid
        fluxes in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the solution data of this
       integration point starts. */
    const su2double *sol = solAndGradInt + nVar*i;

    /*--- Compute the velocities and static energy in this integration point. ---*/
    const su2double DensityInv = 1.0/sol[0];
    su2double vel[3], Velocity2 = 0.0;
    for(unsigned short j=0; j<nDim; ++j) {
      vel[j]     = sol[j+1]*DensityInv;
      Velocity2 += vel[j]*vel[j];
    }

    const su2double TotalEnergy  = sol[nDim+1]*DensityInv;
    const su2double StaticEnergy = TotalEnergy - 0.5*Velocity2;

    /*--- Compute the pressure and the total enthalpy. ---*/
    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
    const su2double Pressure = FluidModel->GetPressure();
    const su2double Htot     = (sol[nDim+1]+Pressure)*DensityInv;

    /* Easier storage of the metric terms in this integration point. The +1
       is present, because the first element of the metric terms is the
       Jacobian in the integration point. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint + 1;

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the
          metric terms in this integration point. Note that these gradients
          must be scaled with the Jacobian. This scaling is already present
          in the metric terms. ---*/
    su2double solGradCart[5][3];
    for(unsigned short k=0; k<nDim; ++k) {
      for(unsigned short j=0; j<nVar; ++j) {
        solGradCart[j][k] = 0.0;
        for(unsigned short l=0; l<nDim; ++l)
          solGradCart[j][k] += sol[j+(l+1)*offDeriv]*metricTerms[k+l*nDim];
      }
    }

    /*--- Compute the Cartesian gradients of the pressure, scaled with
          the Jacobian. ---*/
    su2double pGradCart[] = {0.0, 0.0, 0.0};
    for(unsigned short k=0; k<nDim; ++k) {
      pGradCart[k] = solGradCart[nVar-1][k] + 0.5*Velocity2*solGradCart[0][k];
      for(unsigned short l=0; l<nDim; ++l)
        pGradCart[k] -= vel[l]*solGradCart[l+1][k];
      pGradCart[k] *= Gamma_Minus_One;
    }

    /*--- Abbreviations, which make it easier to compute the divergence term. ---*/
    su2double abv1 = 0.0, abv2 = 0.0, abv3 = 0.0;
    for(unsigned short k=0; k<nDim; ++k) {
      abv1 += solGradCart[k+1][k];
      abv2 += vel[k]*solGradCart[0][k];
      abv3 += vel[k]*(solGradCart[nVar-1][k] + pGradCart[k]);
    }

    /*--- Set the pointer to store the divergence terms for this integration point.
          and compute these terms, multiplied by the integration weight. ---*/
    su2double *divFluxInt = divFlux + nVar*i;

    divFluxInt[0]      = weights[i]*abv1;
    divFluxInt[nVar-1] = weights[i]*(abv3 + Htot*(abv1 - abv2));

    for(unsigned short k=0; k<nDim; ++k) {
      divFluxInt[k+1] = pGradCart[k] + vel[k]*(abv1 - abv2);
      for(unsigned short l=0; l<nDim; ++l)
        divFluxInt[k+1] += vel[l]*solGradCart[k+1][l];

      divFluxInt[k+1] *= weights[i];
    }
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_EulerSolver::ADER_DG_TimeInterpolatePredictorSol(CConfig       *config,
                                                              unsigned short iTime) {

  /* Easier storage of the interpolation coefficients for the current time
     integration point (iTime).       */
  const unsigned short nTimeDOFs    = config->GetnTimeDOFsADER_DG();
  const su2double *DOFToThisTimeInt = timeInterpolDOFToIntegrationADER_DG
                                    + iTime*nTimeDOFs;

  /* Initialize the solution of the owned DOFs to zero. */
  for(unsigned long i=0; i<(nVar*nDOFsLocOwned); ++i)
     VecSolDOFs[i] = 0.0;

  /* Loop over the time DOFs, for which the predictor solution is present. */
  for(unsigned short j=0; j<nTimeDOFs; ++j) {

    /* Add the contribution of this predictor solution to the interpolated solution. */
    const su2double *solPred = VecSolDOFsPredictorADER.data() + j*nVar*nDOFsLocOwned;
    for(unsigned long i=0; i<(nVar*nDOFsLocOwned); ++i)
       VecSolDOFs[i] += DOFToThisTimeInt[j]*solPred[i];
  }
}

void CFEM_DG_EulerSolver::Shock_Capturing_DG(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short iMesh, unsigned short iStep) {

  /*--- Run shock capturing algorithm ---*/
  switch( config->GetKind_FEM_DG_Shock() ) {
    case NONE:
      break;
  }
}

void CFEM_DG_EulerSolver::Volume_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short iMesh, unsigned short iStep) {

  /* Start the MPI communication of the solution in the halo elements. */
  Initiate_MPI_Communication();

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solInt = VecTmpMemory.data();
  su2double *fluxes = solInt + nIntegrationMax*nVar;
  su2double tick = 0.0;

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = ctc::nDim*ctc::nDim + 1;

  /*--- Loop over the owned volume elements to compute the contribution of the
        volume integral in the DG FEM formulation to the residual.       ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Get the data from the corresponding standard element. */
    const unsigned short ind             = volElem[l].indStandardElement;
    const unsigned short nInt            = standardElementsSol[ind].GetNIntegration();
    const unsigned short nDOFs           = volElem[l].nDOFsSol;
    const su2double *matBasisInt         = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
    const su2double *matDerBasisIntTrans = standardElementsSol[ind].GetDerMatBasisFunctionsIntTrans();
    const su2double *weights             = standardElementsSol[ind].GetWeightsIntegration();

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Interpolate the solution to the integration points of    ---*/
    /*---         the element.                                             ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the solution variables for this element. */
    su2double *solDOFs = VecSolDOFs.data() + ctc::nVar*volElem[l].offsetDOFsSolLocal;

    /* Call the general function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nInt, nVar, nDOFs, matBasisInt, solDOFs, solInt);
    config->GEMM_Tock(tick, "Volume_Residual1", nInt, nVar, nDOFs);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the inviscid fluxes, multiplied by minus the     ---*/
    /*---         integration weight, in the integration points.           ---*/
    /*------------------------------------------------------------------------*/

    /* Make a distinction between two and three space dimensions
        in order to have the most efficient code. */
    switch( ctc::nDim ) {

      case 2: {

        /* 2D simulation. Loop over the integration points to compute
           the fluxes. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the metric terms in this integration point and
             compute the inverse of the Jacobian. */
          const su2double *metricTerms = volElem[l].metricTerms.data()
                                       + i*nMetricPerPoint;

          /* Compute the metric terms multiplied by minus the integration weight.
             The minus sign comes from the integration by parts in the weak
             formulation. */
          const su2double wDrdx = -weights[i]*metricTerms[1];
          const su2double wDrdy = -weights[i]*metricTerms[2];

          const su2double wDsdx = -weights[i]*metricTerms[3];
          const su2double wDsdy = -weights[i]*metricTerms[4];

          /* Easier storage of the location where the solution data of this
             integration point starts. */
          const su2double *sol = solInt + ctc::nVar*i;

          /*--- Compute the velocities and static energy in this integration point. ---*/
          const su2double rhoInv       = 1.0/sol[0];
          const su2double u            = sol[1]*rhoInv;
          const su2double v            = sol[2]*rhoInv;
          const su2double TotalEnergy  = sol[3]*rhoInv;
          const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

          /*--- Compute the pressure. ---*/
          FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
          const su2double Pressure = FluidModel->GetPressure();

          /* Set the pointer for the fluxes in this integration point. */
          su2double *flux = fluxes + i*ctc::nDim*ctc::nVar;

          /*--- Fluxes in r-direction. */
          const su2double H     = TotalEnergy + rhoInv*Pressure;
          const su2double rhoUr = sol[1]*wDrdx + sol[2]*wDrdy;

          flux[0] = rhoUr;
          flux[1] = rhoUr*u + Pressure*wDrdx;
          flux[2] = rhoUr*v + Pressure*wDrdy;
          flux[3] = rhoUr*H;

          /*--- Fluxes in s-direction. */
          const su2double rhoUs = sol[1]*wDsdx + sol[2]*wDsdy;

          flux[4] = rhoUs;
          flux[5] = rhoUs*u + Pressure*wDsdx;
          flux[6] = rhoUs*v + Pressure*wDsdy;
          flux[7] = rhoUs*H;
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /* 3D simulation. Loop over the integration points to compute
           the fluxes. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the metric terms in this integration point and
             compute the inverse of the Jacobian. */
          const su2double *metricTerms = volElem[l].metricTerms.data()
                                       + i*nMetricPerPoint;

          /* Compute the metric terms multiplied by minus the integration weight.
             The minus sign comes from the integration by parts in the weak
             formulation. */
          const su2double wDrdx = -weights[i]*metricTerms[1];
          const su2double wDrdy = -weights[i]*metricTerms[2];
          const su2double wDrdz = -weights[i]*metricTerms[3];

          const su2double wDsdx = -weights[i]*metricTerms[4];
          const su2double wDsdy = -weights[i]*metricTerms[5];
          const su2double wDsdz = -weights[i]*metricTerms[6];

          const su2double wDtdx = -weights[i]*metricTerms[7];
          const su2double wDtdy = -weights[i]*metricTerms[8];
          const su2double wDtdz = -weights[i]*metricTerms[9];

          /* Easier storage of the location where the solution data of this
             integration point starts. */
          const su2double *sol    = solInt + ctc::nVar*i;

          /*--- Compute the velocities and static energy in this integration point. ---*/
          const su2double rhoInv       = 1.0/sol[0];
          const su2double u            = sol[1]*rhoInv;
          const su2double v            = sol[2]*rhoInv;
          const su2double w            = sol[3]*rhoInv;
          const su2double TotalEnergy  = sol[4]*rhoInv;
          const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v + w*w);

          /*--- Compute the pressure. ---*/
          FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
          const su2double Pressure = FluidModel->GetPressure();

          /* Set the pointer for the fluxes in this integration point. */
          su2double *flux = fluxes + i*ctc::nDim*ctc::nVar;

          /*--- Fluxes in r-direction. */
          const su2double H     = TotalEnergy + rhoInv*Pressure;
          const su2double rhoUr = sol[1]*wDrdx + sol[2]*wDrdy + sol[3]*wDrdz;

          flux[0] = rhoUr;
          flux[1] = rhoUr*u + Pressure*wDrdx;
          flux[2] = rhoUr*v + Pressure*wDrdy;
          flux[3] = rhoUr*w + Pressure*wDrdz;
          flux[4] = rhoUr*H;

          /*--- Fluxes in s-direction. */
          const su2double rhoUs = sol[1]*wDsdx + sol[2]*wDsdy + sol[3]*wDsdz;

          flux[5] = rhoUs;
          flux[6] = rhoUs*u + Pressure*wDsdx;
          flux[7] = rhoUs*v + Pressure*wDsdy;
          flux[8] = rhoUs*w + Pressure*wDsdz;
          flux[9] = rhoUs*H;

          /*--- Fluxes in t-direction. */
          const su2double rhoUt = sol[1]*wDtdx + sol[2]*wDtdy + sol[3]*wDtdz;

          flux[10] = rhoUt;
          flux[11] = rhoUt*u + Pressure*wDtdx;
          flux[12] = rhoUt*v + Pressure*wDtdy;
          flux[13] = rhoUt*w + Pressure*wDtdz;
          flux[14] = rhoUt*H;
        }

        break;
      }
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
    /*---         integration over the volume element.                     ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the residuals for this volume element. */
    su2double *res = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Call the general function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nDOFs, nVar, nInt*nDim, matDerBasisIntTrans, fluxes, res);
    config->GEMM_Tock(tick, "Volume_Residual2", nDOFs, nVar, nInt*nDim);

  }
}

void CFEM_DG_EulerSolver::Source_Residual(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CNumerics *numerics,
                                          CNumerics *second_numerics,
                                          CConfig *config,
                                          unsigned short iMesh) {

  /*--- Stub for source term integration. ---*/

  bool body_force = config->GetBody_Force();

  if (body_force) {

    su2double *body_force_vector = config->GetBody_Force_Vector();

    /*--- Source term integration goes here... dummy output for now. ---*/

    cout << " Applying a body force of (";
    for( unsigned short iDim = 0; iDim < nDim; iDim++) {
      cout << body_force_vector[iDim];
      if (iDim < nDim-1) cout << ", ";
    }
    cout << ")." << endl;

  }

}

void CFEM_DG_EulerSolver::Surface_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                           CConfig *config, unsigned short iMesh, unsigned short iStep) {

  /* Complete the MPI communication of the solution data. */
  Complete_MPI_Communication();

  /* Compute the residual of the faces that involve a halo element. */
  unsigned long indResFaces = startLocResInternalFacesWithHaloElem[0];
  ResidualFaces(geometry, solver_container, numerics, config, iMesh, iStep,
                nMatchingInternalFacesWithHaloElem[0],
                nMatchingInternalFacesWithHaloElem[1], indResFaces);

  /* Start the communication of the residuals, for which the
     reverse communication must be used. */
  Initiate_MPI_ReverseCommunication();

  /* Compute the residual of the faces that only involve owned elements. */
  indResFaces = startLocResInternalFacesLocalElem[0];
  ResidualFaces(geometry, solver_container, numerics, config, iMesh, iStep,
                nMatchingInternalFacesLocalElem[0],
                nMatchingInternalFacesLocalElem[1], indResFaces);

  /* Complete the communication of the residuals. */
  Complete_MPI_ReverseCommunication();

  /* Create the final residual by summing up all contributions. */
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {

    su2double *resDOF = VecResDOFs.data() + nVar*i;
    for(unsigned long j=nEntriesResFaces[i]; j<nEntriesResFaces[i+1]; ++j) {
      const su2double *resFace = VecResFaces.data() + nVar*entriesResFaces[j];
      for(unsigned short k=0; k<nVar; ++k)
        resDOF[k] += resFace[k];
    }
  }
}

void CFEM_DG_EulerSolver::ResidualFaces(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                        CConfig *config, unsigned short iMesh, unsigned short iStep,
                                        const unsigned long indFaceBeg, const unsigned long indFaceEnd,
                                        unsigned long &indResFaces) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;
  su2double tick = 0.0;

  /*--- Loop over the requested range of matching faces. ---*/
  for(unsigned long l=indFaceBeg; l<indFaceEnd; ++l) {

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Compute the inviscid fluxes in the integration points of ---*/
    /*---         this matching face multiplied by the integration weight. ---*/
    /*------------------------------------------------------------------------*/

    /* Compute the inviscid fluxes in the integration points. */
    InviscidFluxesInternalMatchingFace(config, &matchingInternalFaces[l],
                                       solIntL, solIntR, fluxes, numerics);

    /* Get the number of integration points and integration weights
       from the standard element. */
    const unsigned short ind  = matchingInternalFaces[l].indStandardElement;
    const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();
    const su2double *weights  = standardMatchingFacesSol[ind].GetWeightsIntegration();

    /* Multiply the fluxes with the integration weight of the corresponding
       integration point. */
    for(unsigned short i=0; i<nInt; ++i) {
      su2double *flux = fluxes + i*nVar;

      for(unsigned short j=0; j<nVar; ++j)
        flux[j] *= weights[i];
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the contribution to the residuals from the       ---*/
    /*---         integration over this internal matching face.            ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the position in the residual array for side 0 of
       this face and update the corresponding counter. */
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    su2double *resFace0 = VecResFaces.data() + indResFaces*nVar;
    indResFaces        += nDOFsFace0;

    /* Get the correct form of the basis functions needed for the matrix
       multiplication to compute the residual. */
    const su2double *basisFaceTrans = standardMatchingFacesSol[ind].GetBasisFaceIntegrationTransposeSide0();

    /* Call the general function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nDOFsFace0, nVar, nInt, basisFaceTrans, fluxes, resFace0);
    config->GEMM_Tock(tick,"ResidualFaces1",nDOFsFace0, nVar, nInt);

    /* Easier storage of the position in the residual array for side 1 of
       this face and update the corresponding counter. */
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();
    su2double *resFace1 = VecResFaces.data() + indResFaces*nVar;
    indResFaces        += nDOFsFace1;

    /* Check if the number of DOFs on side 1 is equal to the number of DOFs
       of side 0. In that case the residual for side 1 is obtained by simply
       negating the data from side 0. */
    if(nDOFsFace1 == nDOFsFace0) {
      for(unsigned short i=0; i<(nVar*nDOFsFace1); ++i)
        resFace1[i] = -resFace0[i];
    }
    else {

      /*--- The number of DOFs and hence the polynomial degree of side 1 is
            different from side. Carry out the matrix multiplication to obtain
            the residual. Afterwards the residual is negated, because the
            normal is pointing into the adjacent element. ---*/
      basisFaceTrans = standardMatchingFacesSol[ind].GetBasisFaceIntegrationTransposeSide1();

      config->GEMM_Tick(&tick);
      DenseMatrixProduct(nDOFsFace1, nVar, nInt, basisFaceTrans, fluxes, resFace1);
      config->GEMM_Tock(tick,"ResidualFaces2",nDOFsFace1, nVar, nInt);

      for(unsigned short i=0; i<(nVar*nDOFsFace1); ++i)
        resFace1[i] = -resFace1[i];
    }
  }
}

void CFEM_DG_EulerSolver::InviscidFluxesInternalMatchingFace(
                                          CConfig                       *config,
                                          const CInternalFaceElementFEM *internalFace,
                                          su2double                     *solIntL,
                                          su2double                     *solIntR,
                                          su2double                     *fluxes,
                                          CNumerics                     *numerics) {

  /* Set the pointer solFace to fluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the face and the fluxes. */
  su2double *solFace = fluxes;
  su2double tick = 0.0;

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Interpolate the left state in the integration points of  ---*/
  /*---         the face.                                                ---*/
  /*------------------------------------------------------------------------*/

  /* Get the required information from the corresponding standard face. */
  const unsigned short ind        = internalFace->indStandardElement;
  const unsigned short nInt       = standardMatchingFacesSol[ind].GetNIntegration();
  const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
  const su2double     *basisFace0 = standardMatchingFacesSol[ind].GetBasisFaceIntegrationSide0();

  /*--- Store the solution of the DOFs of side 0 of the face in contiguous memory
        such that the function DenseMatrixProduct can be used to compute the left
        states in the integration points of the face, i.e. side 0. ---*/
  const unsigned long *DOFs = internalFace->DOFsSolFaceSide0.data();
  for(unsigned short i=0; i<nDOFsFace0; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
    su2double       *sol    = solFace + nVar*i;
    for(unsigned short j=0; j<nVar; ++j)
      sol[j] = solDOF[j];
  }

  config->GEMM_Tick(&tick);
  /* Compute the left states. Call the general function to
     carry out the matrix product. */
  DenseMatrixProduct(nInt, nVar, nDOFsFace0, basisFace0, solFace, solIntL);
  config->GEMM_Tock(tick, "InviscidFluxesInternalMatchingFace1", nInt, nVar, nDOFsFace0);

  /*------------------------------------------------------------------------*/
  /*--- Step 2: Interpolate the right state in the integration points of ---*/
  /*---         the face.                                                ---*/
  /*------------------------------------------------------------------------*/

  /* Get the required information from the corresponding standard face. */
  const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();
  const su2double     *basisFace1 = standardMatchingFacesSol[ind].GetBasisFaceIntegrationSide1();

  /*--- Store the solution of the DOFs of side 1 of the face in contiguous memory
        such that the function DenseMatrixProduct can be used to compute the right
        states in the integration points of the face, i.e. side 1. ---*/
  DOFs = internalFace->DOFsSolFaceSide1.data();
  for(unsigned short i=0; i<nDOFsFace1; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
    su2double       *sol    = solFace + nVar*i;
    for(unsigned short j=0; j<nVar; ++j)
      sol[j] = solDOF[j];
  }

  /* Compute the right states. Call the general function to
     carry out the matrix product. */
  config->GEMM_Tick(&tick);
  DenseMatrixProduct(nInt, nVar, nDOFsFace1, basisFace1, solFace, solIntR);
  config->GEMM_Tock(tick, "InviscidFluxesInternalMatchingFace2", nInt, nVar, nDOFsFace0);

  /*------------------------------------------------------------------------*/
  /*--- Step 3: Compute the fluxes in the integration points using the   ---*/
  /*---         approximate Riemann solver.                              ---*/
  /*------------------------------------------------------------------------*/

  /* General function to compute the fluxes in the integration points. */
  ComputeInviscidFluxesFace(config, nInt, internalFace->metricNormalsFace.data(),
                            solIntL, solIntR, fluxes, numerics);
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADER(unsigned short iTime, su2double weight) {

  /* Compute half the integration weight. The reason for doing this is that the
     given integration weight is based on the normalized interval [-1..1], i.e.
     a length of two. */
  const su2double halfWeight = 0.5*weight;

  /* Determine the case we have. */
  if(iTime == 0) {

    /*--- First time integration point. Initialize the total ADER residual
          to the spatial residual multiplied by halfWeight. ---*/
    for(unsigned long i=0; i<VecTotResDOFsADER.size(); ++i)
      VecTotResDOFsADER[i] = halfWeight*VecResDOFs[i];
  }
  else {

    /*--- Not the first integration point. The spatial residual, multiplied by
          halfWeight must be added to the total ADER residual. ---*/
    for(unsigned long i=0; i<VecTotResDOFsADER.size(); ++i)
      VecTotResDOFsADER[i] += halfWeight*VecResDOFs[i];
  }
}

void CFEM_DG_EulerSolver::MultiplyResidualByInverseMassMatrix(CConfig    *config,
                                                              const bool useADER) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *tmpRes = VecTmpMemory.data();
  su2double tick = 0.0;

  /*--- Set the reference to the correct residual. This depends
        whether or not the ADER scheme is used. ---*/
  vector<su2double> &VecRes = useADER ? VecTotResDOFsADER : VecResDOFs;

  /* Loop over the owned volume elements. */
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Easier storage of the residuals for this volume element. */
    su2double *res = VecRes.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Check whether a multiplication must be carried out with the inverse of
       the lumped mass matrix or the full mass matrix. Note that it is crucial
       that the test is performed with the lumpedMassMatrix and not with
       massMatrix. The reason is that for implicit time stepping schemes
       both arrays are in use. */
    if( volElem[l].lumpedMassMatrix.size() ) {

      /* Multiply the residual with the inverse of the lumped mass matrix. */
      for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {

        su2double *resDOF = res + nVar*i;
        su2double lMInv   = 1.0/volElem[l].lumpedMassMatrix[i];
        for(unsigned short k=0; k<nVar; ++k) resDOF[k] *= lMInv;
      }
    }
    else {

      /* Multiply the residual with the inverse of the mass matrix.
         Use the array tmpRes as temporary storage. */
      memcpy(tmpRes, res, nVar*volElem[l].nDOFsSol*sizeof(su2double));
      config->GEMM_Tick(&tick);
      DenseMatrixProduct(volElem[l].nDOFsSol, nVar, volElem[l].nDOFsSol,
                         volElem[l].invMassMatrix.data(), tmpRes, res);
      config->GEMM_Tock(tick, "MultiplyResidualByInverseMassMatrix",
                        volElem[l].nDOFsSol, nVar, volElem[l].nDOFsSol);
    }
  }
}

void CFEM_DG_EulerSolver::Pressure_Forces(CGeometry *geometry, CConfig *config) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solInt  = VecTmpMemory.data();
  su2double *solDOFs = solInt + nIntegrationMax*nVar;

  /*--- Get the information of the angle of attack, reference area, etc. ---*/
  const su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
  const su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
  const su2double RefAreaCoeff    = config->GetRefAreaCoeff();
  const su2double RefLengthMoment = config->GetRefLengthMoment();
  const su2double Gas_Constant    = config->GetGas_ConstantND();
  const su2double *Origin         = config->GetRefOriginMoment(0);
  const bool grid_movement        = config->GetGrid_Movement();

  /*--- Evaluate reference values for non-dimensionalization.
        For dynamic meshes, use the motion Mach number as a reference value
        for computing the force coefficients. Otherwise, use the freestream
        values, which is the standard convention. ---*/
  const su2double RefTemp     = Temperature_Inf;
  const su2double RefDensity  = Density_Inf;

  su2double RefVel2;
  if (grid_movement) {
    const su2double Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    const su2double Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      RefVel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }

  const su2double factor = 1.0/(0.5*RefDensity*RefAreaCoeff*RefVel2);

  /*-- Variables initialization ---*/
  Total_CD = 0.0; Total_CL = 0.0; Total_CSF = 0.0; Total_CEff = 0.0;
  Total_CMx = 0.0;   Total_CMy = 0.0;   Total_CMz = 0.0;
  Total_CFx = 0.0;   Total_CFy = 0.0;   Total_CFz = 0.0;

  AllBound_CD_Inv   = 0.0; AllBound_CL_Inv  = 0.0; AllBound_CSF_Inv = 0.0;
  AllBound_CMx_Inv  = 0.0; AllBound_CMy_Inv = 0.0; AllBound_CMz_Inv = 0.0;
  AllBound_CFx_Inv  = 0.0; AllBound_CFy_Inv = 0.0; AllBound_CFz_Inv = 0.0;
  AllBound_CEff_Inv = 0.0;

  for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring(); ++iMarker_Monitoring) {
    Surface_CL_Inv[iMarker_Monitoring]  = 0.0; Surface_CD_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring] = 0.0; Surface_CEff_Inv[iMarker_Monitoring] = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring] = 0.0; Surface_CFy_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring] = 0.0; Surface_CMx_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring] = 0.0; Surface_CMz_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CL[iMarker_Monitoring]      = 0.0; Surface_CD[iMarker_Monitoring]       = 0.0;
    Surface_CSF[iMarker_Monitoring]     = 0.0; Surface_CEff[iMarker_Monitoring]     = 0.0;
    Surface_CFx[iMarker_Monitoring]     = 0.0; Surface_CFy[iMarker_Monitoring]      = 0.0;
    Surface_CFz[iMarker_Monitoring]     = 0.0; Surface_CMx[iMarker_Monitoring]      = 0.0;
    Surface_CMy[iMarker_Monitoring]     = 0.0; Surface_CMz[iMarker_Monitoring]      = 0.0;
  }

  /*--- Loop over the Euler and Navier-Stokes markers ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /* Check if this boundary must be monitored. */
    const unsigned short Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if(Monitoring == YES) {

      /* Easier storage of the boundary condition. */
      const unsigned short Boundary = config->GetMarker_All_KindBC(iMarker);

      /*--- Obtain the origin for the moment computation for a particular marker ---*/
      for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                       ++iMarker_Monitoring) {
        string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        string Marker_Tag     = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }

      /* Check for a boundary for which the forces must be computed. */
      if((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
         (Boundary == ISOTHERMAL)) {

        /*--- Initialization for this marker ---*/
        CD_Inv[iMarker] = 0.0; CL_Inv[iMarker] = 0.0; CSF_Inv[iMarker] = 0.0;
        CMx_Inv[iMarker] = 0.0;   CMy_Inv[iMarker] = 0.0;   CMz_Inv[iMarker] = 0.0;
        CFx_Inv[iMarker] = 0.0;   CFy_Inv[iMarker] = 0.0;   CFz_Inv[iMarker] = 0.0;
        CEff_Inv[iMarker] = 0.0;

        su2double ForceInviscid[]  = {0.0, 0.0, 0.0};
        su2double MomentInviscid[] = {0.0, 0.0, 0.0};

        /* Easier storage of the boundary faces for this boundary marker. */
        const unsigned long      nSurfElem = boundaries[iMarker].surfElem.size();
        const CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

        /*--- Loop over the faces of this boundary. ---*/
        for(unsigned long l=0; l<nSurfElem; ++l) {

          /* Compute the states in the integration points of the face. */
          LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], solDOFs, solInt);

          /*--- Get the number of integration points and the integration
                weights from the corresponding standard element. ---*/
          const unsigned short ind   = surfElem[l].indStandardElement;
          const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
          const su2double *weights   = standardBoundaryFacesSol[ind].GetWeightsIntegration();

          /* Loop over the integration points of this surface element. */
          for(unsigned short i=0; i<nInt; ++i) {

            /* Easier storage of the solution, the normals and the coordinates
               for this integration point. */
            const su2double *sol     = solInt + i*nVar;
            const su2double *normals = surfElem[l].metricNormalsFace.data()
                                     + i*(nDim+1);
            const su2double *Coord   = surfElem[l].coorIntegrationPoints.data()
                                     + i*nDim;

            /*--- Compute the velocities and pressure in this integration point. ---*/
            const su2double DensityInv = 1.0/sol[0];
            su2double Velocity2 = 0.0;
            for(unsigned short iDim=0; iDim<nDim; ++iDim) {
              const su2double vel = sol[iDim+1]*DensityInv;
              Velocity2          += vel*vel;
            }

            su2double StaticEnergy = sol[nDim+1]*DensityInv - 0.5*Velocity2;

            FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
            const su2double Pressure = FluidModel->GetPressure();

            /*-- Compute the vector from the reference point to the integration
                 point and update the inviscid force. Note that the normal points
                 into the geometry, hence no minus sign. ---*/
            su2double MomentDist[] = {0.0, 0.0, 0.0}, Force[] = {0.0, 0.0, 0.0};
            const su2double forceMag = (Pressure - Pressure_Inf)*weights[i]
                                     * normals[nDim]*factor;

            for(unsigned short iDim=0; iDim<nDim; ++iDim) {
              MomentDist[iDim]     = Coord[iDim] - Origin[iDim];
              Force[iDim]          = forceMag*normals[iDim];
              ForceInviscid[iDim] += Force[iDim];
            }

            /*--- Update the inviscid moment. ---*/
            if (nDim == 3) {
              MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLengthMoment;
              MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
            }
            MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLengthMoment;
          }
        }

        /*--- Project forces and store the non-dimensional coefficients ---*/
        if(nDim == 2) {
          CD_Inv[iMarker]  =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
          CL_Inv[iMarker]  = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
          CEff_Inv[iMarker]   =  CL_Inv[iMarker] / (CD_Inv[iMarker]+EPS);
          CMz_Inv[iMarker]    =  MomentInviscid[2];
          CFx_Inv[iMarker]    =  ForceInviscid[0];
          CFy_Inv[iMarker]    =  ForceInviscid[1];
        }
        if(nDim == 3) {
          CD_Inv[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta)
                                  +  ForceInviscid[1]*sin(Beta)
                                  +  ForceInviscid[2]*sin(Alpha)*cos(Beta);
          CL_Inv[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
          CSF_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha)
                                  +  ForceInviscid[1]*cos(Beta)
                                  -  ForceInviscid[2]*sin(Beta)*sin(Alpha);
          CEff_Inv[iMarker]       =  CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
          CMx_Inv[iMarker]        =  MomentInviscid[0];
          CMy_Inv[iMarker]        =  MomentInviscid[1];
          CMz_Inv[iMarker]        =  MomentInviscid[2];
          CFx_Inv[iMarker]        =  ForceInviscid[0];
          CFy_Inv[iMarker]        =  ForceInviscid[1];
          CFz_Inv[iMarker]        =  ForceInviscid[2];
        }

        AllBound_CD_Inv    += CD_Inv[iMarker];
        AllBound_CL_Inv    += CL_Inv[iMarker];
        AllBound_CSF_Inv   += CSF_Inv[iMarker];
        AllBound_CEff_Inv   = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
        AllBound_CMx_Inv   += CMx_Inv[iMarker];
        AllBound_CMy_Inv   += CMy_Inv[iMarker];
        AllBound_CMz_Inv   += CMz_Inv[iMarker];
        AllBound_CFx_Inv   += CFx_Inv[iMarker];
        AllBound_CFy_Inv   += CFy_Inv[iMarker];
        AllBound_CFz_Inv   += CFz_Inv[iMarker];

        /*--- Compute the coefficients per surface ---*/
        for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                         ++iMarker_Monitoring) {
          string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Inv[iMarker_Monitoring]   += CL_Inv[iMarker];
            Surface_CD_Inv[iMarker_Monitoring]   += CD_Inv[iMarker];
            Surface_CSF_Inv[iMarker_Monitoring]  += CSF_Inv[iMarker];
            Surface_CEff_Inv[iMarker_Monitoring]  = Surface_CL_Inv[iMarker_Monitoring]
                                                  / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
            Surface_CFx_Inv[iMarker_Monitoring]  += CFx_Inv[iMarker];
            Surface_CFy_Inv[iMarker_Monitoring]  += CFy_Inv[iMarker];
            Surface_CFz_Inv[iMarker_Monitoring]  += CFz_Inv[iMarker];
            Surface_CMx_Inv[iMarker_Monitoring]  += CMx_Inv[iMarker];
            Surface_CMy_Inv[iMarker_Monitoring]  += CMy_Inv[iMarker];
            Surface_CMz_Inv[iMarker_Monitoring]  += CMz_Inv[iMarker];
          }
        }
      }
    }
  }

#ifdef HAVE_MPI

  /*--- Parallel mode. The data from all ranks must be gathered.
        Determine the size of the communication buffer. ---*/
  const unsigned long nCommSize = 9*config->GetnMarker_Monitoring() + 9;

  /*--- Define the communication buffers and store to local data in
        the local buffer. ---*/
  vector<su2double> locBuf(nCommSize), globBuf(nCommSize);

  unsigned long ii = 0;
  locBuf[ii++] = AllBound_CD_Inv;  locBuf[ii++] = AllBound_CL_Inv;
  locBuf[ii++] = AllBound_CSF_Inv; locBuf[ii++] = AllBound_CMx_Inv;
  locBuf[ii++] = AllBound_CMy_Inv; locBuf[ii++] = AllBound_CMz_Inv;
  locBuf[ii++] = AllBound_CFx_Inv; locBuf[ii++] = AllBound_CFy_Inv;
  locBuf[ii++] = AllBound_CFz_Inv;

  for(unsigned short i=0; i<config->GetnMarker_Monitoring(); ++i) {
    locBuf[ii++] = Surface_CL_Inv[i];  locBuf[ii++] = Surface_CD_Inv[i];
    locBuf[ii++] = Surface_CSF_Inv[i]; locBuf[ii++] = Surface_CFx_Inv[i];
    locBuf[ii++] = Surface_CFy_Inv[i]; locBuf[ii++] = Surface_CFz_Inv[i];
    locBuf[ii++] = Surface_CMx_Inv[i]; locBuf[ii++] = Surface_CMy_Inv[i];
    locBuf[ii++] = Surface_CMz_Inv[i];
  }

  /* Sum up all the data from all ranks. The result will be available on all ranks. */
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    SU2_MPI::Allreduce(locBuf.data(), globBuf.data(), nCommSize, MPI_DOUBLE,
                       MPI_SUM, MPI_COMM_WORLD);
  }

  /*--- Copy the data back from globBuf into the required variables. ---*/
  ii = 0;
  AllBound_CD_Inv  = globBuf[ii++]; AllBound_CL_Inv  = globBuf[ii++];
  AllBound_CSF_Inv = globBuf[ii++]; AllBound_CMx_Inv = globBuf[ii++];
  AllBound_CMy_Inv = globBuf[ii++]; AllBound_CMz_Inv = globBuf[ii++];
  AllBound_CFx_Inv = globBuf[ii++]; AllBound_CFy_Inv = globBuf[ii++];
  AllBound_CFz_Inv = globBuf[ii++];

  AllBound_CEff_Inv = AllBound_CL_Inv/(AllBound_CD_Inv + EPS);

  for(unsigned short i=0; i<config->GetnMarker_Monitoring(); ++i) {
    Surface_CL_Inv[i]  = globBuf[ii++]; Surface_CD_Inv[i]  = globBuf[ii++];
    Surface_CSF_Inv[i] = globBuf[ii++]; Surface_CFx_Inv[i] = globBuf[ii++];
    Surface_CFy_Inv[i] = globBuf[ii++]; Surface_CFz_Inv[i] = globBuf[ii++];
    Surface_CMx_Inv[i] = globBuf[ii++]; Surface_CMy_Inv[i] = globBuf[ii++];
    Surface_CMz_Inv[i] = globBuf[ii++];

    Surface_CEff_Inv[i] = Surface_CL_Inv[i]/(Surface_CD_Inv[i] + EPS);
  }

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
  Total_CD      = AllBound_CD_Inv;
  Total_CL      = AllBound_CL_Inv;
  Total_CSF = AllBound_CSF_Inv;
  Total_CEff       = Total_CL / (Total_CD + EPS);
  Total_CMx        = AllBound_CMx_Inv;
  Total_CMy        = AllBound_CMy_Inv;
  Total_CMz        = AllBound_CMz_Inv;
  Total_CFx        = AllBound_CFx_Inv;
  Total_CFy        = AllBound_CFy_Inv;
  Total_CFz        = AllBound_CFz_Inv;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                   ++iMarker_Monitoring) {
    Surface_CL[iMarker_Monitoring]      = Surface_CL_Inv[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]      = Surface_CD_Inv[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring] = Surface_CSF_Inv[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CL_Inv[iMarker_Monitoring]
                                           / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];
  }
}

void CFEM_DG_EulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                               CConfig *config, unsigned short iRKStep) {

  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

  for(unsigned short iVar=0; iVar<nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Store the coordinate of the first vertex of this element to give an
       indication for the location of the maximum residual. */
    const unsigned long ind = volElem[l].nodeIDsGrid[0];
    const su2double *coor   = meshPoints[ind].coor;

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar*volElem[l].offsetDOFsSolLocal;
    const su2double *res        = VecResDOFs.data()    + offset;
    const su2double *solDOFsOld = VecSolDOFsOld.data() + offset;
    su2double *solDOFs          = VecSolDOFs.data()    + offset;

    /* Loop over the DOFs for this element and update the solution and the L2 norm. */
    const su2double tmp = RK_AlphaCoeff*VecDeltaTime[l];

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      for(unsigned short iVar=0; iVar<nVar; ++iVar, ++i) {
        solDOFs[i] = solDOFsOld[i] - tmp*res[i];

        AddRes_RMS(iVar, res[i]*res[i]);
        AddRes_Max(iVar, fabs(res[i]), globalIndex, coor);
      }
    }
  }

  /*--- Compute the root mean square residual. Note that the SetResidual_RMS
        function cannot be used, because that is for the FV solver.    ---*/

#ifdef HAVE_MPI
  /*--- Parallel mode. The local L2 norms must be added to obtain the
        global value. Also check for divergence. ---*/
  int nProc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  vector<su2double> rbuf(nVar);

  /*--- Disable the reduce for the residual to avoid overhead if requested. ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    SU2_MPI::Allreduce(Residual_RMS, rbuf.data(), nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(unsigned short iVar=0; iVar<nVar; ++iVar) {

      if (rbuf[iVar] != rbuf[iVar]) {
        if (rank == MASTER_NODE)
          cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      SetRes_RMS(iVar, max(EPS*EPS, sqrt(rbuf[iVar]/nDOFsGlobal)));
    }
  }

#else
  /*--- Sequential mode. Check for a divergence of the solver and compute
        the L2-norm of the residuals. ---*/
  for(unsigned short iVar=0; iVar<nVar; ++iVar) {

    if(GetRes_RMS(iVar) != GetRes_RMS(iVar)) {
      cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
      exit(EXIT_FAILURE);
    }

    SetRes_RMS(iVar, max(EPS*EPS, sqrt(GetRes_RMS(iVar)/nDOFsGlobal)));
  }

#endif
}

void CFEM_DG_EulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                               CConfig *config, unsigned short iRKStep) {

  /*--- Hard-coded classical RK4 coefficients. Will be added to config. ---*/
  su2double RK_FuncCoeff[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
  su2double RK_TimeCoeff[4] = {0.5, 0.5, 1.0, 1.0};

  for(unsigned short iVar=0; iVar<nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Store the coordinate of the first vertex of this element to give an
     indication for the location of the maximum residual. */
    const unsigned long ind = volElem[l].nodeIDsGrid[0];
    const su2double *coor   = meshPoints[ind].coor;

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar*volElem[l].offsetDOFsSolLocal;
    const su2double *res        = VecResDOFs.data()    + offset;
    const su2double *solDOFsOld = VecSolDOFsOld.data() + offset;
    su2double *solDOFsNew       = VecSolDOFsNew.data() + offset;
    su2double *solDOFs          = VecSolDOFs.data()    + offset;

    /* Loop over the DOFs for this element and update the solution and the L2 norm. */
    const su2double tmp_time = -1.0*RK_TimeCoeff[iRKStep]*VecDeltaTime[l];
    const su2double tmp_func = -1.0*RK_FuncCoeff[iRKStep]*VecDeltaTime[l];

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      for(unsigned short iVar=0; iVar<nVar; ++iVar, ++i) {

        if (iRKStep < 3) {
          solDOFsNew[i] += tmp_func*res[i];
          solDOFs[i]     = solDOFsOld[i] + tmp_time*res[i];
        } else {
          solDOFs[i]     = solDOFsNew[i] + tmp_func*res[i];
        }

        AddRes_RMS(iVar, res[i]*res[i]);
        AddRes_Max(iVar, fabs(res[i]), globalIndex, coor);
      }
    }

  }

  /*--- Compute the root mean square residual. Note that the SetResidual_RMS
   function cannot be used, because that is for the FV solver.    ---*/

#ifdef HAVE_MPI
  /*--- Parallel mode. The local L2 norms must be added to obtain the
   global value. Also check for divergence. ---*/
  int nProc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  vector<su2double> rbuf(nVar);

  /*--- Disable the reduce for the residual to avoid overhead if requested. ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    SU2_MPI::Allreduce(Residual_RMS, rbuf.data(), nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(unsigned short iVar=0; iVar<nVar; ++iVar) {

      if (rbuf[iVar] != rbuf[iVar]) {
        if (rank == MASTER_NODE)
          cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      SetRes_RMS(iVar, max(EPS*EPS, sqrt(rbuf[iVar]/nDOFsGlobal)));
    }
  }

#else
  /*--- Sequential mode. Check for a divergence of the solver and compute
   the L2-norm of the residuals. ---*/
  for(unsigned short iVar=0; iVar<nVar; ++iVar) {

    if(GetRes_RMS(iVar) != GetRes_RMS(iVar)) {
      cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
      exit(EXIT_FAILURE);
    }

    SetRes_RMS(iVar, max(EPS*EPS, sqrt(GetRes_RMS(iVar)/nDOFsGlobal)));
  }

#endif
}

void CFEM_DG_EulerSolver::ADER_DG_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                            unsigned short iStep) {

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar*volElem[l].offsetDOFsSolLocal;
    const su2double *res        = VecTotResDOFsADER.data() + offset;
    const su2double *solDOFsOld = VecSolDOFsOld.data()     + offset;
    su2double *solDOFs          = VecSolDOFs.data()        + offset;

    /* Loop over the DOFs for this element and update the solution. */
    for(unsigned short i=0; i<(nVar*volElem[l].nDOFsSol); ++i)
      solDOFs[i] = solDOFsOld[i] - VecDeltaTime[l]*res[i];
  }
}

void CFEM_DG_EulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                        CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /*--- Apply the inviscid wall boundary conditions to compute the right
          state in the integration points. There are two options. Either the
          normal velocity is negated or the normal velocity is set to zero.
          Some experiments are needed to see which formulation gives better
          results. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution and the normals
         for this integration point. */
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;
      const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

      /* Compute the normal component of the momentum variables. */
      su2double rVn = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        rVn += UL[iDim+1]*normals[iDim];

      /* If the normal velocity must be mirrored instead of set to zero,
         the normal component that must be subtracted must be doubled. If the
         normal velocity must be set to zero, simply comment this line. */
      //rVn *= 2.0;

      /* Set the right state. The initial value of the total energy is the
         energy of the left state. */
      UR[0]      = UL[0];
      UR[nDim+1] = UL[nDim+1];
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        UR[iDim+1] = UL[iDim+1] - rVn*normals[iDim];

      /*--- Actually, only the internal energy of UR is equal to UL. If the
            kinetic energy differs for UL and UR, the difference must be
            subtracted from the total energy of UR to obtain the correct
            value. ---*/
      su2double DensityInv = 1.0/UL[0];
      su2double diffKin    = 0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double velL = DensityInv*UL[iDim];
        const su2double velR = DensityInv*UR[iDim];
        diffKin += velL*velL - velR*velR;
      }

      UR[nDim+1] -= 0.5*UL[0]*diffKin;
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /*--- Determine the number of integration points and set the right state
          to the free stream value. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {
      su2double *UR = solIntR + i*nVar;
      for(unsigned short j=0; j<nVar; ++j)
        UR[j] = ConsVarFreeStream[j];
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*ctc::nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*ctc::nVar;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Determine the number of integration points. */
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /* Make a distinction between two and three space dimensions
       in order to have the most efficient code. */
    switch( ctc::nDim ) {

      case 2: {

        /* 2D simulation. Loop over the integration points
           to apply the symmetry condition. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution and the normals
             for this integration point. */
          const su2double *UL      = solIntL + i*ctc::nVar;
                su2double *UR      = solIntR + i*ctc::nVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(ctc::nDim+1);

          /* Compute twice the normal component of the momentum variables. The
             factor 2 comes from the fact that the velocity must be mirrored. */
          const su2double rVn = 2.0*(UL[1]*normals[0] + UL[2]*normals[1]);

          /* Set the right state. Note that the total energy of the right state is
             identical to the left state, because the magnitude of the velocity
             remains the same. */
          UR[0] = UL[0];
          UR[1] = UL[1] - rVn*normals[0];
          UR[2] = UL[2] - rVn*normals[1];
          UR[3] = UL[3];
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /* 3D simulation. Loop over the integration points
           to apply the symmetry condition. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution and the normals
             for this integration point. */
          const su2double *UL      = solIntL + i*ctc::nVar;
                su2double *UR      = solIntR + i*ctc::nVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(ctc::nDim+1);

          /* Compute twice the normal component of the momentum variables. The
             factor 2 comes from the fact that the velocity must be mirrored. */
          const su2double rVn = 2.0*(UL[1]*normals[0] + UL[2]*normals[1] + UL[3]*normals[2]);

          /* Set the right state. Note that the total energy of the right state is
             identical to the left state, because the magnitude of the velocity
             remains the same. */
          UR[0] = UL[0];
          UR[1] = UL[1] - rVn*normals[0];
          UR[2] = UL[2] - rVn*normals[1];
          UR[3] = UL[3] - rVn*normals[2];
          UR[4] = UL[4];
        }

        break;
      }
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                   CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Retrieve the specified total conditions for this inlet. ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double P_Total   = config->GetInlet_Ptotal(Marker_Tag);
  su2double T_Total   = config->GetInlet_Ttotal(Marker_Tag);
  su2double *Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

  /*--- Non-dim. the inputs if necessary, and compute the total enthalpy. ---*/
  P_Total /= config->GetPressure_Ref();
  T_Total /= config->GetTemperature_Ref();

  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double H_Total      = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /*--- Apply the subsonic inlet boundary conditions to compute the right
          state in the integration points. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution and the normals
         for this integration point. */
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;
      const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

      /*--- Compute the normal velocity, the speed of sound squared and
            the pressure in this integration point. ---*/
      const su2double DensityInv = 1.0/UL[0];
      su2double VelocityNormal = 0.0, Velocity2 = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        const su2double vel = UL[iDim+1]*DensityInv;
        VelocityNormal     += vel*normals[iDim];
        Velocity2          += vel*vel;
      }

      su2double StaticEnergy = UL[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(UL[0], StaticEnergy);
      su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
      su2double Pressure    = FluidModel->GetPressure();

      /*--- Compute the Riemann invariant to be extrapolated. ---*/
      const su2double Riemann = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One + VelocityNormal;

      /*--- Total speed of sound. The expression below is also valid for variable cp,
            although a linearization around the value of the left state is performed. ---*/
      const su2double SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - DensityInv*(UL[nDim+1] + Pressure)
                                        +                  0.5*Velocity2) + SoundSpeed2;

      /*--- Dot product of normal and flow direction. This should be
            negative due to outward facing boundary normal convention. ---*/
      su2double alpha = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        alpha += normals[iDim]*Flow_Dir[iDim];

      /*--- Coefficients in the quadratic equation for the magnitude of the velocity. ---*/
      const su2double aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
      const su2double bb = -1.0*Gamma_Minus_One*alpha*Riemann;
      const su2double cc =  0.5*Gamma_Minus_One*Riemann*Riemann
                         -  2.0*SoundSpeed_Total2/Gamma_Minus_One;

      /*--- Solve the equation for the magnitude of the velocity. As this value
            must be positive and both aa and bb are positive (alpha is negative and
            Riemann is positive up till Mach = 5.0 or so, which is not really subsonic
            anymore), it is clear which of the two possible solutions must be taken.
            Some clipping is present, but this is normally not active. ---*/
      su2double dd      = bb*bb - 4.0*aa*cc;   dd      = sqrt(max(0.0, dd));
      su2double Vel_Mag = (-bb + dd)/(2.0*aa); Vel_Mag = max(0.0, Vel_Mag);
      Velocity2 = Vel_Mag*Vel_Mag;

      /*--- Compute speed of sound from total speed of sound eqn. ---*/
      SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

      /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
      su2double Mach2 = Velocity2/SoundSpeed2; Mach2 = min(1.0, Mach2);

      Velocity2 = Mach2*SoundSpeed2;
      Vel_Mag   = sqrt(Velocity2);

      /*--- Static temperature from the speed of sound relation. ---*/
      SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

      const su2double Temperature = SoundSpeed2/(Gamma*Gas_Constant);

      /*--- Static pressure using isentropic relation at a point ---*/
      Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);

      /*--- Density at the inlet from the gas law ---*/
      const su2double Density = Pressure/(Gas_Constant*Temperature);

      /*--- Store the conservative variables in UR. ---*/
      UR[0]      = Density;
      UR[nDim+1] = Pressure/Gamma_Minus_One + 0.5*Density*Velocity2;

      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        UR[iDim+1] = Density*Vel_Mag*Flow_Dir[iDim];
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Retrieve the specified back pressure for this outlet.
        Nondimensionalize, if necessary. ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double P_Exit = config->GetOutlet_Pressure(Marker_Tag);
  P_Exit = P_Exit/config->GetPressure_Ref();

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /*--- Apply the subsonic inlet boundary conditions to compute the right
          state in the integration points. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution and the normals
         for this integration point. */
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;
      const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

      /*--- Compute the normal velocity, the speed of sound squared and
            the pressure in this integration point. ---*/
      const su2double DensityInv = 1.0/UL[0];
      su2double VelocityNormal = 0.0, Velocity2 = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        const su2double vel = UL[iDim+1]*DensityInv;
        VelocityNormal     += vel*normals[iDim];
        Velocity2          += vel*vel;
      }

      su2double StaticEnergy = UL[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(UL[0], StaticEnergy);
      su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
      su2double Pressure    = FluidModel->GetPressure();

      /*--- Subsonic exit flow: there is one incoming characteristic,
            therefore one variable can be specified (back pressure) and is used
            to update the conservative variables. Compute the entropy and the
            acoustic Riemann variable. These invariants, as well as the
            tangential velocity components, are extrapolated. ---*/
      const su2double Entropy = Pressure*pow(DensityInv, Gamma);
      const su2double Riemann = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One + VelocityNormal;

      /*--- Compute the density and normal velocity of the right state. ---*/
      const su2double Density    = pow(P_Exit/Entropy,1.0/Gamma);
      const su2double SoundSpeed = sqrt(Gamma*P_Exit/Density);
      const su2double Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;

      /*--- Store the conservative variables in UR. ---*/
      Velocity2 = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        const su2double vel = UL[iDim+1]*DensityInv
                            + (Vn_Exit-VelocityNormal)*normals[iDim];
        Velocity2 += vel*vel;
        UR[iDim+1] = Density*vel;
      }

      UR[0]      = Density;
      UR[nDim+1] = Pressure/Gamma_Minus_One + 0.5*Density*Velocity2;
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /*--- Apply the subsonic inlet boundary conditions to compute the right
          state in the integration points. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {

#ifdef RINGLEB

      /* Ringleb case. Specify the exact solution for the right solution.
         First determine the pointer to the coordinates of this integration
         point and the pointer to the solution. Afterwards call the function
         RinglebSolution to do the actual job. */
      const su2double *coor = surfElem[l].coorIntegrationPoints + i*nDim;
            su2double *UR   = solIntR + i*nVar;

      RinglebSolution(coor, UR);

#else
      /* No compiler directive specified. Write an error message and exit. */
      int rank = MASTER_NODE;
#ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

      if (rank == MASTER_NODE) {
        cout << endl;
        cout << "In function CFEM_DG_EulerSolver::BC_Custom. " << endl;
        cout << "No or wrong compiler directive specified. This is necessary "
                "for customized boundary conditions." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif

#endif

    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::ResidualInviscidBoundaryFace(
                                      CConfig                  *config,
                                      CNumerics                *conv_numerics,
                                      const CSurfaceElementFEM *surfElem,
                                      const su2double          *solInt0,
                                      const su2double          *solInt1,
                                      su2double                *fluxes,
                                      su2double                *resFaces,
                                      unsigned long            &indResFaces) {

  /*--- Get the information from the standard element. ---*/
  const unsigned short ind   = surfElem->indStandardElement;
  const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();
  const su2double *weights   = standardBoundaryFacesSol[ind].GetWeightsIntegration();
  su2double tick = 0.0;

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Compute the fluxes in the integration points using the   ---*/
  /*---         approximate Riemann solver.                              ---*/
  /*------------------------------------------------------------------------*/

  /* General function to compute the fluxes in the integration points. */
  ComputeInviscidFluxesFace(config, nInt, surfElem->metricNormalsFace.data(),
                            solInt0, solInt1, fluxes, conv_numerics);

  /* Multiply the fluxes with the integration weight of the corresponding
     integration point. */
  for(unsigned short i=0; i<nInt; ++i) {
    su2double *flux = fluxes + i*nVar;

    for(unsigned short j=0; j<nVar; ++j)
      flux[j] *= weights[i];
  }

  /*------------------------------------------------------------------------*/
  /*--- Step 2: Compute the contribution to the residuals from the       ---*/
  /*---         integration over the surface element.                    ---*/
  /*------------------------------------------------------------------------*/

  /* Easier storage of the position in the residual array for this face
     and update the corresponding counter. */
  su2double *resFace = resFaces + indResFaces*nVar;
  indResFaces       += nDOFs;

  /* Get the correct form of the basis functions needed for the matrix
     multiplication to compute the residual. */
  const su2double *basisFaceTrans = standardBoundaryFacesSol[ind].GetBasisFaceIntegrationTranspose();

  /* Call the general function to carry out the matrix product. */
  config->GEMM_Tick(&tick);
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFaceTrans, fluxes, resFace);
  config->GEMM_Tock(tick, "ResidualInviscidBoundaryFace", nDOFs, nVar, nInt);
}

void CFEM_DG_EulerSolver::LeftStatesIntegrationPointsBoundaryFace(CConfig *config,
                                                                  const CSurfaceElementFEM *surfElem,
                                                                  su2double *solFace,
                                                                  su2double *solIntL) {

  /* Get the required information from the corresponding standard face. */
  const unsigned short ind   = surfElem->indStandardElement;
  const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();
  const su2double *basisFace = standardBoundaryFacesSol[ind].GetBasisFaceIntegration();
  su2double tick = 0.0;

  /* Easier storage of the DOFs of the face. */
  const unsigned long *DOFs = surfElem->DOFsSolFace.data();

  /* Copy the solution of the DOFs of the face such that it is contigious
     in memory. */
  for(unsigned short i=0; i<nDOFs; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
    su2double       *sol    = solFace + nVar*i;
    for(unsigned short j=0; j<nVar; ++j)
      sol[j] = solDOF[j];
  }

  /* Call the general function to carry out the matrix product. */
  config->GEMM_Tick(&tick);
  DenseMatrixProduct(nInt, nVar, nDOFs, basisFace, solFace, solIntL);
  config->GEMM_Tock(tick, "LeftStatesIntegrationPointsBoundaryFace", nInt, nVar, nDOFs);

}

void CFEM_DG_EulerSolver::ComputeInviscidFluxesFace(CConfig             *config,
                                                    const unsigned long nPoints,
                                                    const su2double     *normalsFace,
                                                    const su2double     *solL,
                                                    const su2double     *solR,
                                                    su2double           *fluxes,
                                                    CNumerics           *numerics) {

  /* Easier storage of the specific heat ratio. */
  const su2double gm1 = Gamma_Minus_One;

  /* Make a distinction between the several Riemann solvers. */
  switch( config->GetRiemann_Solver_FEM() ) {

    case ROE: {

      /* Roe's approximate Riemann solver. Easier storage of the cut off
         value for the entropy correction. */
      const su2double Delta = config->GetEntropyFix_Coeff();

      /* Make a distinction between two and three space dimensions
         in order to have the most efficient code. */
      switch( ctc::nDim ) {

        case 2: {

          /* Two dimensional simulation. Loop over the number of points. */
          for(unsigned long i=0; i<nPoints; ++i) {
        
            /* Easier storage of the left and right solution, the face normals and
               the flux vector for this point. */
            const su2double *UL   = solL + i*ctc::nVar;
            const su2double *UR   = solR + i*ctc::nVar;
            const su2double *norm = normalsFace + i*(ctc::nDim+1);
                  su2double *flux = fluxes + i*ctc::nVar;

            const su2double nx = norm[0], ny = norm[1], area = norm[2];

            /*--- compute the primitive variables of the left and right state. ---*/
            su2double tmp       = 1.0/UL[0];
            const su2double vxL = tmp*UL[1];
            const su2double vyL = tmp*UL[2];
            const su2double pL  = gm1*(UL[3] - 0.5*(vxL*UL[1] + vyL*UL[2]));

            tmp                 = 1.0/UR[0];
            const su2double vxR = tmp*UR[1];
            const su2double vyR = tmp*UR[2];
            const su2double pR  = gm1*(UR[3] - 0.5*(vxR*UR[1] + vyR*UR[2]));

            /*--- Compute the difference of the conservative mean flow variables. ---*/
            const su2double dr  = UR[0] - UL[0];
            const su2double dru = UR[1] - UL[1];
            const su2double drv = UR[2] - UL[2];
            const su2double drE = UR[3] - UL[3];

            /*--- Compute the Roe average state. ---*/
            const su2double zL = sqrt(UL[0]);
            const su2double zR = sqrt(UR[0]);
            tmp                = 1.0/(zL + zR);

            const su2double rHL = UL[3] + pL;
            const su2double rHR = UR[3] + pR;

            const su2double uAvg = tmp*(zL*vxL + zR*vxR);
            const su2double vAvg = tmp*(zL*vyL + zR*vyR);
            const su2double HAvg = tmp*(rHL/zL + rHR/zR);

            /*--- Compute from the Roe average state some variables, which occur
                  quite often in the matrix vector product to be computed. ---*/
            const su2double alphaAvg = 0.5*(uAvg*uAvg + vAvg*vAvg);
            tmp                      = gm1*(HAvg - alphaAvg);
            const su2double a2Avg    = fabs(tmp);
            const su2double aAvg     = sqrt(a2Avg);
            const su2double vnAvg    = uAvg*nx + vAvg*ny;
            const su2double ovaAvg   = 1.0/aAvg;
            const su2double ova2Avg  = 1.0/a2Avg;

            /*--- Compute the absolute values of the three eigenvalues, apply the
                  entropy correction and multiply by the area to obtain the correct
                  values for the dissipation term. ---*/
            su2double lam1 = fabs(vnAvg + aAvg);
            su2double lam2 = fabs(vnAvg - aAvg);
            su2double lam3 = fabs(vnAvg);

            tmp  = Delta*max(lam1, lam2);
            lam1 = max(lam1, tmp);
            lam2 = max(lam2, tmp);
            lam3 = max(lam3, tmp);

            lam1 *= area;
            lam2 *= area;
            lam3 *= area;

            /* Hack to get a kind of Lax-Friedrichs flux. */
            //lam1 = std::max(lam1, lam2);
            //lam3 = lam2 = lam1;
            /* End hack. */

            /*--- Some abbreviations, which occur quite often in the dissipation terms. ---*/
            const su2double abv1 = 0.5*(lam1 + lam2);
            const su2double abv2 = 0.5*(lam1 - lam2);
            const su2double abv3 = abv1 - lam3;

            const su2double abv4 = gm1*(alphaAvg*dr - uAvg*dru - vAvg*drv + drE);
            const su2double abv5 = nx*dru + ny*drv - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            /*--- Compute the Roe flux vector, which is 0.5*(FL + FR - |A|(UR-UL)). ---*/
            const su2double vnL  = area*(vxL*nx + vyL*ny);
            const su2double vnR  = area*(vxR*nx + vyR*ny);
            const su2double rvnL = UL[0]*vnL;
            const su2double rvnR = UR[0]*vnR;
            const su2double pa   = area*(pL + pR);

            flux[0] = 0.5*(rvnL + rvnR - (lam3*dr + abv6));
            flux[1] = 0.5*(rvnL*vxL + rvnR*vxR + pa*nx - (lam3*dru + uAvg*abv6 + nx*abv7));
            flux[2] = 0.5*(rvnL*vyL + rvnR*vyR + pa*ny - (lam3*drv + vAvg*abv6 + ny*abv7));
            flux[3] = 0.5*( vnL*rHL +  vnR*rHR - (lam3*drE + HAvg*abv6 + vnAvg*abv7));
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /* Three dimensional simulation. Loop over the number of points. */
          for(unsigned long i=0; i<nPoints; ++i) {

            /* Easier storage of the left and right solution, the face normals and
               the flux vector for this point. */
            const su2double *UL   = solL + i*ctc::nVar;
            const su2double *UR   = solR + i*ctc::nVar;
            const su2double *norm = normalsFace + i*(ctc::nDim+1);
                  su2double *flux = fluxes + i*ctc::nVar;

            const su2double nx = norm[0], ny = norm[1], nz = norm[2], area = norm[3];

            /*--- Compute the primitive variables of the left and right state. ---*/
            su2double tmp       = 1.0/UL[0];
            const su2double vxL = tmp*UL[1];
            const su2double vyL = tmp*UL[2];
            const su2double vzL = tmp*UL[3];
            const su2double pL  = gm1*(UL[4] - 0.5*(vxL*UL[1] + vyL*UL[2] + vzL*UL[3]));

            tmp                 = 1.0/UR[0];
            const su2double vxR = tmp*UR[1];
            const su2double vyR = tmp*UR[2];
            const su2double vzR = tmp*UR[3];
            const su2double pR  = gm1*(UR[4] - 0.5*(vxR*UR[1] + vyR*UR[2] + vzR*UR[3]));

            /*--- Compute the difference of the conservative mean flow variables. ---*/
            const su2double dr  = UR[0] - UL[0];
            const su2double dru = UR[1] - UL[1];
            const su2double drv = UR[2] - UL[2];
            const su2double drw = UR[3] - UL[3];
            const su2double drE = UR[4] - UL[4];

            /*--- Compute the Roe average state. ---*/
            const su2double zL = sqrt(UL[0]);
            const su2double zR = sqrt(UR[0]);
            tmp                = 1.0/(zL + zR);

            const su2double rHL = UL[4] + pL;
            const su2double rHR = UR[4] + pR;

            const su2double uAvg = tmp*(zL*vxL + zR*vxR);
            const su2double vAvg = tmp*(zL*vyL + zR*vyR);
            const su2double wAvg = tmp*(zL*vzL + zR*vzR);
            const su2double HAvg = tmp*(rHL/zL + rHR/zR);

            /*--- Compute from the Roe average state some variables, which occur
                  quite often in the matrix vector product to be computed. ---*/
            const su2double alphaAvg = 0.5*(uAvg*uAvg + vAvg*vAvg + wAvg*wAvg);
            tmp                      = gm1*(HAvg - alphaAvg);
            const su2double a2Avg    = fabs(tmp);
            const su2double aAvg     = sqrt(a2Avg);
            const su2double vnAvg    = uAvg*nx + vAvg*ny + wAvg*nz;
            const su2double ovaAvg   = 1.0/aAvg;
            const su2double ova2Avg  = 1.0/a2Avg;

            /*--- Compute the absolute values of the three eigenvalues, apply the
                  entropy correction and multiply by the area to obtain the correct
                  values for the dissipation term. ---*/
            su2double lam1 = fabs(vnAvg + aAvg);
            su2double lam2 = fabs(vnAvg - aAvg);
            su2double lam3 = fabs(vnAvg);

            tmp  = Delta*max(lam1, lam2);
            lam1 = max(lam1, tmp);
            lam2 = max(lam2, tmp);
            lam3 = max(lam3, tmp);

            lam1 *= area;
            lam2 *= area;
            lam3 *= area;

            /* Hack to get a kind of Lax-Friedrichs flux. */
            //lam1 = std::max(lam1, lam2);
            //lam3 = lam2 = lam1;
            /* End hack. */

            // Some abbreviations, which occur quite often in the dissipation terms.
            const su2double abv1 = 0.5*(lam1 + lam2);
            const su2double abv2 = 0.5*(lam1 - lam2);
            const su2double abv3 = abv1 - lam3;

            const su2double abv4 = gm1*(alphaAvg*dr - uAvg*dru - vAvg*drv -wAvg*drw + drE);
            const su2double abv5 = nx*dru + ny*drv + nz*drw - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            // Compute the Roe flux vector, which is 0.5*(FL + FR - |A|(UR-UL)).
            const su2double vnL  = area*(vxL*nx + vyL*ny + vzL*nz);
            const su2double vnR  = area*(vxR*nx + vyR*ny + vzR*nz);
            const su2double rvnL = UL[0]*vnL;
            const su2double rvnR = UR[0]*vnR;
            const su2double pa   = area*(pL + pR);

            flux[0] = 0.5*(rvnL + rvnR - (lam3*dr + abv6));
            flux[1] = 0.5*(rvnL*vxL + rvnR*vxR + pa*nx - (lam3*dru + uAvg*abv6 + nx*abv7));
            flux[2] = 0.5*(rvnL*vyL + rvnR*vyR + pa*ny - (lam3*drv + vAvg*abv6 + ny*abv7));
            flux[3] = 0.5*(rvnL*vzL + rvnR*vzR + pa*nz - (lam3*drw + wAvg*abv6 + nz*abv7));
            flux[4] = 0.5*( vnL*rHL +  vnR*rHR - (lam3*drE + HAvg*abv6 + vnAvg*abv7));
          }

          break;
        }
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    default: {

      /* Riemann solver not explicitly implemented. Fall back to the
         implementation via numerics. This is not efficient. */

      /*--- Data for loading into the CNumerics Riemann solvers.
       This is temporary and not efficient.. just replicating exactly
       the arrays we typically have in order to avoid bugs. We can
       probably be more clever with pointers, etc. ---*/
      su2double Normal[3];
      su2double Prim_L[8];
      su2double Prim_R[8];
  
      Jacobian_i = new su2double*[nVar];
      Jacobian_j = new su2double*[nVar];
      for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
        Jacobian_i[iVar] = new su2double[nVar];
        Jacobian_j[iVar] = new su2double[nVar];
      }
  
      /*--- Loop over the number of points. ---*/
  
      for(unsigned long i=0; i<nPoints; ++i) {
    
        /* Easier storage of the left and right solution, the face normals and
         the flux vector for this point. */
        const su2double *UL   = solL + i*nVar;
        const su2double *UR   = solR + i*nVar;
        const su2double *norm = normalsFace + i*(nDim+1);
        su2double       *flux = fluxes + i*nVar;
    
        /*--- Store and load the normal into numerics. ---*/
        for (unsigned short iDim = 0; iDim < nDim; ++iDim)
          Normal[iDim] = norm[iDim]*norm[nDim];
        numerics->SetNormal(Normal);
    
        /*--- Prepare the primitive states for the numerics class. Note 
         that for the FV solver, we have the following primitive
         variable ordering: Compressible flow, primitive variables nDim+5,
         (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
    
        /*--- Left primitive state ---*/
        Prim_L[0] = 0.0;                                        // Temperature (unused)
        Prim_L[nDim+1] = gm1*UL[nVar-1];
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Prim_L[iDim+1]  = UL[iDim+1]/UL[0];                   // Velocities
          Prim_L[nDim+1] -= gm1*0.5*Prim_L[iDim+1]*UL[iDim+1];  // Pressure
        }
        Prim_L[nDim+2] = UL[0];                                 // Density
        Prim_L[nDim+3] = (UL[nVar-1] + Prim_L[nDim+1]) / UL[0]; // Enthalpy
    
        /*--- Right primitive state ---*/
        Prim_R[0] = 0.0;                                        // Temperature (unused)
        Prim_R[nDim+1] = gm1*UR[nVar-1];
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Prim_R[iDim+1]  = UR[iDim+1]/UR[0];                   // Velocities
          Prim_R[nDim+1] -= gm1*0.5*Prim_R[iDim+1]*UR[iDim+1];  // Pressure
        }
        Prim_R[nDim+2] = UR[0];                                 // Density
        Prim_R[nDim+3] = (UR[nVar-1] + Prim_R[nDim+1]) / UR[0]; // Enthalpy
    
        /*--- Load the primitive states into the numerics class. ---*/
        numerics->SetPrimitive(Prim_L, Prim_R);
    
        /*--- Now simply call the ComputeResidual() function to calculate
         the flux using the chosen approximate Riemann solver. Note that
         the Jacobian arrays here are just dummies for now (no implicit). ---*/
        numerics->ComputeResidual(flux, Jacobian_i, Jacobian_j, config);
      }
  
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        delete [] Jacobian_i[iVar];
        delete [] Jacobian_j[iVar];
      }
      delete [] Jacobian_i;
      delete [] Jacobian_j;

      Jacobian_i = NULL;
      Jacobian_j = NULL;
    }
  }
}

#ifdef RINGLEB

void CFEM_DG_EulerSolver::RinglebSolution(const su2double *coor,
                                                su2double *sol) {

  /* Compute several expononts involving Gamma. */
  const su2double gm1     = Gamma_Minus_One;
  const su2double tovgm1  = 2.0/gm1;
  const su2double tgovgm1 = Gamma*tovgm1;

  /* Easier storage of the coordinates and abbreviate y*y. */
  const su2double x  = coor[0], y = coor[1];
  const su2double y2 = y*y;

  /* Initial guess for q (velocity magnitude) and k (streamline parameter). */
  su2double k = 1.2;
  su2double q = 1.0;

  /* Newton algorithm to solve for the variables q and k for the given x and y. */
  const int iterMax = 500;
  su2double duMaxPrev = 10.0;

  int iter;
  for(iter=0; iter<iterMax; ++iter) {

    /* Compute the speed of sound, the density, the parameter JJ
       and its derivatives w.r.t. q. */
    const su2double a   = sqrt(1.0 - 0.5*gm1*q*q);
    const su2double rho = pow(a,tovgm1);
    const su2double JJ  = 1.0/a + 1.0/(3.0*a*a*a) + 1.0/(5.0*a*a*a*a*a)
                        - 0.5*log((1.0+a)/(1.0-a));

    const su2double dadq   = -0.5*gm1*q/a;
    const su2double drhodq =  2.0*rho*dadq/(gm1*a);
    const su2double dJJdq  =  dadq/(pow(a,6)*(a*a-1.0));

    /* Determine the values of the nonlinear equations to solve
       and its corresponding Jacobian matrix. */
    const su2double y2c = (k*k - q*q)/(k*k*k*k*rho*rho*q*q);
    const su2double f[] = {(2.0/(k*k) - 1.0/(q*q))/(2.0*rho) - 0.5*JJ - x,
                           y2c - y2};
    su2double Jac[2][2];
    Jac[0][0] = -(1.0/(k*k) - 0.50/(q*q))*drhodq/(rho*rho)
              + 1.0/(rho*q*q*q) - 0.5*dJJdq;
    Jac[0][1] = -2.0/(rho*k*k*k);
    Jac[1][0] = -2.0/(k*k*rho*rho*q*q*q) - 2.0*y2c*drhodq/rho;
    Jac[1][1] = (4.0*q*q - 2.0*k*k)/(k*k*k*k*k*rho*rho*q*q);

    /* Determine the update dU. */
    const su2double det  = Jac[0][0]*Jac[1][1] - Jac[0][1]*Jac[1][0];
    const su2double dU[] = {(f[0]*Jac[1][1] - f[1]*Jac[0][1])/det,
                            (f[1]*Jac[0][0] - f[0]*Jac[1][0])/det};

    /* Determine the underrelaxation coefficient alp. */
    const su2double dUMax = max(fabs(dU[0]), fabs(dU[1]));
    su2double alp = 1.0;
    if(     dUMax > 1.0) alp = 0.04;
    else if(dUMax > 0.1) alp = 0.2;

    /* Update q and k. */
    q -= alp*dU[0];
    k -= alp*dU[1];

    /* Convergence check, which is independent of the precision used. */
    if((dUMax < 1.e-3) && (dUMax >= duMaxPrev)) break;
    duMaxPrev = dUMax;
  }

  /* Check if the Newton algorithm actually converged. */
  if(iter == iterMax) {
    cout << "In function CFEM_DG_EulerSolver::RinglebSolution: "
         << "Newton algorithm did not converge." << endl << flush;
    exit(1);
  }

  /* Compute the speed of sound, density and pressure. */
  const su2double a   = sqrt(1.0 - 0.5*gm1*q*q);
  const su2double rho = pow(a,tovgm1);
  const su2double p   = pow(a,tgovgm1)/Gamma;

  /* Initialize the velocity direction to (0,-1), which corresponds to the
     direction on the symmetry line. */
  su2double velDir[] = {0.0, -1.0};

  /* Check for a point away from the symmetry line. */
  if(fabs(y) > 1.e-8) {

    /*--- Determine the derivatives of x and y w.r.t. q, such that the direction
          of the streamline can be determined. This direction is identical to the
          velocity direction. ---*/
    const su2double dadq   = -0.5*gm1*q/a;
    const su2double drhodq =  2.0*rho*dadq/(gm1*a);
    const su2double dJJdq  =  dadq/(pow(a,6)*(a*a-1.0));

    const su2double dxdq = -(1.0/(k*k) - 0.5/(q*q))*drhodq/(rho*rho)
                         +   1.0/(rho*q*q*q) - 0.5*dJJdq;

    const su2double dydq = y > 0.0 ? -y*drhodq/rho - 1.0/(rho*q*q*sqrt(k*k - q*q))
                                   : -y*drhodq/rho + 1.0/(rho*q*q*sqrt(k*k - q*q));

    const su2double vecLen = sqrt(dxdq*dxdq + dydq*dydq);

    velDir[0] = dxdq/vecLen; velDir[1] = dydq/vecLen;
    if(velDir[1] > 0.0){velDir[0] = -velDir[0]; velDir[1] = -velDir[1];}
  }

  /* Compute the conservative variables. Note that both 2D and 3D
     cases are treated correctly. */
  sol[0]      = rho;
  sol[1]      = rho*q*velDir[0];
  sol[2]      = rho*q*velDir[1];
  sol[3]      = 0.0;
  sol[nVar-1] = p/gm1 + 0.5*rho*q*q;
}

#endif

CFEM_DG_NSSolver::CFEM_DG_NSSolver(void) : CFEM_DG_EulerSolver() {

  /*--- Basic array initialization ---*/
  CD_Visc  = NULL; CL_Visc  = NULL; CSF_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL; CMy_Visc = NULL; CMz_Visc = NULL;
  CFx_Visc = NULL; CFy_Visc = NULL; CFz_Visc = NULL;

  ForceViscous = NULL; MomentViscous = NULL; CSkinFriction = NULL;

  /*--- Surface-based array initialization ---*/
  Surface_CL_Visc  = NULL; Surface_CD_Visc  = NULL; Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL; Surface_CFy_Visc = NULL; Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL; Surface_CMy_Visc = NULL; Surface_CMz_Visc = NULL;
  MaxHeatFlux_Visc = NULL; Heat_Visc = NULL;

  /*--- Set the SGS model to NULL and indicate that no SGS model is used. ---*/
  SGSModel     = NULL;
  SGSModelUsed = false;
}

CFEM_DG_NSSolver::CFEM_DG_NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
 : CFEM_DG_EulerSolver(geometry, config, iMesh) {

  /*--- Array initialization ---*/
  CD_Visc = NULL;  CL_Visc = NULL;  CSF_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL; CMy_Visc = NULL; CMz_Visc = NULL;
  CFx_Visc = NULL; CFy_Visc = NULL; CFz_Visc = NULL;

  Surface_CL_Visc  = NULL; Surface_CD_Visc = NULL;  Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL; Surface_CFy_Visc = NULL; Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL; Surface_CMy_Visc = NULL; Surface_CMz_Visc = NULL;
  MaxHeatFlux_Visc = NULL; Heat_Visc = NULL;

  ForceViscous = NULL;  MomentViscous = NULL;
  CSkinFriction = NULL; Cauchy_Serie = NULL;

  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/

  //LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  //LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Non dimensional coefficients ---*/
  ForceViscous  = new su2double[3];
  MomentViscous = new su2double[3];
  CD_Visc       = new su2double[nMarker];
  CL_Visc       = new su2double[nMarker];
  CSF_Visc      = new su2double[nMarker];
  CMx_Visc      = new su2double[nMarker];
  CMy_Visc      = new su2double[nMarker];
  CMz_Visc      = new su2double[nMarker];
  CEff_Visc     = new su2double[nMarker];
  CFx_Visc      = new su2double[nMarker];
  CFy_Visc      = new su2double[nMarker];
  CFz_Visc      = new su2double[nMarker];

  Surface_CL_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc  = new su2double[config->GetnMarker_Monitoring()];

  Heat_Visc        = new su2double[nMarker];
  MaxHeatFlux_Visc = new su2double[nMarker];

  /*--- Init total coefficients ---*/

  Total_CD   = 0.0; Total_CL  = 0.0; Total_CSF = 0.0;
  Total_CMx  = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
  Total_CEff = 0.0;
  Total_CFx  = 0.0; Total_CFy = 0.0; Total_CFz = 0.0;

  /*--- Read farfield conditions from config ---*/

  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  Prandtl_Lam   = config->GetPrandtl_Lam();
  Prandtl_Turb  = config->GetPrandtl_Turb();
  Tke_Inf       = config->GetTke_FreeStreamND();

  /*--- Set the SGS model in case an LES simulation is carried out ---*/

  if(config->GetKind_Solver() == FEM_LES) {

    /* Make a distinction between the SGS models used and set SGSModel and
       SGSModelUsed accordingly. */
    switch( config->GetKind_SGS_Model() ) {

      case IMPLICIT_LES:
        SGSModel     = NULL;
        SGSModelUsed = false;
        break;

      case SMAGORINSKY:
        SGSModel     = new CSmagorinskyModel;
        SGSModelUsed = true;
        break;

      case WALE:
        SGSModel     = new CWALEModel;
        SGSModelUsed = true;
        break;

      default:
        cout << "Unknown SGS model encountered" << endl;

#ifndef HAVE_MPI
        exit(EXIT_FAILURE);
#else
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#endif
    }
  }
  else {

    /* No LES, so no SGS model needed.
       Set the pointer to NULL and the boolean to false. */
    SGSModel     = NULL;
    SGSModelUsed = false;
  }
}

CFEM_DG_NSSolver::~CFEM_DG_NSSolver(void) {
  unsigned short iMarker;

  if (CD_Visc != NULL)       delete [] CD_Visc;
  if (CL_Visc != NULL)       delete [] CL_Visc;
  if (CSF_Visc != NULL)      delete [] CSF_Visc;
  if (CMx_Visc != NULL)      delete [] CMx_Visc;
  if (CMy_Visc != NULL)      delete [] CMy_Visc;
  if (CMz_Visc != NULL)      delete [] CMz_Visc;
  if (CFx_Visc != NULL)      delete [] CFx_Visc;
  if (CFy_Visc != NULL)      delete [] CFy_Visc;
  if (CFz_Visc != NULL)      delete [] CFz_Visc;
  if (CEff_Visc != NULL)     delete [] CEff_Visc;
  if (ForceViscous != NULL)  delete [] ForceViscous;
  if (MomentViscous != NULL) delete [] MomentViscous;


  if (Surface_CL_Visc != NULL)   delete [] Surface_CL_Visc;
  if (Surface_CD_Visc != NULL)   delete [] Surface_CD_Visc;
  if (Surface_CSF_Visc != NULL)  delete [] Surface_CSF_Visc;
  if (Surface_CEff_Visc != NULL) delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != NULL)  delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != NULL)  delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != NULL)  delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != NULL)  delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)  delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)  delete [] Surface_CMz_Visc;

  if (Heat_Visc        != NULL)  delete [] Heat_Visc;
  if (MaxHeatFlux_Visc != NULL)  delete [] MaxHeatFlux_Visc;

  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;

  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }

  if( SGSModel ) delete SGSModel;
}

void CFEM_DG_NSSolver::Friction_Forces(CGeometry *geometry, CConfig *config) {

  /*--------------------------------------------------------------------------*/
  /*--- The assumption is made that the eddy viscosity is zero for viscous ---*/
  /*--- walls. This is true for an integration to the wall, but is this    ---*/
  /*--- also true when wall functions are used?                            ---*/
  /*--------------------------------------------------------------------------*/

  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux_Lam = Gamma/Prandtl_Lam;
  su2double tick = 0.0;

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solInt     = VecTmpMemory.data();
  su2double *gradSolInt = solInt     + nIntegrationMax*nVar;
  su2double *solDOFs    = gradSolInt + nIntegrationMax*nVar*nDim;

  /*--- Get the information of the angle of attack, reference area, etc. ---*/
  const su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
  const su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
  const su2double RefAreaCoeff    = config->GetRefAreaCoeff();
  const su2double RefLengthMoment = config->GetRefLengthMoment();
  const su2double Gas_Constant    = config->GetGas_ConstantND();
  const su2double *Origin         = config->GetRefOriginMoment(0);
  const bool grid_movement        = config->GetGrid_Movement();

  /*--- Evaluate reference values for non-dimensionalization.
        For dynamic meshes, use the motion Mach number as a reference value
        for computing the force coefficients. Otherwise, use the freestream
        values, which is the standard convention. ---*/
  const su2double RefTemp     = Temperature_Inf;
  const su2double RefDensity  = Density_Inf;

  su2double RefVel2;
  if (grid_movement) {
    const su2double Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    const su2double Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      RefVel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }

  const su2double factor = 1.0/(0.5*RefDensity*RefAreaCoeff*RefVel2);

  /*--- Variables initialization ---*/
  AllBound_CD_Visc = 0.0;  AllBound_CL_Visc  = 0.0; AllBound_CSF_Visc = 0.0;
  AllBound_CMx_Visc = 0.0; AllBound_CMy_Visc = 0.0; AllBound_CMz_Visc = 0.0;
  AllBound_CFx_Visc = 0.0; AllBound_CFy_Visc = 0.0; AllBound_CFz_Visc = 0.0;

  AllBound_HeatFlux_Visc = 0.0; AllBound_MaxHeatFlux_Visc = 0.0; AllBound_CEff_Visc = 0.0;

  for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring(); ++iMarker_Monitoring) {
    Surface_CL_Visc[iMarker_Monitoring]  = 0.0; Surface_CD_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring] = 0.0; Surface_CEff_Visc[iMarker_Monitoring] = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring] = 0.0; Surface_CFy_Visc[iMarker_Monitoring]  = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring] = 0.0; Surface_CMx_Visc[iMarker_Monitoring]  = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring] = 0.0; Surface_CMz_Visc[iMarker_Monitoring]  = 0.0;
  }

  /*--- Loop over the Navier-Stokes markers ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /* Check if this boundary must be monitored. */
    const unsigned short Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if(Monitoring == YES) {

      /* Easier storage of the boundary condition. */
      const unsigned short Boundary = config->GetMarker_All_KindBC(iMarker);

      /*--- Obtain the origin for the moment computation for a particular marker ---*/
      for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                       ++iMarker_Monitoring) {
        string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        string Marker_Tag     = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }

      /* Check for a boundary for which the viscous forces must be computed. */
      if((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {

        /*--- Forces initialization at each Marker ---*/
        CD_Visc[iMarker]  = 0.0; CL_Visc[iMarker]  = 0.0; CSF_Visc[iMarker] = 0.0;
        CMx_Visc[iMarker] = 0.0; CMy_Visc[iMarker] = 0.0; CMz_Visc[iMarker] = 0.0;
        CFx_Visc[iMarker] = 0.0; CFy_Visc[iMarker] = 0.0; CFz_Visc[iMarker] = 0.0;

        Heat_Visc[iMarker]  = 0.0; MaxHeatFlux_Visc[iMarker] = 0.0; CEff_Visc[iMarker]       = 0.0;

        for(unsigned short iDim=0; iDim<nDim; ++iDim)
          ForceViscous[iDim] = MomentViscous[iDim] = 0.0;

        /* Easier storage of the boundary faces for this boundary marker. */
        const unsigned long      nSurfElem = boundaries[iMarker].surfElem.size();
        const CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

        /*--- Loop over the faces of this boundary. ---*/
        for(unsigned long l=0; l<nSurfElem; ++l) {

          /* Compute the states in the integration points of the face. */
          LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], solDOFs, solInt);

          /*--- Get the required information from the standard element. ---*/
          const unsigned short ind          = surfElem[l].indStandardElement;
          const unsigned short nInt         = standardBoundaryFacesSol[ind].GetNIntegration();
          const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
          const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();
          const su2double     *weights      = standardBoundaryFacesSol[ind].GetWeightsIntegration();

          /*--- Store the solution of the DOFs of the adjacent element in contiguous
                memory such that the function DenseMatrixProduct can be used to compute
                the gradients solution variables in the integration points of the face. ---*/
          for(unsigned short i=0; i<nDOFsElem; ++i) {
            const su2double *solDOFElem = VecSolDOFs.data() + nVar*surfElem[l].DOFsSolElement[i];
            su2double       *sol        = solDOFs + nVar*i;
           for(unsigned short j=0; j<nVar; ++j)
             sol[j] = solDOFElem[j];
          }

          /* Compute the gradients in the integration points. Call the general function to
             carry out the matrix product. */
          config->GEMM_Tick(&tick);
          DenseMatrixProduct(nInt*nDim, nVar, nDOFsElem, derBasisElem, solDOFs, gradSolInt);
          config->GEMM_Tock(tick, "Friction_Forces", nInt*nDim, nVar, nDOFsElem);

          /* Determine the offset between r- and -s-derivatives, which is also the
             offset between s- and t-derivatives. */
          const unsigned short offDeriv = nVar*nInt;

          /* Loop over the integration points of this surface element. */
          for(unsigned short i=0; i<nInt; ++i) {

            /* Easier storage of the solution, its gradients, the normals, the
               metric terms and the coordinates for this integration point. */
            const su2double *sol         = solInt     + i*nVar;
            const su2double *gradSol     = gradSolInt + nVar*i;
            const su2double *normals     = surfElem[l].metricNormalsFace.data()
                                         + i*(nDim+1);
            const su2double *metricTerms = surfElem[l].metricCoorDerivFace.data()
                                         + i*nDim*nDim;
            const su2double *Coord       = surfElem[l].coorIntegrationPoints.data()
                                         + i*nDim;

            /*--- Compute the Cartesian gradients of the solution. ---*/
            su2double solGradCart[5][3];
            for(unsigned short k=0; k<nDim; ++k) {
              for(unsigned short j=0; j<nVar; ++j) {
                solGradCart[j][k] = 0.0;
                for(unsigned short l=0; l<nDim; ++l)
                  solGradCart[j][k] += gradSol[j+l*offDeriv]*metricTerms[k+l*nDim];
              }
            }

            /*--- Compute the velocities and static energy in this integration point. ---*/
            const su2double DensityInv = 1.0/sol[0];
            su2double vel[3], Velocity2 = 0.0;
            for(unsigned short j=0; j<nDim; ++j) {
              vel[j]     = sol[j+1]*DensityInv;
              Velocity2 += vel[j]*vel[j];
            }

            const su2double TotalEnergy  = sol[nDim+1]*DensityInv;
            const su2double StaticEnergy = TotalEnergy - 0.5*Velocity2;

            /*--- Compute the Cartesian gradients of the velocities and static energy
                  in this integration point and also the divergence of the velocity. ---*/
            su2double velGrad[3][3], StaticEnergyGrad[3], divVel = 0.0;
            for(unsigned short k=0; k<nDim; ++k) {
              StaticEnergyGrad[k] = DensityInv*(solGradCart[nDim+1][k]
                                  -             TotalEnergy*solGradCart[0][k]);
              for(unsigned short j=0; j<nDim; ++j) {
                velGrad[j][k]        = DensityInv*(solGradCart[j+1][k]
                                     -      vel[j]*solGradCart[0][k]);
                StaticEnergyGrad[k] -= vel[j]*velGrad[j][k];
              }
              divVel += velGrad[k][k];
            }

            /*--- Compute the laminar viscosity. ---*/
            FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
            const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

            /*--- Set the value of the second viscosity and compute the divergence
                  term in the viscous normal stresses. ---*/
            const su2double lambda     = -TWO3*ViscosityLam;
            const su2double lamDivTerm =  lambda*divVel;

            /*--- Compute the viscous stress tensor and the normal flux. Note that
                  there is a plus sign for the heat flux, because the normal
                  points into the geometry. ---*/
            su2double tauVis[3][3], qHeatNorm = 0.0;
            for(unsigned short k=0; k<nDim; ++k) {
              tauVis[k][k] = 2.0*ViscosityLam*velGrad[k][k] + lamDivTerm;    // Normal stress
              for(unsigned short j=(k+1); j<nDim; ++j) {
                tauVis[j][k] = ViscosityLam*(velGrad[j][k] + velGrad[k][j]); // Shear stress
                tauVis[k][j] = tauVis[j][k];
              }

              qHeatNorm += ViscosityLam*factHeatFlux_Lam*StaticEnergyGrad[k]*normals[k];
            }

            /*-- Compute the vector from the reference point to the integration
                 point and update the viscous force. Note that the normal points
                 into the geometry, hence the minus sign for the stress. ---*/
            su2double MomentDist[3], Force[3];
            const su2double scaleFac = weights[i]*normals[nDim]*factor;

            for(unsigned short iDim=0; iDim<nDim; ++iDim) {
              Force[iDim] = 0.0;
              for(unsigned short jDim=0; jDim<nDim; ++jDim)
                Force[iDim] -= tauVis[iDim][jDim]*normals[jDim];

              MomentDist[iDim]    = Coord[iDim] - Origin[iDim];
              Force[iDim]        *= scaleFac;
              ForceViscous[iDim] += Force[iDim];
            }

            /*--- Update the viscous moment. ---*/
            if (nDim == 3) {
              MomentViscous[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLengthMoment;
              MomentViscous[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
            }
            MomentViscous[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLengthMoment;

            /* Update the heat flux and maximum heat flux for this marker. */
            Heat_Visc[iMarker] += qHeatNorm*weights[i]*normals[nDim];
            MaxHeatFlux_Visc[iMarker] = max(MaxHeatFlux_Visc[iMarker], fabs(qHeatNorm));
          }
        }

        /*--- Project forces and store the non-dimensional coefficients ---*/
        if (nDim == 2) {
          CD_Visc[iMarker]   =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CL_Visc[iMarker]   = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker] = CL_Visc[iMarker] / (CD_Visc[iMarker]+EPS);
          CMz_Visc[iMarker]  = MomentViscous[2];
          CFx_Visc[iMarker]  = ForceViscous[0];
          CFy_Visc[iMarker]  = ForceViscous[1];
        }
        if (nDim == 3) {
          CD_Visc[iMarker]   =  ForceViscous[0]*cos(Alpha)*cos(Beta)
                             +  ForceViscous[1]*sin(Beta)
                             +  ForceViscous[2]*sin(Alpha)*cos(Beta);
          CL_Visc[iMarker]   = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSF_Visc[iMarker]  = -ForceViscous[0]*sin(Beta)*cos(Alpha)
                             +  ForceViscous[1]*cos(Beta)
                             -  ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker] = CL_Visc[iMarker]/(CD_Visc[iMarker] + EPS);
          CMx_Visc[iMarker]  = MomentViscous[0];
          CMy_Visc[iMarker]  = MomentViscous[1];
          CMz_Visc[iMarker]  = MomentViscous[2];
          CFx_Visc[iMarker]  = ForceViscous[0];
          CFy_Visc[iMarker]  = ForceViscous[1];
          CFz_Visc[iMarker]  = ForceViscous[2];
        }

        AllBound_CD_Visc   += CD_Visc[iMarker];
        AllBound_CL_Visc   += CL_Visc[iMarker];
        AllBound_CSF_Visc  += CSF_Visc[iMarker];
        AllBound_CMx_Visc  += CMx_Visc[iMarker];
        AllBound_CMy_Visc  += CMy_Visc[iMarker];
        AllBound_CMz_Visc  += CMz_Visc[iMarker];
        AllBound_CFx_Visc  += CFx_Visc[iMarker];
        AllBound_CFy_Visc  += CFy_Visc[iMarker];
        AllBound_CFz_Visc  += CFz_Visc[iMarker];

        AllBound_HeatFlux_Visc   += Heat_Visc[iMarker];
        AllBound_MaxHeatFlux_Visc = max(AllBound_MaxHeatFlux_Visc,
                                        MaxHeatFlux_Visc[iMarker]);

        /*--- Compute the coefficients per surface ---*/
        for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                         ++iMarker_Monitoring) {
          string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Visc[iMarker_Monitoring]   += CL_Visc[iMarker];
            Surface_CD_Visc[iMarker_Monitoring]   += CD_Visc[iMarker];
            Surface_CSF_Visc[iMarker_Monitoring]  += CSF_Visc[iMarker];
            Surface_CEff_Visc[iMarker_Monitoring] += CEff_Visc[iMarker];
            Surface_CFx_Visc[iMarker_Monitoring]  += CFx_Visc[iMarker];
            Surface_CFy_Visc[iMarker_Monitoring]  += CFy_Visc[iMarker];
            Surface_CFz_Visc[iMarker_Monitoring]  += CFz_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]  += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]  += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]  += CMz_Visc[iMarker];
          }
        }
      }
    }
  }

#ifdef HAVE_MPI

  /*--- Parallel mode. The data from all ranks must be gathered.
        Determine the size of the communication buffer. ---*/
  const unsigned long nCommSize = 9*config->GetnMarker_Monitoring() + 10;

  /*--- Define the communication buffers and store to local data in
        the local buffer. ---*/
  vector<su2double> locBuf(nCommSize), globBuf(nCommSize);

  unsigned long ii = 0;
  locBuf[ii++] = AllBound_CD_Visc;  locBuf[ii++] = AllBound_CL_Visc;
  locBuf[ii++] = AllBound_CSF_Visc; locBuf[ii++] = AllBound_CMx_Visc;
  locBuf[ii++] = AllBound_CMy_Visc; locBuf[ii++] = AllBound_CMz_Visc;
  locBuf[ii++] = AllBound_CFx_Visc; locBuf[ii++] = AllBound_CFy_Visc;
  locBuf[ii++] = AllBound_CFz_Visc; locBuf[ii++] = AllBound_HeatFlux_Visc;

  for(unsigned short i=0; i<config->GetnMarker_Monitoring(); ++i) {
    locBuf[ii++] = Surface_CL_Visc[i];  locBuf[ii++] = Surface_CD_Visc[i];
    locBuf[ii++] = Surface_CSF_Visc[i]; locBuf[ii++] = Surface_CFx_Visc[i];
    locBuf[ii++] = Surface_CFy_Visc[i]; locBuf[ii++] = Surface_CFz_Visc[i];
    locBuf[ii++] = Surface_CMx_Visc[i]; locBuf[ii++] = Surface_CMy_Visc[i];
    locBuf[ii++] = Surface_CMz_Visc[i];
  }

  /* Sum up all the data from all ranks. The result will be available on all ranks. */
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    SU2_MPI::Allreduce(locBuf.data(), globBuf.data(), nCommSize, MPI_DOUBLE,
                       MPI_SUM, MPI_COMM_WORLD);
  }

  /*--- Copy the data back from globBuf into the required variables. ---*/
  ii = 0;
  AllBound_CD_Visc  = globBuf[ii++]; AllBound_CL_Visc       = globBuf[ii++];
  AllBound_CSF_Visc = globBuf[ii++]; AllBound_CMx_Visc      = globBuf[ii++];
  AllBound_CMy_Visc = globBuf[ii++]; AllBound_CMz_Visc      = globBuf[ii++];
  AllBound_CFx_Visc = globBuf[ii++]; AllBound_CFy_Visc      = globBuf[ii++];
  AllBound_CFz_Visc = globBuf[ii++]; AllBound_HeatFlux_Visc = globBuf[ii++];

  AllBound_CEff_Visc = AllBound_CL_Visc/(AllBound_CD_Visc + EPS);

  for(unsigned short i=0; i<config->GetnMarker_Monitoring(); ++i) {
    Surface_CL_Visc[i]  = globBuf[ii++]; Surface_CD_Visc[i]  = globBuf[ii++];
    Surface_CSF_Visc[i] = globBuf[ii++]; Surface_CFx_Visc[i] = globBuf[ii++];
    Surface_CFy_Visc[i] = globBuf[ii++]; Surface_CFz_Visc[i] = globBuf[ii++];
    Surface_CMx_Visc[i] = globBuf[ii++]; Surface_CMy_Visc[i] = globBuf[ii++];
    Surface_CMz_Visc[i] = globBuf[ii++];

    Surface_CEff_Visc[i] = Surface_CL_Visc[i]/(Surface_CD_Visc[i] + EPS);
  }

  /* Determine the maximum heat flux over all ranks. */
  su2double localMax = AllBound_MaxHeatFlux_Visc;
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    SU2_MPI::Allreduce(&localMax, &AllBound_MaxHeatFlux_Visc, 1, MPI_DOUBLE,
                       MPI_MAX, MPI_COMM_WORLD);
  }
#endif

  /*--- Update the total coefficients (note that all the nodes have the same value)---*/
  Total_CD   += AllBound_CD_Visc;
  Total_CL   += AllBound_CL_Visc;
  Total_CSF  += AllBound_CSF_Visc;
  Total_CEff  = Total_CL / (Total_CD + EPS);
  Total_CMx  += AllBound_CMx_Visc;
  Total_CMy  += AllBound_CMy_Visc;
  Total_CMz  += AllBound_CMz_Visc;
  Total_CFx  += AllBound_CFx_Visc;
  Total_CFy  += AllBound_CFy_Visc;
  Total_CFz  += AllBound_CFz_Visc;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  for (unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                    ++iMarker_Monitoring) {
    Surface_CL[iMarker_Monitoring]   += Surface_CL_Visc[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]   += Surface_CD_Visc[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring]  += Surface_CSF_Visc[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]  = Surface_CL[iMarker_Monitoring] / (Surface_CD[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]  += Surface_CFx_Visc[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]  += Surface_CFy_Visc[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]  += Surface_CFz_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]  += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]  += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]  += Surface_CMz_Visc[iMarker_Monitoring];
  }
}

void CFEM_DG_NSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /* Check whether or not a time stepping scheme is used. */
  const bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  /* Constant factor present in the heat flux vector, namely the ratio of
     thermal conductivity and viscosity. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /* Constant ratio of the second viscosity and the viscosity itself. */
  const su2double lambdaOverMu = -TWO3;

  /* The eigenvalues of the viscous Jacobian, scaled by the kinematic viscosity,
     are 1.0, 2.0 + lambdaOverMu and kOverCv/Mu. The last is variable due to the
     possible presence of an eddy viscosity, but the first two are constant and
     the maximum can be determined. */
  const su2double radOverNuTerm = max(1.0, 2.0+lambdaOverMu);

  /* Store the number of metric points per DOF, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /* Initialize the minimum and maximum time step. */
  Min_Delta_Time = 1.e25; Max_Delta_Time = 0.0;

  /* Easier storage of the CFL number. Note that if we are using explicit
   time stepping, the regular CFL condition has been overwritten with the
   unsteady CFL condition in the config post-processing (if non-zero). */

  const su2double CFL = config->GetCFL(iMesh);

  /*--- Explicit time stepping with imposed time step (eventually will
   allow for local time stepping with this value imposed as the time
   for syncing the cells). If the unsteady CFL is set to zero (default),
   it uses the defined unsteady time step, otherwise it computes the time
   step based on the provided unsteady CFL. Note that the regular CFL
   option in the config is always ignored with time stepping. ---*/
  if (time_stepping && (config->GetUnst_CFL() == 0.0)) {

    /*--- Loop over the owned volume elements and set the fixed dt. ---*/
    for(unsigned long i=0; i<nVolElemOwned; ++i)
      VecDeltaTime[i] = config->GetDelta_UnstTimeND();

  } else {

    /*--- Check for a compressible solver. ---*/
    if(config->GetKind_Regime() == COMPRESSIBLE) {

      /*--- Loop over the owned volume elements. ---*/
      for(unsigned long i=0; i<nVolElemOwned; ++i) {

        /* Determine the offset between the r-derivatives and s-derivatives, which is
           also the offset between s- and t-derivatives, of the solution in the DOFs. */
        const unsigned short offDerivSol = nVar*volElem[i].nDOFsSol;

        /*--- Compute the length scale for this element. Note that in the
              length scale the polynomial degree must be taken into account
              for the high order element. ---*/
        const unsigned short ind = volElem[i].indStandardElement;
        unsigned short nPoly = standardElementsSol[ind].GetNPoly();
        if(nPoly == 0) nPoly = 1;

        const su2double lenScaleInv = nPoly/volElem[i].lenScale;
        const su2double lenScale    = 1.0/lenScaleInv;

        /*--- Compute the gradients of the conserved variables if a subgrid
              scale model for LES is used. ---*/
        su2double *gradSolDOFs = VecTmpMemory.data();
        if( SGSModelUsed ) {
          const unsigned short nDOFs          = volElem[i].nDOFsSol;
          const su2double *matDerBasisSolDOFs = standardElementsSol[ind].GetMatDerBasisFunctionsSolDOFs();
          const su2double *sol                = VecSolDOFs.data() + nVar*volElem[i].offsetDOFsSolLocal;

          DenseMatrixProduct(nDOFs*nDim, nVar, nDOFs, matDerBasisSolDOFs, sol, gradSolDOFs);
        }

        /*--- Loop over the DOFs of this element and determine the maximum wave speed
              and the maximum value of the viscous spectral radius. ---*/
        su2double charVel2Max = 0.0, radViscMax = 0.0;
        for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
          const su2double *solDOF = VecSolDOFs.data() + nVar*(volElem[i].offsetDOFsSolLocal + j);

          /* Compute the velocities. */
          su2double velAbs[3], vel[3];
          const su2double DensityInv = 1.0/solDOF[0];
          su2double Velocity2 = 0.0;
          for(unsigned short iDim=0; iDim<nDim; ++iDim) {
            vel[iDim]    = solDOF[iDim+1]*DensityInv;
            velAbs[iDim] = fabs(vel[iDim]);
            Velocity2 += vel[iDim]*vel[iDim];
          }

          /*--- Compute the maximum value of the wave speed. This is a rather
                conservative estimate. ---*/
          const su2double StaticEnergy = solDOF[nDim+1]*DensityInv - 0.5*Velocity2;
          FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
          const su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
          const su2double SoundSpeed  = sqrt(fabs(SoundSpeed2));

          su2double charVel2 = 0.0;
          for(unsigned short iDim=0; iDim<nDim; ++iDim) {
            const su2double rad = velAbs[iDim] + SoundSpeed;
            charVel2 += rad*rad;
          }

          charVel2Max = max(charVel2Max, charVel2);

          /* Compute the laminar kinematic viscosity and check if an eddy
             viscosity must be determined. */
          const su2double muLam = FluidModel->GetLaminarViscosity();
          su2double muTurb      = 0.0;

          if( SGSModelUsed ) {

            /* Set the pointer gradSolDOF to the location where the gradients
               of this DOF start. */
            const su2double *gradSolDOF = gradSolDOFs + j*nVar;

            /* Easier storage of the metric terms in this DOF. First compute the
               inverse of the Jacobian, the Jacobian is the first entry in the metric
               terms, and afterwards update the metric terms by 1. */
            const su2double *metricTerms = volElem[i].metricTermsSolDOFs.data()
                                         + j*nMetricPerPoint;
            const su2double Jac          = metricTerms[0];
            const su2double JacInv       = 1.0/Jac;
            metricTerms                 += 1;

            /*--- Compute the Cartesian gradients of the independent solution
                  variables from the gradients in parametric coordinates and the metric
                  terms in this DOF. Note that at the end a multiplication with JacInv
                  takes places, because the metric terms are scaled by the Jacobian. ---*/
            su2double solGradCart[5][3];
            for(unsigned short k=0; k<nDim; ++k) {
              for(unsigned short l=0; l<nVar; ++l) {
                solGradCart[l][k] = 0.0;
                for(unsigned short m=0; m<nDim; ++m)
                  solGradCart[l][k] += gradSolDOF[l+m*offDerivSol]*metricTerms[k+m*nDim];
                solGradCart[l][k] *= JacInv;
              }
            }

            /* Compute the Cartesian gradients of the velocities in this DOF. */
            su2double velGrad[3][3];
            for(unsigned short k=0; k<nDim; ++k) {
              for(unsigned short l=0; l<nDim; ++l) {
                velGrad[l][k] = DensityInv*(solGradCart[l+1][k]
                              -      vel[l]*solGradCart[0][k]);
              }
            }

            /* Compute the eddy viscosity. */
            const su2double dist = volElem[i].wallDistanceSolDOFs[j];
            muTurb = SGSModel->ComputeEddyViscosity(nDim, solDOF[0], velGrad,
                                                    lenScale, dist);
          }

          /*--- Determine the viscous spectral radius. ---*/
          const su2double mu           = muLam + muTurb;
          const su2double kOverCv      = muLam*factHeatFlux_Lam
                                       + muTurb*factHeatFlux_Turb;
          const su2double factHeatFlux = kOverCv/mu;

          const su2double radVisc = DensityInv*mu*max(radOverNuTerm, factHeatFlux);

          /* Update the maximum value of the viscous spectral radius. */
          radViscMax = max(radViscMax, radVisc);
        }

        /*--- Compute the time step for the element and update the minimum and
              maximum value. Take the factor for time accurate local time
              stepping into account for the minimum and maximum. ---*/
        const su2double dtInv = lenScaleInv*(sqrt(charVel2Max) + radViscMax*lenScaleInv);

        VecDeltaTime[i] = CFL/dtInv;

        Min_Delta_Time = min(Min_Delta_Time, volElem[i].factTimeLevel*VecDeltaTime[i]);
        Max_Delta_Time = max(Max_Delta_Time, volElem[i].factTimeLevel*VecDeltaTime[i]);
      }
    }
    else {

      /*--- Incompressible solver. ---*/

      int rank = MASTER_NODE;
#ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      if(rank == MASTER_NODE) {
        cout << "In function CFEM_DG_EulerSolver::SetTime_Step" << endl;
        cout << "Incompressible solver not implemented yet" << endl;
      }

#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#else
      exit(EXIT_FAILURE);
#endif

    }

    /*--- Compute the max and the min dt (in parallel). Note that we only
     so this for steady calculations if the high verbosity is set, but we
     always perform the reduction for unsteady calculations where the CFL
     limit is used to set the global time step. ---*/
    if ((config->GetConsole_Output_Verb() == VERB_HIGH) || time_stepping) {
#ifdef HAVE_MPI
      su2double rbuf_time = Min_Delta_Time;
      SU2_MPI::Allreduce(&rbuf_time, &Min_Delta_Time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      rbuf_time = Max_Delta_Time;
      SU2_MPI::Allreduce(&rbuf_time, &Max_Delta_Time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    }

    /*--- For explicit time stepping with an unsteady CFL imposed, use the
          minimum delta time of the entire mesh. As Min_Delta_Time is scaled to
          the time step of the largest time level, a correction must be used
          for the time level when time accurate local time stepping is used. ---*/
    if (time_stepping) {
      for(unsigned long i=0; i<nVolElemOwned; ++i)
        VecDeltaTime[i] = Min_Delta_Time/volElem[i].factTimeLevel;
    }
  }
}

void CFEM_DG_NSSolver::ADER_DG_AliasedPredictorResidual(CConfig           *config,
                                                        CVolumeElementFEM *elem,
                                                        const su2double   *sol,
                                                        su2double         *res,
                                                        su2double         *work) {

  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Get the necessary information from the standard element. ---*/
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *matDerBasisInt         = matBasisInt + nDOFs*nInt;
  const su2double *matDerBasisSolDOFs     = standardElementsSol[ind].GetMatDerBasisFunctionsSolDOFs();
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  unsigned short nPoly = standardElementsSol[ind].GetNPoly();
  if(nPoly == 0) nPoly = 1;

  /* Compute the length scale of the current element for the LES. */
  const su2double lenScale = elem->lenScale/nPoly;

  /* Set the pointers for fluxesDOF, gradFluxesInt, gradSolDOF and divFlux.
     Note that the same array can be used for the fluxesDOF and divFlux and
     gradFluxesInt and gradSolDOFs, because this information is needed after
     each other. */
  su2double *fluxesDOF     = work;
  su2double *gradFluxesInt = fluxesDOF + nDOFs*nVar*nDim;
  su2double *gradSolDOFs   = gradFluxesInt;
  su2double *divFlux       = work;

  /* Determine the offset between the r-derivatives and s-derivatives, which is
     also the offset between s- and t-derivatives, of the fluxes in the
     integration points and the solution in the DOFs. */
  const unsigned short offDerivSol    = nVar*nDOFs;
  const unsigned short offDerivFluxes = nVar*nInt*nDim;

  /* Store the number of metric points per integration point/DOF, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Determine the Cartesian fluxes of the Navier-Stokes        ---*/
  /*---         equations in the DOFs.                                     ---*/
  /*--------------------------------------------------------------------------*/

  /* Compute the derivatives of the solution variables w.r.t. the parametric
     coordinates in the DOFs. */
  DenseMatrixProduct(nDOFs*nDim, nVar, nDOFs, matDerBasisSolDOFs, sol, gradSolDOFs);

  /* Loop over the DOFs. */
  for(unsigned short i=0; i<nDOFs; ++i) {

    /* Set the pointers to the location where the solution, the gradients of the
       solution and the location where the Cartesian fluxes of this DOF
       are stored. */
    const su2double *solDOF     = sol         + i*nVar;
    const su2double *gradSolDOF = gradSolDOFs + i*nVar;
    su2double       *fluxes     = fluxesDOF   + i*nVar*nDim;

    /* Easier storage of the metric terms in this DOF. First compute the
       inverse of the Jacobian, the Jacobian is the first entry in the metric
       terms, and afterwards update the metric terms by 1. */
    const su2double *metricTerms = elem->metricTermsSolDOFs.data()
                                 + i*nMetricPerPoint;
    const su2double Jac          = metricTerms[0];
    const su2double JacInv       = 1.0/Jac;
    metricTerms                 += 1;

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the metric
          terms in this DOF. Note that at the end a multiplication with JacInv
          takes places, because the metric terms are scaled by the Jacobian. ---*/
    su2double solGradCart[5][3];
    for(unsigned short k=0; k<nDim; ++k) {
      for(unsigned short j=0; j<nVar; ++j) {
        solGradCart[j][k] = 0.0;
        for(unsigned short l=0; l<nDim; ++l)
          solGradCart[j][k] += gradSolDOF[j+l*offDerivSol]*metricTerms[k+l*nDim];
        solGradCart[j][k] *= JacInv;
      }
    }

    /*--- Compute the velocities and static energy in this DOF. ---*/
    const su2double DensityInv = 1.0/solDOF[0];
    su2double vel[3], Velocity2 = 0.0;
    for(unsigned short j=0; j<nDim; ++j) {
      vel[j]     = solDOF[j+1]*DensityInv;
      Velocity2 += vel[j]*vel[j];
    }

    const su2double TotalEnergy  = solDOF[nDim+1]*DensityInv;
    const su2double StaticEnergy = TotalEnergy - 0.5*Velocity2;

    /*--- Compute the Cartesian gradients of the velocities and static energy
          in this integration point and also the divergence of the velocity. ---*/
    su2double velGrad[3][3], StaticEnergyGrad[3], divVel = 0.0;
    for(unsigned short k=0; k<nDim; ++k) {
      StaticEnergyGrad[k] = DensityInv*(solGradCart[nDim+1][k]
                          -             TotalEnergy*solGradCart[0][k]);
      for(unsigned short j=0; j<nDim; ++j) {
        velGrad[j][k]        = DensityInv*(solGradCart[j+1][k]
                             -      vel[j]*solGradCart[0][k]);
        StaticEnergyGrad[k] -= vel[j]*velGrad[j][k];
      }
      divVel += velGrad[k][k];
    }

    /*--- Compute the pressure and the laminar viscosity. ---*/
    FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
    const su2double Pressure     = FluidModel->GetPressure();
    const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

    /*--- Compute the eddy viscosity, if needed, and the total viscosity. ---*/
    su2double ViscosityTurb = 0.0;
    if( SGSModelUsed ) {
      const su2double dist = elem->wallDistanceSolDOFs[i];
      ViscosityTurb = SGSModel->ComputeEddyViscosity(nDim, solDOF[0], velGrad,
                                                    lenScale, dist);
    }

    const su2double Viscosity = ViscosityLam + ViscosityTurb;

    /* Compute the total thermal conductivity divided by Cv. */
    const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                            + ViscosityTurb*factHeatFlux_Turb;

    /*--- Set the value of the second viscosity and compute the divergence
          term in the viscous normal stresses. ---*/
    const su2double lambda     = -TWO3*Viscosity;
    const su2double lamDivTerm =  lambda*divVel;

    /*--- Compute the viscous stress tensor. ---*/
    su2double tauVis[3][3];
    for(unsigned short k=0; k<nDim; ++k) {
      tauVis[k][k] = 2.0*Viscosity*velGrad[k][k] + lamDivTerm;    // Normal stress
      for(unsigned short j=(k+1); j<nDim; ++j) {
        tauVis[j][k] = Viscosity*(velGrad[j][k] + velGrad[k][j]); // Shear stress
        tauVis[k][j] = tauVis[j][k];
      }
    }

    /* Loop over the number of dimensions for the number of fluxes. */
    const su2double rH = solDOF[nDim+1] + Pressure;
    unsigned short ll = 0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim) {

      /* Mass flux for the current direction. */
      fluxes[ll++] = solDOF[iDim+1];

      /* Momentum fluxes for the current direction. */
      for(unsigned short jDim=0; jDim<nDim; ++jDim)
        fluxes[ll+jDim] = solDOF[iDim+1]*vel[jDim] - tauVis[jDim][iDim];
      fluxes[ll+iDim]  += Pressure;

      ll += nDim;

      /* Energy flux for the current direction. */
      fluxes[ll] = rH*vel[iDim]                         // Inviscid part
                 - kOverCv*StaticEnergyGrad[iDim];      // Heat flux part
      for(unsigned short jDim=0; jDim<nDim; ++jDim)
        fluxes[ll] -= tauVis[jDim][iDim]*vel[jDim];     // Work of the viscous forces part
      ++ll;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Compute the divergence of the fluxes in the integration    ---*/
  /*---         points and distribute the divergence terms to the DOFs.    ---*/
  /*--------------------------------------------------------------------------*/

  /* Compute the derivatives of the Cartesian fluxes w.r.t. the parametric
     coordinates in the integration points. */
  DenseMatrixProduct(nInt*nDim, nVar*nDim, nDOFs, matDerBasisInt, fluxesDOF,
                     gradFluxesInt);

  /*--- Loop over the integration points to compute the divergence of the fluxes
        in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the data of the derivatives of the
       fluxes of this integration point starts. */
    const su2double *gradFluxes = gradFluxesInt + nVar*nDim*i;

    /* Easier storage of the metric terms in this integration point. The +1
       is present, because the first element of the metric terms is the
       Jacobian in the integration point. Also set the point where the divergence
       of the flux terms is stored for the current integration point. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint + 1;
    su2double       *divFluxInt  = divFlux + nVar*i;

    /* Initialize the divergence to zero. */
    for(unsigned short j=0; j<nVar; ++j) divFluxInt[j] = 0.0;

    /* Loop over the nDim parametric coordinates. */
    unsigned short ll = 0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim) {

      /* Set the pointer to derivatives of this parametric coordinate. */
      const su2double *derFluxes = gradFluxes + iDim*offDerivFluxes;

      /* Loop over the number of Cartesian dimensions. */
      for(unsigned short jDim=0; jDim<nDim; ++jDim, ++ll) {

        /* Set the pointer to the derivatives of the Cartesian flux in
           the jDim direction. */
        const su2double *derFlux = derFluxes + jDim*nVar;

        /* Update the divergence for all equations. */
        for(unsigned short j=0; j<nVar; ++j)
          divFluxInt[j] += derFlux[j]*metricTerms[ll];
      }
    }

    /* Multiply the divergence with the integration weight. */
    for(unsigned short j=0; j<nVar; ++j) divFluxInt[j] *= weights[i];
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_NSSolver::ADER_DG_NonAliasedPredictorResidual(CConfig           *config,
                                                           CVolumeElementFEM *elem,
                                                           const su2double   *sol,
                                                           su2double         *res,
                                                           su2double         *work) {

  /* Constant factor present in the heat flux vector, the inverse of
     the specific heat at constant volume and ratio lambdaOverMu. */
  const su2double factHeatFlux_Lam  =  Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb =  Gamma/Prandtl_Turb;
  const su2double Gas_Constant      =  config->GetGas_ConstantND();
  const su2double CvInv             =  Gamma_Minus_One/Gas_Constant;
  const su2double lambdaOverMu      = -TWO3;

  /*--- Get the necessary information from the standard element. ---*/
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *matDerBasisInt         = matBasisInt + nDOFs*nInt;
  const su2double *matDerBasisSolDOFs     = standardElementsSol[ind].GetMatDerBasisFunctionsSolDOFs();
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  unsigned short nPoly = standardElementsSol[ind].GetNPoly();
  if(nPoly == 0) nPoly = 1;

  /* Compute the length scale of the current element for the LES. */
  const su2double lenScale = elem->lenScale/nPoly;

  /* Set the pointers for solAndGradInt and divFlux to work. The same array
     can be used for both help arrays. Set the pointer for the second
     derivatives such that they are stored after the first derivatives. */
  su2double *solAndGradInt      = work;
  su2double *divFlux            = work;
  const unsigned short offPoint = nVar*(nDim+1)*max(nInt, nDOFs);
  su2double *secDerSolInt       = solAndGradInt + offPoint;

  /* Determine the offset between the solution variables and the r-derivatives,
     which is also the offset between the r- and s-derivatives and the offset
     between s- and t-derivatives. Also determine the offset between the data
     of the second derivatives. */
  const unsigned short offDeriv    = nVar*nInt;
  const unsigned short off2ndDeriv = nDim*offDeriv;

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /* Store the number of additional metric points per integration point, which
     are needed to compute the second derivatives. These terms take the
     non-constant metric into account. */
  const unsigned short nMetric2ndDerPerPoint = nDim*(nDim + nDim*(nDim-1)/2);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Interpolate the conserved solution variables to the        ---*/
  /*---         integration points and also determine the first and second ---*/
  /*---         derivatives of these variables in the integration points.  ---*/
  /*---         All derivatives are w.r.t. the parametric coordinates.     ---*/
  /*--------------------------------------------------------------------------*/

  /*--- The computation of the second derivatives happens in two stages. First
        the first derivatives in the DOFs are determined and these are then
        differentiated again to obtain the second derivatives in the integration
        points. Note that secDerSolInt and solAndGradInt are used for temporary
        storage of the first derivatives in the DOFs. ---*/
  su2double *gradSolDOFs = secDerSolInt;
  DenseMatrixProduct(nDOFs*nDim, nVar, nDOFs, matDerBasisSolDOFs, sol, gradSolDOFs);

  /* Store the gradients in true row major order for each DOF. */
  su2double *firstDerSolDOFs = solAndGradInt;
  unsigned short ll = 0;
  for(unsigned short i=0; i<nDOFs; ++i) {

    for(unsigned short j=0; j<nDim; ++j) {
      const su2double *gradSolDOF = gradSolDOFs + nVar*(i + j*nDOFs);

      for(unsigned short k=0; k<nVar; ++k, ++ll)
        firstDerSolDOFs[ll] = gradSolDOF[k];
    }
  }

  /* Compute the second derivatives w.r.t. the parametric coordinates
     in the integration points. */
  DenseMatrixProduct(nInt*nDim, nVar*nDim, nDOFs, matDerBasisInt,
                     firstDerSolDOFs, secDerSolInt);

  /* Compute the solution and the derivatives w.r.t. the parametric coordinates
     in the integration points. */
  DenseMatrixProduct(nInt*(nDim+1), nVar, nDOFs, matBasisInt, sol, solAndGradInt);

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Determine from the solution and gradients the divergence   ---*/
  /*---         of the fluxes in the integration points.                   ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the integration points to compute the divergence of the fluxes
        in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the solution data of this
       integration point starts. */
    const su2double *sol = solAndGradInt + nVar*i;

    /*--- Compute the velocities and static energy in this integration point. ---*/
    const su2double DensityInv = 1.0/sol[0];
    su2double vel[3], Velocity2 = 0.0;
    for(unsigned short j=0; j<nDim; ++j) {
      vel[j]     = sol[j+1]*DensityInv;
      Velocity2 += vel[j]*vel[j];
    }

    const su2double TotalEnergy  = sol[nDim+1]*DensityInv;
    const su2double StaticEnergy = TotalEnergy - 0.5*Velocity2;

    /*--- Compute the laminar viscosity, its derivative w.r.t. temperature,
          pressure and the total enthalpy. ---*/
    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
    const su2double ViscosityLam = FluidModel->GetLaminarViscosity();
    const su2double dViscLamdT   = FluidModel->GetdmudT_rho();
    const su2double Pressure     = FluidModel->GetPressure();
    const su2double Htot         = (sol[nDim+1]+Pressure)*DensityInv;

    /* Easier storage of the metric terms in this integration point.  Store the
       Jacobian and its inverse for later purposes and update the pointer for
       the metric terms by one, such that it starts at the position where drdx
       of this integration point is stored. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
    const su2double Jac          = metricTerms[0];
    const su2double JacInv       = 1.0/Jac;
    const su2double JacInv2      = JacInv*JacInv;
    metricTerms                 += 1;

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the
          metric terms in this integration point. Note that at the end a
          multiplication with JacInv takes places, because the metric terms
          are scaled by the Jacobian. ---*/
    su2double solGradCart[5][3];
    for(unsigned short k=0; k<nDim; ++k) {
      for(unsigned short j=0; j<nVar; ++j) {
        solGradCart[j][k] = 0.0;
        for(unsigned short l=0; l<nDim; ++l)
          solGradCart[j][k] += sol[j+(l+1)*offDeriv]*metricTerms[k+l*nDim];
        solGradCart[j][k] *= JacInv;
      }
    }

    /* Easier storage of the location where the data of second derivatives of
       this integration point starts as well as the necessary additional metric
       terms to compute the Cartesian second derivatives. */
    const su2double *sol2ndDer         = secDerSolInt + nVar*nDim*i;
    const su2double *metricTerms2ndDer = elem->metricTerms2ndDer.data()
                                       + i*nMetric2ndDerPerPoint;

    /*--- Compute the Cartesian second derivatives of the independent solution
          variables from the gradients and second derivatives in parametric
          coordinates and the metric terms and its derivatives w.r.t. the
          parametric coordinates. ---*/
    su2double sol2ndDerCart[5][3][3];
    for(unsigned short j=0; j<nVar; ++j) {
      unsigned short llDerMet = 0;
      for(unsigned short l=0; l<nDim; ++l) {
        for(unsigned short k=0; k<=l; ++k) {
          sol2ndDerCart[j][k][l] = 0.0;

          /* First the terms coming from the second derivatives w.r.t.
             the parametric coordinates. The result is multiplied by the
             square of the inverse of the Jacobian, because this term is
             not present in the metric terms itself. */
          for(unsigned short m=0; m<nDim; ++m) {
            const su2double *secDerJ = sol2ndDer + m*off2ndDeriv + j;
            for(unsigned short n=0; n<nDim; ++n)
              sol2ndDerCart[j][k][l] += secDerJ[n*nVar]*metricTerms[l+m*nDim]
                                      * metricTerms[k+n*nDim];
          }
          sol2ndDerCart[j][k][l] *= JacInv2;

          /* There is also a contribution to the second derivative from the first
             derivatives and non-constant metric terms. This is added here. */
          for(unsigned short m=0; m<nDim; ++m, ++llDerMet)
            sol2ndDerCart[j][k][l] += metricTerms2ndDer[llDerMet]*sol[j+(m+1)*offDeriv];

          /* The Hessian is symmetric, so copy this derivative to its
             symmetric counterpart. */
          sol2ndDerCart[j][l][k] = sol2ndDerCart[j][k][l];
        }
      }
    }

    /*--- Compute the Cartesian gradients of the pressure, velocity components,
          static energy and dynamic viscosity. ---*/
    su2double pGrad[3], velGrad[3][3], StaticEnergyGrad[3], ViscosityLamGrad[3];
    for(unsigned short k=0; k<nDim; ++k) {
      pGrad[k]            = solGradCart[nVar-1][k] + 0.5*Velocity2*solGradCart[0][k];
      StaticEnergyGrad[k] = DensityInv*(solGradCart[nDim+1][k]
                          -             TotalEnergy*solGradCart[0][k]);
      for(unsigned short l=0; l<nDim; ++l) {
        pGrad[k]            -= vel[l]*solGradCart[l+1][k];
        velGrad[l][k]        = DensityInv*(solGradCart[l+1][k]
                             -             vel[l]*solGradCart[0][k]);
        StaticEnergyGrad[k] -= vel[l]*velGrad[l][k];
      }
      pGrad[k] *= Gamma_Minus_One;

      ViscosityLamGrad[k] = CvInv*StaticEnergyGrad[k]*dViscLamdT;
    }

    /*--- Compute the second derivatives of the velocity components. ---*/
    su2double vel2ndDer[3][3][3];
    for(unsigned short j=0; j<nDim; ++j) {
      for(unsigned short l=0; l<nDim; ++l) {
        for(unsigned short k=0; k<=l; ++k) {
          vel2ndDer[j][k][l] = DensityInv*(sol2ndDerCart[j+1][k][l] - vel[j]*sol2ndDerCart[0][k][l]
                             +     DensityInv*(2.0*vel[j]*solGradCart[0][k]*solGradCart[0][l]
                             -                 solGradCart[j+1][k]*solGradCart[0][l]
                             -                 solGradCart[j+1][l]*solGradCart[0][k]));
          vel2ndDer[j][l][k] = vel2ndDer[j][k][l];
        }
      }
    }

    /*--- Compute the second derivatives of the static energy. Note that this
          term appears in the heat flux and therefore only the pure second
          derivatives are needed. Hence, the cross-derivatives are omitted. ---*/
    su2double StaticEnergy2ndDer[3];
    for(unsigned short l=0; l<nDim; ++l) {
      StaticEnergy2ndDer[l] = DensityInv*(sol2ndDerCart[nDim+1][l][l] - TotalEnergy*sol2ndDerCart[0][0][l]
                            +     2.0*DensityInv*(TotalEnergy*solGradCart[0][l]*solGradCart[0][l]
                            -                     solGradCart[nDim+1][l]*solGradCart[0][l]));
      for(unsigned short k=0; k<nDim; ++k)
        StaticEnergy2ndDer[l] -= vel[k]*vel2ndDer[k][l][l] - velGrad[k][l]*velGrad[k][l];
    }

    /*--- If an SGS model is used the eddy viscosity and its spatial
          derivatives must be computed. ---*/
    su2double ViscosityTurb = 0.0;
    su2double ViscosityTurbGrad[] = {0.0, 0.0, 0.0};

    if( SGSModelUsed ) {
      const su2double dist = elem->wallDistance[i];
      ViscosityTurb = SGSModel->ComputeEddyViscosity(nDim, sol[0], velGrad,
                                                    lenScale, dist);

      su2double densityGrad[3];
      for(unsigned short l=0; l<nDim; ++l) densityGrad[l] = solGradCart[0][l];

      SGSModel->ComputeGradEddyViscosity(nDim, sol[0], densityGrad, velGrad,
                                         vel2ndDer, lenScale, dist,
                                         ViscosityTurbGrad);
    }

    /*--- Compute the total viscosity, the total heat conductivity and their
          gradients. Note that the heat conductivity is divided by the Cv,
          because gradients of internal energy are computed and not temperature. ---*/
    const su2double Viscosity = ViscosityLam + ViscosityTurb;
    const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                            + ViscosityTurb*factHeatFlux_Turb;

    su2double ViscosityGrad[3], kOverCvGrad[3];
    for(unsigned short k=0; k<nDim; ++k) {
      ViscosityGrad[k] = ViscosityLamGrad[k] + ViscosityTurbGrad[k];
      kOverCvGrad[k]   = ViscosityLamGrad[k] *factHeatFlux_Lam
                       + ViscosityTurbGrad[k]*factHeatFlux_Turb;
    }

    /* Abbreviations, which make it easier to compute the divergence term. */
    su2double abv1 = 0.0, abv2 = 0.0, abv3 = 0.0, abv4 = 0.0;
    for(unsigned short k=0; k<nDim; ++k) {
      abv1 += solGradCart[k+1][k];
      abv2 += vel[k]*solGradCart[0][k];
      abv3 += vel[k]*(solGradCart[nVar-1][k] + pGrad[k]);
      abv4 += velGrad[k][k];
    }

    /*--- Set the pointer to store the divergence terms for this integration
          point and compute these terms, multiplied by the integration weight
          and Jacobian. ---*/
    const su2double weightJac = weights[i]*Jac;
    su2double *divFluxInt     = divFlux + nVar*i;

    divFluxInt[0]      = weightJac*abv1;
    divFluxInt[nVar-1] = abv3 + Htot*(abv1 - abv2);

    for(unsigned short k=0; k<nDim; ++k) {
      divFluxInt[k+1]     = pGrad[k] + vel[k]*(abv1 - abv2)
                          - lambdaOverMu*ViscosityGrad[k]*abv4;
      divFluxInt[nVar-1] -= abv4*lambdaOverMu*(Viscosity*velGrad[k][k]
                          +                    vel[k]*ViscosityGrad[k])
                          + kOverCvGrad[k]*StaticEnergyGrad[k]
                          + kOverCv*StaticEnergy2ndDer[k];

      for(unsigned short l=0; l<nDim; ++l) {
        divFluxInt[k+1]    += vel[l]*solGradCart[k+1][l]
                            - lambdaOverMu*Viscosity*vel2ndDer[l][k][l]
                            - Viscosity*(vel2ndDer[k][l][l] + vel2ndDer[l][k][l])
                            - ViscosityGrad[l]*(velGrad[k][l] + velGrad[l][k]);
        divFluxInt[nVar-1] -= (Viscosity*velGrad[k][l] + vel[k]*ViscosityGrad[l])
                            * (velGrad[k][l] + velGrad[l][k])
                            + vel[k]*Viscosity*(vel2ndDer[k][l][l] + vel2ndDer[l][k][l]
                            +                   lambdaOverMu*vel2ndDer[l][k][l]);
      }

      divFluxInt[k+1] *= weightJac;
    }

    divFluxInt[nVar-1] *= weightJac;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Compute the residual in the DOFs, which is the matrix      ---*/
  /*---         product of basisFunctionsIntTrans and divFlux.             ---*/
  /*--------------------------------------------------------------------------*/

  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_NSSolver::Shock_Capturing_DG(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short iMesh, unsigned short iStep) {

  /*--- Run shock capturing algorithm ---*/
  switch( config->GetKind_FEM_DG_Shock() ) {
    case NONE:
      break;
    case PERSSON:
      Shock_Capturing_DG_Persson(geometry, solver_container, numerics, config, iMesh, iStep);
      break;
  }

}
void CFEM_DG_NSSolver::Shock_Capturing_DG_Persson(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short iMesh, unsigned short iStep) {

  /*--- Dummy variable for storing shock sensor value temporarily ---*/
  su2double sensorVal, sensorLowerBound, machNorm, machMax;
  su2double DensityInv, Velocity2, StaticEnergy, SoundSpeed2;

  bool shockExist;
  unsigned short nDOFsPm1;       // Number of DOFs up to polynomial degree p-1

  /*--- Loop over the owned volume elements to sense the shock. If shock exists,
        add artificial viscosity for DG FEM formulation to the residual.  ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Get the data from the corresponding standard element. */
    const unsigned short ind          = volElem[l].indStandardElement;
    const unsigned short nDOFs        = volElem[l].nDOFsSol;
    const unsigned short VTK_TypeElem = volElem[l].VTK_Type;
    const unsigned short nPoly        = standardElementsSol[ind].GetNPoly();
    const su2double *matVanderInv     = standardElementsSol[ind].GetMatVandermondeInv();

    /*----------------------------------------------------------------------------*/
    /*--- Step 1: Calculate the number of DOFs up to polynomial degree p-1.    ---*/
    /*----------------------------------------------------------------------------*/

    switch( VTK_TypeElem ) {
      case TRIANGLE:
        nDOFsPm1 = nPoly*(nPoly+1)/2;
        break;
      case QUADRILATERAL:
        nDOFsPm1 = nPoly*nPoly;
        break;
      case TETRAHEDRON:
        nDOFsPm1 = nPoly*(nPoly+1)*(nPoly+2)/6;
        break;
      case PYRAMID:
        nDOFsPm1 = nPoly*(nPoly+1)*(2*nPoly+1)/6;
        break;
      case PRISM:
        nDOFsPm1 = nPoly*nPoly*(nPoly+1)/2;
        break;
      case HEXAHEDRON:
        nDOFsPm1 = nPoly*nPoly*nPoly;
        break;
    }

    /*---------------------------------------------------------------------*/
    /*--- Step 2: Calculate the shock sensor value for this element.    ---*/
    /*---------------------------------------------------------------------*/

    /* Initialize dummy variable for this volume element */
    sensorVal = 0;
    machMax = -1;
    shockExist = false;
    sensorLowerBound = 1.e15;

    /* Easier storage of the solution variables for this element. */
    const su2double *solDOFs = VecSolDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Temporary storage of mach number for DOFs in this element. */
    su2double *machSolDOFs = VecTmpMemory.data();
    su2double *vecTemp     = machSolDOFs + nDOFs;

    /* Calculate primitive variables and mach number for DOFs in this element.
       Also, track the maximum mach number in this element. */
    for(unsigned short iInd=0; iInd<nDOFs; ++iInd) {

      const su2double *sol = solDOFs + iInd*nVar;
      DensityInv = 1.0/sol[0];
      Velocity2 = 0.0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double vel = sol[iDim]*DensityInv;
        Velocity2 += vel*vel;
      }

      StaticEnergy = sol[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
      SoundSpeed2 = FluidModel->GetSoundSpeed2();
      machSolDOFs[iInd] = sqrt( Velocity2/SoundSpeed2 );
      machMax = max(machSolDOFs[iInd],machMax);
    }

    /* Change the solution coefficients to modal form from nodal form */
    for(unsigned short i=0; i<nDOFs; ++i) {
      vecTemp[i] = 0.0;
      for (unsigned short j=0; j<nDOFs; ++j)
        vecTemp[i] += matVanderInv[i+j*nDOFs]*machSolDOFs[j];
    }

    /* Get the L2 norm of solution coefficients for the highest polynomial order. */
    for(unsigned short i=nDOFsPm1; i<nDOFs; ++i) {
        sensorVal += vecTemp[i]*vecTemp[i];
    }

    /* If the maximum mach number is greater than 1.0, try to calculate the shockSensorValue.
       Otherwise, assign default value. */
    if ( machMax > 1.0) {
      // !!!!!Threshold value for sensorVal should be further investigated
      if(sensorVal > 1.e-15) {
        machNorm = 0.0;

        /*--- Get L2 norm square of vecTemp ---*/
        for (unsigned short i=0; i<nDOFs; ++i) {
          machNorm += vecTemp[i]*vecTemp[i];
        }
        if (machNorm < 1.e-15) {
          // This should not happen
          volElem[l].shockSensorValue = 1000.0;
        }
        else {
          volElem[l].shockSensorValue = log(sensorVal/machNorm);
          shockExist = true;
        }
      }
      else {
        // There is no shock in this element
        volElem[l].shockSensorValue = -1000.0;
      }
    }
    else {
      volElem[l].shockSensorValue = -1000.0;
    }

    /*---------------------------------------------------------------------*/
    /*--- Step 3: Determine artificial viscosity for this element.      ---*/
    /*---------------------------------------------------------------------*/
    if (shockExist) {
      // Following if-else clause is purely empirical from NACA0012 case.
      // Need to develop thorough method for general problems
      if ( nPoly == 1) {
        sensorLowerBound = -6.0;
      }
      else if ( nPoly == 2 ) {
        sensorLowerBound = -12.0;
      }
      else if ( nPoly == 3 ) {
        sensorLowerBound = -12.0;
      }
      else if ( nPoly == 4 ) {
         sensorLowerBound = -17.0;
      }

      // Assign artificial viscosity based on shockSensorValue
      if ( volElem[l].shockSensorValue > sensorLowerBound ) {
        // Following value is initial guess.
        volElem[l].shockArtificialViscosity = 1.e-10;
      }
      else {
        volElem[l].shockArtificialViscosity = 0.0;
      }
    }
    else {
      volElem[l].shockArtificialViscosity = 0.0;
    }
  }

}

void CFEM_DG_NSSolver::Volume_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh, unsigned short iStep) {

  /* Start the MPI communication of the solution in the halo elements. */
  Initiate_MPI_Communication();

  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Set the pointers for the local arrays. ---*/
  su2double tick = 0.0;

  su2double *solAndGradInt = VecTmpMemory.data();
  su2double *fluxes        = solAndGradInt + nIntegrationMax*nVar*(ctc::nDim+1);

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = ctc::nDim*ctc::nDim + 1;

  /*--- Loop over the owned volume elements to compute the contribution of the
        volume integral in the DG FEM formulation to the residual.       ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Get the data from the corresponding standard element. */
    const unsigned short ind             = volElem[l].indStandardElement;
    const unsigned short nInt            = standardElementsSol[ind].GetNIntegration();
    const unsigned short nDOFs           = volElem[l].nDOFsSol;
    const su2double *matBasisInt         = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
    const su2double *matDerBasisIntTrans = standardElementsSol[ind].GetDerMatBasisFunctionsIntTrans();
    const su2double *weights             = standardElementsSol[ind].GetWeightsIntegration();

    unsigned short nPoly = standardElementsSol[ind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    /* Compute the length scale of the current element for the LES. */
    const su2double lenScale = volElem[l].lenScale/nPoly;

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Determine the solution variables and their gradients     ---*/
    /*---         w.r.t. the parametric coordinates in the integration     ---*/
    /*---         points of the element.                                   ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the solution variables for this element. */
    su2double *solDOFs = VecSolDOFs.data() + ctc::nVar*volElem[l].offsetDOFsSolLocal;

    /* Call the general function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nInt*(nDim+1), nVar, nDOFs, matBasisInt, solDOFs, solAndGradInt);
    config->GEMM_Tock(tick, "Volume_Residual1", nInt*(nDim+1), nVar, nDOFs);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the total fluxes (inviscid fluxes minus the      ---*/
    /*---         viscous fluxes), multiplied by minus the integration     ---*/
    /*---         weight, in the integration points.                       ---*/
    /*------------------------------------------------------------------------*/

    config->Tick(&tick);
    /* Determine the offset between the solution variables and the r-derivatives,
       which is also the offset between the r- and s-derivatives and the offset
       between s- and t-derivatives. */
    const unsigned short offDeriv = ctc::nVar*nInt;

    /* Make a distinction between two and three space dimensions
        in order to have the most efficient code. */
    switch( ctc::nDim ) {

      case 2: {

        /* 2D simulation. Loop over the integration points to compute
           the fluxes. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the metric terms in this integration point and
             compute the inverse of the Jacobian. */
          const su2double *metricTerms = volElem[l].metricTerms.data()
                                       + i*nMetricPerPoint;
          const su2double Jac          = metricTerms[0];
          const su2double JacInv       = 1.0/Jac;

          /* Compute the true metric terms in this integration point. */
          const su2double drdx = JacInv*metricTerms[1];
          const su2double drdy = JacInv*metricTerms[2];

          const su2double dsdx = JacInv*metricTerms[3];
          const su2double dsdy = JacInv*metricTerms[4];

          /* Compute the metric terms multiplied by minus the integration weight.
             The minus sign comes from the integration by parts in the weak
             formulation. */
          const su2double wDrdx = -weights[i]*metricTerms[1];
          const su2double wDrdy = -weights[i]*metricTerms[2];

          const su2double wDsdx = -weights[i]*metricTerms[3];
          const su2double wDsdy = -weights[i]*metricTerms[4];

          /* Easier storage of the location where the solution data of this
             integration point starts. */
          const su2double *sol    = solAndGradInt + ctc::nVar*i;
          const su2double *dSolDr = sol    + offDeriv;
          const su2double *dSolDs = dSolDr + offDeriv;

          /*--- Compute the Cartesian gradients of the independent solution
                variables from the gradients in parametric coordinates and the
                metric terms in this integration point. ---*/
          const su2double dRhoDx  = dSolDr[0]*drdx + dSolDs[0]*dsdx;
          const su2double dRhoUDx = dSolDr[1]*drdx + dSolDs[1]*dsdx;
          const su2double dRhoVDx = dSolDr[2]*drdx + dSolDs[2]*dsdx;
          const su2double dRhoEDx = dSolDr[3]*drdx + dSolDs[3]*dsdx;

          const su2double dRhoDy  = dSolDr[0]*drdy + dSolDs[0]*dsdy;
          const su2double dRhoUDy = dSolDr[1]*drdy + dSolDs[1]*dsdy;
          const su2double dRhoVDy = dSolDr[2]*drdy + dSolDs[2]*dsdy;
          const su2double dRhoEDy = dSolDr[3]*drdy + dSolDs[3]*dsdy;

          /*--- Compute the velocities and static energy in this integration point. ---*/
          const su2double rhoInv       = 1.0/sol[0];
          const su2double u            = sol[1]*rhoInv;
          const su2double v            = sol[2]*rhoInv;
          const su2double TotalEnergy  = sol[3]*rhoInv;
          const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

          /*--- Compute the Cartesian gradients of the velocities and static energy
                in this integration point and also the divergence of the velocity. ---*/
          const su2double dudx = rhoInv*(dRhoUDx - u*dRhoDx);
          const su2double dudy = rhoInv*(dRhoUDy - u*dRhoDy);

          const su2double dvdx = rhoInv*(dRhoVDx - v*dRhoDx);
          const su2double dvdy = rhoInv*(dRhoVDy - v*dRhoDy);

          const su2double dStaticEnergyDx = rhoInv*(dRhoEDx - TotalEnergy*dRhoDx)
                                          - u*dudx - v*dvdx;
          const su2double dStaticEnergyDy = rhoInv*(dRhoEDy - TotalEnergy*dRhoDy)
                                          - u*dudy - v*dvdy;

          const su2double divVel = dudx + dvdy;

          /*--- Compute the pressure and the laminar viscosity. ---*/
          FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
          const su2double Pressure     = FluidModel->GetPressure();
          const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

          /*--- If an SGS model is used the eddy viscosity must be computed. ---*/
          su2double ViscosityTurb = 0.0;
          if( SGSModelUsed ) {
            const su2double dist = volElem[l].wallDistance[i];
            su2double velGrad[3][3];
            velGrad[0][0] = dudx; velGrad[0][1] = dudy; velGrad[0][2] = 0.0;
            velGrad[1][0] = dvdx; velGrad[1][1] = dvdy; velGrad[1][2] = 0.0;
            velGrad[2][0] = 0.0;  velGrad[2][1] = 0.0;  velGrad[2][2] = 0.0;

            ViscosityTurb = SGSModel->ComputeEddyViscosity(nDim, sol[0], velGrad,
                                                           lenScale, dist);
          }

          /* Compute the total viscosity and heat conductivity. Note that the heat
             conductivity is divided by the Cv, because gradients of internal energy
             are computed and not temperature. */
          const su2double Viscosity = ViscosityLam + ViscosityTurb;
          const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                                  + ViscosityTurb*factHeatFlux_Turb;

          /*--- Set the value of the second viscosity and compute the divergence
                term in the viscous normal stresses. ---*/
          const su2double lambda     = -TWO3*Viscosity;
          const su2double lamDivTerm =  lambda*divVel;

          /*--- Compute the viscous stress tensor and minus the heatflux vector. ---*/
          const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
          const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
          const su2double tauxy = Viscosity*(dudy + dvdx);

          const su2double qx = kOverCv*dStaticEnergyDx;
          const su2double qy = kOverCv*dStaticEnergyDy;

          /* Set the pointer for the fluxes in this integration point. */
          su2double *flux = fluxes + i*ctc::nDim*ctc::nVar;

          /*--- Fluxes in r-direction. */
          const su2double H     = TotalEnergy + rhoInv*Pressure;
          const su2double rhoUr = sol[1]*wDrdx + sol[2]*wDrdy;

          flux[0] = rhoUr;
          flux[1] = rhoUr*u + (Pressure - tauxx)*wDrdx - tauxy*wDrdy;
          flux[2] = rhoUr*v - tauxy*wDrdx + (Pressure - tauyy)*wDrdy;
          flux[3] = rhoUr*H - (u*tauxx + v*tauxy + qx)*wDrdx
                            - (u*tauxy + v*tauyy + qy)*wDrdy;

          /*--- Fluxes in s-direction. */
          const su2double rhoUs = sol[1]*wDsdx + sol[2]*wDsdy;

          flux[4] = rhoUs;
          flux[5] = rhoUs*u + (Pressure - tauxx)*wDsdx - tauxy*wDsdy;
          flux[6] = rhoUs*v - tauxy*wDsdx + (Pressure - tauyy)*wDsdy;
          flux[7] = rhoUs*H - (u*tauxx + v*tauxy + qx)*wDsdx
                            - (u*tauxy + v*tauyy + qy)*wDsdy;
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /* 3D simulation. Loop over the integration points to compute
           the fluxes. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the metric terms in this integration point and
             compute the inverse of the Jacobian. */
          const su2double *metricTerms = volElem[l].metricTerms.data()
                                       + i*nMetricPerPoint;
          const su2double Jac          = metricTerms[0];
          const su2double JacInv       = 1.0/Jac;

          /* Compute the true metric terms in this integration point. */
          const su2double drdx = JacInv*metricTerms[1];
          const su2double drdy = JacInv*metricTerms[2];
          const su2double drdz = JacInv*metricTerms[3];

          const su2double dsdx = JacInv*metricTerms[4];
          const su2double dsdy = JacInv*metricTerms[5];
          const su2double dsdz = JacInv*metricTerms[6];

          const su2double dtdx = JacInv*metricTerms[7];
          const su2double dtdy = JacInv*metricTerms[8];
          const su2double dtdz = JacInv*metricTerms[9];

          /* Compute the metric terms multiplied by minus the integration weight.
             The minus sign comes from the integration by parts in the weak
             formulation. */
          const su2double wDrdx = -weights[i]*metricTerms[1];
          const su2double wDrdy = -weights[i]*metricTerms[2];
          const su2double wDrdz = -weights[i]*metricTerms[3];

          const su2double wDsdx = -weights[i]*metricTerms[4];
          const su2double wDsdy = -weights[i]*metricTerms[5];
          const su2double wDsdz = -weights[i]*metricTerms[6];

          const su2double wDtdx = -weights[i]*metricTerms[7];
          const su2double wDtdy = -weights[i]*metricTerms[8];
          const su2double wDtdz = -weights[i]*metricTerms[9];

          /* Easier storage of the location where the solution data of this
             integration point starts. */
          const su2double *sol    = solAndGradInt + ctc::nVar*i;
          const su2double *dSolDr = sol    + offDeriv;
          const su2double *dSolDs = dSolDr + offDeriv;
          const su2double *dSolDt = dSolDs + offDeriv;

          /*--- Compute the Cartesian gradients of the independent solution
                variables from the gradients in parametric coordinates and the
                metric terms in this integration point. ---*/
          const su2double dRhoDx  = dSolDr[0]*drdx + dSolDs[0]*dsdx + dSolDt[0]*dtdx;
          const su2double dRhoUDx = dSolDr[1]*drdx + dSolDs[1]*dsdx + dSolDt[1]*dtdx;
          const su2double dRhoVDx = dSolDr[2]*drdx + dSolDs[2]*dsdx + dSolDt[2]*dtdx;
          const su2double dRhoWDx = dSolDr[3]*drdx + dSolDs[3]*dsdx + dSolDt[3]*dtdx;
          const su2double dRhoEDx = dSolDr[4]*drdx + dSolDs[4]*dsdx + dSolDt[4]*dtdx;

          const su2double dRhoDy  = dSolDr[0]*drdy + dSolDs[0]*dsdy + dSolDt[0]*dtdy;
          const su2double dRhoUDy = dSolDr[1]*drdy + dSolDs[1]*dsdy + dSolDt[1]*dtdy;
          const su2double dRhoVDy = dSolDr[2]*drdy + dSolDs[2]*dsdy + dSolDt[2]*dtdy;
          const su2double dRhoWDy = dSolDr[3]*drdy + dSolDs[3]*dsdy + dSolDt[3]*dtdy;
          const su2double dRhoEDy = dSolDr[4]*drdy + dSolDs[4]*dsdy + dSolDt[4]*dtdy;

          const su2double dRhoDz  = dSolDr[0]*drdz + dSolDs[0]*dsdz + dSolDt[0]*dtdz;
          const su2double dRhoUDz = dSolDr[1]*drdz + dSolDs[1]*dsdz + dSolDt[1]*dtdz;
          const su2double dRhoVDz = dSolDr[2]*drdz + dSolDs[2]*dsdz + dSolDt[2]*dtdz;
          const su2double dRhoWDz = dSolDr[3]*drdz + dSolDs[3]*dsdz + dSolDt[3]*dtdz;
          const su2double dRhoEDz = dSolDr[4]*drdz + dSolDs[4]*dsdz + dSolDt[4]*dtdz;

          /*--- Compute the velocities and static energy in this integration point. ---*/
          const su2double rhoInv       = 1.0/sol[0];
          const su2double u            = sol[1]*rhoInv;
          const su2double v            = sol[2]*rhoInv;
          const su2double w            = sol[3]*rhoInv;
          const su2double TotalEnergy  = sol[4]*rhoInv;
          const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v + w*w);

          /*--- Compute the Cartesian gradients of the velocities and static energy
                in this integration point and also the divergence of the velocity. ---*/
          const su2double dudx = rhoInv*(dRhoUDx - u*dRhoDx);
          const su2double dudy = rhoInv*(dRhoUDy - u*dRhoDy);
          const su2double dudz = rhoInv*(dRhoUDz - u*dRhoDz);

          const su2double dvdx = rhoInv*(dRhoVDx - v*dRhoDx);
          const su2double dvdy = rhoInv*(dRhoVDy - v*dRhoDy);
          const su2double dvdz = rhoInv*(dRhoVDz - v*dRhoDz);

          const su2double dwdx = rhoInv*(dRhoWDx - w*dRhoDx);
          const su2double dwdy = rhoInv*(dRhoWDy - w*dRhoDy);
          const su2double dwdz = rhoInv*(dRhoWDz - w*dRhoDz);

          const su2double dStaticEnergyDx = rhoInv*(dRhoEDx - TotalEnergy*dRhoDx)
                                          - u*dudx - v*dvdx - w*dwdx;
          const su2double dStaticEnergyDy = rhoInv*(dRhoEDy - TotalEnergy*dRhoDy)
                                          - u*dudy - v*dvdy - w*dwdy;
          const su2double dStaticEnergyDz = rhoInv*(dRhoEDz - TotalEnergy*dRhoDz)
                                          - u*dudz - v*dvdz - w*dwdz;

          const su2double divVel = dudx + dvdy + dwdz;

          /*--- Compute the pressure and the laminar viscosity. ---*/
          FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
          const su2double Pressure     = FluidModel->GetPressure();
          const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

          /*--- If an SGS model is used the eddy viscosity must be computed. ---*/
          su2double ViscosityTurb = 0.0;
          if( SGSModelUsed ) {
            const su2double dist = volElem[l].wallDistance[i];
            su2double velGrad[3][3];
            velGrad[0][0] = dudx; velGrad[0][1] = dudy; velGrad[0][2] = dudz;
            velGrad[1][0] = dvdx; velGrad[1][1] = dvdy; velGrad[1][2] = dvdz;
            velGrad[2][0] = dwdx; velGrad[2][1] = dwdy; velGrad[2][2] = dwdz;

            ViscosityTurb = SGSModel->ComputeEddyViscosity(nDim, sol[0], velGrad,
                                                           lenScale, dist);
          }

          /* Compute the total viscosity and heat conductivity. Note that the heat
             conductivity is divided by the Cv, because gradients of internal energy
             are computed and not temperature. */
          const su2double Viscosity = ViscosityLam + ViscosityTurb;
          const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                                  + ViscosityTurb*factHeatFlux_Turb;

          /*--- Set the value of the second viscosity and compute the divergence
                term in the viscous normal stresses. ---*/
          const su2double lambda     = -TWO3*Viscosity;
          const su2double lamDivTerm =  lambda*divVel;

          /*--- Compute the viscous stress tensor and minus the heatflux vector. ---*/
          const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
          const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
          const su2double tauzz = 2.0*Viscosity*dwdz + lamDivTerm;

          const su2double tauxy = Viscosity*(dudy + dvdx);
          const su2double tauxz = Viscosity*(dudz + dwdx);
          const su2double tauyz = Viscosity*(dvdz + dwdy);

          const su2double qx = kOverCv*dStaticEnergyDx;
          const su2double qy = kOverCv*dStaticEnergyDy;
          const su2double qz = kOverCv*dStaticEnergyDz;

          /* Set the pointer for the fluxes in this integration point. */
          su2double *flux = fluxes + i*ctc::nDim*ctc::nVar;

          /*--- Fluxes in r-direction. */
          const su2double H     = TotalEnergy + rhoInv*Pressure;
          const su2double rhoUr = sol[1]*wDrdx + sol[2]*wDrdy + sol[3]*wDrdz;

          flux[0] = rhoUr;
          flux[1] = rhoUr*u + (Pressure - tauxx)*wDrdx - tauxy*wDrdy - tauxz*wDrdz;
          flux[2] = rhoUr*v - tauxy*wDrdx + (Pressure - tauyy)*wDrdy - tauyz*wDrdz;
          flux[3] = rhoUr*w - tauxz*wDrdx - tauyz*wDrdy + (Pressure - tauzz)*wDrdz;
          flux[4] = rhoUr*H - (u*tauxx + v*tauxy + w*tauxz + qx)*wDrdx
                            - (u*tauxy + v*tauyy + w*tauyz + qy)*wDrdy
                            - (u*tauxz + v*tauyz + w*tauzz + qz)*wDrdz;

          /*--- Fluxes in s-direction. */
          const su2double rhoUs = sol[1]*wDsdx + sol[2]*wDsdy + sol[3]*wDsdz;

          flux[5] = rhoUs;
          flux[6] = rhoUs*u + (Pressure - tauxx)*wDsdx - tauxy*wDsdy - tauxz*wDsdz;
          flux[7] = rhoUs*v - tauxy*wDsdx + (Pressure - tauyy)*wDsdy - tauyz*wDsdz;
          flux[8] = rhoUs*w - tauxz*wDsdx - tauyz*wDsdy + (Pressure - tauzz)*wDsdz;
          flux[9] = rhoUs*H - (u*tauxx + v*tauxy + w*tauxz + qx)*wDsdx
                            - (u*tauxy + v*tauyy + w*tauyz + qy)*wDsdy
                            - (u*tauxz + v*tauyz + w*tauzz + qz)*wDsdz;

          /*--- Fluxes in t-direction. */
          const su2double rhoUt = sol[1]*wDtdx + sol[2]*wDtdy + sol[3]*wDtdz;

          flux[10] = rhoUt;
          flux[11] = rhoUt*u + (Pressure - tauxx)*wDtdx - tauxy*wDtdy - tauxz*wDtdz;
          flux[12] = rhoUt*v - tauxy*wDtdx + (Pressure - tauyy)*wDtdy - tauyz*wDtdz;
          flux[13] = rhoUt*w - tauxz*wDtdx - tauyz*wDtdy + (Pressure - tauzz)*wDtdz;
          flux[14] = rhoUt*H - (u*tauxx + v*tauxy + w*tauxz + qx)*wDtdx
                             - (u*tauxy + v*tauyy + w*tauyz + qy)*wDtdy
                             - (u*tauxz + v*tauyz + w*tauzz + qz)*wDtdz;
        }

        break;
      }
    }

    config->Tock(tick, "IR_2_1", 4);

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
    /*---         integration over the volume element.                     ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the residuals for this volume element. */
    su2double *res = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Call the general function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nDOFs, nVar, nInt*nDim, matDerBasisIntTrans, fluxes, res);
    config->GEMM_Tock(tick, "Volume_Residual2", nDOFs, nVar, nInt*nDim);

  }
}

void CFEM_DG_NSSolver::ResidualFaces(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iStep,
                                     const unsigned long indFaceBeg, const unsigned long indFaceEnd,
                                     unsigned long &indResFaces) {

  /* Determine whether or not the Cartesian gradients of the basis functions
     are stored. */
  const bool CartGradBasisFunctionsStored = config->GetStore_Cart_Grad_BasisFunctions_DGFEM();

  /*--- Set the pointers for the local arrays. ---*/
  su2double tick = 0.0;
  su2double tick2 = 0.0;

  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

  su2double *solIntL       = VecTmpMemory.data();
  su2double *solIntR       = solIntL       + nIntegrationMax*nVar;
  su2double *viscosityIntL = solIntR       + nIntegrationMax*nVar;
  su2double *kOverCvIntL   = viscosityIntL + nIntegrationMax;
  su2double *viscosityIntR = kOverCvIntL   + nIntegrationMax;
  su2double *kOverCvIntR   = viscosityIntR + nIntegrationMax;
  su2double *gradSolInt    = kOverCvIntR   + nIntegrationMax;
  su2double *fluxes        = gradSolInt    + sizeGradSolInt;
  su2double *viscFluxes    = fluxes        + sizeFluxes;

  /*--- Loop over the requested range of matching faces. ---*/
  for(unsigned long l=indFaceBeg; l<indFaceEnd; ++l) {

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Compute the inviscid fluxes in the integration points of ---*/
    /*---         this matching face.                                      ---*/
    /*------------------------------------------------------------------------*/

    config->Tick(&tick);
    /* Compute the inviscid fluxes in the integration points. */
    InviscidFluxesInternalMatchingFace(config, &matchingInternalFaces[l],
                                       solIntL, solIntR, fluxes, numerics);
    config->Tock(tick, "ER_1_1", 4);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the viscous fluxes in the integration points of  ---*/
    /*---         this matching face and subtract them from the already    ---*/
    /*---         computed inviscid fluxes.                                ---*/
    /*------------------------------------------------------------------------*/

    config->Tick(&tick);
    /*--- Get the information from the standard face to compute the viscous
          fluxes in the integration points on the left side, i.e. side 0. ---*/
    const unsigned short ind          = matchingInternalFaces[l].indStandardElement;
    const unsigned short nInt         = standardMatchingFacesSol[ind].GetNIntegration();
          unsigned short nDOFsElem    = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
    const su2double     *derBasisElem = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationSide0();

    /* Get the length scales of the adjacent elements. */
    const su2double lenScale0 = volElem[matchingInternalFaces[l].elemID0].lenScale;
    const su2double lenScale1 = volElem[matchingInternalFaces[l].elemID1].lenScale;

    /* Compute the length scale for the LES of the element of side 0. */
    unsigned short iind  = volElem[matchingInternalFaces[l].elemID0].indStandardElement;
    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale0_LES = lenScale0/nPoly;

    /* Call the general function to compute the viscous flux in normal
       direction for side 0. */
    ViscousNormalFluxFace(config, nInt, nDOFsElem, 0.0, false, derBasisElem, solIntL,
                          matchingInternalFaces[l].DOFsSolElementSide0.data(),
                          matchingInternalFaces[l].metricCoorDerivFace0.data(),
                          matchingInternalFaces[l].metricNormalsFace.data(),
                          matchingInternalFaces[l].wallDistance.data(), lenScale0_LES,
                          gradSolInt, viscFluxes, viscosityIntL, kOverCvIntL);

    /*--- Subtract half of the viscous fluxes from the inviscid fluxes. The
          factor 0.5 comes from the fact that the average of the viscous fluxes
          of side 0 and side 1 must be taken in the DG-FEM formulation. ---*/
    for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] -= 0.5*viscFluxes[j];

    /*--- Get the information from the standard face to compute the viscous
          fluxes in the integration points on the right side, i.e. side 1. ---*/
    nDOFsElem    = standardMatchingFacesSol[ind].GetNDOFsElemSide1();
    derBasisElem = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationSide1();

    /* Compute the length scale for the LES of the element of side 1. */
    iind  = volElem[matchingInternalFaces[l].elemID1].indStandardElement;
    nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale1_LES = lenScale1/nPoly;

    /* Call the general function to compute the viscous flux in normal
       direction for side 1. */
    ViscousNormalFluxFace(config, nInt, nDOFsElem, 0.0, false, derBasisElem, solIntR,
                          matchingInternalFaces[l].DOFsSolElementSide1.data(),
                          matchingInternalFaces[l].metricCoorDerivFace1.data(),
                          matchingInternalFaces[l].metricNormalsFace.data(),
                          matchingInternalFaces[l].wallDistance.data(), lenScale1_LES,
                          gradSolInt, viscFluxes, viscosityIntR, kOverCvIntR);

    /*--- Subtract half of the viscous fluxes from the inviscid fluxes. ---*/
    for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] -= 0.5*viscFluxes[j];

    config->Tock(tick, "ER_1_2", 4);

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the penalty terms in the integration points of   ---*/
    /*---         this matching face and them to the already stored        ---*/
    /*---         inviscid and viscous fluxes.                             ---*/
    /*------------------------------------------------------------------------*/

    config->Tick(&tick);
    /* Get the required constant needed for the penalty terms. */
    const su2double ConstPenFace = standardMatchingFacesSol[ind].GetPenaltyConstant();

    /* Call the function PenaltyTermsFluxFace to compute the actual penalty
       terms. Use the array viscFluxes as storage. */
    PenaltyTermsFluxFace(nInt, solIntL, solIntR, viscosityIntL, viscosityIntR,
                         kOverCvIntL, kOverCvIntR, ConstPenFace, lenScale0, lenScale1,
                         matchingInternalFaces[l].metricNormalsFace.data(),
                         viscFluxes);

    /* Add the penalty fluxes to the earlier computed fluxes. */
    for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] += viscFluxes[j];

    /* Multiply the fluxes with the integration weight of the corresponding
       integration point. */
    const su2double *weights = standardMatchingFacesSol[ind].GetWeightsIntegration();
    for(unsigned short i=0; i<nInt; ++i) {
      su2double *flux = fluxes + i*ctc::nVar;

      for(unsigned short j=0; j<ctc::nVar; ++j)
        flux[j] *= weights[i];
    }
    config->Tock(tick, "ER_1_3", 4);

    /*------------------------------------------------------------------------*/
    /*--- Step 4: Compute the contribution to the residuals from the       ---*/
    /*---         integration over this internal matching face.            ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the position in the residual array for side 0 of
       this face and update the corresponding counter. */
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    su2double *resFace0 = VecResFaces.data() + indResFaces*ctc::nVar;
    indResFaces        += nDOFsFace0;

    /* Get the correct form of the basis functions needed for the matrix
       multiplication to compute the residual. */
    const su2double *basisFaceTrans = standardMatchingFacesSol[ind].GetBasisFaceIntegrationTransposeSide0();

    /* Call the general function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nDOFsFace0, nVar, nInt, basisFaceTrans, fluxes, resFace0);
    config->GEMM_Tock(tick, "ResidualFaces1", nDOFsFace0, nVar, nInt);

    /* Easier storage of the position in the residual array for side 1 of
       this face and update the corresponding counter. */
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();
    su2double *resFace1 = VecResFaces.data() + indResFaces*nVar;
    indResFaces        += nDOFsFace1;

    /* Check if the number of DOFs on side 1 is equal to the number of DOFs
       of side 0. In that case the residual for side 1 is obtained by simply
       negating the data from side 0. */
    if(nDOFsFace1 == nDOFsFace0) {
      for(unsigned short i=0; i<(nVar*nDOFsFace1); ++i)
        resFace1[i] = -resFace0[i];
    }
    else {

      /*--- The number of DOFs and hence the polynomial degree of side 1 is
            different from side. Carry out the matrix multiplication to obtain
            the residual. Afterwards the residual is negated, because the
            normal is pointing into the adjacent element. ---*/
      basisFaceTrans = standardMatchingFacesSol[ind].GetBasisFaceIntegrationTransposeSide1();
      config->GEMM_Tick(&tick);
      DenseMatrixProduct(nDOFsFace1, nVar, nInt, basisFaceTrans, fluxes, resFace1);
      config->GEMM_Tock(tick, "ResidualFaces2", nDOFsFace1, nVar, nInt);

      for(unsigned short i=0; i<(nVar*nDOFsFace1); ++i)
        resFace1[i] = -resFace1[i];
    }
    config->Tock(tick, "ER_1_4", 4);

    /*------------------------------------------------------------------------*/
    /*--- Step 5: Compute the symmetrizing terms, if present, in the       ---*/
    /*---         integration points of this matching face.                ---*/
    /*------------------------------------------------------------------------*/

    if( symmetrizingTermsPresent ) {

      config->Tick(&tick);
      /* Compute the symmetrizing fluxes in the nDim directions. */
      SymmetrizingFluxesFace(nInt, solIntL, solIntR, viscosityIntL, viscosityIntR,
                             kOverCvIntL, kOverCvIntR,
                             matchingInternalFaces[l].metricNormalsFace.data(),
                             fluxes);

      /*--- Multiply the fluxes just computed by their integration weights and
            -theta/2. The parameter theta is the parameter in the Interior Penalty
            formulation, the factor 1/2 comes in from the averaging and the minus
            sign is from the convention that the viscous fluxes come with a minus
            sign in this code. ---*/
      const su2double halfTheta = 0.5*config->GetTheta_Interior_Penalty_DGFEM();

      for(unsigned short i=0; i<nInt; ++i) {
        su2double *flux        = fluxes + i*ctc::nVar*ctc::nDim;
        const su2double wTheta = -halfTheta*weights[i];

        for(unsigned short j=0; j<(ctc::nVar*ctc::nDim); ++j)
          flux[j] *= wTheta;
      }
      config->Tock(tick, "ER_1_5", 4);

      /*------------------------------------------------------------------------*/
      /*--- Step 6: Distribute the symmetrizing terms to the DOFs. Note that ---*/
      /*---         these terms must be distributed to all the DOFs of the   ---*/
      /*---         adjacent elements, not only to the DOFs of the face.     ---*/
      /*------------------------------------------------------------------------*/

      /* Easier storage of the number of DOFs for the adjacent elements. */
      const unsigned short nDOFsElem0 = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
      const unsigned short nDOFsElem1 = standardMatchingFacesSol[ind].GetNDOFsElemSide1();

      config->Tick(&tick);

      /* The Cartesian gradients of the element basis functions of side 0
         are needed. Check if this data is stored. */
      const su2double *cartGrad;
      if( CartGradBasisFunctionsStored )
        cartGrad = matchingInternalFaces[l].metricElemSide0.data();
      else {

        /* The gradients are not stored. Use the array gradSolInt to store
           these derivatives. Get the derivatives w.r.t. the parametric
           coordinates of the element on side 0 of the face. */
        const su2double *derBasisElemTrans = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationTransposeSide0();

        /*--- Create the Cartesian derivatives of the basis functions
              in the integration points. ---*/
        unsigned int ii = 0;
        for(unsigned short j=0; j<nDOFsElem0; ++j) {
          for(unsigned short i=0; i<nInt; ++i, ii+=ctc::nDim) {

            /* Easier storage of the derivatives of the basis function w.r.t. the
               parametric coordinates, the location where to store the Cartesian
               derivatives of the basis functions, and the metric terms in this
               integration point. */
            const su2double *derParam    = derBasisElemTrans + ii;
            const su2double *metricTerms = matchingInternalFaces[l].metricCoorDerivFace0.data()
                                         + i*ctc::nDim*ctc::nDim;
                  su2double *derCar      = gradSolInt + ii;

            /*--- Loop over the dimensions to compute the Cartesian derivatives
                  of the basis functions. ---*/
#pragma simd 
            for(unsigned short k=0; k<ctc::nDim; ++k) {
              derCar[k] = 0.0;
              for(unsigned short l=0; l<ctc::nDim; ++l)
                derCar[k] += derParam[l]*metricTerms[k+l*ctc::nDim];
            }
          }
        }

        /* Set the pointer of gradSolInt to cartGrad, such that the latter can
           be used in the call to the matrix multiplication. */
        cartGrad = gradSolInt;
      }

      /* Set the pointer where to store the current residual and update the
         counter indResFaces. */
      su2double *resElem0 = VecResFaces.data() + indResFaces*nVar;
      indResFaces        += nDOFsElem0;

      /* Call the general function to carry out the matrix product to compute
         the residual for side 0. */
      config->GEMM_Tick(&tick2);
      DenseMatrixProduct(nDOFsElem0, nVar, nInt*nDim, cartGrad, fluxes, resElem0);
      config->GEMM_Tock(tick2, "ResidualFaces3", nDOFsElem0, nVar, nInt*nDim);

      /* The Cartesian gradients of the element basis functions of side 1
         are needed. Check if this data is stored. */
      if( CartGradBasisFunctionsStored )
        cartGrad = matchingInternalFaces[l].metricElemSide1.data();
      else {

        /* The Cartesian gradients must be computed. Get the derivatives w.r.t.
           the parametric coordinates of the element on side 1 of the face. */
        const su2double *derBasisElemTrans = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationTransposeSide1();

        /*--- Create the Cartesian derivatives of the basis functions in the
              integration points. Use gradSolInt for storage. ---*/
        unsigned int ii = 0;
        for(unsigned short j=0; j<nDOFsElem1; ++j) {
          for(unsigned short i=0; i<nInt; ++i, ii+=ctc::nDim) {

            /* Easier storage of the derivatives of the basis function w.r.t. the
               parametric coordinates, the location where to store the Cartesian
               derivatives of the basis functions, and the metric terms in this
               integration point. */
            const su2double *derParam    = derBasisElemTrans + ii;
            const su2double *metricTerms = matchingInternalFaces[l].metricCoorDerivFace1.data()
                                         + i*ctc::nDim*ctc::nDim;
                  su2double *derCar      = gradSolInt + ii;

            /*--- Loop over the dimensions to compute the Cartesian derivatives
                  of the basis functions. ---*/
#pragma simd
            for(unsigned short k=0; k<ctc::nDim; ++k) {
              derCar[k] = 0.0;
              for(unsigned short l=0; l<ctc::nDim; ++l)
                derCar[k] += derParam[l]*metricTerms[k+l*ctc::nDim];
            }
          }
        }

        /* Set the pointer of gradSolInt to cartGrad, such that the latter can
           be used in the call to the matrix multiplication. */
        cartGrad = gradSolInt;
      }

      /* Set the pointer where to store the current residual and update the
         counter indResFaces. */
      su2double *resElem1 = VecResFaces.data() + indResFaces*nVar;
      indResFaces        += nDOFsElem1;

      /* Call the general function to carry out the matrix product to compute
         the residual for side 1. Note that the symmetrizing residual should not
         be negated, because two minus signs enter the formulation for side 1,
         which cancel each other. */
      config->GEMM_Tick(&tick2);
      DenseMatrixProduct(nDOFsElem1, nVar, nInt*nDim, cartGrad, fluxes, resElem1);
      config->GEMM_Tock(tick2, "ResidualFaces4", nDOFsElem1, nVar, nInt*nDim);
    }
    config->Tock(tick, "ER_1_6", 4);
  }
}

void CFEM_DG_NSSolver::ViscousNormalFluxFace(CConfig              *config,
                                             const unsigned short nInt,
                                             const unsigned short nDOFsElem,
                                             const su2double      Wall_HeatFlux,
                                             const bool           HeatFlux_Prescribed,
                                             const su2double      *derBasisElem,
                                             const su2double      *solInt,
                                             const unsigned long  *DOFsElem,
                                             const su2double      *metricCoorDerivFace,
                                             const su2double      *metricNormalsFace,
                                             const su2double      *wallDistanceInt,
                                             const su2double      lenScale_LES,
                                                   su2double      *gradSolInt,
                                                   su2double      *viscNormFluxes,
                                                   su2double      *viscosityInt,
                                                   su2double      *kOverCvInt) {

  su2double tick = 0.0;

  /* Constant factor present in the heat flux vector. Set it to zero if the heat
     flux is prescribed, such that no if statements are needed in the loop. */
  const su2double factHeatFlux_Lam  = HeatFlux_Prescribed ? 0.0: Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = HeatFlux_Prescribed ? 0.0: Gamma/Prandtl_Turb;

  /* Set the value of the prescribed heat flux for the same reason. */
  const su2double HeatFlux = HeatFlux_Prescribed ? Wall_HeatFlux : 0.0;

  /* Set the pointer solElem to viscNormFluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the element and the fluxes to be computed. */
  su2double *solElem = viscNormFluxes;

  /*--- Store the solution of the DOFs of the adjacent element in contiguous
        memory such that the function DenseMatrixProduct can be used to compute
        the gradients solution variables in the integration points of the face. ---*/
  for(unsigned short i=0; i<nDOFsElem; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + ctc::nVar*DOFsElem[i];
    su2double       *sol    = solElem + ctc::nVar*i;
    for(unsigned short j=0; j<ctc::nVar; ++j)
      sol[j] = solDOF[j];
  }

  /* Compute the gradients in the integration points. Call the general function to
     carry out the matrix product. */
  config->GEMM_Tick(&tick);
  DenseMatrixProduct(nInt*nDim, nVar, nDOFsElem, derBasisElem, solElem, gradSolInt);
  config->GEMM_Tock(tick, "ViscousNormalFluxFace", nInt*nDim, nVar, nDOFsElem);

  /* Determine the offset between r- and -s-derivatives, which is also the
     offset between s- and t-derivatives. */
  const unsigned short offDeriv = ctc::nVar*nInt;

  /* Make a distinction between two and three space dimensions
     in order to have the most efficient code. */
  switch( ctc::nDim ) {

    case 2: {

      /* 2D simulation. Loop over the integration points to
         compute the viscous fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the metric terms needed to compute the Cartesian
           gradients in this integration point and the starting locations of
           the solution and the gradients, w.r.t. the parametric coordinates
           of this solution. */
        const su2double *metricTerms = metricCoorDerivFace + i*ctc::nDim*ctc::nDim;
        const su2double *sol         = solInt     + ctc::nVar*i;
        const su2double *dSolDr      = gradSolInt + ctc::nVar*i;
        const su2double *dSolDs      = dSolDr     + offDeriv;

        /* Easier storage of the metric terms in this integration point. */
        const su2double drdx = metricTerms[0];
        const su2double drdy = metricTerms[1];

        const su2double dsdx = metricTerms[2];
        const su2double dsdy = metricTerms[3];

        /*--- Compute the Cartesian gradients of the solution. ---*/
        su2double solGradCart[4][2];

        solGradCart[0][0] = dSolDr[0]*drdx + dSolDs[0]*dsdx;
        solGradCart[1][0] = dSolDr[1]*drdx + dSolDs[1]*dsdx;
        solGradCart[2][0] = dSolDr[2]*drdx + dSolDs[2]*dsdx;
        solGradCart[3][0] = dSolDr[3]*drdx + dSolDs[3]*dsdx;

        solGradCart[0][1] = dSolDr[0]*drdy + dSolDs[0]*dsdy;
        solGradCart[1][1] = dSolDr[1]*drdy + dSolDs[1]*dsdy;
        solGradCart[2][1] = dSolDr[2]*drdy + dSolDs[2]*dsdy;
        solGradCart[3][1] = dSolDr[3]*drdy + dSolDs[3]*dsdy;

        /*--- Call the function ViscousNormalFluxIntegrationPoint to compute the
              actual normal viscous flux. The viscosity and thermal conductivity
              are stored for later use. ---*/
        const su2double *normal  = metricNormalsFace + i*(ctc::nDim+1);
        su2double *normalFlux    = viscNormFluxes + i*ctc::nVar;
        const su2double wallDist = wallDistanceInt ? wallDistanceInt[i] : 0.0;

        su2double Viscosity, kOverCv;

        ViscousNormalFluxIntegrationPoint_2D(sol, solGradCart, normal, HeatFlux,
                                             factHeatFlux_Lam, factHeatFlux_Turb,
                                             wallDist, lenScale_LES,
                                             Viscosity, kOverCv, normalFlux);
        viscosityInt[i] = Viscosity;
        kOverCvInt[i]   = kOverCv;
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case 3: {

      /* 3D simulation. Loop over the integration points to
         compute the viscous fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the metric terms needed to compute the Cartesian
           gradients in this integration point and the starting locations of
           the solution and the gradients, w.r.t. the parametric coordinates
           of this solution. */
        const su2double *metricTerms = metricCoorDerivFace + i*ctc::nDim*ctc::nDim;
        const su2double *sol         = solInt     + ctc::nVar*i;
        const su2double *dSolDr      = gradSolInt + ctc::nVar*i;
        const su2double *dSolDs      = dSolDr     + offDeriv;
        const su2double *dSolDt      = dSolDs     + offDeriv;

        /* Easier storage of the metric terms in this integration point. */
        const su2double drdx = metricTerms[0];
        const su2double drdy = metricTerms[1];
        const su2double drdz = metricTerms[2];

        const su2double dsdx = metricTerms[3];
        const su2double dsdy = metricTerms[4];
        const su2double dsdz = metricTerms[5];

        const su2double dtdx = metricTerms[6];
        const su2double dtdy = metricTerms[7];
        const su2double dtdz = metricTerms[8];

        /*--- Compute the Cartesian gradients of the solution. ---*/
        su2double solGradCart[5][3];

        solGradCart[0][0] = dSolDr[0]*drdx + dSolDs[0]*dsdx + dSolDt[0]*dtdx;
        solGradCart[1][0] = dSolDr[1]*drdx + dSolDs[1]*dsdx + dSolDt[1]*dtdx;
        solGradCart[2][0] = dSolDr[2]*drdx + dSolDs[2]*dsdx + dSolDt[2]*dtdx;
        solGradCart[3][0] = dSolDr[3]*drdx + dSolDs[3]*dsdx + dSolDt[3]*dtdx;
        solGradCart[4][0] = dSolDr[4]*drdx + dSolDs[4]*dsdx + dSolDt[4]*dtdx;

        solGradCart[0][1] = dSolDr[0]*drdy + dSolDs[0]*dsdy + dSolDt[0]*dtdy;
        solGradCart[1][1] = dSolDr[1]*drdy + dSolDs[1]*dsdy + dSolDt[1]*dtdy;
        solGradCart[2][1] = dSolDr[2]*drdy + dSolDs[2]*dsdy + dSolDt[2]*dtdy;
        solGradCart[3][1] = dSolDr[3]*drdy + dSolDs[3]*dsdy + dSolDt[3]*dtdy;
        solGradCart[4][1] = dSolDr[4]*drdy + dSolDs[4]*dsdy + dSolDt[4]*dtdy;

        solGradCart[0][2] = dSolDr[0]*drdz + dSolDs[0]*dsdz + dSolDt[0]*dtdz;
        solGradCart[1][2] = dSolDr[1]*drdz + dSolDs[1]*dsdz + dSolDt[1]*dtdz;
        solGradCart[2][2] = dSolDr[2]*drdz + dSolDs[2]*dsdz + dSolDt[2]*dtdz;
        solGradCart[3][2] = dSolDr[3]*drdz + dSolDs[3]*dsdz + dSolDt[3]*dtdz;
        solGradCart[4][2] = dSolDr[4]*drdz + dSolDs[4]*dsdz + dSolDt[4]*dtdz;

        /*--- Call the function ViscousNormalFluxIntegrationPoint to compute the
              actual normal viscous flux. The viscosity and thermal conductivity
              are stored for later use. ---*/
        const su2double *normal  = metricNormalsFace + i*(ctc::nDim+1);
        su2double *normalFlux    = viscNormFluxes + i*ctc::nVar;
        const su2double wallDist = wallDistanceInt ? wallDistanceInt[i] : 0.0;

        su2double Viscosity, kOverCv;

        ViscousNormalFluxIntegrationPoint_3D(sol, solGradCart, normal, HeatFlux,
                                             factHeatFlux_Lam, factHeatFlux_Turb,
                                             wallDist, lenScale_LES,
                                             Viscosity, kOverCv, normalFlux);
        viscosityInt[i] = Viscosity;
        kOverCvInt[i]   = kOverCv;
      }

      break;
    }
  }
}

void CFEM_DG_NSSolver::ViscousNormalFluxIntegrationPoint_2D(const su2double *sol,
                                                            const su2double solGradCart[4][2],
                                                            const su2double *normal,
                                                            const su2double HeatFlux,
                                                            const su2double factHeatFlux_Lam,
                                                            const su2double factHeatFlux_Turb,
                                                            const su2double wallDist,
                                                            const su2double lenScale_LES,
                                                                  su2double &Viscosity,
                                                                  su2double &kOverCv,
                                                                  su2double *normalFlux) {

  /*--- Compute the velocities and static energy in this integration point. ---*/
  const su2double rhoInv = 1.0/sol[0];
  const su2double u = rhoInv*sol[1];
  const su2double v = rhoInv*sol[2];

  const su2double TotalEnergy  = rhoInv*sol[3];
  const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

  /*--- Compute the Cartesian gradients of the velocities and static energy
        in this integration point and also the divergence of the velocity. ---*/
  const su2double dudx = rhoInv*(solGradCart[1][0] - u*solGradCart[0][0]);
  const su2double dudy = rhoInv*(solGradCart[1][1] - u*solGradCart[0][1]);

  const su2double dvdx = rhoInv*(solGradCart[2][0] - v*solGradCart[0][0]);
  const su2double dvdy = rhoInv*(solGradCart[2][1] - v*solGradCart[0][1]);

  const su2double dStaticEnergyDx = rhoInv*(solGradCart[3][0]
                                  -         TotalEnergy*solGradCart[0][0])
                                  - u*dudx - v*dvdx;
  const su2double dStaticEnergyDy = rhoInv*(solGradCart[3][1]
                                  -         TotalEnergy*solGradCart[0][1])
                                  - u*dudy - v*dvdy;

  const su2double divVel = dudx + dvdy;

  /*--- Compute the laminar viscosity. ---*/
  FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
  const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

  /*--- Compute the eddy viscosity, if needed. ---*/
  su2double ViscosityTurb = 0.0;
  if( SGSModelUsed ) {
    su2double velGrad[3][3];
    velGrad[0][0] = dudx; velGrad[0][1] = dudy; velGrad[0][2] = 0.0;
    velGrad[1][0] = dvdx; velGrad[1][1] = dvdy; velGrad[1][2] = 0.0;
    velGrad[2][0] = 0.0;  velGrad[2][1] = 0.0;  velGrad[2][2] = 0.0;

    ViscosityTurb = SGSModel->ComputeEddyViscosity(nDim, sol[0], velGrad,
                                                   lenScale_LES, wallDist);
  }

  /* Compute the total viscosity and heat conductivity. Note that the heat
     conductivity is divided by the Cv, because gradients of internal energy
     are computed and not temperature. */
  Viscosity = ViscosityLam + ViscosityTurb;
  kOverCv   = ViscosityLam*factHeatFlux_Lam + ViscosityTurb*factHeatFlux_Turb;

  /*--- Set the value of the second viscosity and compute the divergence
        term in the viscous normal stresses. ---*/
  const su2double lambda     = -TWO3*Viscosity;
  const su2double lamDivTerm =  lambda*divVel;

  /*--- Compute the viscous stress tensor and minus the heatflux vector. ---*/
  const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
  const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
  const su2double tauxy = Viscosity*(dudy + dvdx);

  const su2double qx = kOverCv*dStaticEnergyDx;
  const su2double qy = kOverCv*dStaticEnergyDy;

  /* Compute the unscaled normal vector. */
  const su2double nx = normal[0]*normal[2];
  const su2double ny = normal[1]*normal[2];

  /*--- Compute the viscous normal flux. Note that the energy flux get a 
        contribution from both the prescribed and the computed heat flux.
        At least one of these terms is zero. ---*/
  normalFlux[0] = 0.0;
  normalFlux[1] = tauxx*nx + tauxy*ny;
  normalFlux[2] = tauxy*nx + tauyy*ny;
  normalFlux[3] = normal[2]*HeatFlux
                + (u*tauxx + v*tauxy + qx)*nx + (u*tauxy + v*tauyy + qy)*ny;
}

void CFEM_DG_NSSolver::ViscousNormalFluxIntegrationPoint_3D(const su2double *sol,
                                                            const su2double solGradCart[5][3],
                                                            const su2double *normal,
                                                            const su2double HeatFlux,
                                                            const su2double factHeatFlux_Lam,
                                                            const su2double factHeatFlux_Turb,
                                                            const su2double wallDist,
                                                            const su2double lenScale_LES,
                                                                  su2double &Viscosity,
                                                                  su2double &kOverCv,
                                                                  su2double *normalFlux) {

  /*--- Compute the velocities and static energy in this integration point. ---*/
  const su2double rhoInv = 1.0/sol[0];
  const su2double u = rhoInv*sol[1];
  const su2double v = rhoInv*sol[2];
  const su2double w = rhoInv*sol[3];

  const su2double TotalEnergy  = rhoInv*sol[4];
  const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v + w*w);

  /*--- Compute the Cartesian gradients of the velocities and static energy
        in this integration point and also the divergence of the velocity. ---*/
  const su2double dudx = rhoInv*(solGradCart[1][0] - u*solGradCart[0][0]);
  const su2double dudy = rhoInv*(solGradCart[1][1] - u*solGradCart[0][1]);
  const su2double dudz = rhoInv*(solGradCart[1][2] - u*solGradCart[0][2]);

  const su2double dvdx = rhoInv*(solGradCart[2][0] - v*solGradCart[0][0]);
  const su2double dvdy = rhoInv*(solGradCart[2][1] - v*solGradCart[0][1]);
  const su2double dvdz = rhoInv*(solGradCart[2][2] - v*solGradCart[0][2]);

  const su2double dwdx = rhoInv*(solGradCart[3][0] - w*solGradCart[0][0]);
  const su2double dwdy = rhoInv*(solGradCart[3][1] - w*solGradCart[0][1]);
  const su2double dwdz = rhoInv*(solGradCart[3][2] - w*solGradCart[0][2]);

  const su2double dStaticEnergyDx = rhoInv*(solGradCart[4][0]
                                  -         TotalEnergy*solGradCart[0][0])
                                  - u*dudx - v*dvdx - w*dwdx;
  const su2double dStaticEnergyDy = rhoInv*(solGradCart[4][1] 
                                  -         TotalEnergy*solGradCart[0][1])
                                  - u*dudy - v*dvdy - w*dwdy;
  const su2double dStaticEnergyDz = rhoInv*(solGradCart[4][2] 
                                  -         TotalEnergy*solGradCart[0][2])
                                  - u*dudz - v*dvdz - w*dwdz;

  const su2double divVel = dudx + dvdy + dwdz;

  /*--- Compute the laminar viscosity. ---*/
  FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
  const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

  /*--- Compute the eddy viscosity, if needed. ---*/
  su2double ViscosityTurb = 0.0;
  if( SGSModelUsed ) {
    su2double velGrad[3][3];
    velGrad[0][0] = dudx; velGrad[0][1] = dudy; velGrad[0][2] = dudz;
    velGrad[1][0] = dvdx; velGrad[1][1] = dvdy; velGrad[1][2] = dvdz;
    velGrad[2][0] = dwdx; velGrad[2][1] = dwdy; velGrad[2][2] = dwdz;

    ViscosityTurb = SGSModel->ComputeEddyViscosity(nDim, sol[0], velGrad,
                                                   lenScale_LES, wallDist);
  }

  /* Compute the total viscosity and heat conductivity. Note that the heat
     conductivity is divided by the Cv, because gradients of internal energy
     are computed and not temperature. */
  Viscosity = ViscosityLam + ViscosityTurb;
  kOverCv   = ViscosityLam*factHeatFlux_Lam + ViscosityTurb*factHeatFlux_Turb;

  /*--- Set the value of the second viscosity and compute the divergence
        term in the viscous normal stresses. ---*/
  const su2double lambda     = -TWO3*Viscosity;
  const su2double lamDivTerm =  lambda*divVel;

  /*--- Compute the viscous stress tensor and minus the heatflux vector. ---*/
  const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
  const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
  const su2double tauzz = 2.0*Viscosity*dwdz + lamDivTerm;

  const su2double tauxy = Viscosity*(dudy + dvdx);
  const su2double tauxz = Viscosity*(dudz + dwdx);
  const su2double tauyz = Viscosity*(dvdz + dwdy);

  const su2double qx = kOverCv*dStaticEnergyDx;
  const su2double qy = kOverCv*dStaticEnergyDy;
  const su2double qz = kOverCv*dStaticEnergyDz;

  /* Compute the unscaled normal vector. */
  const su2double nx = normal[0]*normal[3];
  const su2double ny = normal[1]*normal[3];
  const su2double nz = normal[2]*normal[3];

  /*--- Compute the viscous normal flux. Note that the energy flux get a 
        contribution from both the prescribed and the computed heat flux.
        At least one of these terms is zero. ---*/
  normalFlux[0] = 0.0;
  normalFlux[1] = tauxx*nx + tauxy*ny + tauxz*nz;
  normalFlux[2] = tauxy*nx + tauyy*ny + tauyz*nz;
  normalFlux[3] = tauxz*nx + tauyz*ny + tauzz*nz;
  normalFlux[4] = normal[3]*HeatFlux
                + (u*tauxx + v*tauxy + w*tauxz + qx)*nx
                + (u*tauxy + v*tauyy + w*tauyz + qy)*ny
                + (u*tauxz + v*tauyz + w*tauzz + qz)*nz;
}

void CFEM_DG_NSSolver::PenaltyTermsFluxFace(const unsigned short nInt,
                                            const su2double      *solInt0,
                                            const su2double      *solInt1,
                                            const su2double      *viscosityInt0,
                                            const su2double      *viscosityInt1,
                                            const su2double      *kOverCvInt0,
                                            const su2double      *kOverCvInt1,
                                            const su2double      ConstPenFace,
                                            const su2double      lenScale0,
                                            const su2double      lenScale1,
                                            const su2double      *metricNormalsFace,
                                                  su2double      *penaltyFluxes) {

  /* Constant ratio of the second viscosity and the viscosity itself. */
  const su2double lambdaOverMu = -TWO3;

  /* The eigenvalues of the viscous Jacobian, scaled by the kinematic viscosity,
     are 1.0, 2.0 + lambdaOverMu and kOverCv/Mu. The last is variable due to the
     possible presence of an eddy viscosity, but the first two are constant and
     the maximum can be determined. */
  const su2double radOverNuTerm = max(1.0, 2.0+lambdaOverMu);

  /*--- Loop over the integration points to compute the penalty fluxes. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the variables for this integration point. */
    const su2double *sol0   = solInt0 + nVar*i;
    const su2double *sol1   = solInt1 + nVar*i;
    const su2double *normal = metricNormalsFace + i*(nDim+1);
    su2double       *flux   = penaltyFluxes + nVar*i;

    /* Determine the ratio of kOverCv and mu for both sides and compute the
       spectral radius of the viscous terms, scaled by kinematic viscosity. */
    const su2double factHeatFlux0 = kOverCvInt0[i]/viscosityInt0[i];
    const su2double factHeatFlux1 = kOverCvInt1[i]/viscosityInt1[i];

    const su2double radOverNu0 = max(radOverNuTerm, factHeatFlux0);
    const su2double radOverNu1 = max(radOverNuTerm, factHeatFlux1);

    /* Compute the kinematic viscosities of both sides. Multiply it by
       ConstPenFace and divide by the length scale. */
    const su2double nu0 = ConstPenFace*viscosityInt0[i]/(lenScale0*sol0[0]);
    const su2double nu1 = ConstPenFace*viscosityInt1[i]/(lenScale1*sol1[0]);

    /* Compute the penalty parameter of this face as the maximum of the
       penalty parameter from both sides. Multiply by the area to obtain the
       correct expression. */
    const su2double pen0 = radOverNu0*nu0;
    const su2double pen1 = radOverNu1*nu1;

    const su2double penFace = normal[nDim]*max(pen0, pen1);

    /* Compute the penalty flux, where it is assumed that the normal points from
       side 0 to side 1. */
    flux[0] = 0;
    for(unsigned short j=0; j<nVar; ++j)
      flux[j] = penFace*(sol0[j] - sol1[j]);
  }
}

void CFEM_DG_NSSolver::SymmetrizingFluxesFace(const unsigned short nInt,
                                              const su2double      *solInt0,
                                              const su2double      *solInt1,
                                              const su2double      *viscosityInt0,
                                              const su2double      *viscosityInt1,
                                              const su2double      *kOverCvInt0,
                                              const su2double      *kOverCvInt1,
                                              const su2double      *metricNormalsFace,
                                                    su2double      *symmFluxes) {

  /* Constant ratio of the second viscosity and the viscosity itself. */
  const su2double lambdaOverMu = -TWO3;

  /*--- Set two factors such that either the original or the transposed diffusion
        tensor is taken in the symmetrizing fluxes. ---*/

  /* Use the following line for the original formulation. */
  const su2double alpha = lambdaOverMu;

  /* Use the following line for the transposed formulation. */
  //const su2double alpha = 1.0;

  /* Other constants, which appear in the symmetrizing fluxes. */
  const su2double beta     = lambdaOverMu + 1.0 - alpha;
  const su2double alphaP1  = alpha + 1.0;
  const su2double lambdaP1 = lambdaOverMu + 1.0;

  /* Make a distinction between two and three space dimensions
     in order to have the most efficient code. */
  switch( ctc::nDim ) {

    case 2: {

      /*--- 2D simulation. Loop over the number of integration points
            of the face. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const su2double *sol0   = solInt0 + ctc::nVar*i;
        const su2double *sol1   = solInt1 + ctc::nVar*i;
        const su2double *normal = metricNormalsFace + i*(ctc::nDim+1);
        su2double       *flux   = symmFluxes + i*ctc::nDim*ctc::nVar;

        /* Determine the difference in conservative variables. Multiply these
           differences by the length of the normal vector to obtain the correct
           dimensions for the symmetrizing fluxes. */
        const su2double dSol[] = {normal[2]*(sol0[0] - sol1[0]),
                                  normal[2]*(sol0[1] - sol1[1]),
                                  normal[2]*(sol0[2] - sol1[2]),
                                  normal[2]*(sol0[3] - sol1[3])};

        /*--- Compute the terms that occur in the symmetrizing fluxes
              for state 0 and state 1. ---*/
        const su2double DensityInv0 = 1.0/sol0[0],           DensityInv1 = 1.0/sol1[0];
        const su2double Etot0 = DensityInv0*sol0[3],         Etot1 = DensityInv1*sol1[3];
        const su2double nu0 = DensityInv0*viscosityInt0[i],  nu1 = DensityInv1*viscosityInt1[i];
        const su2double kScal0 = DensityInv0*kOverCvInt0[i], kScal1 = DensityInv1*kOverCvInt1[i];

        const su2double vel0[] = {DensityInv0*sol0[1], DensityInv0*sol0[2]};
        const su2double vel1[] = {DensityInv1*sol1[1], DensityInv1*sol1[2]};

        const su2double velNorm0 = vel0[0]*normal[0] + vel0[1]*normal[1];
        const su2double velNorm1 = vel1[0]*normal[0] + vel1[1]*normal[1];

        const su2double velSquared0 = vel0[0]*vel0[0] + vel0[1]*vel0[1];
        const su2double velSquared1 = vel1[0]*vel1[0] + vel1[1]*vel1[1];

        /*--- Compute the average of the terms that occur in the symmetrizing
              fluxes. The average of the left and right terms is taken, rather
              than the terms evaluated at the average state, because the viscous
              fluxes are also computed as the average of the fluxes and not the
              fluxes of the averaged state. ---*/
        const su2double nuAvg    = 0.5*(nu0    + nu1);
        const su2double kScalAvg = 0.5*(kScal0 + kScal1);
        const su2double nuVelSquaredAvg     = 0.5*(nu0*velSquared0 + nu1*velSquared1);
        const su2double nuVelNormAve        = 0.5*(nu0*velNorm0    + nu1*velNorm1);
        const su2double kScalEminVelSquaredAve = 0.5*(kScal0*(Etot0-velSquared0)
                                               +      kScal1*(Etot1-velSquared1));

        const su2double nuVelAvg[] = {0.5*(nu0*vel0[0] + nu1*vel1[0]),
                                      0.5*(nu0*vel0[1] + nu1*vel1[1])};
        const su2double kScalVelAvg[] = {0.5*(kScal0*vel0[0] + kScal1*vel1[0]),
                                         0.5*(kScal0*vel0[1] + kScal1*vel1[1])};
        const su2double nuVelVelAvg[] = {0.5*(nu0*vel0[0]*velNorm0 + nu1*vel1[0]*velNorm1),
                                         0.5*(nu0*vel0[1]*velNorm0 + nu1*vel1[1]*velNorm1)};

        /*--- Abbreviations to make the flux computations a bit more efficient. ---*/
        const su2double abv1 = normal[0]  *dSol[1] + normal[1]  *dSol[2];
        const su2double abv2 = nuVelAvg[0]*dSol[1] + nuVelAvg[1]*dSol[2];

        const su2double abv2kScal = kScalVelAvg[0]*dSol[1] + kScalVelAvg[1]*dSol[2];

        const su2double abv3 = beta*(nuAvg*abv1 - nuVelNormAve*dSol[0]);
        const su2double abv4 = kScalAvg*dSol[3] - abv2kScal
                             - kScalEminVelSquaredAve*dSol[0] + abv2;

        /*--- Compute the symmetrizing fluxes. ---*/
        flux[0] = 0.0;
        flux[1] = abv3 + alphaP1*normal[0]*(nuAvg*dSol[1] - nuVelAvg[0]*dSol[0]);
        flux[2] = nuAvg*(normal[0]*dSol[2] + alpha*normal[1]*dSol[1])
                - (normal[0]*nuVelAvg[1] + alpha*normal[1]*nuVelAvg[0])*dSol[0];
        flux[3] = normal[0]*abv4
                - (lambdaP1*nuVelVelAvg[0] + nuVelSquaredAvg*normal[0])*dSol[0]
                + alpha*nuVelNormAve*dSol[1] + beta*nuVelAvg[0]*abv1;

        flux[4] = 0.0;
        flux[5] = nuAvg*(normal[1]*dSol[1] + alpha*normal[0]*dSol[2])
                - (normal[1]*nuVelAvg[0] + alpha*normal[0]*nuVelAvg[1])*dSol[0];
        flux[6] = abv3 + alphaP1*normal[1]*(nuAvg*dSol[2] - nuVelAvg[1]*dSol[0]);
        flux[7] = normal[1]*abv4
                - (lambdaP1*nuVelVelAvg[1] + nuVelSquaredAvg*normal[1])*dSol[0]
                + alpha*nuVelNormAve*dSol[2] + beta*nuVelAvg[1]*abv1;
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case 3: {

      /*--- 3D simulation. Loop over the number of integration points
            of the face. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const su2double *sol0   = solInt0 + ctc::nVar*i;
        const su2double *sol1   = solInt1 + ctc::nVar*i;
        const su2double *normal = metricNormalsFace + i*(ctc::nDim+1);
        su2double       *flux   = symmFluxes + i*ctc::nDim*ctc::nVar;

        /* Determine the difference in conservative variables. Multiply these
           differences by the length of the normal vector to obtain the correct
           dimensions for the symmetrizing fluxes. */
        const su2double dSol[] = {normal[3]*(sol0[0] - sol1[0]),
                                  normal[3]*(sol0[1] - sol1[1]),
                                  normal[3]*(sol0[2] - sol1[2]),
                                  normal[3]*(sol0[3] - sol1[3]),
                                  normal[3]*(sol0[4] - sol1[4])};

        /*--- Compute the terms that occur in the symmetrizing fluxes
              for state 0 and state 1. ---*/
        const su2double DensityInv0 = 1.0/sol0[0],           DensityInv1 = 1.0/sol1[0];
        const su2double Etot0 = DensityInv0*sol0[4],         Etot1 = DensityInv1*sol1[4];
        const su2double nu0 = DensityInv0*viscosityInt0[i],  nu1 = DensityInv1*viscosityInt1[i];
        const su2double kScal0 = DensityInv0*kOverCvInt0[i], kScal1 = DensityInv1*kOverCvInt1[i];

        const su2double vel0[] = {DensityInv0*sol0[1], DensityInv0*sol0[2], DensityInv0*sol0[3]};
        const su2double vel1[] = {DensityInv1*sol1[1], DensityInv1*sol1[2], DensityInv1*sol1[3]};

        const su2double velNorm0 = vel0[0]*normal[0] + vel0[1]*normal[1] + vel0[2]*normal[2];
        const su2double velNorm1 = vel1[0]*normal[0] + vel1[1]*normal[1] + vel1[2]*normal[2];

        const su2double velSquared0 = vel0[0]*vel0[0] + vel0[1]*vel0[1] + vel0[2]*vel0[2];
        const su2double velSquared1 = vel1[0]*vel1[0] + vel1[1]*vel1[1] + vel1[2]*vel1[2];

        /*--- Compute the average of the terms that occur in the symmetrizing
              fluxes. The average of the left and right terms is taken, rather
              than the terms evaluated at the average state, because the viscous
              fluxes are also computed as the average of the fluxes and not the
              fluxes of the averaged state. ---*/
        const su2double nuAvg    = 0.5*(nu0    + nu1);
        const su2double kScalAvg = 0.5*(kScal0 + kScal1);
        const su2double nuVelSquaredAvg     = 0.5*(nu0*velSquared0 + nu1*velSquared1);
        const su2double nuVelNormAve        = 0.5*(nu0*velNorm0    + nu1*velNorm1);
        const su2double kScalEminVelSquaredAve = 0.5*(kScal0*(Etot0-velSquared0)
                                               +      kScal1*(Etot1-velSquared1));

        const su2double nuVelAvg[] = {0.5*(nu0*vel0[0] + nu1*vel1[0]),
                                      0.5*(nu0*vel0[1] + nu1*vel1[1]),
                                      0.5*(nu0*vel0[2] + nu1*vel1[2])};
        const su2double kScalVelAvg[] = {0.5*(kScal0*vel0[0] + kScal1*vel1[0]),
                                         0.5*(kScal0*vel0[1] + kScal1*vel1[1]),
                                         0.5*(kScal0*vel0[2] + kScal1*vel1[2])};
        const su2double nuVelVelAvg[] = {0.5*(nu0*vel0[0]*velNorm0 + nu1*vel1[0]*velNorm1),
                                         0.5*(nu0*vel0[1]*velNorm0 + nu1*vel1[1]*velNorm1),
                                         0.5*(nu0*vel0[2]*velNorm0 + nu1*vel1[2]*velNorm1)};

        /*--- Abbreviations to make the flux computations a bit more efficient. ---*/
        const su2double abv1 = normal[0]  *dSol[1] + normal[1]  *dSol[2] + normal[2]  *dSol[3];
        const su2double abv2 = nuVelAvg[0]*dSol[1] + nuVelAvg[1]*dSol[2] + nuVelAvg[2]*dSol[3];

        const su2double abv2kScal = kScalVelAvg[0]*dSol[1] + kScalVelAvg[1]*dSol[2]
                                  + kScalVelAvg[2]*dSol[3];

        const su2double abv3 = beta*(nuAvg*abv1 - nuVelNormAve*dSol[0]);
        const su2double abv4 = kScalAvg*dSol[4] - abv2kScal
                             - kScalEminVelSquaredAve*dSol[0] + abv2;

        /*--- Compute the symmetrizing fluxes. ---*/
        flux[0] = 0.0;
        flux[1] = abv3 + alphaP1*normal[0]*(nuAvg*dSol[1] - nuVelAvg[0]*dSol[0]);
        flux[2] = nuAvg*(normal[0]*dSol[2] + alpha*normal[1]*dSol[1])
                - (normal[0]*nuVelAvg[1] + alpha*normal[1]*nuVelAvg[0])*dSol[0];
        flux[3] = nuAvg*(normal[0]*dSol[3] + alpha*normal[2]*dSol[1])
                - (normal[0]*nuVelAvg[2] + alpha*normal[2]*nuVelAvg[0])*dSol[0];
        flux[4] = normal[0]*abv4
                - (lambdaP1*nuVelVelAvg[0] + nuVelSquaredAvg*normal[0])*dSol[0]
                + alpha*nuVelNormAve*dSol[1] + beta*nuVelAvg[0]*abv1;

        flux[5] = 0.0;
        flux[6] = nuAvg*(normal[1]*dSol[1] + alpha*normal[0]*dSol[2])
                - (normal[1]*nuVelAvg[0] + alpha*normal[0]*nuVelAvg[1])*dSol[0];
        flux[7] = abv3 + alphaP1*normal[1]*(nuAvg*dSol[2] - nuVelAvg[1]*dSol[0]);
        flux[8] = nuAvg*(normal[1]*dSol[3] + alpha*normal[2]*dSol[2])
                - (normal[1]*nuVelAvg[2] + alpha*normal[2]*nuVelAvg[1])*dSol[0];
        flux[9] = normal[1]*abv4
                - (lambdaP1*nuVelVelAvg[1] + nuVelSquaredAvg*normal[1])*dSol[0]
                + alpha*nuVelNormAve*dSol[2] + beta*nuVelAvg[1]*abv1;

        flux[10] = 0.0;
        flux[11] = nuAvg*(normal[2]*dSol[1] + alpha*normal[0]*dSol[3])
                 - (normal[2]*nuVelAvg[0] + alpha*normal[0]*nuVelAvg[2])*dSol[0];
        flux[12] = nuAvg*(normal[2]*dSol[2] + alpha*normal[1]*dSol[3])
                 - (normal[2]*nuVelAvg[1] + alpha*normal[1]*nuVelAvg[2])*dSol[0];
        flux[13] = abv3 + alphaP1*normal[2]*(nuAvg*dSol[3] - nuVelAvg[2]*dSol[0]);
        flux[14] = normal[2]*abv4
                 - (lambdaP1*nuVelVelAvg[2] + nuVelSquaredAvg*normal[2])*dSol[0]
                 + alpha*nuVelNormAve*dSol[3] + beta*nuVelAvg[2]*abv1;
      }

      break;
    }
  }
}

void CFEM_DG_NSSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

  su2double *solIntL      = VecTmpMemory.data();
  su2double *solIntR      = solIntL      + nIntegrationMax*nVar;
  su2double *viscosityInt = solIntR      + nIntegrationMax*nVar;
  su2double *kOverCvInt   = viscosityInt + nIntegrationMax;
  su2double *gradSolInt   = kOverCvInt   + nIntegrationMax;
  su2double *fluxes       = gradSolInt   + sizeGradSolInt;
  su2double *viscFluxes   = fluxes       + sizeFluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /*--- Determine the number of integration points and set the right state
          to the free stream value. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {
      su2double *UR = solIntR + i*nVar;
      for(unsigned short j=0; j<nVar; ++j)
        UR[j] = ConsVarFreeStream[j];
    }

    /* Compute the length scale for the LES of the adjacent element. */
    const unsigned long  elemID = surfElem[l].volElemID;
    const unsigned short iind   = volElem[elemID].indStandardElement;

    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, nInt, nDOFsElem, 0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(), lenScale_LES,
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Sym_Plane(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *conv_numerics,
                                    CNumerics *visc_numerics,
                                    CConfig *config,
                                    unsigned short val_marker) {
  su2double tick = 0.0;

  /* Constant factor present in the heat flux vector, namely the ratio of
     thermal conductivity and viscosity. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Set the pointers for the local arrays. ---*/
  unsigned int sizeFluxes = nIntegrationMax*ctc::nDim;
  sizeFluxes = ctc::nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*ctc::nDim*max(nVar,nDOFsMax);

  su2double *solIntL      = VecTmpMemory.data();
  su2double *solIntR      = solIntL      + nIntegrationMax*nVar;
  su2double *viscosityInt = solIntR      + nIntegrationMax*nVar;
  su2double *kOverCvInt   = viscosityInt + nIntegrationMax;
  su2double *gradSolInt   = kOverCvInt   + nIntegrationMax;
  su2double *fluxes       = gradSolInt   + sizeGradSolInt;
  su2double *viscFluxes   = fluxes       + sizeFluxes;

  /* Set the pointer solElem to fluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the element and the fluxes to be computed. */
  su2double *solElem = fluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Get the required information from the standard element. */
    const unsigned short ind          = surfElem[l].indStandardElement;
    const unsigned short nInt         = standardBoundaryFacesSol[ind].GetNIntegration();
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    /* Compute the length scale for the LES of the adjacent element. */
    const unsigned long  elemID = surfElem[l].volElemID;
    const unsigned short iind   = volElem[elemID].indStandardElement;

    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

    /* Easier storage of the wall distance array for this surface element. */
    const su2double *wallDistance = surfElem[l].wallDistance.data();

    /* Determine the offset between r- and -s-derivatives, which is also the
       offset between s- and t-derivatives. */
    const unsigned short offDeriv = ctc::nVar*nInt;

    /* Store the solution of the DOFs of the adjacent element in contiguous
       memory such that the function DenseMatrixProduct can be used to compute
       the gradients solution variables in the integration points of the face. */
    for(unsigned short i=0; i<nDOFsElem; ++i) {
      const su2double *solDOF = VecSolDOFs.data() + ctc::nVar*surfElem[l].DOFsSolElement[i];
      su2double       *sol    = solElem + ctc::nVar*i;
      for(unsigned short j=0; j<ctc::nVar; ++j)
        sol[j] = solDOF[j];
    }

    /* Compute the left gradients in the integration points. Call the general
       function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nInt*nDim, nVar, nDOFsElem, derBasisElem, solElem, gradSolInt);
    config->GEMM_Tock(tick, "BC_Sym_Plane", nInt*nDim, nVar, nDOFsElem);

    /* Make a distinction between two and three space dimensions
       in order to have the most efficient code. */
    switch( ctc::nDim ) {

      case 2: {

        /*--- 2D simulation. Loop over the integration points to apply the
              symmetry condition and to compute the viscous normal fluxes. This
              is done in the same loop, because the gradients in the right
              integration points are constructed using the symmetry boundary
              condition. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution and the normals
             for this integration point. */
          const su2double *UL      = solIntL + i*ctc::nVar;
                su2double *UR      = solIntR + i*ctc::nVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(ctc::nDim+1);

          /* Compute twice the normal component of the momentum variables. The
             factor 2 comes from the fact that the velocity must be mirrored. */
          const su2double rVn = 2.0*(UL[1]*normals[0] + UL[2]*normals[1]);

          /* Set the right state. Note that the total energy of the right state is
             identical to the left state, because the magnitude of the velocity
             remains the same. */
          UR[0] = UL[0];
          UR[1] = UL[1] - rVn*normals[0];
          UR[2] = UL[2] - rVn*normals[1];
          UR[3] = UL[3];

          /* Easier storage of the metric terms and left gradients of this
             integration point. */
          const su2double *metricTerms = surfElem[l].metricCoorDerivFace.data()
                                       + i*ctc::nDim*ctc::nDim;
          const su2double *dULDr       = gradSolInt + ctc::nVar*i;
          const su2double *dULDs       = dULDr + offDeriv;

          /* Easier storage of the metric terms in this integration point. */
          const su2double drdx = metricTerms[0];
          const su2double drdy = metricTerms[1];
          const su2double dsdx = metricTerms[2];
          const su2double dsdy = metricTerms[3];

          /* Compute the Cartesian gradients of the left solution. */
          su2double ULGradCart[4][2];
          ULGradCart[0][0] = dULDr[0]*drdx + dULDs[0]*dsdx;
          ULGradCart[1][0] = dULDr[1]*drdx + dULDs[1]*dsdx;
          ULGradCart[2][0] = dULDr[2]*drdx + dULDs[2]*dsdx;
          ULGradCart[3][0] = dULDr[3]*drdx + dULDs[3]*dsdx;

          ULGradCart[0][1] = dULDr[0]*drdy + dULDs[0]*dsdy;
          ULGradCart[1][1] = dULDr[1]*drdy + dULDs[1]*dsdy;
          ULGradCart[2][1] = dULDr[2]*drdy + dULDs[2]*dsdy;
          ULGradCart[3][1] = dULDr[3]*drdy + dULDs[3]*dsdy;

          /* Determine the normal gradients of all variables. */
          su2double ULGradNorm[4];
          ULGradNorm[0] = ULGradCart[0][0]*normals[0] + ULGradCart[0][1]*normals[1];
          ULGradNorm[1] = ULGradCart[1][0]*normals[0] + ULGradCart[1][1]*normals[1];
          ULGradNorm[2] = ULGradCart[2][0]*normals[0] + ULGradCart[2][1]*normals[1];
          ULGradNorm[3] = ULGradCart[3][0]*normals[0] + ULGradCart[3][1]*normals[1];

          /* For the construction of the gradients of the right solution, also
             the Cartesian gradients and the normal gradient of the normal
             momentum is needed. This is computed below. */
          su2double GradCartNormMomL[2];
          GradCartNormMomL[0] = ULGradCart[1][0]*normals[0] + ULGradCart[2][0]*normals[1];
          GradCartNormMomL[1] = ULGradCart[1][1]*normals[0] + ULGradCart[2][1]*normals[1];

          const su2double GradNormNormMomL = ULGradNorm[1]*normals[0] + ULGradNorm[2]*normals[1];

          /* Abbreviate twice the normal vector. */
          const su2double tnx = 2.0*normals[0], tny = 2.0*normals[1];

          /*--- Construct the gradients of the right solution. The tangential
                gradients of the normal momentum and normal gradients of the other
                variables, density, energy and tangential momentum, must be negated.
                For the common situation that the symmetry plane coincides with an
                x- or y-plane this boils down to a multiplication by -1 or +1,
                depending on the variable and gradient. However, the implementation
                below is valid for an arbitrary orientation of the symmetry plane. ---*/
          su2double URGradCart[4][2];
          URGradCart[0][0] = ULGradCart[0][0] - tnx*ULGradNorm[0];
          URGradCart[1][0] = ULGradCart[1][0] - tnx*ULGradNorm[1] - tnx*GradCartNormMomL[0] + tnx*tnx*GradNormNormMomL;
          URGradCart[2][0] = ULGradCart[2][0] - tnx*ULGradNorm[2] - tny*GradCartNormMomL[0] + tnx*tny*GradNormNormMomL;
          URGradCart[3][0] = ULGradCart[3][0] - tnx*ULGradNorm[3];

          URGradCart[0][1] = ULGradCart[0][1] - tny*ULGradNorm[0];
          URGradCart[1][1] = ULGradCart[1][1] - tny*ULGradNorm[1] - tnx*GradCartNormMomL[1] + tny*tnx*GradNormNormMomL;
          URGradCart[2][1] = ULGradCart[2][1] - tny*ULGradNorm[2] - tny*GradCartNormMomL[1] + tny*tny*GradNormNormMomL;
          URGradCart[3][1] = ULGradCart[3][1] - tny*ULGradNorm[3];

          /*--- Compute the viscous fluxes of the left and right state and average
                them to compute the correct viscous flux for a symmetry boundary. ---*/
          const su2double wallDist = wallDistance ? wallDistance[i] : 0.0;
          su2double viscFluxL[4], viscFluxR[4];
          su2double Viscosity, kOverCv;

          ViscousNormalFluxIntegrationPoint_2D(UR, URGradCart, normals, 0.0,
                                               factHeatFlux_Lam, factHeatFlux_Turb,
                                               wallDist, lenScale_LES, Viscosity,
                                               kOverCv, viscFluxR);
          ViscousNormalFluxIntegrationPoint_2D(UL, ULGradCart, normals, 0.0,
                                               factHeatFlux_Lam, factHeatFlux_Turb,
                                               wallDist, lenScale_LES, Viscosity,
                                               kOverCv, viscFluxL);
          viscosityInt[i] = Viscosity;
          kOverCvInt[i]   = kOverCv;

          su2double *viscNormalFlux = viscFluxes + i*ctc::nVar;
          viscNormalFlux[0] = 0.5*(viscFluxL[0] + viscFluxR[0]);
          viscNormalFlux[1] = 0.5*(viscFluxL[1] + viscFluxR[1]);
          viscNormalFlux[2] = 0.5*(viscFluxL[2] + viscFluxR[2]);
          viscNormalFlux[3] = 0.5*(viscFluxL[3] + viscFluxR[3]);
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /*--- 3D simulation. Loop over the integration points to apply the
              symmetry condition and to compute the viscous normal fluxes. This
              is done in the same loop, because the gradients in the right
              integration points are constructed using the symmetry boundary
              condition. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution and the normals
             for this integration point. */
          const su2double *UL      = solIntL + i*ctc::nVar;
                su2double *UR      = solIntR + i*ctc::nVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(ctc::nDim+1);

          /* Compute twice the normal component of the momentum variables. The
             factor 2 comes from the fact that the velocity must be mirrored. */
          const su2double rVn = 2.0*(UL[1]*normals[0] + UL[2]*normals[1] + UL[3]*normals[2]);

          /* Set the right state. Note that the total energy of the right state is
             identical to the left state, because the magnitude of the velocity
             remains the same. */
          UR[0] = UL[0];
          UR[1] = UL[1] - rVn*normals[0];
          UR[2] = UL[2] - rVn*normals[1];
          UR[3] = UL[3] - rVn*normals[2];
          UR[4] = UL[4];

          /* Easier storage of the metric terms and left gradients of this
             integration point. */
          const su2double *metricTerms = surfElem[l].metricCoorDerivFace.data()
                                       + i*ctc::nDim*ctc::nDim;
          const su2double *dULDr       = gradSolInt + ctc::nVar*i;
          const su2double *dULDs       = dULDr + offDeriv;
          const su2double *dULDt       = dULDs + offDeriv;

          /* Easier storage of the metric terms in this integration point. */
          const su2double drdx = metricTerms[0];
          const su2double drdy = metricTerms[1];
          const su2double drdz = metricTerms[2];

          const su2double dsdx = metricTerms[3];
          const su2double dsdy = metricTerms[4];
          const su2double dsdz = metricTerms[5];

          const su2double dtdx = metricTerms[6];
          const su2double dtdy = metricTerms[7];
          const su2double dtdz = metricTerms[8];

          /* Compute the Cartesian gradients of the left solution. */
          su2double ULGradCart[5][3];
          ULGradCart[0][0] = dULDr[0]*drdx + dULDs[0]*dsdx + dULDt[0]*dtdx;
          ULGradCart[1][0] = dULDr[1]*drdx + dULDs[1]*dsdx + dULDt[1]*dtdx;
          ULGradCart[2][0] = dULDr[2]*drdx + dULDs[2]*dsdx + dULDt[2]*dtdx;
          ULGradCart[3][0] = dULDr[3]*drdx + dULDs[3]*dsdx + dULDt[3]*dtdx;
          ULGradCart[4][0] = dULDr[4]*drdx + dULDs[4]*dsdx + dULDt[4]*dtdx;

          ULGradCart[0][1] = dULDr[0]*drdy + dULDs[0]*dsdy + dULDt[0]*dtdy;
          ULGradCart[1][1] = dULDr[1]*drdy + dULDs[1]*dsdy + dULDt[1]*dtdy;
          ULGradCart[2][1] = dULDr[2]*drdy + dULDs[2]*dsdy + dULDt[2]*dtdy;
          ULGradCart[3][1] = dULDr[3]*drdy + dULDs[3]*dsdy + dULDt[3]*dtdy;
          ULGradCart[4][1] = dULDr[4]*drdy + dULDs[4]*dsdy + dULDt[4]*dtdy;

          ULGradCart[0][2] = dULDr[0]*drdz + dULDs[0]*dsdz + dULDt[0]*dtdz;
          ULGradCart[1][2] = dULDr[1]*drdz + dULDs[1]*dsdz + dULDt[1]*dtdz;
          ULGradCart[2][2] = dULDr[2]*drdz + dULDs[2]*dsdz + dULDt[2]*dtdz;
          ULGradCart[3][2] = dULDr[3]*drdz + dULDs[3]*dsdz + dULDt[3]*dtdz;
          ULGradCart[4][2] = dULDr[4]*drdz + dULDs[4]*dsdz + dULDt[4]*dtdz;

          /* Determine the normal gradients of all variables. */
          su2double ULGradNorm[5];
          ULGradNorm[0] = ULGradCart[0][0]*normals[0] + ULGradCart[0][1]*normals[1] + ULGradCart[0][2]*normals[2];
          ULGradNorm[1] = ULGradCart[1][0]*normals[0] + ULGradCart[1][1]*normals[1] + ULGradCart[1][2]*normals[2];
          ULGradNorm[2] = ULGradCart[2][0]*normals[0] + ULGradCart[2][1]*normals[1] + ULGradCart[2][2]*normals[2];
          ULGradNorm[3] = ULGradCart[3][0]*normals[0] + ULGradCart[3][1]*normals[1] + ULGradCart[3][2]*normals[2];
          ULGradNorm[4] = ULGradCart[4][0]*normals[0] + ULGradCart[4][1]*normals[1] + ULGradCart[4][2]*normals[2];

          /* For the construction of the gradients of the right solution, also
             the Cartesian gradients and the normal gradient of the normal
             momentum is needed. This is computed below. */
          su2double GradCartNormMomL[3];
          GradCartNormMomL[0] = ULGradCart[1][0]*normals[0] + ULGradCart[2][0]*normals[1] + ULGradCart[3][0]*normals[2];
          GradCartNormMomL[1] = ULGradCart[1][1]*normals[0] + ULGradCart[2][1]*normals[1] + ULGradCart[3][1]*normals[2];
          GradCartNormMomL[2] = ULGradCart[1][2]*normals[0] + ULGradCart[2][2]*normals[1] + ULGradCart[3][2]*normals[2];

          const su2double GradNormNormMomL = ULGradNorm[1]*normals[0]
                                           + ULGradNorm[2]*normals[1]
                                           + ULGradNorm[3]*normals[2];

          /* Abbreviate twice the normal vector. */
          const su2double tnx = 2.0*normals[0], tny = 2.0*normals[1], tnz = 2.0*normals[2];

          /*--- Construct the gradients of the right solution. The tangential
                gradients of the normal momentum and normal gradients of the other
                variables, density, energy and tangential momentum, must be negated.
                For the common situation that the symmetry plane coincides with an
                x-, y- or z-plane this boils down to a multiplication by -1 or +1,
                depending on the variable and gradient. However, the implementation
                below is valid for an arbitrary orientation of the symmetry plane. ---*/
          su2double URGradCart[5][3];
          URGradCart[0][0] = ULGradCart[0][0] - tnx*ULGradNorm[0];
          URGradCart[1][0] = ULGradCart[1][0] - tnx*ULGradNorm[1] - tnx*GradCartNormMomL[0] + tnx*tnx*GradNormNormMomL;
          URGradCart[2][0] = ULGradCart[2][0] - tnx*ULGradNorm[2] - tny*GradCartNormMomL[0] + tnx*tny*GradNormNormMomL;
          URGradCart[3][0] = ULGradCart[3][0] - tnx*ULGradNorm[3] - tnz*GradCartNormMomL[0] + tnx*tnz*GradNormNormMomL;
          URGradCart[4][0] = ULGradCart[4][0] - tnx*ULGradNorm[4];

          URGradCart[0][1] = ULGradCart[0][1] - tny*ULGradNorm[0];
          URGradCart[1][1] = ULGradCart[1][1] - tny*ULGradNorm[1] - tnx*GradCartNormMomL[1] + tny*tnx*GradNormNormMomL;
          URGradCart[2][1] = ULGradCart[2][1] - tny*ULGradNorm[2] - tny*GradCartNormMomL[1] + tny*tny*GradNormNormMomL;
          URGradCart[3][1] = ULGradCart[3][1] - tny*ULGradNorm[3] - tnz*GradCartNormMomL[1] + tny*tnz*GradNormNormMomL;
          URGradCart[4][1] = ULGradCart[4][1] - tny*ULGradNorm[4];

          URGradCart[0][2] = ULGradCart[0][2] - tnz*ULGradNorm[0];
          URGradCart[1][2] = ULGradCart[1][2] - tnz*ULGradNorm[1] - tnx*GradCartNormMomL[2] + tnz*tnx*GradNormNormMomL;
          URGradCart[2][2] = ULGradCart[2][2] - tnz*ULGradNorm[2] - tny*GradCartNormMomL[2] + tnz*tny*GradNormNormMomL;
          URGradCart[3][2] = ULGradCart[3][2] - tnz*ULGradNorm[3] - tnz*GradCartNormMomL[2] + tnz*tnz*GradNormNormMomL;
          URGradCart[4][2] = ULGradCart[4][2] - tnz*ULGradNorm[4];

          /*--- Compute the viscous fluxes of the left and right state and average
                them to compute the correct viscous flux for a symmetry boundary. ---*/
          const su2double wallDist = wallDistance ? wallDistance[i] : 0.0;
          su2double viscFluxL[5], viscFluxR[5];
          su2double Viscosity, kOverCv;

          ViscousNormalFluxIntegrationPoint_3D(UR, URGradCart, normals, 0.0,
                                               factHeatFlux_Lam, factHeatFlux_Turb,
                                               wallDist, lenScale_LES, Viscosity,
                                               kOverCv, viscFluxR);
          ViscousNormalFluxIntegrationPoint_3D(UL, ULGradCart, normals, 0.0,
                                               factHeatFlux_Lam, factHeatFlux_Turb,
                                               wallDist, lenScale_LES, Viscosity,
                                               kOverCv, viscFluxL);
          viscosityInt[i] = Viscosity;
          kOverCvInt[i]   = kOverCv;

          su2double *viscNormalFlux = viscFluxes + i*ctc::nVar;
          viscNormalFlux[0] = 0.5*(viscFluxL[0] + viscFluxR[0]);
          viscNormalFlux[1] = 0.5*(viscFluxL[1] + viscFluxR[1]);
          viscNormalFlux[2] = 0.5*(viscFluxL[2] + viscFluxR[2]);
          viscNormalFlux[3] = 0.5*(viscFluxL[3] + viscFluxR[3]);
          viscNormalFlux[4] = 0.5*(viscFluxL[4] + viscFluxR[4]);
        }

        break;
      }
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

  su2double *solIntL      = VecTmpMemory.data();
  su2double *solIntR      = solIntL      + nIntegrationMax*nVar;
  su2double *viscosityInt = solIntR      + nIntegrationMax*nVar;
  su2double *kOverCvInt   = viscosityInt + nIntegrationMax;
  su2double *gradSolInt   = kOverCvInt   + nIntegrationMax;
  su2double *fluxes       = gradSolInt   + sizeGradSolInt;
  su2double *viscFluxes   = fluxes       + sizeFluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the length scale for the LES of the adjacent element. */
    const unsigned long  elemID = surfElem[l].volElemID;
    const unsigned short iind   = volElem[elemID].indStandardElement;

    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

    /*--- Apply the inviscid wall boundary conditions to compute the right
          state in the integration points. There are two options. Either the
          normal velocity is negated or the normal velocity is set to zero.
          Some experiments are needed to see which formulation gives better
          results. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution and the normals
         for this integration point. */
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;
      const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

      /* Compute the normal component of the momentum variables. */
      su2double rVn = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        rVn += UL[iDim+1]*normals[iDim];

      /* If the normal velocity must be mirrored instead of set to zero,
         the normal component that must be subtracted must be doubled. If the
         normal velocity must be set to zero, simply comment this line. */
      //rVn *= 2.0;

      /* Set the right state. The initial value of the total energy is the
         energy of the left state. */
      UR[0]      = UL[0];
      UR[nDim+1] = UL[nDim+1];
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        UR[iDim+1] = UL[iDim+1] - rVn*normals[iDim];

      /*--- Actually, only the internal energy of UR is equal to UL. If the
            kinetic energy differs for UL and UR, the difference must be
            subtracted from the total energy of UR to obtain the correct
            value. ---*/
      su2double DensityInv = 1.0/UL[0];
      su2double diffKin    = 0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double velL = DensityInv*UL[iDim];
        const su2double velR = DensityInv*UR[iDim];
        diffKin += velL*velL - velR*velR;
      }

      UR[nDim+1] -= 0.5*UL[0]*diffKin;
    }

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, nInt, nDOFsElem, 0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(), lenScale_LES,
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, numerics, &surfElem[l], solIntL, solIntR,
                                gradSolInt, fluxes, viscFluxes, viscosityInt,
                                kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Retrieve the specified total conditions for this inlet. ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double P_Total   = config->GetInlet_Ptotal(Marker_Tag);
  su2double T_Total   = config->GetInlet_Ttotal(Marker_Tag);
  su2double *Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

  /*--- Non-dim. the inputs if necessary, and compute the total enthalpy. ---*/
  P_Total /= config->GetPressure_Ref();
  T_Total /= config->GetTemperature_Ref();

  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double H_Total      = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;

  /*--- Set the pointers for the local arrays. ---*/
  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

  su2double *solIntL      = VecTmpMemory.data();
  su2double *solIntR      = solIntL      + nIntegrationMax*nVar;
  su2double *viscosityInt = solIntR      + nIntegrationMax*nVar;
  su2double *kOverCvInt   = viscosityInt + nIntegrationMax;
  su2double *gradSolInt   = kOverCvInt   + nIntegrationMax;
  su2double *fluxes       = gradSolInt   + sizeGradSolInt;
  su2double *viscFluxes   = fluxes       + sizeFluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the length scale for the LES of the adjacent element. */
    const unsigned long  elemID = surfElem[l].volElemID;
    const unsigned short iind   = volElem[elemID].indStandardElement;

    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

    /*--- Apply the subsonic inlet boundary conditions to compute the right
          state in the integration points. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution and the normals
         for this integration point. */
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;
      const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

      /*--- Compute the normal velocity, the speed of sound squared and
            the pressure in this integration point. ---*/
      const su2double DensityInv = 1.0/UL[0];
      su2double VelocityNormal = 0.0, Velocity2 = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        const su2double vel = UL[iDim+1]*DensityInv;
        VelocityNormal     += vel*normals[iDim];
        Velocity2          += vel*vel;
      }

      su2double StaticEnergy = UL[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(UL[0], StaticEnergy);
      su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
      su2double Pressure    = FluidModel->GetPressure();

      /*--- Compute the Riemann invariant to be extrapolated. ---*/
      const su2double Riemann = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One + VelocityNormal;

      /*--- Total speed of sound. The expression below is also valid for variable cp,
            although a linearization around the value of the left state is performed. ---*/
      const su2double SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - DensityInv*(UL[nDim+1] + Pressure)
                                        +                  0.5*Velocity2) + SoundSpeed2;

      /*--- Dot product of normal and flow direction. This should be
            negative due to outward facing boundary normal convention. ---*/
      su2double alpha = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        alpha += normals[iDim]*Flow_Dir[iDim];

      /*--- Coefficients in the quadratic equation for the magnitude of the velocity. ---*/
      const su2double aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
      const su2double bb = -1.0*Gamma_Minus_One*alpha*Riemann;
      const su2double cc =  0.5*Gamma_Minus_One*Riemann*Riemann
                         -  2.0*SoundSpeed_Total2/Gamma_Minus_One;

      /*--- Solve the equation for the magnitude of the velocity. As this value
            must be positive and both aa and bb are positive (alpha is negative and
            Riemann is positive up till Mach = 5.0 or so, which is not really subsonic
            anymore), it is clear which of the two possible solutions must be taken.
            Some clipping is present, but this is normally not active. ---*/
      su2double dd      = bb*bb - 4.0*aa*cc;   dd      = sqrt(max(0.0, dd));
      su2double Vel_Mag = (-bb + dd)/(2.0*aa); Vel_Mag = max(0.0, Vel_Mag);
      Velocity2 = Vel_Mag*Vel_Mag;

      /*--- Compute speed of sound from total speed of sound eqn. ---*/
      SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

      /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
      su2double Mach2 = Velocity2/SoundSpeed2; Mach2 = min(1.0, Mach2);

      Velocity2 = Mach2*SoundSpeed2;
      Vel_Mag   = sqrt(Velocity2);

      /*--- Static temperature from the speed of sound relation. ---*/
      SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

      const su2double Temperature = SoundSpeed2/(Gamma*Gas_Constant);

      /*--- Static pressure using isentropic relation at a point ---*/
      Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);

      /*--- Density at the inlet from the gas law ---*/
      const su2double Density = Pressure/(Gas_Constant*Temperature);

      /*--- Store the conservative variables in UR. ---*/
      UR[0]      = Density;
      UR[nDim+1] = Pressure/Gamma_Minus_One + 0.5*Density*Velocity2;

      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        UR[iDim+1] = Density*Vel_Mag*Flow_Dir[iDim];
    }

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, nInt, nDOFsElem, 0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(), lenScale_LES,
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                 CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Retrieve the specified back pressure for this outlet.
        Nondimensionalize, if necessary. ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double P_Exit = config->GetOutlet_Pressure(Marker_Tag);
  P_Exit = P_Exit/config->GetPressure_Ref();

  /*--- Set the pointers for the local arrays. ---*/
  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

  su2double *solIntL      = VecTmpMemory.data();
  su2double *solIntR      = solIntL      + nIntegrationMax*nVar;
  su2double *viscosityInt = solIntR      + nIntegrationMax*nVar;
  su2double *kOverCvInt   = viscosityInt + nIntegrationMax;
  su2double *gradSolInt   = kOverCvInt   + nIntegrationMax;
  su2double *fluxes       = gradSolInt   + sizeGradSolInt;
  su2double *viscFluxes   = fluxes       + sizeFluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the length scale for the LES of the adjacent element. */
    const unsigned long  elemID = surfElem[l].volElemID;
    const unsigned short iind   = volElem[elemID].indStandardElement;

    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

    /*--- Apply the subsonic inlet boundary conditions to compute the right
          state in the integration points. ---*/
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution and the normals
         for this integration point. */
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;
      const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

      /*--- Compute the normal velocity, the speed of sound squared and
            the pressure in this integration point. ---*/
      const su2double DensityInv = 1.0/UL[0];
      su2double VelocityNormal = 0.0, Velocity2 = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        const su2double vel = UL[iDim+1]*DensityInv;
        VelocityNormal     += vel*normals[iDim];
        Velocity2          += vel*vel;
      }

      su2double StaticEnergy = UL[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(UL[0], StaticEnergy);
      su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
      su2double Pressure    = FluidModel->GetPressure();

      /*--- Subsonic exit flow: there is one incoming characteristic,
            therefore one variable can be specified (back pressure) and is used
            to update the conservative variables. Compute the entropy and the
            acoustic Riemann variable. These invariants, as well as the
            tangential velocity components, are extrapolated. ---*/
      const su2double Entropy = Pressure*pow(DensityInv, Gamma);
      const su2double Riemann = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One + VelocityNormal;

      /*--- Compute the density and normal velocity of the right state. ---*/
      const su2double Density    = pow(P_Exit/Entropy,1.0/Gamma);
      const su2double SoundSpeed = sqrt(Gamma*P_Exit/Density);
      const su2double Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;

      /*--- Store the conservative variables in UR. ---*/
      Velocity2 = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        const su2double vel = UL[iDim+1]*DensityInv
                            + (Vn_Exit-VelocityNormal)*normals[iDim];
        Velocity2 += vel*vel;
        UR[iDim+1] = Density*vel;
      }

      UR[0]      = Density;
      UR[nDim+1] = Pressure/Gamma_Minus_One + 0.5*Density*Velocity2;
    }

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, nInt, nDOFsElem, 0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(), lenScale_LES,
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /* Set the factor for the wall velocity. For factWallVel = 0, the right state
     contains the wall velocity. For factWallVel = 1.0, the velocity of the
     right state is obtained by negating the interior velocity w.r.t. the
     velocity of the wall. */
  const su2double factWallVel = 0.0;
  // const su2double factWallVel = 1.0;

  /* Get the wall heat flux. */
  const string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  const su2double Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  /*--- Set the pointers for the local arrays. ---*/
  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

  su2double *solIntL      = VecTmpMemory.data();
  su2double *solIntR      = solIntL      + nIntegrationMax*nVar;
  su2double *viscosityInt = solIntR      + nIntegrationMax*nVar;
  su2double *kOverCvInt   = viscosityInt + nIntegrationMax;
  su2double *gradSolInt   = kOverCvInt   + nIntegrationMax;
  su2double *fluxes       = gradSolInt   + sizeGradSolInt;
  su2double *viscFluxes   = fluxes       + sizeFluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the length scale for the LES of the adjacent element. */
    const unsigned long  elemID = surfElem[l].volElemID;
    const unsigned short iind   = volElem[elemID].indStandardElement;

    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

    /* Determine the number of integration points. */
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Apply the heat flux wall boundary conditions to compute the right
          state. There are two options. Either the velocity is negated or it
          is set to zero. Some experiments are needed to see which formulation
          gives better results. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution for this integration point. */
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;

      /* Set the right state. The initial value of the total energy is the
         energy of the left state. Also compute the difference in kinetic
         energy between the left and right state. */
      UR[0]      = UL[0];
      UR[nDim+1] = UL[nDim+1];

      su2double DensityInv = 1.0/UL[0];
      su2double diffKin    = 0.0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        UR[iDim] = -factWallVel*UL[iDim];
        const su2double velL = DensityInv*UL[iDim];
        const su2double velR = DensityInv*UR[iDim];
        diffKin += velL*velL - velR*velR;
      }

      /* As only the internal energy of UR is equal to UL, the difference
         in kinetic energy must be subtracted. */
      UR[nDim+1] -= 0.5*UR[0]*diffKin;
    }

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, nInt, nDOFsElem, Wall_HeatFlux, true, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(), lenScale_LES,
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /* Set the factor for the wall velocity. For factWallVel = 0, the right state
     contains the wall velocity. For factWallVel = 1.0, the velocity of the
     right state is obtained by negating the interior velocity w.r.t. the
     velocity of the wall. */
  const su2double factWallVel = 0.0;
  // const su2double factWallVel = 1.0;

  /* Get the wall temperature. */
  const string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  const su2double TWall   = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

  /* Compute the prescribed value of the energy (per unit mass). */
  const su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double Cv           = Gas_Constant/Gamma_Minus_One;
  const su2double StaticEnergy = Cv*TWall;

  /*--- Set the pointers for the local arrays. ---*/
  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

  su2double *solIntL      = VecTmpMemory.data();
  su2double *solIntR      = solIntL      + nIntegrationMax*nVar;
  su2double *viscosityInt = solIntR      + nIntegrationMax*nVar;
  su2double *kOverCvInt   = viscosityInt + nIntegrationMax;
  su2double *gradSolInt   = kOverCvInt   + nIntegrationMax;
  su2double *fluxes       = gradSolInt   + sizeGradSolInt;
  su2double *viscFluxes   = fluxes       + sizeFluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the length scale for the LES of the adjacent element. */
    const unsigned long  elemID = surfElem[l].volElemID;
    const unsigned short iind   = volElem[elemID].indStandardElement;

    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

    /* Determine the number of integration points. */
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Apply the heat flux wall boundary conditions to compute the right
          state. There are two options. Either the velocity is negated or it
          is set to zero. Some experiments are needed to see which formulation
          gives better results. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution for this integration point. */
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;

      /* Set the right state for the density and the momentum variables of the
         right state. Compute twice the possible kinetic energy. */
      UR[0] = UL[0];
      su2double DensityInv = 1.0/UL[0];
      su2double kinEner    = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        UR[iDim+1] = -factWallVel*UL[iDim+1];
        const su2double velR = DensityInv*UR[iDim];
        kinEner += velR*velR;
      }

      /* Compute the total energy of the right state. */
      UR[nDim+1] = UR[0]*(StaticEnergy + 0.5*kinEner);
    }

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, nInt, nDOFsElem, 0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(), lenScale_LES,
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {

#ifdef CUSTOM_BC_NSUNITQUAD
  /* Get the flow angle, which is stored in the angle of attack and the
     viscosity coefficient. */
  const su2double flowAngle = config->GetAoA()*PI_NUMBER/180.0;
  const su2double mu        = config->GetViscosity_FreeStreamND();

  const su2double cosFlowAngle = cos(flowAngle);
  const su2double sinFlowAngle = sin(flowAngle);
#endif

  /*--- Set the pointers for the local arrays. ---*/
  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

  su2double *solIntL      = VecTmpMemory.data();
  su2double *solIntR      = solIntL      + nIntegrationMax*nVar;
  su2double *viscosityInt = solIntR      + nIntegrationMax*nVar;
  su2double *kOverCvInt   = viscosityInt + nIntegrationMax;
  su2double *gradSolInt   = kOverCvInt   + nIntegrationMax;
  su2double *fluxes       = gradSolInt   + sizeGradSolInt;
  su2double *viscFluxes   = fluxes       + sizeFluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker][0];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the length scale for the LES of the adjacent element. */
    const unsigned long  elemID = surfElem[l].volElemID;
    const unsigned short iind   = volElem[elemID].indStandardElement;

    unsigned short nPoly = standardElementsSol[iind].GetNPoly();
    if(nPoly == 0) nPoly = 1;

    const su2double lenScale_LES = volElem[elemID].lenScale/nPoly;

    /* Determine the number of integration points. */
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Loop over the integration points to compute the right state via the
          customized boundary conditions. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

#ifdef CUSTOM_BC_NSUNITQUAD

      /* Easier storage of the right solution for this integration point and
         the coordinates of this integration point. */
      const su2double *coor = surfElem[l].coorIntegrationPoints + i*nDim;
            su2double *UR   = solIntR + i*nVar;

      /*--- Set the exact solution in this integration point. ---*/
      const double xTilde = coor[0]*cosFlowAngle - coor[1]*sinFlowAngle;
      const double yTilde = coor[0]*sinFlowAngle + coor[1]*cosFlowAngle;

      UR[0] =  1.0;
      UR[1] =  cosFlowAngle*yTilde*yTilde;
      UR[2] = -sinFlowAngle*yTilde*yTilde;
      UR[3] =  (2.0*mu*xTilde + 10)/Gamma_Minus_One
            +  0.5*yTilde*yTilde*yTilde*yTilde;
#else

      /* No compiler directive specified. Write an error message and exit. */
      int rank = MASTER_NODE;
#ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

      if (rank == MASTER_NODE) {
        cout << endl;
        cout << "In function CFEM_DG_NSSolver::BC_Custom. " << endl;
        cout << "No or wrong compiler directive specified. This is necessary "
                "for customized boundary conditions." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif

#endif
    }

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
    const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, nInt, nDOFsElem, 0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(), lenScale_LES,
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, numerics, &surfElem[l], solIntL, solIntR,
                                gradSolInt, fluxes, viscFluxes, viscosityInt,
                                kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::ResidualViscousBoundaryFace(
                                      CConfig                  *config,
                                      CNumerics                *conv_numerics,
                                      const CSurfaceElementFEM *surfElem,
                                      const su2double          *solInt0,
                                      const su2double          *solInt1,
                                      su2double                *gradSolInt,
                                      su2double                *fluxes,
                                      su2double                *viscFluxes,
                                      const su2double          *viscosityInt,
                                      const su2double          *kOverCvInt,
                                      su2double                *resFaces,
                                      unsigned long            &indResFaces) {

  su2double tick = 0.0;

  /* Determine whether or not the Cartesian gradients of the basis functions
     are stored. */
  const bool CartGradBasisFunctionsStored = config->GetStore_Cart_Grad_BasisFunctions_DGFEM();
  
  /*--- Get the required information from the standard element. ---*/
  const unsigned short ind          = surfElem->indStandardElement;
  const unsigned short nInt         = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFs        = standardBoundaryFacesSol[ind].GetNDOFsFace();
  const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
  const su2double      ConstPenFace = standardBoundaryFacesSol[ind].GetPenaltyConstant();
  const su2double     *weights      = standardBoundaryFacesSol[ind].GetWeightsIntegration();

  /* General function to compute the inviscid fluxes, using an approximate
     Riemann solver in the integration points. */
  ComputeInviscidFluxesFace(config, nInt, surfElem->metricNormalsFace.data(),
                            solInt0, solInt1, fluxes, conv_numerics);

  /* Subtract the viscous fluxes from the inviscid fluxes. */
  for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] -= viscFluxes[j];

  /* Get the length scale for the adjacent element. */
  const su2double lenScale = volElem[surfElem->volElemID].lenScale;

  /* Call the function PenaltyTermsFluxFace to compute the actual penalty
     terms. Use the array viscFluxes as storage. */
  PenaltyTermsFluxFace(nInt, solInt0, solInt1, viscosityInt, viscosityInt,
                       kOverCvInt, kOverCvInt, ConstPenFace, lenScale, lenScale,
                       surfElem->metricNormalsFace.data(), viscFluxes);

  /* Add the penalty fluxes to the earlier computed fluxes. */
  for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] += viscFluxes[j];

  /* Multiply the fluxes with the integration weight of the corresponding
     integration point. */
  for(unsigned short i=0; i<nInt; ++i) {
    su2double *flux = fluxes + i*nVar;

    for(unsigned short j=0; j<nVar; ++j)
      flux[j] *= weights[i];
  }

  /*------------------------------------------------------------------------*/
  /*--- Step 2: Compute the contribution to the residuals from the       ---*/
  /*---         integration over the surface element of the invisid      ---*/
  /*---         fluxes, viscous fluxes and penalty terms.                ---*/
  /*------------------------------------------------------------------------*/

  /* Easier storage of the position in the residual array for this face
     and update the corresponding counter. */
  su2double *resFace = resFaces + indResFaces*nVar;
  indResFaces       += nDOFs;

  /* Get the correct form of the basis functions needed for the matrix
     multiplication to compute the residual. */
  const su2double *basisFaceTrans = standardBoundaryFacesSol[ind].GetBasisFaceIntegrationTranspose();

  /* Call the general function to carry out the matrix product. */
  config->GEMM_Tick(&tick);
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFaceTrans, fluxes, resFace);
  config->GEMM_Tock(tick, "ResidualViscousBoundaryFace1", nDOFs, nVar, nInt);

  /*------------------------------------------------------------------------*/
  /*--- Step 3: Compute the symmetrizing terms, if present, in the       ---*/
  /*---         integration points of this boundary face.                ---*/
  /*------------------------------------------------------------------------*/

  if( symmetrizingTermsPresent ) {

    /* Compute the symmetrizing fluxes in the nDim directions. */
    SymmetrizingFluxesFace(nInt, solInt0, solInt1, viscosityInt,
                           viscosityInt, kOverCvInt, kOverCvInt,
                           surfElem->metricNormalsFace.data(), fluxes);

    /*--- Multiply the fluxes just computed by their integration weights and
          -theta/2. The parameter theta is the parameter in the Interior Penalty
          formulation, the factor 1/2 comes in from the averaging and the minus
          sign is from the convention that the viscous fluxes comes with a minus
          sign in this code. ---*/
    const su2double halfTheta = 0.5*config->GetTheta_Interior_Penalty_DGFEM();

    for(unsigned short i=0; i<nInt; ++i) {
      su2double *flux        = fluxes + i*nVar*nDim;
      const su2double wTheta = -halfTheta*weights[i];

      for(unsigned short j=0; j<(nVar*nDim); ++j)
        flux[j] *= wTheta;
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 4: Distribute the symmetrizing terms to the DOFs. Note that ---*/
    /*---         these terms must be distributed to all the DOFs of the   ---*/
    /*---         adjacent element, not only to the DOFs of the face.      ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the position in the residual array for this face
       and update the corresponding counter. */
    su2double *resElem = resFaces + indResFaces*nVar;
    indResFaces       += nDOFsElem;

    /* The Cartesian gradients of the basis functions of the adjacent element
         are needed. Check if this data is stored. */
    const su2double *cartGrad;
    if( CartGradBasisFunctionsStored )
      cartGrad = surfElem->metricElem.data();
    else {

      /* Get the derivatives of the basis functions w.r.t. the parametric
         coordinates of the element. */
      const su2double *derBasisElemTrans = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegrationTranspose();

      /*--- Create the Cartesian derivatives of the basis functions in the
            integration points. Use gradSolInt for storage. ---*/
      unsigned int ii = 0;
      for(unsigned short j=0; j<nDOFsElem; ++j) {
        for(unsigned short i=0; i<nInt; ++i, ii+=ctc::nDim) {

          /* Easier storage of the derivatives of the basis function w.r.t. the
             parametric coordinates, the location where to store the Cartesian
             derivatives of the basis functions, and the metric terms in this
             integration point. */
          const su2double *derParam    = derBasisElemTrans + ii;
          const su2double *metricTerms = surfElem->metricCoorDerivFace.data()
                                       + i*ctc::nDim*ctc::nDim;
                su2double *derCar      = gradSolInt + ii;

          /*--- Loop over the dimensions to compute the Cartesian derivatives
                of the basis functions. ---*/
          for(unsigned short k=0; k<ctc::nDim; ++k) {
            derCar[k] = 0.0;
            for(unsigned short l=0; l<ctc::nDim; ++l)
              derCar[k] += derParam[l]*metricTerms[k+l*ctc::nDim];
          }
        }
      }

      /* Set the pointer of gradSolInt to cartGrad, such that the latter can
           be used in the call to the matrix multiplication. */
        cartGrad = gradSolInt;
    }

    /* Call the general function to carry out the matrix product to compute
       the residual. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nDOFsElem, nVar, nInt*nDim, cartGrad, fluxes, resElem);
    config->GEMM_Tock(tick, "ResidualViscousBoundaryFace2", nDOFsElem, nVar, nInt*nDim);

  }
}
