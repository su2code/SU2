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
}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CConfig *config, unsigned short val_nDim, unsigned short iMesh) : CSolver() {

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

  /*--- Dummy solver constructor that calls the SetNondim. routine in
        order to load the flow non-dim. information into the config class.
        This is needed to complete a partitioning for time-accurate local
        time stepping that depends on the flow state. ---*/
  nDim = val_nDim;
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  SetNondimensionalization(config, iMesh, false);
}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

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
  nVolElemHaloPerTimeLevel     = DGGeometry->GetNVolElemHaloPerTimeLevel();

  ownedElemAdjLowTimeLevel = DGGeometry->GetOwnedElemAdjLowTimeLevel();
  haloElemAdjLowTimeLevel  = DGGeometry->GetHaloElemAdjLowTimeLevel();

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

  LagrangianBeginTimeIntervalADER_DG     = DGGeometry->GetLagrangianBeginTimeIntervalADER_DG();
  timeInterpolDOFToIntegrationADER_DG    = DGGeometry->GetTimeInterpolDOFToIntegrationADER_DG();
  timeInterpolAdjDOFToIntegrationADER_DG = DGGeometry->GetTimeInterpolAdjDOFToIntegrationADER_DG();

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

    sizeVecTmp = nIntegrationMax*(4 + 3*nVar) + sizeFluxes + sizeGradSolInt
               + max(nIntegrationMax,nDOFsMax)*nVar;
  }
  else {

    /* Inviscid simulation. */
    unsigned int sizeVol = nVar*nIntegrationMax*(nDim+2);
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

  VecSolDOFs.resize(nVar*nDOFsLocOwned);

  /*--- Allocate the memory for the working vectors for the solution variables
        for all the time levels used. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
  VecWorkSolDOFs.resize(nTimeLevels);

  vector<unsigned long> nDOFsTimeLevels(nTimeLevels, 0);
  for(unsigned long i=0; i<nVolElemTot; ++i) {
    nDOFsTimeLevels[volElem[i].timeLevel] += volElem[i].nDOFsSol;
    if(volElem[i].offsetDOFsSolPrevTimeLevel != ULONG_MAX)
      nDOFsTimeLevels[volElem[i].timeLevel-1] += volElem[i].nDOFsSol;
  }

  for(unsigned short i=0; i<nTimeLevels; ++i)
    VecWorkSolDOFs[i].resize(nVar*nDOFsTimeLevels[i]);

  /*--- Check for the ADER-DG time integration scheme and allocate the memory
        for the additional vectors. ---*/
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {

    const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();

    VecTotResDOFsADER.assign(nVar*nDOFsLocTot, 0.0);
    VecSolDOFsPredictorADER.resize(nTimeDOFs*nVar*nDOFsLocTot);
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

  /*--- Store the number of DOFs in the geometry class in case of restart. ---*/
  geometry->SetnPointDomain(nDOFsLocOwned);
  geometry->SetGlobal_nPointDomain(nDOFsGlobal);

  /*--- Allocate the memory to store the time steps, residuals, etc. ---*/
  VecDeltaTime.resize(nVolElemOwned);

  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG)
    VecResDOFs.resize(nVar*nDOFsLocOwned);
  else
    VecResDOFs.resize(nVar*nDOFsLocTot);

  nEntriesResFaces.assign(nDOFsLocTot+1, 0);
  nEntriesResAdjFaces.assign(nDOFsLocTot+1, 0);
  startLocResFacesMarkers.resize(nMarker);

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

    /* Determine the time level of the face. */
    const unsigned long  elem0     = matchingInternalFaces[i].elemID0;
    const unsigned long  elem1     = matchingInternalFaces[i].elemID1;
    const unsigned short timeLevel = min(volElem[elem0].timeLevel,
                                         volElem[elem1].timeLevel);

    /* The terms that only contribute to the DOFs located on the face.
       Check the time level of the element compared to the face. */
    const unsigned short ind = matchingInternalFaces[i].indStandardElement;
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    sizeVecResFaces += nDOFsFace0;
    if(timeLevel == volElem[elem0].timeLevel) {
      for(unsigned short j=0; j<nDOFsFace0; ++j)
        ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolFaceSide0[j]+1];
    }
    else {
      for(unsigned short j=0; j<nDOFsFace0; ++j)
        ++nEntriesResAdjFaces[matchingInternalFaces[i].DOFsSolFaceSide0[j]+1];
    }

    sizeVecResFaces += nDOFsFace1;
    if(timeLevel == volElem[elem1].timeLevel) {
      for(unsigned short j=0; j<nDOFsFace1; ++j)
        ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolFaceSide1[j]+1];
    }
    else {
      for(unsigned short j=0; j<nDOFsFace1; ++j)
        ++nEntriesResAdjFaces[matchingInternalFaces[i].DOFsSolFaceSide1[j]+1];
    }

    /* The symmetrizing terms, if present, contribute to all
       the DOFs of the adjacent elements. */
    if( symmetrizingTermsPresent ) {
      const unsigned short nDOFsElem0 = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
      const unsigned short nDOFsElem1 = standardMatchingFacesSol[ind].GetNDOFsElemSide1();

      sizeVecResFaces += nDOFsElem0;
      if(timeLevel == volElem[elem0].timeLevel) {
        for(unsigned short j=0; j<nDOFsElem0; ++j)
          ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolElementSide0[j]+1];
      }
      else {
        for(unsigned short j=0; j<nDOFsElem0; ++j)
          ++nEntriesResAdjFaces[matchingInternalFaces[i].DOFsSolElementSide0[j]+1];
      }

      sizeVecResFaces += nDOFsElem1;
      if(timeLevel == volElem[elem1].timeLevel) {
        for(unsigned short j=0; j<nDOFsElem1; ++j)
          ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolElementSide1[j]+1];
      }
      else {
        for(unsigned short j=0; j<nDOFsElem1; ++j)
          ++nEntriesResAdjFaces[matchingInternalFaces[i].DOFsSolElementSide1[j]+1];
      }
    }

    /* Store the position of the residual in the appropriate entry. */
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

  /*--- Put nEntriesResFaces and nEntriesResAdjFaces in cumulative storage
        format and allocate the memory for entriesResFaces, entriesResAdjFaces
        and VecResFaces. ---*/
  for(unsigned long i=0; i<nDOFsLocTot; ++i) {
    nEntriesResFaces[i+1]    += nEntriesResFaces[i];
    nEntriesResAdjFaces[i+1] += nEntriesResAdjFaces[i];
  }

  entriesResFaces.resize(nEntriesResFaces[nDOFsLocTot]);
  entriesResAdjFaces.resize(nEntriesResAdjFaces[nDOFsLocTot]);
  VecResFaces.resize(nVar*sizeVecResFaces);

  /*--- Repeat the loops over the internal and boundary faces, but now store
        the entries in entriesResFaces and entriesResAdjFaces. A counter
        variable is needed to keep track of the appropriate location in
        entriesResFaces and entriesResAdjFaces. ---*/
  vector<unsigned long> counterEntries    = nEntriesResFaces;
  vector<unsigned long> counterEntriesAdj = nEntriesResAdjFaces;

  /* First the loop over the internal matching faces. */
  sizeVecResFaces = 0;
  for(unsigned long i=0; i<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++i) {

    /* Determine the time level of the face. */
    const unsigned long  elem0     = matchingInternalFaces[i].elemID0;
    const unsigned long  elem1     = matchingInternalFaces[i].elemID1;
    const unsigned short timeLevel = min(volElem[elem0].timeLevel,
                                         volElem[elem1].timeLevel);

    /* The terms that only contribute to the DOFs located on the face. */
    const unsigned short ind = matchingInternalFaces[i].indStandardElement;
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    if(timeLevel == volElem[elem0].timeLevel) {
      for(unsigned short j=0; j<nDOFsFace0; ++j) {
        const unsigned long jj = counterEntries[matchingInternalFaces[i].DOFsSolFaceSide0[j]]++;
        entriesResFaces[jj]    = sizeVecResFaces++;
      }
    }
    else {
      for(unsigned short j=0; j<nDOFsFace0; ++j) {
        const unsigned long jj = counterEntriesAdj[matchingInternalFaces[i].DOFsSolFaceSide0[j]]++;
        entriesResAdjFaces[jj] = sizeVecResFaces++;
      }
    }

    if(timeLevel == volElem[elem1].timeLevel) {
      for(unsigned short j=0; j<nDOFsFace1; ++j) {
        const unsigned long jj = counterEntries[matchingInternalFaces[i].DOFsSolFaceSide1[j]]++;
        entriesResFaces[jj]    = sizeVecResFaces++;
      }
    }
    else {
      for(unsigned short j=0; j<nDOFsFace1; ++j) {
        const unsigned long jj = counterEntriesAdj[matchingInternalFaces[i].DOFsSolFaceSide1[j]]++;
        entriesResAdjFaces[jj] = sizeVecResFaces++;
      }
    }

    /* The symmetrizing terms, if present, contribute to all
       the DOFs of the adjacent elements. */
    if( symmetrizingTermsPresent ) {
      const unsigned short nDOFsElem0 = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
      const unsigned short nDOFsElem1 = standardMatchingFacesSol[ind].GetNDOFsElemSide1();

      if(timeLevel == volElem[elem0].timeLevel) {
        for(unsigned short j=0; j<nDOFsElem0; ++j) {
          const unsigned long jj = counterEntries[matchingInternalFaces[i].DOFsSolElementSide0[j]]++;
          entriesResFaces[jj]    = sizeVecResFaces++;
        }
      }
      else {
        for(unsigned short j=0; j<nDOFsElem0; ++j) {
          const unsigned long jj = counterEntriesAdj[matchingInternalFaces[i].DOFsSolElementSide0[j]]++;
          entriesResAdjFaces[jj] = sizeVecResFaces++;
        }
      }

      if(timeLevel == volElem[elem1].timeLevel) {
        for(unsigned short j=0; j<nDOFsElem1; ++j) {
          const unsigned long jj = counterEntries[matchingInternalFaces[i].DOFsSolElementSide1[j]]++;
          entriesResFaces[jj]    = sizeVecResFaces++;
        }
      }
      else {
        for(unsigned short j=0; j<nDOFsElem1; ++j) {
          const unsigned long jj = counterEntriesAdj[matchingInternalFaces[i].DOFsSolElementSide1[j]]++;
          entriesResAdjFaces[jj] = sizeVecResFaces++;
        }
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

  /*--- Start the solution from the free-stream state. This is overruled
        when a restart takes place. ---*/
  unsigned long ii = 0;
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
    for(unsigned short j=0; j<nVar; ++j, ++ii) {
      VecSolDOFs[ii] = ConsVarFreeStream[j];
    }
  }

  /* Check if the exact Jacobian of the spatial discretization must be
     determined. If so, the color of each DOF must be determined, which
     is converted to the DOFs for each color. */
  if( config->GetJacobian_Spatial_Discretization_Only() ) {

    /* Write a message that the graph coloring is performed. */
    if(rank == MASTER_NODE)
      cout << "Creating the vertex colors of the graph. " << std::endl;

    /* Determine the graph of the stencil of every local DOF. */
    DetermineGraphDOFs(DGGeometry, config);

    /* Carry out the vertex coloring of the graph. */
    vector<int> colorLocalDOFs;
    GraphVertexColoring(config, nDOFsPerRank, nonZeroEntriesJacobian,
                        nGlobalColors, colorLocalDOFs);

    /* Write a message that the all volume DOFs have been colored. */
    if(rank == MASTER_NODE)
      cout << "There are " << nGlobalColors
           << " colors present in the graph for " << nDOFsPerRank.back()
           << " DOFs." << std::endl;

    /* Determine the meta data needed for the computation of the Jacobian. */
    MetaDataJacobianComputation(DGGeometry, colorLocalDOFs);
  }

  /* Set up the persistent communication for the conservative variables and
     the reverse communication for the residuals of the halo elements. */
  Prepare_MPI_Communication(DGGeometry, config);

  /* Set up the list of tasks to be carried out in the computational expensive
     part of the code. For the Runge-Kutta schemes this is typically the
     computation of the spatial residual, while for ADER this list contains
     the tasks to be done for one space time step. */
  SetUpTaskList(config);
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
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);

      FluidModel = new CIdealGas(1.4, config->GetGas_Constant(), config->GetCompute_Entropy());
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

      FluidModel = new CIdealGas(Gamma, config->GetGas_Constant(), config->GetCompute_Entropy());
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

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
          To accomplish this, simply set the non-dimensional coefficients to the
          dimensional ones. This will be overruled later.---*/
    config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
    config->SetMu_SND(config->GetMu_S());

    config->SetMu_ConstantND(config->GetMu_Constant());

    /*--- Reynolds based initialization ---*/

    if (reynolds_init) {

      /*--- First, check if there is mesh motion. If yes, use the Mach
       number relative to the body to initialize the flow. ---*/

      if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
      else Velocity_Reynolds = ModVel_FreeStream;

      /*--- For viscous flows, pressure will be computed from a density
            that is found from the Reynolds number. The viscosity is computed
            from the dimensional version of Sutherland's law or the constant
            viscosity, depending on the input option.---*/

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

  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/max((Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream()), 1.e-25);
  config->SetOmega_FreeStream(Omega_FreeStream);

  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/max((Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream()), 1.e-25);
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);

  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/

  /*--- Delete the original (dimensional) FluidModel object before replacing. ---*/

  delete FluidModel;

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:
      FluidModel = new CIdealGas(1.4, Gas_ConstantND, config->GetCompute_Entropy());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case IDEAL_GAS:
      FluidModel = new CIdealGas(Gamma, Gas_ConstantND, config->GetCompute_Entropy());
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
    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

    /*--- Sutherland's model ---*/

    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());

    /* constant thermal conductivity model */
    config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);

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
          cout << "Laminar Viscosity: " << config->GetMu_Constant();
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          break;

        case SUTHERLAND:
          cout << "Viscosity Model: SUTHERLAND "<< endl;
          cout << "Ref. Laminar Viscosity: " << config->GetMu_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Ref. Temperature: " << config->GetMu_Temperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Sutherland Constant: "<< config->GetMu_S();
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
          cout << "Molecular Conductivity: " << config->GetKt_Constant()<< " W/m^2.K." << endl;
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

void CFEM_DG_EulerSolver::DetermineGraphDOFs(const CMeshFEM *FEMGeometry,
                                             CConfig        *config) {

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

  /* Determine the global values of the locally owned DOFs. */
  vector<unsigned long> DOFsGlobalID(nDOFsLocTot);
  for(unsigned long i=0; i<nDOFsLocOwned; ++i)
    DOFsGlobalID[i] = nDOFsPerRank[rank] + i;

  /*--- Get the communication information from DG_Geometry. Note that for a
        FEM DG discretization the communication entities of FEMGeometry contain
        the volume elements. ---*/
  const vector<int>                    &ranksSend    = FEMGeometry->GetRanksSend();
  const vector<int>                    &ranksRecv    = FEMGeometry->GetRanksRecv();
  const vector<vector<unsigned long> > &elementsSend = FEMGeometry->GetEntitiesSend();
  const vector<vector<unsigned long> > &elementsRecv = FEMGeometry->GetEntitiesRecv();

#ifdef HAVE_MPI

  /* Parallel implementation. Allocate the memory for the send buffers and
     loop over the number of ranks to which I have to send data. Note that
     self communication is not excluded here. */
  vector<vector<unsigned long> > sendBuf;
  sendBuf.resize(ranksSend.size());

  vector<SU2_MPI::Request> sendReqs(ranksSend.size());

  for(unsigned i=0; i<ranksSend.size(); ++i) {

    /* Fill the send buffer for this rank. */
    for(unsigned long j=0; j<elementsSend[i].size(); ++j) {
      const unsigned long jj = elementsSend[i][j];
      for(unsigned short k=0; k<volElem[jj].nDOFsSol; ++k) {
        const unsigned long kk = volElem[jj].offsetDOFsSolLocal + k;
        sendBuf[i].push_back(DOFsGlobalID[kk]);
      }
    }

    /* Send the data using non-blocking sends to avoid deadlock. */
    int dest = ranksSend[i];
    SU2_MPI::Isend(sendBuf[i].data(), sendBuf[i].size(), MPI_UNSIGNED_LONG,
                   dest, dest, MPI_COMM_WORLD, &sendReqs[i]);
  }

  /* Create a map of the receive rank to the index in ranksRecv. */
  map<int,int> rankToIndRecvBuf;
  for(int i=0; i<(int) ranksRecv.size(); ++i)
    rankToIndRecvBuf[ranksRecv[i]] = i;

  /* Loop over the number of ranks from which I receive data. */
  for(unsigned long i=0; i<ranksRecv.size(); ++i) {

    /* Block until a message arrives and determine the source and size
       of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for the receive buffer, receive the message
       and determine the actual index of this rank in ranksRecv. */
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank, MPI_COMM_WORLD, &status);

    map<int,int>::const_iterator MI = rankToIndRecvBuf.find(source);
    source = MI->second;

    /* Loop over the receive elements and copy the data for its DOFs in
       DOFsGlobalID. */
    unsigned long ii = 0;
    for(unsigned long j=0; j<elementsRecv[source].size(); ++j) {
      const unsigned long jj = elementsRecv[source][j];
      for(unsigned short k=0; k<volElem[jj].nDOFsSol; ++k, ++ii) {
        const unsigned long kk = volElem[jj].offsetDOFsSolLocal + k;
        DOFsGlobalID[kk] = recvBuf[ii];
      }
    }
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(ranksSend.size(), sendReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(MPI_COMM_WORLD);

#else

  /*--- Sequential implementation. Halo's may be present due to periodic
        boundary conditions. A loop over the receiving ranks is carried out
        to cover both situations, i.e. halo's present and halo's not present. ---*/
  for(unsigned long i=0; i<ranksRecv.size(); ++i) {
    for(unsigned long j=0; j<elementsRecv[i].size(); ++j) {
      const unsigned long elemR = elementsRecv[i][j];
      const unsigned long elemS = elementsSend[i][j];

      for(unsigned short k=0; k<volElem[elemR].nDOFsSol; ++k) {
        const unsigned long kR = volElem[elemR].offsetDOFsSolLocal + k;
        const unsigned long kS = volElem[elemS].offsetDOFsSolLocal + k;

        DOFsGlobalID[kR] = DOFsGlobalID[kS];
      }
    }
  }

#endif

  /*-------------------------------------------------------------------*/
  /* Step 3: Determine the entries of the graph for the locally stored */
  /*         data. Note that this data must also be determined for the */
  /*         halo DOFs, because the faces between partitions are       */
  /*         stored only once.                                         */
  /*********************************************************************/

  /* Allocate the memory for the first index of nonZeroEntriesJacobian.
     During the construction of this data, also the halo DOFs are needed,
     hence it is allocate for all DOFs stored on this rank. */
  nonZeroEntriesJacobian.resize(nDOFsLocTot);

  /* The residual of the DOFs of a volume element depends on all the
     DOFs of that volume element. Set these dependencies. Only needed
     for the owned elements. */
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
      const unsigned long jj = volElem[i].offsetDOFsSolLocal + j;
      for(unsigned short k=0; k<volElem[i].nDOFsSol; ++k) {
        const unsigned long kk = volElem[i].offsetDOFsSolLocal + k;
        nonZeroEntriesJacobian[jj].push_back(DOFsGlobalID[kk]);
      }
    }
  }

  /* Loop over the internal matching faces to set the dependencies of the
     DOFs w.r.t. the DOFs of the neighbors. It is assumed that the residual
     of the DOFs depends on all the DOFs of the neighboring elements. This
     is always correct for the Navier-Stokes equations. For the Euler equations
     in combination with a lumped mass matrix, this is a conservative
     estimate. However, it is also correct when the full mass matrix is used.
     Therefore, this assumption is used here. */
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
  for(unsigned long i=0; i<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++i) {

    /* Easier storage of the elements on both sides. */
    const unsigned long elem0 = matchingInternalFaces[i].elemID0;
    const unsigned long elem1 = matchingInternalFaces[i].elemID1;

    /* Double loop over the DOFs of both elements and set the dependencies
       accordingly. */
    for(unsigned short j=0; j<volElem[elem0].nDOFsSol; ++j) {
      const unsigned long jj = volElem[elem0].offsetDOFsSolLocal + j;
      for(unsigned short k=0; k<volElem[elem1].nDOFsSol; ++k) {
        const unsigned long kk = volElem[elem1].offsetDOFsSolLocal + k;

        nonZeroEntriesJacobian[jj].push_back(DOFsGlobalID[kk]);
        nonZeroEntriesJacobian[kk].push_back(DOFsGlobalID[jj]);
      }
    }
  }

  /*-------------------------------------------------------------------*/
  /* Step 4: Carry out the communication for the halo DOFs, such that  */
  /*         all the graph information is obtained for the owned DOFs. */
  /*********************************************************************/

#ifdef HAVE_MPI

  /* Parallel implementation. Allocate the memory for the inverse send buffers
     and loop over the number of ranks to which I have to send data for the
     inverse communication pattern. Note that self communication is not
     excluded here, because efficiency is not really important for this function. */
  vector<vector<unsigned long> > invSendBuf;
  invSendBuf.resize(ranksRecv.size());

  vector<SU2_MPI::Request> invSendReqs(ranksRecv.size());

  for(unsigned i=0; i<ranksRecv.size(); ++i) {

    /* Fill the inverse send buffer for this rank. */
    for(unsigned long j=0; j<elementsRecv[i].size(); ++j) {
      const unsigned long jj = elementsRecv[i][j];
      for(unsigned short k=0; k<volElem[jj].nDOFsSol; ++k) {
        const unsigned long kk = volElem[jj].offsetDOFsSolLocal + k;
        invSendBuf[i].push_back(nonZeroEntriesJacobian[kk].size());
        invSendBuf[i].insert(invSendBuf[i].end(),
                             nonZeroEntriesJacobian[kk].begin(),
                             nonZeroEntriesJacobian[kk].end());
      }
    }

    /* Send the data using non-blocking sends to avoid deadlock. */
    int dest = ranksRecv[i];
    SU2_MPI::Isend(invSendBuf[i].data(), invSendBuf[i].size(), MPI_UNSIGNED_LONG,
                   dest, dest+1, MPI_COMM_WORLD, &invSendReqs[i]);
  }

  /* Create a map of the inverse receive (i.e. the original send) rank
     to the index in ranksSend. */
  map<int,int> rankToIndSendBuf;
  for(int i=0; i<(int) ranksSend.size(); ++i)
    rankToIndSendBuf[ranksSend[i]] = i;

  /* Loop over the number of ranks from which I receive data
     in the inverse communication. */
  for(unsigned long i=0; i<ranksSend.size(); ++i) {

    /* Block until a message arrives and determine the source and size
       of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, MPI_COMM_WORLD, &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for the receive buffer, receive the message
       and determine the actual index of this rank in ranksSend. */
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank+1, MPI_COMM_WORLD, &status);

    map<int,int>::const_iterator MI = rankToIndSendBuf.find(source);
    source = MI->second;

    /* Loop over the elements for which the data was just received. */
    int ii = 0;
    for(unsigned long j=0; j<elementsSend[source].size(); ++j) {
      const unsigned long jj = elementsSend[source][j];

      for(unsigned short k=0; k<volElem[jj].nDOFsSol; ++k) {
        const unsigned long kk       = volElem[jj].offsetDOFsSolLocal + k;
        const unsigned long nEntries = recvBuf[ii++];

        nonZeroEntriesJacobian[kk].insert(nonZeroEntriesJacobian[kk].end(),
                                          recvBuf.data() + ii,
                                          recvBuf.data() + ii + nEntries);
        ii += nEntries;
      }
    }
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(ranksRecv.size(), invSendReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(MPI_COMM_WORLD);

#else
  /*--- Sequential implementation. Just add the data of the halo DOFs
        to the corresponding owned DOFs. Note that only when periodic
        boundary conditions are present, something will happen. ---*/
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    for(unsigned long j=0; j<elementsSend[i].size(); ++j) {
      const unsigned long elemR = elementsRecv[i][j];
      const unsigned long elemS = elementsSend[i][j];

      for(unsigned short k=0; k<volElem[elemR].nDOFsSol; ++k) {
        const unsigned long kR = volElem[elemR].offsetDOFsSolLocal + k;
        const unsigned long kS = volElem[elemS].offsetDOFsSolLocal + k;

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

  /* Loop over the locally owned DOFs and remove the possible double
     entries. This is accomplished by first sorting the data followed
     by the removal of duplicate entries using the unique function. */
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
    sort(nonZeroEntriesJacobian[i].begin(), nonZeroEntriesJacobian[i].end());
    vector<unsigned long>::iterator lastEntry;
    lastEntry = unique(nonZeroEntriesJacobian[i].begin(), nonZeroEntriesJacobian[i].end());
    nonZeroEntriesJacobian[i].erase(lastEntry, nonZeroEntriesJacobian[i].end());
  }

  /* Delete the memory of the halo data. To make sure that the memory
     is physically deleted, use the swap function. */
  for(unsigned long i=nDOFsLocOwned; i<nDOFsLocTot; ++i)
    vector<unsigned long>().swap(nonZeroEntriesJacobian[i]);
}

void CFEM_DG_EulerSolver::MetaDataJacobianComputation(const CMeshFEM    *FEMGeometry,
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

  /* Get the send and receive ranks for FEMGeometry. */
  const vector<int> &ranksSend = FEMGeometry->GetRanksSend();
  const vector<int> &ranksRecv = FEMGeometry->GetRanksRecv();

  /* Store the ranks this rank communicates with in a map. These are both
     sending and receiving ranks. Make sure to exclude my own rank. */
  map<int,int> rankCommToInd;
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    if(ranksSend[i] != rank) {
      const int ind = rankCommToInd.size();
      rankCommToInd[ranksSend[i]] = ind;
    }
  }

  for(unsigned long i=0; i<ranksRecv.size(); ++i) {
    map<int,int>::const_iterator MI = rankCommToInd.find(ranksRecv[i]);
    if((MI == rankCommToInd.end()) &&(ranksRecv[i] != rank)) {
      const int ind = rankCommToInd.size();
      rankCommToInd[ranksRecv[i]] = ind;
    }
  }

  /* Allocate the memory for the first index of the send buffers. */
  vector<vector<unsigned long> > sendBuf(rankCommToInd.size(), vector<unsigned long>(0));

  /* Define the map from the global matrix index to the color and the set to
     keep track whether an index has already been treated. */
  map<unsigned long, int> mapMatrixIndToColor;
  set<unsigned long> setTreatedIndices;

  /* Loop over the owned DOFs and its matrix entries to create the contents of
     either mapMatrixIndToColor (owned DOFs) or sendBuf (unowned DOFs). */
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
    for(unsigned long j=0; j<nonZeroEntriesJacobian[i].size(); ++j) {
      const unsigned long jj = nonZeroEntriesJacobian[i][j];

      /* Check if jj has not been stored yet. */
      pair<set<unsigned long>::iterator,bool> attempt = setTreatedIndices.insert(jj);
      if( attempt.second ) {

        /* Determine the rank where this DOF is stored. */
        vector<unsigned long>::iterator low;
        low = lower_bound(nDOFsPerRank.begin(), nDOFsPerRank.end(), jj);
        int rankDOF = low - nDOFsPerRank.begin();
        if(*low > jj) --rankDOF;

        /* If rankDOF is the current rank, create the entry in
           mapMatrixIndToColor. Otherwise store jj in the appropriate send
           buffer. */
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
                   dest, dest+2, MPI_COMM_WORLD, &sendReqs[i]);
  }

  /* Loop over the ranks from which I receive data to be processed. The number
     of ranks is equal to the number of ranks to which I just sent data. */
  vector<SU2_MPI::Request> sendReturnReqs(rankCommToInd.size());
  vector<vector<int> > sendReturnBuf(rankCommToInd.size(), vector<int>(0));
  for(unsigned long i=0; i<rankCommToInd.size(); ++i) {

    /* Block until a message arrives and determine the source and size
       of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+2, MPI_COMM_WORLD, &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for the receive buffer as well as for the
       return send buffer. Receive the message afterwards. */
    vector<unsigned long> recvBuf(sizeMess);
    sendReturnBuf[i].resize(sizeMess);

    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank+2, MPI_COMM_WORLD, &status);

    /* Loop over the data just received and fill the return send buffer
       with the color of the DOFs. */
    for(int j=0; j<sizeMess; ++j) {
      const unsigned long jj = recvBuf[j] - nDOFsPerRank[rank];
      if(jj >= nDOFsLocOwned) {
        cout << "This DOF should be owned, but it is not. This should not happen." << endl;
        exit(1);
      }
      sendReturnBuf[i][j] = colorLocalDOFs[jj];
    }

    /* Send the return buffer back to the calling rank. Again use non-blocking
       sends to avoid deadlock. */
    SU2_MPI::Isend(sendReturnBuf[i].data(), sendReturnBuf[i].size(), MPI_INT,
                   source, source+3, MPI_COMM_WORLD, &sendReturnReqs[i]);
  }

  /* Complete the first round of non-blocking sends. */
  SU2_MPI::Waitall(sendReqs.size(), sendReqs.data(), MPI_STATUSES_IGNORE);

  /* Loop over the ranks from which I receive the return information. */
  for(unsigned long i=0; i<rankCommToInd.size(); ++i) {

    /* Block until a message arrives and determine the source of the message
       and its index in the original send buffers. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+3, MPI_COMM_WORLD, &status);
    int source = status.MPI_SOURCE;

    MI = rankCommToInd.find(source);
    const int ind = MI->second;

    /* Allocate the memory for the receive buffer and receive the message using
       a blocking receive. */
    vector<int> recvBuf(sendBuf[ind].size());
    SU2_MPI::Recv(recvBuf.data(), recvBuf.size(), MPI_INT,
                  source, rank+3, MPI_COMM_WORLD, &status);

    /* Loop over the data just received and add them to the map
       mapMatrixIndToColor .*/
    for(unsigned long j=0; j<sendBuf[ind].size(); ++j)
      mapMatrixIndToColor[sendBuf[ind][j]] = recvBuf[j];
  }

  /* Complete the second round of non-blocking sends. */
  SU2_MPI::Waitall(sendReturnReqs.size(), sendReturnReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(MPI_COMM_WORLD);

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Part 3: For each locally owned DOF it is determined which matrix   ---*/
  /*---         entry is computed for which color. This information is     ---*/
  /*---         stored in the vector of vectors colorToIndEntriesJacobian. ---*/
  /*--------------------------------------------------------------------------*/

  /* Loop over the ownded DOFs. */
  colorToIndEntriesJacobian.resize(nDOFsLocOwned);
  for(unsigned long i=0; i<nDOFsLocOwned; ++i) {

    /* Initialize the elements of colorToIndEntriesJacobian to -1 to indicate
       that the color does not have a matrix entity . */
    colorToIndEntriesJacobian[i].assign(nGlobalColors, -1);

    /* Loop over the matrix entries of this DOF. */
    for(int j=0; j<(int) nonZeroEntriesJacobian[i].size(); ++j) {

      /* Search for the entry in mapMatrixIndToColor and determine the
         color of this entry. */
      map<unsigned long, int>::const_iterator MMI;
      MMI = mapMatrixIndToColor.find(nonZeroEntriesJacobian[i][j]);
      if(MMI == mapMatrixIndToColor.end()) {
        cout << "Matrix entry not found in mapMatrixIndToColor." << endl;
        cout << "This should not happen." << endl;
        exit(1);
      }
      const int color = MMI->second;

      /* Set the correct information in colorToIndEntriesJacobian. */
      colorToIndEntriesJacobian[i][color] = j;
    }
  }
}

void CFEM_DG_EulerSolver::SetUpTaskList(CConfig *config) {

  /* Check whether an ADER space-time step must be carried out.
     When only a spatial Jacobian this is false per definition.  */
  if( (config->GetKind_TimeIntScheme_Flow() == ADER_DG) &&
     !(config->GetJacobian_Spatial_Discretization_Only()) ) {

    /*------------------------------------------------------------------------*/
    /* ADER time integration with local time stepping. There are 2^(M-1)      */
    /* subtime steps, where M is the number of time levels. For each of       */
    /* the subtime steps a number of tasks must be carried out, hence the     */
    /* total list can become rather lengthy.                                  */
    /*------------------------------------------------------------------------*/

    /* Determine the number of subtime steps and abbreviate the number of
       time integration points a bit easier.  */
    const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
    const unsigned int nSubTimeSteps = pow(2,nTimeLevels-1);

    const unsigned short nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();

    /* Define the two dimensional vector to store the latest index for a
       certain task for every time level. */
    vector<vector<int> > indexInList(CTaskDefinition::ADER_UPDATE_SOLUTION+1,
                                     vector<int>(nTimeLevels, -1));

    /* Loop over the number of subtime steps in the algorithm. */
    for(unsigned int step=0; step<nSubTimeSteps; ++step) {

      /* Determine the time level for which an update must be carried out
         for this step. */
      unsigned int ii = step+1;
      unsigned short timeLevel = 0;
      while( !(ii%2) ) {ii/=2; ++timeLevel;}

      /* Definition of the variable to store a previous indices of
         tasks that must have been completed. */
      int prevInd[4];

      /* Carry out the predictor step of the communication elements of level 0
         if these elements are present on this rank. */
      unsigned long elemBeg = nVolElemOwnedPerTimeLevel[0]
                            + nVolElemInternalPerTimeLevel[0];
      unsigned long elemEnd = nVolElemOwnedPerTimeLevel[1];
      if(elemEnd > elemBeg) {
        prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][0];
        indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0] = tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS, 0, prevInd[0]));
      }

      /* Initiate the communication of elements of level 0, if there is
         something to be communicated. */
#ifdef HAVE_MPI
      if( commRequests[0].size() ) {
        if(elemEnd > elemBeg)   // Data needs to be computed before sending.
          prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0];
        else                    // Data only needs to be received.
          prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][0];

        /* Make sure that any previous communication has been completed. */
        prevInd[1] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][0];

        /* Create the task. */
        indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][0] = tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION, 0,
                                            prevInd[0], prevInd[1]));
      }
#endif

      /* Check if the time level is not nTimeLevels-1. In that case carry out
         the predictor step of the communication elements for the next time
         level and initiate their communication, if appropriate on this rank. */
      if(timeLevel < (nTimeLevels-1)) {
        const unsigned short nL = timeLevel+1;

        /* Carry out the predictor step of the communication elements of level nL
           if these elements are present on this rank. */
        elemBeg = nVolElemOwnedPerTimeLevel[nL] + nVolElemInternalPerTimeLevel[nL];
        elemEnd = nVolElemOwnedPerTimeLevel[nL+1];

        if(elemEnd > elemBeg) {
          prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][nL];
          indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL] = tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS, nL, prevInd[0]));
        }

        /* Initiate the communication of elements of level nL, if there is
           something to be communicated. */
#ifdef HAVE_MPI
        if( commRequests[nL].size() ) {
          if(elemEnd > elemBeg)   // Data needs to be computed before sending.
            prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL];
          else                    // Data only needs to be received.
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][nL];

          /* Make sure that any previous communication has been completed. */
          prevInd[1] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][nL];

          /* Create the actual task. */
          indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][nL] = tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION, nL,
                                              prevInd[0], prevInd[1]));
        }
#endif
      }

      /* Carry out the predictor step of the internal elements of time level 0,
         if these elements are present on this rank. */
      elemBeg = nVolElemOwnedPerTimeLevel[0];
      elemEnd = nVolElemOwnedPerTimeLevel[0] + nVolElemInternalPerTimeLevel[0];

      if(elemEnd > elemBeg) {
        prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][0];
        indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][0] = tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS, 0, prevInd[0]));
      }

      /* Determine the tasks to be completed before the communication of time
         level 0 can be completed. */
      prevInd[0] = prevInd[1] = prevInd[2] = -1;
#ifdef HAVE_MPI
      if( commRequests[0].size() )
        prevInd[0] = indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][0];
#endif

      if( elementsSendSelfComm[0].size() ) {
        prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0];
        prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][0];
      }

      /* Make sure that the -1 in prevInd, if any, are numbered last. */
      sort(prevInd, prevInd+3, greater<int>());

      /* Complete the communication of time level 0, if there is something
         to be completed. */
      if(prevInd[0] > -1) {
        indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][0] = tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION, 0,
                                            prevInd[0], prevInd[1], prevInd[2]));
      }

      /* Check if the time level is not nTimeLevels-1. In that case carry out
         the predictor step of the internal elements for the next time
         level and complete the communication for this time level,
         if appropriate on this rank. */
      if(timeLevel < (nTimeLevels-1)) {
        const unsigned short nL = timeLevel+1;

        /* Carry out the predictor step of the internal elements of level nL
           if these elements are present on this rank. */
        elemBeg = nVolElemOwnedPerTimeLevel[nL];
        elemEnd = nVolElemOwnedPerTimeLevel[nL] + nVolElemInternalPerTimeLevel[nL];

        if(elemEnd > elemBeg) {
          prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][nL];
          indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][nL] = tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS, nL, prevInd[0]));
        }

        /* Determine the tasks to be completed before the communication of time
           level nL can be completed. */
        prevInd[0] = prevInd[1] = prevInd[2] = -1;
#ifdef HAVE_MPI
        if( commRequests[nL].size() )
          prevInd[0] = indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][nL];
#endif

        if( elementsSendSelfComm[nL].size() ) {
          prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL];
          prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][nL];
        }

        /* Make sure that the -1 in prevInd, if any, are numbered last. */
        sort(prevInd, prevInd+3, greater<int>());

        /* Complete the communication of time level nL, if there is something
           to be completed. */
        if(prevInd[0] > -1) {
          indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][nL] = tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION, nL,
                                              prevInd[0], prevInd[1], prevInd[2]));
        }
      }

      /* Loop over the time integration points to compute the space time
         integral for the corrector step. */
      for(unsigned short intPoint=0; intPoint<nTimeIntegrationPoints; ++intPoint) {

        /* Loop over the time levels to be treated. */
        for(unsigned short level=0; level<=timeLevel; ++level) {

          /* The solution of the owned elements of the current time level and
             the adjacent owned elements of the next time level must be
             interpolated to the integration point intPoint. Check if these
             elements are present. */
          unsigned long nAdjOwnedElem = 0;
          if(level < (nTimeLevels-1))
            nAdjOwnedElem = ownedElemAdjLowTimeLevel[level+1].size();

          unsigned long nOwnedElem = nVolElemOwnedPerTimeLevel[level+1]
                                   - nVolElemOwnedPerTimeLevel[level];
          if(nAdjOwnedElem || nOwnedElem) {

            /* Determine the tasks that should have been completed before this
               task can be carried out. This depends on the time integration
               point. */
            if(intPoint == 0) {

              /* First time integration point. The predictor steps of the
                 owned elements must have been completed before this task
                 can be carried out. */
              if( nOwnedElem ) {
                prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][level];
                prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][level];
              }
              else {
                prevInd[0] = prevInd[1] = -1;
              }

              if( nAdjOwnedElem ) {
                prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][level+1];
                prevInd[3] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][level+1];
              }
              else {
                prevInd[2] = prevInd[3] = -1;
              }
            }
            else {

              /* Not the first integration point. The residual computations of
                 the previous integration point must have been completed, because
                 the solution in the work vectors is overwritten. */
              prevInd[0] = indexInList[CTaskDefinition::VOLUME_RESIDUAL][level];
              prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS][level];
              prevInd[2] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
              prevInd[3] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level];
            }

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2], prevInd[3]));

            /* The info on the integration point and whether or not this
               time integration corresponds to the second part for the
               adjacent elements must be added to the task just created. */
            tasksList.back().intPointADER          = intPoint;
            tasksList.back().secondPartTimeIntADER = level < timeLevel;

            /* If artificial viscosity is used for the shock capturing
               terms, these terms must be computed for the owned elements,
               including the adjacent ones. */
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS][level];
            indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS,
                                                level, prevInd[0]));
          }

          /* The solution of the halo elements of the current time level and
             the adjacent halo elements of the next time level must be
             interpolated to the integration point intPoint. Check if these
             elements are present. */
          unsigned long nAdjHaloElem = 0;
          if(level < (nTimeLevels-1))
            nAdjHaloElem = haloElemAdjLowTimeLevel[level+1].size();

          unsigned long nHaloElem = nVolElemHaloPerTimeLevel[level+1]
                                  - nVolElemHaloPerTimeLevel[level];
          if(nAdjHaloElem || nHaloElem) {

            /* Determine the tasks that should have been completed before this
               task can be carried out. This depends on the time integration
               point. */
            if(intPoint == 0) {

              /* First time integration point. The communication of the
                 halo elements must have been completed before this task
                 can be carried out. */
             if( nHaloElem )
               prevInd[0] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level];
             else
               prevInd[0] = -1;

             if( nAdjHaloElem )
               prevInd[1] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level+1];
             else
               prevInd[1] = -1;
            }
            else {

              /* Not the first integration point. The residual computation of
                 the previous integration point must have been completed, because
                 the solution in the work vectors is overwritten. */
              prevInd[0] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
              prevInd[1] = -1;
            }

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS,
                                                level, prevInd[0], prevInd[1]));

            /* The info on the integration point and whether or not this
               time integration corresponds to the second part for the
               adjacent elements must be added to the task just created. */
            tasksList.back().intPointADER          = intPoint;
            tasksList.back().secondPartTimeIntADER = level < timeLevel;

            /* If artificial viscosity is used for the shock capturing
               terms, these terms must be computed for the halo elements,
               including the adjacent ones. */
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][level];
            indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS,
                                                level, prevInd[0]));
          }

          /* Compute the surface residuals for this time level that involve
             halo elements, if present. Afterwards, accumulate the space time
             residuals for the halo elements. */
          if(nMatchingInternalFacesWithHaloElem[level+1] >
             nMatchingInternalFacesWithHaloElem[level]) {

            /* Create the dependencies for the surface residual part that involve
               halo elements. For all but the first integration point, make sure
               that the previous residual is already accumulated, because this
               task will overwrite that residual. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            prevInd[1] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level];

            if(intPoint == 0)
              prevInd[2] = -1;
            else
              prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level];

            /* Create the task for the surface residual. */
            indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2]));

            /* Create the task to accumulate the surface residuals of the halo
               elements. Make sure to set the integration point for this task. */
            prevInd[0] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
            indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS,
                                                level, prevInd[0]));
            tasksList.back().intPointADER = intPoint;
          }

          /* If this is the last integration point, initiate the reverse
             communication for the residuals of this time level, if there
             is something to communicate. */
          if(intPoint == (nTimeIntegrationPoints-1)) {

            /* Check if there is something to communicate. Despite the fact
               that self communication takes place when the reverse communication
               is completed, a check for self communication is still needed,
               because of the periodic transformations. */
            bool commData = false;

#ifdef HAVE_MPI
            if( commRequests[level].size() ) commData = true;
#endif
            if(commData || elementsSendSelfComm[level].size()) {

              /* Determine the dependencies before the reverse communication
                 can start. */
              prevInd[0] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level];
              if( haloElemAdjLowTimeLevel[level].size() )
                prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level-1];
              else
                prevInd[1] = -1;

              prevInd[2] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][level];

              /* Create the task. */
              indexInList[CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION][level] = tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION,
                                                  level, prevInd[0], prevInd[1], prevInd[2]));
            }
          }

          /* Compute the contribution of the volume integral, boundary conditions and
             surface integral between owned elements for this time level, if present. */
          if( nOwnedElem ) {

            /* Create the dependencies for this task. The shock capturing viscosity of
               the owned elements must be completed and for all integration points but
               the first, the computed residuals must have been accumulated in the space
               time residual, because the tasks below will overwrite the residuals. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            if(intPoint == 0)
              prevInd[1] = -1;
            else
              prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];

            /* Create the tasks. */
            indexInList[CTaskDefinition::VOLUME_RESIDUAL][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::VOLUME_RESIDUAL, level,
                                                prevInd[0], prevInd[1]));

            indexInList[CTaskDefinition::BOUNDARY_CONDITIONS][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS, level,
                                                prevInd[0], prevInd[1]));

            if(nMatchingInternalFacesLocalElem[level+1] >
               nMatchingInternalFacesLocalElem[level]) {
              indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level] = tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS,
                                                  level, prevInd[0], prevInd[1]));
            }
          }

          /* Accumulate the space time residuals for the owned elements of this
             level, if these elements are present. */
          if(nAdjOwnedElem || nOwnedElem) {

            /* Create the dependencies for this task. */
            prevInd[0] = indexInList[CTaskDefinition::VOLUME_RESIDUAL][level];
            prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS][level];
            prevInd[2] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level];
            prevInd[3] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level] = tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2], prevInd[3]));
            tasksList.back().intPointADER = intPoint;
          }

          /* If this is the last integration point, complete the reverse
             communication for the residuals of this time level, if there
             is something to communicate. */
          if(intPoint == (nTimeIntegrationPoints-1)) {

            bool commData = false;

#ifdef HAVE_MPI
            if( commRequests[level].size() ) commData = true;
#endif
            if(commData || elementsSendSelfComm[level].size()) {

              /* Determine the dependencies before the reverse communication
                 can be completed. */
              prevInd[0] = indexInList[CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION][level];
              prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];
              if( ownedElemAdjLowTimeLevel[level].size() )
                prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level-1];
              else
                prevInd[2] = -1;

              /* Create the task. */
              indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][level] = tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION,
                                                  level, prevInd[0], prevInd[1], prevInd[2]));
            }
          }

        } /* End loop time level. */

      } /* End loop time integration points. */

      /* Loop again over the number of active time levels for this subtime
         step and compute the update for the DOFs of the owned elements. */
      for(unsigned short level=0; level<=timeLevel; ++level) {
        if(nVolElemOwnedPerTimeLevel[level+1] > nVolElemOwnedPerTimeLevel[level]) {

          /* Updates must be carried out for this time level. First multiply the
             residuals by the inverse of the mass matrix. This task can be
             completed if the communication of this level has been completed.
             However, it is possible that no communication is needed and
             therefore the accumulation of the space time residuals is also
             added to the dependency list. */
          prevInd[0] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][level];
          prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];
          if( ownedElemAdjLowTimeLevel[level].size() )
            prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level-1];
          else
            prevInd[2] = -1;

          indexInList[CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX][level] = tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX,
                                               level, prevInd[0], prevInd[1], prevInd[2]));

          /* Compute the new state vector for this time level. */
          prevInd[0] = indexInList[CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX][level];
          indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][level] = tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_UPDATE_SOLUTION,
                                               level, prevInd[0]));
        }
      }

    } /* End loop subtime steps. */


    // EXTRA FOR DEBUGGING
/*  for(int i=0; i<size; ++i) {
      if(i == rank) {
        cout << endl;
        cout << "Task list for rank " << rank << endl;
        cout << "------------------------------------------------" << endl;
        cout << "Number of tasks: " << tasksList.size() << endl;
        for(unsigned long j=0; j<tasksList.size(); ++j) {

          cout << "Task " << j << ": ";
          switch( tasksList[j].task ) {
            case CTaskDefinition::NO_TASK: cout << "NO_TASK" << endl; break;
            case CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS: cout << "ADER_PREDICTOR_STEP_COMM_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS: cout << "ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS" << endl; break;
            case CTaskDefinition::INITIATE_MPI_COMMUNICATION: cout << "INITIATE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::COMPLETE_MPI_COMMUNICATION: cout << "COMPLETE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION: cout << "INITIATE_REVERSE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION: cout << "COMPLETE_REVERSE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS: cout << "ADER_TIME_INTERPOLATE_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS: cout << "ADER_TIME_INTERPOLATE_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS: cout << "SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS: cout << "SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::VOLUME_RESIDUAL: cout << "VOLUME_RESIDUAL" << endl; break;
            case CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS: cout << "SURFACE_RESIDUAL_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS: cout << "SURFACE_RESIDUAL_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::BOUNDARY_CONDITIONS: cout << "BOUNDARY_CONDITIONS" << endl; break;
            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS: cout << "SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS: cout << "SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS: cout << "ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS: cout << "ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX: cout << "MULTIPLY_INVERSE_MASS_MATRIX" << endl; break;
            case CTaskDefinition::ADER_UPDATE_SOLUTION: cout << "ADER_UPDATE_SOLUTION" << endl; break;
            default: cout << "This cannot happen" << endl;
          }
          cout << " Time level: " << tasksList[j].timeLevel
               << " Integration point: " << tasksList[j].intPointADER
               << " Second part: " << tasksList[j].secondPartTimeIntADER << endl;
          cout << " Depends on tasks:";
          for(unsigned short k=0; k<tasksList[j].nIndMustBeCompleted; ++k)
            cout << " " << tasksList[j].indMustBeCompleted[k];
          cout << endl << endl;
        }

      }

#ifdef HAVE_MPI
      SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
    } */

//  cout << "CFEM_DG_EulerSolver::SetUpTaskList: ADER tasklist printed";
//  exit(1);

    // END EXTRA FOR DEBUGGING.
  }
  else {

    /*------------------------------------------------------------------------*/
    /* Standard time integration scheme for which the spatial residual must   */
    /* be computed for the DOFS of the owned elements. This results in a      */
    /* relatively short tasks list, which can be set easily.                  */
    /*------------------------------------------------------------------------*/

    tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION,                   0));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS,     0));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::VOLUME_RESIDUAL,                              0,  1));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION,                   0,  0));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS,      0,  3));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS,               0,  1, 4));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS,  0,  5));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION,           0,  6));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS,              0,  1));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS,                          0,  1));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION,           0,  2, 7));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS, 0,  2, 8, 9, 10));
    tasksList.push_back(CTaskDefinition(CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX,                 0, 11));
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
  /*--- Step 1. Create the triple vector, which contains elements to be    ---*/
  /*---         sent and received per time level. Note that when time      ---*/
  /*---         accurate local time stepping is not used, this is just a   ---*/
  /*---         copy of elementsSend and elementsRecv.                     ---*/
  /*--------------------------------------------------------------------------*/

  /* Allocate the first two indices of the triple vectors. */
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

  vector<vector<vector<unsigned long> > > elemSendPerTimeLevel, elemRecvPerTimeLevel;
  elemSendPerTimeLevel.resize(nTimeLevels);
  elemRecvPerTimeLevel.resize(nTimeLevels);

  for(unsigned short i=0; i<nTimeLevels; ++i) {
    elemSendPerTimeLevel[i].resize(elementsSend.size());
    elemRecvPerTimeLevel[i].resize(elementsRecv.size());
  }

  /* Loop over the send data and set the corresponding data
     in elemSendPerTimeLevel. */
  for(unsigned long i=0; i<elementsSend.size(); ++i) {
    for(unsigned long j=0; j<elementsSend[i].size(); ++j) {
      const unsigned long ii = elementsSend[i][j];
      elemSendPerTimeLevel[volElem[ii].timeLevel][i].push_back(ii);
    }
  }

  /* Loop over the receive data and set the corresponding data
     in recvSendPerTimeLevel. */
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

  /* Allocate the first index of elementsSendSelfComm and elementsRecvSelfComm. */
  elementsSendSelfComm.resize(nTimeLevels);
  elementsRecvSelfComm.resize(nTimeLevels);

  /* Loop over the ranks to which data is sent and copy the elements for self
     communication for all time levels. */
  for(unsigned long i=0; i<ranksSend.size(); ++i) {
    if(ranksSend[i] == rank) {
      for(unsigned short j=0; j<nTimeLevels; ++j)
        elementsSendSelfComm[j] = elemSendPerTimeLevel[j][i];
    }
  }

  /* Loop over the ranks from which data is received and copy the elements
     for self communication for all time levels. */
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

  /* Determine the number of time DOFs that must be communicated.
     For non-ADER schemes this is simply set to 1. */
  unsigned short nTimeDOFs = 1;
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG)
    nTimeDOFs = config->GetnTimeDOFsADER_DG();

  /* Determine the number of items per DOF that must be communicated. */
  const unsigned short nItemsPerDOF = nTimeDOFs*nVar;

  /* Allocate the memory for the first index of the vectors that
     determine the MPI communication patterns. */
  commRequests.resize(nTimeLevels);
  elementsRecvMPIComm.resize(nTimeLevels);
  elementsSendMPIComm.resize(nTimeLevels);
  ranksRecvMPI.resize(nTimeLevels);
  ranksSendMPI.resize(nTimeLevels);
  commRecvBuf.resize(nTimeLevels);
  commSendBuf.resize(nTimeLevels);

  /* Loop over the time levels. */
  for(unsigned short level=0; level<nTimeLevels; ++level) {

    /* Determine the number of ranks from which data will be received
       for this time level. Self communication is excluded. */
    int nRankRecv = 0;
    for(unsigned long i=0; i<ranksRecv.size(); ++i) {
      if((ranksRecv[i] != rank) && elemRecvPerTimeLevel[level][i].size()) ++nRankRecv;
    }

    /* Determine the number of ranks to which data will be send
       for this time level. Self communication is excluded. */
    int nRankSend = 0;
    for(unsigned long i=0; i<ranksSend.size(); ++i) {
      if((ranksSend[i] != rank) && elemSendPerTimeLevel[level][i].size()) ++nRankSend;
    }

    /* Allocate the memory for the second index of the vectors that
       determine the MPI communication. */
    commRequests[level].resize(nRankRecv+nRankSend);
    elementsRecvMPIComm[level].resize(nRankRecv);
    elementsSendMPIComm[level].resize(nRankSend);
    ranksRecvMPI[level].resize(nRankRecv);
    ranksSendMPI[level].resize(nRankSend);
    commRecvBuf[level].resize(nRankRecv);
    commSendBuf[level].resize(nRankSend);

    /* Determine the receive information. */
    nRankRecv = 0;
    for(unsigned long i=0; i<ranksRecv.size(); ++i) {
      if((ranksRecv[i] != rank) && elemRecvPerTimeLevel[level][i].size()) {

        /* Copy the elements to be received and the rank from where they come. */
        elementsRecvMPIComm[level][nRankRecv] = elemRecvPerTimeLevel[level][i];
        ranksRecvMPI[level][nRankRecv] = ranksRecv[i];

        /* Determine the size of the receive buffer and allocate the memory. */
        unsigned long sizeBuf = 0;
        for(unsigned long j=0; j<elementsRecvMPIComm[level][nRankRecv].size(); ++j) {
          const unsigned long jj = elementsRecvMPIComm[level][nRankRecv][j];
          sizeBuf += volElem[jj].nDOFsSol;
        }

        sizeBuf *= nItemsPerDOF;
        commRecvBuf[level][nRankRecv].resize(sizeBuf);

        /* Update nRankRecv. */
        ++nRankRecv;
      }
    }

    /* Determine the send information. */
    nRankSend = 0;
    for(unsigned long i=0; i<ranksSend.size(); ++i) {
      if((ranksSend[i] != rank) && elemSendPerTimeLevel[level][i].size()) {

        /* Copy the elements to be sent and the rank to where they are sent. */
        elementsSendMPIComm[level][nRankSend] = elemSendPerTimeLevel[level][i];
        ranksSendMPI[level][nRankSend] = ranksSend[i];

        /* Determine the size of the send buffer and allocate the memory. */
        unsigned long sizeBuf = 0;
        for(unsigned long j=0; j<elementsSendMPIComm[level][nRankSend].size(); ++j) {
          const unsigned long jj = elementsSendMPIComm[level][nRankSend][j];
          sizeBuf += volElem[jj].nDOFsSol;
        }

        sizeBuf *= nItemsPerDOF;
        commSendBuf[level][nRankSend].resize(sizeBuf);

        /* Update nRankSend. */
        ++nRankSend;
      }
    }
  }

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 4. Store the information for the rotational periodic          ---*/
  /*---         corrections for the halo elements.                         ---*/
  /*--------------------------------------------------------------------------*/

  /* Get the data for the rotational periodic halos from FEMGeometry. */
  const vector<unsigned short> &markersRotPer = FEMGeometry->GetRotPerMarkers();
  const vector<vector<unsigned long> > &rotPerHalos = FEMGeometry->GetRotPerHalos();

  /* Allocate the memory for the first and second index of
     halosRotationalPeriodicity. Also allocate the memory to store the
     rotation matrices of the periodic transformations. */
  const unsigned short nRotPerMarkers = markersRotPer.size();
  halosRotationalPeriodicity.resize(nTimeLevels);
  for(unsigned short i=0; i<nTimeLevels; ++i)
    halosRotationalPeriodicity[i].resize(nRotPerMarkers);

  rotationMatricesPeriodicity.resize(9*nRotPerMarkers);

  /* Loop over the rotational periodic transformations. */
  unsigned int ii = 0;
  for(unsigned short i=0; i<nRotPerMarkers; ++i) {

   /* Get the rotation angles from config for this marker. */
   const unsigned short pInd = markersRotPer[i];
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

    /* Loop over the elements of this periodic transformation and store them
       in the appropriate location of halosRotationalPeriodicity. */
    for(unsigned long j=0; j<rotPerHalos[i].size(); ++j) {
      const unsigned long jj = rotPerHalos[i][j];
      halosRotationalPeriodicity[volElem[jj].timeLevel][i].push_back(jj);
    }
  }
}

void CFEM_DG_EulerSolver::Initiate_MPI_Communication(CConfig *config,
                                                     const unsigned short timeLevel) {
#ifdef HAVE_MPI

  /* Check if there is anything to communicate. */
  if( commRequests[timeLevel].size() ) {

    /* Set the pointer to the memory, whose data must be communicated.
       This depends on the time integration scheme used. For ADER the data of
       the predictor part is communicated, while for the other time integration
       schemes typically the working solution is communicated. Note that if the
       working solution is communicated, there is only one time level. */
    unsigned short nTimeDOFs;
    su2double *commData;
    if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {
      nTimeDOFs = config->GetnTimeDOFsADER_DG();
      commData  = VecSolDOFsPredictorADER.data();
    }
    else {
      nTimeDOFs = 1;
      commData  = VecWorkSolDOFs[0].data();
    }

    /* Loop over the number of ranks to which this rank must send data. */
    int indComm = 0;
    for(unsigned long i=0; i<ranksSendMPI[timeLevel].size(); ++i, ++indComm) {
      unsigned long ii = 0;

      /* Loop over the elements to be sent and copy the solution data
         into the send buffer. */
      su2double *sendBuf = commSendBuf[timeLevel][i].data();
      for(unsigned long j=0; j<elementsSendMPIComm[timeLevel][i].size(); ++j) {
        const unsigned long jj = elementsSendMPIComm[timeLevel][i][j];
        const unsigned long nItems = volElem[jj].nDOFsSol * nVar;
        const unsigned long nBytes = nItems * sizeof(su2double);

        for(unsigned short k=0; k<nTimeDOFs; ++k) {
          const unsigned long indS = nVar*(volElem[jj].offsetDOFsSolLocal + k*nDOFsLocTot);
          memcpy(sendBuf+ii, commData+indS, nBytes);
          ii += nItems;
        }
      }

      /* Send the data using non-blocking sends. */
      int dest = ranksSendMPI[timeLevel][i];
      int tag  = dest + timeLevel;
      SU2_MPI::Isend(sendBuf, ii, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD,
                     &commRequests[timeLevel][indComm]);
    }

    /* Loop over the number of ranks from which data is received. */
    for(unsigned long i=0; i<ranksRecvMPI[timeLevel].size(); ++i, ++indComm) {

      /* Post the non-blocking receive. */
      int source = ranksRecvMPI[timeLevel][i];
      int tag    = rank + timeLevel;
      SU2_MPI::Irecv(commRecvBuf[timeLevel][i].data(),
                     commRecvBuf[timeLevel][i].size(),
                     MPI_DOUBLE, source, tag, MPI_COMM_WORLD,
                     &commRequests[timeLevel][indComm]);
    }
  }

#endif
}

bool CFEM_DG_EulerSolver::Complete_MPI_Communication(CConfig *config,
                                                     const unsigned short timeLevel,
                                                     const bool commMustBeCompleted) {

  /* Set the pointer to the memory, whose data must be communicated.
     This depends on the time integration scheme used. For ADER the data of
     the predictor part is communicated, while for the other time integration
     schemes typically the working solution is communicated. Note that if the
     working solution is communicated, there is only one time level. */
  unsigned short nTimeDOFs;
  su2double *commData;
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {
    nTimeDOFs = config->GetnTimeDOFsADER_DG();
    commData  = VecSolDOFsPredictorADER.data();
  }
  else {
    nTimeDOFs = 1;
    commData  = VecWorkSolDOFs[0].data();
  }

  /*-----------------------------------------------------------------------*/
  /*--- Complete the MPI communication, if needed and if possible, and  ---*/
  /*--- copy the data from the receive buffers into the correct         ---*/
  /*--- location in commData.                                           ---*/
  /*-----------------------------------------------------------------------*/

#ifdef HAVE_MPI
  if( commRequests[timeLevel].size() ) {

    /*--- There are communication requests to be completed. Check if these
          requests must be completed. In that case Waitall is used.
          Otherwise, Testall is used to check if all the requests have
          been completed. If not, false is returned. ---*/
    if( commMustBeCompleted ) {
      SU2_MPI::Waitall(commRequests[timeLevel].size(),
                       commRequests[timeLevel].data(), MPI_STATUSES_IGNORE);
    }
    else {
      int flag;
      SU2_MPI::Testall(commRequests[timeLevel].size(),
                       commRequests[timeLevel].data(), &flag, MPI_STATUSES_IGNORE);
      if( !flag ) return false;
    }

    /* Loop over the number of ranks from which this rank has received data. */
    for(unsigned long i=0; i<ranksRecvMPI[timeLevel].size(); ++i) {
      unsigned long ii = 0;

      /* Loop over the elements to be received and copy the solution data
         into commData. */
      su2double *recvBuf = commRecvBuf[timeLevel][i].data();
      for(unsigned long j=0; j<elementsRecvMPIComm[timeLevel][i].size(); ++j) {
        const unsigned long jj = elementsRecvMPIComm[timeLevel][i][j];
        const unsigned long nItems = volElem[jj].nDOFsSol * nVar;
        const unsigned long nBytes = nItems * sizeof(su2double);

        for(unsigned short k=0; k<nTimeDOFs; ++k) {
          const unsigned long indR = nVar*(volElem[jj].offsetDOFsSolLocal + k*nDOFsLocTot);
          memcpy(commData+indR, recvBuf+ii, nBytes);
          ii += nItems;
        }
      }
    }
  }
#endif

  /*-----------------------------------------------------------------------*/
  /*---               Carry out the self communication.                 ---*/
  /*-----------------------------------------------------------------------*/

  /* Loop over the number of elements involved and copy the data of the DOFs. */
  for(unsigned long i=0; i<elementsSendSelfComm[timeLevel].size(); ++i) {
    const unsigned long elemS  = elementsSendSelfComm[timeLevel][i];
    const unsigned long elemR  = elementsRecvSelfComm[timeLevel][i];
    const unsigned long nBytes = volElem[elemS].nDOFsSol * nVar*sizeof(su2double);

    for(unsigned short j=0; j<nTimeDOFs; ++j) {
      const unsigned long indS = nVar*(volElem[elemS].offsetDOFsSolLocal + j*nDOFsLocTot);
      const unsigned long indR = nVar*(volElem[elemR].offsetDOFsSolLocal + j*nDOFsLocTot);

      memcpy(commData+indR, commData+indS, nBytes);
    }
  }

  /*------------------------------------------------------------------------*/
  /*--- Correct the vector quantities in the rotational periodic halo's. ---*/
  /*------------------------------------------------------------------------*/

  /*--- Loop over the markers for which a rotational periodic
        correction must be applied to the momentum variables. ---*/
  unsigned int ii = 0;
  const unsigned long nRotPerMarkers = halosRotationalPeriodicity[timeLevel].size();
  for(unsigned long k=0; k<nRotPerMarkers; ++k) {

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
    const unsigned long nHaloElem = halosRotationalPeriodicity[timeLevel][k].size();
    for(unsigned long j=0; j<nHaloElem; ++j) {

      /* Easier storage of the halo index and loop over its DOFs, including
         the multiple time DOFs for ADER. */
      for(unsigned short tInd=0; tInd<nTimeDOFs; ++tInd) {
        const unsigned long tIndOff = tInd*nDOFsLocTot;
        const unsigned long ind     = halosRotationalPeriodicity[timeLevel][k][j];
        for(unsigned short i=0; i<volElem[ind].nDOFsSol; ++i) {

          /* Determine the pointer in commData where the solution of this DOF
             is stored and copy the momentum variables. Note that a rotational
             correction can only take place for a 3D simulation. */
          su2double *sol = commData + nVar*(volElem[ind].offsetDOFsSolLocal + i + tIndOff);

          const su2double ru = sol[1], rv = sol[2], rw = sol[3];

          /* Correct the momentum variables. */
          sol[1] = rotMatrix[0][0]*ru + rotMatrix[0][1]*rv + rotMatrix[0][2]*rw;
          sol[2] = rotMatrix[1][0]*ru + rotMatrix[1][1]*rv + rotMatrix[1][2]*rw;
          sol[3] = rotMatrix[2][0]*ru + rotMatrix[2][1]*rv + rotMatrix[2][2]*rw;
        }
      }
    }
  }

  /* Return true to indicate that the communication has been completed. */
  return true;
}

void CFEM_DG_EulerSolver::Initiate_MPI_ReverseCommunication(CConfig *config,
                                                            const unsigned short timeLevel) {

  /* Set the pointer to the residual to be communicated. */
  su2double *resComm;
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) resComm = VecTotResDOFsADER.data();
  else                                                resComm = VecResDOFs.data();

  /*------------------------------------------------------------------------*/
  /*--- Correct the vector residuals in the rotational periodic halo's.  ---*/
  /*------------------------------------------------------------------------*/

  /*--- Loop over the markers for which a rotational periodic
        correction must be applied to the momentum variables. ---*/
  unsigned int ii = 0;
  const unsigned short nRotPerMarkers = halosRotationalPeriodicity[0].size();
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
    const unsigned long nHaloElem = halosRotationalPeriodicity[timeLevel][k].size();
    for(unsigned long j=0; j<nHaloElem; ++j) {

      /* Easier storage of the halo index and loop over its DOFs. */
      const unsigned long ind = halosRotationalPeriodicity[timeLevel][k][j];
      for(unsigned short i=0; i<volElem[ind].nDOFsSol; ++i) {

        /* Determine the pointer in VecResDOFs where the residual of this DOF
           is stored and copy the momentum residuals. Note that a rotational
           correction can only take place for a 3D simulation. */
        su2double *res = resComm + nVar*(volElem[ind].offsetDOFsSolLocal + i);

        const su2double ru = res[1], rv = res[2], rw = res[3];

        /* Correct the momentum variables. */
        res[1] = rotMatrix[0][0]*ru + rotMatrix[0][1]*rv + rotMatrix[0][2]*rw;
        res[2] = rotMatrix[1][0]*ru + rotMatrix[1][1]*rv + rotMatrix[1][2]*rw;
        res[3] = rotMatrix[2][0]*ru + rotMatrix[2][1]*rv + rotMatrix[2][2]*rw;
      }
    }
  }

#ifdef HAVE_MPI

  /* Check if there is anything to communicate. */
  if( commRequests[timeLevel].size() ) {

    /* Loop over the number of ranks from which this rank receives data in
       the original communication pattern. In the reverse pattern, data
       has to be sent. */
    int indComm = 0;
    for(unsigned long i=0; i<ranksRecvMPI[timeLevel].size(); ++i, ++indComm) {
      unsigned long ii = 0;

      /* Loop over the elements copy the residual data into recvBuf. */
      su2double *recvBuf = commRecvBuf[timeLevel][i].data();
      for(unsigned long j=0; j<elementsRecvMPIComm[timeLevel][i].size(); ++j) {
        const unsigned long jj = elementsRecvMPIComm[timeLevel][i][j];
        const unsigned long nItems = volElem[jj].nDOFsSol * nVar;
        const unsigned long nBytes = nItems * sizeof(su2double);

        const unsigned long indS = nVar*volElem[jj].offsetDOFsSolLocal;
        memcpy(recvBuf+ii, resComm+indS, nBytes);
        ii += nItems;
      }

      /* Send the data using non-blocking sends. */
      int dest = ranksRecvMPI[timeLevel][i];
      int tag  = dest + timeLevel + 20;
      SU2_MPI::Isend(recvBuf, ii, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD,
                     &commRequests[timeLevel][indComm]);
    }

    /* Post the non-blocking receives. As this is the reverse communication,
       a loop over the sending ranks must be carried out. */
    for(unsigned long i=0; i<ranksSendMPI[timeLevel].size(); ++i, ++indComm) {

      int source = ranksSendMPI[timeLevel][i];
      int tag    = rank + timeLevel + 20;
      SU2_MPI::Irecv(commSendBuf[timeLevel][i].data(),
                     commSendBuf[timeLevel][i].size(),
                     MPI_DOUBLE, source, tag, MPI_COMM_WORLD,
                     &commRequests[timeLevel][indComm]);
    }
  }

#endif

}

bool CFEM_DG_EulerSolver::Complete_MPI_ReverseCommunication(CConfig *config,
                                                            const unsigned short timeLevel,
                                                            const bool commMustBeCompleted) {
  /* Set the pointer to the residual to be communicated. */
  su2double *resComm;
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) resComm = VecTotResDOFsADER.data();
  else                                                resComm = VecResDOFs.data();

#ifdef HAVE_MPI

  /*-----------------------------------------------------------------------*/
  /*---   Complete the MPI communication, if needed and if possible.    ---*/
  /*-----------------------------------------------------------------------*/

  /* Check if there are any requests to complete for this time level. */
  if( commRequests[timeLevel].size() ) {

    /*--- There are communication requests to be completed. Check if these
          requests must be completed. In that case Waitall is used.
          Otherwise, Testall is used to check if all the requests have
          been completed. If not, return false. ---*/
    if( commMustBeCompleted ) {
      SU2_MPI::Waitall(commRequests[timeLevel].size(),
                       commRequests[timeLevel].data(), MPI_STATUSES_IGNORE);
    }
    else {
      int flag;
      SU2_MPI::Testall(commRequests[timeLevel].size(),
                       commRequests[timeLevel].data(), &flag, MPI_STATUSES_IGNORE);
      if( !flag ) return false;
    }

    /*-------------------------------------------------------------------------*/
    /*---    Update the residuals of the owned DOFs with the data received. ---*/
    /*-------------------------------------------------------------------------*/

    /*--- Loop over the received residual data from all ranks and update the
          residual of the DOFs of the corresponding elements. Note that in
          reverse mode the send communication data must be used. ---*/
    for(unsigned long i=0; i<ranksSendMPI[timeLevel].size(); ++i) {

      /* Loop over the elements that must be updated. */
      unsigned long nn = 0;
      for(unsigned long j=0; j<elementsSendMPIComm[timeLevel][i].size(); ++j) {
        const unsigned long jj = elementsSendMPIComm[timeLevel][i][j];

        /* Easier storage of the starting residual of this element and the
           data in the buffer and update the residuals of the DOFs. */
        const unsigned short nRes = nVar*volElem[jj].nDOFsSol;
        su2double            *res = resComm + nVar*volElem[jj].offsetDOFsSolLocal;
        const su2double      *buf = commSendBuf[timeLevel][i].data() + nn;

        nn += nRes;
        for(unsigned short k=0; k<nRes; ++k)
          res[k] += buf[k];
      }
    }
  }

#endif

  /*-----------------------------------------------------------------------*/
  /*---               Carry out the self communication.                 ---*/
  /*-----------------------------------------------------------------------*/

  /*--- Loop over the number of elements involved and update the residuals
        of the owned DOFs. ---*/
  for(unsigned long i=0; i<elementsSendSelfComm[timeLevel].size(); ++i) {

    const unsigned long volOwned = elementsSendSelfComm[timeLevel][i];
    const unsigned long volHalo  = elementsRecvSelfComm[timeLevel][i];
    su2double *resOwned = resComm + nVar*volElem[volOwned].offsetDOFsSolLocal;
    su2double *resHalo  = resComm + nVar*volElem[volHalo].offsetDOFsSolLocal;

    const unsigned short nRes = nVar*volElem[volOwned].nDOFsSol;
    for(unsigned short j=0; j<nRes; ++j)
      resOwned[j] += resHalo[j];
  }

  /* Initialize the halo residuals of the just completed communication pattern
     to zero if ADER-DG is used. */
  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {
    const unsigned long iBeg      = nVar*volElem[nVolElemHaloPerTimeLevel[timeLevel]].offsetDOFsSolLocal;
    const unsigned long volIndEnd = nVolElemHaloPerTimeLevel[timeLevel+1]-1;
    const unsigned long iEnd      = nVar*(volElem[volIndEnd].offsetDOFsSolLocal
                                  +       volElem[volIndEnd].nDOFsSol);

    for(unsigned long i=iBeg; i<iEnd; ++i)
      VecTotResDOFsADER[i] = 0.0;
  }

  /* Return true to indicate that the communication has been completed. */
  return true;
}

void CFEM_DG_EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

#ifdef INVISCID_VORTEX

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

}

void CFEM_DG_EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long ErrorCounter = 0;

  /*-----------------------------------------------------------------------------*/
  /*--- Check for non-physical points. Only needed for a compressible solver. ---*/
  /*-----------------------------------------------------------------------------*/

  if(config->GetKind_Regime() == COMPRESSIBLE) {

    /*--- Make a distinction between 2D and 3D for optimal performance. ---*/
    switch( nDim ) {

      case 2: {
        /*--- 2D simulation. Loop over the owned DOFs and check for
              non-physical solutions. If found, the state is overwritten with
              the free-stream values. This usually does not work. ---*/
        for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
          su2double *solDOF = VecSolDOFs.data() + nVar*i;
          const su2double DensityInv   = 1.0/solDOF[0];
          const su2double Mom2         = solDOF[1]*solDOF[1] + solDOF[2]*solDOF[2];
          const su2double StaticEnergy = DensityInv*(solDOF[3] - 0.5*DensityInv*Mom2);

          FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
          const su2double Pressure    = FluidModel->GetPressure();
          const su2double Temperature = FluidModel->GetTemperature();

          if((Pressure < 0.0) || (solDOF[0] < 0.0) || (Temperature < 0.0)) {
            ++ErrorCounter;

            solDOF[0] = ConsVarFreeStream[0]; solDOF[1] = ConsVarFreeStream[1];
            solDOF[2] = ConsVarFreeStream[2]; solDOF[3] = ConsVarFreeStream[3];
          }
        }

        break;
      }

      case 3: {
        /*--- 3D simulation. Loop over the owned DOFs and check for
              non-physical solutions. If found, the state is overwritten with
              the free-stream values. This usually does not work. ---*/
        for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
          su2double *solDOF = VecSolDOFs.data() + nVar*i;
          const su2double DensityInv   = 1.0/solDOF[0];
          const su2double Mom2         = solDOF[1]*solDOF[1] + solDOF[2]*solDOF[2]
                                       + solDOF[3]*solDOF[3];
          const su2double StaticEnergy = DensityInv*(solDOF[4] - 0.5*DensityInv*Mom2);

          FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
          const su2double Pressure    = FluidModel->GetPressure();
          const su2double Temperature = FluidModel->GetTemperature();

          if((Pressure < 0.0) || (solDOF[0] < 0.0) || (Temperature < 0.0)) {
            ++ErrorCounter;

            solDOF[0] = ConsVarFreeStream[0]; solDOF[1] = ConsVarFreeStream[1];
            solDOF[2] = ConsVarFreeStream[2]; solDOF[3] = ConsVarFreeStream[3];
            solDOF[4] = ConsVarFreeStream[4];
          }
        }

        break;
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

  /*-----------------------------------------------------------------------------*/
  /*                       Check for grid motion.                                */
  /*-----------------------------------------------------------------------------*/

  const bool harmonic_balance = config->GetUnsteady_Simulation() == HARMONIC_BALANCE;
  if(config->GetGrid_Movement() && !harmonic_balance) {

    /*--- Determine the type of grid motion. ---*/
    switch( config->GetKind_GridMovement(0) ) {

      case RIGID_MOTION: {

        /* Rigid body motion described. At the moment this is only
           possible for the Classical Runge Kutta scheme. */
        if(config->GetKind_TimeIntScheme() != CLASSICAL_RK4_EXPLICIT)
          SU2_MPI::Error("Rigid body motion only possible for CLASSICAL_RK4_EXPLICIT.",
                         CURRENT_FUNCTION);

        /* Determine whether or not it is needed to compute the motion data. */
        const unsigned long ExtIter = config->GetExtIter();

        bool computeMotion = false, firstTime = false;
        if(ExtIter == 0 && iStep == 0) computeMotion = firstTime = true;
        if(iStep == 1 || iStep == 3)   computeMotion = true;

        if( computeMotion ) {

          /* Determine the number of time levels. For the classical Runge-Kutta
             scheme this should be 1. */
          const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

          /* Determine the time for which the motion data must be determined. */
          const su2double deltaT = config->GetDelta_UnstTimeND();

          su2double tNew = ExtIter*deltaT;
          if( iStep ) tNew += 0.25*(iStep+1)*deltaT;

          /* Determine the time for which the currently stored position was
             calculated. */
          const su2double tOld = tNew - 0.5*deltaT;

          /* Hard code the plunging motion in y-direction. */
          const su2double b2New    = 0.25*tNew*tNew*(3.0-tNew);
          const su2double db2Newdt = 0.75*tNew*(2.0-tNew);

          const su2double b2Old = 0.25*tOld*tOld*(3.0-tOld);
          const su2double dB2   = b2New - b2Old;

          /*-------------------------------------------------------------------*/
          /*--- Update the coordinates of the DOFs (both grid and solution) ---*/
          /*--- and of the integration points (volume and surfaces).        ---*/
          /*-------------------------------------------------------------------*/

          /* The coordinates must be updated if this is not the first
             first stage of the first time step. */
          if( !firstTime ) {

            /* Adapt the coordinates of the grid points. */
            for(unsigned long l=0; l<nMeshPoints; ++l)
              meshPoints[l].coor[1] += dB2;

            /* Loop over the owned volume elements. */
            for(unsigned long l=0; l<nVolElemOwned; ++l) {

              /* Get the required data from the corresponding standard element. */
              const unsigned short ind  = volElem[l].indStandardElement;
              const unsigned short nInt = standardElementsSol[ind].GetNIntegration();

              /* Adapt the coordinates of the volume integration points. */
              for(unsigned short i=0; i<nInt; ++i) {
                su2double *coor = volElem[l].coorIntegrationPoints.data() + i*nDim;
                coor[1] += dB2;
              }

              /* Adapt the coordinates of the solution DOFs. */
              for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {
                su2double *coor = volElem[l].coorSolDOFs.data() + i*nDim;
                coor[1] += dB2;
              }
            }

            /* Loop over the internal matching faces. */
            for(unsigned long l=0; l<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++l) {

              /* Get the required information from the standard element. */
              const unsigned short ind  = matchingInternalFaces[l].indStandardElement;
              const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();

              /* Adapt the coordinates of the surface integration points. */
              for(unsigned short i=0; i<nInt; ++i) {
                su2double *coor = matchingInternalFaces[l].coorIntegrationPoints.data() + i*nDim;
                coor[1] += dB2;
              }
            }

            /* The physical boundary faces. Exclude the periodic boundaries,
               because these are not physical boundaries. */
            for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

              if( !(boundaries[iMarker].periodicBoundary) ) {

                /* Easier storage of the variables for this boundary. */
                const unsigned long nSurfElem = boundaries[iMarker].surfElem.size();
                CSurfaceElementFEM  *surfElem = boundaries[iMarker].surfElem.data();

                /*--- Loop over the surface elements and update the coordinates of
                      the integration points. ---*/
                for(unsigned long l=0; l<nSurfElem; ++l) {
                  const unsigned short ind  = surfElem[l].indStandardElement;
                  const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

                  for(unsigned short i=0; i<nInt; ++i) {
                    su2double *coor = surfElem[l].coorIntegrationPoints.data() + i*nDim;
                    coor[1] += dB2;
                  }
                }
              }
            }
          }

          /*-------------------------------------------------------------------*/
          /*--- Compute the grid velocities of the volume solution DOFs and ---*/
          /*--- and of the integration points (volume and surfaces).        ---*/
          /*-------------------------------------------------------------------*/

          /* Loop over the owned volume elements. */
          for(unsigned long l=0; l<nVolElemOwned; ++l) {

            /* Get the required data from the corresponding standard element. */
            const unsigned short ind  = volElem[l].indStandardElement;
            const unsigned short nInt = standardElementsSol[ind].GetNIntegration();

            /* Compute the grid velocity in the integration points. */
            for(unsigned short i=0; i<nInt; ++i) {
              su2double *gridVel = volElem[l].gridVelocities.data() + i*nDim;
              gridVel[1] = db2Newdt;
            }

            /* Compute the grid velocities for the solution DOFs. */
            for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {
              su2double *gridVel = volElem[l].gridVelocitiesSolDOFs.data() + i*nDim;
              gridVel[1] = db2Newdt;
            }
          }

          /* Loop over the internal matching faces. */
          for(unsigned long l=0; l<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++l) {

            /* Get the required information from the standard element. */
            const unsigned short ind  = matchingInternalFaces[l].indStandardElement;
            const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();

            /* Compute the grid velocity in the integration points. */
            for(unsigned short i=0; i<nInt; ++i) {
              su2double *gridVel = matchingInternalFaces[l].gridVelocities.data() + i*nDim;
              gridVel[1] = db2Newdt;
            }
          }

          /* The physical boundary faces. Exclude the periodic boundaries,
             because these are not physical boundaries. */
          for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

            if( !(boundaries[iMarker].periodicBoundary) ) {

              /* Easier storage of the variables for this boundary. */
              const unsigned long nSurfElem = boundaries[iMarker].surfElem.size();
              CSurfaceElementFEM  *surfElem = boundaries[iMarker].surfElem.data();

              /*--- Loop over the surface elements and compute the grid velocities
                    of the integration points. ---*/
              for(unsigned long l=0; l<nSurfElem; ++l) {
                const unsigned short ind  = surfElem[l].indStandardElement;
                const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

                for(unsigned short i=0; i<nInt; ++i) {
                  su2double *gridVel = surfElem[l].gridVelocities.data() + i*nDim;
                  gridVel[1] = db2Newdt;
                }
              }
            }
          }
        }

        break;
      }

      default: {
        SU2_MPI::Error("Only rigid body motion possible for DG-FEM solver at the moment.",
                       CURRENT_FUNCTION);
      }
    }
  }
}

void CFEM_DG_EulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iMesh) { }

void CFEM_DG_EulerSolver::ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                                                 CNumerics **numerics, CConfig *config,
                                                 unsigned short iMesh, unsigned short RunTime_EqSystem) {

  /* Write a message that the Jacobian is being computed. */
  if(rank == MASTER_NODE) {
    cout << endl;
    cout << "Computing the Jacobian via coloring." << endl;
    cout << endl << flush;
  }

  /* Easier storage of the number of variables squared. */
  const unsigned short nVar2 = nVar*nVar;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Computation of the actual Jacobian by looping over the     ---*/
  /*---         of colors and disturbing the appropriate DOFs for each     ---*/
  /*---         color. Note that for each color an additional loop over    ---*/
  /*---         the number of variables is needed to create the Jacobian.  ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the number of non-zeros in the Jacobian matrix for the owned
     DOFs in cumulative storage format. */
  vector<unsigned long> nNonZeroEntries(nDOFsLocOwned+1);

  nNonZeroEntries[0] = 0;
  for(unsigned i=0; i<nDOFsLocOwned; ++i)
    nNonZeroEntries[i+1] = nNonZeroEntries[i] + nonZeroEntriesJacobian[i].size();

  /* Copy the solution into the working variables. */
  Set_OldSolution(geometry);

  /* Allocate the memory for local part of the Jacobian. Note that passivedouble
     must be used for the Jacobian matrix. */
  vector<passivedouble> Jacobian(nVar2*nNonZeroEntries[nDOFsLocOwned]);

  /* Loop over the colors and write a message, if needed. */
  for(int color=0; color<nGlobalColors; ++color) {

    if(rank == MASTER_NODE && !(color%10)) {
      cout << "Color " << color+1 << " out of " << nGlobalColors
           << " is computed" << endl;
    }

    /* Loop over the number of variables. */
    for(short var=0; var<nVar; ++var) {

      /* Loop over the DOFs that must be disturbed for this color. */
      for(unsigned long j=0; j<localDOFsPerColor[color].size(); ++j) {
        const unsigned long jj = localDOFsPerColor[color][j];

        su2double *solDOF = VecWorkSolDOFs[0].data() + jj*nVar;

#ifdef CODI_FORWARD_TYPE
        solDOF[var].setGradient(1.0);
#else
        solDOF[var] += 0.001;   /* This is to avoid a compiler warning. */
#endif
      }

      /* Carry out all the tasks to compute the residual. */
      ProcessTaskList_DG(geometry, solver_container, numerics, config, iMesh);

      /* Loop over all the locally owned DOFs. */
      for(unsigned long i=0; i<nDOFsLocOwned; ++i) {

        /* Check if a Jacobian entry is computed for this color. */
        if(colorToIndEntriesJacobian[i][color] >= 0) {

          /* Set the pointer to entry in Jacobian where the derivative data
             must be stored. */
          const int ind = colorToIndEntriesJacobian[i][color];
          passivedouble *Jac = Jacobian.data() + nVar2*(nNonZeroEntries[i] + ind);

          /* Set the pointer to the residual of the current DOF. */
          const su2double *resDOF = VecResDOFs.data() + i*nVar;

          /* Store the matrix entries. */
          for(unsigned short j=0; j<nVar; ++j) {
#ifdef CODI_FORWARD_TYPE
            Jac[var+j*nVar] = resDOF[j].getGradient();
#else
            Jac[var+j*nVar] = resDOF[j];   /* This is to avoid a compiler warning. */
#endif
          }
        }
      }

      /* Reset the disturbances for this color. */
      for(unsigned long j=0; j<localDOFsPerColor[color].size(); ++j) {
        const unsigned long jj = localDOFsPerColor[color][j];

        su2double *solDOF = VecWorkSolDOFs[0].data() + jj*nVar;

#ifdef CODI_FORWARD_TYPE
        solDOF[var].setGradient(0.0);
#else
        solDOF[var] -= 0.001;   /* This is to avoid a compiler warning. */
#endif
      }

    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: The writing of the Jacobian in Block Compressed Row        ---*/
  /*---         Storage to file.                                           ---*/
  /*--------------------------------------------------------------------------*/

  /* Write a message that the Jacobian is written to file. */
  if(rank == MASTER_NODE) {
    cout << endl;
    cout << "Writing the Jacobian Block Compressed Row Storage." << flush;
  }

#ifdef HAVE_MPI
  /* Parallel mode. The parallel IO functionality is used to Write
     the Jacobian matrix. */
  cout << "CFEM_DG_EulerSolver::ComputeSpatialJacobian: Not implemented yet" << endl;
  exit(1);

#else
  /* Sequential mode. The entire file is written by one rank,
     Open the file for writing. */
  FILE *fJac = fopen("Jacobian_DG.bin", "wb");
  if( !fJac) {
    cout << "Could not open the file Jacobian_DG for binary writing." << endl;
    exit(1);
  }

  /* Write the number of DOFs (Block rows in the matrix), the number of
     variables (the dimension of each block) and the number of neighboring
     blocks per DOF in cumulative storage format. */
  fwrite(&nDOFsPerRank[1], 1, sizeof(unsigned long), fJac);
  fwrite(&nVar, 1, sizeof(unsigned short), fJac);
  fwrite(nNonZeroEntries.data(), nNonZeroEntries.size(),
         sizeof(unsigned long), fJac);

  /* Write the block column indices of the neighboring blocks. */
  for(unsigned long i=0; i<nDOFsLocOwned; ++i)
    fwrite(nonZeroEntriesJacobian[i].data(), nonZeroEntriesJacobian[i].size(),
           sizeof(unsigned long), fJac);

  /* Write the actual matrix elements. */
  fwrite(Jacobian.data(), Jacobian.size(), sizeof(passivedouble), fJac);

  /* Close the file again. */
  fclose(fJac);

#endif

  /* Write a message that the writing is done.. */
  if(rank == MASTER_NODE) {
    cout << " Done." << endl;
    cout << endl << flush;
  }
}

void CFEM_DG_EulerSolver::Set_OldSolution(CGeometry *geometry) {

  memcpy(VecWorkSolDOFs[0].data(), VecSolDOFs.data(), VecSolDOFs.size()*sizeof(su2double));
}

void CFEM_DG_EulerSolver::Set_NewSolution(CGeometry *geometry) {

  memcpy(VecSolDOFsNew.data(), VecSolDOFs.data(), VecSolDOFs.size()*sizeof(su2double));
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

  /*--- Explicit time stepping with imposed time step. If the unsteady CFL is
        set to zero (default), it uses the defined unsteady time step,
        otherwise it computes the time step based on the provided unsteady CFL.
        Note that the regular CFL option in the config is always ignored with
        time stepping. ---*/
  if (time_stepping && (config->GetUnst_CFL() == 0.0)) {

    /*--- Loop over the owned volume elements and set the fixed dt. ---*/
    for(unsigned long i=0; i<nVolElemOwned; ++i)
      VecDeltaTime[i] = config->GetDelta_UnstTimeND();

  } else {

    /*--- Check for a compressible solver. ---*/
    if(config->GetKind_Regime() == COMPRESSIBLE) {

      /*--- Loop over the owned volume elements. ---*/
      for(unsigned long i=0; i<nVolElemOwned; ++i) {

        su2double charVel2Max = 0.0;

        /*--- Make a distinction between 2D and 3D for optimal performance. ---*/
        switch( nDim ) {

          case 2: {
            /*--- 2D simulation. Loop over the DOFs of this element
                  and determine the maximum wave speed. ---*/
            for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
              const su2double *solDOF  = VecSolDOFs.data()
                                       + nVar*(volElem[i].offsetDOFsSolLocal + j);
              const su2double *gridVel = volElem[i].gridVelocitiesSolDOFs.data()
                                       + j*nDim;

              /* Compute the velocities and the internal energy per unit mass. */
              const su2double DensityInv   = 1.0/solDOF[0];
              const su2double u            = DensityInv*solDOF[1];
              const su2double v            = DensityInv*solDOF[2];
              const su2double StaticEnergy = DensityInv*solDOF[3] - 0.5*(u*u + v*v);

              /*--- Compute the maximum value of the wave speed. This is a rather
                    conservative estimate. ---*/
              FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
              const su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
              const su2double SoundSpeed  = sqrt(fabs(SoundSpeed2));

              const su2double radx     = fabs(u-gridVel[0]) + SoundSpeed;
              const su2double rady     = fabs(v-gridVel[1]) + SoundSpeed;
              const su2double charVel2 = radx*radx + rady*rady;

              charVel2Max = max(charVel2Max, charVel2);
            }

            break;
          }

          case 3: {
            /*--- 3D simulation. Loop over the DOFs of this element
                  and determine the maximum wave speed. ---*/
            for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
              const su2double *solDOF = VecSolDOFs.data()
                                      + nVar*(volElem[i].offsetDOFsSolLocal + j);
              const su2double *gridVel = volElem[i].gridVelocitiesSolDOFs.data()
                                       + j*nDim;

              /* Compute the velocities and the internal energy per unit mass. */
              const su2double DensityInv   = 1.0/solDOF[0];
              const su2double u            = DensityInv*solDOF[1];
              const su2double v            = DensityInv*solDOF[2];
              const su2double w            = DensityInv*solDOF[3];
              const su2double StaticEnergy = DensityInv*solDOF[4] - 0.5*(u*u + v*v + w*w);

              /*--- Compute the maximum value of the wave speed. This is a rather
                    conservative estimate. ---*/
              FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
              const su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
              const su2double SoundSpeed  = sqrt(fabs(SoundSpeed2));

              const su2double radx     = fabs(u-gridVel[0]) + SoundSpeed;
              const su2double rady     = fabs(v-gridVel[1]) + SoundSpeed;
              const su2double radz     = fabs(w-gridVel[2]) + SoundSpeed;
              const su2double charVel2 = radx*radx + rady*rady + radz*radz;

              charVel2Max = max(charVel2Max, charVel2);
            }

            break;
          }
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

      SU2_MPI::Error("Incompressible solver not implemented yet", CURRENT_FUNCTION);
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

void CFEM_DG_EulerSolver::ProcessTaskList_DG(CGeometry *geometry,  CSolver **solver_container,
                                             CNumerics **numerics, CConfig *config,
                                             unsigned short iMesh) {
  /* Variable for the internal timing and store the number of time levels.. */
  double tick = 0.0;
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

  /* Define and initialize the bool vector, that indicates whether or
     not the tasks from the list have been completed. */
  vector<bool> taskCompleted(tasksList.size(), false);

  /* While loop to carry out all the tasks in tasksList. */
  unsigned long lowestIndexInList = 0;
  while(lowestIndexInList < tasksList.size()) {

    /* Find the next task that can be carried out. The outer loop is there
       to make sure that a communication is completed in case there are no
       other tasks */
    for(unsigned short j=0; j<2; ++j) {
      bool taskCarriedOut = false;
      for(unsigned long i=lowestIndexInList; i<tasksList.size(); ++i) {

        /* Determine whether or not it can be attempted to carry out
           this task. */
        bool taskCanBeCarriedOut = !taskCompleted[i];
        for(unsigned short ind=0; ind<tasksList[i].nIndMustBeCompleted; ++ind) {
          if( !taskCompleted[tasksList[i].indMustBeCompleted[ind]] )
            taskCanBeCarriedOut = false;
        }

        if( taskCanBeCarriedOut ) {

          /*--- Determine the actual task to be carried out and do so. The
                only tasks that may fail are the completion of the non-blocking
                communication. If that is the case the next task needs to be
                found. ---*/
          switch( tasksList[i].task ) {

            case CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS: {

              /* Carry out the ADER predictor step for the elements whose
                 solution must be communicated for this time level. */
              const unsigned short level   = tasksList[i].timeLevel;
              const unsigned long  elemBeg = nVolElemOwnedPerTimeLevel[level]
                                           + nVolElemInternalPerTimeLevel[level];
              const unsigned long  elemEnd = nVolElemOwnedPerTimeLevel[level+1];

              config->Tick(&tick);
              ADER_DG_PredictorStep(config, elemBeg, elemEnd);
              config->Tock(tick,"ADER_DG_PredictorStep",2);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS: {

              /* Carry out the ADER predictor step for the elements whose
                 solution must not be communicated for this time level. */
              const unsigned short level   = tasksList[i].timeLevel;
              const unsigned long  elemBeg = nVolElemOwnedPerTimeLevel[level];
              const unsigned long  elemEnd = nVolElemOwnedPerTimeLevel[level]
                                           + nVolElemInternalPerTimeLevel[level];
              config->Tick(&tick);
              ADER_DG_PredictorStep(config, elemBeg, elemEnd);
              config->Tock(tick,"ADER_DG_PredictorStep",2);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::INITIATE_MPI_COMMUNICATION: {

              /* Start the MPI communication of the solution in the halo elements. */
              Initiate_MPI_Communication(config, tasksList[i].timeLevel);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::COMPLETE_MPI_COMMUNICATION: {

              /* Attempt to complete the MPI communication of the solution data.
                 For j==0, SU2_MPI::Testall will be used, which returns false if
                 not all requests can be completed. In that case the next task on
                 the list is carried out. If j==1, this means that the next
                 tasks are waiting for this communication to be completed and
                 hence MPI_Waitall is used. */
              if( Complete_MPI_Communication(config, tasksList[i].timeLevel,
                                             j==1) )
                taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION: {

              /* Start the communication of the residuals, for which the
                 reverse communication must be used. */
              Initiate_MPI_ReverseCommunication(config, tasksList[i].timeLevel);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION: {

              /* Attempt to complete the MPI communication of the residual data.
                 For j==0, SU2_MPI::Testall will be used, which returns false if
                 not all requests can be completed. In that case the next task on
                 the list is carried out. If j==1, this means that the next
                 tasks are waiting for this communication to be completed and
                 hence MPI_Waitall is used. */
              if( Complete_MPI_ReverseCommunication(config, tasksList[i].timeLevel,
                                                    j==1) )
                taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS: {

              /* Interpolate the predictor solution of the owned elements
                 in time to the given time integration point for the
                 given time level. */
              const unsigned short level = tasksList[i].timeLevel;
              unsigned long nAdjElem = 0, *adjElem = NULL;
              if(level < (nTimeLevels-1)) {
                nAdjElem = ownedElemAdjLowTimeLevel[level+1].size();
                adjElem  = ownedElemAdjLowTimeLevel[level+1].data();
              }

              config->Tick(&tick);
              ADER_DG_TimeInterpolatePredictorSol(config, tasksList[i].intPointADER,
                                                  nVolElemOwnedPerTimeLevel[level],
                                                  nVolElemOwnedPerTimeLevel[level+1],
                                                  nAdjElem, adjElem,
                                                  tasksList[i].secondPartTimeIntADER,
                                                  VecWorkSolDOFs[level].data());
              config->Tock(tick,"ADER_DG_TimeInterpolatePredictorSol",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS: {

              /* Interpolate the predictor solution of the halo elements
                 in time to the given time integration point for the
                 given time level. */
              const unsigned short level = tasksList[i].timeLevel;
              unsigned long nAdjElem = 0, *adjElem = NULL;
              if(level < (nTimeLevels-1)) {
                nAdjElem = haloElemAdjLowTimeLevel[level+1].size();
                adjElem  = haloElemAdjLowTimeLevel[level+1].data();
              }

              config->Tick(&tick);
              ADER_DG_TimeInterpolatePredictorSol(config, tasksList[i].intPointADER,
                                                  nVolElemHaloPerTimeLevel[level],
                                                  nVolElemHaloPerTimeLevel[level+1],
                                                  nAdjElem, adjElem,
                                                  tasksList[i].secondPartTimeIntADER,
                                                  VecWorkSolDOFs[level].data());
              config->Tock(tick,"ADER_DG_TimeInterpolatePredictorSol",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS: {

              /*--- Compute the artificial viscosity for shock capturing in DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              config->Tick(&tick);
              Shock_Capturing_DG(config, nVolElemOwnedPerTimeLevel[level],
                                 nVolElemOwnedPerTimeLevel[level+1]);
              config->Tock(tick,"Shock_Capturing",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS: {

              /*--- Compute the artificial viscosity for shock capturing in DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              config->Tick(&tick);
              Shock_Capturing_DG(config, nVolElemHaloPerTimeLevel[level],
                                 nVolElemHaloPerTimeLevel[level+1]);
              config->Tock(tick,"Shock_Capturing",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::VOLUME_RESIDUAL: {

              /*--- Compute the volume portion of the residual. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              config->Tick(&tick);
              Volume_Residual(config, nVolElemOwnedPerTimeLevel[level],
                              nVolElemOwnedPerTimeLevel[level+1]);
              config->Tock(tick,"Volume_Residual",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS: {

              /* Compute the residual of the faces that only involve owned elements. */
              const unsigned short level = tasksList[i].timeLevel;
              config->Tick(&tick);
              unsigned long indResFaces = startLocResInternalFacesLocalElem[level];
              ResidualFaces(config, nMatchingInternalFacesLocalElem[level],
                            nMatchingInternalFacesLocalElem[level+1],
                            indResFaces, numerics[CONV_TERM]);
              config->Tock(tick,"ResidualFaces",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS: {

              /* Compute the residual of the faces that involve a halo element. */
              const unsigned short level = tasksList[i].timeLevel;
              config->Tick(&tick);
              unsigned long indResFaces = startLocResInternalFacesWithHaloElem[level];
              ResidualFaces(config, nMatchingInternalFacesWithHaloElem[level],
                            nMatchingInternalFacesWithHaloElem[level+1],
                            indResFaces, numerics[CONV_TERM]);
              config->Tock(tick,"ResidualFaces",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::BOUNDARY_CONDITIONS: {

              /*--- Apply the boundary conditions ---*/
              config->Tick(&tick);
              Boundary_Conditions(tasksList[i].timeLevel, config, numerics);
              config->Tock(tick,"Boundary_Conditions",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS: {

              /* Create the final residual by summing up all contributions. */
              config->Tick(&tick);
              CreateFinalResidual(tasksList[i].timeLevel, true);
              config->Tock(tick,"CreateFinalResidual",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS: {

              /* Create the final residual by summing up all contributions. */
              config->Tick(&tick);
              CreateFinalResidual(tasksList[i].timeLevel, false);
              config->Tock(tick,"CreateFinalResidual",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS: {

              /* Accumulate the space time residuals for the owned elements
                 for ADER-DG. */
              config->Tick(&tick);
              AccumulateSpaceTimeResidualADEROwnedElem(config, tasksList[i].timeLevel,
                                                       tasksList[i].intPointADER);
              config->Tock(tick,"AccumulateSpaceTimeResidualADEROwnedElem",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS: {

              /* Accumulate the space time residuals for the halo elements
                 for ADER-DG. */
              config->Tick(&tick);
              AccumulateSpaceTimeResidualADERHaloElem(config, tasksList[i].timeLevel,
                                                      tasksList[i].intPointADER);
              config->Tock(tick,"AccumulateSpaceTimeResidualADERHaloElem",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX: {

              /*--- Multiply the residual by the (lumped) mass matrix, to obtain the final value. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              config->Tick(&tick);
              const bool useADER = config->GetKind_TimeIntScheme() == ADER_DG;
              MultiplyResidualByInverseMassMatrix(config, useADER,
                                                  nVolElemOwnedPerTimeLevel[level],
                                                  nVolElemOwnedPerTimeLevel[level+1]);
              config->Tock(tick,"MultiplyResidualByInverseMassMatrix",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_UPDATE_SOLUTION: {

              /*--- Perform the update step for ADER-DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              config->Tick(&tick);
              ADER_DG_Iteration(nVolElemOwnedPerTimeLevel[level],
                                nVolElemOwnedPerTimeLevel[level+1]);
              config->Tock(tick,"ADER_DG_Iteration",3);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            default: {

              cout << "Task not defined. This should not happen." << endl;
              exit(1);
            }
          }
        }

        /* Break the inner loop if a task has been carried out. */
        if( taskCarriedOut ) break;
      }

      /* Break the outer loop if a task has been carried out. */
      if( taskCarriedOut ) break;
    }

    /* Update the value of lowestIndexInList. */
    for(; lowestIndexInList < tasksList.size(); ++lowestIndexInList)
      if( !taskCompleted[lowestIndexInList] ) break;
  }
}

void CFEM_DG_EulerSolver::ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                                    CNumerics **numerics, CConfig *config,
                                                    unsigned short iMesh, unsigned short RunTime_EqSystem) {
  double tick = 0.0;

  /* Preprocessing. */
  config->Tick(&tick);
  Preprocessing(geometry, solver_container, config, iMesh, 0, RunTime_EqSystem, false);
  TolerancesADERPredictorStep();
  config->Tock(tick,"Preprocessing",2);

  /* Process the tasks list to carry out one ADER space time integration step. */
  config->Tick(&tick);
  ProcessTaskList_DG(geometry, solver_container, numerics, config, iMesh);
  config->Tock(tick,"ProcessTaskList_DG",3);

  /* Postprocessing. */
  config->Tick(&tick);
  Postprocessing(geometry, solver_container, config, iMesh);
  config->Tock(tick,"Postprocessing",2);
}

void CFEM_DG_EulerSolver::TolerancesADERPredictorStep(void) {

  /* Determine the maximum values of the conservative variables of the
     locally stored DOFs. Make a distinction between 2D and 3D for
     performance reasons. */
  su2double URef[] = {0.0, 0.0, 0.0, 0.0, 0.0};

  switch( nDim ) {
    case 2: {
      for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
        const su2double *solDOF = VecSolDOFs.data() + i*nVar;
        URef[0] = max(URef[0], fabs(solDOF[0]));
        URef[1] = max(URef[1], fabs(solDOF[1]));
        URef[2] = max(URef[2], fabs(solDOF[2]));
        URef[3] = max(URef[3], fabs(solDOF[3]));
      }
      break;
    }

    case 3: {
      for(unsigned long i=0; i<nDOFsLocOwned; ++i) {
        const su2double *solDOF = VecSolDOFs.data() + i*nVar;
        URef[0] = max(URef[0], fabs(solDOF[0]));
        URef[1] = max(URef[1], fabs(solDOF[1]));
        URef[2] = max(URef[2], fabs(solDOF[2]));
        URef[3] = max(URef[3], fabs(solDOF[3]));
        URef[4] = max(URef[3], fabs(solDOF[4]));
      }
      break;
    }
  }

  /* Determine the maximum values on all ranks. */
  TolSolADER.resize(nVar);

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(URef, TolSolADER.data(), nVar, MPI_DOUBLE, MPI_MAX,
                     MPI_COMM_WORLD);
#else
  for(unsigned short i=0; i<nVar; ++i) TolSolADER[i] = URef[i];
#endif

  /* Determine the maximum scale of the momentum variables and adapt the
     corresponding values of TolSolADER accordingly. */
  su2double momRef = TolSolADER[1];
  for(unsigned short i=2; i<=nDim; ++i) momRef = max(momRef, TolSolADER[i]);
  for(unsigned short i=1; i<=nDim; ++i) TolSolADER[i] = momRef;

  /* Currently the maximum values of the conserved variables are stored.
     Multiply by the relative tolerance to obtain the true tolerance values. */
  for(unsigned short i=0; i<nVar; ++i) TolSolADER[i] *= 1.e-6;
}

void CFEM_DG_EulerSolver::ADER_DG_PredictorStep(CConfig             *config,
                                                const unsigned long elemBeg,
                                                const unsigned long elemEnd) {

  /*----------------------------------------------------------------------*/
  /*---        Get the data of the ADER time integration scheme.       ---*/
  /*----------------------------------------------------------------------*/

  const unsigned short nTimeDOFs              = config->GetnTimeDOFsADER_DG();
  const unsigned short nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();
  const su2double     *timeIntegrationWeights = config->GetWeightsIntegrationADER_DG();
  const bool          useAliasedPredictor     = config->GetKind_ADER_Predictor() == ADER_ALIASED_PREDICTOR;

  /*--------------------------------------------------------------------------*/
  /*--- Loop over the given elemen range to compute the predictor solution.---*/
  /*--- For the predictor solution only an integration over the element is ---*/
  /*--- performed to obtain the weak formulation, i.e. no integration by   ---*/
  /*--- parts. As a consequence there is no surface term and also no       ---*/
  /*--- coupling with neighboring elements. This is also the reason why    ---*/
  /*--- the ADER scheme is only conditionally stable.                      ---*/
  /*--------------------------------------------------------------------------*/

  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* Easier storage of the number of spatial DOFs, the total
       number of variables for this element per time level and the
       number of space time DOFs. */
    const unsigned short nDOFs          = volElem[l].nDOFsSol;
    const unsigned short nVarNDOFs      = nVar*nDOFs;
    const unsigned short nDOFsNTimeDOFs = nDOFs*nTimeDOFs;

    /* Set the pointers for the working variables. */
    su2double *resInt = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    su2double *solPred = VecTmpMemory.data();
    su2double *solOld  = solPred + nVarNDOFs*nTimeDOFs;
    su2double *resSol  = solOld  + nVarNDOFs*nTimeDOFs;
    su2double *resTot  = resSol  + nVarNDOFs*nTimeDOFs;
    su2double *solInt  = resTot  + nVarNDOFs*nTimeDOFs;
    su2double *work    = solInt  + nVarNDOFs;

    /* Initialize the predictor solution to the current solution. */
    const su2double    *solCur = VecSolDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;
    const unsigned long nBytes = nVarNDOFs*sizeof(su2double);

    for(unsigned short j=0; j<nTimeDOFs; ++j) {
      su2double *solPredTimeInd = solPred + j*nVarNDOFs;
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

      su2double *resSolTimeInd = resSol + j*nVarNDOFs;
      for(unsigned short i=0; i<nVarNDOFs; ++i)
        resSolTimeInd[i] = resSol[i]*LagrangianBeginTimeIntervalADER_DG[j];
    }

    /* Multiply the values of resSol for the first time DOF with the */
    /* value of its corresponding Lagrangian interpolation function. */
    for(unsigned short i=0; i<nVarNDOFs; ++i)
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
        /*--- integration point. It is likely that this integration point  ---*/
        /*--- coincides with the location of one of the time DOFs. When    ---*/
        /*--- this is the case the interpolation boils down to a copy.     ---*/
        /*--------------------------------------------------------------------*/

        /* Store the interpolation data for this time integration point. */
        const su2double *DOFToThisTimeInt = timeInterpolDOFToIntegrationADER_DG
                                          + intPoint*nTimeDOFs;

        /* Initialize the interpolated solution to zero. */
        for(unsigned short i=0; i<nVarNDOFs; ++i) solInt[i] = 0.0;

        /* Carry out the actual interpolation. */
        for(unsigned short j=0; j<nTimeDOFs; ++j) {

          const su2double *solPredTimeInd = solPred + j*nVarNDOFs;
          for(unsigned short i=0; i<nVarNDOFs; ++i)
            solInt[i] += DOFToThisTimeInt[j]*solPredTimeInd[i];
        }

        /*--------------------------------------------------------------------*/
        /*--- Compute the spatial residual of the predictor step for the   ---*/
        /*--- current time integration point. A distinction is             ---*/
        /*--- made between an aliased and a non-aliased evaluation of the  ---*/
        /*--- predictor residual and between 2D and 3D. The latter is for  ---*/
        /*--- performance reasons.                                         ---*/
        /*--------------------------------------------------------------------*/

        switch( nDim ) {
          case 2: {
            if( useAliasedPredictor )
              ADER_DG_AliasedPredictorResidual_2D(config, &volElem[l], solInt, resInt, work);
            else
              ADER_DG_NonAliasedPredictorResidual_2D(config, &volElem[l], solInt, resInt, work);

            break;
          }

          case 3: {
            if( useAliasedPredictor )
              ADER_DG_AliasedPredictorResidual_3D(config, &volElem[l], solInt, resInt, work);
            else
              ADER_DG_NonAliasedPredictorResidual_3D(config, &volElem[l], solInt, resInt, work);

            break;
          }
        }

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
          su2double *res = resTot + j*nVarNDOFs;
          for(unsigned short i=0; i<nVarNDOFs; ++i)
            res[i] -= w*resInt[i];
        }
      }

      /* Solve for the new values of solPred, which are obtained by
         carrying out the matrix product iterMat X resTot. */
      DenseMatrixProduct(nDOFsNTimeDOFs, nVar, nDOFsNTimeDOFs,
                         volElem[l].ADERIterationMatrix.data(),
                         resTot, solPred);

      /*--- Determine whether or not the iteration process is considered
            converged. Make a distinction between 2D and 3D for performance
            reasons. ---*/
      bool converged = true;

      switch( nDim ) {
        case 2: {
          /* 2D simulation. Compute the L2 norm of the updates. */
          su2double L2[] = {0.0, 0.0, 0.0, 0.0};

          for(unsigned short j=0; j<nDOFsNTimeDOFs; ++j) {
            const su2double *solP = solPred + j*nVar;
            const su2double *solO = solOld  + j*nVar;

            su2double diff;
            diff = solP[0] - solO[0]; L2[0] += diff*diff;
            diff = solP[1] - solO[1]; L2[1] += diff*diff;
            diff = solP[2] - solO[2]; L2[2] += diff*diff;
            diff = solP[3] - solO[3]; L2[3] += diff*diff;
          }

          L2[0] = sqrt(L2[0]/nDOFsNTimeDOFs);
          L2[1] = sqrt(L2[1]/nDOFsNTimeDOFs);
          L2[2] = sqrt(L2[2]/nDOFsNTimeDOFs);
          L2[3] = sqrt(L2[3]/nDOFsNTimeDOFs);

          if(L2[0] > TolSolADER[0]) converged = false;
          if(L2[1] > TolSolADER[1]) converged = false;
          if(L2[2] > TolSolADER[2]) converged = false;
          if(L2[3] > TolSolADER[3]) converged = false;

          break;
        }

        case 3: {
          /* 3D simulation. Compute the L2 norm of the updates. */
          su2double L2[] = {0.0, 0.0, 0.0, 0.0, 0.0};

          for(unsigned short j=0; j<nDOFsNTimeDOFs; ++j) {
            const su2double *solP = solPred + j*nVar;
            const su2double *solO = solOld  + j*nVar;

            su2double diff;
            diff = solP[0] - solO[0]; L2[0] += diff*diff;
            diff = solP[1] - solO[1]; L2[1] += diff*diff;
            diff = solP[2] - solO[2]; L2[2] += diff*diff;
            diff = solP[3] - solO[3]; L2[3] += diff*diff;
            diff = solP[4] - solO[4]; L2[4] += diff*diff;
          }

          L2[0] = sqrt(L2[0]/nDOFsNTimeDOFs);
          L2[1] = sqrt(L2[1]/nDOFsNTimeDOFs);
          L2[2] = sqrt(L2[2]/nDOFsNTimeDOFs);
          L2[3] = sqrt(L2[3]/nDOFsNTimeDOFs);
          L2[4] = sqrt(L2[4]/nDOFsNTimeDOFs);

          if(L2[0] > TolSolADER[0]) converged = false;
          if(L2[1] > TolSolADER[1]) converged = false;
          if(L2[2] > TolSolADER[2]) converged = false;
          if(L2[3] > TolSolADER[3]) converged = false;
          if(L2[4] > TolSolADER[4]) converged = false;

          break;
        }
      }

      /* Break the iteration loop if the solution is considered converged. */
      if( converged ) break;
    }

    /* Store the predictor solution in the correct location of
       VecSolDOFsPredictorADER. */
    for(unsigned short j=0; j<nTimeDOFs; ++j) {
      su2double *solCur      = VecSolDOFsPredictorADER.data()
                             + nVar*(j*nDOFsLocTot + volElem[l].offsetDOFsSolLocal);
      su2double *solPredTime = solPred + j*nVarNDOFs;

      memcpy(solCur, solPredTime, nBytes);
    }
  }
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_2D(CConfig           *config,
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
  su2double *gradFluxesInt = fluxesDOF + 8*nDOFs; /* nDOFs*nVar*nDim. */
  su2double *divFlux       = work;

  /* Determine the offset between the r-derivatives and s-derivatives of
     the fluxes. */
  const unsigned short offDeriv = 8*nInt;  /* nVar*nInt*nDim. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 5;  /* nDim*nDim + 1. */

  /*-- Compute the Cartesian fluxes in the DOFs. ---*/
  for(unsigned short i=0; i<nDOFs; ++i) {

    /* Set the pointers to the location where the solution of this DOF is
       stored and the location where the Cartesian fluxes are stored. */
    const su2double *solDOF = sol       + 4*i;    /* nVar*i */
    su2double       *fluxes = fluxesDOF + 8*i;    /* nVar*nDim*i */

    /* Set the pointer to the grid velocities at the location of the
       solution DOFS. */
    const su2double *gridVel = elem->gridVelocitiesSolDOFs.data() + 2*i; /* nDim*i. */

    /*--- Compute the velocities and pressure in this DOF. ---*/
    const su2double DensityInv   = 1.0/solDOF[0];
    const su2double u            = DensityInv*solDOF[1];
    const su2double v            = DensityInv*solDOF[2];
    const su2double StaticEnergy = DensityInv*solDOF[3] - 0.5*(u*u + v*v);

    FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
    const su2double Pressure = FluidModel->GetPressure();

    /*--- Compute the Cartesian fluxes in this DOF. First the fluxes in
          x-direction are stored, followed by the fluxes in y-direction. ---*/
    const su2double uRel = u - gridVel[0];
    const su2double vRel = v - gridVel[1];

    fluxes[0] = solDOF[0]*uRel;
    fluxes[1] = solDOF[1]*uRel + Pressure;
    fluxes[2] = solDOF[2]*uRel;
    fluxes[3] = solDOF[3]*uRel + Pressure*u;

    fluxes[4] = solDOF[0]*vRel;
    fluxes[5] = solDOF[1]*vRel;
    fluxes[6] = solDOF[2]*vRel + Pressure;
    fluxes[7] = solDOF[3]*vRel + Pressure*v;
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
    const su2double *gradFluxesDr = gradFluxesInt + 8*i;       /* nVar*nDim*i. */
    const su2double *gradFluxesDs = gradFluxesDr  + offDeriv;

    /* Easier storage of the metric terms in this integration point.
       Also set the point where the divergence of the flux terms is
       stored for the current integration point. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
    su2double       *divFluxInt  = divFlux + 4*i;  /* nVar*i. */

    /* Compute the metric terms multiplied by the integration weight. Note that the
       first term in the metric terms is the Jacobian. */
    const su2double wDrdx = weights[i]*metricTerms[1];
    const su2double wDrdy = weights[i]*metricTerms[2];

    const su2double wDsdx = weights[i]*metricTerms[3];
    const su2double wDsdy = weights[i]*metricTerms[4];

    /* Compute the divergence of the fluxes, multiplied by the integration weight. */
    divFluxInt[0] = gradFluxesDr[0]*wDrdx + gradFluxesDs[0]*wDsdx
                  + gradFluxesDr[4]*wDrdy + gradFluxesDs[4]*wDsdy;
    divFluxInt[1] = gradFluxesDr[1]*wDrdx + gradFluxesDs[1]*wDsdx
                  + gradFluxesDr[5]*wDrdy + gradFluxesDs[5]*wDsdy;
    divFluxInt[2] = gradFluxesDr[2]*wDrdx + gradFluxesDs[2]*wDsdx
                  + gradFluxesDr[6]*wDrdy + gradFluxesDs[6]*wDsdy;
    divFluxInt[3] = gradFluxesDr[3]*wDrdx + gradFluxesDs[3]*wDsdx
                  + gradFluxesDr[7]*wDrdy + gradFluxesDs[7]*wDsdy;
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_3D(CConfig           *config,
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
  su2double *gradFluxesInt = fluxesDOF + 15*nDOFs; /* nDOFs*nVar*nDim. */
  su2double *divFlux       = work;

  /* Determine the offset between the r-derivatives and s-derivatives, which is
     also the offset between s- and t-derivatives, of the fluxes. */
  const unsigned short offDeriv = 15*nInt; /* nVar*nInt*nDim. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 10;  /* nDim*nDim + 1. */

  /*-- Compute the Cartesian fluxes in the DOFs. ---*/
  for(unsigned short i=0; i<nDOFs; ++i) {

    /* Set the pointers to the location where the solution of this DOF is
       stored and the location where the Cartesian fluxes are stored. */
    const su2double *solDOF = sol       +  5*i;  /* i*nVar. */
    su2double       *fluxes = fluxesDOF + 15*i;  /* i*nVar*nDim. */

    /* Set the pointer to the grid velocities at the location of the
       solution DOFS. */
    const su2double *gridVel = elem->gridVelocitiesSolDOFs.data() + 3*i; /* nDim*i. */

    /*--- Compute the velocities and pressure in this DOF. ---*/
    const su2double DensityInv   = 1.0/solDOF[0];
    const su2double u            = DensityInv*solDOF[1];
    const su2double v            = DensityInv*solDOF[2];
    const su2double w            = DensityInv*solDOF[3];
    const su2double StaticEnergy = DensityInv*solDOF[4] - 0.5*(u*u + v*v + w*w);

    FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
    const su2double Pressure = FluidModel->GetPressure();

    /*--- Compute the Cartesian fluxes in this DOF. First the fluxes in
          x-direction are stored, followed by the fluxes in y-direction
          and the fluxes in z-direction. ---*/
    const su2double uRel = u - gridVel[0];
    const su2double vRel = v - gridVel[1];
    const su2double wRel = w - gridVel[2];

    fluxes[0] = solDOF[0]*uRel;
    fluxes[1] = solDOF[1]*uRel + Pressure;
    fluxes[2] = solDOF[2]*uRel;
    fluxes[3] = solDOF[3]*uRel;
    fluxes[4] = solDOF[4]*uRel + Pressure*u;

    fluxes[5] = solDOF[0]*vRel;
    fluxes[6] = solDOF[1]*vRel;
    fluxes[7] = solDOF[2]*vRel + Pressure;
    fluxes[8] = solDOF[3]*vRel;
    fluxes[9] = solDOF[4]*vRel + Pressure*v;

    fluxes[10] = solDOF[0]*wRel;
    fluxes[11] = solDOF[1]*wRel;
    fluxes[12] = solDOF[2]*wRel;
    fluxes[13] = solDOF[3]*wRel + Pressure;
    fluxes[14] = solDOF[4]*wRel + Pressure*w;
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
    const su2double *gradFluxesDr = gradFluxesInt + 15*i;       /* nVar*nDim*i. */
    const su2double *gradFluxesDs = gradFluxesDr  + offDeriv;
    const su2double *gradFluxesDt = gradFluxesDs  + offDeriv;

    /* Easier storage of the metric terms in this integration point.
       Also set the point where the divergence of the flux terms is
       stored for the current integration point. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
    su2double       *divFluxInt  = divFlux + 5*i;  /* nVar*i. */

    /* Compute the metric terms multiplied by the integration weight. Note that the
       first term in the metric terms is the Jacobian. */
    const su2double wDrdx = weights[i]*metricTerms[1];
    const su2double wDrdy = weights[i]*metricTerms[2];
    const su2double wDrdz = weights[i]*metricTerms[3];

    const su2double wDsdx = weights[i]*metricTerms[4];
    const su2double wDsdy = weights[i]*metricTerms[5];
    const su2double wDsdz = weights[i]*metricTerms[6];

    const su2double wDtdx = weights[i]*metricTerms[7];
    const su2double wDtdy = weights[i]*metricTerms[8];
    const su2double wDtdz = weights[i]*metricTerms[9];

    /* Compute the divergence of the fluxes, multiplied by the integration weight. */
    divFluxInt[0] = gradFluxesDr[0] *wDrdx + gradFluxesDs[0] *wDsdx + gradFluxesDt[0] *wDtdx
                  + gradFluxesDr[5] *wDrdy + gradFluxesDs[5] *wDsdy + gradFluxesDt[5] *wDtdy
                  + gradFluxesDr[10]*wDrdz + gradFluxesDs[10]*wDsdz + gradFluxesDt[10]*wDtdz;
    divFluxInt[1] = gradFluxesDr[1] *wDrdx + gradFluxesDs[1] *wDsdx + gradFluxesDt[1] *wDtdx
                  + gradFluxesDr[6] *wDrdy + gradFluxesDs[6] *wDsdy + gradFluxesDt[6] *wDtdy
                  + gradFluxesDr[11]*wDrdz + gradFluxesDs[11]*wDsdz + gradFluxesDt[11]*wDtdz;
    divFluxInt[2] = gradFluxesDr[2] *wDrdx + gradFluxesDs[2] *wDsdx + gradFluxesDt[2] *wDtdx
                  + gradFluxesDr[7] *wDrdy + gradFluxesDs[7] *wDsdy + gradFluxesDt[7] *wDtdy
                  + gradFluxesDr[12]*wDrdz + gradFluxesDs[12]*wDsdz + gradFluxesDt[12]*wDtdz;
    divFluxInt[3] = gradFluxesDr[3] *wDrdx + gradFluxesDs[3] *wDsdx + gradFluxesDt[3] *wDtdx
                  + gradFluxesDr[8] *wDrdy + gradFluxesDs[8] *wDsdy + gradFluxesDt[8] *wDtdy
                  + gradFluxesDr[13]*wDrdz + gradFluxesDs[13]*wDsdz + gradFluxesDt[13]*wDtdz;
    divFluxInt[4] = gradFluxesDr[4] *wDrdx + gradFluxesDs[4] *wDsdx + gradFluxesDt[4] *wDtdx
                  + gradFluxesDr[9] *wDrdy + gradFluxesDs[9] *wDsdy + gradFluxesDt[9] *wDtdy
                  + gradFluxesDr[14]*wDrdz + gradFluxesDs[14]*wDsdz + gradFluxesDt[14]*wDtdz;
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_2D(CConfig           *config,
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
     which is also the offset between the r- and s-derivatives. */
  const unsigned short offDeriv = 4*nInt;  /* nVar*nInt. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 5;  /* nDim*nDim + 1. */

  /* Compute the solution and the derivatives w.r.t. the parametric coordinates
     in the integration points. The first argument is nInt*(nDim+1). */
  DenseMatrixProduct(nInt*3, nVar, nDOFs, matBasisInt, sol, solAndGradInt);

  /*--- Loop over the integration points to compute the divergence of the inviscid
        fluxes in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the solution and gradient data
       of this integration point starts. */
    const su2double *sol   = solAndGradInt + 4*i; /* nVar*i. */
    const su2double *solDr = sol   + offDeriv;
    const su2double *solDs = solDr + offDeriv;

    /*--- Compute the velocities, pressure and total enthalpy
          in this integration point. ---*/
    const su2double DensityInv   = 1.0/sol[0];
    const su2double u            = DensityInv*sol[1];
    const su2double v            = DensityInv*sol[2];
    const su2double kinEnergy    = 0.5*(u*u + v*v);
    const su2double StaticEnergy = DensityInv*sol[3] - kinEnergy;

    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
    const su2double Pressure = FluidModel->GetPressure();
    const su2double Htot     = DensityInv*(sol[3] + Pressure);

    /* Easier storage of the grid velocities in this integration point. */
    const su2double *gridVel = elem->gridVelocities.data() + 2*i; /* nDim*i. */

    /* Easier storage of the metric terms. Note that the first term is the Jacobian. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

    const su2double drdx = metricTerms[1];
    const su2double drdy = metricTerms[2];
    const su2double dsdx = metricTerms[3];
    const su2double dsdy = metricTerms[4];

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the
          metric terms in this integration point. Note that these gradients
          must be scaled with the Jacobian. This scaling is already present
          in the metric terms. ---*/
    su2double solGradCart[4][2];

    solGradCart[0][0] = solDr[0]*drdx + solDs[0]*dsdx;
    solGradCart[1][0] = solDr[1]*drdx + solDs[1]*dsdx;
    solGradCart[2][0] = solDr[2]*drdx + solDs[2]*dsdx;
    solGradCart[3][0] = solDr[3]*drdx + solDs[3]*dsdx;

    solGradCart[0][1] = solDr[0]*drdy + solDs[0]*dsdy;
    solGradCart[1][1] = solDr[1]*drdy + solDs[1]*dsdy;
    solGradCart[2][1] = solDr[2]*drdy + solDs[2]*dsdy;
    solGradCart[3][1] = solDr[3]*drdy + solDs[3]*dsdy;

    /*--- Compute the divergence of the grid velocity.
          SET TO ZERO FOR NOW. THIS IS NOT CORRECT!!!!. ---*/
    const su2double divGridVel = 0.0;  // Must be multiplied by the Jacobian.

    /*--- Compute the Cartesian gradients of the pressure, scaled with
          the Jacobian. ---*/
    const su2double dpdx = Gamma_Minus_One*(solGradCart[3][0] + kinEnergy*solGradCart[0][0]
                         -                  u*solGradCart[1][0] - v*solGradCart[2][0]);
    const su2double dpdy = Gamma_Minus_One*(solGradCart[3][1] + kinEnergy*solGradCart[0][1]
                         -                  u*solGradCart[1][1] - v*solGradCart[2][1]);

    /*--- Abbreviations, which make it easier to compute the divergence term. ---*/
    const su2double abv1 =    solGradCart[1][0] +   solGradCart[2][1];
    const su2double abv2 = u* solGradCart[0][0] + v*solGradCart[0][1];
    const su2double abv3 = u*(solGradCart[3][0] + dpdx)
                         + v*(solGradCart[3][1] + dpdy);

    /*--- Set the pointer to store the divergence terms for this integration point.
          and compute these terms, multiplied by the integration weight. ---*/
    su2double *divFluxInt = divFlux + 4*i;  /* nVar*i. */

    divFluxInt[0] = weights[i]*(abv1 - sol[0]*divGridVel
                  -             gridVel[0]*solGradCart[0][0] - gridVel[1]*solGradCart[0][1]);
    divFluxInt[1] = weights[i]*(dpdx + u*(abv1 - abv2 + solGradCart[1][0]) + v*solGradCart[1][1]
                  -             sol[1]*divGridVel
                  -             gridVel[0]*solGradCart[1][0] - gridVel[1]*solGradCart[1][1]);
    divFluxInt[2] = weights[i]*(dpdy + v*(abv1 - abv2 + solGradCart[2][1]) + u*solGradCart[2][0]
                  -             sol[2]*divGridVel
                  -             gridVel[0]*solGradCart[2][0] - gridVel[1]*solGradCart[2][1]);
    divFluxInt[3] = weights[i]*(abv3 + Htot*(abv1 - abv2) - sol[3]*divGridVel
                  -             gridVel[0]*solGradCart[3][0] - gridVel[1]*solGradCart[3][1]);
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_3D(CConfig           *config,
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
  const unsigned short offDeriv = 5*nInt; /* nVar*nInt. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 10;  /* nDim*nDim + 1. */

  /* Compute the solution and the derivatives w.r.t. the parametric coordinates
     in the integration points. The first argument is nInt*(nDim+1). */
  DenseMatrixProduct(nInt*4, nVar, nDOFs, matBasisInt, sol, solAndGradInt);

  /*--- Loop over the integration points to compute the divergence of the inviscid
        fluxes in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the solution and gradient data
       of this integration point starts. */
    const su2double *sol   = solAndGradInt + 5*i; /* nVar*i. */
    const su2double *solDr = sol   + offDeriv;
    const su2double *solDs = solDr + offDeriv;
    const su2double *solDt = solDs + offDeriv;

    /*--- Compute the velocities, pressure and total enthalpy
          in this integration point. ---*/
    const su2double DensityInv   = 1.0/sol[0];
    const su2double u            = DensityInv*sol[1];
    const su2double v            = DensityInv*sol[2];
    const su2double w            = DensityInv*sol[3];
    const su2double kinEnergy    = 0.5*(u*u + v*v + w*w);
    const su2double StaticEnergy = DensityInv*sol[4] - kinEnergy;

    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
    const su2double Pressure = FluidModel->GetPressure();
    const su2double Htot     = DensityInv*(sol[4] + Pressure);

    /* Easier storage of the grid velocities in this integration point. */
    const su2double *gridVel = elem->gridVelocities.data() + 3*i; /* nDim*i. */

    /* Easier storage of the metric terms. Note that the first term is the Jacobian. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

    const su2double drdx = metricTerms[1];
    const su2double drdy = metricTerms[2];
    const su2double drdz = metricTerms[3];

    const su2double dsdx = metricTerms[4];
    const su2double dsdy = metricTerms[5];
    const su2double dsdz = metricTerms[6];

    const su2double dtdx = metricTerms[7];
    const su2double dtdy = metricTerms[8];
    const su2double dtdz = metricTerms[9];

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the
          metric terms in this integration point. Note that these gradients
          must be scaled with the Jacobian. This scaling is already present
          in the metric terms. ---*/
    su2double solGradCart[5][3];

    solGradCart[0][0] = solDr[0]*drdx + solDs[0]*dsdx + solDt[0]*dtdx;
    solGradCart[1][0] = solDr[1]*drdx + solDs[1]*dsdx + solDt[1]*dtdx;
    solGradCart[2][0] = solDr[2]*drdx + solDs[2]*dsdx + solDt[2]*dtdx;
    solGradCart[3][0] = solDr[3]*drdx + solDs[3]*dsdx + solDt[3]*dtdx;
    solGradCart[4][0] = solDr[4]*drdx + solDs[4]*dsdx + solDt[4]*dtdx;

    solGradCart[0][1] = solDr[0]*drdy + solDs[0]*dsdy + solDs[0]*dtdy;
    solGradCart[1][1] = solDr[1]*drdy + solDs[1]*dsdy + solDs[1]*dtdy;
    solGradCart[2][1] = solDr[2]*drdy + solDs[2]*dsdy + solDs[2]*dtdy;
    solGradCart[3][1] = solDr[3]*drdy + solDs[3]*dsdy + solDs[3]*dtdy;
    solGradCart[4][1] = solDr[4]*drdy + solDs[4]*dsdy + solDs[4]*dtdy;

    solGradCart[0][2] = solDr[0]*drdz + solDs[0]*dsdz + solDs[0]*dtdz;
    solGradCart[1][2] = solDr[1]*drdz + solDs[1]*dsdz + solDs[1]*dtdz;
    solGradCart[2][2] = solDr[2]*drdz + solDs[2]*dsdz + solDs[2]*dtdz;
    solGradCart[3][2] = solDr[3]*drdz + solDs[3]*dsdz + solDs[3]*dtdz;
    solGradCart[4][2] = solDr[4]*drdz + solDs[4]*dsdz + solDs[4]*dtdz;

    /*--- Compute the divergence of the grid velocity.
          SET TO ZERO FOR NOW. THIS IS NOT CORRECT!!!!. ---*/
    const su2double divGridVel = 0.0;   // Must be multiplied by the Jacobian.

     /*--- Compute the Cartesian gradients of the pressure, scaled with
          the Jacobian. ---*/
    su2double pGradCart[3];

    pGradCart[0] = Gamma_Minus_One*(solGradCart[4][0] + kinEnergy*solGradCart[0][0]
                 -                u*solGradCart[1][0] - v*solGradCart[2][0]
                 -                w*solGradCart[3][0]);
    pGradCart[1] = Gamma_Minus_One*(solGradCart[4][1] + kinEnergy*solGradCart[0][1]
                 -                u*solGradCart[1][1] - v*solGradCart[2][1]
                 -                w*solGradCart[3][1]);
    pGradCart[2] = Gamma_Minus_One*(solGradCart[4][2] + kinEnergy*solGradCart[0][2]
                 -                u*solGradCart[1][2] - v*solGradCart[2][2]
                 -                w*solGradCart[3][2]);

    /*--- Abbreviations, which make it easier to compute the divergence term. ---*/
    const su2double abv1 =    solGradCart[1][0] +   solGradCart[2][1] +   solGradCart[3][2];
    const su2double abv2 = u* solGradCart[0][0] + v*solGradCart[0][1] + w*solGradCart[0][2];
    const su2double abv3 = u*(solGradCart[4][0] +   pGradCart[0])
                         + v*(solGradCart[4][1] +   pGradCart[1])
                         + w*(solGradCart[4][2] +   pGradCart[2]);

    /*--- Set the pointer to store the divergence terms for this integration point.
          and compute these terms, multiplied by the integration weight. ---*/
    su2double *divFluxInt = divFlux + 5*i;  /* nVar*i. */

    divFluxInt[0] = weights[i]*(abv1
                  -             sol[0]*divGridVel            - gridVel[0]*solGradCart[0][0]
                  -             gridVel[1]*solGradCart[0][1] - gridVel[2]*solGradCart[0][2]);
    divFluxInt[1] = weights[i]*(pGradCart[0] + u*(abv1 - abv2 + solGradCart[1][0])
                  +             v*solGradCart[1][1] + w*solGradCart[1][2]
                  -             sol[1]*divGridVel            - gridVel[0]*solGradCart[1][0]
                  -             gridVel[1]*solGradCart[1][1] - gridVel[2]*solGradCart[1][2]);
    divFluxInt[2] = weights[i]*(pGradCart[1] + v*(abv1 - abv2 + solGradCart[2][1])
                  +             u*solGradCart[2][0] + w*solGradCart[2][2]
                  -             sol[2]*divGridVel            - gridVel[0]*solGradCart[2][0]
                  -             gridVel[1]*solGradCart[2][1] - gridVel[2]*solGradCart[2][2]);
    divFluxInt[3] = weights[i]*(pGradCart[2] + w*(abv1 - abv2 + solGradCart[3][2])
                  +             u*solGradCart[3][0] + v*solGradCart[3][1]
                  -             sol[3]*divGridVel            - gridVel[0]*solGradCart[3][0]
                  -             gridVel[1]*solGradCart[3][1] - gridVel[2]*solGradCart[3][2]);
    divFluxInt[4] = weights[i]*(abv3 + Htot*(abv1 - abv2)
                  -             sol[4]*divGridVel            - gridVel[0]*solGradCart[4][0]
                  -             gridVel[1]*solGradCart[4][1] - gridVel[2]*solGradCart[4][2]);
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_EulerSolver::ADER_DG_TimeInterpolatePredictorSol(CConfig             *config,
                                                              const unsigned short iTime,
                                                              const unsigned long  elemBeg,
                                                              const unsigned long  elemEnd,
                                                              const unsigned long  nAdjElem,
                                                              const unsigned long  *adjElem,
                                                              const bool           secondPartTimeInt,
                                                              su2double            *solTimeLevel) {

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Interpolate the solution to the given integration point    ---*/
  /*---         for the element range elemBeg to elemEnd. These elements   ---*/
  /*---         belong to the time level considered (not needed to know    ---*/
  /*---         the actual value in this function) and are therefore       ---*/
  /*---         continuous in memory.                                      ---*/
  /*--------------------------------------------------------------------------*/

  /* Easier storage of the interpolation coefficients for the current time
     integration point (iTime) for the elements of this time interval. */
  const unsigned short nTimeDOFs    = config->GetnTimeDOFsADER_DG();
  const su2double *DOFToThisTimeInt = timeInterpolDOFToIntegrationADER_DG
                                    + iTime*nTimeDOFs;

  /* Loop over the element range of this time level. */
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* Determine the number of solution variables for this element and
       set the pointer where the solution variables for this element must be
       stored in solTimeLevel. */
    const unsigned short nSolVar = nVar*volElem[l].nDOFsSol;
    su2double           *solDOFs = solTimeLevel + nVar*volElem[l].offsetDOFsSolThisTimeLevel;

    /* Initialize the solution to zero. */
    for(unsigned short i=0; i<nSolVar; ++i) solDOFs[i] = 0.0;

    /* Loop over the time DOFs, for which the predictor solution is present. */
    for(unsigned short j=0; j<nTimeDOFs; ++j) {

      /* Add the contribution of this predictor solution to the interpolated solution. */
      const su2double *solPred = VecSolDOFsPredictorADER.data()
                               + nVar*(j*nDOFsLocTot + volElem[l].offsetDOFsSolLocal);
      for(unsigned short i=0; i<nSolVar; ++i)
        solDOFs[i] += DOFToThisTimeInt[j]*solPred[i];
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Interpolate the solution to the given integration point    ---*/
  /*---         for the elements of the next time level, which are         ---*/
  /*---         adjacent to elements of the current time level. Note that  ---*/
  /*---         these elements are not contiguous in memory.               ---*/
  /*--------------------------------------------------------------------------*/

  /* For the faces with adjacent elements of a higher time level, the time
     integration takes place with twice the number of time integration points
     from the perspective of the higher time level. Hence the integration points
     iTime must be corrected if the state to be interpolated corresponds to the
     second part of the time integration interval. Perform this correction and
     determine the corresponding interpolation coefficients. */
  unsigned short iiTime = iTime;
  if( secondPartTimeInt ) iiTime += config->GetnTimeIntegrationADER_DG();

  DOFToThisTimeInt = timeInterpolAdjDOFToIntegrationADER_DG + iiTime*nTimeDOFs;

  /* Loop over the adjacent elements. */
  for(unsigned long l=0; l<nAdjElem; ++l) {
    const unsigned long ll = adjElem[l];

    /* Determine the number of solution variables for this element and
       set the pointer where the solution variables for this element must be
       stored in solTimeLevel. */
    const unsigned short nSolVar = nVar*volElem[ll].nDOFsSol;
    su2double           *solDOFs = solTimeLevel + nVar*volElem[ll].offsetDOFsSolPrevTimeLevel;

    /* Initialize the solution to zero. */
    for(unsigned short i=0; i<nSolVar; ++i) solDOFs[i] = 0.0;

    /* Loop over the time DOFs, for which the predictor solution is present. */
    for(unsigned short j=0; j<nTimeDOFs; ++j) {

      /* Add the contribution of this predictor solution to the interpolated solution. */
      const su2double *solPred = VecSolDOFsPredictorADER.data()
                               + nVar*(j*nDOFsLocTot + volElem[ll].offsetDOFsSolLocal);
      for(unsigned short i=0; i<nSolVar; ++i)
        solDOFs[i] += DOFToThisTimeInt[j]*solPred[i];
    }
  }
}

void CFEM_DG_EulerSolver::Shock_Capturing_DG(CConfig             *config,
                                             const unsigned long elemBeg,
                                             const unsigned long elemEnd) {

  /*--- Run shock capturing algorithm ---*/
  switch( config->GetKind_FEM_DG_Shock() ) {
    case NONE:
      break;
  }
}

void CFEM_DG_EulerSolver::Volume_Residual(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd) {

  /*--- Determine whether body force term is present in configuration. ---*/
  bool body_force = config->GetBody_Force();
  const su2double *body_force_vector = body_force ? config->GetBody_Force_Vector() : NULL;

  /*--- Set the pointers for the local arrays. ---*/
  su2double *sources = VecTmpMemory.data();
  su2double *solInt  = sources + nIntegrationMax*nVar;
  su2double *fluxes  = solInt  + nIntegrationMax*nVar;
  double tick = 0.0;

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /*--- Loop over the given element range to compute the contribution of the
        volume integral in the DG FEM formulation to the residual.  ---*/
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* Get the data from the corresponding standard element. */
    const unsigned short ind             = volElem[l].indStandardElement;
    const unsigned short nInt            = standardElementsSol[ind].GetNIntegration();
    const unsigned short nDOFs           = volElem[l].nDOFsSol;
    const su2double *matBasisInt         = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
    const su2double *matBasisIntTrans    = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
    const su2double *matDerBasisIntTrans = standardElementsSol[ind].GetDerMatBasisFunctionsIntTrans();
    const su2double *weights             = standardElementsSol[ind].GetWeightsIntegration();

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Interpolate the solution to the integration points of    ---*/
    /*---         the element.                                             ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the solution variables for this element. */
    const unsigned short timeLevel = volElem[l].timeLevel;
    const su2double *solDOFs = VecWorkSolDOFs[timeLevel].data()
                             + nVar*volElem[l].offsetDOFsSolThisTimeLevel;

    /* Call the general function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nInt, nVar, nDOFs, matBasisInt, solDOFs, solInt);
    config->GEMM_Tock(tick, "Volume_Residual1", nInt, nVar, nDOFs);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the inviscid fluxes, multiplied by minus the     ---*/
    /*---         integration weight, in the integration points.           ---*/
    /*---         If needed, also the source terms are computed.           ---*/
    /*------------------------------------------------------------------------*/

    /* Make a distinction between two and three space dimensions
        in order to have the most efficient code. */
    switch( nDim ) {

      case 2: {

        /* 2D simulation. Loop over the integration points to compute
           the fluxes. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the metric terms and grid velocities
             in this integration point. */
          const su2double *metricTerms = volElem[l].metricTerms.data()
                                       + i*nMetricPerPoint;
          const su2double Jac          = metricTerms[0];
          const su2double *gridVel     = volElem[l].gridVelocities.data() + 2*i;  /* nDim*i. */

          /* Compute the metric terms multiplied by minus the integration weight.
             The minus sign comes from the integration by parts in the weak
             formulation. */
          const su2double wDrdx = -weights[i]*metricTerms[1];
          const su2double wDrdy = -weights[i]*metricTerms[2];

          const su2double wDsdx = -weights[i]*metricTerms[3];
          const su2double wDsdy = -weights[i]*metricTerms[4];

          /* Easier storage of the location where the solution data of this
             integration point starts. */
          const su2double *sol = solInt + nVar*i;

          /*--- Compute the velocities and static energy in this integration point. ---*/
          const su2double rhoInv       = 1.0/sol[0];
          const su2double u            = sol[1]*rhoInv;
          const su2double v            = sol[2]*rhoInv;
          const su2double TotalEnergy  = sol[3]*rhoInv;
          const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

          /*--- Compute the pressure. ---*/
          FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
          const su2double Pressure = FluidModel->GetPressure();

          /* Compute the relative velocities w.r.t. the grid. */
          const su2double uRel = u - gridVel[0];
          const su2double vRel = v - gridVel[1];

          /* Set the pointer for the fluxes in this integration point. */
          su2double *flux = fluxes + i*nDim*nVar;

          /*--- Fluxes in r-direction. */
          const su2double Ur = uRel*wDrdx + vRel*wDrdy;

          flux[0] = sol[0]*Ur;
          flux[1] = sol[1]*Ur + Pressure*wDrdx;
          flux[2] = sol[2]*Ur + Pressure*wDrdy;
          flux[3] = sol[3]*Ur + Pressure*(u*wDrdx + v*wDrdy);

          /*--- Fluxes in s-direction. */
          const su2double Us = uRel*wDsdx + vRel*wDsdy;

          flux[4] = sol[0]*Us;
          flux[5] = sol[1]*Us + Pressure*wDsdx;
          flux[6] = sol[2]*Us + Pressure*wDsdy;
          flux[7] = sol[3]*Us + Pressure*(u*wDsdx + v*wDsdy);

          /*--- If needed, compute the body forces in this integration point.
                Note that the source terms are multiplied with minus the
                integration weight in order to be consistent with the
                formulation of the residual. ---*/
          if( body_force ) {
            su2double *source         = sources + i*nVar;
            const su2double weightJac = weights[i]*Jac;

            source[0] =  0.0;
            source[1] = -weightJac*body_force_vector[0];
            source[2] = -weightJac*body_force_vector[1];
            source[3] = -weightJac*(u*body_force_vector[0] + v*body_force_vector[1]);
          }
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /* 3D simulation. Loop over the integration points to compute
           the fluxes. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the metric terms and grid velocities
             in this integration point. */
          const su2double *metricTerms = volElem[l].metricTerms.data()
                                       + i*nMetricPerPoint;
          const su2double Jac          = metricTerms[0];
          const su2double *gridVel     = volElem[l].gridVelocities.data() + 3*i;  /* nDim*i. */

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
          const su2double *sol    = solInt + nVar*i;

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

          /* Compute the relative velocities w.r.t. the grid. */
          const su2double uRel = u - gridVel[0];
          const su2double vRel = v - gridVel[1];
          const su2double wRel = w - gridVel[2];

          /* Set the pointer for the fluxes in this integration point. */
          su2double *flux = fluxes + i*nDim*nVar;

          /*--- Fluxes in r-direction. */
          const su2double Ur = uRel*wDrdx + vRel*wDrdy + wRel*wDrdz;

          flux[0] = sol[0]*Ur;
          flux[1] = sol[1]*Ur + Pressure*wDrdx;
          flux[2] = sol[2]*Ur + Pressure*wDrdy;
          flux[3] = sol[3]*Ur + Pressure*wDrdz;
          flux[4] = sol[4]*Ur + Pressure*(u*wDrdx + v*wDrdy + w*wDrdz);

          /*--- Fluxes in s-direction. */
          const su2double Us = uRel*wDsdx + vRel*wDsdy + wRel*wDsdz;

          flux[5] = sol[0]*Us;
          flux[6] = sol[1]*Us + Pressure*wDsdx;
          flux[7] = sol[2]*Us + Pressure*wDsdy;
          flux[8] = sol[3]*Us + Pressure*wDsdz;
          flux[9] = sol[4]*Us + Pressure*(u*wDsdx + v*wDsdy + w*wDsdz);

          /*--- Fluxes in t-direction. */
          const su2double Ut = uRel*wDtdx + vRel*wDtdy + wRel*wDtdz;

          flux[10] = sol[0]*Ut;
          flux[11] = sol[1]*Ut + Pressure*wDtdx;
          flux[12] = sol[2]*Ut + Pressure*wDtdy;
          flux[13] = sol[3]*Ut + Pressure*wDtdz;
          flux[14] = sol[4]*Ut + Pressure*(u*wDtdx + v*wDtdy + w*wDtdz);

          /*--- If needed, compute the body forces in this integration point.
                Note that the source terms are multiplied with minus the
                integration weight in order to be consistent with the
                formulation of the residual. ---*/
          if( body_force ) {
            su2double *source         = sources + i*nVar;
            const su2double weightJac = weights[i]*Jac;

            source[0] =  0.0;
            source[1] = -weightJac*body_force_vector[0];
            source[2] = -weightJac*body_force_vector[1];
            source[3] = -weightJac*body_force_vector[2];
            source[4] = -weightJac*(u*body_force_vector[0] + v*body_force_vector[1]
                      +             w*body_force_vector[2]);
          }
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

    /* Add the contribution from the source terms, if needed. Use solInt
       as temporary storage for the matrix product. */
    if( body_force ) {

      /* Call the general function to carry out the matrix product. */
      config->GEMM_Tick(&tick);
      DenseMatrixProduct(nDOFs, nVar, nInt, matBasisIntTrans, sources, solInt);
      config->GEMM_Tock(tick, "Volume_Residual3", nDOFs, nVar, nInt);

      /* Add the residuals due to source terms to the volume residuals */
      for(unsigned short i=0; i<(nDOFs*nVar); ++i)
        res[i] += solInt[i];
    }
  }
}

void CFEM_DG_EulerSolver::Boundary_Conditions(const unsigned short timeLevel,
                                              CConfig              *config,
                                              CNumerics            **numerics){

  /* Loop over all boundaries. */
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /* Determine the range of faces for this time level and test if any
       surface element for this marker must be treated at all. */
    const unsigned long surfElemBeg = boundaries[iMarker].nSurfElem[timeLevel];
    const unsigned long surfElemEnd = boundaries[iMarker].nSurfElem[timeLevel+1];

    if(surfElemEnd > surfElemBeg) {

      /* Set the starting position in the vector for the face residuals for
         this boundary marker and time level and set the pointer to the boundary
         faces for this boundary marker. */
      su2double *resFaces = VecResFaces.data()
                          + nVar*startLocResFacesMarkers[iMarker][timeLevel];

      const CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

      /* Apply the appropriate boundary condition. */
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case EULER_WALL:
          BC_Euler_Wall(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                        numerics[CONV_BOUND_TERM]);
          break;
        case FAR_FIELD:
          BC_Far_Field(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                       numerics[CONV_BOUND_TERM]);
          break;
        case SYMMETRY_PLANE:
          BC_Sym_Plane(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                       numerics[CONV_BOUND_TERM]);
          break;
        case INLET_FLOW:
          BC_Inlet(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                   numerics[CONV_BOUND_TERM], iMarker);
          break;
        case OUTLET_FLOW:
          BC_Outlet(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                    numerics[CONV_BOUND_TERM], iMarker);
          break;
        case ISOTHERMAL:
          BC_Isothermal_Wall(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                             numerics[CONV_BOUND_TERM], iMarker);
          break;
        case HEAT_FLUX:
          BC_HeatFlux_Wall(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                           numerics[CONV_BOUND_TERM], iMarker);
          break;
        case RIEMANN_BOUNDARY:
          BC_Riemann(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                     numerics[CONV_BOUND_TERM], iMarker);
          break;
        case CUSTOM_BOUNDARY:
          BC_Custom(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                    numerics[CONV_BOUND_TERM]);
          break;
        case PERIODIC_BOUNDARY:  // Nothing to be done for a periodic boundary.
          break;
        default:
          SU2_MPI::Error("BC not implemented.", CURRENT_FUNCTION);
      }
    }
  }
}

void CFEM_DG_EulerSolver::ResidualFaces(CConfig             *config,
                                        const unsigned long indFaceBeg,
                                        const unsigned long indFaceEnd,
                                        unsigned long       &indResFaces,
                                        CNumerics           *numerics) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;
  double tick = 0.0;

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
            different from side 0. Carry out the matrix multiplication to obtain
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
  double tick = 0.0;

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Interpolate the left state in the integration points of  ---*/
  /*---         the face.                                                ---*/
  /*------------------------------------------------------------------------*/

  /* Get the required information from the corresponding standard face. */
  const unsigned short ind        = internalFace->indStandardElement;
  const unsigned short nInt       = standardMatchingFacesSol[ind].GetNIntegration();
  const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
  const su2double     *basisFace0 = standardMatchingFacesSol[ind].GetBasisFaceIntegrationSide0();

  /* Determine the time level of the face. */
  const unsigned long  elem0         = internalFace->elemID0;
  const unsigned long  elem1         = internalFace->elemID1;
  const unsigned short timeLevelFace = min(volElem[elem0].timeLevel,
                                           volElem[elem1].timeLevel);

  /* The numbering of the DOFs is given for the entire grid. However, in the
     working vectors for the solution a numbering per time level is used.
     Therefore an offset must be applied, which can be obtained from the
     corresponding volume element. Note that this offset is non-negative
     and depends whether or not the element belongs to the same time
     level as the face. */
  unsigned long offset;
  if(timeLevelFace == volElem[elem0].timeLevel)
    offset = volElem[elem0].offsetDOFsSolLocal - volElem[elem0].offsetDOFsSolThisTimeLevel;
  else
    offset = volElem[elem0].offsetDOFsSolLocal - volElem[elem0].offsetDOFsSolPrevTimeLevel;

  /*--- Store the solution of the DOFs of side 0 of the face in contiguous memory
        such that the function DenseMatrixProduct can be used to compute the left
        states in the integration points of the face, i.e. side 0. ---*/
  const unsigned long nBytes = nVar*sizeof(su2double);
  const unsigned long *DOFs = internalFace->DOFsSolFaceSide0.data();
  for(unsigned short i=0; i<nDOFsFace0; ++i) {
    const su2double *solDOF = VecWorkSolDOFs[timeLevelFace].data()
                            + nVar*(DOFs[i] - offset);
    su2double       *sol    = solFace + nVar*i;
    memcpy(sol, solDOF, nBytes);
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

  /* Determine the offset between the numbering used in the DOFs of the face
     and the numbering used in the working solution vectors, see above for a
     more detailed explanation. */
  if(timeLevelFace == volElem[elem1].timeLevel)
    offset = volElem[elem1].offsetDOFsSolLocal - volElem[elem1].offsetDOFsSolThisTimeLevel;
  else
    offset = volElem[elem1].offsetDOFsSolLocal - volElem[elem1].offsetDOFsSolPrevTimeLevel;

  /*--- Store the solution of the DOFs of side 1 of the face in contiguous memory
        such that the function DenseMatrixProduct can be used to compute the right
        states in the integration points of the face, i.e. side 1. ---*/
  DOFs = internalFace->DOFsSolFaceSide1.data();
  for(unsigned short i=0; i<nDOFsFace1; ++i) {
    const su2double *solDOF = VecWorkSolDOFs[timeLevelFace].data()
                            + nVar*(DOFs[i] - offset);
    su2double       *sol    = solFace + nVar*i;
    memcpy(sol, solDOF, nBytes);
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
                            internalFace->gridVelocities.data(),
                            solIntL, solIntR, fluxes, numerics);
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADEROwnedElem(
                                                     CConfig             *config,
                                                     const unsigned short timeLevel,
                                                     const unsigned short intPoint) {

  /* Compute half the integration weight. The reason for doing this is that the
     given integration weight is based on the normalized interval [-1..1], i.e.
     a length of two. Also compute a quarter of the weight, which is necessary
     to accumulate the face residuals of the DOFs of elements adjacent to
     element of a lower time level. */
  const su2double *timeIntegrationWeights = config->GetWeightsIntegrationADER_DG();

  const su2double halfWeight  = 0.5 *timeIntegrationWeights[intPoint];
  const su2double quartWeight = 0.25*timeIntegrationWeights[intPoint];

  /* Determine the owned element range that must be considered. */
  const unsigned long elemBegOwned = nVolElemOwnedPerTimeLevel[timeLevel];
  const unsigned long elemEndOwned = nVolElemOwnedPerTimeLevel[timeLevel+1];

  /* Add the residuals coming from the volume integral to VecTotResDOFsADER. */
  for(unsigned long l=elemBegOwned; l<elemEndOwned; ++l) {
    const unsigned long offset  = nVar*volElem[l].offsetDOFsSolLocal;
    const su2double    *res     = VecResDOFs.data() + offset;
    su2double          *resADER = VecTotResDOFsADER.data() + offset;

    for(unsigned short i=0; i<(nVar*volElem[l].nDOFsSol); ++i)
      resADER[i] += halfWeight*res[i];
  }

  /* Add the residuals coming from the surface integral to VecTotResDOFsADER.
     This part is from faces with the same time level as the element. */
  for(unsigned long l=elemBegOwned; l<elemEndOwned; ++l) {
    for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {
      const unsigned long ii = volElem[l].offsetDOFsSolLocal + i;
      su2double *resADER = VecTotResDOFsADER.data() + nVar*ii;

      for(unsigned long j=nEntriesResFaces[ii]; j<nEntriesResFaces[ii+1]; ++j) {
        const su2double *resFace = VecResFaces.data() + nVar*entriesResFaces[j];
        for(unsigned short k=0; k<nVar; ++k)
          resADER[k] += halfWeight*resFace[k];
      }
    }
  }

  /* Check if this is not the last time level. */
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
  if(timeLevel < (nTimeLevels-1)) {

    /* There may exist faces of this time level, which have a neighboring
       element of the next time level. The residuals of the DOFs of such
       elements must be updated with the face residuals of the current
       time level. However, on these elements the time step is twice as
       large. This is taken into account by multiplying the face residual
       with quartWeight when accumulating. */
    const unsigned long nAdjElem = ownedElemAdjLowTimeLevel[timeLevel+1].size();
    const unsigned long *adjElem = ownedElemAdjLowTimeLevel[timeLevel+1].data();

    for(unsigned l=0; l<nAdjElem; ++l) {
      const unsigned long ll = adjElem[l];
      for(unsigned short i=0; i<volElem[ll].nDOFsSol; ++i) {
        const unsigned long ii = volElem[ll].offsetDOFsSolLocal + i;
        su2double *resADER = VecTotResDOFsADER.data() + nVar*ii;

        for(unsigned long j=nEntriesResAdjFaces[ii]; j<nEntriesResAdjFaces[ii+1]; ++j) {
          const su2double *resFace = VecResFaces.data() + nVar*entriesResAdjFaces[j];
          for(unsigned short k=0; k<nVar; ++k)
            resADER[k] += quartWeight*resFace[k];
        }
      }
    }
  }
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADERHaloElem(
                                                     CConfig             *config,
                                                     const unsigned short timeLevel,
                                                     const unsigned short intPoint) {

  /* Compute half the integration weight. The reason for doing this is that the
     given integration weight is based on the normalized interval [-1..1], i.e.
     a length of two. Also compute a quarter of the weight, which is necessary
     to accumulate the face residuals of the DOFs of elements adjacent to
     element of a lower time level. */
  const su2double *timeIntegrationWeights = config->GetWeightsIntegrationADER_DG();

  const su2double halfWeight  = 0.5 *timeIntegrationWeights[intPoint];
  const su2double quartWeight = 0.25*timeIntegrationWeights[intPoint];

  /* Determine the halo element range that must be considered. */
  const unsigned long elemBegHalo = nVolElemHaloPerTimeLevel[timeLevel];
  const unsigned long elemEndHalo = nVolElemHaloPerTimeLevel[timeLevel+1];

  /* Add the residuals coming from the surface integral to VecTotResDOFsADER.
     This part is from faces with the same time level as the element. */
  for(unsigned long l=elemBegHalo; l<elemEndHalo; ++l) {
    for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {
      const unsigned long ii = volElem[l].offsetDOFsSolLocal + i;
      su2double *resADER = VecTotResDOFsADER.data() + nVar*ii;

      for(unsigned long j=nEntriesResFaces[ii]; j<nEntriesResFaces[ii+1]; ++j) {
        const su2double *resFace = VecResFaces.data() + nVar*entriesResFaces[j];
        for(unsigned short k=0; k<nVar; ++k)
          resADER[k] += halfWeight*resFace[k];
      }
    }
  }

  /* Check if this is not the last time level. */
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
  if(timeLevel < (nTimeLevels-1)) {

    /* There may exist faces of this time level, which have a neighboring
       element of the next time level. The residuals of the DOFs of such
       elements must be updated with the face residuals of the current
       time level. However, on these elements the time step is twice as
       large. This is taken into account by multiplying the face residual
       with quartWeight when accumulating. */
    const unsigned long nAdjElem = haloElemAdjLowTimeLevel[timeLevel+1].size();
    const unsigned long *adjElem = haloElemAdjLowTimeLevel[timeLevel+1].data();

    for(unsigned l=0; l<nAdjElem; ++l) {
      const unsigned long ll = adjElem[l];
      for(unsigned short i=0; i<volElem[ll].nDOFsSol; ++i) {
        const unsigned long ii = volElem[ll].offsetDOFsSolLocal + i;
        su2double *resADER = VecTotResDOFsADER.data() + nVar*ii;

        for(unsigned long j=nEntriesResAdjFaces[ii]; j<nEntriesResAdjFaces[ii+1]; ++j) {
          const su2double *resFace = VecResFaces.data() + nVar*entriesResAdjFaces[j];
          for(unsigned short k=0; k<nVar; ++k)
            resADER[k] += quartWeight*resFace[k];
        }
      }
    }
  }
}

void CFEM_DG_EulerSolver::CreateFinalResidual(const unsigned short timeLevel,
                                              const bool ownedElements) {

  /* Determine the element range for which the final residual must
     be created. */
  unsigned long elemStart, elemEnd;
  if( ownedElements ) {
    elemStart = nVolElemOwnedPerTimeLevel[0];
    elemEnd   = nVolElemOwnedPerTimeLevel[timeLevel+1];
  }
  else {
    elemStart = nVolElemHaloPerTimeLevel[0];
    elemEnd   = nVolElemHaloPerTimeLevel[timeLevel+1];
  }

  /* For the halo elements the residual is initialized to zero. */
  if( !ownedElements ) {

    for(unsigned long l=elemStart; l<elemEnd; ++l) {
      su2double *resDOFsElem = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;
      for(unsigned short i=0; i<(nVar*volElem[l].nDOFsSol); ++i)
        resDOFsElem[i] = 0.0;
    }
  }

  /* Loop over the required element range. */
  for(unsigned long l=elemStart; l<elemEnd; ++l) {

    /* Loop over the DOFs of this element. */
    for(unsigned long i=volElem[l].offsetDOFsSolLocal;
                      i<(volElem[l].offsetDOFsSolLocal+volElem[l].nDOFsSol); ++i) {

      /* Create the final residual by summing up all contributions. */
      su2double *resDOF = VecResDOFs.data() + nVar*i;
      for(unsigned long j=nEntriesResFaces[i]; j<nEntriesResFaces[i+1]; ++j) {
        const su2double *resFace = VecResFaces.data() + nVar*entriesResFaces[j];
        for(unsigned short k=0; k<nVar; ++k)
          resDOF[k] += resFace[k];
      }
    }
  }
}

void CFEM_DG_EulerSolver::MultiplyResidualByInverseMassMatrix(
                                              CConfig            *config,
                                              const bool          useADER,
                                              const unsigned long elemBeg,
                                              const unsigned long elemEnd) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *tmpRes = VecTmpMemory.data();
  double tick = 0.0;

  /*--- Set the reference to the correct residual. This depends
        whether or not the ADER scheme is used. ---*/
  vector<su2double> &VecRes = useADER ? VecTotResDOFsADER : VecResDOFs;

  /* Loop over the owned volume elements. */
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* Easier storage of the residuals for this volume element. */
    su2double *res = VecRes.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Check whether a multiplication must be carried out with the inverse of
       the lumped mass matrix or the full mass matrix. Note that it is crucial
       that the test is performed with the lumpedMassMatrix and not with
       massMatrix. The reason is that for implicit time stepping schemes
       both arrays may be in use. */
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

  double tick = 0.0;

  /*--- Set the pointers for the local arrays and determine the number of bytes
        for the call to memcpy later in this function. ---*/
  su2double *solInt  = VecTmpMemory.data();
  su2double *solCopy = solInt + nIntegrationMax*nVar;
  const unsigned long nBytes = nVar*sizeof(su2double);

  /*--- Get the information of the angle of attack, reference area, etc. ---*/
  const su2double Alpha        = config->GetAoA()*PI_NUMBER/180.0;
  const su2double Beta         = config->GetAoS()*PI_NUMBER/180.0;
  const su2double RefArea      = config->GetRefArea();
  const su2double RefLength    = config->GetRefLength();
  const su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double *Origin      = config->GetRefOriginMoment(0);
  const bool grid_movement     = config->GetGrid_Movement();

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

  const su2double factor = 1.0/(0.5*RefDensity*RefArea*RefVel2);

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

          /* Get the required information from the corresponding standard face. */
          const unsigned short ind   = surfElem[l].indStandardElement;
          const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
          const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();
          const su2double *basisFace = standardBoundaryFacesSol[ind].GetBasisFaceIntegration();
          const su2double *weights   = standardBoundaryFacesSol[ind].GetWeightsIntegration();

          /* Easier storage of the DOFs of the face. */
          const unsigned long *DOFs = surfElem[l].DOFsSolFace.data();

          /* Copy the solution of the DOFs of the face such that it is
             contiguous in memory. */
          for(unsigned short i=0; i<nDOFs; ++i) {
            const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
            su2double       *sol    = solCopy + nVar*i;
            memcpy(sol, solDOF, nBytes);
          }

          /* Call the general function to carry out the matrix product to determine
             the solution in the integration points. */
          config->GEMM_Tick(&tick);
          DenseMatrixProduct(nInt, nVar, nDOFs, basisFace, solCopy, solInt);
          config->GEMM_Tock(tick, "Pressure_Forces", nInt, nVar, nDOFs);

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
              MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLength;
              MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLength;
            }
            MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLength;
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
          CD_Inv[iMarker]   =  ForceInviscid[0]*cos(Alpha)*cos(Beta)
                            +  ForceInviscid[1]*sin(Beta)
                            +  ForceInviscid[2]*sin(Alpha)*cos(Beta);
          CL_Inv[iMarker]   = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
          CSF_Inv[iMarker]  = -ForceInviscid[0]*sin(Beta)*cos(Alpha)
                            +  ForceInviscid[1]*cos(Beta)
                            -  ForceInviscid[2]*sin(Beta)*sin(Alpha);
          CEff_Inv[iMarker] =  CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
          CMx_Inv[iMarker]  =  MomentInviscid[0];
          CMy_Inv[iMarker]  =  MomentInviscid[1];
          CMz_Inv[iMarker]  =  MomentInviscid[2];
          CFx_Inv[iMarker]  =  ForceInviscid[0];
          CFy_Inv[iMarker]  =  ForceInviscid[1];
          CFz_Inv[iMarker]  =  ForceInviscid[2];
        }

        AllBound_CD_Inv  += CD_Inv[iMarker];
        AllBound_CL_Inv  += CL_Inv[iMarker];
        AllBound_CSF_Inv += CSF_Inv[iMarker];
        AllBound_CEff_Inv = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
        AllBound_CMx_Inv += CMx_Inv[iMarker];
        AllBound_CMy_Inv += CMy_Inv[iMarker];
        AllBound_CMz_Inv += CMz_Inv[iMarker];
        AllBound_CFx_Inv += CFx_Inv[iMarker];
        AllBound_CFy_Inv += CFy_Inv[iMarker];
        AllBound_CFz_Inv += CFz_Inv[iMarker];

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

  const su2double      RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  const unsigned short nRKStages     = config->GetnRKStep();

  for(unsigned short iVar=0; iVar<nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /* Set the pointer to the array where the new state vector must be stored.
     If this is the last RK step, the solution should be stored in
     VecSolDOFs, such that a new time step can be taken. Otherwise, the
     solution should be stored in the first time level of the working vectors
     for the solution. Note that only one working vector will be present for a
     RK scheme, because time accurate local time stepping is not possible yet
     with RK schemes. */
  su2double *solNew;
  if(iRKStep == (nRKStages-1)) solNew = VecSolDOFs.data();
  else                         solNew = VecWorkSolDOFs[0].data();

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Store the coordinate of the first vertex of this element to give an
       indication for the location of the maximum residual. */
    const unsigned long ind = volElem[l].nodeIDsGrid[0];
    const su2double *coor   = meshPoints[ind].coor;

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar*volElem[l].offsetDOFsSolLocal;
    const su2double *res        = VecResDOFs.data() + offset;
    const su2double *solDOFsOld = VecSolDOFs.data() + offset;
    su2double *solDOFs          = solNew            + offset;

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
  vector<su2double> rbuf(nVar);

  /*--- Disable the reduce for the residual to avoid overhead if requested. ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    SU2_MPI::Allreduce(Residual_RMS, rbuf.data(), nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(unsigned short iVar=0; iVar<nVar; ++iVar) {

      if (rbuf[iVar] != rbuf[iVar])
        SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

      SetRes_RMS(iVar, max(EPS*EPS, sqrt(rbuf[iVar]/nDOFsGlobal)));
    }
  }

#else
  /*--- Sequential mode. Check for a divergence of the solver and compute
        the L2-norm of the residuals. ---*/
  for(unsigned short iVar=0; iVar<nVar; ++iVar) {

    if(GetRes_RMS(iVar) != GetRes_RMS(iVar))
      SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

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
    const su2double *solDOFsOld = VecSolDOFs.data()    + offset;
    su2double *solDOFsNew       = VecSolDOFsNew.data() + offset;

    su2double *solDOFs;
    if(iRKStep < 3) solDOFs = VecWorkSolDOFs[0].data() + offset;
    else            solDOFs = VecSolDOFs.data()        + offset;

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
  vector<su2double> rbuf(nVar);

  /*--- Disable the reduce for the residual to avoid overhead if requested. ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    SU2_MPI::Allreduce(Residual_RMS, rbuf.data(), nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(unsigned short iVar=0; iVar<nVar; ++iVar) {

      if (rbuf[iVar] != rbuf[iVar])
        SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

      SetRes_RMS(iVar, max(EPS*EPS, sqrt(rbuf[iVar]/nDOFsGlobal)));
    }
  }

#else
  /*--- Sequential mode. Check for a divergence of the solver and compute
   the L2-norm of the residuals. ---*/
  for(unsigned short iVar=0; iVar<nVar; ++iVar) {

    if(GetRes_RMS(iVar) != GetRes_RMS(iVar))
      SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

    SetRes_RMS(iVar, max(EPS*EPS, sqrt(GetRes_RMS(iVar)/nDOFsGlobal)));
  }

#endif
}

void CFEM_DG_EulerSolver::ADER_DG_Iteration(const unsigned long elemBeg,
                                            const unsigned long elemEnd) {

  /*--- Update the solution by looping over the given range
        of volume elements. ---*/
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset = nVar*volElem[l].offsetDOFsSolLocal;
    su2double *res     = VecTotResDOFsADER.data() + offset;
    su2double *solDOFs = VecSolDOFs.data()        + offset;

    /* Loop over the DOFs for this element and update the solution.
       Initialize the residual to zero afterwards. */
    for(unsigned short i=0; i<(nVar*volElem[l].nDOFsSol); ++i) {
      solDOFs[i] -= VecDeltaTime[l]*res[i];
      res[i]      = 0.0;
    }
  }
}

void CFEM_DG_EulerSolver::BoundaryStates_Euler_Wall(CConfig                  *config,
                                                    const CSurfaceElementFEM *surfElem,
                                                    const su2double          *solIntL,
                                                    su2double                *solIntR) {

  /*--- Apply the inviscid wall boundary conditions to compute the right
        state in the integration points. There are two options. Either the
        normal velocity is negated or the normal velocity is set to zero,
        relative to the grid motion. Some experiments are needed to see which
        formulation gives better results. ---*/
  const unsigned short ind  = surfElem->indStandardElement;
  const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the left and right solution, the normals and the
       grid velocity for this integration point. */
    const su2double *UL      = solIntL + i*nVar;
          su2double *UR      = solIntR + i*nVar;
    const su2double *normals = surfElem->metricNormalsFace.data() + i*(nDim+1);
    const su2double *gridVel = surfElem->gridVelocities.data() + i*nDim;

    /* Compute the normal component of the momentum variables, relative
       to the grid motion. */
    su2double rVn = 0.0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      rVn += (UL[iDim+1]-UL[0]*gridVel[iDim])*normals[iDim];

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
}

void CFEM_DG_EulerSolver::BoundaryStates_Inlet(CConfig                  *config,
                                               const CSurfaceElementFEM *surfElem,
                                               unsigned short           val_marker,
                                               const su2double          *solIntL,
                                               su2double                *solIntR) {

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

  /*--- Apply the subsonic inlet boundary conditions to compute the right
        state in the integration points. ---*/
  const unsigned short ind  = surfElem->indStandardElement;
  const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the left and right solution and the normals
       for this integration point. */
    const su2double *UL      = solIntL + i*nVar;
          su2double *UR      = solIntR + i*nVar;
    const su2double *normals = surfElem->metricNormalsFace.data() + i*(nDim+1);

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
}

void CFEM_DG_EulerSolver::BoundaryStates_Outlet(CConfig                  *config,
                                                const CSurfaceElementFEM *surfElem,
                                                unsigned short           val_marker,
                                                const su2double          *solIntL,
                                                su2double                *solIntR) {

  /*--- Retrieve the specified back pressure for this outlet.
        Nondimensionalize, if necessary. ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double P_Exit = config->GetOutlet_Pressure(Marker_Tag);
  P_Exit = P_Exit/config->GetPressure_Ref();

  /*--- Apply the subsonic inlet boundary conditions to compute the right
        state in the integration points. ---*/
  const unsigned short ind  = surfElem->indStandardElement;
  const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the left and right solution and the normals
       for this integration point. */
    const su2double *UL      = solIntL + i*nVar;
          su2double *UR      = solIntR + i*nVar;
    const su2double *normals = surfElem->metricNormalsFace.data() + i*(nDim+1);

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
}

void CFEM_DG_EulerSolver::BoundaryStates_Riemann(CConfig                  *config,
                                                 const CSurfaceElementFEM *surfElem,
                                                 unsigned short           val_marker,
                                                 const su2double          *solIntL,
                                                 su2double                *solIntR) {

  /* Retrieve the corresponding string for this marker. */
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /* Determine the number of integration points for the face. */
  const unsigned short ind  = surfElem->indStandardElement;
  const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

  /*--- Make a distinction between the different type of boundary conditions
        and set the initial guess of the right states in the integration
        points of the face. ---*/
  switch( config->GetKind_Data_Riemann(Marker_Tag) ) {
    case TOTAL_CONDITIONS_PT: {

      /* Total conditions specified. Retrieve the non-dimensional total
         pressure and temperature as well as the flow direction. */
      su2double P_Total   = config->GetRiemann_Var1(Marker_Tag);
      su2double T_Total   = config->GetRiemann_Var2(Marker_Tag);
      su2double *Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);

      P_Total /= config->GetPressure_Ref();
      T_Total /= config->GetTemperature_Ref();

      /* Compute the total enthalpy and entropy from these values. */
      FluidModel->SetTDState_PT(P_Total, T_Total);
      const su2double Enthalpy_e = FluidModel->GetStaticEnergy()
                                 + FluidModel->GetPressure()/FluidModel->GetDensity();
      const su2double Entropy_e  = FluidModel->GetEntropy();

      /* Loop over the integration points of the face. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the left and right solution
           for this integration point. */
        const su2double *UL = solIntL + i*nVar;
        su2double       *UR = solIntR + i*nVar;

        /* Compute the velocity squared of the left point. This value is
           extrapolated and is therefore also the velocity of the right point. */
        const su2double DensityInv = 1.0/UL[0];
        su2double Velocity2_e = 0.0;
        for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
          const su2double vel = UL[iDim]*DensityInv;
          Velocity2_e        += vel*vel;
        }

        /* Determine the static enthalpy, the density and the internal
           and total energy per unit mass for the right state. */
        const su2double StaticEnthalpy_e = Enthalpy_e - 0.5*Velocity2_e;

        FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
        const su2double Density_e = FluidModel->GetDensity();
        const su2double StaticEnergy_e = FluidModel->GetStaticEnergy();
        const su2double Energy_e       = StaticEnergy_e + 0.5*Velocity2_e;

        /* Set the conservative variables of the right state. */
        UR[0]      = Density_e;
        UR[nDim+1] = Density_e*Energy_e;

        const su2double momTot = Density_e*sqrt(Velocity2_e);
        for(unsigned short iDim=0; iDim<nDim; ++iDim)
          UR[iDim+1] = momTot*Flow_Dir[iDim];
      }

      break;
    }

    case STATIC_SUPERSONIC_INFLOW_PT: {

      /* Supersonic inlet with specified pressure, temperature and Mach
         number. Retrieve the non-dimensional static pressure and
         temperature as well as the three components of the Mach number. */
      su2double P_static = config->GetRiemann_Var1(Marker_Tag);
      su2double T_static = config->GetRiemann_Var2(Marker_Tag);
      su2double *Mach    = config->GetRiemann_FlowDir(Marker_Tag);

      P_static /= config->GetPressure_Ref();
      T_static /= config->GetTemperature_Ref();

      /* Compute the prescribed density, static energy per unit mass
         and speed of sound. */
      FluidModel->SetTDState_PT(P_static, T_static);
      const su2double Density_e      = FluidModel->GetDensity();
      const su2double StaticEnergy_e = FluidModel->GetStaticEnergy();
      const su2double SoundSpeed     = FluidModel->GetSoundSpeed();

      /* Determine the magnitude of the Mach number. */
      su2double MachMag = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        MachMag += Mach[iDim]*Mach[iDim];
      MachMag = sqrt(MachMag);

      /* Determine the magnitude of the velocity and the prescribed
         conservative variables. */
      const su2double VelMag = MachMag*SoundSpeed;

      su2double UCons[5];
      UCons[0]      = Density_e;
      UCons[nDim+1] = Density_e*(StaticEnergy_e + 0.5*VelMag*VelMag);

      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        UCons[iDim+1] = Density_e*Mach[iDim]*SoundSpeed;

      /* Loop over the integration points of the face and set the right
         state in the integration points. */
      const unsigned long nBytes = nVar*sizeof(su2double);
      for(unsigned short i=0; i<nInt; ++i) {
        su2double *UR = solIntR + i*nVar;
        memcpy(UR, UCons, nBytes);
      }

      break;
    }

    case STATIC_SUPERSONIC_INFLOW_PD: {

      /* Supersonic inlet with specified pressure, temperature and Mach
         number. Retrieve the non-dimensional static pressure and
         temperature as well as the three components of the Mach number. */
      su2double P_static   = config->GetRiemann_Var1(Marker_Tag);
      su2double Rho_static = config->GetRiemann_Var2(Marker_Tag);
      su2double *Mach      = config->GetRiemann_FlowDir(Marker_Tag);

      P_static /= config->GetPressure_Ref();
      Rho_static /= config->GetDensity_Ref();

      /* Compute the prescribed pressure, static energy per unit mass
         and speed of sound. */
      FluidModel->SetTDState_Prho(P_static, Rho_static);
      const su2double Density_e      = FluidModel->GetDensity();
      const su2double StaticEnergy_e = FluidModel->GetStaticEnergy();
      const su2double SoundSpeed     = FluidModel->GetSoundSpeed();

      /* Determine the magnitude of the Mach number. */
      su2double MachMag = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        MachMag += Mach[iDim]*Mach[iDim];
      MachMag = sqrt(MachMag);

      /* Determine the magnitude of the velocity and the prescribed
         conservative variables. */
      const su2double VelMag = MachMag*SoundSpeed;

      su2double UCons[5];
      UCons[0]      = Density_e;
      UCons[nDim+1] = Density_e*(StaticEnergy_e + 0.5*VelMag*VelMag);

      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        UCons[iDim+1] = Density_e*Mach[iDim]*SoundSpeed;

      /* Loop over the integration points of the face and set the right
         state in the integration points. */
      const unsigned long nBytes = nVar*sizeof(su2double);
      for(unsigned short i=0; i<nInt; ++i) {
        su2double *UR = solIntR + i*nVar;
        memcpy(UR, UCons, nBytes);
      }

      break;
    }

    case DENSITY_VELOCITY: {

      /* Subsonic inlet with specified density, velocity magnitude and
         flow direction. Retrieve the non-dimensional data. */
      su2double Density_e = config->GetRiemann_Var1(Marker_Tag);
      su2double VelMag_e  = config->GetRiemann_Var2(Marker_Tag);
      su2double *Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);

      Density_e /= config->GetDensity_Ref();
      VelMag_e  /= config->GetVelocity_Ref();

      /* Loop over the integration points of the face. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the left and right solution
           for this integration point. */
        const su2double *UL = solIntL + i*nVar;
        su2double       *UR = solIntR + i*nVar;

        /* Compute the total energy per unit mass, which is extrapolated. */
        const su2double Energy_e = UL[nDim+1]/UL[0];

        /* Set the conservative variables of the right state. */
        UR[0]      = Density_e;
        UR[nDim+1] = Density_e*Energy_e;

        const su2double momTot = Density_e*VelMag_e;
        for(unsigned short iDim=0; iDim<nDim; ++iDim)
          UR[iDim+1] = momTot*Flow_Dir[iDim];
      }

      break;
    }

    case STATIC_PRESSURE: {

      /* Subsonic outlet with specified pressure.
         Retrieve the non-dimensional data. */
      const su2double Pressure_e = config->GetRiemann_Var1(Marker_Tag)
                                 / config->GetPressure_Ref();

      /* Loop over the integration points of the face. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the left and right solution
           for this integration point. */
        const su2double *UL = solIntL + i*nVar;
        su2double       *UR = solIntR + i*nVar;

        /* Extrapolate the density and set the thermodynamic state. */
        UR[0] = UL[0];
        FluidModel->SetTDState_Prho(Pressure_e, UR[0]);

        /* Extrapolate the velocity. As the density is also extrapolated,
           this means that the momentum variables are identical for UL and UR.
           Also determine the velocity squared. */
        const su2double DensityInv = 1.0/UL[0];
        su2double Velocity2_e = 0.0;
        for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
          UR[iDim] = UL[iDim];

          const su2double vel = UL[iDim]*DensityInv;
          Velocity2_e        += vel*vel;
        }

        /* Compute the total energy per unit volume. */
        UR[nDim+1] = UR[0]*(FluidModel->GetStaticEnergy() + 0.5*Velocity2_e);
      }

      break;
    }

    default: {
      SU2_MPI::Error("Invalid Riemann input!", CURRENT_FUNCTION);
    }
  }

  /***************************************************************************/
  /* Correction of the boundary state using a characteristic reconstruction. */
  /***************************************************************************/

  /* Make a distinction between 2D and 3D. */
  switch( nDim ) {
    case 2: {

      /* Two-dimensional simulation. Loop over the number of integration points. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the left and right solution and the normals
           for this integration point. */
        const su2double *UL      = solIntL + i*nVar;
              su2double *UR      = solIntR + i*nVar;
        const su2double *normals = surfElem->metricNormalsFace.data() + i*(nDim+1);

        /* Compute the difference in conservative variables. */
        const su2double dr  = UR[0] - UL[0];
        const su2double dru = UR[1] - UL[1];
        const su2double drv = UR[2] - UL[2];
        const su2double drE = UR[3] - UL[3];

        /* Compute the primitive variables of the left state. */
        const su2double tmp    = 1.0/UL[0];
        const su2double vxL    = tmp*UL[1];
        const su2double vyL    = tmp*UL[2];
        const su2double alphaL = 0.5*(vxL*vxL + vyL*vyL);
        const su2double eL     = tmp*UL[3] - alphaL;

        const su2double nx  = normals[0];
        const su2double ny  = normals[1];
        const su2double vnL = vxL*nx + vyL*ny;

        FluidModel->SetTDState_rhoe(UL[0], eL);

        const su2double aL  = FluidModel->GetSoundSpeed();
        const su2double a2L = aL*aL;
        const su2double pL  = FluidModel->GetPressure();
        const su2double HL  = (UL[3] + pL)*tmp;

        const su2double ovaL  = 1.0/aL;
        const su2double ova2L = 1.0/a2L;

        /* Compute the eigenvalues. */
        su2double lam1 = vnL + aL;
        su2double lam2 = vnL - aL;
        su2double lam3 = vnL;

        /* For the characteristic reconstruction, only the outgoing
           characteristics (negative eigenvalues) must be taken into account.
           Set the eigenvalues accordingly. */
        lam1 = (lam1 < 0.0) ? 1.0 : 0.0;
        lam2 = (lam2 < 0.0) ? 1.0 : 0.0;
        lam3 = (lam3 < 0.0) ? 1.0 : 0.0;

        /* Some abbreviations, which occur quite often in the reconstruction. */
        const su2double abv1 = 0.5*(lam1 + lam2);
        const su2double abv2 = 0.5*(lam1 - lam2);
        const su2double abv3 = abv1 - lam3;

        const su2double abv4 = Gamma_Minus_One*(alphaL*dr - vxL*dru - vyL*drv + drE);
        const su2double abv5 = nx*dru + ny*drv - vnL*dr;
        const su2double abv6 = abv3*abv4*ova2L + abv2*abv5*ovaL;
        const su2double abv7 = abv2*abv4*ovaL  + abv3*abv5;

        /* Compute the right state. */
        UR[0] = UL[0] + lam3*dr  + abv6;
        UR[1] = UL[1] + lam3*dru + vxL*abv6 + nx*abv7;
        UR[2] = UL[2] + lam3*drv + vyL*abv6 + ny*abv7;
        UR[3] = UL[3] + lam3*drE + HL*abv6  + vnL*abv7;
      }

      break;
    }

    case 3: {

      /* Three-dimensional simulation. Loop over the number of integration points. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the left and right solution and the normals
           for this integration point. */
        const su2double *UL      = solIntL + i*nVar;
              su2double *UR      = solIntR + i*nVar;
        const su2double *normals = surfElem->metricNormalsFace.data() + i*(nDim+1);

        /* Compute the difference in conservative variables. */
        const su2double dr  = UR[0] - UL[0];
        const su2double dru = UR[1] - UL[1];
        const su2double drv = UR[2] - UL[2];
        const su2double drw = UR[3] - UL[3];
        const su2double drE = UR[4] - UL[4];

        /* Compute the primitive variables of the left state. */
        const su2double tmp    = 1.0/UL[0];
        const su2double vxL    = tmp*UL[1];
        const su2double vyL    = tmp*UL[2];
        const su2double vzL    = tmp*UL[3];
        const su2double alphaL = 0.5*(vxL*vxL + vyL*vyL + vzL*vzL);
        const su2double eL     = tmp*UL[4] - alphaL;

        const su2double nx  = normals[0];
        const su2double ny  = normals[1];
        const su2double nz  = normals[2];
        const su2double vnL = vxL*nx + vyL*ny + vzL*nz;

        FluidModel->SetTDState_rhoe(UL[0], eL);

        const su2double aL  = FluidModel->GetSoundSpeed();
        const su2double a2L = aL*aL;
        const su2double pL  = FluidModel->GetPressure();
        const su2double HL  = (UL[4] + pL)*tmp;

        const su2double ovaL  = 1.0/aL;
        const su2double ova2L = 1.0/a2L;

        /* Compute the eigenvalues. */
        su2double lam1 = vnL + aL;
        su2double lam2 = vnL - aL;
        su2double lam3 = vnL;

        /* For the characteristic reconstruction, only the outgoing
           characteristics (negative eigenvalues) must be taken into account.
           Set the eigenvalues accordingly. */
        lam1 = (lam1 < 0.0) ? 1.0 : 0.0;
        lam2 = (lam2 < 0.0) ? 1.0 : 0.0;
        lam3 = (lam3 < 0.0) ? 1.0 : 0.0;

        /* Some abbreviations, which occur quite often in the reconstruction. */
        const su2double abv1 = 0.5*(lam1 + lam2);
        const su2double abv2 = 0.5*(lam1 - lam2);
        const su2double abv3 = abv1 - lam3;

        const su2double abv4 = Gamma_Minus_One*(alphaL*dr - vxL*dru - vyL*drv - vzL*drw + drE);
        const su2double abv5 = nx*dru + ny*drv + nz*drw - vnL*dr;
        const su2double abv6 = abv3*abv4*ova2L + abv2*abv5*ovaL;
        const su2double abv7 = abv2*abv4*ovaL  + abv3*abv5;

        /* Compute the right state. */
        UR[0] = UL[0] + lam3*dr  + abv6;
        UR[1] = UL[1] + lam3*dru + vxL*abv6 + nx*abv7;
        UR[2] = UL[2] + lam3*drv + vyL*abv6 + ny*abv7;
        UR[3] = UL[3] + lam3*drw + vzL*abv6 + nz*abv7;
        UR[4] = UL[4] + lam3*drE + HL*abv6  + vnL*abv7;
      }

      break;
    }
  }
}

void CFEM_DG_EulerSolver::BC_Euler_Wall(CConfig                  *config,
                                        const unsigned long      surfElemBeg,
                                        const unsigned long      surfElemEnd,
                                        const CSurfaceElementFEM *surfElem,
                                        su2double                *resFaces,
                                        CNumerics                *conv_numerics) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the right state by applying the inviscid wall BC's. */
    BoundaryStates_Euler_Wall(config, &surfElem[l], solIntL, solIntR);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Far_Field(CConfig                  *config,
                                       const unsigned long      surfElemBeg,
                                       const unsigned long      surfElemEnd,
                                       const CSurfaceElementFEM *surfElem,
                                       su2double                *resFaces,
                                       CNumerics                *conv_numerics) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

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

void CFEM_DG_EulerSolver::BC_Sym_Plane(CConfig                  *config,
                                       const unsigned long      surfElemBeg,
                                       const unsigned long      surfElemEnd,
                                       const CSurfaceElementFEM *surfElem,
                                       su2double                *resFaces,
                                       CNumerics                *conv_numerics) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Determine the number of integration points. */
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /* Make a distinction between two and three space dimensions
       in order to have the most efficient code. */
    switch( nDim ) {

      case 2: {

        /* 2D simulation. Loop over the integration points
           to apply the symmetry condition. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution and the normals
             for this integration point. */
          const su2double *UL      = solIntL + i*nVar;
                su2double *UR      = solIntR + i*nVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

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
          const su2double *UL      = solIntL + i*nVar;
                su2double *UR      = solIntR + i*nVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

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

void CFEM_DG_EulerSolver::BC_Inlet(CConfig                  *config,
                                   const unsigned long      surfElemBeg,
                                   const unsigned long      surfElemEnd,
                                   const CSurfaceElementFEM *surfElem,
                                   su2double                *resFaces,
                                   CNumerics                *conv_numerics,
                                   unsigned short           val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the right state by applying the subsonic inlet BC's. */
    BoundaryStates_Inlet(config, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Outlet(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    unsigned short           val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the right state by applying the subsonic outlet BC's. */
    BoundaryStates_Outlet(config, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Riemann(CConfig                  *config,
                                     const unsigned long      surfElemBeg,
                                     const unsigned long      surfElemEnd,
                                     const CSurfaceElementFEM *surfElem,
                                     su2double                *resFaces,
                                     CNumerics                *conv_numerics,
                                     unsigned short           val_marker) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the right state by applying the Riemann BC's. */
    BoundaryStates_Riemann(config, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                 solIntR, fluxes, resFaces, indResFaces);
  }
}

void CFEM_DG_EulerSolver::BC_Custom(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics) {

  /*--- Set the pointers for the local arrays. ---*/
  su2double *solIntL = VecTmpMemory.data();
  su2double *solIntR = solIntL + nIntegrationMax*nVar;
  su2double *fluxes  = solIntR + nIntegrationMax*nVar;

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /*--- Apply the boundary conditions to compute the right
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

      SU2_MPI::Error("No or wrong compiler directive specified. This is necessary for customized boundary conditions.",
                     CURRENT_FUNCTION);
#endif

    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
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
  double tick = 0.0;

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Compute the fluxes in the integration points using the   ---*/
  /*---         approximate Riemann solver.                              ---*/
  /*------------------------------------------------------------------------*/

  /* General function to compute the fluxes in the integration points. */
  ComputeInviscidFluxesFace(config, nInt, surfElem->metricNormalsFace.data(),
                            surfElem->gridVelocities.data(),
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
  double tick = 0.0;

  /* Easier storage of the DOFs of the face. */
  const unsigned long *DOFs = surfElem->DOFsSolFace.data();

  /* The numbering of the DOFs is given for the entire grid. However, in the
     working vectors for the solution a numbering per time level is used.
     Therefore an offset must be applied, which can be obtained from the
     corresponding volume element. Note that this offset is non-negative. */
  const unsigned long  volID     = surfElem->volElemID;
  const unsigned short timeLevel = volElem[volID].timeLevel;
  const unsigned long  offset    = volElem[volID].offsetDOFsSolLocal
                                 - volElem[volID].offsetDOFsSolThisTimeLevel;

  /* Copy the solution of the DOFs of the face such that it is contiguous
     in memory. */
  const unsigned long nBytes = nVar*sizeof(su2double);
  for(unsigned short i=0; i<nDOFs; ++i) {

    const su2double *solDOF = VecWorkSolDOFs[timeLevel].data() + nVar*(DOFs[i]-offset);
    su2double       *sol    = solFace + nVar*i;
    memcpy(sol, solDOF, nBytes);
  }

  /* Call the general function to carry out the matrix product. */
  config->GEMM_Tick(&tick);
  DenseMatrixProduct(nInt, nVar, nDOFs, basisFace, solFace, solIntL);
  config->GEMM_Tock(tick, "LeftStatesIntegrationPointsBoundaryFace", nInt, nVar, nDOFs);
}

void CFEM_DG_EulerSolver::ComputeInviscidFluxesFace(CConfig             *config,
                                                    const unsigned long nPoints,
                                                    const su2double     *normalsFace,
                                                    const su2double     *gridVelFace,
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
      switch( nDim ) {

        case 2: {

          /* Two dimensional simulation. Loop over the number of points. */
          for(unsigned long i=0; i<nPoints; ++i) {

            /* Easier storage of the left and right solution, the face normals,
               the grid velocities and the flux vector for this point. */
            const su2double *UL      = solL + i*nVar;
            const su2double *UR      = solR + i*nVar;
            const su2double *norm    = normalsFace + i*(nDim+1);
            const su2double *gridVel = gridVelFace + i*nDim;
                  su2double *flux    = fluxes + i*nVar;

            const su2double nx = norm[0], ny = norm[1], halfArea = 0.5*norm[2];

            /* Compute the normal grid velocity. */
            const su2double gridVelNorm = gridVel[0]*nx + gridVel[1]*ny;

            /*--- Compute the primitive variables of the left and right state. ---*/
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
            const su2double unAvg    = vnAvg - gridVelNorm;
            const su2double ovaAvg   = 1.0/aAvg;
            const su2double ova2Avg  = 1.0/a2Avg;

            /*--- Compute the absolute values of the three eigenvalues and
                  apply the entropy correction. ---*/
            su2double lam1 = fabs(unAvg + aAvg);
            su2double lam2 = fabs(unAvg - aAvg);
            su2double lam3 = fabs(unAvg);

            tmp  = Delta*max(lam1, lam2);
            lam1 = max(lam1, tmp);
            lam2 = max(lam2, tmp);
            lam3 = max(lam3, tmp);

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
            const su2double vnL = vxL*nx + vyL*ny;
            const su2double vnR = vxR*nx + vyR*ny;
            const su2double unL = vnL - gridVelNorm;
            const su2double unR = vnR - gridVelNorm;
            const su2double pa  = pL + pR;

            flux[0] = halfArea*(UL[0]*unL + UR[0]*unR - (lam3*dr + abv6));
            flux[1] = halfArea*(UL[1]*unL + UR[1]*unR + pa*nx
                    -           (lam3*dru + uAvg*abv6 + nx*abv7));
            flux[2] = halfArea*(UL[2]*unL + UR[2]*unR + pa*ny
                    -           (lam3*drv + vAvg*abv6 + ny*abv7));
            flux[3] = halfArea*(UL[3]*unL + UR[3]*unR + pL*vnL + pR*vnR
                    -           (lam3*drE + HAvg*abv6 + vnAvg*abv7));
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /* Three dimensional simulation. Loop over the number of points. */
          for(unsigned long i=0; i<nPoints; ++i) {

            /* Easier storage of the left and right solution, the face normals,
               the grid velocities and the flux vector for this point. */
            const su2double *UL   = solL + i*nVar;
            const su2double *UR   = solR + i*nVar;
            const su2double *norm = normalsFace + i*(nDim+1);
            const su2double *gridVel = gridVelFace + i*nDim;
                  su2double *flux = fluxes + i*nVar;

            const su2double nx = norm[0], ny = norm[1], nz = norm[2];
            const su2double halfArea = 0.5*norm[3];

            /* Compute the normal grid velocity. */
            const su2double gridVelNorm = gridVel[0]*nx + gridVel[1]*ny + gridVel[2]*nz;

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
            const su2double unAvg    = vnAvg - gridVelNorm;
            const su2double ovaAvg   = 1.0/aAvg;
            const su2double ova2Avg  = 1.0/a2Avg;

            /*--- Compute the absolute values of the three eigenvalues and
                  apply the entropy correction. ---*/
            su2double lam1 = fabs(unAvg + aAvg);
            su2double lam2 = fabs(unAvg - aAvg);
            su2double lam3 = fabs(unAvg);

            tmp  = Delta*max(lam1, lam2);
            lam1 = max(lam1, tmp);
            lam2 = max(lam2, tmp);
            lam3 = max(lam3, tmp);

            /* Hack to get a kind of Lax-Friedrichs flux. */
            //lam1 = std::max(lam1, lam2);
            //lam3 = lam2 = lam1;
            /* End hack. */

            /*--- Some abbreviations, which occur quite often in the dissipation terms. ---*/
            const su2double abv1 = 0.5*(lam1 + lam2);
            const su2double abv2 = 0.5*(lam1 - lam2);
            const su2double abv3 = abv1 - lam3;

            const su2double abv4 = gm1*(alphaAvg*dr - uAvg*dru - vAvg*drv -wAvg*drw + drE);
            const su2double abv5 = nx*dru + ny*drv + nz*drw - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            /*--- Compute the Roe flux vector, which is 0.5*(FL + FR - |A|(UR-UL)). ---*/
            const su2double vnL = vxL*nx + vyL*ny + vzL*nz;
            const su2double vnR = vxR*nx + vyR*ny + vzR*nz;
            const su2double unL = vnL - gridVelNorm;
            const su2double unR = vnR - gridVelNorm;
            const su2double pa  = pL + pR;

            flux[0] = halfArea*(UL[0]*unL + UR[0]*unR - (lam3*dr + abv6));
            flux[1] = halfArea*(UL[1]*unL + UR[1]*unR + pa*nx
                    -           (lam3*dru + uAvg*abv6 + nx*abv7));
            flux[2] = halfArea*(UL[2]*unL + UR[2]*unR + pa*ny
                    -           (lam3*drv + vAvg*abv6 + ny*abv7));
            flux[3] = halfArea*(UL[3]*unL + UR[3]*unR + pa*nz
                    -           (lam3*drw + wAvg*abv6 + nz*abv7));
            flux[4] = halfArea*(UL[4]*unL + UR[4]*unR + pL*vnL + pR*vnR
                    -           (lam3*drE + HAvg*abv6 + vnAvg*abv7));
          }

          break;
        }
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case LAX_FRIEDRICH: {

      /* Local Lax-Friedrich (Rusanov) flux Make a distinction between two and
         three space dimensions in order to have the most efficient code. */
      switch( nDim ) {

        case 2: {

          /* Two dimensional simulation. Loop over the number of points. */
          for(unsigned long i=0; i<nPoints; ++i) {

            /* Easier storage of the left and right solution, the face normals,
               the grid velocities and the flux vector for this point. */
            const su2double *UL      = solL + i*nVar;
            const su2double *UR      = solR + i*nVar;
            const su2double *norm    = normalsFace + i*(nDim+1);
            const su2double *gridVel = gridVelFace + i*nDim;
                  su2double *flux    = fluxes + i*nVar;

            const su2double nx = norm[0], ny = norm[1], halfArea = 0.5*norm[2];

            /* Compute the normal grid velocity. */
            const su2double gridVelNorm = gridVel[0]*nx + gridVel[1]*ny;

            /*--- Compute the primitive variables of the left and right state. ---*/
            su2double tmp       = 1.0/UL[0];
            const su2double vxL = tmp*UL[1];
            const su2double vyL = tmp*UL[2];
            const su2double pL  = gm1*(UL[3] - 0.5*(vxL*UL[1] + vyL*UL[2]));
            const su2double a2L = Gamma*pL*tmp;

            tmp                 = 1.0/UR[0];
            const su2double vxR = tmp*UR[1];
            const su2double vyR = tmp*UR[2];
            const su2double pR  = gm1*(UR[3] - 0.5*(vxR*UR[1] + vyR*UR[2]));
            const su2double a2R = Gamma*pR*tmp;

            /*--- Compute the difference of the conservative mean flow variables. ---*/
            const su2double dr  = UR[0] - UL[0];
            const su2double dru = UR[1] - UL[1];
            const su2double drv = UR[2] - UL[2];
            const su2double drE = UR[3] - UL[3];

            /*--- Compute the spectral radii of the left and right state
                  and take the maximum for the dissipation terms. ---*/
            const su2double vnL = vxL*nx + vyL*ny;
            const su2double vnR = vxR*nx + vyR*ny;
            const su2double unL = vnL - gridVelNorm;
            const su2double unR = vnR - gridVelNorm;

            const su2double radL = fabs(unL) + sqrt(fabs(a2L));
            const su2double radR = fabs(unR) + sqrt(fabs(a2R));
            const su2double rad  = max(radL, radR);

            /*--- Compute the flux vector, which is 0.5*(FL + FR - rad(UR-UL)). ---*/
            const su2double pa = pL + pR;

            flux[0] = halfArea*(UL[0]*unL + UR[0]*unR - rad*dr);
            flux[1] = halfArea*(UL[1]*unL + UR[1]*unR + pa*nx - rad*dru);
            flux[2] = halfArea*(UL[2]*unL + UR[2]*unR + pa*ny - rad*drv);
            flux[3] = halfArea*(UL[3]*unL + UR[3]*unR + pL*vnL + pR*vnR - rad*drE);
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /* Three dimensional simulation. Loop over the number of points. */
          for(unsigned long i=0; i<nPoints; ++i) {

            /* Easier storage of the left and right solution, the face normals,
               the grid velocities and the flux vector for this point. */
            const su2double *UL   = solL + i*nVar;
            const su2double *UR   = solR + i*nVar;
            const su2double *norm = normalsFace + i*(nDim+1);
            const su2double *gridVel = gridVelFace + i*nDim;
                  su2double *flux = fluxes + i*nVar;

            const su2double nx = norm[0], ny = norm[1], nz = norm[2];
            const su2double halfArea = 0.5*norm[3];

            /* Compute the normal grid velocity. */
            const su2double gridVelNorm = gridVel[0]*nx + gridVel[1]*ny + gridVel[2]*nz;

            /*--- Compute the primitive variables of the left and right state. ---*/
            su2double tmp       = 1.0/UL[0];
            const su2double vxL = tmp*UL[1];
            const su2double vyL = tmp*UL[2];
            const su2double vzL = tmp*UL[3];
            const su2double pL  = gm1*(UL[4] - 0.5*(vxL*UL[1] + vyL*UL[2] + vzL*UL[3]));
            const su2double a2L = Gamma*pL*tmp;

            tmp                 = 1.0/UR[0];
            const su2double vxR = tmp*UR[1];
            const su2double vyR = tmp*UR[2];
            const su2double vzR = tmp*UR[3];
            const su2double pR  = gm1*(UR[4] - 0.5*(vxR*UR[1] + vyR*UR[2] + vzR*UR[3]));
            const su2double a2R = Gamma*pR*tmp;

            /*--- Compute the difference of the conservative mean flow variables. ---*/
            const su2double dr  = UR[0] - UL[0];
            const su2double dru = UR[1] - UL[1];
            const su2double drv = UR[2] - UL[2];
            const su2double drw = UR[3] - UL[3];
            const su2double drE = UR[4] - UL[4];

            /*--- Compute the spectral radii of the left and right state
                  and take the maximum for the dissipation terms. ---*/
            const su2double vnL = vxL*nx + vyL*ny + vzL*nz;
            const su2double vnR = vxR*nx + vyR*ny + vzR*nz;
            const su2double unL = vnL - gridVelNorm;
            const su2double unR = vnR - gridVelNorm;

            const su2double radL = fabs(unL) + sqrt(fabs(a2L));
            const su2double radR = fabs(unR) + sqrt(fabs(a2R));
            const su2double rad  = max(radL, radR);

            /*--- Compute the flux vector, which is 0.5*(FL + FR - rad(UR-UL)). ---*/
            const su2double pa = pL + pR;

            flux[0] = halfArea*(UL[0]*unL + UR[0]*unR - rad*dr);
            flux[1] = halfArea*(UL[1]*unL + UR[1]*unR + pa*nx - rad*dru);
            flux[2] = halfArea*(UL[2]*unL + UR[2]*unR + pa*ny - rad*drv);
            flux[3] = halfArea*(UL[3]*unL + UR[3]*unR + pa*nz - rad*drw);
            flux[4] = halfArea*(UL[4]*unL + UR[4]*unR + pL*vnL + pR*vnR - rad*drE);
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

        /* Easier storage of the left and right solution, the face normals,
           the grid velocities and the flux vector for this point. */
        const su2double *UL      = solL + i*nVar;
        const su2double *UR      = solR + i*nVar;
        const su2double *norm    = normalsFace + i*(nDim+1);
        const su2double *gridVel = gridVelFace + i*nDim;
        su2double       *flux    = fluxes + i*nVar;

        /*--- Store and load the normal into numerics. ---*/
        for (unsigned short iDim = 0; iDim < nDim; ++iDim)
          Normal[iDim] = norm[iDim]*norm[nDim];
        numerics->SetNormal(Normal);

        /*--- Load the grid velocities into numerics. ---*/
        su2double vGrid[] = {0.0, 0.0, 0.0};
        for(unsigned short iDim=0; iDim<nDim; ++iDim)
          vGrid[iDim] = gridVel[iDim];
        numerics->SetGridVel(vGrid, vGrid);

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

void CFEM_DG_EulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iVar;
  unsigned long index;

  string UnstExt, text_line;
  ifstream restart_file;

  const bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);

  string restart_filename = config->GetSolution_FlowFileName();

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned short rbuf_NotMatching = 0;
  unsigned long nDOF_Read = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) {
        VecSolDOFs[nVar*iPoint_Local+iVar] = Restart_Data[index+iVar];
      }
      /*--- Update the local counter nDOF_Read. ---*/
      ++nDOF_Read;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/
  if(nDOF_Read < nDOFsLocOwned) rbuf_NotMatching = 1;

#ifdef HAVE_MPI
  unsigned short sbuf_NotMatching = rbuf_NotMatching;
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, MPI_COMM_WORLD);
#endif

  if (rbuf_NotMatching != 0)
    SU2_MPI::Error(string("The solution file ") + restart_filename.data() +
                   string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."),
                   CURRENT_FUNCTION);

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

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

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
        SU2_MPI::Error("Unknown SGS model encountered", CURRENT_FUNCTION);
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
  double tick = 0.0;

  /*--- Set the pointers for the local arrays and determine the number
        of bytes for the call to memcpy later in this function. ---*/
  su2double *solInt     = VecTmpMemory.data();
  su2double *gradSolInt = solInt     + nIntegrationMax*nVar;
  su2double *solCopy    = gradSolInt + nIntegrationMax*nVar*nDim;
  const unsigned long nBytes = nVar*sizeof(su2double);

  /*--- Get the information of the angle of attack, reference area, etc. ---*/
  const su2double Alpha        = config->GetAoA()*PI_NUMBER/180.0;
  const su2double Beta         = config->GetAoS()*PI_NUMBER/180.0;
  const su2double RefArea      = config->GetRefArea();
  const su2double RefLength    = config->GetRefLength();
  const su2double Gas_Constant = config->GetGas_ConstantND();
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

  const su2double factor = 1.0/(0.5*RefDensity*RefArea*RefVel2);

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

          /* Get the required information from the corresponding standard face. */
          const unsigned short ind       = surfElem[l].indStandardElement;
          const unsigned short nInt      = standardBoundaryFacesSol[ind].GetNIntegration();
          const unsigned short nDOFsFace = standardBoundaryFacesSol[ind].GetNDOFsFace();
          const unsigned short nDOFsElem = standardBoundaryFacesSol[ind].GetNDOFsElem();
          const su2double *basisFace     = standardBoundaryFacesSol[ind].GetBasisFaceIntegration();
          const su2double *derBasisElem  = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();
          const su2double *weights       = standardBoundaryFacesSol[ind].GetWeightsIntegration();

          /* Copy the solution of the DOFs of the face such that it is
             contiguous in memory. */
          for(unsigned short i=0; i<nDOFsFace; ++i) {
            const su2double *solDOF = VecSolDOFs.data() + nVar*surfElem[l].DOFsSolFace[i];
            su2double       *sol    = solCopy + nVar*i;
            memcpy(sol, solDOF, nBytes);
          }

          /* Call the general function to carry out the matrix product to determine
             the solution in the integration points. */
          config->GEMM_Tick(&tick);
          DenseMatrixProduct(nInt, nVar, nDOFsFace, basisFace, solCopy, solInt);
          config->GEMM_Tock(tick, "Friction_Forces", nInt, nVar, nDOFsFace);

          /*--- Store the solution of the DOFs of the adjacent element in contiguous
                memory such that the function DenseMatrixProduct can be used to compute
                the gradients solution variables in the integration points of the face. ---*/
          for(unsigned short i=0; i<nDOFsElem; ++i) {
            const su2double *solDOF = VecSolDOFs.data() + nVar*surfElem[l].DOFsSolElement[i];
            su2double       *sol    = solCopy + nVar*i;
            memcpy(sol, solDOF, nBytes);
          }

          /* Compute the gradients in the integration points. Call the general function to
             carry out the matrix product. */
          config->GEMM_Tick(&tick);
          DenseMatrixProduct(nInt*nDim, nVar, nDOFsElem, derBasisElem, solCopy, gradSolInt);
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
              MomentViscous[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLength;
              MomentViscous[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLength;
            }
            MomentViscous[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLength;

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

        su2double charVel2Max = 0.0, radViscMax = 0.0;

        /*--- Make a distinction between 2D and 3D for optimal performance. ---*/
        switch( nDim ) {

          case 2: {
            /*--- 2D simulation. Loop over the DOFs of this element
                  and determine the maximum wave speed and the maximum
                  value of the viscous spectral radius. ---*/
            for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
              const su2double *solDOF  = VecSolDOFs.data()
                                       + nVar*(volElem[i].offsetDOFsSolLocal + j);
              const su2double *gridVel = volElem[i].gridVelocitiesSolDOFs.data()
                                       + j*nDim;

              /* Compute the velocities and the internal energy per unit mass. */
              const su2double DensityInv   = 1.0/solDOF[0];
              const su2double u            = DensityInv*solDOF[1];
              const su2double v            = DensityInv*solDOF[2];
              const su2double StaticEnergy = DensityInv*solDOF[3] - 0.5*(u*u + v*v);

              /*--- Compute the maximum value of the wave speed. This is a rather
                    conservative estimate. ---*/
              FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
              const su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
              const su2double SoundSpeed  = sqrt(fabs(SoundSpeed2));

              const su2double radx     = fabs(u-gridVel[0]) + SoundSpeed;
              const su2double rady     = fabs(v-gridVel[1]) + SoundSpeed;
              const su2double charVel2 = radx*radx + rady*rady;

              charVel2Max = max(charVel2Max, charVel2);

              /* Compute the laminar kinematic viscosity and check if an eddy
                 viscosity must be determined. */
              const su2double muLam = FluidModel->GetLaminarViscosity();
              su2double muTurb      = 0.0;

              if( SGSModelUsed ) {

                /* Set the pointers to the locations where the gradients
                   of this DOF start. */
                const su2double *solDOFDr = gradSolDOFs + j*nVar;
                const su2double *solDOFDs = solDOFDr    + offDerivSol;

                /* Compute the true value of the metric terms in this DOF. Note that in
                   metricTerms the metric terms scaled by the Jacobian are stored. */
                const su2double *metricTerms = volElem[i].metricTermsSolDOFs.data()
                                             + j*nMetricPerPoint;
                const su2double JacInv       = 1.0/metricTerms[0];

                const su2double drdx = JacInv*metricTerms[1];
                const su2double drdy = JacInv*metricTerms[2];

                const su2double dsdx = JacInv*metricTerms[3];
                const su2double dsdy = JacInv*metricTerms[4];

                /*--- Compute the Cartesian gradients of the independent solution
                      variables from the gradients in parametric coordinates and the metric
                      terms in this DOF. ---*/
                su2double solGradCart[3][2];

                solGradCart[0][0] = solDOFDr[0]*drdx + solDOFDs[0]*dsdx;
                solGradCart[1][0] = solDOFDr[1]*drdx + solDOFDs[1]*dsdx;
                solGradCart[2][0] = solDOFDr[2]*drdx + solDOFDs[2]*dsdx;

                solGradCart[0][1] = solDOFDr[0]*drdy + solDOFDs[0]*dsdy;
                solGradCart[1][1] = solDOFDr[1]*drdy + solDOFDs[1]*dsdy;
                solGradCart[2][1] = solDOFDr[2]*drdy + solDOFDs[2]*dsdy;

                /*--- Compute the Cartesian gradients of the velocities. ---*/
                const su2double dudx = DensityInv*(solGradCart[1][0] - u*solGradCart[0][0]);
                const su2double dvdx = DensityInv*(solGradCart[2][0] - v*solGradCart[0][0]);
                const su2double dudy = DensityInv*(solGradCart[1][1] - u*solGradCart[0][1]);
                const su2double dvdy = DensityInv*(solGradCart[2][1] - v*solGradCart[0][1]);

                /* Compute the eddy viscosity. */
                const su2double dist = volElem[i].wallDistanceSolDOFs[j];
                muTurb = SGSModel->ComputeEddyViscosity_2D(solDOF[0], dudx, dudy,
                                                           dvdx, dvdy, lenScale,
                                                           dist);
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

            break;
          }

          case 3: {
            /*--- 3D simulation. Loop over the DOFs of this element
                  and determine the maximum wave speed and the maximum
                  value of the viscous spectral radius. ---*/
            for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
              const su2double *solDOF  = VecSolDOFs.data()
                                       + nVar*(volElem[i].offsetDOFsSolLocal + j);
              const su2double *gridVel = volElem[i].gridVelocitiesSolDOFs.data()
                                       + j*nDim;

              /* Compute the velocities and the internal energy per unit mass. */
              const su2double DensityInv   = 1.0/solDOF[0];
              const su2double u            = DensityInv*solDOF[1];
              const su2double v            = DensityInv*solDOF[2];
              const su2double w            = DensityInv*solDOF[3];
              const su2double StaticEnergy = DensityInv*solDOF[4] - 0.5*(u*u + v*v + w*w);

              /*--- Compute the maximum value of the wave speed. This is a rather
                    conservative estimate. ---*/
              FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
              const su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();
              const su2double SoundSpeed  = sqrt(fabs(SoundSpeed2));

              const su2double radx     = fabs(u-gridVel[0]) + SoundSpeed;
              const su2double rady     = fabs(v-gridVel[1]) + SoundSpeed;
              const su2double radz     = fabs(w-gridVel[2]) + SoundSpeed;
              const su2double charVel2 = radx*radx + rady*rady + radz*radz;

              charVel2Max = max(charVel2Max, charVel2);

              /* Compute the laminar kinematic viscosity and check if an eddy
                 viscosity must be determined. */
              const su2double muLam = FluidModel->GetLaminarViscosity();
              su2double muTurb      = 0.0;

              if( SGSModelUsed ) {

                /* Set the pointers to the locations where the gradients
                   of this DOF start. */
                const su2double *solDOFDr = gradSolDOFs + j*nVar;
                const su2double *solDOFDs = solDOFDr    + offDerivSol;
                const su2double *solDOFDt = solDOFDs    + offDerivSol;

                /* Compute the true value of the metric terms in this DOF. Note that in
                   metricTerms the metric terms scaled by the Jacobian are stored. */
                const su2double *metricTerms = volElem[i].metricTermsSolDOFs.data()
                                             + j*nMetricPerPoint;
                const su2double JacInv       = 1.0/metricTerms[0];

                const su2double drdx = JacInv*metricTerms[1];
                const su2double drdy = JacInv*metricTerms[2];
                const su2double drdz = JacInv*metricTerms[3];

                const su2double dsdx = JacInv*metricTerms[4];
                const su2double dsdy = JacInv*metricTerms[5];
                const su2double dsdz = JacInv*metricTerms[6];

                const su2double dtdx = JacInv*metricTerms[7];
                const su2double dtdy = JacInv*metricTerms[8];
                const su2double dtdz = JacInv*metricTerms[9];

                /*--- Compute the Cartesian gradients of the independent solution
                      variables from the gradients in parametric coordinates and the metric
                      terms in this DOF. ---*/
                su2double solGradCart[4][3];

                solGradCart[0][0] = solDOFDr[0]*drdx + solDOFDs[0]*dsdx + solDOFDt[0]*dtdx;
                solGradCart[1][0] = solDOFDr[1]*drdx + solDOFDs[1]*dsdx + solDOFDt[1]*dtdx;
                solGradCart[2][0] = solDOFDr[2]*drdx + solDOFDs[2]*dsdx + solDOFDt[2]*dtdx;
                solGradCart[3][0] = solDOFDr[3]*drdx + solDOFDs[3]*dsdx + solDOFDt[3]*dtdx;

                solGradCart[0][1] = solDOFDr[0]*drdy + solDOFDs[0]*dsdy + solDOFDt[0]*dtdy;
                solGradCart[1][1] = solDOFDr[1]*drdy + solDOFDs[1]*dsdy + solDOFDt[1]*dtdy;
                solGradCart[2][1] = solDOFDr[2]*drdy + solDOFDs[2]*dsdy + solDOFDt[2]*dtdy;
                solGradCart[3][1] = solDOFDr[3]*drdy + solDOFDs[3]*dsdy + solDOFDt[3]*dtdy;

                solGradCart[0][2] = solDOFDr[0]*drdz + solDOFDs[0]*dsdz + solDOFDt[0]*dtdz;
                solGradCart[1][2] = solDOFDr[1]*drdz + solDOFDs[1]*dsdz + solDOFDt[1]*dtdz;
                solGradCart[2][2] = solDOFDr[2]*drdz + solDOFDs[2]*dsdz + solDOFDt[2]*dtdz;
                solGradCart[3][2] = solDOFDr[3]*drdz + solDOFDs[3]*dsdz + solDOFDt[3]*dtdz;

                /*--- Compute the Cartesian gradients of the velocities. ---*/
                const su2double dudx = DensityInv*(solGradCart[1][0] - u*solGradCart[0][0]);
                const su2double dudy = DensityInv*(solGradCart[1][1] - u*solGradCart[0][1]);
                const su2double dudz = DensityInv*(solGradCart[1][2] - u*solGradCart[0][2]);

                const su2double dvdx = DensityInv*(solGradCart[2][0] - v*solGradCart[0][0]);
                const su2double dvdy = DensityInv*(solGradCart[2][1] - v*solGradCart[0][1]);
                const su2double dvdz = DensityInv*(solGradCart[2][2] - v*solGradCart[0][2]);

                const su2double dwdx = DensityInv*(solGradCart[3][0] - w*solGradCart[0][0]);
                const su2double dwdy = DensityInv*(solGradCart[3][1] - w*solGradCart[0][1]);
                const su2double dwdz = DensityInv*(solGradCart[3][2] - w*solGradCart[0][2]);

                /* Compute the eddy viscosity. */
                const su2double dist = volElem[i].wallDistanceSolDOFs[j];
                muTurb = SGSModel->ComputeEddyViscosity_3D(solDOF[0], dudx, dudy, dudz,
                                                           dvdx, dvdy, dvdz, dwdx, dwdy,
                                                           dwdz, lenScale, dist);
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

            break;
          }
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

      SU2_MPI::Error("Incompressible solver not implemented yet", CURRENT_FUNCTION);
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

void CFEM_DG_NSSolver::ADER_DG_AliasedPredictorResidual_2D(CConfig           *config,
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
  su2double *gradFluxesInt = fluxesDOF + 8*nDOFs; /* nDOFs*nVar*nDim. */
  su2double *gradSolDOFs   = gradFluxesInt;
  su2double *divFlux       = work;

  /* Determine the offset between the r-derivatives and s-derivatives of the
     fluxes in the integration points and the offset between the r-derivatives
     and s-derivatives of the solution in the DOFs. */
  const unsigned short offDerivSol    = 4*nDOFs;  /* nVar*nDOFs. */
  const unsigned short offDerivFluxes = 8*nInt;   /* nVar*nInt*nDim. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 5;  /* nDim*nDim + 1. */

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
    const su2double *solDOF   = sol         + 4*i;     /* i*nVar. */
    const su2double *solDOFDr = gradSolDOFs + 4*i;     /* i*nVar. */
    const su2double *solDOFDs = solDOFDr    + offDerivSol;
    su2double       *fluxes   = fluxesDOF   + 8*i;     /* i*nVar*nDim. */

    /* Set the pointer to the grid velocities at the location of the
       solution DOFS. */
    const su2double *gridVel = elem->gridVelocitiesSolDOFs.data() + 2*i; /* nDim*i. */

    /* Compute the true value of the metric terms in this DOF. Note that in
       metricTerms the metric terms scaled by the Jacobian are stored. */
    const su2double *metricTerms = elem->metricTermsSolDOFs.data()
                                 + i*nMetricPerPoint;
    const su2double JacInv       = 1.0/metricTerms[0];

    const su2double drdx = JacInv*metricTerms[1];
    const su2double drdy = JacInv*metricTerms[2];

    const su2double dsdx = JacInv*metricTerms[3];
    const su2double dsdy = JacInv*metricTerms[4];

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the metric
          terms in this DOF. ---*/
    su2double solGradCart[4][2];

    solGradCart[0][0] = solDOFDr[0]*drdx + solDOFDs[0]*dsdx;
    solGradCart[1][0] = solDOFDr[1]*drdx + solDOFDs[1]*dsdx;
    solGradCart[2][0] = solDOFDr[2]*drdx + solDOFDs[2]*dsdx;
    solGradCart[3][0] = solDOFDr[3]*drdx + solDOFDs[3]*dsdx;

    solGradCart[0][1] = solDOFDr[0]*drdy + solDOFDs[0]*dsdy;
    solGradCart[1][1] = solDOFDr[1]*drdy + solDOFDs[1]*dsdy;
    solGradCart[2][1] = solDOFDr[2]*drdy + solDOFDs[2]*dsdy;
    solGradCart[3][1] = solDOFDr[3]*drdy + solDOFDs[3]*dsdy;

    /*--- Compute the velocities, pressure and laminar viscosity in this DOF. ---*/
    const su2double DensityInv   = 1.0/solDOF[0];
    const su2double u            = DensityInv*solDOF[1];
    const su2double v            = DensityInv*solDOF[2];
    const su2double TotalEnergy  = DensityInv*solDOF[3];
    const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v);

    FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
    const su2double Pressure     = FluidModel->GetPressure();
    const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

    /*--- Compute the Cartesian gradients of the velocities and
          static energy in this DOF. ---*/
    const su2double dudx = DensityInv*(solGradCart[1][0] - u*solGradCart[0][0]);
    const su2double dvdx = DensityInv*(solGradCart[2][0] - v*solGradCart[0][0]);
    const su2double dudy = DensityInv*(solGradCart[1][1] - u*solGradCart[0][1]);
    const su2double dvdy = DensityInv*(solGradCart[2][1] - v*solGradCart[0][1]);

    const su2double dedx = DensityInv*(solGradCart[3][0] - TotalEnergy*solGradCart[0][0])
                         - u*dudx - v*dvdx;
    const su2double dedy = DensityInv*(solGradCart[3][1] - TotalEnergy*solGradCart[0][1])
                         - u*dudy - v*dvdy;

    /*--- Compute the eddy viscosity, if needed, and the total viscosity. ---*/
    su2double ViscosityTurb = 0.0;
    if( SGSModelUsed )
      ViscosityTurb = SGSModel->ComputeEddyViscosity_2D(solDOF[0], dudx, dudy,
                                                        dvdx, dvdy, lenScale,
                                                        elem->wallDistanceSolDOFs[i]);
    const su2double Viscosity = ViscosityLam + ViscosityTurb;

    /* Compute the total thermal conductivity divided by Cv. */
    const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                            + ViscosityTurb*factHeatFlux_Turb;

    /*--- Set the value of the second viscosity and compute the divergence
          term in the viscous normal stresses. ---*/
    const su2double lambda     = -TWO3*Viscosity;
    const su2double lamDivTerm =  lambda*(dudx + dvdy);

    /*--- Compute the viscous stress tensor. ---*/
    const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
    const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
    const su2double tauxy = Viscosity*(dudy + dvdx);

    /*--- Compute the Cartesian fluxes in this DOF. First the fluxes in
          x-direction are stored, followed by the fluxes in y-direction. ---*/
    const su2double uRel = u - gridVel[0];
    const su2double vRel = v - gridVel[1];

    fluxes[0] = solDOF[0]*uRel;
    fluxes[1] = solDOF[1]*uRel + Pressure - tauxx;
    fluxes[2] = solDOF[2]*uRel - tauxy;
    fluxes[3] = solDOF[3]*uRel + Pressure*u - kOverCv*dedx - u*tauxx - v*tauxy;;

    fluxes[4] = solDOF[0]*vRel;
    fluxes[5] = solDOF[1]*vRel - tauxy;
    fluxes[6] = solDOF[2]*vRel + Pressure - tauyy;
    fluxes[7] = solDOF[3]*vRel + Pressure*v - kOverCv*dedy - u*tauxy - v*tauyy;
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
    const su2double *gradFluxesDr = gradFluxesInt + 8*i;       /* nVar*nDim*i. */
    const su2double *gradFluxesDs = gradFluxesDr  + offDerivFluxes;

    /* Easier storage of the metric terms in this integration point.
       Also set the point where the divergence of the flux terms is
       stored for the current integration point. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
    su2double       *divFluxInt  = divFlux + 4*i;  /* nVar*i. */

    /* Compute the metric terms multiplied by the integration weight. Note that the
       first term in the metric terms is the Jacobian. */
    const su2double wDrdx = weights[i]*metricTerms[1];
    const su2double wDrdy = weights[i]*metricTerms[2];

    const su2double wDsdx = weights[i]*metricTerms[3];
    const su2double wDsdy = weights[i]*metricTerms[4];

    /* Compute the divergence of the fluxes, multiplied by the integration weight. */
    divFluxInt[0] = gradFluxesDr[0]*wDrdx + gradFluxesDs[0]*wDsdx
                  + gradFluxesDr[4]*wDrdy + gradFluxesDs[4]*wDsdy;
    divFluxInt[1] = gradFluxesDr[1]*wDrdx + gradFluxesDs[1]*wDsdx
                  + gradFluxesDr[5]*wDrdy + gradFluxesDs[5]*wDsdy;
    divFluxInt[2] = gradFluxesDr[2]*wDrdx + gradFluxesDs[2]*wDsdx
                  + gradFluxesDr[6]*wDrdy + gradFluxesDs[6]*wDsdy;
    divFluxInt[3] = gradFluxesDr[3]*wDrdx + gradFluxesDs[3]*wDsdx
                  + gradFluxesDr[7]*wDrdy + gradFluxesDs[7]*wDsdy;
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_NSSolver::ADER_DG_AliasedPredictorResidual_3D(CConfig           *config,
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
  su2double *gradFluxesInt = fluxesDOF + 15*nDOFs; /* nDOFs*nVar*nDim. */
  su2double *gradSolDOFs   = gradFluxesInt;
  su2double *divFlux       = work;

  /* Determine the offset between the r-derivatives and s-derivatives, which is
     also the offset between s- and t-derivatives, of the fluxes in the
     integration points and the solution in the DOFs. */
  const unsigned short offDerivSol    =  5*nDOFs;    /* nVar*nDOFs.     */
  const unsigned short offDerivFluxes = 15*nInt;     /* nVar*nInt*nDim. */

  /* Store the number of metric points per integration point/DOF for readability. */
  const unsigned short nMetricPerPoint = 10;  /* nDim*nDim + 1. */

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
    const su2double *solDOF   = sol         + 5*i;    /* i*nVar. */
    const su2double *solDOFDr = gradSolDOFs + 5*i;    /* i*nVar. */
    const su2double *solDOFDs = solDOFDr    + offDerivSol;
    const su2double *solDOFDt = solDOFDs    + offDerivSol;
    su2double       *fluxes   = fluxesDOF   + 15*i;   /* i*nVar*nDim. */

    /* Set the pointer to the grid velocities at the location of the
       solution DOFS. */
    const su2double *gridVel = elem->gridVelocitiesSolDOFs.data() + 3*i; /* nDim*i. */

    /* Compute the true value of the metric terms in this DOF. Note that in
       metricTerms the metric terms scaled by the Jacobian are stored. */
    const su2double *metricTerms = elem->metricTermsSolDOFs.data()
                                 + i*nMetricPerPoint;
    const su2double JacInv       = 1.0/metricTerms[0];

    const su2double drdx = JacInv*metricTerms[1];
    const su2double drdy = JacInv*metricTerms[2];
    const su2double drdz = JacInv*metricTerms[3];

    const su2double dsdx = JacInv*metricTerms[4];
    const su2double dsdy = JacInv*metricTerms[5];
    const su2double dsdz = JacInv*metricTerms[6];

    const su2double dtdx = JacInv*metricTerms[7];
    const su2double dtdy = JacInv*metricTerms[8];
    const su2double dtdz = JacInv*metricTerms[9];

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the metric
          terms in this DOF. ---*/
    su2double solGradCart[5][3];

    solGradCart[0][0] = solDOFDr[0]*drdx + solDOFDs[0]*dsdx + solDOFDt[0]*dtdx;
    solGradCart[1][0] = solDOFDr[1]*drdx + solDOFDs[1]*dsdx + solDOFDt[1]*dtdx;
    solGradCart[2][0] = solDOFDr[2]*drdx + solDOFDs[2]*dsdx + solDOFDt[2]*dtdx;
    solGradCart[3][0] = solDOFDr[3]*drdx + solDOFDs[3]*dsdx + solDOFDt[3]*dtdx;
    solGradCart[4][0] = solDOFDr[4]*drdx + solDOFDs[4]*dsdx + solDOFDt[4]*dtdx;

    solGradCart[0][1] = solDOFDr[0]*drdy + solDOFDs[0]*dsdy + solDOFDt[0]*dtdy;
    solGradCart[1][1] = solDOFDr[1]*drdy + solDOFDs[1]*dsdy + solDOFDt[1]*dtdy;
    solGradCart[2][1] = solDOFDr[2]*drdy + solDOFDs[2]*dsdy + solDOFDt[2]*dtdy;
    solGradCart[3][1] = solDOFDr[3]*drdy + solDOFDs[3]*dsdy + solDOFDt[3]*dtdy;
    solGradCart[4][1] = solDOFDr[4]*drdy + solDOFDs[4]*dsdy + solDOFDt[4]*dtdy;

    solGradCart[0][2] = solDOFDr[0]*drdz + solDOFDs[0]*dsdz + solDOFDt[0]*dtdz;
    solGradCart[1][2] = solDOFDr[1]*drdz + solDOFDs[1]*dsdz + solDOFDt[1]*dtdz;
    solGradCart[2][2] = solDOFDr[2]*drdz + solDOFDs[2]*dsdz + solDOFDt[2]*dtdz;
    solGradCart[3][2] = solDOFDr[3]*drdz + solDOFDs[3]*dsdz + solDOFDt[3]*dtdz;
    solGradCart[4][2] = solDOFDr[4]*drdz + solDOFDs[4]*dsdz + solDOFDt[4]*dtdz;

    /*--- Compute the velocities, pressure and laminar viscosity in this DOF. ---*/
    const su2double DensityInv   = 1.0/solDOF[0];
    const su2double u            = DensityInv*solDOF[1];
    const su2double v            = DensityInv*solDOF[2];
    const su2double w            = DensityInv*solDOF[3];
    const su2double TotalEnergy  = DensityInv*solDOF[4];
    const su2double StaticEnergy = TotalEnergy - 0.5*(u*u + v*v + w*w);

    FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
    const su2double Pressure     = FluidModel->GetPressure();
    const su2double ViscosityLam = FluidModel->GetLaminarViscosity();

    /*--- Compute the Cartesian gradients of the velocities and
          static energy in this DOF. ---*/
    const su2double dudx = DensityInv*(solGradCart[1][0] - u*solGradCart[0][0]);
    const su2double dudy = DensityInv*(solGradCart[1][1] - u*solGradCart[0][1]);
    const su2double dudz = DensityInv*(solGradCart[1][2] - u*solGradCart[0][2]);

    const su2double dvdx = DensityInv*(solGradCart[2][0] - v*solGradCart[0][0]);
    const su2double dvdy = DensityInv*(solGradCart[2][1] - v*solGradCart[0][1]);
    const su2double dvdz = DensityInv*(solGradCart[2][2] - v*solGradCart[0][2]);

    const su2double dwdx = DensityInv*(solGradCart[3][0] - w*solGradCart[0][0]);
    const su2double dwdy = DensityInv*(solGradCart[3][1] - w*solGradCart[0][1]);
    const su2double dwdz = DensityInv*(solGradCart[3][2] - w*solGradCart[0][2]);

    const su2double dedx = DensityInv*(solGradCart[4][0] - TotalEnergy*solGradCart[0][0])
                         - u*dudx - v*dvdx - w*dwdx;
    const su2double dedy = DensityInv*(solGradCart[4][1] - TotalEnergy*solGradCart[0][1])
                         - u*dudy - v*dvdy - w*dwdy;
    const su2double dedz = DensityInv*(solGradCart[4][2] - TotalEnergy*solGradCart[0][2])
                         - u*dudz - v*dvdz - w*dwdz;

    /*--- Compute the eddy viscosity, if needed, and the total viscosity. ---*/
    su2double ViscosityTurb = 0.0;
    if( SGSModelUsed )
      ViscosityTurb = SGSModel->ComputeEddyViscosity_3D(solDOF[0], dudx, dudy, dudz,
                                                        dvdx, dvdy, dvdz, dwdx,
                                                        dwdy, dwdz, lenScale,
                                                        elem->wallDistanceSolDOFs[i]);
    const su2double Viscosity = ViscosityLam + ViscosityTurb;

    /* Compute the total thermal conductivity divided by Cv. */
    const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                            + ViscosityTurb*factHeatFlux_Turb;

    /*--- Set the value of the second viscosity and compute the divergence
          term in the viscous normal stresses. ---*/
    const su2double lambda     = -TWO3*Viscosity;
    const su2double lamDivTerm =  lambda*(dudx + dvdy + dwdz);

    /*--- Compute the viscous stress tensor. ---*/
    const su2double tauxx = 2.0*Viscosity*dudx + lamDivTerm;
    const su2double tauyy = 2.0*Viscosity*dvdy + lamDivTerm;
    const su2double tauzz = 2.0*Viscosity*dwdz + lamDivTerm;

    const su2double tauxy = Viscosity*(dudy + dvdx);
    const su2double tauxz = Viscosity*(dudz + dwdx);
    const su2double tauyz = Viscosity*(dvdz + dwdy);

    /*--- Compute the Cartesian fluxes in this DOF. First the fluxes in
          x-direction are stored, followed by the fluxes in y-direction
          and the fluxes in z-direction. ---*/
    const su2double uRel = u - gridVel[0];
    const su2double vRel = v - gridVel[1];
    const su2double wRel = w - gridVel[2];

    fluxes[0] = solDOF[0]*uRel;
    fluxes[1] = solDOF[1]*uRel + Pressure - tauxx;
    fluxes[2] = solDOF[2]*uRel - tauxy;
    fluxes[3] = solDOF[3]*uRel - tauxz;
    fluxes[4] = solDOF[4]*uRel + Pressure*u - kOverCv*dedx - u*tauxx - v*tauxy - w*tauxz;

    fluxes[5] = solDOF[0]*vRel;
    fluxes[6] = solDOF[1]*vRel - tauxy;
    fluxes[7] = solDOF[2]*vRel + Pressure - tauyy;
    fluxes[8] = solDOF[3]*vRel - tauyz;
    fluxes[9] = solDOF[4]*vRel + Pressure*v - kOverCv*dedy - u*tauxy - v*tauyy - w*tauyz;

    fluxes[10] = solDOF[0]*wRel;
    fluxes[11] = solDOF[1]*wRel - tauxz;
    fluxes[12] = solDOF[2]*wRel - tauyz;
    fluxes[13] = solDOF[3]*wRel + Pressure - tauzz;
    fluxes[14] = solDOF[4]*wRel + Pressure*w - kOverCv*dedz - u*tauxz - v*tauyz - w*tauzz;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Compute the divergence of the fluxes in the integration    ---*/
  /*---         points and distribute the divergence terms to the DOFs.    ---*/
  /*--------------------------------------------------------------------------*/

  /* Compute the derivatives of the Cartesian fluxes w.r.t. the parametric
     coordinates in the integration points. */
  DenseMatrixProduct(nInt*nDim, nVar*nDim, nDOFs, matDerBasisInt, fluxesDOF,
                     gradFluxesInt);

  /*--- Loop over the integration points to compute the divergence of the inviscid
        fluxes in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the data of the derivatives of the
       fluxes of this integration point starts. */
    const su2double *gradFluxesDr = gradFluxesInt + 15*i;       /* nVar*nDim*i. */
    const su2double *gradFluxesDs = gradFluxesDr  + offDerivFluxes;
    const su2double *gradFluxesDt = gradFluxesDs  + offDerivFluxes;

    /* Easier storage of the metric terms in this integration point.
       Also set the point where the divergence of the flux terms is
       stored for the current integration point. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
    su2double       *divFluxInt  = divFlux + 5*i;  /* nVar*i. */

    /* Compute the metric terms multiplied by the integration weight. Note that the
       first term in the metric terms is the Jacobian. */
    const su2double wDrdx = weights[i]*metricTerms[1];
    const su2double wDrdy = weights[i]*metricTerms[2];
    const su2double wDrdz = weights[i]*metricTerms[3];

    const su2double wDsdx = weights[i]*metricTerms[4];
    const su2double wDsdy = weights[i]*metricTerms[5];
    const su2double wDsdz = weights[i]*metricTerms[6];

    const su2double wDtdx = weights[i]*metricTerms[7];
    const su2double wDtdy = weights[i]*metricTerms[8];
    const su2double wDtdz = weights[i]*metricTerms[9];

    /* Compute the divergence of the fluxes, multiplied by the integration weight. */
    divFluxInt[0] = gradFluxesDr[0] *wDrdx + gradFluxesDs[0] *wDsdx + gradFluxesDt[0] *wDtdx
                  + gradFluxesDr[5] *wDrdy + gradFluxesDs[5] *wDsdy + gradFluxesDt[5] *wDtdy
                  + gradFluxesDr[10]*wDrdz + gradFluxesDs[10]*wDsdz + gradFluxesDt[10]*wDtdz;
    divFluxInt[1] = gradFluxesDr[1] *wDrdx + gradFluxesDs[1] *wDsdx + gradFluxesDt[1] *wDtdx
                  + gradFluxesDr[6] *wDrdy + gradFluxesDs[6] *wDsdy + gradFluxesDt[6] *wDtdy
                  + gradFluxesDr[11]*wDrdz + gradFluxesDs[11]*wDsdz + gradFluxesDt[11]*wDtdz;
    divFluxInt[2] = gradFluxesDr[2] *wDrdx + gradFluxesDs[2] *wDsdx + gradFluxesDt[2] *wDtdx
                  + gradFluxesDr[7] *wDrdy + gradFluxesDs[7] *wDsdy + gradFluxesDt[7] *wDtdy
                  + gradFluxesDr[12]*wDrdz + gradFluxesDs[12]*wDsdz + gradFluxesDt[12]*wDtdz;
    divFluxInt[3] = gradFluxesDr[3] *wDrdx + gradFluxesDs[3] *wDsdx + gradFluxesDt[3] *wDtdx
                  + gradFluxesDr[8] *wDrdy + gradFluxesDs[8] *wDsdy + gradFluxesDt[8] *wDtdy
                  + gradFluxesDr[13]*wDrdz + gradFluxesDs[13]*wDsdz + gradFluxesDt[13]*wDtdz;
    divFluxInt[4] = gradFluxesDr[4] *wDrdx + gradFluxesDs[4] *wDsdx + gradFluxesDt[4] *wDtdx
                  + gradFluxesDr[9] *wDrdy + gradFluxesDs[9] *wDsdy + gradFluxesDt[9] *wDtdy
                  + gradFluxesDr[14]*wDrdz + gradFluxesDs[14]*wDsdz + gradFluxesDt[14]*wDtdz;
  }

  /* Compute the residual in the DOFs, which is the matrix product of
     basisFunctionsIntTrans and divFlux. */
  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_NSSolver::ADER_DG_NonAliasedPredictorResidual_2D(CConfig           *config,
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
  const unsigned short offPoint = 12*max(nInt, nDOFs);  /*nVar*(nDim+1)*max(..). */
  su2double *secDerSolInt       = solAndGradInt + offPoint;

  /* Determine the offset between the solution variables and the r-derivatives,
     which is also the offset between the r- and s-derivatives. Also determine
     the offset between the data of the second derivatives. */
  const unsigned short offDeriv    = 4*nInt;       /* nVar*nInt. */
  const unsigned short off2ndDeriv = 2*offDeriv;   /* nDim*offDeriv. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 5;  /* nDim*nDim + 1. */

  /* Store the number of additional metric points per integration point, which
     are needed to compute the second derivatives. These terms take the
     non-constant metric into account. */
  const unsigned short nMetric2ndDerPerPoint = 6; /*nDim*(nDim + nDim*(nDim-1)/2). */

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
  for(unsigned short i=0; i<nDOFs; ++i) {
    const su2double *dSolDr = gradSolDOFs + 4*i;     /* nVar*i. */
    const su2double *dSolDs = dSolDr      + 4*nDOFs; /* nVar*nDOFs. */
    su2double *firstDer = firstDerSolDOFs + 8*i;     /* nVar*nDim*i. */

    memcpy(firstDer,   dSolDr, 4*sizeof(su2double));
    memcpy(firstDer+4, dSolDs, 4*sizeof(su2double));
  }

  /* Compute the second derivatives w.r.t. the parametric coordinates
     in the integration points. */
  DenseMatrixProduct(nInt*nDim, nVar*nDim, nDOFs, matDerBasisInt,
                     firstDerSolDOFs, secDerSolInt);

  /* Compute the solution and the derivatives w.r.t. the parametric coordinates
     in the integration points. The first argument is nInt*(nDim+1). */
  DenseMatrixProduct(nInt*3, nVar, nDOFs, matBasisInt, sol, solAndGradInt);

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Determine from the solution and gradients the divergence   ---*/
  /*---         of the fluxes in the integration points.                   ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the integration points to compute the divergence of the fluxes
        in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the solution and gradient data
       of this integration point starts. */
    const su2double *sol   = solAndGradInt + 4*i; /* nVar*i. */
    const su2double *solDr = sol   + offDeriv;
    const su2double *solDs = solDr + offDeriv;

    /*--- Compute the velocities, pressure and total enthalpy
          in this integration point. ---*/
    const su2double DensityInv   = 1.0/sol[0];
    const su2double u            = DensityInv*sol[1];
    const su2double v            = DensityInv*sol[2];
    const su2double kinEnergy    = 0.5*(u*u + v*v);
    const su2double TotalEnergy  = DensityInv*sol[3];
    const su2double StaticEnergy = TotalEnergy - kinEnergy;

    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
    const su2double Pressure = FluidModel->GetPressure();
    const su2double Htot     = DensityInv*(sol[3] + Pressure);

    /* Easier storage of the grid velocities in this integration point. */
    const su2double *gridVel = elem->gridVelocities.data() + 2*i; /* nDim*i. */

    /* Compute the laminar viscosity and its derivative w.r.t. temperature. */
    const su2double ViscosityLam = FluidModel->GetLaminarViscosity();
    const su2double dViscLamdT   = FluidModel->GetdmudT_rho();

    /* Compute the true metric terms. Note in metricTerms the actual metric
       terms multiplied by the Jacobian are stored. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
    const su2double Jac          = metricTerms[0];
    const su2double JacInv       = 1.0/Jac;

    const su2double drdx = JacInv*metricTerms[1];
    const su2double drdy = JacInv*metricTerms[2];
    const su2double dsdx = JacInv*metricTerms[3];
    const su2double dsdy = JacInv*metricTerms[4];

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the
          metric terms in this integration point. ---*/
    su2double solGradCart[4][2];

    solGradCart[0][0] = solDr[0]*drdx + solDs[0]*dsdx;
    solGradCart[1][0] = solDr[1]*drdx + solDs[1]*dsdx;
    solGradCart[2][0] = solDr[2]*drdx + solDs[2]*dsdx;
    solGradCart[3][0] = solDr[3]*drdx + solDs[3]*dsdx;

    solGradCart[0][1] = solDr[0]*drdy + solDs[0]*dsdy;
    solGradCart[1][1] = solDr[1]*drdy + solDs[1]*dsdy;
    solGradCart[2][1] = solDr[2]*drdy + solDs[2]*dsdy;
    solGradCart[3][1] = solDr[3]*drdy + solDs[3]*dsdy;

    /* Easier storage of the location where the data of second derivatives of
       this integration point starts as well as the necessary additional metric
       terms to compute the Cartesian second derivatives. */
    const su2double *sol2ndDerDr       = secDerSolInt + 8*i;  /* nVar*nDim*i. */
    const su2double *sol2ndDerDs       = sol2ndDerDr + off2ndDeriv;
    const su2double *metricTerms2ndDer = elem->metricTerms2ndDer.data()
                                       + i*nMetric2ndDerPerPoint;

    /*--- Compute the Cartesian second derivatives of the independent solution
          variables from the gradients and second derivatives in parametric
          coordinates and the metric terms and its derivatives w.r.t. the
          parametric coordinates. ---*/
    su2double sol2ndDerCart[4][2][2];

    sol2ndDerCart[0][0][0] = sol2ndDerDr[0]*drdx*drdx + sol2ndDerDs[4]*dsdx*dsdx
                           + (sol2ndDerDr[4] + sol2ndDerDs[0])*drdx*dsdx
                           + solDr[0]*metricTerms2ndDer[0] + solDs[0]*metricTerms2ndDer[1];
    sol2ndDerCart[1][0][0] = sol2ndDerDr[1]*drdx*drdx + sol2ndDerDs[5]*dsdx*dsdx
                           + (sol2ndDerDr[5] + sol2ndDerDs[1])*drdx*dsdx
                           + solDr[1]*metricTerms2ndDer[0] + solDs[1]*metricTerms2ndDer[1];
    sol2ndDerCart[2][0][0] = sol2ndDerDr[2]*drdx*drdx + sol2ndDerDs[6]*dsdx*dsdx
                           + (sol2ndDerDr[6] + sol2ndDerDs[2])*drdx*dsdx
                           + solDr[2]*metricTerms2ndDer[0] + solDs[2]*metricTerms2ndDer[1];
    sol2ndDerCart[3][0][0] = sol2ndDerDr[3]*drdx*drdx + sol2ndDerDs[7]*dsdx*dsdx
                           + (sol2ndDerDr[7] + sol2ndDerDs[3])*drdx*dsdx
                           + solDr[3]*metricTerms2ndDer[0] + solDs[3]*metricTerms2ndDer[1];

    sol2ndDerCart[0][1][0] = sol2ndDerDr[0]*drdx*drdy + sol2ndDerDs[4]*dsdx*dsdy
                           + 0.5*(sol2ndDerDr[4] + sol2ndDerDs[0])*(drdx*dsdy + dsdx*drdy)
                           + solDr[0]*metricTerms2ndDer[2] + solDs[0]*metricTerms2ndDer[3];
    sol2ndDerCart[1][1][0] = sol2ndDerDr[1]*drdx*drdy + sol2ndDerDs[5]*dsdx*dsdy
                           + 0.5*(sol2ndDerDr[5] + sol2ndDerDs[1])*(drdx*dsdy + dsdx*drdy)
                           + solDr[1]*metricTerms2ndDer[2] + solDs[1]*metricTerms2ndDer[3];
    sol2ndDerCart[2][1][0] = sol2ndDerDr[2]*drdx*drdy + sol2ndDerDs[6]*dsdx*dsdy
                           + 0.5*(sol2ndDerDr[6] + sol2ndDerDs[2])*(drdx*dsdy + dsdx*drdy)
                           + solDr[2]*metricTerms2ndDer[2] + solDs[2]*metricTerms2ndDer[3];
    sol2ndDerCart[3][1][0] = sol2ndDerDr[3]*drdx*drdy + sol2ndDerDs[7]*dsdx*dsdy
                           + 0.5*(sol2ndDerDr[7] + sol2ndDerDs[3])*(drdx*dsdy + dsdx*drdy)
                           + solDr[3]*metricTerms2ndDer[2] + solDs[3]*metricTerms2ndDer[3];

    sol2ndDerCart[0][0][1] = sol2ndDerCart[0][1][0];
    sol2ndDerCart[1][0][1] = sol2ndDerCart[1][1][0];
    sol2ndDerCart[2][0][1] = sol2ndDerCart[2][1][0];
    sol2ndDerCart[3][0][1] = sol2ndDerCart[3][1][0];

    sol2ndDerCart[0][1][1] = sol2ndDerDr[0]*drdy*drdy + sol2ndDerDs[4]*dsdy*dsdy
                           + (sol2ndDerDr[4] + sol2ndDerDs[0])*drdy*dsdy
                           + solDr[0]*metricTerms2ndDer[4] + solDs[0]*metricTerms2ndDer[5];
    sol2ndDerCart[1][1][1] = sol2ndDerDr[1]*drdy*drdy + sol2ndDerDs[5]*dsdy*dsdy
                           + (sol2ndDerDr[5] + sol2ndDerDs[1])*drdy*dsdy
                           + solDr[1]*metricTerms2ndDer[4] + solDs[1]*metricTerms2ndDer[5];
    sol2ndDerCart[2][1][1] = sol2ndDerDr[2]*drdy*drdy + sol2ndDerDs[6]*dsdy*dsdy
                           + (sol2ndDerDr[6] + sol2ndDerDs[2])*drdy*dsdy
                           + solDr[2]*metricTerms2ndDer[4] + solDs[2]*metricTerms2ndDer[5];
    sol2ndDerCart[3][1][1] = sol2ndDerDr[3]*drdy*drdy + sol2ndDerDs[7]*dsdy*dsdy
                           + (sol2ndDerDr[7] + sol2ndDerDs[3])*drdy*dsdy
                           + solDr[3]*metricTerms2ndDer[4] + solDs[3]*metricTerms2ndDer[5];

    /*--- Compute the Cartesian gradients of the pressure, velocity components,
          static energy and dynamic viscosity. ---*/
    const su2double dpdx = Gamma_Minus_One*(solGradCart[3][0] + kinEnergy*solGradCart[0][0]
                         -                  u*solGradCart[1][0] - v*solGradCart[2][0]);
    const su2double dpdy = Gamma_Minus_One*(solGradCart[3][1] + kinEnergy*solGradCart[0][1]
                         -                  u*solGradCart[1][1] - v*solGradCart[2][1]);

    const su2double dudx = DensityInv*(solGradCart[1][0] - u*solGradCart[0][0]);
    const su2double dudy = DensityInv*(solGradCart[1][1] - u*solGradCart[0][1]);
    const su2double dvdx = DensityInv*(solGradCart[2][0] - v*solGradCart[0][0]);
    const su2double dvdy = DensityInv*(solGradCart[2][1] - v*solGradCart[0][1]);

    const su2double dedx = DensityInv*(solGradCart[3][0] - TotalEnergy*solGradCart[0][0])
                         - u*dudx - v*dvdx;
    const su2double dedy = DensityInv*(solGradCart[3][1] - TotalEnergy*solGradCart[0][1])
                         - u*dudy - v*dvdy;

    const su2double dViscLamdx = CvInv*dedx*dViscLamdT;
    const su2double dViscLamdy = CvInv*dedy*dViscLamdT;

    /*--- Compute the second derivatives of the velocity components. ---*/
    const su2double d2udxx = DensityInv*(sol2ndDerCart[1][0][0] - u*sol2ndDerCart[0][0][0]
                           +    2.0*DensityInv*solGradCart[0][0]*(u*solGradCart[0][0] - solGradCart[1][0]));
    const su2double d2udyy = DensityInv*(sol2ndDerCart[1][1][1] - u*sol2ndDerCart[0][1][1]
                           +    2.0*DensityInv*solGradCart[0][1]*(u*solGradCart[0][1] - solGradCart[1][1]));
    const su2double d2udxy = DensityInv*(sol2ndDerCart[1][0][1] - u*sol2ndDerCart[0][0][1]
                           +       DensityInv*(solGradCart[0][0]*(u*solGradCart[0][1] - solGradCart[1][1])
                           +                   solGradCart[0][1]*(u*solGradCart[0][0] - solGradCart[1][0])));

    const su2double d2vdxx = DensityInv*(sol2ndDerCart[2][0][0] - v*sol2ndDerCart[0][0][0]
                           +    2.0*DensityInv*solGradCart[0][0]*(v*solGradCart[0][0] - solGradCart[2][0]));
    const su2double d2vdyy = DensityInv*(sol2ndDerCart[2][1][1] - v*sol2ndDerCart[0][1][1]
                           +    2.0*DensityInv*solGradCart[0][1]*(v*solGradCart[0][1] - solGradCart[2][1]));
    const su2double d2vdxy = DensityInv*(sol2ndDerCart[2][0][1] - v*sol2ndDerCart[0][0][1]
                           +       DensityInv*(solGradCart[0][0]*(v*solGradCart[0][1] - solGradCart[2][1])
                           +                   solGradCart[0][1]*(v*solGradCart[0][0] - solGradCart[2][0])));

    /*--- Compute the second derivatives of the static energy. Note that this
          term appears in the heat flux and therefore only the pure second
          derivatives are needed. Hence, the cross-derivatives are omitted. ---*/
    const su2double d2edxx = DensityInv*(sol2ndDerCart[3][0][0] - TotalEnergy*sol2ndDerCart[0][0][0]
                           +    2.0*DensityInv*solGradCart[0][0]*(TotalEnergy*solGradCart[0][0] - solGradCart[3][0]))
                           -    u*d2udxx - dudx*dudx - v*d2vdxx - dvdx*dvdx;
    const su2double d2edyy = DensityInv*(sol2ndDerCart[3][1][1] - TotalEnergy*sol2ndDerCart[0][1][1]
                           +    2.0*DensityInv*solGradCart[0][1]*(TotalEnergy*solGradCart[0][1] - solGradCart[3][1]))
                           -    u*d2udyy - dudy*dudy - v*d2vdyy - dvdy*dvdy;

    /*--- If an SGS model is used the eddy viscosity and its spatial
          derivatives must be computed. ---*/
    su2double ViscosityTurb = 0.0;
    su2double dViscTurbdx = 0.0, dViscTurbdy = 0.0;

    if( SGSModelUsed ) {
      const su2double dist = elem->wallDistance[i];
      ViscosityTurb = SGSModel->ComputeEddyViscosity_2D(sol[0], dudx, dudy, dvdx,
                                                        dvdy, lenScale, dist);

      SGSModel->ComputeGradEddyViscosity_2D(sol[0], solGradCart[0][0], solGradCart[0][1],
                                            dudx, dudy, dvdx, dvdy, d2udxx, d2udyy,
                                            d2udxy, d2vdxx, d2vdyy, d2vdxy, lenScale,
                                            dist, dViscTurbdx, dViscTurbdy);
    }

    /*--- Compute the total viscosity, the total heat conductivity and their
          gradients. Note that the heat conductivity is divided by the Cv,
          because gradients of internal energy are computed and not temperature. ---*/
    const su2double Viscosity = ViscosityLam + ViscosityTurb;
    const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                            + ViscosityTurb*factHeatFlux_Turb;

    const su2double dViscDx = dViscLamdx + dViscTurbdx;
    const su2double dViscDy = dViscLamdy + dViscTurbdy;

    const su2double dkOverCvdx = dViscLamdx *factHeatFlux_Lam
                               + dViscTurbdx*factHeatFlux_Turb;
    const su2double dkOverCvdy = dViscLamdy *factHeatFlux_Lam
                               + dViscTurbdy*factHeatFlux_Turb;

    /* Abbreviations, which make it easier to compute the divergence term. */
    const su2double abv1 = solGradCart[1][0] + solGradCart[2][1];
    const su2double abv2 = u*solGradCart[0][0] + v*solGradCart[0][1];
    const su2double abv3 = u*(solGradCart[3][0] + dpdx) + v*(solGradCart[3][1] + dpdy);
    const su2double abv4 = dudx + dvdy;

    /*--- Compute the divergence of the grid velocity.
          SET TO ZERO FOR NOW. THIS IS NOT CORRECT!!!!. ---*/
    const su2double divGridVel = 0.0;

    /*--- Set the pointer to store the divergence terms for this integration
          point and compute these terms, multiplied by the integration weight
          and Jacobian. ---*/
    const su2double weightJac = weights[i]*Jac;
    su2double *divFluxInt     = divFlux + 4*i;   /* nVar*i. */

    divFluxInt[0] = weightJac*(abv1 - sol[0]*divGridVel
                  -            gridVel[0]*solGradCart[0][0] - gridVel[1]*solGradCart[0][1]);
    divFluxInt[1] = weightJac*(dpdx + u*(abv1-abv2) - lambdaOverMu*abv4*dViscDx
                  +            u*solGradCart[1][0] + v*solGradCart[1][1]
                  -            lambdaOverMu*Viscosity*(d2udxx + d2vdxy)
                  -            Viscosity*(2.0*d2udxx + d2udyy + d2vdxy)
                  -            2.0*dViscDx*dudx - dViscDy*(dudy+dvdx)
                  -            sol[1]*divGridVel
                  -            gridVel[0]*solGradCart[1][0] - gridVel[1]*solGradCart[1][1]);
    divFluxInt[2] = weightJac*(dpdy + v*(abv1-abv2) - lambdaOverMu*abv4*dViscDy
                  +            u*solGradCart[2][0] + v*solGradCart[2][1]
                  -            lambdaOverMu*Viscosity*(d2udxy + d2vdyy)
                  -            Viscosity*(2.0*d2vdyy + d2vdxx + d2udxy)
                  -            dViscDx*(dudy + dvdx) - 2.0*dViscDy*dvdy
                  -            sol[2]*divGridVel
                  -            gridVel[0]*solGradCart[2][0] - gridVel[1]*solGradCart[2][1]);
    divFluxInt[3] = weightJac*(abv3 + Htot*(abv1 - abv2)
                  -            abv4*lambdaOverMu*(Viscosity*abv4 + u*dViscDx + v*dViscDy)
                  -            dkOverCvdx*dedx - dkOverCvdy*dedy - kOverCv*(d2edxx + d2edyy)
                  -            (Viscosity*dudx + u*dViscDx)*2.0*dudx
                  -            (Viscosity*dvdy + v*dViscDy)*2.0*dvdy
                  -            (Viscosity*dudy + u*dViscDy + Viscosity*dvdx + v*dViscDx)*(dudy + dvdx)
                  -            Viscosity*u*(d2udxx+d2udyy + (1.0+lambdaOverMu)*(d2udxx+d2vdxy))
                  -            Viscosity*v*(d2vdxx+d2vdyy + (1.0+lambdaOverMu)*(d2udxy+d2vdyy))
                  -            sol[3]*divGridVel
                  -            gridVel[0]*solGradCart[3][0] - gridVel[1]*solGradCart[3][1]);
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Compute the residual in the DOFs, which is the matrix      ---*/
  /*---         product of basisFunctionsIntTrans and divFlux.             ---*/
  /*--------------------------------------------------------------------------*/

  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_NSSolver::ADER_DG_NonAliasedPredictorResidual_3D(CConfig           *config,
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
  const unsigned short offPoint = 20*max(nInt, nDOFs);  /*nVar*(nDim+1)*max(..). */
  su2double *secDerSolInt       = solAndGradInt + offPoint;

  /* Determine the offset between the solution variables and the r-derivatives,
     which is also the offset between the r- and s-derivatives and the offset
     between s- and t-derivatives. Also determine the offset between the data
     of the second derivatives. */
  const unsigned short offDeriv    = 5*nInt;     /* nVar*nInt. */
  const unsigned short off2ndDeriv = 3*offDeriv; /* nDim*offDeriv. */

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 10;  /* nDim*nDim + 1. */

  /* Store the number of additional metric points per integration point, which
     are needed to compute the second derivatives. These terms take the
     non-constant metric into account. */
  const unsigned short nMetric2ndDerPerPoint = 18; /*nDim*(nDim + nDim*(nDim-1)/2). */

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
  for(unsigned short i=0; i<nDOFs; ++i) {
    const su2double *dSolDr = gradSolDOFs +  5*i;     /* nVar*i. */
    const su2double *dSolDs = dSolDr      +  5*nDOFs; /* nVar*nDOFs. */
    const su2double *dSolDt = dSolDs      +  5*nDOFs; /* nVar*nDOFs. */
    su2double *firstDer = firstDerSolDOFs + 15*i;     /* nVar*nDim*i. */

    memcpy(firstDer,    dSolDr, 5*sizeof(su2double));
    memcpy(firstDer+ 5, dSolDs, 5*sizeof(su2double));
    memcpy(firstDer+10, dSolDt, 5*sizeof(su2double));
  }

  /* Compute the second derivatives w.r.t. the parametric coordinates
     in the integration points. */
  DenseMatrixProduct(nInt*nDim, nVar*nDim, nDOFs, matDerBasisInt,
                     firstDerSolDOFs, secDerSolInt);

  /* Compute the solution and the derivatives w.r.t. the parametric coordinates
     in the integration points. The first argument is nInt*(nDim+1). */
  DenseMatrixProduct(nInt*4, nVar, nDOFs, matBasisInt, sol, solAndGradInt);

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Determine from the solution and gradients the divergence   ---*/
  /*---         of the fluxes in the integration points.                   ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the integration points to compute the divergence of the fluxes
        in these integration points, multiplied by the integration weight. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the location where the solution and gradient data
       of this integration point starts. */
    const su2double *sol   = solAndGradInt + 5*i; /* nVar*i. */
    const su2double *solDr = sol   + offDeriv;
    const su2double *solDs = solDr + offDeriv;
    const su2double *solDt = solDs + offDeriv;

    /*--- Compute the velocities, pressure and total enthalpy
          in this integration point. ---*/
    const su2double DensityInv   = 1.0/sol[0];
    const su2double u            = DensityInv*sol[1];
    const su2double v            = DensityInv*sol[2];
    const su2double w            = DensityInv*sol[3];
    const su2double kinEnergy    = 0.5*(u*u + v*v + w*w);
    const su2double TotalEnergy  = DensityInv*sol[4];
    const su2double StaticEnergy = TotalEnergy - kinEnergy;

    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
    const su2double Pressure = FluidModel->GetPressure();
    const su2double Htot     = DensityInv*(sol[4] + Pressure);

    /* Easier storage of the grid velocities in this integration point. */
    const su2double *gridVel = elem->gridVelocities.data() + 3*i; /* nDim*i. */

    /* Compute the laminar viscosity and its derivative w.r.t. temperature. */
    const su2double ViscosityLam = FluidModel->GetLaminarViscosity();
    const su2double dViscLamdT   = FluidModel->GetdmudT_rho();

    /* Compute the true metric terms. Note in metricTerms the actual metric
       terms multiplied by the Jacobian are stored. */
    const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
    const su2double Jac          = metricTerms[0];
    const su2double JacInv       = 1.0/Jac;

    const su2double drdx = JacInv*metricTerms[1];
    const su2double drdy = JacInv*metricTerms[2];
    const su2double drdz = JacInv*metricTerms[3];

    const su2double dsdx = JacInv*metricTerms[4];
    const su2double dsdy = JacInv*metricTerms[5];
    const su2double dsdz = JacInv*metricTerms[6];

    const su2double dtdx = JacInv*metricTerms[7];
    const su2double dtdy = JacInv*metricTerms[8];
    const su2double dtdz = JacInv*metricTerms[9];

    /*--- Compute the Cartesian gradients of the independent solution
          variables from the gradients in parametric coordinates and the
          metric terms in this integration point. ---*/
    su2double solGradCart[5][3];

    solGradCart[0][0] = solDr[0]*drdx + solDs[0]*dsdx + solDt[0]*dtdx;
    solGradCart[1][0] = solDr[1]*drdx + solDs[1]*dsdx + solDt[1]*dtdx;
    solGradCart[2][0] = solDr[2]*drdx + solDs[2]*dsdx + solDt[2]*dtdx;
    solGradCart[3][0] = solDr[3]*drdx + solDs[3]*dsdx + solDt[3]*dtdx;
    solGradCart[4][0] = solDr[4]*drdx + solDs[4]*dsdx + solDt[4]*dtdx;

    solGradCart[0][1] = solDr[0]*drdy + solDs[0]*dsdy + solDt[0]*dtdy;
    solGradCart[1][1] = solDr[1]*drdy + solDs[1]*dsdy + solDt[1]*dtdy;
    solGradCart[2][1] = solDr[2]*drdy + solDs[2]*dsdy + solDt[2]*dtdy;
    solGradCart[3][1] = solDr[3]*drdy + solDs[3]*dsdy + solDt[3]*dtdy;
    solGradCart[4][1] = solDr[4]*drdy + solDs[4]*dsdy + solDt[4]*dtdy;

    solGradCart[0][2] = solDr[0]*drdz + solDs[0]*dsdz + solDt[0]*dtdz;
    solGradCart[1][2] = solDr[1]*drdz + solDs[1]*dsdz + solDt[1]*dtdz;
    solGradCart[2][2] = solDr[2]*drdz + solDs[2]*dsdz + solDt[2]*dtdz;
    solGradCart[3][2] = solDr[3]*drdz + solDs[3]*dsdz + solDt[3]*dtdz;
    solGradCart[4][2] = solDr[4]*drdz + solDs[4]*dsdz + solDt[4]*dtdz;

    /* Easier storage of the location where the data of second derivatives of
       this integration point starts as well as the necessary additional metric
       terms to compute the Cartesian second derivatives. */
    const su2double *sol2ndDerDr       = secDerSolInt + 15*i;  /* nVar*nDim*i. */
    const su2double *sol2ndDerDs       = sol2ndDerDr + off2ndDeriv;
    const su2double *sol2ndDerDt       = sol2ndDerDs + off2ndDeriv;
    const su2double *metricTerms2ndDer = elem->metricTerms2ndDer.data()
                                       + i*nMetric2ndDerPerPoint;

    /*--- Compute the Cartesian second derivatives of the independent solution
          variables from the gradients and second derivatives in parametric
          coordinates and the metric terms and its derivatives w.r.t. the
          parametric coordinates. ---*/
    su2double sol2ndDerCart[5][3][3];

    sol2ndDerCart[0][0][0] = sol2ndDerDr[0]*drdx*drdx + sol2ndDerDs[5]*dsdx*dsdx + sol2ndDerDs[10]*dtdx*dtdx
                           + (sol2ndDerDr[5]  + sol2ndDerDs[0])*drdx*dsdx + (sol2ndDerDr[10] + sol2ndDerDt[0])*drdx*dtdx
                           + (sol2ndDerDs[10] + sol2ndDerDt[5])*dsdx*dtdx
                           + solDr[0]*metricTerms2ndDer[0] + solDs[0]*metricTerms2ndDer[1] + solDt[0]*metricTerms2ndDer[2];
    sol2ndDerCart[1][0][0] = sol2ndDerDr[1]*drdx*drdx + sol2ndDerDs[6]*dsdx*dsdx + sol2ndDerDs[11]*dtdx*dtdx
                           + (sol2ndDerDr[6]  + sol2ndDerDs[1])*drdx*dsdx + (sol2ndDerDr[11] + sol2ndDerDt[1])*drdx*dtdx
                           + (sol2ndDerDs[11] + sol2ndDerDt[6])*dsdx*dtdx
                           + solDr[1]*metricTerms2ndDer[0] + solDs[1]*metricTerms2ndDer[1] + solDt[1]*metricTerms2ndDer[2];
    sol2ndDerCart[2][0][0] = sol2ndDerDr[2]*drdx*drdx + sol2ndDerDs[7]*dsdx*dsdx + sol2ndDerDs[12]*dtdx*dtdx
                           + (sol2ndDerDr[7]  + sol2ndDerDs[2])*drdx*dsdx + (sol2ndDerDr[12] + sol2ndDerDt[2])*drdx*dtdx
                           + (sol2ndDerDs[12] + sol2ndDerDt[7])*dsdx*dtdx
                           + solDr[2]*metricTerms2ndDer[0] + solDs[2]*metricTerms2ndDer[1] + solDt[2]*metricTerms2ndDer[2];
    sol2ndDerCart[3][0][0] = sol2ndDerDr[3]*drdx*drdx + sol2ndDerDs[8]*dsdx*dsdx + sol2ndDerDs[13]*dtdx*dtdx
                           + (sol2ndDerDr[8]  + sol2ndDerDs[3])*drdx*dsdx + (sol2ndDerDr[13] + sol2ndDerDt[3])*drdx*dtdx
                           + (sol2ndDerDs[13] + sol2ndDerDt[8])*dsdx*dtdx
                           + solDr[3]*metricTerms2ndDer[0] + solDs[3]*metricTerms2ndDer[1] + solDt[3]*metricTerms2ndDer[2];
    sol2ndDerCart[4][0][0] = sol2ndDerDr[4]*drdx*drdx + sol2ndDerDs[9]*dsdx*dsdx + sol2ndDerDs[14]*dtdx*dtdx
                           + (sol2ndDerDr[9]  + sol2ndDerDs[4])*drdx*dsdx + (sol2ndDerDr[14] + sol2ndDerDt[4])*drdx*dtdx
                           + (sol2ndDerDs[14] + sol2ndDerDt[9])*dsdx*dtdx
                           + solDr[4]*metricTerms2ndDer[0] + solDs[4]*metricTerms2ndDer[1] + solDt[4]*metricTerms2ndDer[2];

    sol2ndDerCart[0][1][0] = sol2ndDerDr[0]*drdx*drdy + sol2ndDerDs[5]*dsdx*dsdy + sol2ndDerDs[10]*dtdx*dtdy
                           + 0.5*(sol2ndDerDr[5]  + sol2ndDerDs[0])*(drdx*dsdy + dsdx*drdy)
                           + 0.5*(sol2ndDerDr[10] + sol2ndDerDt[0])*(drdx*dtdy + dtdx*drdy)
                           + 0.5*(sol2ndDerDs[10] + sol2ndDerDt[5])*(dsdx*dtdy + dtdx*dsdy)
                           + solDr[0]*metricTerms2ndDer[3] + solDs[0]*metricTerms2ndDer[4] + solDt[0]*metricTerms2ndDer[5];
    sol2ndDerCart[1][1][0] = sol2ndDerDr[1]*drdx*drdy + sol2ndDerDs[6]*dsdx*dsdy + sol2ndDerDs[11]*dtdx*dtdy
                           + 0.5*(sol2ndDerDr[6]  + sol2ndDerDs[1])*(drdx*dsdy + dsdx*drdy)
                           + 0.5*(sol2ndDerDr[11] + sol2ndDerDt[1])*(drdx*dtdy + dtdx*drdy)
                           + 0.5*(sol2ndDerDs[11] + sol2ndDerDt[6])*(dsdx*dtdy + dtdx*dsdy)
                           + solDr[1]*metricTerms2ndDer[3] + solDs[1]*metricTerms2ndDer[4] + solDt[1]*metricTerms2ndDer[5];
    sol2ndDerCart[2][1][0] = sol2ndDerDr[2]*drdx*drdy + sol2ndDerDs[7]*dsdx*dsdy + sol2ndDerDs[12]*dtdx*dtdy
                           + 0.5*(sol2ndDerDr[7]  + sol2ndDerDs[2])*(drdx*dsdy + dsdx*drdy)
                           + 0.5*(sol2ndDerDr[12] + sol2ndDerDt[2])*(drdx*dtdy + dtdx*drdy)
                           + 0.5*(sol2ndDerDs[12] + sol2ndDerDt[7])*(dsdx*dtdy + dtdx*dsdy)
                           + solDr[2]*metricTerms2ndDer[3] + solDs[2]*metricTerms2ndDer[4] + solDt[2]*metricTerms2ndDer[5];
    sol2ndDerCart[3][1][0] = sol2ndDerDr[3]*drdx*drdy + sol2ndDerDs[8]*dsdx*dsdy + sol2ndDerDs[13]*dtdx*dtdy
                           + 0.5*(sol2ndDerDr[8]  + sol2ndDerDs[3])*(drdx*dsdy + dsdx*drdy)
                           + 0.5*(sol2ndDerDr[13] + sol2ndDerDt[3])*(drdx*dtdy + dtdx*drdy)
                           + 0.5*(sol2ndDerDs[13] + sol2ndDerDt[8])*(dsdx*dtdy + dtdx*dsdy)
                           + solDr[3]*metricTerms2ndDer[3] + solDs[3]*metricTerms2ndDer[4] + solDt[3]*metricTerms2ndDer[5];
    sol2ndDerCart[4][1][0] = sol2ndDerDr[4]*drdx*drdy + sol2ndDerDs[9]*dsdx*dsdy + sol2ndDerDs[14]*dtdx*dtdy
                           + 0.5*(sol2ndDerDr[9]  + sol2ndDerDs[4])*(drdx*dsdy + dsdx*drdy)
                           + 0.5*(sol2ndDerDr[14] + sol2ndDerDt[4])*(drdx*dtdy + dtdx*drdy)
                           + 0.5*(sol2ndDerDs[14] + sol2ndDerDt[9])*(dsdx*dtdy + dtdx*dsdy)
                           + solDr[4]*metricTerms2ndDer[3] + solDs[4]*metricTerms2ndDer[4] + solDt[4]*metricTerms2ndDer[5];

    sol2ndDerCart[0][0][1] = sol2ndDerCart[0][1][0];
    sol2ndDerCart[1][0][1] = sol2ndDerCart[1][1][0];
    sol2ndDerCart[2][0][1] = sol2ndDerCart[2][1][0];
    sol2ndDerCart[3][0][1] = sol2ndDerCart[3][1][0];
    sol2ndDerCart[4][0][1] = sol2ndDerCart[4][1][0];

    sol2ndDerCart[0][1][1] = sol2ndDerDr[0]*drdy*drdy + sol2ndDerDs[5]*dsdy*dsdy + sol2ndDerDs[10]*dtdy*dtdy
                           + (sol2ndDerDr[5]  + sol2ndDerDs[0])*drdy*dsdy + (sol2ndDerDr[10] + sol2ndDerDt[0])*drdy*dtdy
                           + (sol2ndDerDs[10] + sol2ndDerDt[5])*dsdy*dtdy
                           + solDr[0]*metricTerms2ndDer[6] + solDs[0]*metricTerms2ndDer[7] + solDt[0]*metricTerms2ndDer[8];
    sol2ndDerCart[1][1][1] = sol2ndDerDr[1]*drdy*drdy + sol2ndDerDs[6]*dsdy*dsdy + sol2ndDerDs[11]*dtdy*dtdy
                           + (sol2ndDerDr[6]  + sol2ndDerDs[1])*drdy*dsdy + (sol2ndDerDr[11] + sol2ndDerDt[1])*drdy*dtdy
                           + (sol2ndDerDs[11] + sol2ndDerDt[6])*dsdy*dtdy
                           + solDr[1]*metricTerms2ndDer[6] + solDs[1]*metricTerms2ndDer[7] + solDt[1]*metricTerms2ndDer[8];
    sol2ndDerCart[2][1][1] = sol2ndDerDr[2]*drdy*drdy + sol2ndDerDs[7]*dsdy*dsdy + sol2ndDerDs[12]*dtdy*dtdy
                           + (sol2ndDerDr[7]  + sol2ndDerDs[2])*drdy*dsdy + (sol2ndDerDr[12] + sol2ndDerDt[2])*drdy*dtdy
                           + (sol2ndDerDs[12] + sol2ndDerDt[7])*dsdy*dtdy
                           + solDr[2]*metricTerms2ndDer[6] + solDs[2]*metricTerms2ndDer[7] + solDt[2]*metricTerms2ndDer[8];
    sol2ndDerCart[3][1][1] = sol2ndDerDr[3]*drdy*drdy + sol2ndDerDs[8]*dsdy*dsdy + sol2ndDerDs[13]*dtdy*dtdy
                           + (sol2ndDerDr[8]  + sol2ndDerDs[3])*drdy*dsdy + (sol2ndDerDr[13] + sol2ndDerDt[3])*drdy*dtdy
                           + (sol2ndDerDs[13] + sol2ndDerDt[8])*dsdy*dtdy
                           + solDr[3]*metricTerms2ndDer[6] + solDs[3]*metricTerms2ndDer[7] + solDt[3]*metricTerms2ndDer[8];
    sol2ndDerCart[4][1][1] = sol2ndDerDr[4]*drdy*drdy + sol2ndDerDs[9]*dsdy*dsdy + sol2ndDerDs[14]*dtdy*dtdy
                           + (sol2ndDerDr[9]  + sol2ndDerDs[4])*drdy*dsdy + (sol2ndDerDr[14] + sol2ndDerDt[4])*drdy*dtdy
                           + (sol2ndDerDs[14] + sol2ndDerDt[9])*dsdy*dtdy
                           + solDr[4]*metricTerms2ndDer[6] + solDs[4]*metricTerms2ndDer[7] + solDt[4]*metricTerms2ndDer[8];

    sol2ndDerCart[0][2][0] = sol2ndDerDr[0]*drdx*drdz + sol2ndDerDs[5]*dsdx*dsdz + sol2ndDerDs[10]*dtdx*dtdz
                           + 0.5*(sol2ndDerDr[5]  + sol2ndDerDs[0])*(drdx*dsdz + dsdx*drdz)
                           + 0.5*(sol2ndDerDr[10] + sol2ndDerDt[0])*(drdx*dtdz + dtdx*drdz)
                           + 0.5*(sol2ndDerDs[10] + sol2ndDerDt[5])*(dsdx*dtdz + dtdx*dsdz)
                           + solDr[0]*metricTerms2ndDer[9] + solDs[0]*metricTerms2ndDer[10] + solDt[0]*metricTerms2ndDer[11];
    sol2ndDerCart[1][2][0] = sol2ndDerDr[1]*drdx*drdz + sol2ndDerDs[6]*dsdx*dsdz + sol2ndDerDs[11]*dtdx*dtdz
                           + 0.5*(sol2ndDerDr[6]  + sol2ndDerDs[1])*(drdx*dsdz + dsdx*drdz)
                           + 0.5*(sol2ndDerDr[11] + sol2ndDerDt[1])*(drdx*dtdz + dtdx*drdz)
                           + 0.5*(sol2ndDerDs[11] + sol2ndDerDt[6])*(dsdx*dtdz + dtdx*dsdz)
                           + solDr[1]*metricTerms2ndDer[9] + solDs[1]*metricTerms2ndDer[10] + solDt[1]*metricTerms2ndDer[11];
    sol2ndDerCart[2][2][0] = sol2ndDerDr[2]*drdx*drdz + sol2ndDerDs[7]*dsdx*dsdz + sol2ndDerDs[12]*dtdx*dtdz
                           + 0.5*(sol2ndDerDr[7]  + sol2ndDerDs[2])*(drdx*dsdz + dsdx*drdz)
                           + 0.5*(sol2ndDerDr[12] + sol2ndDerDt[2])*(drdx*dtdz + dtdx*drdz)
                           + 0.5*(sol2ndDerDs[12] + sol2ndDerDt[7])*(dsdx*dtdz + dtdx*dsdz)
                           + solDr[2]*metricTerms2ndDer[9] + solDs[2]*metricTerms2ndDer[10] + solDt[2]*metricTerms2ndDer[11];
    sol2ndDerCart[3][2][0] = sol2ndDerDr[3]*drdx*drdz + sol2ndDerDs[8]*dsdx*dsdz + sol2ndDerDs[13]*dtdx*dtdz
                           + 0.5*(sol2ndDerDr[8]  + sol2ndDerDs[3])*(drdx*dsdz + dsdx*drdz)
                           + 0.5*(sol2ndDerDr[13] + sol2ndDerDt[3])*(drdx*dtdz + dtdx*drdz)
                           + 0.5*(sol2ndDerDs[13] + sol2ndDerDt[8])*(dsdx*dtdz + dtdx*dsdz)
                           + solDr[3]*metricTerms2ndDer[9] + solDs[3]*metricTerms2ndDer[10] + solDt[3]*metricTerms2ndDer[11];
    sol2ndDerCart[4][2][0] = sol2ndDerDr[4]*drdx*drdz + sol2ndDerDs[9]*dsdx*dsdz + sol2ndDerDs[14]*dtdx*dtdz
                           + 0.5*(sol2ndDerDr[9]  + sol2ndDerDs[4])*(drdx*dsdz + dsdx*drdz)
                           + 0.5*(sol2ndDerDr[14] + sol2ndDerDt[4])*(drdx*dtdz + dtdx*drdz)
                           + 0.5*(sol2ndDerDs[14] + sol2ndDerDt[9])*(dsdx*dtdz + dtdx*dsdz)
                           + solDr[4]*metricTerms2ndDer[9] + solDs[4]*metricTerms2ndDer[10] + solDt[4]*metricTerms2ndDer[11];

    sol2ndDerCart[0][0][2] = sol2ndDerCart[0][2][0];
    sol2ndDerCart[1][0][2] = sol2ndDerCart[1][2][0];
    sol2ndDerCart[2][0][2] = sol2ndDerCart[2][2][0];
    sol2ndDerCart[3][0][2] = sol2ndDerCart[3][2][0];
    sol2ndDerCart[4][0][2] = sol2ndDerCart[4][2][0];

    sol2ndDerCart[0][2][1] = sol2ndDerDr[0]*drdy*drdz + sol2ndDerDs[5]*dsdy*dsdz + sol2ndDerDs[10]*dtdy*dtdz
                           + 0.5*(sol2ndDerDr[5]  + sol2ndDerDs[0])*(drdy*dsdz + dsdy*drdz)
                           + 0.5*(sol2ndDerDr[10] + sol2ndDerDt[0])*(drdy*dtdz + dtdy*drdz)
                           + 0.5*(sol2ndDerDs[10] + sol2ndDerDt[5])*(dsdy*dtdz + dtdy*dsdz)
                           + solDr[0]*metricTerms2ndDer[12] + solDs[0]*metricTerms2ndDer[13] + solDt[0]*metricTerms2ndDer[14];
    sol2ndDerCart[1][2][1] = sol2ndDerDr[1]*drdy*drdz + sol2ndDerDs[6]*dsdy*dsdz + sol2ndDerDs[11]*dtdy*dtdz
                           + 0.5*(sol2ndDerDr[6]  + sol2ndDerDs[1])*(drdy*dsdz + dsdy*drdz)
                           + 0.5*(sol2ndDerDr[11] + sol2ndDerDt[1])*(drdy*dtdz + dtdy*drdz)
                           + 0.5*(sol2ndDerDs[11] + sol2ndDerDt[6])*(dsdy*dtdz + dtdy*dsdz)
                           + solDr[1]*metricTerms2ndDer[12] + solDs[1]*metricTerms2ndDer[13] + solDt[1]*metricTerms2ndDer[14];
    sol2ndDerCart[2][2][1] = sol2ndDerDr[2]*drdy*drdz + sol2ndDerDs[7]*dsdy*dsdz + sol2ndDerDs[12]*dtdy*dtdz
                           + 0.5*(sol2ndDerDr[7]  + sol2ndDerDs[2])*(drdy*dsdz + dsdy*drdz)
                           + 0.5*(sol2ndDerDr[12] + sol2ndDerDt[2])*(drdy*dtdz + dtdy*drdz)
                           + 0.5*(sol2ndDerDs[12] + sol2ndDerDt[7])*(dsdy*dtdz + dtdy*dsdz)
                           + solDr[2]*metricTerms2ndDer[12] + solDs[2]*metricTerms2ndDer[13] + solDt[2]*metricTerms2ndDer[14];
    sol2ndDerCart[3][2][1] = sol2ndDerDr[3]*drdy*drdz + sol2ndDerDs[8]*dsdy*dsdz + sol2ndDerDs[13]*dtdy*dtdz
                           + 0.5*(sol2ndDerDr[8]  + sol2ndDerDs[3])*(drdy*dsdz + dsdy*drdz)
                           + 0.5*(sol2ndDerDr[13] + sol2ndDerDt[3])*(drdy*dtdz + dtdy*drdz)
                           + 0.5*(sol2ndDerDs[13] + sol2ndDerDt[8])*(dsdy*dtdz + dtdy*dsdz)
                           + solDr[3]*metricTerms2ndDer[12] + solDs[3]*metricTerms2ndDer[13] + solDt[3]*metricTerms2ndDer[14];
    sol2ndDerCart[4][2][1] = sol2ndDerDr[4]*drdy*drdz + sol2ndDerDs[9]*dsdy*dsdz + sol2ndDerDs[14]*dtdy*dtdz
                           + 0.5*(sol2ndDerDr[9]  + sol2ndDerDs[4])*(drdy*dsdz + dsdy*drdz)
                           + 0.5*(sol2ndDerDr[14] + sol2ndDerDt[4])*(drdy*dtdz + dtdy*drdz)
                           + 0.5*(sol2ndDerDs[14] + sol2ndDerDt[9])*(dsdy*dtdz + dtdy*dsdz)
                           + solDr[4]*metricTerms2ndDer[12] + solDs[4]*metricTerms2ndDer[13] + solDt[4]*metricTerms2ndDer[14];

    sol2ndDerCart[0][1][2] = sol2ndDerCart[0][2][1];
    sol2ndDerCart[1][1][2] = sol2ndDerCart[1][2][1];
    sol2ndDerCart[2][1][2] = sol2ndDerCart[2][2][1];
    sol2ndDerCart[3][1][2] = sol2ndDerCart[3][2][1];
    sol2ndDerCart[4][1][2] = sol2ndDerCart[4][2][1];

    sol2ndDerCart[0][2][2] = sol2ndDerDr[0]*drdz*drdz + sol2ndDerDs[5]*dsdz*dsdz + sol2ndDerDs[10]*dtdz*dtdz
                           + (sol2ndDerDr[5]  + sol2ndDerDs[0])*drdz*dsdz + (sol2ndDerDr[10] + sol2ndDerDt[0])*drdz*dtdz
                           + (sol2ndDerDs[10] + sol2ndDerDt[5])*dsdz*dtdz
                           + solDr[0]*metricTerms2ndDer[15] + solDs[0]*metricTerms2ndDer[16] + solDt[0]*metricTerms2ndDer[17];
    sol2ndDerCart[1][2][2] = sol2ndDerDr[1]*drdz*drdz + sol2ndDerDs[6]*dsdz*dsdz + sol2ndDerDs[11]*dtdz*dtdz
                           + (sol2ndDerDr[6]  + sol2ndDerDs[1])*drdz*dsdz + (sol2ndDerDr[11] + sol2ndDerDt[1])*drdz*dtdz
                           + (sol2ndDerDs[11] + sol2ndDerDt[6])*dsdz*dtdz
                           + solDr[1]*metricTerms2ndDer[15] + solDs[1]*metricTerms2ndDer[16] + solDt[1]*metricTerms2ndDer[17];
    sol2ndDerCart[2][2][2] = sol2ndDerDr[2]*drdz*drdz + sol2ndDerDs[7]*dsdz*dsdz + sol2ndDerDs[12]*dtdz*dtdz
                           + (sol2ndDerDr[7]  + sol2ndDerDs[2])*drdz*dsdz + (sol2ndDerDr[12] + sol2ndDerDt[2])*drdz*dtdz
                           + (sol2ndDerDs[12] + sol2ndDerDt[7])*dsdz*dtdz
                           + solDr[2]*metricTerms2ndDer[15] + solDs[2]*metricTerms2ndDer[16] + solDt[2]*metricTerms2ndDer[17];
    sol2ndDerCart[3][2][2] = sol2ndDerDr[3]*drdz*drdz + sol2ndDerDs[8]*dsdz*dsdz + sol2ndDerDs[13]*dtdz*dtdz
                           + (sol2ndDerDr[8]  + sol2ndDerDs[3])*drdz*dsdz + (sol2ndDerDr[13] + sol2ndDerDt[3])*drdz*dtdz
                           + (sol2ndDerDs[13] + sol2ndDerDt[8])*dsdz*dtdz
                           + solDr[3]*metricTerms2ndDer[15] + solDs[3]*metricTerms2ndDer[16] + solDt[3]*metricTerms2ndDer[17];
    sol2ndDerCart[4][2][2] = sol2ndDerDr[4]*drdz*drdz + sol2ndDerDs[9]*dsdz*dsdz + sol2ndDerDs[14]*dtdz*dtdz
                           + (sol2ndDerDr[9]  + sol2ndDerDs[4])*drdz*dsdz + (sol2ndDerDr[14] + sol2ndDerDt[4])*drdz*dtdz
                           + (sol2ndDerDs[14] + sol2ndDerDt[9])*dsdz*dtdz
                           + solDr[4]*metricTerms2ndDer[15] + solDs[4]*metricTerms2ndDer[16] + solDt[4]*metricTerms2ndDer[17];

    /*--- Compute the Cartesian gradients of the pressure, velocity components,
          static energy and dynamic viscosity. ---*/
    const su2double dpdx = Gamma_Minus_One*(solGradCart[4][0] + kinEnergy*solGradCart[0][0]
                         -                  u*solGradCart[1][0] - v*solGradCart[2][0] - w*solGradCart[3][0]);
    const su2double dpdy = Gamma_Minus_One*(solGradCart[3][1] + kinEnergy*solGradCart[0][1]
                         -                  u*solGradCart[1][1] - v*solGradCart[2][1] - w*solGradCart[3][1]);
    const su2double dpdz = Gamma_Minus_One*(solGradCart[3][2] + kinEnergy*solGradCart[0][2]
                         -                  u*solGradCart[1][2] - v*solGradCart[2][2] - w*solGradCart[3][2]);

    const su2double dudx = DensityInv*(solGradCart[1][0] - u*solGradCart[0][0]);
    const su2double dudy = DensityInv*(solGradCart[1][1] - u*solGradCart[0][1]);
    const su2double dudz = DensityInv*(solGradCart[1][2] - u*solGradCart[0][2]);

    const su2double dvdx = DensityInv*(solGradCart[2][0] - v*solGradCart[0][0]);
    const su2double dvdy = DensityInv*(solGradCart[2][1] - v*solGradCart[0][1]);
    const su2double dvdz = DensityInv*(solGradCart[2][2] - v*solGradCart[0][2]);

    const su2double dwdx = DensityInv*(solGradCart[3][0] - w*solGradCart[0][0]);
    const su2double dwdy = DensityInv*(solGradCart[3][1] - w*solGradCart[0][1]);
    const su2double dwdz = DensityInv*(solGradCart[3][2] - w*solGradCart[0][2]);

    const su2double dedx = DensityInv*(solGradCart[4][0] - TotalEnergy*solGradCart[0][0])
                         - u*dudx - v*dvdx - w*dwdx;
    const su2double dedy = DensityInv*(solGradCart[4][1] - TotalEnergy*solGradCart[0][1])
                         - u*dudy - v*dvdy - w*dwdy;
    const su2double dedz = DensityInv*(solGradCart[4][2] - TotalEnergy*solGradCart[0][2])
                         - u*dudz - v*dvdz - w*dwdz;

    const su2double dViscLamdx = CvInv*dedx*dViscLamdT;
    const su2double dViscLamdy = CvInv*dedy*dViscLamdT;
    const su2double dViscLamdz = CvInv*dedz*dViscLamdT;

    /*--- Compute the second derivatives of the velocity components. ---*/
    const su2double d2udxx = DensityInv*(sol2ndDerCart[1][0][0] - u*sol2ndDerCart[0][0][0]
                           +    2.0*DensityInv*solGradCart[0][0]*(u*solGradCart[0][0] - solGradCart[1][0]));
    const su2double d2udyy = DensityInv*(sol2ndDerCart[1][1][1] - u*sol2ndDerCart[0][1][1]
                           +    2.0*DensityInv*solGradCart[0][1]*(u*solGradCart[0][1] - solGradCart[1][1]));
    const su2double d2udzz = DensityInv*(sol2ndDerCart[1][2][2] - u*sol2ndDerCart[0][2][2]
                           +    2.0*DensityInv*solGradCart[0][2]*(u*solGradCart[0][2] - solGradCart[1][2]));
    const su2double d2udxy = DensityInv*(sol2ndDerCart[1][0][1] - u*sol2ndDerCart[0][0][1]
                           +       DensityInv*(solGradCart[0][0]*(u*solGradCart[0][1] - solGradCart[1][1])
                           +                   solGradCart[0][1]*(u*solGradCart[0][0] - solGradCart[1][0])));
    const su2double d2udxz = DensityInv*(sol2ndDerCart[1][0][2] - u*sol2ndDerCart[0][0][2]
                           +       DensityInv*(solGradCart[0][0]*(u*solGradCart[0][2] - solGradCart[1][2])
                           +                   solGradCart[0][2]*(u*solGradCart[0][0] - solGradCart[1][0])));
    const su2double d2udyz = DensityInv*(sol2ndDerCart[1][1][2] - u*sol2ndDerCart[0][1][2]
                           +       DensityInv*(solGradCart[0][1]*(u*solGradCart[0][2] - solGradCart[1][2])
                           +                   solGradCart[0][2]*(u*solGradCart[0][1] - solGradCart[1][1])));

    const su2double d2vdxx = DensityInv*(sol2ndDerCart[2][0][0] - v*sol2ndDerCart[0][0][0]
                           +    2.0*DensityInv*solGradCart[0][0]*(v*solGradCart[0][0] - solGradCart[2][0]));
    const su2double d2vdyy = DensityInv*(sol2ndDerCart[2][1][1] - v*sol2ndDerCart[0][1][1]
                           +    2.0*DensityInv*solGradCart[0][1]*(v*solGradCart[0][1] - solGradCart[2][1]));
    const su2double d2vdzz = DensityInv*(sol2ndDerCart[2][2][2] - v*sol2ndDerCart[0][2][2]
                           +    2.0*DensityInv*solGradCart[0][2]*(v*solGradCart[0][2] - solGradCart[2][2]));
    const su2double d2vdxy = DensityInv*(sol2ndDerCart[2][0][1] - v*sol2ndDerCart[0][0][1]
                           +       DensityInv*(solGradCart[0][0]*(v*solGradCart[0][1] - solGradCart[2][1])
                           +                   solGradCart[0][1]*(v*solGradCart[0][0] - solGradCart[2][0])));
    const su2double d2vdxz = DensityInv*(sol2ndDerCart[2][0][2] - v*sol2ndDerCart[0][0][2]
                           +       DensityInv*(solGradCart[0][0]*(v*solGradCart[0][2] - solGradCart[2][2])
                           +                   solGradCart[0][2]*(v*solGradCart[0][0] - solGradCart[2][0])));
    const su2double d2vdyz = DensityInv*(sol2ndDerCart[2][1][2] - v*sol2ndDerCart[0][1][2]
                           +       DensityInv*(solGradCart[0][1]*(v*solGradCart[0][2] - solGradCart[2][2])
                           +                   solGradCart[0][2]*(v*solGradCart[0][1] - solGradCart[2][1])));

    const su2double d2wdxx = DensityInv*(sol2ndDerCart[3][0][0] - w*sol2ndDerCart[0][0][0]
                           +    2.0*DensityInv*solGradCart[0][0]*(w*solGradCart[0][0] - solGradCart[3][0]));
    const su2double d2wdyy = DensityInv*(sol2ndDerCart[3][1][1] - w*sol2ndDerCart[0][1][1]
                           +    2.0*DensityInv*solGradCart[0][1]*(w*solGradCart[0][1] - solGradCart[3][1]));
    const su2double d2wdzz = DensityInv*(sol2ndDerCart[3][2][2] - w*sol2ndDerCart[0][2][2]
                           +    2.0*DensityInv*solGradCart[0][2]*(w*solGradCart[0][2] - solGradCart[3][2]));
    const su2double d2wdxy = DensityInv*(sol2ndDerCart[3][0][1] - w*sol2ndDerCart[0][0][1]
                           +       DensityInv*(solGradCart[0][0]*(w*solGradCart[0][1] - solGradCart[3][1])
                           +                   solGradCart[0][1]*(w*solGradCart[0][0] - solGradCart[3][0])));
    const su2double d2wdxz = DensityInv*(sol2ndDerCart[3][0][2] - w*sol2ndDerCart[0][0][2]
                           +       DensityInv*(solGradCart[0][0]*(w*solGradCart[0][2] - solGradCart[3][2])
                           +                   solGradCart[0][2]*(w*solGradCart[0][0] - solGradCart[3][0])));
    const su2double d2wdyz = DensityInv*(sol2ndDerCart[3][1][2] - w*sol2ndDerCart[0][1][2]
                           +       DensityInv*(solGradCart[0][1]*(w*solGradCart[0][2] - solGradCart[3][2])
                           +                   solGradCart[0][2]*(w*solGradCart[0][1] - solGradCart[3][1])));

    /*--- Compute the second derivatives of the static energy. Note that this
          term appears in the heat flux and therefore only the pure second
          derivatives are needed. Hence, the cross-derivatives are omitted. ---*/
    const su2double d2edxx = DensityInv*(sol2ndDerCart[4][0][0] - TotalEnergy*sol2ndDerCart[0][0][0]
                           +    2.0*DensityInv*solGradCart[0][0]*(TotalEnergy*solGradCart[0][0] - solGradCart[4][0]))
                           -    u*d2udxx - dudx*dudx - v*d2vdxx - dvdx*dvdx - w*d2wdxx - dwdx*dwdx;
    const su2double d2edyy = DensityInv*(sol2ndDerCart[4][1][1] - TotalEnergy*sol2ndDerCart[0][1][1]
                           +    2.0*DensityInv*solGradCart[0][1]*(TotalEnergy*solGradCart[0][1] - solGradCart[4][1]))
                           -    u*d2udyy - dudy*dudy - v*d2vdyy - dvdy*dvdy - w*d2wdyy - dwdy*dwdy;
    const su2double d2edzz = DensityInv*(sol2ndDerCart[4][2][2] - TotalEnergy*sol2ndDerCart[0][2][2]
                           +    2.0*DensityInv*solGradCart[0][2]*(TotalEnergy*solGradCart[0][2] - solGradCart[4][2]))
                           -    u*d2udzz - dudz*dudz - v*d2vdzz - dvdz*dvdz - w*d2wdzz - dwdz*dwdz;

    /*--- If an SGS model is used the eddy viscosity and its spatial
          derivatives must be computed. ---*/
    su2double ViscosityTurb = 0.0;
    su2double dViscTurbdx = 0.0, dViscTurbdy = 0.0, dViscTurbdz = 0.0;

    if( SGSModelUsed ) {
      const su2double dist = elem->wallDistance[i];
      ViscosityTurb = SGSModel->ComputeEddyViscosity_3D(sol[0], dudx, dudy, dudz,
                                                        dvdx, dvdy, dvdz, dwdx,
                                                        dwdy, dwdz, lenScale, dist);

      SGSModel->ComputeGradEddyViscosity_3D(sol[0], solGradCart[0][0], solGradCart[0][1],
                                            solGradCart[0][2], dudx, dudy, dudz, dvdx,
                                            dvdy, dvdz, dwdx, dwdy, dwdz, d2udxx, d2udyy,
                                            d2udzz, d2udxy, d2udxz, d2udyz, d2vdxx, d2vdyy,
                                            d2vdzz, d2vdxy, d2vdxz, d2vdyz, d2wdxx, d2wdyy,
                                            d2wdzz, d2wdxy, d2wdxz, d2wdyz, lenScale,
                                            dist, dViscTurbdx, dViscTurbdy, dViscTurbdz);
    }

    /*--- Compute the total viscosity, the total heat conductivity and their
          gradients. Note that the heat conductivity is divided by the Cv,
          because gradients of internal energy are computed and not temperature. ---*/
    const su2double Viscosity = ViscosityLam + ViscosityTurb;
    const su2double kOverCv = ViscosityLam *factHeatFlux_Lam
                            + ViscosityTurb*factHeatFlux_Turb;

    const su2double dViscDx = dViscLamdx + dViscTurbdx;
    const su2double dViscDy = dViscLamdy + dViscTurbdy;
    const su2double dViscDz = dViscLamdz + dViscTurbdz;

    const su2double dkOverCvdx = dViscLamdx *factHeatFlux_Lam
                               + dViscTurbdx*factHeatFlux_Turb;
    const su2double dkOverCvdy = dViscLamdy *factHeatFlux_Lam
                               + dViscTurbdy*factHeatFlux_Turb;
    const su2double dkOverCvdz = dViscLamdz *factHeatFlux_Lam
                               + dViscTurbdz*factHeatFlux_Turb;

    /* Abbreviations, which make it easier to compute the divergence term. */
    const su2double abv1 = solGradCart[1][0] + solGradCart[2][1] + solGradCart[3][2];
    const su2double abv2 = u*solGradCart[0][0] + v*solGradCart[0][1] + w*solGradCart[0][2];
    const su2double abv3 = u*(solGradCart[4][0] + dpdx) + v*(solGradCart[4][1] + dpdy)
                         + w*(solGradCart[4][2] + dpdz);
    const su2double abv4 = dudx + dvdy + dwdz;

    /*--- Compute the divergence of the grid velocity.
          SET TO ZERO FOR NOW. THIS IS NOT CORRECT!!!!. ---*/
    const su2double divGridVel = 0.0;

    /*--- Set the pointer to store the divergence terms for this integration
          point and compute these terms, multiplied by the integration weight
          and Jacobian. ---*/
    const su2double weightJac = weights[i]*Jac;
    su2double *divFluxInt     = divFlux + 5*i;   /* nVar*i. */

    divFluxInt[0] = weightJac*(abv1
                  -            sol[0]*divGridVel            - gridVel[0]*solGradCart[0][0]
                  -            gridVel[1]*solGradCart[0][1] - gridVel[2]*solGradCart[0][2]);
    divFluxInt[1] = weightJac*(dpdx + u*(abv1-abv2) - lambdaOverMu*abv4*dViscDx
                  +            u*solGradCart[1][0] + v*solGradCart[1][1] + w*solGradCart[1][2]
                  -            lambdaOverMu*Viscosity*(d2udxx + d2vdxy + d2wdxz)
                  -            Viscosity*(2.0*d2udxx + d2udyy + d2vdxy + d2udzz + d2wdxz)
                  -            2.0*dViscDx*dudx - dViscDy*(dudy+dvdx) - dViscDz*(dudz+dwdx)
                  -            sol[1]*divGridVel            - gridVel[0]*solGradCart[1][0]
                  -            gridVel[1]*solGradCart[1][1] - gridVel[2]*solGradCart[1][2]);
    divFluxInt[2] = weightJac*(dpdy + v*(abv1-abv2) - lambdaOverMu*abv4*dViscDy
                  +            u*solGradCart[2][0] + v*solGradCart[2][1] + w*solGradCart[2][2]
                  -            lambdaOverMu*Viscosity*(d2udxy + d2vdxy + d2wdyz)
                  -            Viscosity*(d2udxy + d2vdxx + 2.0*d2vdyy + d2vdzz + d2wdyz)
                  -            dViscDx*(dudy+dvdx) - 2.0*dViscDy*dvdy - dViscDz*(dvdz+dwdy)
                  -            sol[2]*divGridVel            - gridVel[0]*solGradCart[2][0]
                  -            gridVel[1]*solGradCart[2][1] - gridVel[2]*solGradCart[2][2]);
    divFluxInt[3] = weightJac*(dpdz + w*(abv1-abv2) - lambdaOverMu*abv4*dViscDz
                  +            u*solGradCart[3][0] + v*solGradCart[3][1] + w*solGradCart[3][2]
                  -            lambdaOverMu*Viscosity*(d2udxz + d2vdyz + d2wdzz)
                  -            Viscosity*(d2udxz + d2wdxx + d2vdyz + d2wdyy + 2.0*d2wdzz)
                  -            dViscDx*(dudz+dwdx) - dViscDy*(dvdz+dwdy) - 2.0*dViscDz*dwdz
                  -            sol[3]*divGridVel            - gridVel[0]*solGradCart[3][0]
                  -            gridVel[1]*solGradCart[3][1] - gridVel[2]*solGradCart[3][2]);
    divFluxInt[4] = weightJac*(abv3 + Htot*(abv1 - abv2)
                  -            abv4*lambdaOverMu*(Viscosity*abv4 + u*dViscDx + v*dViscDy + w*dViscDz)
                  -            dkOverCvdx*dedx - dkOverCvdy*dedy - dkOverCvdz*dedz
                  -            kOverCv*(d2edxx + d2edyy + d2edzz)
                  -            (Viscosity*dudx + u*dViscDx)*2.0*dudx
                  -            (Viscosity*dvdy + v*dViscDy)*2.0*dvdy
                  -            (Viscosity*dwdz + w*dViscDz)*2.0*dwdz
                  -            (Viscosity*dudy + u*dViscDy + Viscosity*dvdx + v*dViscDx)*(dudy + dvdx)
                  -            (Viscosity*dudz + u*dViscDz + Viscosity*dwdx + w*dViscDx)*(dudz + dwdx)
                  -            (Viscosity*dvdz + v*dViscDz + Viscosity*dwdy + w*dViscDy)*(dvdz + dwdy)
                  -            Viscosity*u*(d2udxx+d2udyy+d2udzz + (1.0+lambdaOverMu)*(d2udxx+d2vdxy+d2wdxz))
                  -            Viscosity*v*(d2vdxx+d2vdyy+d2vdzz + (1.0+lambdaOverMu)*(d2udxy+d2vdyy+d2wdyz))
                  -            Viscosity*w*(d2wdxx+d2wdyy+d2wdzz + (1.0+lambdaOverMu)*(d2udxz+d2vdyz+d2wdzz))
                  -            sol[4]*divGridVel            - gridVel[0]*solGradCart[4][0]
                  -            gridVel[1]*solGradCart[4][1] - gridVel[2]*solGradCart[4][2]);
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Compute the residual in the DOFs, which is the matrix      ---*/
  /*---         product of basisFunctionsIntTrans and divFlux.             ---*/
  /*--------------------------------------------------------------------------*/

  DenseMatrixProduct(nDOFs, nVar, nInt, basisFunctionsIntTrans, divFlux, res);
}

void CFEM_DG_NSSolver::Shock_Capturing_DG(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd) {

  /*--- Run shock capturing algorithm ---*/
  switch( config->GetKind_FEM_DG_Shock() ) {
    case NONE:
      break;
    case PERSSON:
      Shock_Capturing_DG_Persson(elemBeg, elemEnd);
      break;
  }

}
void CFEM_DG_NSSolver::Shock_Capturing_DG_Persson(const unsigned long elemBeg,
                                                  const unsigned long elemEnd) {

  /*--- Dummy variable for storing shock sensor value temporarily ---*/
  su2double sensorVal, sensorLowerBound, machNorm, machMax;
  su2double DensityInv, Velocity2, StaticEnergy, SoundSpeed2, Velocity2Rel;

  bool shockExist;
  unsigned short nDOFsPm1;       // Number of DOFs up to polynomial degree p-1

  /*--- Loop over the given range of elements to sense the shock. If shock exists,
        add artificial viscosity for DG FEM formulation to the residual.  ---*/
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

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
    const unsigned short timeLevel = volElem[l].timeLevel;
    const su2double *solDOFs = VecWorkSolDOFs[timeLevel].data()
                             + nVar*volElem[l].offsetDOFsSolThisTimeLevel;

    /* Temporary storage of mach number for DOFs in this element. */
    su2double *machSolDOFs = VecTmpMemory.data();
    su2double *vecTemp     = machSolDOFs + nDOFs;

    /* Calculate primitive variables and mach number for DOFs in this element.
       Also, track the maximum mach number in this element. */
    for(unsigned short iInd=0; iInd<nDOFs; ++iInd) {

      const su2double *sol     = solDOFs + iInd*nVar;
      const su2double *gridVel = volElem[l].gridVelocitiesSolDOFs.data() + iInd*nDim;
      DensityInv = 1.0/sol[0];
      Velocity2 = 0.0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double vel    = sol[iDim]*DensityInv;
        const su2double velRel = vel - gridVel[iDim-1];
        Velocity2    += vel*vel;
        Velocity2Rel += velRel*velRel;
      }

      StaticEnergy = sol[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
      SoundSpeed2 = FluidModel->GetSoundSpeed2();
      machSolDOFs[iInd] = sqrt( Velocity2Rel/SoundSpeed2 );
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
      // Need to develop thorough method for general problems.
      switch ( nPoly ) {
        case 1:  sensorLowerBound =  -6.0; break;
        case 2:  sensorLowerBound = -12.0; break;
        case 3:  sensorLowerBound = -12.0; break;
        case 4:  sensorLowerBound = -17.0; break;
        default: sensorLowerBound = -17.0; break;
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

void CFEM_DG_NSSolver::Volume_Residual(CConfig             *config,
                                       const unsigned long elemBeg,
                                       const unsigned long elemEnd) {

  /*--- Determine whether body force term is present in configuration. ---*/
  bool body_force = config->GetBody_Force();
  const su2double *body_force_vector = body_force ? config->GetBody_Force_Vector() : NULL;

  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

  /*--- Set the pointers for the local arrays. ---*/
  double tick = 0.0;

  su2double *sources       = VecTmpMemory.data();
  su2double *solAndGradInt = sources       + nIntegrationMax*nVar;
  su2double *fluxes        = solAndGradInt + nIntegrationMax*nVar*(nDim+1);

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /*--- Loop over the given element range to compute the contribution of the
        volume integral in the DG FEM formulation to the residual. ---*/
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* Get the data from the corresponding standard element. */
    const unsigned short ind             = volElem[l].indStandardElement;
    const unsigned short nInt            = standardElementsSol[ind].GetNIntegration();
    const unsigned short nDOFs           = volElem[l].nDOFsSol;
    const su2double *matBasisInt         = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
    const su2double *matDerBasisIntTrans = standardElementsSol[ind].GetDerMatBasisFunctionsIntTrans();
    const su2double *matBasisIntTrans    = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
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
    const unsigned short timeLevel = volElem[l].timeLevel;
    const su2double *solDOFs = VecWorkSolDOFs[timeLevel].data()
                             + nVar*volElem[l].offsetDOFsSolThisTimeLevel;

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
    const unsigned short offDeriv = nVar*nInt;

    /* Make a distinction between two and three space dimensions
        in order to have the most efficient code. */
    switch( nDim ) {

      case 2: {

        /* 2D simulation. Loop over the integration points to compute
           the fluxes. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the metric terms and grid velocities in this
             integration point and compute the inverse of the Jacobian. */
          const su2double *metricTerms = volElem[l].metricTerms.data()
                                       + i*nMetricPerPoint;
          const su2double *gridVel     = volElem[l].gridVelocities.data() + 2*i;  /* nDim*i. */
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
          const su2double *sol    = solAndGradInt + nVar*i;
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
          if( SGSModelUsed )
            ViscosityTurb = SGSModel->ComputeEddyViscosity_2D(sol[0], dudx, dudy,
                                                              dvdx, dvdy, lenScale,
                                                              volElem[l].wallDistance[i]);

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

          /* Compute the relative velocities w.r.t. the grid. */
          const su2double uRel = u - gridVel[0];
          const su2double vRel = v - gridVel[1];

          /* Compute the viscous normal stress minus the pressure. */
          const su2double tauxxMP = tauxx - Pressure;
          const su2double tauyyMP = tauyy - Pressure;

          /* Set the pointer for the fluxes in this integration point. */
          su2double *flux = fluxes + i*nDim*nVar;

          /*--- Fluxes in r-direction. */
          const su2double Ur = uRel*wDrdx + vRel*wDrdy;

          flux[0] = sol[0]*Ur;
          flux[1] = sol[1]*Ur - tauxxMP*wDrdx - tauxy*wDrdy;
          flux[2] = sol[2]*Ur - tauxy*wDrdx - tauyyMP*wDrdy;
          flux[3] = sol[3]*Ur - (u*tauxxMP + v*tauxy + qx)*wDrdx
                              - (u*tauxy + v*tauyyMP + qy)*wDrdy;

          /*--- Fluxes in s-direction. */
          const su2double Us = uRel*wDsdx + vRel*wDsdy;

          flux[4] = sol[0]*Us;
          flux[5] = sol[1]*Us - tauxxMP*wDsdx - tauxy*wDsdy;
          flux[6] = sol[2]*Us - tauxy*wDsdx - tauyyMP*wDsdy;
          flux[7] = sol[3]*Us - (u*tauxxMP + v*tauxy + qx)*wDsdx
                              - (u*tauxy + v*tauyyMP + qy)*wDsdy;

          /*--- If needed, compute the body forces in this integration point.
                Note that the source terms are multiplied with minus the
                integration weight in order to be consistent with the
                formulation of the residual. ---*/
          if( body_force ) {
            su2double *source         = sources + i*nVar;
            const su2double weightJac = weights[i]*Jac;

            source[0] =  0.0;
            source[1] = -weightJac*body_force_vector[0];
            source[2] = -weightJac*body_force_vector[1];
            source[3] = -weightJac*(u*body_force_vector[0] + v*body_force_vector[1]);
          }
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /* 3D simulation. Loop over the integration points to compute
           the fluxes. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the metric terms and grid velocities in this
             integration point and compute the inverse of the Jacobian. */
          const su2double *metricTerms = volElem[l].metricTerms.data()
                                       + i*nMetricPerPoint;
          const su2double *gridVel     = volElem[l].gridVelocities.data() + 3*i;  /* nDim*i. */
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
          const su2double *sol    = solAndGradInt + nVar*i;
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
          if( SGSModelUsed )
            ViscosityTurb = SGSModel->ComputeEddyViscosity_3D(sol[0], dudx, dudy, dudz,
                                                              dvdx, dvdy, dvdz, dwdx,
                                                              dwdy, dwdz, lenScale,
                                                              volElem[l].wallDistance[i]);

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

          /* Compute the relative velocities w.r.t. the grid. */
          const su2double uRel = u - gridVel[0];
          const su2double vRel = v - gridVel[1];
          const su2double wRel = w - gridVel[2];

          /* Compute the viscous normal stress minus the pressure. */
          const su2double tauxxMP = tauxx - Pressure;
          const su2double tauyyMP = tauyy - Pressure;
          const su2double tauzzMP = tauzz - Pressure;

          /* Set the pointer for the fluxes in this integration point. */
          su2double *flux = fluxes + i*nDim*nVar;

          /*--- Fluxes in r-direction. */
          const su2double Ur = uRel*wDrdx + vRel*wDrdy + wRel*wDrdz;

          flux[0] = sol[0]*Ur;
          flux[1] = sol[1]*Ur - tauxxMP*wDrdx - tauxy*wDrdy - tauxz*wDrdz;
          flux[2] = sol[2]*Ur - tauxy*wDrdx - tauyyMP*wDrdy - tauyz*wDrdz;
          flux[3] = sol[3]*Ur - tauxz*wDrdx - tauyz*wDrdy - tauzzMP*wDrdz;
          flux[4] = sol[4]*Ur - (u*tauxxMP + v*tauxy + w*tauxz + qx)*wDrdx
                              - (u*tauxy + v*tauyyMP + w*tauyz + qy)*wDrdy
                              - (u*tauxz + v*tauyz + w*tauzzMP + qz)*wDrdz;

          /*--- Fluxes in s-direction. */
          const su2double Us = uRel*wDsdx + vRel*wDsdy + wRel*wDsdz;

          flux[5] = sol[0]*Us;
          flux[6] = sol[1]*Us - tauxxMP*wDsdx - tauxy*wDsdy - tauxz*wDsdz;
          flux[7] = sol[2]*Us - tauxy*wDsdx - tauyyMP*wDsdy - tauyz*wDsdz;
          flux[8] = sol[3]*Us - tauxz*wDsdx - tauyz*wDsdy - tauzzMP*wDsdz;
          flux[9] = sol[4]*Us - (u*tauxxMP + v*tauxy + w*tauxz + qx)*wDsdx
                              - (u*tauxy + v*tauyyMP + w*tauyz + qy)*wDsdy
                              - (u*tauxz + v*tauyz + w*tauzzMP + qz)*wDsdz;

          /*--- Fluxes in t-direction. */
          const su2double Ut = uRel*wDtdx + vRel*wDtdy + wRel*wDtdz;

          flux[10] = sol[0]*Ut;
          flux[11] = sol[1]*Ut - tauxxMP*wDtdx - tauxy*wDtdy - tauxz*wDtdz;
          flux[12] = sol[2]*Ut - tauxy*wDtdx - tauyyMP*wDtdy - tauyz*wDtdz;
          flux[13] = sol[3]*Ut - tauxz*wDtdx - tauyz*wDtdy - tauzzMP*wDtdz;
          flux[14] = sol[4]*Ut - (u*tauxxMP + v*tauxy + w*tauxz + qx)*wDtdx
                               - (u*tauxy + v*tauyyMP + w*tauyz + qy)*wDtdy
                               - (u*tauxz + v*tauyz + w*tauzzMP + qz)*wDtdz;

          /*--- If needed, compute the body forces in this integration point.
                Note that the source terms are multiplied with minus the
                integration weight in order to be consistent with the
                formulation of the residual. ---*/
          if( body_force ) {
            su2double *source         = sources + i*nVar;
            const su2double weightJac = weights[i]*Jac;

            source[0] =  0.0;
            source[1] = -weightJac*body_force_vector[0];
            source[2] = -weightJac*body_force_vector[1];
            source[3] = -weightJac*body_force_vector[2];
            source[4] = -weightJac*(u*body_force_vector[0] + v*body_force_vector[1]
                      +             w*body_force_vector[2]);
          }
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

    /* Add the contribution from the source terms, if needed. Use solAndGradInt
       as temporary storage for the matrix product. */
    if( body_force ) {

      /* Call the general function to carry out the matrix product. */
      config->GEMM_Tick(&tick);
      DenseMatrixProduct(nDOFs, nVar, nInt, matBasisIntTrans, sources, solAndGradInt);
      config->GEMM_Tock(tick, "Volume_Residual3", nDOFs, nVar, nInt);

      /* Add the residuals due to source terms to the volume residuals */
      for(unsigned short i=0; i<(nDOFs*nVar); ++i)
        res[i] += solAndGradInt[i];
    }
  }
}

void CFEM_DG_NSSolver::ResidualFaces(CConfig             *config,
                                     const unsigned long indFaceBeg,
                                     const unsigned long indFaceEnd,
                                     unsigned long       &indResFaces,
                                     CNumerics           *numerics) {

  /*--- Set the pointers for the local arrays. ---*/
  double tick = 0.0;
  double tick2 = 0.0;

  unsigned int sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, (unsigned int) nDOFsMax);

  const unsigned int sizeGradSolInt = nIntegrationMax*nDim*nVar;

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
    const su2double     *derBasisElem = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationSide0();

    /* Get the length scales of the adjacent elements. */
    const unsigned long elemID0 = matchingInternalFaces[l].elemID0;
    const unsigned long elemID1 = matchingInternalFaces[l].elemID1;

    const su2double lenScale0 = volElem[elemID0].lenScale;
    const su2double lenScale1 = volElem[elemID1].lenScale;

    /* Determine the time level of the face, which is the minimum
       of the levels of the adjacent elements. */
    const unsigned short timeLevelFace = min(volElem[elemID0].timeLevel,
                                             volElem[elemID1].timeLevel);

    /* Call the general function to compute the viscous flux in normal
       direction for side 0. */
    ViscousNormalFluxFace(config, &volElem[elemID0], timeLevelFace, nInt,
                          0.0, false, derBasisElem, solIntL,
                          matchingInternalFaces[l].DOFsSolElementSide0.data(),
                          matchingInternalFaces[l].metricCoorDerivFace0.data(),
                          matchingInternalFaces[l].metricNormalsFace.data(),
                          matchingInternalFaces[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityIntL, kOverCvIntL);

    /*--- Subtract half of the viscous fluxes from the inviscid fluxes. The
          factor 0.5 comes from the fact that the average of the viscous fluxes
          of side 0 and side 1 must be taken in the DG-FEM formulation. ---*/
    for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] -= 0.5*viscFluxes[j];

    /*--- Get the information from the standard face to compute the viscous
          fluxes in the integration points on the right side, i.e. side 1. ---*/
    derBasisElem = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationSide1();

    /* Call the general function to compute the viscous flux in normal
       direction for side 1. */
    ViscousNormalFluxFace(config, &volElem[elemID1], timeLevelFace, nInt,
                          0.0, false, derBasisElem, solIntR,
                          matchingInternalFaces[l].DOFsSolElementSide1.data(),
                          matchingInternalFaces[l].metricCoorDerivFace1.data(),
                          matchingInternalFaces[l].metricNormalsFace.data(),
                          matchingInternalFaces[l].wallDistance.data(),
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
      su2double *flux = fluxes + i*nVar;

      for(unsigned short j=0; j<nVar; ++j)
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
    su2double *resFace0 = VecResFaces.data() + indResFaces*nVar;
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
    /*--- Step 5: Compute and distribute the symmetrizing terms, if        ---*/
    /*---         present. Note that these terms must be distributed to    ---*/
    /*---         DOFs of the adjacent element, not only of the face.      ---*/
    /*------------------------------------------------------------------------*/

    if( symmetrizingTermsPresent ) {

      config->Tick(&tick);
      /* Compute the symmetrizing fluxes in the nDim directions in the
         integration points of the face. */
      SymmetrizingFluxesFace(nInt, solIntL, solIntR, viscosityIntL, viscosityIntR,
                             kOverCvIntL, kOverCvIntR,
                             matchingInternalFaces[l].metricNormalsFace.data(),
                             fluxes);

      /* Transform the fluxes, such that they must be multiplied with the
         gradients w.r.t. the parametric coordinates rather than the
         Cartesian coordinates of the basis functions. Also a multiplication
         with the integration weight and the theta parameter is carried out.
         Use gradSolInt as a buffer to store the parametric fluxes.
         First side 0 of the face.  */
      su2double *paramFluxes = gradSolInt;
      const su2double halfTheta = 0.5*config->GetTheta_Interior_Penalty_DGFEM();

      TransformSymmetrizingFluxes(nInt, halfTheta, fluxes, weights,
                                  matchingInternalFaces[l].metricCoorDerivFace0.data(),
                                  paramFluxes);

      /* Easier storage of the number of DOFs for the adjacent elements. */
      const unsigned short nDOFsElem0 = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
      const unsigned short nDOFsElem1 = standardMatchingFacesSol[ind].GetNDOFsElemSide1();

      /* Set the pointer where to store the current residual and update the
         counter indResFaces. */
      su2double *resElem0 = VecResFaces.data() + indResFaces*nVar;
      indResFaces        += nDOFsElem0;

      /* Call the general function to carry out the matrix product to compute
         the residual for side 0. */
      const su2double *derBasisElemTrans = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationTransposeSide0();
      config->GEMM_Tick(&tick2);
      DenseMatrixProduct(nDOFsElem0, nVar, nInt*nDim, derBasisElemTrans, paramFluxes, resElem0);
      config->GEMM_Tock(tick2, "ResidualFaces3", nDOFsElem0, nVar, nInt*nDim);

      /* Transform the fluxes, such that they must be multiplied with the
         gradients w.r.t. the parametric coordinates rather than the
         Cartesian coordinates of the basis functions. Also a multiplication
         with the integration weight and the theta parameter is carried out.
         Use gradSolInt as a buffer to store the parametric fluxes.
         Now side 1 of the face.  */
      TransformSymmetrizingFluxes(nInt, halfTheta, fluxes, weights,
                                  matchingInternalFaces[l].metricCoorDerivFace1.data(),
                                  paramFluxes);

      /* Set the pointer where to store the current residual and update the
         counter indResFaces. */
      su2double *resElem1 = VecResFaces.data() + indResFaces*nVar;
      indResFaces        += nDOFsElem1;

      /* Call the general function to carry out the matrix product to compute
         the residual for side 1. Note that the symmetrizing residual should not
         be negated, because two minus signs enter the formulation for side 1,
         which cancel each other. */
      derBasisElemTrans = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationTransposeSide1();
      config->GEMM_Tick(&tick2);
      DenseMatrixProduct(nDOFsElem1, nVar, nInt*nDim, derBasisElemTrans, paramFluxes, resElem1);
      config->GEMM_Tock(tick2, "ResidualFaces4", nDOFsElem1, nVar, nInt*nDim);
    }
    config->Tock(tick, "ER_1_5", 4);
  }
}

void CFEM_DG_NSSolver::ViscousNormalFluxFace(CConfig                 *config,
                                             const CVolumeElementFEM *adjVolElem,
                                             const unsigned short    timeLevelFace,
                                             const unsigned short    nInt,
                                             const su2double         Wall_HeatFlux,
                                             const bool              HeatFlux_Prescribed,
                                             const su2double         *derBasisElem,
                                             const su2double         *solInt,
                                             const unsigned long     *DOFsElem,
                                             const su2double         *metricCoorDerivFace,
                                             const su2double         *metricNormalsFace,
                                             const su2double         *wallDistanceInt,
                                                   su2double         *gradSolInt,
                                                   su2double         *viscNormFluxes,
                                                   su2double         *viscosityInt,
                                                   su2double         *kOverCvInt) {

  double tick = 0.0;

  /* Constant factor present in the heat flux vector. Set it to zero if the heat
     flux is prescribed, such that no if statements are needed in the loop. */
  const su2double factHeatFlux_Lam  = HeatFlux_Prescribed ? su2double(0.0): Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = HeatFlux_Prescribed ? su2double(0.0): Gamma/Prandtl_Turb;

  /* Set the value of the prescribed heat flux for the same reason. */
  const su2double HeatFlux = HeatFlux_Prescribed ? Wall_HeatFlux : su2double(0.0);

  /* Determine the number of DOFs of the adjacent element and the offset that
     must be applied to access the correct data for this element in the working
     vector of the solution. If the element has the same time level as the face
     offsetDOFsSolThisTimeLevel is used, otherwise offsetDOFsSolPrevTimeLevel
     is the correct value. */
  const unsigned short nDOFsElem = adjVolElem->nDOFsSol;
  unsigned long offset;
  if(adjVolElem->timeLevel == timeLevelFace)
    offset = adjVolElem->offsetDOFsSolLocal - adjVolElem->offsetDOFsSolThisTimeLevel;
  else
    offset = adjVolElem->offsetDOFsSolLocal - adjVolElem->offsetDOFsSolPrevTimeLevel;

  /* Compute the length scale for the LES of the adjacent element. */
  const unsigned short iind  = adjVolElem->indStandardElement;
  unsigned short       nPoly = standardElementsSol[iind].GetNPoly();
  if(nPoly == 0) nPoly = 1;

  const su2double lenScale_LES = adjVolElem->lenScale/nPoly;

  /* Set the pointer solElem to viscNormFluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the element and the fluxes to be computed. */
  su2double *solElem = viscNormFluxes;

  /*--- Store the solution of the DOFs of the adjacent element in contiguous
        memory such that the function DenseMatrixProduct can be used to compute
        the gradients solution variables in the integration points of the face. ---*/
  const unsigned long nBytes = nVar*sizeof(su2double);
  for(unsigned short i=0; i<nDOFsElem; ++i) {
    const su2double *solDOF = VecWorkSolDOFs[timeLevelFace].data()
                            + nVar*(DOFsElem[i] - offset);
    su2double       *sol    = solElem + nVar*i;
    memcpy(sol, solDOF, nBytes);
  }

  /* Compute the gradients in the integration points. Call the general function to
     carry out the matrix product. */
  config->GEMM_Tick(&tick);
  DenseMatrixProduct(nInt*nDim, nVar, nDOFsElem, derBasisElem, solElem, gradSolInt);
  config->GEMM_Tock(tick, "ViscousNormalFluxFace", nInt*nDim, nVar, nDOFsElem);

  /* Determine the offset between r- and -s-derivatives, which is also the
     offset between s- and t-derivatives. */
  const unsigned short offDeriv = nVar*nInt;

  /* Make a distinction between two and three space dimensions
     in order to have the most efficient code. */
  switch( nDim ) {

    case 2: {

      /* 2D simulation. Loop over the integration points to
         compute the viscous fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the metric terms needed to compute the Cartesian
           gradients in this integration point and the starting locations of
           the solution and the gradients, w.r.t. the parametric coordinates
           of this solution. */
        const su2double *metricTerms = metricCoorDerivFace + i*nDim*nDim;
        const su2double *sol         = solInt     + nVar*i;
        const su2double *dSolDr      = gradSolInt + nVar*i;
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
        const su2double *normal  = metricNormalsFace + i*(nDim+1);
        su2double *normalFlux    = viscNormFluxes + i*nVar;
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
        const su2double *metricTerms = metricCoorDerivFace + i*nDim*nDim;
        const su2double *sol         = solInt     + nVar*i;
        const su2double *dSolDr      = gradSolInt + nVar*i;
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
        const su2double *normal  = metricNormalsFace + i*(nDim+1);
        su2double *normalFlux    = viscNormFluxes + i*nVar;
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
  if( SGSModelUsed )
    ViscosityTurb = SGSModel->ComputeEddyViscosity_2D(sol[0], dudx, dudy, dvdx,
                                                      dvdy, lenScale_LES, wallDist);

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
  if( SGSModelUsed )
    ViscosityTurb = SGSModel->ComputeEddyViscosity_3D(sol[0], dudx, dudy, dudz,
                                                      dvdx, dvdy, dvdz, dwdx,
                                                      dwdy, dwdz, lenScale_LES,
                                                      wallDist);

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

  /*--- Make a distinction between 2D and 3D for efficiency. ---*/
  switch ( nDim ) {
    case 2: {

      /* 2D simulation. Loop over the integration points to compute
         the penalty fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const su2double *sol0   = solInt0 + 4*i;           // nVar*i
        const su2double *sol1   = solInt1 + 4*i;           // nVar*i
        const su2double *normal = metricNormalsFace + 3*i; // (nDim+1)*i
        su2double       *flux   = penaltyFluxes + 4*i;     // nVar*i

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

        const su2double penFace = normal[2]*max(pen0, pen1);

        /* Compute the penalty flux, where it is assumed that the normal points from
           side 0 to side 1. */
        flux[0] = penFace*(sol0[0] - sol1[0]);
        flux[1] = penFace*(sol0[1] - sol1[1]);
        flux[2] = penFace*(sol0[2] - sol1[2]);
        flux[3] = penFace*(sol0[3] - sol1[3]);
      }

      break;
    }

    /*----------------------------------------------------------------------*/

    case 3: {

      /* 3D simulation. Loop over the integration points to compute
         the penalty fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const su2double *sol0   = solInt0 + 5*i;           // nVar*i
        const su2double *sol1   = solInt1 + 5*i;           // nVar*i
        const su2double *normal = metricNormalsFace + 4*i; // (nDim+1)*i
        su2double       *flux   = penaltyFluxes + 5*i;     // nVar*i

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

        const su2double penFace = normal[3]*max(pen0, pen1);

        /* Compute the penalty flux, where it is assumed that the normal points from
           side 0 to side 1. */
        flux[0] = penFace*(sol0[0] - sol1[0]);
        flux[1] = penFace*(sol0[1] - sol1[1]);
        flux[2] = penFace*(sol0[2] - sol1[2]);
        flux[3] = penFace*(sol0[3] - sol1[3]);
        flux[4] = penFace*(sol0[4] - sol1[4]);
      }

      break;
    }
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
  switch( nDim ) {

    case 2: {

      /*--- 2D simulation. Loop over the number of integration points
            of the face. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the variables for this integration point. */
        const su2double *sol0   = solInt0 + nVar*i;
        const su2double *sol1   = solInt1 + nVar*i;
        const su2double *normal = metricNormalsFace + i*(nDim+1);
        su2double       *flux   = symmFluxes + i*nDim*nVar;

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
        const su2double *sol0   = solInt0 + nVar*i;
        const su2double *sol1   = solInt1 + nVar*i;
        const su2double *normal = metricNormalsFace + i*(nDim+1);
        su2double       *flux   = symmFluxes + i*nDim*nVar;

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

void CFEM_DG_NSSolver::TransformSymmetrizingFluxes(const unsigned short nInt,
                                                   const su2double      halfTheta,
                                                   const su2double      *symmFluxes,
                                                   const su2double      *weights,
                                                   const su2double      *metricCoorFace,
                                                         su2double      *paramFluxes) {

  /*--- Transform the fluxes, such that they must be multiplied with the
        gradients w.r.t. the parametric coordinates rather than the
        Cartesian coordinates of the basis functions. This involves the
        multiplication with the metric terms. Also multiply the fluxes with
        their integration weights and -theta/2. The parameter theta is the
        parameter in the Interior Penalty formulation, the factor 1/2 comes in
        from the averaging and the minus sign is from the convention that the
        viscous fluxes come with a minus sign in this code. ---*/

  /* Make a distinction between two and three space dimensions
     in order to have the most efficient code. */
  switch( nDim ) {

    case 2: {
      /* Two-dimensional computation. Loop over the integration
         points to correct the fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the location where the the fluxes
           are stored for this integration point as well as an
           abbreviation for the integration weights time halfTheta. */
        const su2double *flux  =  symmFluxes + i*nVar*nDim;
        const su2double wTheta = -halfTheta*weights[i];
        su2double *paramFlux   =  paramFluxes + i*nVar*nDim;

        /* Compute the modified metric terms. */
        const su2double *metricTerms = metricCoorFace + 4*i;   // The 4 is nDim*nDim;
        const su2double drdx = wTheta*metricTerms[0];
        const su2double drdy = wTheta*metricTerms[1];
        const su2double dsdx = wTheta*metricTerms[2];
        const su2double dsdy = wTheta*metricTerms[3];

        /* Parametric fluxes in r-direction. */
        paramFlux[0] = flux[0]*drdx + flux[4]*drdy;
        paramFlux[1] = flux[1]*drdx + flux[5]*drdy;
        paramFlux[2] = flux[2]*drdx + flux[6]*drdy;
        paramFlux[3] = flux[3]*drdx + flux[7]*drdy;

        /* Parametric fluxes in s-direction. */
        paramFlux[4] = flux[0]*dsdx + flux[4]*dsdy;
        paramFlux[5] = flux[1]*dsdx + flux[5]*dsdy;
        paramFlux[6] = flux[2]*dsdx + flux[6]*dsdy;
        paramFlux[7] = flux[3]*dsdx + flux[7]*dsdy;
      }

      break;
    }

    case 3: {
      /* Three-dimensional computation. Loop over the integration
         points to correct the fluxes. */
      for(unsigned short i=0; i<nInt; ++i) {

        /* Easier storage of the location where the the fluxes
           are stored for this integration point as well as an
           abbreviation for the integration weights time halfTheta. */
        const su2double *flux  =  symmFluxes + i*nVar*nDim;
        const su2double wTheta = -halfTheta*weights[i];
        su2double *paramFlux   =  paramFluxes + i*nVar*nDim;

        /* Compute the modified metric terms. */
        const su2double *metricTerms = metricCoorFace + 9*i;   // The 9 is nDim*nDim;
        const su2double drdx = wTheta*metricTerms[0];
        const su2double drdy = wTheta*metricTerms[1];
        const su2double drdz = wTheta*metricTerms[2];

        const su2double dsdx = wTheta*metricTerms[3];
        const su2double dsdy = wTheta*metricTerms[4];
        const su2double dsdz = wTheta*metricTerms[5];

        const su2double dtdx = wTheta*metricTerms[6];
        const su2double dtdy = wTheta*metricTerms[7];
        const su2double dtdz = wTheta*metricTerms[8];

        /* Parametric fluxes in r-direction. */
        paramFlux[0] = flux[0]*drdx + flux[5]*drdy + flux[10]*drdz;
        paramFlux[1] = flux[1]*drdx + flux[6]*drdy + flux[11]*drdz;
        paramFlux[2] = flux[2]*drdx + flux[7]*drdy + flux[12]*drdz;
        paramFlux[3] = flux[3]*drdx + flux[8]*drdy + flux[13]*drdz;
        paramFlux[4] = flux[4]*drdx + flux[9]*drdy + flux[14]*drdz;

        /* Parametric fluxes in s-direction. */
        paramFlux[5] = flux[0]*dsdx + flux[5]*dsdy + flux[10]*dsdz;
        paramFlux[6] = flux[1]*dsdx + flux[6]*dsdy + flux[11]*dsdz;
        paramFlux[7] = flux[2]*dsdx + flux[7]*dsdy + flux[12]*dsdz;
        paramFlux[8] = flux[3]*dsdx + flux[8]*dsdy + flux[13]*dsdz;
        paramFlux[9] = flux[4]*dsdx + flux[9]*dsdy + flux[14]*dsdz;

        /* Parametric fluxes in t-direction. */
        paramFlux[10] = flux[0]*dtdx + flux[5]*dtdy + flux[10]*dtdz;
        paramFlux[11] = flux[1]*dtdx + flux[6]*dtdy + flux[11]*dtdz;
        paramFlux[12] = flux[2]*dtdx + flux[7]*dtdy + flux[12]*dtdz;
        paramFlux[13] = flux[3]*dtdx + flux[8]*dtdy + flux[13]*dtdz;
        paramFlux[14] = flux[4]*dtdx + flux[9]*dtdy + flux[14]*dtdz;
      }

      break;
    }
  }
}


void CFEM_DG_NSSolver::BC_Euler_Wall(CConfig                  *config,
                                     const unsigned long      surfElemBeg,
                                     const unsigned long      surfElemEnd,
                                     const CSurfaceElementFEM *surfElem,
                                     su2double                *resFaces,
                                     CNumerics                *conv_numerics){

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

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the right state by applying the inviscid wall BC's. */
    BoundaryStates_Euler_Wall(config, &surfElem[l], solIntL, solIntR);

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned long  elemID    = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short ind      = surfElem[l].indStandardElement;
    const unsigned short nInt     = standardBoundaryFacesSol[ind].GetNIntegration();
    const su2double *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, &volElem[elemID], timeLevel, nInt,
                          0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Far_Field(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics){

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

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

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

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned long  elemID    = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const su2double *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, &volElem[elemID], timeLevel, nInt,
                          0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Sym_Plane(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics){
  double tick = 0.0;

  /* Constant factor present in the heat flux vector, namely the ratio of
     thermal conductivity and viscosity. */
  const su2double factHeatFlux_Lam  = Gamma/Prandtl_Lam;
  const su2double factHeatFlux_Turb = Gamma/Prandtl_Turb;

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

  /* Set the pointer solElem to fluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the element and the fluxes to be computed. */
  su2double *solElem = fluxes;

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

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

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Determine the offset that must be applied to access the correct data for
       the adjacent element in the working vector of the solution. This is a
       boundary face, which has the same time level as the adjacent element. */
    const unsigned long offset = volElem[elemID].offsetDOFsSolLocal
                               - volElem[elemID].offsetDOFsSolThisTimeLevel;

    /* Easier storage of the wall distance array for this surface element. */
    const su2double *wallDistance = surfElem[l].wallDistance.data();

    /* Determine the offset between r- and -s-derivatives, which is also the
       offset between s- and t-derivatives. */
    const unsigned short offDeriv = nVar*nInt;

    /* Store the solution of the DOFs of the adjacent element in contiguous
       memory such that the function DenseMatrixProduct can be used to compute
       the gradients solution variables in the integration points of the face. */
    const unsigned long nBytes = nVar*sizeof(su2double);
    for(unsigned short i=0; i<nDOFsElem; ++i) {
      const su2double *solDOF = VecWorkSolDOFs[timeLevel].data()
                              + nVar*(surfElem[l].DOFsSolElement[i] - offset);
      su2double       *sol    = solElem + nVar*i;
      memcpy(sol, solDOF, nBytes);
    }

    /* Compute the left gradients in the integration points. Call the general
       function to carry out the matrix product. */
    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nInt*nDim, nVar, nDOFsElem, derBasisElem, solElem, gradSolInt);
    config->GEMM_Tock(tick, "BC_Sym_Plane", nInt*nDim, nVar, nDOFsElem);

    /* Make a distinction between two and three space dimensions
       in order to have the most efficient code. */
    switch( nDim ) {

      case 2: {

        /*--- 2D simulation. Loop over the integration points to apply the
              symmetry condition and to compute the viscous normal fluxes. This
              is done in the same loop, because the gradients in the right
              integration points are constructed using the symmetry boundary
              condition. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution and the normals
             for this integration point. */
          const su2double *UL      = solIntL + i*nVar;
                su2double *UR      = solIntR + i*nVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

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
                                       + i*nDim*nDim;
          const su2double *dULDr       = gradSolInt + nVar*i;
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

          su2double *viscNormalFlux = viscFluxes + i*nVar;
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
          const su2double *UL      = solIntL + i*nVar;
                su2double *UR      = solIntR + i*nVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

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
                                       + i*nDim*nDim;
          const su2double *dULDr       = gradSolInt + nVar*i;
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

          su2double *viscNormalFlux = viscFluxes + i*nVar;
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

void CFEM_DG_NSSolver::BC_Inlet(CConfig                  *config,
                                const unsigned long      surfElemBeg,
                                const unsigned long      surfElemEnd,
                                const CSurfaceElementFEM *surfElem,
                                su2double                *resFaces,
                                CNumerics                *conv_numerics,
                                unsigned short           val_marker) {

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

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the right state by applying the subsonic inlet BC's. */
    BoundaryStates_Inlet(config, &surfElem[l], val_marker, solIntL, solIntR);

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned long  elemID    = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short ind      = surfElem[l].indStandardElement;
    const unsigned short nInt     = standardBoundaryFacesSol[ind].GetNIntegration();
    const su2double *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, &volElem[elemID], timeLevel, nInt,
                          0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Outlet(CConfig                  *config,
                                 const unsigned long      surfElemBeg,
                                 const unsigned long      surfElemEnd,
                                 const CSurfaceElementFEM *surfElem,
                                 su2double                *resFaces,
                                 CNumerics                *conv_numerics,
                                 unsigned short           val_marker) {

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

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the right state by applying the subsonic outlet BC's. */
    BoundaryStates_Outlet(config, &surfElem[l], val_marker, solIntL, solIntR);

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned long  elemID    = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short ind      = surfElem[l].indStandardElement;
    const unsigned short nInt     = standardBoundaryFacesSol[ind].GetNIntegration();
    const su2double *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, &volElem[elemID], timeLevel, nInt,
                          0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_HeatFlux_Wall(CConfig                  *config,
                                        const unsigned long      surfElemBeg,
                                        const unsigned long      surfElemEnd,
                                        const CSurfaceElementFEM *surfElem,
                                        su2double                *resFaces,
                                        CNumerics                *conv_numerics,
                                        unsigned short           val_marker) {

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

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Determine the number of integration points. */
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Apply the heat flux wall boundary conditions to compute the right
          state. There are two options. Either the velocity is negated or it
          is set to zero. Some experiments are needed to see which formulation
          gives better results. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the grid velocity and the left and right solution
         for this integration point. */
      const su2double *gridVel = surfElem[l].gridVelocities.data() + i*nDim;
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;

      /* Set the right state. The initial value of the total energy is the
         energy of the left state. Also compute the difference in kinetic
         energy between the left and right state. */
      UR[0]      = UL[0];
      UR[nDim+1] = UL[nDim+1];

      su2double DensityInv = 1.0/UL[0];
      su2double diffKin    = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        const su2double velL = DensityInv*UL[iDim+1];
        const su2double dV   = factWallVel*(velL-gridVel[iDim]);
        const su2double velR = gridVel[iDim] - dV;

        UR[iDim+1] = UR[0]*velR;
        diffKin   += velL*velL - velR*velR;
      }

      /* As only the internal energy of UR is equal to UL, the difference
         in kinetic energy must be subtracted. */
      UR[nDim+1] -= 0.5*UR[0]*diffKin;
    }

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned long  elemID    = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const su2double *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, &volElem[elemID], timeLevel, nInt,
                          Wall_HeatFlux, true, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Isothermal_Wall(CConfig                  *config,
                                          const unsigned long      surfElemBeg,
                                          const unsigned long      surfElemEnd,
                                          const CSurfaceElementFEM *surfElem,
                                          su2double                *resFaces,
                                          CNumerics                *conv_numerics,
                                          unsigned short           val_marker) {

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

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Determine the number of integration points. */
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Apply the heat flux wall boundary conditions to compute the right
          state. There are two options. Either the velocity is negated or it
          is set to zero. Some experiments are needed to see which formulation
          gives better results. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the grid velocity and the left and right solution
         for this integration point. */
      const su2double *gridVel = surfElem[l].gridVelocities.data() + i*nDim;
      const su2double *UL      = solIntL + i*nVar;
            su2double *UR      = solIntR + i*nVar;

      /* Set the right state for the density and the momentum variables of the
         right state. Also compute twice the kinetic energy of the right state. */
      UR[0] = UL[0];
      su2double DensityInv = 1.0/UL[0];
      su2double kinEner    = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        const su2double velL = DensityInv*UL[iDim+1];
        const su2double dV   = factWallVel*(velL-gridVel[iDim]);
        const su2double velR = gridVel[iDim] - dV;

        UR[iDim+1] = UR[0]*velR;
        kinEner   += velR*velR;
      }

      /* Compute the total energy of the right state. */
      UR[nDim+1] = UR[0]*(StaticEnergy + 0.5*kinEner);
    }

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned long  elemID    = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const su2double *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, &volElem[elemID], timeLevel, nInt,
                          0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Riemann(CConfig                  *config,
                                  const unsigned long      surfElemBeg,
                                  const unsigned long      surfElemEnd,
                                  const CSurfaceElementFEM *surfElem,
                                  su2double                *resFaces,
                                  CNumerics                *conv_numerics,
                                  unsigned short           val_marker) {

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

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

    /* Compute the right state by applying the Riemann BC's. */
    BoundaryStates_Riemann(config, &surfElem[l], val_marker, solIntL, solIntR);

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned long  elemID    = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const unsigned short ind      = surfElem[l].indStandardElement;
    const unsigned short nInt     = standardBoundaryFacesSol[ind].GetNIntegration();
    const su2double *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, &volElem[elemID], timeLevel, nInt,
                          0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
  }
}

void CFEM_DG_NSSolver::BC_Custom(CConfig                  *config,
                                 const unsigned long      surfElemBeg,
                                 const unsigned long      surfElemEnd,
                                 const CSurfaceElementFEM *surfElem,
                                 su2double                *resFaces,
                                 CNumerics                *conv_numerics) {

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

  /*--- Loop over the given range of boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, &surfElem[l], fluxes, solIntL);

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

      SU2_MPI::Error("No or wrong compiler directive specified. This is necessary for customized boundary conditions.",
                     CURRENT_FUNCTION);

#endif
    }

    /* Determine the time level of the boundary face, which is equal
       to the time level of the adjacent element. */
    const unsigned long  elemID    = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[elemID].timeLevel;

    /* Call the general function to compute the viscous flux in normal
       direction for the face. */
    const su2double *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();

    ViscousNormalFluxFace(config, &volElem[elemID], timeLevel, nInt,
                          0.0, false, derBasisElem, solIntL,
                          surfElem[l].DOFsSolElement.data(),
                          surfElem[l].metricCoorDerivFace.data(),
                          surfElem[l].metricNormalsFace.data(),
                          surfElem[l].wallDistance.data(),
                          gradSolInt, viscFluxes, viscosityInt, kOverCvInt);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], solIntL,
                                solIntR, gradSolInt, fluxes, viscFluxes,
                                viscosityInt, kOverCvInt, resFaces, indResFaces);
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

  double tick = 0.0;

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
                            surfElem->gridVelocities.data(),
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
  /*--- Step 3: Compute and distribute the symmetrizing terms, if        ---*/
  /*---         present. Note that these terms must be distributed to    ---*/
  /*---         DOFs of the adjacent element, not only of the face.      ---*/
  /*------------------------------------------------------------------------*/

  if( symmetrizingTermsPresent ) {

    /* Compute the symmetrizing fluxes in the nDim directions. */
    SymmetrizingFluxesFace(nInt, solInt0, solInt1, viscosityInt,
                           viscosityInt, kOverCvInt, kOverCvInt,
                           surfElem->metricNormalsFace.data(), fluxes);

    /* Transform the fluxes, such that they must be multiplied with the
       gradients w.r.t. the parametric coordinates rather than the
       Cartesian coordinates of the basis functions. Also a multiplication
       with the integration weight and the theta parameter is carried out.
       Use gradSolInt as a buffer to store the parametric fluxes. */
    su2double *paramFluxes = gradSolInt;
    const su2double halfTheta = 0.5*config->GetTheta_Interior_Penalty_DGFEM();

    TransformSymmetrizingFluxes(nInt, halfTheta, fluxes, weights,
                                surfElem->metricCoorDerivFace.data(),
                                paramFluxes);

    /* Easier storage of the position in the residual array for this face
       and update the corresponding counter. */
    su2double *resElem = resFaces + indResFaces*nVar;
    indResFaces       += nDOFsElem;

    /* Call the general function to carry out the matrix product to compute
       the residual. */
    const su2double *derBasisElemTrans = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegrationTranspose();

    config->GEMM_Tick(&tick);
    DenseMatrixProduct(nDOFsElem, nVar, nInt*nDim, derBasisElemTrans, paramFluxes, resElem);
    config->GEMM_Tock(tick, "ResidualViscousBoundaryFace2", nDOFsElem, nVar, nInt*nDim);
  }
}
