/*!
 * \file CFEM_DG_EulerSolver.cpp
 * \brief Main subroutines for solving finite element Euler flow problems
 * \author J. Alonso, E. van der Weide, T. Economon
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


#include "../../include/solvers/CFEM_DG_EulerSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CIdealGas.hpp"
#include "../../include/fluid/CVanDerWaalsGas.hpp"
#include "../../include/fluid/CPengRobinson.hpp"
#include "../../include/fluid/CCoolProp.hpp"
#include "../../include/fluid/CDataDrivenFluid.hpp"

enum {
SIZE_ARR_NORM = 8
};

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver() : CSolver() {

  /*--- Basic array initialization ---*/

  FluidModel = nullptr;

  CD_Inv = nullptr; CL_Inv = nullptr; CSF_Inv = nullptr;  CEff_Inv = nullptr;
  CMx_Inv = nullptr; CMy_Inv = nullptr; CMz_Inv = nullptr;
  CFx_Inv = nullptr; CFy_Inv = nullptr; CFz_Inv = nullptr;

  /*--- Surface-based array initialization ---*/
  Surface_CL_Inv = nullptr; Surface_CD_Inv = nullptr; Surface_CSF_Inv = nullptr; Surface_CEff_Inv = nullptr;
  Surface_CFx_Inv = nullptr; Surface_CFy_Inv = nullptr; Surface_CFz_Inv = nullptr;
  Surface_CMx_Inv = nullptr; Surface_CMy_Inv = nullptr; Surface_CMz_Inv = nullptr;

  Surface_CL = nullptr; Surface_CD = nullptr; Surface_CSF = nullptr; Surface_CEff = nullptr;
  Surface_CFx = nullptr; Surface_CFy = nullptr; Surface_CFz = nullptr;
  Surface_CMx = nullptr; Surface_CMy = nullptr; Surface_CMz = nullptr;

  /*--- Initialization of the boolean symmetrizingTermsPresent. ---*/
  symmetrizingTermsPresent = true;

  /*--- Initialize the variables for the monitoring data to avoid
        possible valgrind warnings. ---*/
  Total_CL   = Total_CD  = Total_CSF = 0.0;
  Total_CFx  = Total_CFy = Total_CFz = 0.0;
  Total_CMx  = Total_CMy = Total_CMz = 0.0;
  Total_CEff = 0.0;

  /*--- Initialize the pointer for performing the BLAS functionalities. ---*/
  blasFunctions = nullptr;

}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CConfig *config, unsigned short val_nDim, unsigned short iMesh) : CSolver() {

  /*--- Basic array initialization ---*/

  FluidModel = nullptr;

  CD_Inv = nullptr; CL_Inv = nullptr; CSF_Inv = nullptr;  CEff_Inv = nullptr;
  CMx_Inv = nullptr; CMy_Inv = nullptr; CMz_Inv = nullptr;
  CFx_Inv = nullptr; CFy_Inv = nullptr; CFz_Inv = nullptr;

  /*--- Surface-based array initialization ---*/
  Surface_CL_Inv = nullptr; Surface_CD_Inv = nullptr; Surface_CSF_Inv = nullptr; Surface_CEff_Inv = nullptr;
  Surface_CFx_Inv = nullptr; Surface_CFy_Inv = nullptr; Surface_CFz_Inv = nullptr;
  Surface_CMx_Inv = nullptr; Surface_CMy_Inv = nullptr; Surface_CMz_Inv = nullptr;

  Surface_CL = nullptr; Surface_CD = nullptr; Surface_CSF = nullptr; Surface_CEff = nullptr;
  Surface_CFx = nullptr; Surface_CFy = nullptr; Surface_CFz = nullptr;
  Surface_CMx = nullptr; Surface_CMy = nullptr; Surface_CMz = nullptr;

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

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

  /*--- Initialize the variables for the monitoring data to avoid
        possible valgrind warnings. ---*/
  Total_CL   = Total_CD  = Total_CSF = 0.0;
  Total_CFx  = Total_CFy = Total_CFz = 0.0;
  Total_CMx  = Total_CMy = Total_CMz = 0.0;
  Total_CEff = 0.0;

  /*--- Initialize the pointer for performing the BLAS functionalities. ---*/
  blasFunctions = nullptr;

}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

  /*--- Array initialization ---*/
  FluidModel = nullptr;

  CD_Inv = nullptr; CL_Inv = nullptr; CSF_Inv = nullptr; CEff_Inv = nullptr;
  CMx_Inv = nullptr;   CMy_Inv = nullptr;   CMz_Inv = nullptr;
  CFx_Inv = nullptr;   CFy_Inv = nullptr;   CFz_Inv = nullptr;

  Surface_CL_Inv = nullptr; Surface_CD_Inv = nullptr; Surface_CSF_Inv = nullptr; Surface_CEff_Inv = nullptr;
  Surface_CFx_Inv = nullptr;   Surface_CFy_Inv = nullptr;   Surface_CFz_Inv = nullptr;
  Surface_CMx_Inv = nullptr;   Surface_CMy_Inv = nullptr;   Surface_CMz_Inv = nullptr;

  Surface_CL = nullptr; Surface_CD = nullptr; Surface_CSF = nullptr; Surface_CEff = nullptr;
  Surface_CFx = nullptr;   Surface_CFy = nullptr;   Surface_CFz = nullptr;
  Surface_CMx = nullptr;   Surface_CMy = nullptr;   Surface_CMz = nullptr;

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Allocate the memory for blasFunctions. ---*/
  blasFunctions = new CBlasStructure;

  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure. ---*/
  nDim    = geometry->GetnDim();
  nMarker = config->GetnMarker_All();

  const bool compressible = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  if( compressible ) nVar = nDim + 2;
  else               nVar = nDim + 1;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  auto *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

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

  timeCoefADER_DG                        = DGGeometry->GetTimeCoefADER_DG();
  timeInterpolDOFToIntegrationADER_DG    = DGGeometry->GetTimeInterpolDOFToIntegrationADER_DG();
  timeInterpolAdjDOFToIntegrationADER_DG = DGGeometry->GetTimeInterpolAdjDOFToIntegrationADER_DG();

  /*--- Determine the maximum number of integration points used. Usually this
        is for the volume integral, but to avoid problems the faces are also
        taken into account.  ---*/
  unsigned short nIntegrationMax = 0;
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
  unsigned short nDOFsMax = 0;
  for(unsigned short i=0; i<nStandardElementsSol; ++i) {
    const unsigned short nDOFs = standardElementsSol[i].GetNDOFs();
    nDOFsMax = max(nDOFsMax, nDOFs);
  }

  for(unsigned short i=0; i<nStandardMatchingFacesSol; ++i) {
    const unsigned short nDOFs = standardMatchingFacesSol[i].GetNDOFsElemSide1();
    nDOFsMax = max(nDOFsMax, nDOFs);
  }

  /*--- Determine the size of the work array. ---*/
  const unsigned short nPadGemm = config->GetSizeMatMulPadding();
  if( config->GetViscous() ) {

    /* Viscous simulation. */
    unsigned int sizeFluxes = nIntegrationMax*nDim;
    sizeFluxes = nPadGemm*max(sizeFluxes, (unsigned int) nDOFsMax);

    const unsigned int sizeGradSolInt = nIntegrationMax*nDim*max(nPadGemm,nDOFsMax);

    sizeWorkArray = nIntegrationMax*(4 + 3*nPadGemm) + sizeFluxes + sizeGradSolInt
                  + max(nIntegrationMax,nDOFsMax)*nPadGemm + nPadGemm*nDOFsMax;
  }
  else {

    /* Inviscid simulation. */
    unsigned int sizeVol = nPadGemm*nIntegrationMax*(nDim+2) + nPadGemm*nDOFsMax;
    unsigned int sizeSur = nPadGemm*(2*nIntegrationMax + max(nIntegrationMax,nDOFsMax));

    sizeWorkArray = max(sizeVol, sizeSur);
  }

  if(config->GetKind_TimeIntScheme_Flow() == ADER_DG) {

    /*--- ADER-DG scheme. Determine the size needed for the predictor step
          and make sure that sizeWorkArray is big enough. This size depends
          whether an aliased or a non-aliased predictor step is used. Note
          that the size estimates are for viscous computations. ---*/
    const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();

    unsigned int sizePredictorADER = 4*nPadGemm*nDOFsMax*nTimeDOFs
                                   +   nPadGemm*nDOFsMax;

    if(config->GetKind_ADER_Predictor() == ADER_ALIASED_PREDICTOR)
      sizePredictorADER += nDim*nPadGemm*nDOFsMax + nDim*nDim*nPadGemm*nIntegrationMax;
    else
      sizePredictorADER += (nDim+1)*nPadGemm*max(nIntegrationMax, nDOFsMax)
                         + nDim*nDim*nPadGemm*max(nIntegrationMax,nDOFsMax);

    sizeWorkArray = max(sizeWorkArray, sizePredictorADER);
  }

  /*--- Perform the non-dimensionalization for the flow equations using the
        specified reference values. ---*/
  SetNondimensionalization(config, iMesh, true);

  /*--- Check if we are executing a verification case. If so, the
   VerificationSolution object will be instantiated for a particular
   option from the available library of verification solutions. Note
   that this is done after SetNondim(), as problem-specific initial
   parameters are needed by the solution constructors. ---*/

  SetVerificationSolution(nDim, nVar, config);

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual_RMS.resize(nVar,1.e-35);
  Residual_Max.resize(nVar,1.e-35);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

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
  SU2_MPI::Allreduce(&nDOFsLocOwned, &nDOFsGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
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
    CGraphColoringStructure graphColoring;
    vector<int> colorLocalDOFs;
    graphColoring.GraphVertexColoring(config, nDOFsPerRank, nonZeroEntriesJacobian,
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

  /*--- Add the solver name. ---*/
  SolverName = "DG.FLOW";
}

CFEM_DG_EulerSolver::~CFEM_DG_EulerSolver() {

  delete FluidModel;
  delete blasFunctions;

  /*--- Array deallocation ---*/
            delete [] CD_Inv;
            delete [] CL_Inv;
           delete [] CSF_Inv;
           delete [] CMx_Inv;
           delete [] CMy_Inv;
           delete [] CMz_Inv;
           delete [] CFx_Inv;
           delete [] CFy_Inv;
           delete [] CFz_Inv;
    delete [] Surface_CL_Inv;
    delete [] Surface_CD_Inv;
   delete [] Surface_CSF_Inv;
  delete [] Surface_CEff_Inv;
   delete [] Surface_CFx_Inv;
   delete [] Surface_CFy_Inv;
   delete [] Surface_CFz_Inv;
   delete [] Surface_CMx_Inv;
   delete [] Surface_CMy_Inv;
   delete [] Surface_CMz_Inv;
        delete [] Surface_CL;
        delete [] Surface_CD;
       delete [] Surface_CSF;
      delete [] Surface_CEff;
       delete [] Surface_CFx;
       delete [] Surface_CFy;
       delete [] Surface_CFz;
       delete [] Surface_CMx;
       delete [] Surface_CMy;
       delete [] Surface_CMz;
          delete [] CEff_Inv;

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
  bool unsteady           = (config->GetTime_Marching() != TIME_MARCHING::STEADY);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool turbulent          = (config->GetKind_Solver() == MAIN_SOLVER::FEM_RANS) || (config->GetKind_Solver() == MAIN_SOLVER::FEM_LES);
  bool tkeNeeded          = ((turbulent) && ((config->GetKind_Turb_Model() == TURB_MODEL::SST) ));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == FREESTREAM_OPTION::TEMPERATURE_FS);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  /*--- The dimensional viscosity is needed to determine the free-stream conditions.
        To accomplish this, simply set the non-dimensional coefficients to the
        dimensional ones. This will be overruled later.---*/

  config->SetTemperature_Ref(1.0);
  config->SetViscosity_Ref(1.0);
  config->SetConductivity_Ref(1.0);

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

    case COOLPROP:

      FluidModel = new CCoolProp(config->GetFluid_Name());
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

    case DATADRIVEN_FLUID:
      FluidModel = new CDataDrivenFluid(config, false);
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
  Froude            = ModVel_FreeStream/sqrt(STANDARD_GRAVITY*Length_Ref);         config->SetFroude(Froude);

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

    case COOLPROP:
      FluidModel = new CCoolProp(config->GetFluid_Name());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case DATADRIVEN_FLUID:
      FluidModel = new CDataDrivenFluid(config);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
  }

  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

  if (viscous) {
    FluidModel->SetLaminarViscosityModel(config);
    FluidModel->SetThermalConductivityModel(config);
    FluidModel->SetMassDiffusivityModel(config); // nijso: TODO, needs to be tested
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

    stringstream NonDimTableOut, ModelTableOut;
    stringstream Unit;

    cout << endl;
    PrintingToolbox::CTablePrinter ModelTable(&ModelTableOut);
    ModelTableOut <<"-- Models:"<< endl;

    ModelTable.AddColumn("Viscosity Model", 25);
    ModelTable.AddColumn("Conductivity Model", 26);
    ModelTable.AddColumn("Fluid Model", 25);
    ModelTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
    ModelTable.PrintHeader();

    PrintingToolbox::CTablePrinter NonDimTable(&NonDimTableOut);
    NonDimTable.AddColumn("Name", 22);
    NonDimTable.AddColumn("Dim. value", 14);
    NonDimTable.AddColumn("Ref. value", 14);
    NonDimTable.AddColumn("Unit", 10);
    NonDimTable.AddColumn("Non-dim. value", 14);
    NonDimTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);

    NonDimTableOut <<"-- Fluid properties:"<< endl;

    NonDimTable.PrintHeader();

    if (viscous) {

      switch(config->GetKind_ViscosityModel()){
      case VISCOSITYMODEL::CONSTANT:
        ModelTable << "CONSTANT_VISCOSITY";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << config->GetMu_Constant() << config->GetMu_Constant()/config->GetMu_ConstantND() << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case VISCOSITYMODEL::COOLPROP:
        ModelTable << "COOLPROP";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << "--" << "--" << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case VISCOSITYMODEL::SUTHERLAND:
        ModelTable << "SUTHERLAND";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Ref. Viscosity" <<  config->GetMu_Ref() <<  config->GetViscosity_Ref() << Unit.str() << config->GetMu_RefND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "K";
        else if (config->GetSystemMeasurements() == US) Unit << "R";
        NonDimTable << "Sutherland Temp." << config->GetMu_Temperature_Ref() <<  config->GetTemperature_Ref() << Unit.str() << config->GetMu_Temperature_RefND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "K";
        else if (config->GetSystemMeasurements() == US) Unit << "R";
        NonDimTable << "Sutherland Const." << config->GetMu_S() << config->GetTemperature_Ref() << Unit.str() << config->GetMu_SND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      default:
        break;

      }
      switch(config->GetKind_ConductivityModel()){
      case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
        ModelTable << "CONSTANT_PRANDTL";
        NonDimTable << "Prandtl (Lam.)"  << "-" << "-" << "-" << config->GetPrandtl_Lam();
        Unit.str("");
        NonDimTable << "Prandtl (Turb.)" << "-" << "-" << "-" << config->GetPrandtl_Turb();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case CONDUCTIVITYMODEL::CONSTANT:
        ModelTable << "CONSTANT";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << config->GetThermal_Conductivity_Constant() << config->GetThermal_Conductivity_Constant()/config->GetThermal_Conductivity_ConstantND() << Unit.str() << config->GetThermal_Conductivity_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case CONDUCTIVITYMODEL::COOLPROP:
        ModelTable << "COOLPROP";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << "--" << "--" << Unit.str() << config->GetThermal_Conductivity_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      default:
        break;

      }
    } else {
      ModelTable << "-" << "-";
    }

    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    if (config->GetKind_FluidModel() == COOLPROP) {
      CCoolProp auxFluidModel(config->GetFluid_Name());
      NonDimTable << "Gas Constant" << auxFluidModel.GetGas_Constant() << config->GetGas_Constant_Ref()
                  << Unit.str() << auxFluidModel.GetGas_Constant()/config->GetGas_Constant_Ref();
    }
    else {
        NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
    }
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    if (config->GetKind_FluidModel() == COOLPROP) {
      NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << "-";
    }
    else {
      NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << Gamma;
    }
    Unit.str("");

    switch(config->GetKind_FluidModel()){
    case STANDARD_AIR:
      ModelTable << "STANDARD_AIR";
      break;
    case IDEAL_GAS:
      ModelTable << "IDEAL_GAS";
      break;
    case VW_GAS:
      ModelTable << "VW_GAS";
      break;
    case PR_GAS:
      ModelTable << "PR_GAS";
      break;
    case COOLPROP:
      ModelTable << "CoolProp library";
      break;
    }

    if (config->GetKind_FluidModel() == VW_GAS || config->GetKind_FluidModel() == PR_GAS){
        NonDimTable << "Critical Pressure" << config->GetPressure_Critical() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_Critical() /config->GetPressure_Ref();
        Unit.str("");
        Unit << "K";
        NonDimTable << "Critical Temperature" << config->GetTemperature_Critical() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_Critical() /config->GetTemperature_Ref();
        Unit.str("");
    }
    if (config->GetKind_FluidModel() == COOLPROP) {
        CCoolProp auxFluidModel(config->GetFluid_Name());
        NonDimTable << "Critical Pressure" << auxFluidModel.GetPressure_Critical() << config->GetPressure_Ref() << Unit.str() << auxFluidModel.GetPressure_Critical() /config->GetPressure_Ref();
        Unit.str("");
        Unit << "K";
        NonDimTable << "Critical Temperature" << auxFluidModel.GetTemperature_Critical() << config->GetTemperature_Ref() << Unit.str() << auxFluidModel.GetTemperature_Critical() /config->GetTemperature_Ref();
        Unit.str("");
    }
    NonDimTable.PrintFooter();

    NonDimTableOut <<"-- Initial and free-stream conditions:"<< endl;

    NonDimTable.PrintHeader();

    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Static Pressure" << config->GetPressure_FreeStream() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "kg/m^3";
    else if (config->GetSystemMeasurements() == US) Unit << "slug/ft^3";
    NonDimTable << "Density" << config->GetDensity_FreeStream() << config->GetDensity_Ref() << Unit.str() << config->GetDensity_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "K";
    else if (config->GetSystemMeasurements() == US) Unit << "R";
    NonDimTable << "Temperature" << config->GetTemperature_FreeStream() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
    else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
    NonDimTable << "Total Energy" << config->GetEnergy_FreeStream() << config->GetEnergy_Ref() << Unit.str() << config->GetEnergy_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m/s";
    else if (config->GetSystemMeasurements() == US) Unit << "ft/s";
    NonDimTable << "Velocity-X" << config->GetVelocity_FreeStream()[0] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[0];
    NonDimTable << "Velocity-Y" << config->GetVelocity_FreeStream()[1] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[1];
    if (nDim == 3) {
      NonDimTable << "Velocity-Z" << config->GetVelocity_FreeStream()[2] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[2];
    }
    NonDimTable << "Velocity Magnitude" << config->GetModVel_FreeStream() << config->GetVelocity_Ref() << Unit.str() << config->GetModVel_FreeStreamND();
    Unit.str("");

    if (viscous) {
      NonDimTable.PrintFooter();
      if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
      NonDimTable << "Viscosity" << config->GetViscosity_FreeStream() << config->GetViscosity_Ref() << Unit.str() << config->GetViscosity_FreeStreamND();
      Unit.str("");
      if      (config->GetSystemMeasurements() == SI) Unit << "W/m^2.K";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf/ft.s.R";
      NonDimTable << "Conductivity" << "-" << config->GetThermal_Conductivity_Ref() << Unit.str() << "-";
      Unit.str("");
      if (turbulent) {
        if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
        else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
        NonDimTable << "Turb. Kin. Energy" << config->GetTke_FreeStream() << config->GetTke_FreeStream()/config->GetTke_FreeStreamND() << Unit.str() << config->GetTke_FreeStreamND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "1/s";
        else if (config->GetSystemMeasurements() == US) Unit << "1/s";
        NonDimTable << "Spec. Dissipation" << config->GetOmega_FreeStream() << config->GetOmega_FreeStream()/config->GetOmega_FreeStreamND() << Unit.str() << config->GetOmega_FreeStreamND();
        Unit.str("");
      }
    }

    NonDimTable.PrintFooter();
    NonDimTable << "Mach Number" << "-" << "-" << "-" << config->GetMach();
    if (viscous) {
      NonDimTable << "Reynolds Number" << "-" << "-" << "-" << config->GetReynolds();
    }
    NonDimTable.PrintFooter();
    ModelTable.PrintFooter();

    if (unsteady){
      NonDimTableOut << "-- Unsteady conditions" << endl;
      NonDimTable.PrintHeader();
      NonDimTable << "Total Time" << config->GetMax_Time() << config->GetTime_Ref() << "s" << config->GetMax_Time()/config->GetTime_Ref();
      Unit.str("");
      NonDimTable << "Time Step" << config->GetTime_Step() << config->GetTime_Ref() << "s" << config->GetDelta_UnstTimeND();
      Unit.str("");
      NonDimTable.PrintFooter();
    }

    cout << ModelTableOut.str();
    cout << NonDimTableOut.str();

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
                   dest, dest, SU2_MPI::GetComm(), &sendReqs[i]);
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
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for the receive buffer, receive the message
       and determine the actual index of this rank in ranksRecv. */
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank, SU2_MPI::GetComm(), &status);

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
  SU2_MPI::Barrier(SU2_MPI::GetComm());

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
                   dest, dest+1, SU2_MPI::GetComm(), &invSendReqs[i]);
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
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for the receive buffer, receive the message
       and determine the actual index of this rank in ranksSend. */
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank+1, SU2_MPI::GetComm(), &status);

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
  SU2_MPI::Barrier(SU2_MPI::GetComm());

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
        int rankDOF = (int)(low - nDOFsPerRank.begin());
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
                   dest, dest+2, SU2_MPI::GetComm(), &sendReqs[i]);
  }

  /* Loop over the ranks from which I receive data to be processed. The number
     of ranks is equal to the number of ranks to which I just sent data. */
  vector<SU2_MPI::Request> sendReturnReqs(rankCommToInd.size());
  vector<vector<int> > sendReturnBuf(rankCommToInd.size(), vector<int>(0));
  for(unsigned long i=0; i<rankCommToInd.size(); ++i) {

    /* Block until a message arrives and determine the source and size
       of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+2, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for the receive buffer as well as for the
       return send buffer. Receive the message afterwards. */
    vector<unsigned long> recvBuf(sizeMess);
    sendReturnBuf[i].resize(sizeMess);

    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank+2, SU2_MPI::GetComm(), &status);

    /* Loop over the data just received and fill the return send buffer
       with the color of the DOFs. */
    for(int j=0; j<sizeMess; ++j) {
      const unsigned long jj = recvBuf[j] - nDOFsPerRank[rank];
      if(jj >= nDOFsLocOwned) {
        SU2_MPI::Error("This DOF should be owned, but it is not. This should not happen.",CURRENT_FUNCTION);
      }
      sendReturnBuf[i][j] = colorLocalDOFs[jj];
    }

    /* Send the return buffer back to the calling rank. Again use non-blocking
       sends to avoid deadlock. */
    SU2_MPI::Isend(sendReturnBuf[i].data(), sendReturnBuf[i].size(), MPI_INT,
                   source, source+3, SU2_MPI::GetComm(), &sendReturnReqs[i]);
  }

  /* Complete the first round of non-blocking sends. */
  SU2_MPI::Waitall(sendReqs.size(), sendReqs.data(), MPI_STATUSES_IGNORE);

  /* Loop over the ranks from which I receive the return information. */
  for(unsigned long i=0; i<rankCommToInd.size(); ++i) {

    /* Block until a message arrives and determine the source of the message
       and its index in the original send buffers. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+3, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    MI = rankCommToInd.find(source);
    const int ind = MI->second;

    /* Allocate the memory for the receive buffer and receive the message using
       a blocking receive. */
    vector<int> recvBuf(sendBuf[ind].size());
    SU2_MPI::Recv(recvBuf.data(), recvBuf.size(), MPI_INT,
                  source, rank+3, SU2_MPI::GetComm(), &status);

    /* Loop over the data just received and add them to the map
       mapMatrixIndToColor .*/
    for(unsigned long j=0; j<sendBuf[ind].size(); ++j)
      mapMatrixIndToColor[sendBuf[ind][j]] = recvBuf[j];
  }

  /* Complete the second round of non-blocking sends. */
  SU2_MPI::Waitall(sendReturnReqs.size(), sendReturnReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(SU2_MPI::GetComm());

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
     When only a spatial Jacobian is computed this is false per definition.  */
  if( (config->GetKind_TimeIntScheme_Flow() == ADER_DG) &&
     !(config->GetJacobian_Spatial_Discretization_Only()) ) {

    /*------------------------------------------------------------------------*/
    /* ADER time integration with local time stepping. There are 2^(M-1)      */
    /* subtime steps, where M is the number of time levels. For each of       */
    /* the subtime steps a number of tasks must be carried out, hence the     */
    /* total list can become rather lengthy.                                  */
    /*------------------------------------------------------------------------*/

    /*--- Determine whether or not the predictor solution must be interpolated
          for the time levels on this rank. Make a distinction between owned
          and halo elements. Also determine whether or not boundary conditions
          are present and whether or not these depend on halo data. ---*/
    const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
    vector<bool> interpolOwnedElem(nTimeLevels, false);
    vector<bool> interpolHaloElem(nTimeLevels, false);
    vector<bool> BCPresent(nTimeLevels, false);
    vector<bool> BCDependOnHalos(nTimeLevels, false);

    for(unsigned short level=0; level<nTimeLevels; ++level) {

      /* First check the boundary conditions. */
      for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); iMarker++) {

        const unsigned long surfElemBeg = boundaries[iMarker].nSurfElem[level];
        const unsigned long surfElemEnd = boundaries[iMarker].nSurfElem[level+1];
        if(surfElemEnd > surfElemBeg) {
          BCPresent[level] = true;
          if( boundaries[iMarker].haloInfoNeededForBC ) BCDependOnHalos[level] = true;
        }
      }

      /* Determine whether or not the owned elements must be interpolated. */
      if(nVolElemOwnedPerTimeLevel[level+1] > nVolElemOwnedPerTimeLevel[level])
        interpolOwnedElem[level] = true;
      if(nMatchingInternalFacesLocalElem[level+1] > nMatchingInternalFacesLocalElem[level])
       interpolOwnedElem[level] = true;
      if(nMatchingInternalFacesWithHaloElem[level+1] > nMatchingInternalFacesWithHaloElem[level])
        interpolOwnedElem[level] = true;
      if( BCPresent[level] )
        interpolOwnedElem[level] = true;

      /* Determine whether or not the halo elements must be interpolated. */
      if(nMatchingInternalFacesWithHaloElem[level+1] > nMatchingInternalFacesWithHaloElem[level])
        interpolHaloElem[level] = true;
      if( BCDependOnHalos[level] )
        interpolHaloElem[level] = true;
    }

    /* Determine the number of subtime steps and abbreviate the number of
       time integration points a bit easier.  */
    const unsigned int   nSubTimeSteps          = pow(2,nTimeLevels-1);
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

      /* Definition of the variable to store previous indices of
         tasks that must have been completed. */
      int prevInd[5];

      /* Carry out the predictor step of the communication elements of level 0
         if these elements are present on this rank. */
      unsigned long elemBeg = nVolElemOwnedPerTimeLevel[0]
                            + nVolElemInternalPerTimeLevel[0];
      unsigned long elemEnd = nVolElemOwnedPerTimeLevel[1];
      if(elemEnd > elemBeg) {
        prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][0];
        indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0] = (int)tasksList.size();
        tasksList.emplace_back(CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS, 0, prevInd[0]);
      }

      /* Initiate the communication of elements of level 0, if there is
         something to be communicated. */
#ifdef HAVE_MPI
      if( !commRequests[0].empty() ) {
        if(elemEnd > elemBeg)   // Data needs to be computed before sending.
          prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0];
        else                    // Data only needs to be received.
          prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][0];

        /* Make sure that any previous communication has been completed. */
        prevInd[1] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][0];

        /* Create the task. */
        indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][0] = tasksList.size();
        tasksList.emplace_back(CTaskDefinition::INITIATE_MPI_COMMUNICATION, 0,
                                            prevInd[0], prevInd[1]);
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
          indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL] = (int)tasksList.size();
          tasksList.emplace_back(CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS, nL, prevInd[0]);
        }

        /* Initiate the communication of elements of level nL, if there is
           something to be communicated. */
#ifdef HAVE_MPI
        if( !commRequests[nL].empty() ) {
          if(elemEnd > elemBeg)   // Data needs to be computed before sending.
            prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL];
          else                    // Data only needs to be received.
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][nL];

          /* Make sure that any previous communication has been completed. */
          prevInd[1] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][nL];

          /* Create the actual task. */
          indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][nL] = tasksList.size();
          tasksList.emplace_back(CTaskDefinition::INITIATE_MPI_COMMUNICATION, nL,
                                              prevInd[0], prevInd[1]);
        }
#endif
      }

      /* Carry out the predictor step of the internal elements of time level 0,
         if these elements are present on this rank. */
      elemBeg = nVolElemOwnedPerTimeLevel[0];
      elemEnd = nVolElemOwnedPerTimeLevel[0] + nVolElemInternalPerTimeLevel[0];

      if(elemEnd > elemBeg) {
        prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][0];
        indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][0] = (int)tasksList.size();
        tasksList.emplace_back(CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS, 0, prevInd[0]);
      }

      /* Determine the tasks to be completed before the communication of time
         level 0 can be completed. */
      prevInd[0] = prevInd[1] = prevInd[2] = -1;
#ifdef HAVE_MPI
      if( !commRequests[0].empty() )
        prevInd[0] = indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][0];
#endif

      if( !elementsSendSelfComm[0].empty() ) {
        prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0];
        prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][0];
      }

      /* Make sure that the -1 in prevInd, if any, are numbered last. */
      sort(prevInd, prevInd+3, greater<int>());

      /* Complete the communication of time level 0, if there is something
         to be completed. */
      if(prevInd[0] > -1) {
        indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][0] = (int)tasksList.size();
        tasksList.emplace_back(CTaskDefinition::COMPLETE_MPI_COMMUNICATION, 0,
                                            prevInd[0], prevInd[1], prevInd[2]);
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
          indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][nL] = (int)tasksList.size();
          tasksList.emplace_back(CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS, nL, prevInd[0]);
        }

        /* Determine the tasks to be completed before the communication of time
           level nL can be completed. */
        prevInd[0] = prevInd[1] = prevInd[2] = -1;
#ifdef HAVE_MPI
        if( !commRequests[nL].empty() )
          prevInd[0] = indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][nL];
#endif

        if( !elementsSendSelfComm[nL].empty() ) {
          prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL];
          prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][nL];
        }

        /* Make sure that the -1 in prevInd, if any, are numbered last. */
        sort(prevInd, prevInd+3, greater<int>());

        /* Complete the communication of time level nL, if there is something
           to be completed. */
        if(prevInd[0] > -1) {
          indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][nL] = (int)tasksList.size();
          tasksList.emplace_back(CTaskDefinition::COMPLETE_MPI_COMMUNICATION, nL,
                                              prevInd[0], prevInd[1], prevInd[2]);
        }
      }

      /* Loop over the time integration points to compute the space time
         integral for the corrector step. */
      for(unsigned short intPoint=0; intPoint<nTimeIntegrationPoints; ++intPoint) {

        /* Loop over the time levels to be treated. */
        for(unsigned short level=0; level<=timeLevel; ++level) {

          /* Easier storage of the number of owned elements for this level
             and the number of owned elements of the adjacent time level
             that are needed for the computation of the surface integral. */
          unsigned long nOwnedElem = nVolElemOwnedPerTimeLevel[level+1]
                                   - nVolElemOwnedPerTimeLevel[level];

          unsigned long nAdjOwnedElem = 0;
          if(level < (nTimeLevels-1))
            nAdjOwnedElem = ownedElemAdjLowTimeLevel[level+1].size();

          /* Check if the solution of the owned elements of the current time
             level and the adjacent owned elements of the next time level must
             be interpolated to the integration point intPoint. */
          if( interpolOwnedElem[level] ) {

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

              prevInd[4] = -1;
            }
            else {

              /* Not the first integration point. The residual computations of
                 the previous integration point must have been completed, because
                 the solution in the work vectors is overwritten. */
              prevInd[0] = indexInList[CTaskDefinition::VOLUME_RESIDUAL][level];
              prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level];
              prevInd[2] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
              prevInd[3] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
              prevInd[4] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level];
            }

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2], prevInd[3], prevInd[4]);

            /* The info on the integration point and whether or not this
               time integration corresponds to the second part for the
               adjacent elements must be added to the task just created. */
            tasksList.back().intPointADER          = intPoint;
            tasksList.back().secondPartTimeIntADER = level < timeLevel;

            /* If artificial viscosity is used for the shock capturing
               terms, these terms must be computed for the owned elements,
               including the adjacent ones. */
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS][level];
            indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS,
                                                level, prevInd[0]);
          }

          /* The solution of the halo elements of the current time level and
             the adjacent halo elements of the next time level must be
             interpolated to the integration point intPoint. Check if these
             elements are present. */
          if( interpolHaloElem[level] ) {

            unsigned long nAdjHaloElem = 0;
            if(level < (nTimeLevels-1))
              nAdjHaloElem = haloElemAdjLowTimeLevel[level+1].size();

            unsigned long nHaloElem = nVolElemHaloPerTimeLevel[level+1]
                                    - nVolElemHaloPerTimeLevel[level];

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
              prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
            }

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS,
                                                level, prevInd[0], prevInd[1]);

            /* The info on the integration point and whether or not this
               time integration corresponds to the second part for the
               adjacent elements must be added to the task just created. */
            tasksList.back().intPointADER          = intPoint;
            tasksList.back().secondPartTimeIntADER = level < timeLevel;

            /* If artificial viscosity is used for the shock capturing
               terms, these terms must be computed for the halo elements,
               including the adjacent ones. */
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][level];
            indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS,
                                                level, prevInd[0]);
          }

          /* Check whether there are boundary conditions that involve halo elements. */
          if( BCDependOnHalos[level] ) {

            /* Create the dependency list for this task. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            prevInd[1] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level];

            /* Create the task for the boundary conditions that involve halo elements. */
            indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO,
                                                level, prevInd[0], prevInd[1]);
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
            indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2]);

            /* Create the task to accumulate the surface residuals of the halo
               elements. Make sure to set the integration point for this task. */
            prevInd[0] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
            indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS,
                                                level, prevInd[0]);
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
            if( !commRequests[level].empty() ) commData = true;
#endif
            if(commData || !elementsSendSelfComm[level].empty()) {

              /* Determine the dependencies before the reverse communication
                 can start. */
              prevInd[0] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level];
              if( !haloElemAdjLowTimeLevel[level].empty() )
                prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level-1];
              else
                prevInd[1] = -1;

              prevInd[2] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level];

              /* Create the task. */
              indexInList[CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION][level] = (int)tasksList.size();
              tasksList.emplace_back(CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION,
                                                  level, prevInd[0], prevInd[1], prevInd[2]);
            }
          }

          /* Compute the contribution of the volume integral, boundary conditions that only
             involve owned elements and surface integral between owned elements for this
             time level, if present. */
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
            indexInList[CTaskDefinition::VOLUME_RESIDUAL][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::VOLUME_RESIDUAL, level,
                                                prevInd[0], prevInd[1]);

            indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED, level,
                                                prevInd[0], prevInd[1]);

            if(nMatchingInternalFacesLocalElem[level+1] >
               nMatchingInternalFacesLocalElem[level]) {
              indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level] = (int)tasksList.size();
              tasksList.emplace_back(CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS,
                                                  level, prevInd[0], prevInd[1]);
            }
          }

          /* Accumulate the space time residuals for the owned elements of this
             level, if these elements are present. */
          if(nAdjOwnedElem || nOwnedElem) {

            /* Create the dependencies for this task. */
            prevInd[0] = indexInList[CTaskDefinition::VOLUME_RESIDUAL][level];
            prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level];
            prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
            prevInd[3] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level];
            prevInd[4] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.emplace_back(CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2], prevInd[3], prevInd[4]);
            tasksList.back().intPointADER = intPoint;
          }

          /* If this is the last integration point, complete the reverse
             communication for the residuals of this time level, if there
             is something to communicate. */
          if(intPoint == (nTimeIntegrationPoints-1)) {

            bool commData = false;

#ifdef HAVE_MPI
            if( !commRequests[level].empty() ) commData = true;
#endif
            if(commData || !elementsSendSelfComm[level].empty()) {

              /* Determine the dependencies before the reverse communication
                 can be completed. */
              prevInd[0] = indexInList[CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION][level];
              prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];
              if( !ownedElemAdjLowTimeLevel[level].empty() )
                prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level-1];
              else
                prevInd[2] = -1;

              /* Create the task. */
              indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][level] = (int)tasksList.size();
              tasksList.emplace_back(CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION,
                                                  level, prevInd[0], prevInd[1], prevInd[2]);
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
          if( !ownedElemAdjLowTimeLevel[level].empty() )
            prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level-1];
          else
            prevInd[2] = -1;

          indexInList[CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX][level] = (int)tasksList.size();
          tasksList.emplace_back(CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX,
                                               level, prevInd[0], prevInd[1], prevInd[2]);

          /* Compute the new state vector for this time level. */
          prevInd[0] = indexInList[CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX][level];
          indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][level] = (int)tasksList.size();
          tasksList.emplace_back(CTaskDefinition::ADER_UPDATE_SOLUTION,
                                               level, prevInd[0]);
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
            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED: cout << "BOUNDARY_CONDITIONS_DEPEND_ON_OWNED" << endl; break;
            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO: cout << "BOUNDARY_CONDITIONS_DEPEND_ON_HALO" << endl; break;
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
      SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
    }

    cout << "CFEM_DG_EulerSolver::SetUpTaskList: ADER tasklist printed" << endl;
    exit(1); */

    // END EXTRA FOR DEBUGGING.
  }
  else {

    /*------------------------------------------------------------------------*/
    /* Standard time integration scheme for which the spatial residual must   */
    /* be computed for the DOFS of the owned elements. This results in a      */
    /* relatively short tasks list, which can be set easily.                  */
    /*------------------------------------------------------------------------*/

    tasksList.emplace_back(CTaskDefinition::INITIATE_MPI_COMMUNICATION,                   0);                   // Task  0
    tasksList.emplace_back(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS,     0);                   // Task  1
    tasksList.emplace_back(CTaskDefinition::VOLUME_RESIDUAL,                              0,  1);               // Task  2
    tasksList.emplace_back(CTaskDefinition::COMPLETE_MPI_COMMUNICATION,                   0,  0);               // Task  3
    tasksList.emplace_back(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS,      0,  3);               // Task  4
    tasksList.emplace_back(CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS,               0,  1, 4);            // Task  5
    tasksList.emplace_back(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO,           0,  1, 3);            // Task  6
    tasksList.emplace_back(CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS,  0,  5);               // Task  7
    tasksList.emplace_back(CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION,           0,  7);               // Task  8
    tasksList.emplace_back(CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS,              0,  1);               // Task  9
    tasksList.emplace_back(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED,          0,  1);               // Task 10
    tasksList.emplace_back(CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION,           0,  2, 8);            // Task 11
    tasksList.emplace_back(CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS, 0,  2, 6, 9, 10, 11); // Task 12
    tasksList.emplace_back(CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX,                 0, 12);               // Task 13
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
      if((ranksRecv[i] != rank) && !elemRecvPerTimeLevel[level][i].empty()) ++nRankRecv;
    }

    /* Determine the number of ranks to which data will be send
       for this time level. Self communication is excluded. */
    int nRankSend = 0;
    for(unsigned long i=0; i<ranksSend.size(); ++i) {
      if((ranksSend[i] != rank) && !elemSendPerTimeLevel[level][i].empty()) ++nRankSend;
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
      if((ranksRecv[i] != rank) && !elemRecvPerTimeLevel[level][i].empty()) {

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
      if((ranksSend[i] != rank) && !elemSendPerTimeLevel[level][i].empty()) {

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
  if( !commRequests[timeLevel].empty() ) {

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

        for(unsigned short k=0; k<nTimeDOFs; ++k) {
          const unsigned long indS = nVar*(volElem[jj].offsetDOFsSolLocal + k*nDOFsLocTot);
          for(unsigned long mm=0; mm<nItems; ++mm)
            sendBuf[ii+mm] = commData[indS+mm];
          ii += nItems;
        }
      }

      /* Send the data using non-blocking sends. */
      int dest = ranksSendMPI[timeLevel][i];
      int tag  = dest + timeLevel;
      SU2_MPI::Isend(sendBuf, ii, MPI_DOUBLE, dest, tag, SU2_MPI::GetComm(),
                     &commRequests[timeLevel][indComm]);
    }

    /* Loop over the number of ranks from which data is received. */
    for(unsigned long i=0; i<ranksRecvMPI[timeLevel].size(); ++i, ++indComm) {

      /* Post the non-blocking receive. */
      int source = ranksRecvMPI[timeLevel][i];
      int tag    = rank + timeLevel;
      SU2_MPI::Irecv(commRecvBuf[timeLevel][i].data(),
                     commRecvBuf[timeLevel][i].size(),
                     MPI_DOUBLE, source, tag, SU2_MPI::GetComm(),
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
  if( !commRequests[timeLevel].empty() ) {

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

        for(unsigned short k=0; k<nTimeDOFs; ++k) {
          const unsigned long indR = nVar*(volElem[jj].offsetDOFsSolLocal + k*nDOFsLocTot);
          for(unsigned long mm=0; mm<nItems; ++mm)
            commData[indR+mm] = recvBuf[ii+mm];
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
    const unsigned long nItems = volElem[elemS].nDOFsSol * nVar;

    for(unsigned short j=0; j<nTimeDOFs; ++j) {
      const unsigned long indS = nVar*(volElem[elemS].offsetDOFsSolLocal + j*nDOFsLocTot);
      const unsigned long indR = nVar*(volElem[elemR].offsetDOFsSolLocal + j*nDOFsLocTot);

      for(unsigned long mm=0; mm<nItems; ++mm)
        commData[indR+mm] = commData[indS+mm];
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
  if( !commRequests[timeLevel].empty() ) {

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

        const unsigned long indS = nVar*volElem[jj].offsetDOFsSolLocal;
        for(unsigned long mm=0; mm<nItems; ++mm)
          recvBuf[ii+mm] = resComm[indS+mm];
        ii += nItems;
      }

      /* Send the data using non-blocking sends. */
      int dest = ranksRecvMPI[timeLevel][i];
      int tag  = dest + timeLevel + 20;
      SU2_MPI::Isend(recvBuf, ii, MPI_DOUBLE, dest, tag, SU2_MPI::GetComm(),
                     &commRequests[timeLevel][indComm]);
    }

    /* Post the non-blocking receives. As this is the reverse communication,
       a loop over the sending ranks must be carried out. */
    for(unsigned long i=0; i<ranksSendMPI[timeLevel].size(); ++i, ++indComm) {

      int source = ranksSendMPI[timeLevel][i];
      int tag    = rank + timeLevel + 20;
      SU2_MPI::Irecv(commSendBuf[timeLevel][i].data(),
                     commSendBuf[timeLevel][i].size(),
                     MPI_DOUBLE, source, tag, SU2_MPI::GetComm(),
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
  if( !commRequests[timeLevel].empty() ) {

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

void CFEM_DG_EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  /*--- Check if a verification solution is to be computed. ---*/
  if ((VerificationSolution)  && (TimeIter == 0)) {

    /* Loop over the owned elements. */
    for(unsigned long i=0; i<nVolElemOwned; ++i) {

      /* Loop over the DOFs of this element. */
      for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {

        /* Set the pointers to the coordinates and solution of this DOF. */
        const su2double *coor = volElem[i].coorSolDOFs.data() + j*nDim;
        su2double *solDOF     = VecSolDOFs.data() + nVar*(volElem[i].offsetDOFsSolLocal + j);

        /* Set the solution in this DOF to the initial condition provided by
           the verification solution class. This can be the exact solution,
           but this is not necessary. */
        VerificationSolution->GetInitialCondition(coor, solDOF);
      }
    }
  }
}

void CFEM_DG_EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long ErrorCounter = 0;

  /*-----------------------------------------------------------------------------*/
  /*--- Check for non-physical points. Only needed for a compressible solver. ---*/
  /*-----------------------------------------------------------------------------*/

  if(config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {

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
    if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
      unsigned long MyErrorCounter = ErrorCounter;
      SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#endif
      if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
    }
  }

  /*-----------------------------------------------------------------------------*/
  /*                       Check for grid motion.                                */
  /*-----------------------------------------------------------------------------*/

  const bool harmonic_balance = config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE;
  if(config->GetGrid_Movement() && !harmonic_balance) {

    /*--- Determine the type of grid motion. ---*/
    switch( config->GetKind_GridMovement() ) {

      case ROTATING_FRAME:
      case STEADY_TRANSLATION: {

        /* Steady motion. This is accounted for in the general preprocessing,
           where grid velocities are computed. */
        break;
      }

      case RIGID_MOTION: {

        /* Rigid body motion described. At the moment this is only
           possible for the Classical Runge Kutta scheme. */
        if(config->GetKind_TimeIntScheme() != CLASSICAL_RK4_EXPLICIT)
          SU2_MPI::Error("Rigid body motion only possible for CLASSICAL_RK4_EXPLICIT.",
                         CURRENT_FUNCTION);

        /* Determine whether or not it is needed to compute the motion data. */
        const unsigned long TimeIter = config->GetTimeIter();

        bool computeMotion = false, firstTime = false;
        if(TimeIter == 0 && iStep == 0) computeMotion = firstTime = true;
        if(iStep == 1 || iStep == 3)   computeMotion = true;

        if( computeMotion ) {

          /* Determine the number of time levels. For the classical Runge-Kutta
             scheme this should be 1. */
          const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

          /* Determine the time for which the motion data must be determined. */
          const su2double deltaT = config->GetDelta_UnstTimeND();

          su2double tNew = TimeIter*deltaT;
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
  Set_OldSolution();

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

#ifdef CODI_FORWARD_TYPE
          /* Set the pointer to the residual of the current DOF. */
          const su2double *resDOF = VecResDOFs.data() + i*nVar;
#endif

          /* Store the matrix entries. */
          for(unsigned short j=0; j<nVar; ++j) {
#ifdef CODI_FORWARD_TYPE
            Jac[var+j*nVar] = resDOF[j].getGradient();
#else
            Jac[var+j*nVar] = 0.0;   /* This is to avoid a compiler warning. */
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

void CFEM_DG_EulerSolver::Set_OldSolution() {

  for(unsigned long i=0; i<VecSolDOFs.size(); ++i)
    VecWorkSolDOFs[0][i] = VecSolDOFs[i];
}

void CFEM_DG_EulerSolver::Set_NewSolution() {

  for(unsigned long i=0; i<VecSolDOFs.size(); ++i)
    VecSolDOFsNew[i] = VecSolDOFs[i];
}

void CFEM_DG_EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /* Check whether or not a time stepping scheme is used. */
  const bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;

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
    if(config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {

      /*--- Loop over the owned volume elements. ---*/
      for(unsigned long l=0; l<nVolElemOwned; ++l) {

        su2double charVel2Max = 0.0;

        /*--- Make a distinction between 2D and 3D for optimal performance. ---*/
        switch( nDim ) {

          case 2: {
            /*--- 2D simulation. Loop over the DOFs of this element
                  and determine the maximum wave speed. ---*/
            for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {
              const su2double *solDOF  = VecSolDOFs.data()
                                       + nVar*(volElem[l].offsetDOFsSolLocal + i);
              const su2double *gridVel = volElem[l].gridVelocitiesSolDOFs.data()
                                       + i*nDim;

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
            for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {
              const su2double *solDOF = VecSolDOFs.data()
                                      + nVar*(volElem[l].offsetDOFsSolLocal + i);
              const su2double *gridVel = volElem[l].gridVelocitiesSolDOFs.data()
                                       + i*nDim;

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
        const unsigned short ind = volElem[l].indStandardElement;
        unsigned short nPoly = standardElementsSol[ind].GetNPoly();
        if(nPoly == 0) nPoly = 1;

        const su2double lenScaleInv = nPoly/volElem[l].lenScale;
        const su2double dtInv       = lenScaleInv*sqrt(charVel2Max);

        VecDeltaTime[l] = CFL/dtInv;

        const su2double dtEff = volElem[l].factTimeLevel*VecDeltaTime[l];
        Min_Delta_Time = min(Min_Delta_Time, dtEff);
        Max_Delta_Time = max(Max_Delta_Time, dtEff);
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
    if ((config->GetComm_Level() == COMM_FULL) || time_stepping) {
#ifdef HAVE_MPI
      su2double rbuf_time = Min_Delta_Time;
      SU2_MPI::Allreduce(&rbuf_time, &Min_Delta_Time, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());

      rbuf_time = Max_Delta_Time;
      SU2_MPI::Allreduce(&rbuf_time, &Max_Delta_Time, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
#endif
    }

    /*--- For explicit time stepping with an unsteady CFL imposed, use the
          minimum delta time of the entire mesh. As Min_Delta_Time is scaled to
          the time step of the largest time level, a correction must be used
          for the time level when time accurate local time stepping is used. ---*/
    if (time_stepping) {
      for(unsigned long l=0; l<nVolElemOwned; ++l)
        VecDeltaTime[l] = Min_Delta_Time/volElem[l].factTimeLevel;

      config->SetDelta_UnstTimeND(Min_Delta_Time);
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
  /* Easier storage of the number of time levels.. */
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

  /* Define and initialize the bool vector, that indicates whether or
     not the tasks from the list have been completed. */
  vector<bool> taskCompleted(tasksList.size(), false);

  /* Allocate the memory for the work array and initialize it to zero to avoid
     warnings in debug mode  about uninitialized memory when padding is applied. */
  vector<su2double> workArrayVec(sizeWorkArray, 0.0);
  su2double *workArray = workArrayVec.data();

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

              ADER_DG_PredictorStep(config, elemBeg, elemEnd, workArray);
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
              ADER_DG_PredictorStep(config, elemBeg, elemEnd, workArray);
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
              unsigned long nAdjElem = 0, *adjElem = nullptr;
              if(level < (nTimeLevels-1)) {
                nAdjElem = ownedElemAdjLowTimeLevel[level+1].size();
                adjElem  = ownedElemAdjLowTimeLevel[level+1].data();
              }

              ADER_DG_TimeInterpolatePredictorSol(config, tasksList[i].intPointADER,
                                                  nVolElemOwnedPerTimeLevel[level],
                                                  nVolElemOwnedPerTimeLevel[level+1],
                                                  nAdjElem, adjElem,
                                                  tasksList[i].secondPartTimeIntADER,
                                                  VecWorkSolDOFs[level].data());
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS: {

              /* Interpolate the predictor solution of the halo elements
                 in time to the given time integration point for the
                 given time level. */
              const unsigned short level = tasksList[i].timeLevel;
              unsigned long nAdjElem = 0, *adjElem = nullptr;
              if(level < (nTimeLevels-1)) {
                nAdjElem = haloElemAdjLowTimeLevel[level+1].size();
                adjElem  = haloElemAdjLowTimeLevel[level+1].data();
              }

              ADER_DG_TimeInterpolatePredictorSol(config, tasksList[i].intPointADER,
                                                  nVolElemHaloPerTimeLevel[level],
                                                  nVolElemHaloPerTimeLevel[level+1],
                                                  nAdjElem, adjElem,
                                                  tasksList[i].secondPartTimeIntADER,
                                                  VecWorkSolDOFs[level].data());
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS: {

              /*--- Compute the artificial viscosity for shock capturing in DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              Shock_Capturing_DG(config, nVolElemOwnedPerTimeLevel[level],
                                 nVolElemOwnedPerTimeLevel[level+1], workArray);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS: {

              /*--- Compute the artificial viscosity for shock capturing in DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              Shock_Capturing_DG(config, nVolElemHaloPerTimeLevel[level],
                                 nVolElemHaloPerTimeLevel[level+1], workArray);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::VOLUME_RESIDUAL: {

              /*--- Compute the volume portion of the residual. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              Volume_Residual(config, nVolElemOwnedPerTimeLevel[level],
                              nVolElemOwnedPerTimeLevel[level+1], workArray);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS: {

              /* Compute the residual of the faces that only involve owned elements. */
              const unsigned short level = tasksList[i].timeLevel;
              unsigned long indResFaces = startLocResInternalFacesLocalElem[level];
              ResidualFaces(config, nMatchingInternalFacesLocalElem[level],
                            nMatchingInternalFacesLocalElem[level+1],
                            indResFaces, numerics[CONV_TERM], workArray);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS: {

              /* Compute the residual of the faces that involve a halo element. */
              const unsigned short level = tasksList[i].timeLevel;
              unsigned long indResFaces = startLocResInternalFacesWithHaloElem[level];
              ResidualFaces(config, nMatchingInternalFacesWithHaloElem[level],
                            nMatchingInternalFacesWithHaloElem[level+1],
                            indResFaces, numerics[CONV_TERM], workArray);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED: {

              /*--- Apply the boundary conditions that only depend on data
                    of owned elements. ---*/
              Boundary_Conditions(tasksList[i].timeLevel, config, numerics, false,
                                  workArray);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO: {

              /*--- Apply the boundary conditions that also depend on data
                    of halo elements. ---*/
              Boundary_Conditions(tasksList[i].timeLevel, config, numerics, true,
                                  workArray);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS: {

              /* Create the final residual by summing up all contributions. */
              CreateFinalResidual(tasksList[i].timeLevel, true);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS: {

              /* Create the final residual by summing up all contributions. */
              CreateFinalResidual(tasksList[i].timeLevel, false);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS: {

              /* Accumulate the space time residuals for the owned elements
                 for ADER-DG. */
              AccumulateSpaceTimeResidualADEROwnedElem(config, tasksList[i].timeLevel,
                                                       tasksList[i].intPointADER);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS: {

              /* Accumulate the space time residuals for the halo elements
                 for ADER-DG. */
              AccumulateSpaceTimeResidualADERHaloElem(config, tasksList[i].timeLevel,
                                                      tasksList[i].intPointADER);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX: {

              /*--- Multiply the residual by the (lumped) mass matrix, to obtain the final value. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              const bool useADER = config->GetKind_TimeIntScheme() == ADER_DG;
              MultiplyResidualByInverseMassMatrix(config, useADER,
                                                  nVolElemOwnedPerTimeLevel[level],
                                                  nVolElemOwnedPerTimeLevel[level+1],
                                                  workArray);
              taskCarriedOut = taskCompleted[i] = true;
              break;
            }

            case CTaskDefinition::ADER_UPDATE_SOLUTION: {

              /*--- Perform the update step for ADER-DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              ADER_DG_Iteration(nVolElemOwnedPerTimeLevel[level],
                                nVolElemOwnedPerTimeLevel[level+1]);
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
  /* Preprocessing. */
  Preprocessing(geometry, solver_container, config, iMesh, 0, RunTime_EqSystem, false);
  TolerancesADERPredictorStep();

  /* Process the tasks list to carry out one ADER space time integration step. */
  ProcessTaskList_DG(geometry, solver_container, numerics, config, iMesh);

  /* Postprocessing. */
  Postprocessing(geometry, solver_container, config, iMesh);
}

void CFEM_DG_EulerSolver::TolerancesADERPredictorStep() {

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
                     SU2_MPI::GetComm());
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
                                                const unsigned long elemEnd,
                                                su2double           *workArray) {

  /* Get the data of the ADER time integration scheme. */
  const unsigned short nTimeDOFs              = config->GetnTimeDOFsADER_DG();
  const unsigned short nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();
  const su2double     *timeIntegrationWeights = config->GetWeightsIntegrationADER_DG();
  const bool          useAliasedPredictor     = config->GetKind_ADER_Predictor() == ADER_ALIASED_PREDICTOR;

   /* Determine the number of solution entities that are treated simultaneously
      in the matrix products to obtain good gemm performance. A solution entity
      can be a time integration point (to compute the fluxes) or a time DOF
      in the actual iteration process. */
  const unsigned short NPad       = config->GetSizeMatMulPadding();
  const unsigned short nTimeSimul = NPad/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short NPadMin = 64/sizeof(passivedouble);

  /* Determine the number of time integration points that are treated
     simultaneously for better gemm performance as well as the number of
     loops to treat all time integration points. */
  unsigned short nTimeSimulLastInts = nTimeIntegrationPoints%nTimeSimul;
  unsigned short nTimeIntLoop       = nTimeIntegrationPoints/nTimeSimul;
  if(nTimeSimulLastInts == 0) nTimeSimulLastInts = nTimeSimul;
  else                        ++nTimeIntLoop;

  unsigned short NPadLastInts = nTimeSimulLastInts*nVar;
  if( NPadLastInts%NPadMin ) NPadLastInts += NPadMin - (NPadLastInts%NPadMin);

  /* Idem for the time DOFs. */
  unsigned short nTimeSimulLastDOFs = nTimeDOFs%nTimeSimul;
  unsigned short nTimeDOFLoop       = nTimeDOFs/nTimeSimul;
  if(nTimeSimulLastDOFs == 0) nTimeSimulLastDOFs = nTimeSimul;
  else                        ++nTimeDOFLoop;

  unsigned short NPadLastDOFs = nTimeSimulLastDOFs*nVar;
  if( NPadLastDOFs%NPadMin ) NPadLastDOFs += NPadMin - (NPadLastDOFs%NPadMin);

  /* Determine the total padded size to store the residual for all the
     spatial and time DOFs. The residual is stored in such a way that a
     memcpy is avoided when the multiplication with iteration matrix, i.e.
     mass matrix, is carried out. */
  const unsigned short NPadResTot = (nTimeDOFLoop-1)*NPad + NPadLastDOFs;

  /*--------------------------------------------------------------------------*/
  /*--- Loop over the given elemen range to compute the predictor solution.---*/
  /*--- For the predictor solution only an integration over the element is ---*/
  /*--- performed to obtain the weak formulation, i.e. no integration by   ---*/
  /*--- parts. As a consequence there is no surface term and also no       ---*/
  /*--- coupling with neighboring elements. This is also the reason why    ---*/
  /*--- the ADER scheme is only conditionally stable.                      ---*/
  /*--------------------------------------------------------------------------*/

  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* Easier storage of the number of spatial DOFs and the total
       number of variables for this element per time level. */
    const unsigned short nDOFs     = volElem[l].nDOFsSol;
    const unsigned short nVarNDOFs = nVar*nDOFs;

    /* Set the pointers for the working variables. */
    su2double *resInt  = workArray;
    su2double *solPred = resInt  + NPad*nDOFs;
    su2double *dSol    = solPred + nVarNDOFs*nTimeDOFs;
    su2double *resTot  = dSol    + nVarNDOFs*nTimeDOFs;
    su2double *solInt  = resTot  + NPadResTot*nDOFs;
    su2double *work    = solInt  + NPad*nDOFs;

    /* Initialize the predictor solution to the current solution. */
    const su2double *solCur = VecSolDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    for(unsigned short j=0; j<nTimeDOFs; ++j) {
      su2double *solPredTimeInd = solPred + j*nVarNDOFs;
      for(unsigned short mm=0; mm<nVarNDOFs; ++mm)
        solPredTimeInd[mm] = solCur[mm];
    }

    /*-------------------------------------------------------------------------*/
    /*--- Iterative algorithm to compute the predictor solution for all     ---*/
    /*--- the time DOFs of this element simultaneously.                     ---*/
    /*-------------------------------------------------------------------------*/

    for(;;) {

      /* Initialize the total residual to zero. */
      for(unsigned short i=0; i<(NPadResTot*nDOFs); ++i) resTot[i] = 0.0;

      /* Loop, such that all time integration points are treated. Note that
         inside this loop, several time integration points are treated
         simultaneously. */
      for(unsigned short timeIntLoop=0; timeIntLoop<nTimeIntLoop; ++timeIntLoop) {

        /* Determine the starting and ending time integration point, as well as
           the value of the padding to be applied. */
        const unsigned short intStart = timeIntLoop*nTimeSimul;
        const unsigned short intEnd   = (timeIntLoop == (nTimeIntLoop-1)) ?
                                        intStart + nTimeSimulLastInts : intStart + nTimeSimul;

        const unsigned short NPadInt = (timeIntLoop == (nTimeIntLoop-1)) ? NPadLastInts : NPad;

        /*--------------------------------------------------------------------*/
        /*--- Interpolate the predictor solution to the time integration   ---*/
        /*--- points to be computed. It is possible, or even likely, that  ---*/
        /*--- these integration points coincide with the locations of the  ---*/
        /*--- time DOFs. When this is the case the interpolation boils     ---*/
        /*--- down to a copy. Note that the solution of the integration    ---*/
        /*--- points are stored such that a simultaneous treatment of the  ---*/
        /*--- integration points is possible.                              ---*/
        /*--------------------------------------------------------------------*/

        /* Initialize the interpolated solution to zero. */
        for(unsigned short i=0; i<(NPadInt*nDOFs); ++i) solInt[i] = 0.0;

        for(unsigned short intPoint=intStart; intPoint<intEnd; ++intPoint) {

          /* Store the interpolation data for this time integration point. */
          const su2double *DOFToThisTimeInt = timeInterpolDOFToIntegrationADER_DG
                                            + intPoint*nTimeDOFs;

          /* Set the pointer to the address of the 1st variable of the 1st DOF
             for this integration point. In this way there is no need to
             account for this offset in the actual interpolation loop. */
          su2double *solThisInt = solInt + nVar*(intPoint - intStart);

          /* Carry out the actual interpolation. Make a distinction between
             2D and 3D for performance reasons. */
          switch( nDim ) {
            case 2: {
              for(unsigned short j=0; j<nTimeDOFs; ++j) {

                if(DOFToThisTimeInt[j] != 0.0) {
                  const unsigned int indPredTime = j*nVarNDOFs;

                  for(unsigned short i=0; i<nDOFs; ++i) {
                    const unsigned int indSol  = i*NPadInt;
                    const unsigned int indPred = indPredTime + i*nVar;

                    solThisInt[indSol]   += DOFToThisTimeInt[j]*solPred[indPred];
                    solThisInt[indSol+1] += DOFToThisTimeInt[j]*solPred[indPred+1];
                    solThisInt[indSol+2] += DOFToThisTimeInt[j]*solPred[indPred+2];
                    solThisInt[indSol+3] += DOFToThisTimeInt[j]*solPred[indPred+3];
                  }
                }
              }

              break;
            }

            case 3: {
              for(unsigned short j=0; j<nTimeDOFs; ++j) {

                if(DOFToThisTimeInt[j] != 0.0) {
                  const unsigned int indPredTime = j*nVarNDOFs;

                  for(unsigned short i=0; i<nDOFs; ++i) {
                    const unsigned int indSol  = i*NPadInt;
                    const unsigned int indPred = indPredTime + i*nVar;

                    solThisInt[indSol]   += DOFToThisTimeInt[j]*solPred[indPred];
                    solThisInt[indSol+1] += DOFToThisTimeInt[j]*solPred[indPred+1];
                    solThisInt[indSol+2] += DOFToThisTimeInt[j]*solPred[indPred+2];
                    solThisInt[indSol+3] += DOFToThisTimeInt[j]*solPred[indPred+3];
                    solThisInt[indSol+4] += DOFToThisTimeInt[j]*solPred[indPred+4];
                  }
                }
              }

              break;
            }
          }
        }

        /*--------------------------------------------------------------------*/
        /*--- Compute the spatial residual of the predictor step for the   ---*/
        /*--- time integration points considered. A distinction is made    ---*/
        /*--- between an aliased and a non-aliased evaluation of the       ---*/
        /*--- predictor residual and between 2D and 3D. The latter is for  ---*/
        /*--- performance reasons.                                         ---*/
        /*--------------------------------------------------------------------*/

        const unsigned short nIntSimul = intEnd - intStart;
        switch( nDim ) {
          case 2: {
            if( useAliasedPredictor )
              ADER_DG_AliasedPredictorResidual_2D(config, &volElem[l], solInt,
                                                  nIntSimul, NPadInt, resInt, work);
            else
              ADER_DG_NonAliasedPredictorResidual_2D(config, &volElem[l], solInt,
                                                     nIntSimul, NPadInt, resInt, work);
            break;
          }

          case 3: {
            if( useAliasedPredictor )
              ADER_DG_AliasedPredictorResidual_3D(config, &volElem[l], solInt,
                                                  nIntSimul, NPadInt, resInt, work);
            else
              ADER_DG_NonAliasedPredictorResidual_3D(config, &volElem[l], solInt,
                                                     nIntSimul, NPadInt, resInt, work);
            break;
          }
        }

        /*--------------------------------------------------------------------*/
        /*--- Update the total residual with the residual of the time      ---*/
        /*--- integration points considered. Note the minus sign, because  ---*/
        /*--- the residual is put on the RHS of the equation. Also note    ---*/
        /*--- that the residuals are stored in such a way that a memcpy    ---*/
        /*--- is avoided when the multiplication takes place with the      ---*/
        /*--- iteration matrix (i.e. mass matrix) later on.                ---*/
        /*--------------------------------------------------------------------*/

        /* Loop such that all time DOFs are treated. Note that in this loop
           several time DOFs are treated simultaneously. */
        for(unsigned short timeDOFLoop=0; timeDOFLoop<nTimeDOFLoop; ++timeDOFLoop) {

          /* Set the pointer where the residual for this chunk of time DOFs starts. */
          su2double *res = resTot + timeDOFLoop*NPad*nDOFs;

          /* Determine the starting and ending time DOF, as well as
             the value of the padding to be applied. */
          const unsigned short DOFStart = timeDOFLoop*nTimeSimul;
          const unsigned short DOFEnd   = (timeDOFLoop == (nTimeDOFLoop-1)) ?
                                          DOFStart + nTimeSimulLastDOFs : DOFStart + nTimeSimul;

          const unsigned short NPadDOF = (timeDOFLoop == (nTimeDOFLoop-1)) ? NPadLastDOFs : NPad;

          /* Loop over the time integration points that are treated simultaneously. */
          for(unsigned short intPoint=intStart; intPoint<intEnd; ++intPoint) {

            /* Set the pointer to the address of the 1st residual of the 1st DOF
               for this integration point. In this way there is no need to
               account for this offset in the loop below. */
            const su2double *resThisInt = resInt + nVar*(intPoint - intStart);

            /* Store the interpolation data for this time integration point. */
            const su2double *DOFToThisTimeInt = timeInterpolDOFToIntegrationADER_DG
                                              + intPoint*nTimeDOFs;

            /* Update the actual interpolation. Make a distinction between
               2D and 3D for performance reasons. */
            switch( nDim ) {
              case 2: {

                /* Loop over the time DOFs for this time integration point and check
                   if this time integration point contributes to the residual of
                   this time DOF. */
                for(unsigned short iTimeDOF=DOFStart; iTimeDOF<DOFEnd; ++iTimeDOF) {
                  if(DOFToThisTimeInt[iTimeDOF] != 0.0) {

                    /* Determine the index where the 1st residual of the 1st spatial
                       DOF starts for this time DOF. Also determine the multiplication
                       factor for this time DOF. The factor 0.5 is present, because
                       the length of the interval in parameter space, [-1..1], is 2. */
                    const unsigned int indResThisTimeDOF = nVar*(iTimeDOF-DOFStart);
                    const su2double w = 0.5*timeIntegrationWeights[intPoint]
                                      * VecDeltaTime[l]*DOFToThisTimeInt[iTimeDOF];

                    /* Loop over the number of spatial DOFs to update their residuals. Note
                       the minus sign, see the comment in the larger comment section above. */
                    for(unsigned short i=0; i<nDOFs; ++i) {

                      const unsigned int indResThisInt = i*NPadInt;
                      const unsigned int indResDOF     = indResThisTimeDOF + i*NPadDOF;

                      res[indResDOF]   -= w*resThisInt[indResThisInt];
                      res[indResDOF+1] -= w*resThisInt[indResThisInt+1];
                      res[indResDOF+2] -= w*resThisInt[indResThisInt+2];
                      res[indResDOF+3] -= w*resThisInt[indResThisInt+3];
                    }
                  }
                }

                break;
              }

              case 3: {

                /* Loop over the time DOFs for this time integration point and check
                   if this time integration point contributes to the residual of
                   this time DOF. */
                for(unsigned short iTimeDOF=DOFStart; iTimeDOF<DOFEnd; ++iTimeDOF) {
                  if(DOFToThisTimeInt[iTimeDOF] != 0.0) {

                    /* Determine the index where the 1st residual of the 1st spatial
                       DOF starts for this time DOF. Also determine the multiplication
                       factor for this time DOF. The factor 0.5 is present, because
                       the length of the interval in parameter space, [-1..1], is 2. */
                    const unsigned int indResThisTimeDOF = nVar*(iTimeDOF-DOFStart);
                    const su2double w = 0.5*timeIntegrationWeights[intPoint]
                                      * VecDeltaTime[l]*DOFToThisTimeInt[iTimeDOF];

                    /* Loop over the number of spatial DOFs to update their residuals. Note
                       the minus sign, see the comment in the larger comment section above. */
                    for(unsigned short i=0; i<nDOFs; ++i) {

                      const unsigned int indResThisInt = i*NPadInt;
                      const unsigned int indResDOF     = indResThisTimeDOF + i*NPadDOF;

                      res[indResDOF]   -= w*resThisInt[indResThisInt];
                      res[indResDOF+1] -= w*resThisInt[indResThisInt+1];
                      res[indResDOF+2] -= w*resThisInt[indResThisInt+2];
                      res[indResDOF+3] -= w*resThisInt[indResThisInt+3];
                      res[indResDOF+4] -= w*resThisInt[indResThisInt+4];
                    }
                  }
                }

                break;
              }
            }
          }
        }
      }

      /*----------------------------------------------------------------------*/
      /*--- Determine the updates of the solution of the time DOFs         ---*/
      /*--- relative to the solution of the previous time step. This is    ---*/
      /*--- accomplished by multiplying the total residual of the spatial  ---*/
      /*--- DOFs in all time DOFs by the inverse of the mass matrix and an ---*/
      /*--- appropriate weight (timeCoefADER_DG). The residual of the time ---*/
      /*--- DOFs is stored such that multiple time DOFs can be treated     ---*/
      /*--- simultaneously without the need for a memcpy. In this way the  ---*/
      /*--- performance of the matrix multiplication is optimal.           ---*/
      /*----------------------------------------------------------------------*/

      /* Initialize the update to zero. */
      for(unsigned short i=0; i<(nVarNDOFs*nTimeDOFs); ++i) dSol[i] = 0.0;

      /* Loop such that all time DOFs are treated. Note that in this loop
         several time DOFs are treated simultaneously. */
      for(unsigned short timeDOFLoop=0; timeDOFLoop<nTimeDOFLoop; ++timeDOFLoop) {

        /* Set the pointer where the residual for this chunk of time DOFs starts. */
        su2double *res = resTot + timeDOFLoop*NPad*nDOFs;

        /* Determine the starting and ending time DOF, as well as
           the value of the padding to be applied. */
        const unsigned short DOFStart = timeDOFLoop*nTimeSimul;
        const unsigned short DOFEnd   = (timeDOFLoop == (nTimeDOFLoop-1)) ?
                                        DOFStart + nTimeSimulLastDOFs : DOFStart + nTimeSimul;

        const unsigned short NPadDOF = (timeDOFLoop == (nTimeDOFLoop-1)) ? NPadLastDOFs : NPad;

        /* Carry out the multiplication with the inverse of the mass matrix.
           The result is stored in resInt as temporary storage. */
        blasFunctions->gemm(nDOFs, NPadDOF, nDOFs, volElem[l].invMassMatrix.data(),
                            res, resInt, config);

        /* Make a distinction between 2D and 3D for performance reasons. */
        switch( nDim ) {
          case 2: {

            /* Loop over the time DOFs for this chunk of time DOFs. */
            for(unsigned short iTimeDOF=DOFStart; iTimeDOF<DOFEnd; ++iTimeDOF) {

              /* Index of the 1st entry of the 1st spatial DOF starts for this time DOF. */
              const unsigned int indResThisTimeDOF = nVar*(iTimeDOF-DOFStart);

              /* Loop over all the time DOFs to compute the contribution of
                 iTimeDOF to the update. */
              for(unsigned short jTimeDOF=0; jTimeDOF<nTimeDOFs; ++jTimeDOF) {

                /* Easier storage of the ADER coefficient of timeDOF jTimeDOF
                   to the residual of iTimeDOF and determine the starting
                   index for the update for timeDOF jTimeDOF. */
                const su2double coefADER = timeCoefADER_DG[jTimeDOF*nTimeDOFs+iTimeDOF];
                const unsigned int indDSolTimeInd = jTimeDOF*nVarNDOFs;

                /* Loop over the spatial DOFs for timeDOF jTimeDOF and update dSol. */
                for(unsigned short i=0; i<nDOFs; ++i) {

                  const unsigned int indResInt = indResThisTimeDOF + i*NPadDOF;
                  const unsigned int indDSol   = indDSolTimeInd    + i*nVar;

                  dSol[indDSol]   += coefADER*resInt[indResInt];
                  dSol[indDSol+1] += coefADER*resInt[indResInt+1];
                  dSol[indDSol+2] += coefADER*resInt[indResInt+2];
                  dSol[indDSol+3] += coefADER*resInt[indResInt+3];
                }
              }
            }

            break;
          }

          case 3: {

            /* Loop over the time DOFs for this chunk of time DOFs. */
            for(unsigned short iTimeDOF=DOFStart; iTimeDOF<DOFEnd; ++iTimeDOF) {

              /* Index of the 1st entry of the 1st spatial DOF starts for this time DOF. */
              const unsigned int indResThisTimeDOF = nVar*(iTimeDOF-DOFStart);

              /* Loop over all the time DOFs to compute the contribution of
                 iTimeDOF to the update. */
              for(unsigned short jTimeDOF=0; jTimeDOF<nTimeDOFs; ++jTimeDOF) {

                /* Easier storage of the ADER coefficient of timeDOF jTimeDOF
                   to the residual of iTimeDOF and determine the starting
                   index for the update for timeDOF jTimeDOF. */
                const su2double coefADER = timeCoefADER_DG[jTimeDOF*nTimeDOFs+iTimeDOF];
                const unsigned int indDSolTimeInd = jTimeDOF*nVarNDOFs;

                /* Loop over the spatial DOFs for timeDOF jTimeDOF and update dSol. */
                for(unsigned short i=0; i<nDOFs; ++i) {

                  const unsigned int indResInt = indResThisTimeDOF + i*NPadDOF;
                  const unsigned int indDSol   = indDSolTimeInd    + i*nVar;

                  dSol[indDSol]   += coefADER*resInt[indResInt];
                  dSol[indDSol+1] += coefADER*resInt[indResInt+1];
                  dSol[indDSol+2] += coefADER*resInt[indResInt+2];
                  dSol[indDSol+3] += coefADER*resInt[indResInt+3];
                  dSol[indDSol+4] += coefADER*resInt[indResInt+4];
                }
              }
            }

            break;
          }
        }
      }

      /*----------------------------------------------------------------------*/
      /*--- Determine whether or not the iteration process is considered   ---*/
      /*--- converged. Make a distinction between 2D and 3D for            ---*/
      /*--- performance reasons.                                           ---*/
      /*----------------------------------------------------------------------*/

      /* Initialize converged to true. It will be overwritten below when one
         or more DOFs in both space and time are not considered converged. */
      bool converged = true;

      switch( nDim ) {
        case 2: {
          /* 2D simulation. Loop over the time DOFs to determine the new
             solution and to determine whether or not the iteration procedure
             is considered converged. */
          for(unsigned short j=0; j<nTimeDOFs; ++j) {
            su2double       *solTimeInd  = solPred + j*nVarNDOFs;
            const su2double *dSolTimeInd = dSol    + j*nVarNDOFs;

            for(unsigned short i=0; i<nDOFs; ++i) {
              const unsigned short off = i*nVar;
              const su2double *ssolCur = solCur + off;
              const su2double *ddSol   = dSolTimeInd + off;
              su2double       *ssol    = solTimeInd  + off;

              su2double solOld;
              solOld  = ssol[0]; ssol[0] = ssolCur[0] + ddSol[0];
              solOld -= ssol[0]; if(fabs(solOld) > TolSolADER[0]) converged = false;

              solOld  = ssol[1]; ssol[1] = ssolCur[1] + ddSol[1];
              solOld -= ssol[1]; if(fabs(solOld) > TolSolADER[1]) converged = false;

              solOld  = ssol[2]; ssol[2] = ssolCur[2] + ddSol[2];
              solOld -= ssol[2]; if(fabs(solOld) > TolSolADER[2]) converged = false;

              solOld  = ssol[3]; ssol[3] = ssolCur[3] + ddSol[3];
              solOld -= ssol[3]; if(fabs(solOld) > TolSolADER[3]) converged = false;
            }
          }

          break;
        }

        case 3: {
          /* 3D simulation. Loop over the time DOFs to determine the new
             solution and to determine whether or not the iteration procedure
             is considered converged. */
          for(unsigned short j=0; j<nTimeDOFs; ++j) {
            su2double       *solTimeInd  = solPred + j*nVarNDOFs;
            const su2double *dSolTimeInd = dSol    + j*nVarNDOFs;

            for(unsigned short i=0; i<nDOFs; ++i) {
              const unsigned short off = i*nVar;
              const su2double *ssolCur = solCur + off;
              const su2double *ddSol   = dSolTimeInd + off;
              su2double       *ssol    = solTimeInd  + off;

              su2double solOld;
              solOld  = ssol[0]; ssol[0] = ssolCur[0] + ddSol[0];
              solOld -= ssol[0]; if(fabs(solOld) > TolSolADER[0]) converged = false;

              solOld  = ssol[1]; ssol[1] = ssolCur[1] + ddSol[1];
              solOld -= ssol[1]; if(fabs(solOld) > TolSolADER[1]) converged = false;

              solOld  = ssol[2]; ssol[2] = ssolCur[2] + ddSol[2];
              solOld -= ssol[2]; if(fabs(solOld) > TolSolADER[2]) converged = false;

              solOld  = ssol[3]; ssol[3] = ssolCur[3] + ddSol[3];
              solOld -= ssol[3]; if(fabs(solOld) > TolSolADER[3]) converged = false;

              solOld  = ssol[4]; ssol[4] = ssolCur[4] + ddSol[4];
              solOld -= ssol[4]; if(fabs(solOld) > TolSolADER[4]) converged = false;
            }
          }

          break;
        }
      }

      /* Break the iteration loop if the solution is considered converged. */
      if( converged ) break;
    }

    /* Store the predictor solution in the correct location of
       VecSolDOFsPredictorADER. */
    for(unsigned short j=0; j<nTimeDOFs; ++j) {
      su2double *solADERPred = VecSolDOFsPredictorADER.data()
                             + nVar*(j*nDOFsLocTot + volElem[l].offsetDOFsSolLocal);
      su2double *solPredTime = solPred + j*nVarNDOFs;

      for(unsigned short mm=0; mm<nVarNDOFs; ++mm)
        solADERPred[mm] = solPredTime[mm];
    }
  }
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                                              CVolumeElementFEM    *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  /* Get the necessary information from the standard element. */
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *matDerBasisInt         = matBasisInt + nDOFs*nInt;
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  /* Set the pointers for fluxes in the DOFs, the gradient of the fluxes in
     the integration points and the divergence of the fluxes in the integration
     points. Note that the same array can be used for the first and last array,
     because divFlux is only needed after the fluxes in the DOFs. */
  su2double *fluxXDOF     = work;
  su2double *fluxYDOF     = fluxXDOF + NPad*nDOFs;
  su2double *gradFluxXInt = fluxYDOF + NPad*nDOFs;
  su2double *gradFluxYInt = gradFluxXInt + nDim*NPad*nInt;
  su2double *divFlux      = work;

  /* Determine the offset between the r-derivatives and s-derivatives of
     the fluxes. */
  const unsigned short offDeriv = NPad*nInt;

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 5;  /* nDim*nDim + 1. */

  /*--------------------------------------------------------------------------*/
  /*---          Construct the Cartesian fluxes in the DOFs.               ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the DOFs to compute the Cartesian fluxes in the DOFs. ---*/
    for(unsigned short i=0; i<nDOFs; ++i) {

      /* Set the pointers to the locations where the solution and the Cartesian
         fluxes are stored for this DOF. */
      const unsigned short offDOF = i*NPad + simul*nVar;
      const su2double *solDOF     = sol      + offDOF;
      su2double *fluxX            = fluxXDOF + offDOF;
      su2double *fluxY            = fluxYDOF + offDOF;

      /* Set the pointer to the grid velocities at the location of the
         solution DOFS. THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL
         MOTION IS SPECIFIED, THE DATA FOR THIS DOF FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *gridVel = elem->gridVelocitiesSolDOFs.data() + 2*i; /* nDim*i. */

      /* Compute the velocities and pressure in this DOF. */
      const su2double DensityInv   = 1.0/solDOF[0];
      const su2double u            = DensityInv*solDOF[1];
      const su2double v            = DensityInv*solDOF[2];
      const su2double StaticEnergy = DensityInv*solDOF[3] - 0.5*(u*u + v*v);

      FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
      const su2double Pressure = FluidModel->GetPressure();

      /* The Cartesian fluxes in the x-direction. */
      const su2double uRel = u - gridVel[0];
      fluxX[0] = solDOF[0]*uRel;
      fluxX[1] = solDOF[1]*uRel + Pressure;
      fluxX[2] = solDOF[2]*uRel;
      fluxX[3] = solDOF[3]*uRel + Pressure*u;

      /* The Cartesian fluxes in the y-direction. */
      const su2double vRel = v - gridVel[1];
      fluxY[0] = solDOF[0]*vRel;
      fluxY[1] = solDOF[1]*vRel;
      fluxY[2] = solDOF[2]*vRel + Pressure;
      fluxY[3] = solDOF[3]*vRel + Pressure*v;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the derivatives of the Cartesian fluxes w.r.t. the         ---*/
  /*--- parametric coordinates in the integration points.                  ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxXDOF, gradFluxXInt, config);
  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxYDOF, gradFluxYInt, config);

  /*--------------------------------------------------------------------------*/
  /*--- Compute the divergence of the fluxes in the integration points,    ---*/
  /*--- multiplied by the integration weight.                              ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the integration points of the element. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Set the pointers for this integration point. */
      const unsigned short offInt  = i*NPad + simul*nVar;
      const su2double *gradFluxXDr = gradFluxXInt + offInt;
      const su2double *gradFluxXDs = gradFluxXDr  + offDeriv;
      const su2double *gradFluxYDr = gradFluxYInt + offInt;
      const su2double *gradFluxYDs = gradFluxYDr  + offDeriv;
      su2double       *divFluxInt  = divFlux      + offInt;

      /* Easier storage of the metric terms in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
         BE TAKEN. */
      const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

      /* Compute the metric terms multiplied by the integration weight. Note that the
         first term in the metric terms is the Jacobian. */
      const su2double wDrdx = weights[i]*metricTerms[1];
      const su2double wDrdy = weights[i]*metricTerms[2];

      const su2double wDsdx = weights[i]*metricTerms[3];
      const su2double wDsdy = weights[i]*metricTerms[4];

      /* Compute the divergence of the fluxes, multiplied by the
         integration weight. */
      divFluxInt[0] = gradFluxXDr[0]*wDrdx + gradFluxXDs[0]*wDsdx
                    + gradFluxYDr[0]*wDrdy + gradFluxYDs[0]*wDsdy;
      divFluxInt[1] = gradFluxXDr[1]*wDrdx + gradFluxXDs[1]*wDsdx
                    + gradFluxYDr[1]*wDrdy + gradFluxYDs[1]*wDsdy;
      divFluxInt[2] = gradFluxXDr[2]*wDrdx + gradFluxXDs[2]*wDsdx
                    + gradFluxYDr[2]*wDrdy + gradFluxYDs[2]*wDsdy;
      divFluxInt[3] = gradFluxXDr[3]*wDrdx + gradFluxXDs[3]*wDsdx
                    + gradFluxYDr[3]*wDrdy + gradFluxYDs[3]*wDsdy;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the body force to the divergence of the fluxes, if such a      ---*/
  /*--- force is present.                                                  ---*/
  /*--------------------------------------------------------------------------*/

  if( config->GetBody_Force() ) {

    /* Easier storage of the body force. */
    const su2double *body_force_vector = config->GetBody_Force_Vector();

    /* Compute the solution in the integration points of the element.
       Use gradFluxYInt to store this solution. */
    su2double *solInt = gradFluxYInt;

    blasFunctions->gemm(nInt, NPad, nDOFs, matBasisInt, sol, solInt, config);

    /*--- Loop over the number of entities that are treated simultaneously. */
    for(unsigned short simul=0; simul<nSimul; ++simul) {

      /*--- Loop over the integration points of the element. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Set the pointers for this integration point. */
        const unsigned short offInt  = i*NPad + simul*nVar;
        const su2double *solThisInt = solInt  + offInt;
        su2double       *divFluxInt = divFlux + offInt;

        /* Easier storage of the metric terms in this integration point.
           THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
           THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
           BE TAKEN. */
        const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

        /* Compute the velocities. */
        const su2double rhoInv = 1.0/solThisInt[0];
        const su2double u      = solThisInt[1]*rhoInv;
        const su2double v      = solThisInt[2]*rhoInv;

        /* Add the body force to the flux divergence for the momentum and energy
           equation. Note that the source terms are multiplied with minus the
           integration weight in order to be consistent with the formulation of
           the residual. Also note that for the energy source term the absolute
           velocity must be taken and not the relative. */
        const su2double weightJac = weights[i]*metricTerms[0];

        divFluxInt[1] -= weightJac*body_force_vector[0];
        divFluxInt[2] -= weightJac*body_force_vector[1];
        divFluxInt[3] -= weightJac*(u*body_force_vector[0] + v*body_force_vector[1]);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the source terms of the manufactured solution to the divergence---*/
  /*--- of the fluxes, if a manufactured solution is used.                 ---*/
  /*--------------------------------------------------------------------------*/

  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Loop over the number of entities that are treated simultaneously. */
      for(unsigned short simul=0; simul<nSimul; ++simul) {

        /*--- Loop over the integration points of the element. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Set the pointers for this integration point. */
          const unsigned short offInt  = i*NPad + simul*nVar;
          su2double       *divFluxInt = divFlux + offInt;

          /* Set the pointer to the coordinates in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *coor = elem->coorIntegrationPoints.data() + i*nDim;

          /* Easier storage of the metric terms in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
          const su2double weightJac    = weights[i]*metricTerms[0];

          /* Compute the source terms of the manufactured solution.
             THIS IS A TEMPORARY IMPLEMENTATION. FOR AN ACTUAL TIME ACCURATE
             SIMULATION THE CORRECT TIME MUST BE GIVEN TO THIS FUNCTION. */
          su2double sourceMan[4];
          VerificationSolution->GetMMSSourceTerm(coor, 0.0, sourceMan);

          /* Add the source terms to the flux divergence. Note that the source
             terms are multiplied with minus the integration weight in order
             to be consistent with the formulation of the residual. */
          divFluxInt[0] -= weightJac*sourceMan[0];
          divFluxInt[1] -= weightJac*sourceMan[1];
          divFluxInt[2] -= weightJac*sourceMan[2];
          divFluxInt[3] -= weightJac*sourceMan[3];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the residual in the DOFs, which is the matrix product of   ---*/
  /*--- basisFunctionsIntTrans and divFlux.                                ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nDOFs, NPad, nInt, basisFunctionsIntTrans, divFlux, res, config);
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                                              CVolumeElementFEM    *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  /* Get the necessary information from the standard element. */
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *matDerBasisInt         = matBasisInt + nDOFs*nInt;
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  /* Set the pointers for fluxes in the DOFs, the gradient of the fluxes in
     the integration points and the divergence of the fluxes in the integration
     points. Note that the same array can be used for the first and last array,
     because divFlux is only needed after the fluxes in the DOFs. */
  su2double *fluxXDOF     = work;
  su2double *fluxYDOF     = fluxXDOF + NPad*nDOFs;
  su2double *fluxZDOF     = fluxYDOF + NPad*nDOFs;
  su2double *gradFluxXInt = fluxZDOF + NPad*nDOFs;
  su2double *gradFluxYInt = gradFluxXInt + nDim*NPad*nInt;
  su2double *gradFluxZInt = gradFluxYInt + nDim*NPad*nInt;
  su2double *divFlux      = work;

  /* Determine the offset between the r-derivatives and s-derivatives, which is
     also the offset between s- and t-derivatives, of the fluxes. */
  const unsigned short offDeriv = NPad*nInt;

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 10;  /* nDim*nDim + 1. */

  /*--------------------------------------------------------------------------*/
  /*---          Construct the Cartesian fluxes in the DOFs.               ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the DOFs to compute the Cartesian fluxes in the DOFs. ---*/
    for(unsigned short i=0; i<nDOFs; ++i) {

      /* Set the pointers to the locations where the solution and the Cartesian
         fluxes are stored for this DOF. */
      const unsigned short offDOF = i*NPad + simul*nVar;
      const su2double *solDOF     = sol      + offDOF;
      su2double *fluxX            = fluxXDOF + offDOF;
      su2double *fluxY            = fluxYDOF + offDOF;
      su2double *fluxZ            = fluxZDOF + offDOF;

      /* Set the pointer to the grid velocities at the location of the
         solution DOFS. THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL
         MOTION IS SPECIFIED, THE DATA FOR THIS DOF FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *gridVel = elem->gridVelocitiesSolDOFs.data() + 3*i; /* nDim*i. */

      /* Compute the velocities and pressure in this DOF. */
      const su2double DensityInv   = 1.0/solDOF[0];
      const su2double u            = DensityInv*solDOF[1];
      const su2double v            = DensityInv*solDOF[2];
      const su2double w            = DensityInv*solDOF[3];
      const su2double StaticEnergy = DensityInv*solDOF[4] - 0.5*(u*u + v*v + w*w);

      FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
      const su2double Pressure = FluidModel->GetPressure();

      /* The Cartesian fluxes in the x-direction. */
      const su2double uRel = u - gridVel[0];
      fluxX[0] = solDOF[0]*uRel;
      fluxX[1] = solDOF[1]*uRel + Pressure;
      fluxX[2] = solDOF[2]*uRel;
      fluxX[3] = solDOF[3]*uRel;
      fluxX[4] = solDOF[4]*uRel + Pressure*u;

      /* The Cartesian fluxes in the y-direction. */
      const su2double vRel = v - gridVel[1];
      fluxY[0] = solDOF[0]*vRel;
      fluxY[1] = solDOF[1]*vRel;
      fluxY[2] = solDOF[2]*vRel + Pressure;
      fluxY[3] = solDOF[3]*vRel;
      fluxY[4] = solDOF[4]*vRel + Pressure*v;

      /* The Cartesian fluxes in the z-direction. */
      const su2double wRel = w - gridVel[2];
      fluxZ[0] = solDOF[0]*wRel;
      fluxZ[1] = solDOF[1]*wRel;
      fluxZ[2] = solDOF[2]*wRel;
      fluxZ[3] = solDOF[3]*wRel + Pressure;
      fluxZ[4] = solDOF[4]*wRel + Pressure*w;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the derivatives of the Cartesian fluxes w.r.t. the         ---*/
  /*--- parametric coordinates in the integration points.                  ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxXDOF, gradFluxXInt, config);
  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxYDOF, gradFluxYInt, config);
  blasFunctions->gemm(nInt*nDim, NPad, nDOFs, matDerBasisInt, fluxZDOF, gradFluxZInt, config);

  /*--------------------------------------------------------------------------*/
  /*--- Compute the divergence of the fluxes in the integration points,    ---*/
  /*--- multiplied by the integration weight.                              ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the integration points of the element. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Set the pointers for this integration point. */
      const unsigned short offInt  = i*NPad + simul*nVar;
      const su2double *gradFluxXDr = gradFluxXInt + offInt;
      const su2double *gradFluxXDs = gradFluxXDr  + offDeriv;
      const su2double *gradFluxXDt = gradFluxXDs  + offDeriv;
      const su2double *gradFluxYDr = gradFluxYInt + offInt;
      const su2double *gradFluxYDs = gradFluxYDr  + offDeriv;
      const su2double *gradFluxYDt = gradFluxYDs  + offDeriv;
      const su2double *gradFluxZDr = gradFluxZInt + offInt;
      const su2double *gradFluxZDs = gradFluxZDr  + offDeriv;
      const su2double *gradFluxZDt = gradFluxZDs  + offDeriv;
      su2double       *divFluxInt  = divFlux      + offInt;

      /* Easier storage of the metric terms in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

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
      divFluxInt[0] = gradFluxXDr[0]*wDrdx + gradFluxXDs[0]*wDsdx + gradFluxXDt[0]*wDtdx
                    + gradFluxYDr[0]*wDrdy + gradFluxYDs[0]*wDsdy + gradFluxYDt[0]*wDtdy
                    + gradFluxZDr[0]*wDrdz + gradFluxZDs[0]*wDsdz + gradFluxZDt[0]*wDtdz;
      divFluxInt[1] = gradFluxXDr[1]*wDrdx + gradFluxXDs[1]*wDsdx + gradFluxXDt[1]*wDtdx
                    + gradFluxYDr[1]*wDrdy + gradFluxYDs[1]*wDsdy + gradFluxYDt[1]*wDtdy
                    + gradFluxZDr[1]*wDrdz + gradFluxZDs[1]*wDsdz + gradFluxZDt[1]*wDtdz;
      divFluxInt[2] = gradFluxXDr[2]*wDrdx + gradFluxXDs[2]*wDsdx + gradFluxXDt[2]*wDtdx
                    + gradFluxYDr[2]*wDrdy + gradFluxYDs[2]*wDsdy + gradFluxYDt[2]*wDtdy
                    + gradFluxZDr[2]*wDrdz + gradFluxZDs[2]*wDsdz + gradFluxZDt[2]*wDtdz;
      divFluxInt[3] = gradFluxXDr[3]*wDrdx + gradFluxXDs[3]*wDsdx + gradFluxXDt[3]*wDtdx
                    + gradFluxYDr[3]*wDrdy + gradFluxYDs[3]*wDsdy + gradFluxYDt[3]*wDtdy
                    + gradFluxZDr[3]*wDrdz + gradFluxZDs[3]*wDsdz + gradFluxZDt[3]*wDtdz;
      divFluxInt[4] = gradFluxXDr[4]*wDrdx + gradFluxXDs[4]*wDsdx + gradFluxXDt[4]*wDtdx
                    + gradFluxYDr[4]*wDrdy + gradFluxYDs[4]*wDsdy + gradFluxYDt[4]*wDtdy
                    + gradFluxZDr[4]*wDrdz + gradFluxZDs[4]*wDsdz + gradFluxZDt[4]*wDtdz;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the body force to the divergence of the fluxes, if such a      ---*/
  /*--- force is present.                                                  ---*/
  /*--------------------------------------------------------------------------*/

  if( config->GetBody_Force() ) {

    /* Easier storage of the body force. */
    const su2double *body_force_vector = config->GetBody_Force_Vector();

    /* Compute the solution in the integration points of the element.
       Use gradFluxYInt to store this solution. */
    su2double *solInt = gradFluxYInt;

    blasFunctions->gemm(nInt, NPad, nDOFs, matBasisInt, sol, solInt, config);

    /*--- Loop over the number of entities that are treated simultaneously. */
    for(unsigned short simul=0; simul<nSimul; ++simul) {

      /*--- Loop over the integration points of the element. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /* Set the pointers for this integration point. */
        const unsigned short offInt = i*NPad + simul*nVar;
        const su2double *solThisInt = solInt  + offInt;
        su2double       *divFluxInt = divFlux + offInt;

        /* Easier storage of the metric terms in this integration point.
           THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
           THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
           BE TAKEN. */
        const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

        /* Compute the velocities. */
        const su2double rhoInv = 1.0/solThisInt[0];
        const su2double u      = solThisInt[1]*rhoInv;
        const su2double v      = solThisInt[2]*rhoInv;
        const su2double w      = solThisInt[3]*rhoInv;

        /* Add the body force to the flux divergence for the momentum and energy
           equation. Note that the source terms are multiplied with minus the
           integration weight in order to be consistent with the formulation of
           the residual. Also note that for the energy source term the absolute
           velocity must be taken and not the relative. */
        const su2double weightJac = weights[i]*metricTerms[0];

        divFluxInt[1] -= weightJac*body_force_vector[0];
        divFluxInt[2] -= weightJac*body_force_vector[1];
        divFluxInt[3] -= weightJac*body_force_vector[2];
        divFluxInt[4] -= weightJac*(u*body_force_vector[0] + v*body_force_vector[1]
                       +            w*body_force_vector[2]);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the source terms of the manufactured solution to the divergence---*/
  /*--- of the fluxes, if a manufactured solution is used.                 ---*/
  /*--------------------------------------------------------------------------*/

  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Loop over the number of entities that are treated simultaneously. */
      for(unsigned short simul=0; simul<nSimul; ++simul) {

        /*--- Loop over the integration points of the element. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Set the pointers for this integration point. */
          const unsigned short offInt  = i*NPad + simul*nVar;
          su2double       *divFluxInt = divFlux + offInt;

          /* Set the pointer to the coordinates in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *coor = elem->coorIntegrationPoints.data() + i*nDim;

          /* Easier storage of the metric terms in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
          const su2double weightJac    = weights[i]*metricTerms[0];

          /* Compute the source terms of the manufactured solution.
             THIS IS A TEMPORARY IMPLEMENTATION. FOR AN ACTUAL TIME ACCURATE
             SIMULATION THE CORRECT TIME MUST BE GIVEN TO THIS FUNCTION. */
          su2double sourceMan[5];
          VerificationSolution->GetMMSSourceTerm(coor, 0.0, sourceMan);

          /* Add the source terms to the flux divergence. Note that the source
             terms are multiplied with minus the integration weight in order
             to be consistent with the formulation of the residual. */
          divFluxInt[0] -= weightJac*sourceMan[0];
          divFluxInt[1] -= weightJac*sourceMan[1];
          divFluxInt[2] -= weightJac*sourceMan[2];
          divFluxInt[3] -= weightJac*sourceMan[3];
          divFluxInt[4] -= weightJac*sourceMan[4];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the residual in the DOFs, which is the matrix product of   ---*/
  /*--- basisFunctionsIntTrans and divFlux.                                ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nDOFs, NPad, nInt, basisFunctionsIntTrans, divFlux, res, config);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                                                 CVolumeElementFEM    *elem,
                                                                 const su2double      *sol,
                                                                 const unsigned short nSimul,
                                                                 const unsigned short NPad,
                                                                 su2double            *res,
                                                                 su2double            *work) {

  /* Set the pointers for solAndGradInt and divFlux to work. The same array
     can be used for both help arrays. */
  su2double *solAndGradInt = work;
  su2double *divFlux       = work;

  /* Get the necessary information from the standard element. */
  const unsigned short ind                = elem->indStandardElement;
  const unsigned short nInt               = standardElementsSol[ind].GetNIntegration();
  const unsigned short nDOFs              = elem->nDOFsSol;
  const su2double *matBasisInt            = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
  const su2double *basisFunctionsIntTrans = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
  const su2double *weights                = standardElementsSol[ind].GetWeightsIntegration();

  /* Check if a body force is present and set it accordingly. */
  su2double bodyForceX = 0.0, bodyForceY = 0.0;
  if( config->GetBody_Force() ) {
    const su2double *body_force_vector = config->GetBody_Force_Vector();
    bodyForceX = body_force_vector[0];
    bodyForceY = body_force_vector[1];
  }

  /* Determine the offset between the solution variables and the r-derivatives,
     which is also the offset between the r- and s-derivatives. */
  const unsigned short offDeriv = NPad*nInt;

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 5;  /* nDim*nDim + 1. */

  /*--------------------------------------------------------------------------*/
  /*--- Compute the solution and the derivatives w.r.t. the parametric     ---*/
  /*--- coordinates in the integration points. The first argument in       ---*/
  /*--- the call to blasFunctions->gemm is nInt*(nDim+1).                  ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nInt*3, NPad, nDOFs, matBasisInt, sol, solAndGradInt, config);

  /*--------------------------------------------------------------------------*/
  /*--- Compute the divergence of the inviscid fluxes, multiplied by the   ---*/
  /*--- integration weight in the integration points of the element.       ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /*--- Loop over the integration points. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the location where the solution, gradient data and
         the divergence of this integration point starts. */
      const unsigned short offInt = NPad*i + simul*nVar;
      const su2double *sol   = solAndGradInt + offInt;
      const su2double *solDr = sol     + offDeriv;
      const su2double *solDs = solDr   + offDeriv;
      su2double *divFluxInt  = divFlux + offInt;

      /*--- Compute the velocities, pressure and total enthalpy
            in this integration point. ---*/
      const su2double rho = sol[0];
      const su2double ru  = sol[1];
      const su2double rv  = sol[2];
      const su2double rE  = sol[3];

      const su2double rhoInv       = 1.0/rho;
      const su2double u            = rhoInv*ru;
      const su2double v            = rhoInv*rv;
      const su2double kinEnergy    = 0.5*(u*u + v*v);
      const su2double StaticEnergy = rhoInv*rE - kinEnergy;

      FluidModel->SetTDState_rhoe(rho, StaticEnergy);
      const su2double Pressure = FluidModel->GetPressure();
      const su2double Htot     = rhoInv*(rE + Pressure);

      /* Set the pointer to the grid velocities in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *gridVel = elem->gridVelocities.data() + 2*i; /* nDim*i. */

      /* Easier storage of the metric terms in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;

      const su2double drdx = metricTerms[1];
      const su2double drdy = metricTerms[2];
      const su2double dsdx = metricTerms[3];
      const su2double dsdy = metricTerms[4];

      /* Compute the Cartesian gradients of the independent solution
         variables from the gradients in parametric coordinates and the
         metric terms in this integration point. Note that these gradients
         must be scaled with the Jacobian. This scaling is already present
         in the metric terms. */
      const su2double drhodx = solDr[0]*drdx + solDs[0]*dsdx;
      const su2double drudx  = solDr[1]*drdx + solDs[1]*dsdx;
      const su2double drvdx  = solDr[2]*drdx + solDs[2]*dsdx;
      const su2double drEdx  = solDr[3]*drdx + solDs[3]*dsdx;

      const su2double drhody = solDr[0]*drdy + solDs[0]*dsdy;
      const su2double drudy  = solDr[1]*drdy + solDs[1]*dsdy;
      const su2double drvdy  = solDr[2]*drdy + solDs[2]*dsdy;
      const su2double drEdy  = solDr[3]*drdy + solDs[3]*dsdy;

      /*--- Compute the divergence of the grid velocity.
            SET TO ZERO FOR NOW. THIS IS NOT CORRECT!!!!. ---*/
      const su2double divGridVel = 0.0;  // Must be multiplied by the Jacobian.

      /* Compute the Cartesian gradients of the pressure, scaled with
         the Jacobian. */
      const su2double dpdx = Gamma_Minus_One*(drEdx + kinEnergy*drhodx
                           -                  u*drudx - v*drvdx);
      const su2double dpdy = Gamma_Minus_One*(drEdy + kinEnergy*drhody
                           -                  u*drudy - v*drvdy);

      /* Abbreviations, which make it easier to compute the divergence term. */
      const su2double abv1 =    drudx  +   drvdy;
      const su2double abv2 = u* drhodx + v*drhody;
      const su2double abv3 = u*(drEdx  + dpdx) + v*(drEdy + dpdy);

      /* Compute the divergence terms for this integration point,
         multiplied by the integration weight. */
      divFluxInt[0] = weights[i]*(abv1 - rho*divGridVel
                    -             gridVel[0]*drhodx - gridVel[1]*drhody);
      divFluxInt[1] = weights[i]*(dpdx + u*(abv1 - abv2 + drudx) + v*drudy
                    -             ru*divGridVel - gridVel[0]*drudx - gridVel[1]*drudy);
      divFluxInt[2] = weights[i]*(dpdy + v*(abv1 - abv2 + drvdy) + u*drvdx
                    -             rv*divGridVel - gridVel[0]*drvdx - gridVel[1]*drvdy);
      divFluxInt[3] = weights[i]*(abv3 + Htot*(abv1 - abv2) - rE*divGridVel
                    -             gridVel[0]*drEdx - gridVel[1]*drEdy);

      /* Add the body force to the flux divergence for the momentum and energy
         equation. Note that the source terms are multiplied with minus the
         integration weight in order to be consistent with the formulation of
         the residual. Also note that for the energy source term the absolute
         velocity must be taken and not the relative. */
      const su2double weightJac = weights[i]*metricTerms[0];

      divFluxInt[1] -= weightJac*bodyForceX;
      divFluxInt[2] -= weightJac*bodyForceY;
      divFluxInt[3] -= weightJac*(u*bodyForceX + v*bodyForceY);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the source terms of the manufactured solution to the divergence---*/
  /*--- of the fluxes, if a manufactured solution is used.                 ---*/
  /*--------------------------------------------------------------------------*/

  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Loop over the number of entities that are treated simultaneously. */
      for(unsigned short simul=0; simul<nSimul; ++simul) {

        /*--- Loop over the integration points of the element. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Set the pointers for this integration point. */
          const unsigned short offInt  = i*NPad + simul*nVar;
          su2double       *divFluxInt = divFlux + offInt;

          /* Set the pointer to the coordinates in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *coor = elem->coorIntegrationPoints.data() + i*nDim;

          /* Easier storage of the metric terms in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
          const su2double weightJac    = weights[i]*metricTerms[0];

          /* Compute the source terms of the manufactured solution.
             THIS IS A TEMPORARY IMPLEMENTATION. FOR AN ACTUAL TIME ACCURATE
             SIMULATION THE CORRECT TIME MUST BE GIVEN TO THIS FUNCTION. */
          su2double sourceMan[4];
          VerificationSolution->GetMMSSourceTerm(coor, 0.0, sourceMan);

          /* Add the source terms to the flux divergence. Note that the source
             terms are multiplied with minus the integration weight in order
             to be consistent with the formulation of the residual. */
          divFluxInt[0] -= weightJac*sourceMan[0];
          divFluxInt[1] -= weightJac*sourceMan[1];
          divFluxInt[2] -= weightJac*sourceMan[2];
          divFluxInt[3] -= weightJac*sourceMan[3];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the residual in the DOFs, which is the matrix product of   ---*/
  /*--- basisFunctionsIntTrans and divFlux.                                ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nDOFs, NPad, nInt, basisFunctionsIntTrans, divFlux, res, config);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                                                 CVolumeElementFEM    *elem,
                                                                 const su2double      *sol,
                                                                 const unsigned short nSimul,
                                                                 const unsigned short NPad,
                                                                 su2double            *res,
                                                                 su2double            *work) {

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

  /* Check if a body force is present and set it accordingly. */
  su2double bodyForceX = 0.0, bodyForceY = 0.0, bodyForceZ = 0.0;
  if( config->GetBody_Force() ) {
    const su2double *body_force_vector = config->GetBody_Force_Vector();
    bodyForceX = body_force_vector[0];
    bodyForceY = body_force_vector[1];
    bodyForceZ = body_force_vector[2];
  }

  /* Determine the offset between the solution variables and the r-derivatives,
     which is also the offset between the r- and s-derivatives and the offset
     between s- and t-derivatives. */
  const unsigned short offDeriv = NPad*nInt;

  /* Store the number of metric points per integration point for readability. */
  const unsigned short nMetricPerPoint = 10;  /* nDim*nDim + 1. */

  /*--------------------------------------------------------------------------*/
  /*--- Compute the solution and the derivatives w.r.t. the parametric     ---*/
  /*--- coordinates in the integration points. The first argument in       ---*/
  /*--- the call to blasFunctions->gemm is nInt*(nDim+1).                  ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nInt*4, NPad, nDOFs, matBasisInt, sol, solAndGradInt, config);

  /*--- Loop over the number of entities that are treated simultaneously. */
  for(unsigned short simul=0; simul<nSimul; ++simul) {

    /* Set the pointers to the location where the solution and the divergence
       of the 1st DOF of this entity is stored, such that this offset does not
       need to be taken into account anymore in the loop over the integration
       points below. */
    const su2double *solAndGradEntity = solAndGradInt + simul*nVar;
    su2double       *divFluxEntity    = divFlux       + simul*nVar;

    /*--- Loop over the integration points. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the location where the solution and gradient data
         of this integration point starts. */
      const su2double *sol   = solAndGradEntity + NPad*i;
      const su2double *solDr = sol   + offDeriv;
      const su2double *solDs = solDr + offDeriv;
      const su2double *solDt = solDs + offDeriv;
      su2double *divFluxInt  = divFluxEntity + NPad*i;

        /*--- Compute the velocities, pressure and total enthalpy
            in this integration point. ---*/
      const su2double rho = sol[0];
      const su2double ru  = sol[1];
      const su2double rv  = sol[2];
      const su2double rw  = sol[3];
      const su2double rE  = sol[4];

      const su2double rhoInv       = 1.0/rho;
      const su2double u            = rhoInv*ru;
      const su2double v            = rhoInv*rv;
      const su2double w            = rhoInv*rw;
      const su2double kinEnergy    = 0.5*(u*u + v*v + w*w);
      const su2double StaticEnergy = rhoInv*rE - kinEnergy;

      FluidModel->SetTDState_rhoe(rho, StaticEnergy);
      const su2double Pressure = FluidModel->GetPressure();
      const su2double Htot     = rhoInv*(rE + Pressure);

      /* Set the pointer to the grid velocities in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
      const su2double *gridVel = elem->gridVelocities.data() + 3*i; /* nDim*i. */

      /* Easier storage of the metric terms in this integration point.
         THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
         THE DATA FOR THIS SPATIAL INTEGRATION POINT FOR THE CURRENT TIME
         INTEGRATION POINT MUST BE TAKEN. */
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

      /* Compute the Cartesian gradients of the independent solution
         variables from the gradients in parametric coordinates and the
         metric terms in this integration point. Note that these gradients
         must be scaled with the Jacobian. This scaling is already present
         in the metric terms. */
      const su2double drhodx = solDr[0]*drdx + solDs[0]*dsdx + solDt[0]*dtdx;
      const su2double drudx  = solDr[1]*drdx + solDs[1]*dsdx + solDt[1]*dtdx;
      const su2double drvdx  = solDr[2]*drdx + solDs[2]*dsdx + solDt[2]*dtdx;
      const su2double drwdx  = solDr[3]*drdx + solDs[3]*dsdx + solDt[3]*dtdx;
      const su2double drEdx  = solDr[4]*drdx + solDs[4]*dsdx + solDt[4]*dtdx;

      const su2double drhody = solDr[0]*drdy + solDs[0]*dsdy + solDs[0]*dtdy;
      const su2double drudy  = solDr[1]*drdy + solDs[1]*dsdy + solDs[1]*dtdy;
      const su2double drvdy  = solDr[2]*drdy + solDs[2]*dsdy + solDs[2]*dtdy;
      const su2double drwdy  = solDr[3]*drdy + solDs[3]*dsdy + solDs[3]*dtdy;
      const su2double drEdy  = solDr[4]*drdy + solDs[4]*dsdy + solDs[4]*dtdy;

      const su2double drhodz = solDr[0]*drdz + solDs[0]*dsdz + solDs[0]*dtdz;
      const su2double drudz  = solDr[1]*drdz + solDs[1]*dsdz + solDs[1]*dtdz;
      const su2double drvdz  = solDr[2]*drdz + solDs[2]*dsdz + solDs[2]*dtdz;
      const su2double drwdz  = solDr[3]*drdz + solDs[3]*dsdz + solDs[3]*dtdz;
      const su2double drEdz  = solDr[4]*drdz + solDs[4]*dsdz + solDs[4]*dtdz;

      /*--- Compute the divergence of the grid velocity.
            SET TO ZERO FOR NOW. THIS IS NOT CORRECT!!!!. ---*/
      const su2double divGridVel = 0.0;   // Must be multiplied by the Jacobian.

      /* Compute the Cartesian gradients of the pressure, scaled with
         the Jacobian. */
      const su2double dpdx = Gamma_Minus_One*(drEdx + kinEnergy*drhodx
                           -                  u*drudx - v*drvdx - w*drwdx);
      const su2double dpdy = Gamma_Minus_One*(drEdy + kinEnergy*drhody
                           -                  u*drudy - v*drvdy - w*drwdy);
      const su2double dpdz = Gamma_Minus_One*(drEdz + kinEnergy*drhodz
                           -                  u*drudz - v*drvdz - w*drwdz);

      /* Abbreviations, which make it easier to compute the divergence term. */
      const su2double abv1 =    drudx  +   drvdy  +   drwdz;
      const su2double abv2 = u* drhodx + v*drhody + w*drhodz;
      const su2double abv3 = u*(drEdx + dpdx) + v*(drEdy + dpdy) + w*(drEdz + dpdz);

      /* Compute the divergence terms for this integration point,
         multiplied by the integration weight. */
      divFluxInt[0] = weights[i]*(abv1 - rho*divGridVel
                    -             gridVel[0]*drhodx - gridVel[1]*drhody - gridVel[2]*drhodz);
      divFluxInt[1] = weights[i]*(dpdx + u*(abv1 - abv2 + drudx) + v*drudy + w*drudz
                    -             ru*divGridVel - gridVel[0]*drudx - gridVel[1]*drudy - gridVel[2]*drudz);
      divFluxInt[2] = weights[i]*(dpdy + v*(abv1 - abv2 + drvdy) + u*drvdx + w*drvdz
                    -             rv*divGridVel - gridVel[0]*drvdx - gridVel[1]*drvdy - gridVel[2]*drvdz);
      divFluxInt[3] = weights[i]*(dpdz + w*(abv1 - abv2 + drwdz) + u*drwdx + v*drwdy
                    -             rw*divGridVel - gridVel[0]*drwdx - gridVel[1]*drwdy - gridVel[2]*drwdz);
      divFluxInt[4] = weights[i]*(abv3 + Htot*(abv1 - abv2) - rE*divGridVel
                    -             gridVel[0]*drEdx - gridVel[1]*drEdy - gridVel[2]*drEdz);

      /* Add the body force to the flux divergence for the momentum and energy
         equation. Note that the source terms are multiplied with minus the
         integration weight in order to be consistent with the formulation of
         the residual. Also note that for the energy source term the absolute
         velocity must be taken and not the relative. */
      const su2double weightJac = weights[i]*metricTerms[0];

      divFluxInt[1] -= weightJac*bodyForceX;
      divFluxInt[2] -= weightJac*bodyForceY;
      divFluxInt[3] -= weightJac*bodyForceZ;
      divFluxInt[4] -= weightJac*(u*bodyForceX + v*bodyForceY + w*bodyForceZ);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Add the source terms of the manufactured solution to the divergence---*/
  /*--- of the fluxes, if a manufactured solution is used.                 ---*/
  /*--------------------------------------------------------------------------*/

  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Loop over the number of entities that are treated simultaneously. */
      for(unsigned short simul=0; simul<nSimul; ++simul) {

        /*--- Loop over the integration points of the element. ---*/
        for(unsigned short i=0; i<nInt; ++i) {

          /* Set the pointers for this integration point. */
          const unsigned short offInt  = i*NPad + simul*nVar;
          su2double       *divFluxInt = divFlux + offInt;

          /* Set the pointer to the coordinates in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *coor = elem->coorIntegrationPoints.data() + i*nDim;

          /* Easier storage of the metric terms in this integration point.
             THIS IS A TEMPORARY IMPLEMENTATION. WHEN AN ACTUAL MOTION IS SPECIFIED,
             THE DATA FOR THIS DOF FOR THE CURRENT TIME INTEGRATION POINT MUST
             BE TAKEN. */
          const su2double *metricTerms = elem->metricTerms.data() + i*nMetricPerPoint;
          const su2double weightJac    = weights[i]*metricTerms[0];

          /* Compute the source terms of the manufactured solution.
             THIS IS A TEMPORARY IMPLEMENTATION. FOR AN ACTUAL TIME ACCURATE
             SIMULATION THE CORRECT TIME MUST BE GIVEN TO THIS FUNCTION. */
          su2double sourceMan[5];
          VerificationSolution->GetMMSSourceTerm(coor, 0.0, sourceMan);

          /* Add the source terms to the flux divergence. Note that the source
             terms are multiplied with minus the integration weight in order
             to be consistent with the formulation of the residual. */
          divFluxInt[0] -= weightJac*sourceMan[0];
          divFluxInt[1] -= weightJac*sourceMan[1];
          divFluxInt[2] -= weightJac*sourceMan[2];
          divFluxInt[3] -= weightJac*sourceMan[3];
          divFluxInt[4] -= weightJac*sourceMan[4];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Compute the residual in the DOFs, which is the matrix product of   ---*/
  /*--- basisFunctionsIntTrans and divFlux.                                ---*/
  /*--------------------------------------------------------------------------*/

  blasFunctions->gemm(nDOFs, NPad, nInt, basisFunctionsIntTrans, divFlux, res, config);
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
                                             const unsigned long elemEnd,
                                             su2double           *workArray) {

  /*--- Run shock capturing algorithm ---*/
  switch( config->GetKind_FEM_DG_Shock() ) {
    case FEM_SHOCK_CAPTURING_DG::NONE:
      break;
    case FEM_SHOCK_CAPTURING_DG::PERSSON:
      // to be done
      //Shock_Capturing_DG_Persson(elemBeg, elemEnd, workArray);
      break;
  }
}

void CFEM_DG_EulerSolver::Volume_Residual(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd,
                                          su2double           *workArray) {

  /*--- Determine whether a body force term is present. ---*/
  bool body_force = config->GetBody_Force();
  const su2double *body_force_vector = body_force ? config->GetBody_Force_Vector() : nullptr;

  /*--- Get the physical time for MMS if necessary. ---*/
  su2double time = 0.0;
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

  /* Determine the number of elements that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nElemSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

  /*--- Loop over the given element range to compute the contribution of the
        volume integral in the DG FEM formulation to the residual. Multiple
        elements are treated simultaneously to improve the performance
        of the matrix multiplications. As a consequence, the update of the
        counter l happens at the end of this loop section. ---*/
  for(unsigned long l=elemBeg; l<elemEnd;) {

    /* Determine the end index for this chunk of elements and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(volElem, l, elemEnd, nElemSimul, nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the required data from the corresponding standard element. */
    const unsigned short nInt            = standardElementsSol[ind].GetNIntegration();
    const unsigned short nDOFs           = volElem[l].nDOFsSol;
    const su2double *matBasisInt         = standardElementsSol[ind].GetMatBasisFunctionsIntegration();
    const su2double *matBasisIntTrans    = standardElementsSol[ind].GetBasisFunctionsIntegrationTrans();
    const su2double *matDerBasisIntTrans = standardElementsSol[ind].GetDerMatBasisFunctionsIntTrans();
    const su2double *weights             = standardElementsSol[ind].GetWeightsIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solDOFs = workArray;
    su2double *sources = solDOFs + nDOFs*NPad;
    su2double *solInt  = sources + nInt *NPad;
    su2double *fluxes  = solInt  + nInt *NPad;

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Interpolate the solution to the integration points of    ---*/
    /*---         the element.                                             ---*/
    /*------------------------------------------------------------------------*/

    /* Copy the solution of the DOFs into solDOFs for this chunk of elements. */
    for(unsigned short ll=0; ll<llEnd; ++ll) {

      /* Easier storage of the solution variables for this element. */
      const unsigned long lInd = l + ll;
      const unsigned short timeLevel = volElem[lInd].timeLevel;
      const su2double *solDOFsElem = VecWorkSolDOFs[timeLevel].data()
                                   + nVar*volElem[lInd].offsetDOFsSolThisTimeLevel;

      /* Loop over the DOFs and copy the data. */
      const unsigned short llNVar = ll*nVar;
      for(unsigned short i=0; i<nDOFs; ++i)
        for(unsigned short mm=0; mm<nVar; ++mm)
          solDOFs[i*NPad+llNVar+mm] = solDOFsElem[i*nVar+mm];
    }

    /* Call the general function to carry out the matrix product to determine
       the solution in the integration points of the chunk of elements. */
    blasFunctions->gemm(nInt, NPad, nDOFs, matBasisInt, solDOFs, solInt, config);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the inviscid fluxes, multiplied by minus the     ---*/
    /*---         integration weight, in the integration points.           ---*/
    /*---         If needed, also the source terms are computed.           ---*/
    /*------------------------------------------------------------------------*/

    /* Make a distinction between two and three space dimensions
       in order to have the most efficient code. */
    switch( nDim ) {

      case 2: {

        /* 2D simulation. Loop over the chunk of elements and loop over the
           integration points of the elements to compute the fluxes. */
        for(unsigned short ll=0; ll<llEnd; ++ll) {
          const unsigned short llNVar = ll*nVar;
          const unsigned long  lInd   = l + ll;

          for(unsigned short i=0; i<nInt; ++i) {
            const unsigned short iNPad = i*NPad;

            /* Easier storage of the metric terms and grid velocities
               in this integration point. */
            const su2double *metricTerms = volElem[lInd].metricTerms.data()
                                         + i*nMetricPerPoint;
            const su2double Jac          = metricTerms[0];
            const su2double *gridVel     = volElem[lInd].gridVelocities.data() + i*nDim;

            /* Compute the metric terms multiplied by minus the integration weight.
               The minus sign comes from the integration by parts in the weak
               formulation. */
            const su2double wDrdx = -weights[i]*metricTerms[1];
            const su2double wDrdy = -weights[i]*metricTerms[2];

            const su2double wDsdx = -weights[i]*metricTerms[3];
            const su2double wDsdy = -weights[i]*metricTerms[4];

            /* Easier storage of the location where the solution data of this
               integration point starts. */
            const su2double *sol = solInt + iNPad + llNVar;

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
            su2double *flux = fluxes + nDim*iNPad + llNVar;

            /*--- Fluxes in r-direction. */
            const su2double Ur = uRel*wDrdx + vRel*wDrdy;

            flux[0] = sol[0]*Ur;
            flux[1] = sol[1]*Ur + Pressure*wDrdx;
            flux[2] = sol[2]*Ur + Pressure*wDrdy;
            flux[3] = sol[3]*Ur + Pressure*(u*wDrdx + v*wDrdy);

            /*--- Fluxes in s-direction. */
            flux = flux + NPad;
            const su2double Us = uRel*wDsdx + vRel*wDsdy;

            flux[0] = sol[0]*Us;
            flux[1] = sol[1]*Us + Pressure*wDsdx;
            flux[2] = sol[2]*Us + Pressure*wDsdy;
            flux[3] = sol[3]*Us + Pressure*(u*wDsdx + v*wDsdy);

            /*--- If needed, compute the body forces in this integration point.
                  Note that the source terms are multiplied with minus the
                  integration weight in order to be consistent with the
                  formulation of the residual. Note that for the energy source
                  term the absolute velocity must be taken and not the
                  relative. ---*/
            if( body_force ) {
              su2double *source         = sources + iNPad + llNVar;
              const su2double weightJac = weights[i]*Jac;

              source[0] =  0.0;
              source[1] = -weightJac*body_force_vector[0];
              source[2] = -weightJac*body_force_vector[1];
              source[3] = -weightJac*(u*body_force_vector[0] + v*body_force_vector[1]);
            }
          }
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /* 3D simulation. Loop over the chunk of elements and loop over the
           integration points of the elements to compute the fluxes. */
        for(unsigned short ll=0; ll<llEnd; ++ll) {
          const unsigned short llNVar = ll*nVar;
          const unsigned long  lInd   = l + ll;
          for(unsigned short i=0; i<nInt; ++i) {
            const unsigned short iNPad = i*NPad;

            /* Easier storage of the metric terms and grid velocities
               in this integration point. */
            const su2double *metricTerms = volElem[lInd].metricTerms.data()
                                         + i*nMetricPerPoint;
            const su2double Jac          = metricTerms[0];
            const su2double *gridVel     = volElem[lInd].gridVelocities.data() + i*nDim;

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
            const su2double *sol = solInt + iNPad + llNVar;

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
            su2double *flux = fluxes + nDim*iNPad + llNVar;

            /*--- Fluxes in r-direction. */
            const su2double Ur = uRel*wDrdx + vRel*wDrdy + wRel*wDrdz;

            flux[0] = sol[0]*Ur;
            flux[1] = sol[1]*Ur + Pressure*wDrdx;
            flux[2] = sol[2]*Ur + Pressure*wDrdy;
            flux[3] = sol[3]*Ur + Pressure*wDrdz;
            flux[4] = sol[4]*Ur + Pressure*(u*wDrdx + v*wDrdy + w*wDrdz);

            /*--- Fluxes in s-direction. */
            flux = flux + NPad;
            const su2double Us = uRel*wDsdx + vRel*wDsdy + wRel*wDsdz;

            flux[0] = sol[0]*Us;
            flux[1] = sol[1]*Us + Pressure*wDsdx;
            flux[2] = sol[2]*Us + Pressure*wDsdy;
            flux[3] = sol[3]*Us + Pressure*wDsdz;
            flux[4] = sol[4]*Us + Pressure*(u*wDsdx + v*wDsdy + w*wDsdz);

            /*--- Fluxes in t-direction. */
            flux = flux + NPad;
            const su2double Ut = uRel*wDtdx + vRel*wDtdy + wRel*wDtdz;

            flux[0] = sol[0]*Ut;
            flux[1] = sol[1]*Ut + Pressure*wDtdx;
            flux[2] = sol[2]*Ut + Pressure*wDtdy;
            flux[3] = sol[3]*Ut + Pressure*wDtdz;
            flux[4] = sol[4]*Ut + Pressure*(u*wDtdx + v*wDtdy + w*wDtdz);

            /*--- If needed, compute the body forces in this integration point.
                  Note that the source terms are multiplied with minus the
                  integration weight in order to be consistent with the
                  formulation of the residual. Note that for the energy source
                  term the absolute velocity must be taken and not the
                  relative. ---*/
            if( body_force ) {
              su2double *source         = sources + iNPad + llNVar;
              const su2double weightJac = weights[i]*Jac;

              source[0] =  0.0;
              source[1] = -weightJac*body_force_vector[0];
              source[2] = -weightJac*body_force_vector[1];
              source[3] = -weightJac*body_force_vector[2];
              source[4] = -weightJac*(u*body_force_vector[0] + v*body_force_vector[1]
                        +             w*body_force_vector[2]);
            }
          }
        }

        break;
      }
    }

    /* Initialize addSourceTerms to body_force. The value of addSourceTerms
       is set to true when a manufactured solution is computed. */
    bool addSourceTerms = body_force;

    /* Check whether or not a manufactured solution is used. */
    if( VerificationSolution ) {
      if( VerificationSolution->IsManufacturedSolution() ) {

        /*--- For the manufactured solutions a source term must be added. If a
              standard source term has not been specified, initialize the source
              terms to zero and set addSourceTerms to true. ---*/
        addSourceTerms = true;
        if( !body_force ) {
          for(unsigned short i=0; i<(nInt*NPad); ++i)
            sources[i] = 0.0;
        }

        /*--- Loop over the chunk of elements and its integration points. ---*/
        for(unsigned short ll=0; ll<llEnd; ++ll) {
          const unsigned short llNVar = ll*nVar;
          const unsigned long  lInd   = l + ll;
          for(unsigned short i=0; i<nInt; ++i) {
            const unsigned short iNPad = i*NPad;

            /* Determine the integration weight multiplied by the Jacobian. */
            const su2double *metricTerms = volElem[lInd].metricTerms.data()
                                         + i*nMetricPerPoint;
            const su2double weightJac    = weights[i]*metricTerms[0];

            /* Set the pointer to the coordinates in this integration point and
               call the function to compute the source terms for the manufactured
               solution. */
            const su2double *coor = volElem[lInd].coorIntegrationPoints.data() + i*nDim;

            su2double sourceMan[5];

            VerificationSolution->GetMMSSourceTerm(coor, time, sourceMan);

            /*--- Subtract the source term of the manufactured solution, multiplied
                  by the appropriate weight, from the possibly earlier computed
                  source term. It is subtracted in order to be consistent with
                  the definition of the residual used in this code. ---*/
            su2double *source = sources + iNPad + llNVar;
            for(unsigned short k=0; k<nVar; ++k)
              source[k] -= weightJac*sourceMan[k];
          }
        }
      }
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
    /*---         integration over the volume element.                     ---*/
    /*------------------------------------------------------------------------*/

    /* Call the general function to carry out the matrix product.
       Use solDOFs as a temporary storage for the matrix product. */
    blasFunctions->gemm(nDOFs, NPad, nInt*nDim, matDerBasisIntTrans, fluxes, solDOFs, config);

    /* Add the contribution from the source terms, if needed. Use solInt
       as temporary storage for the matrix product. */
    if( addSourceTerms ) {

      /* Call the general function to carry out the matrix product. */
      blasFunctions->gemm(nDOFs, NPad, nInt, matBasisIntTrans, sources, solInt, config);

      /* Add the residuals due to source terms to the volume residuals */
      for(unsigned short i=0; i<(nDOFs*NPad); ++i)
        solDOFs[i] += solInt[i];
    }

    /* Loop over the elements in this chunk to store the residuals
       in the appropriate locations. */
    for(unsigned short ll=0; ll<llEnd; ++ll) {
      const unsigned short llNVar = ll*nVar;
      const unsigned long  lInd   = l + ll;

      /* Easier storage of the residuals for this volume element. */
      su2double *res = VecResDOFs.data() + nVar*volElem[lInd].offsetDOFsSolLocal;

      /* Loop over the DOFs and copy the data. */
      for(unsigned short i=0; i<nDOFs; ++i)
        for(unsigned short mm=0; mm<nVar; ++mm)
          res[i*nVar+mm] = solDOFs[i*NPad+llNVar+mm];
    }

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::Boundary_Conditions(const unsigned short timeLevel,
                                              CConfig              *config,
                                              CNumerics            **numerics,
                                              const bool           haloInfoNeededForBC,
                                              su2double            *workArray){

  /* Loop over all boundaries. */
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /* Check if this boundary marker must be treated at all. */
    if(boundaries[iMarker].haloInfoNeededForBC == haloInfoNeededForBC) {

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
                          numerics[CONV_BOUND_TERM], workArray);
            break;
          case FAR_FIELD:
            BC_Far_Field(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                         numerics[CONV_BOUND_TERM], workArray);
            break;
          case SYMMETRY_PLANE:
            BC_Sym_Plane(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                         numerics[CONV_BOUND_TERM], workArray);
            break;
          case SUPERSONIC_INLET: /* Use far field for this. When a more detailed state
                                    needs to be specified, use a Riemann boundary. */
            BC_Far_Field(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                         numerics[CONV_BOUND_TERM], workArray);
            break;
          case SUPERSONIC_OUTLET:
            BC_Supersonic_Outlet(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                                 numerics[CONV_BOUND_TERM], workArray);
            break;
          case INLET_FLOW:
            BC_Inlet(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                     numerics[CONV_BOUND_TERM], iMarker, workArray);
            break;
          case OUTLET_FLOW:
            BC_Outlet(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                      numerics[CONV_BOUND_TERM], iMarker, workArray);
            break;
          case ISOTHERMAL:
            BC_Isothermal_Wall(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                               numerics[CONV_BOUND_TERM], iMarker, workArray);
            break;
          case HEAT_FLUX:
            BC_HeatFlux_Wall(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                             numerics[CONV_BOUND_TERM], iMarker, workArray);
            break;
          case RIEMANN_BOUNDARY:
            BC_Riemann(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                       numerics[CONV_BOUND_TERM], iMarker, workArray);
            break;
          case CUSTOM_BOUNDARY:
            BC_Custom(config, surfElemBeg, surfElemEnd, surfElem, resFaces,
                      numerics[CONV_BOUND_TERM], workArray);
            break;
          case PERIODIC_BOUNDARY:  // Nothing to be done for a periodic boundary.
            break;
          default:
            SU2_MPI::Error("BC not implemented.", CURRENT_FUNCTION);
        }
      }
    }
  }
}

void CFEM_DG_EulerSolver::ResidualFaces(CConfig             *config,
                                        const unsigned long indFaceBeg,
                                        const unsigned long indFaceEnd,
                                        unsigned long       &indResFaces,
                                        CNumerics           *numerics,
                                        su2double           *workArray) {

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of matching faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=indFaceBeg; l<indFaceEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(matchingInternalFaces, l, indFaceEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardMatchingFacesSol[ind].GetNIntegration();
    const su2double *weights  = standardMatchingFacesSol[ind].GetWeightsIntegration();

    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*max(nInt, nDOFsFace0);
    su2double *fluxes  = solIntR + NPad*max(nInt, nDOFsFace1);

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Compute the inviscid fluxes in the integration points of ---*/
    /*---         this matching face multiplied by the integration weight. ---*/
    /*------------------------------------------------------------------------*/

    /* Compute the inviscid fluxes in the integration points. */
    InviscidFluxesInternalMatchingFace(config, l, lEnd, NPad,
                                       solIntL, solIntR, fluxes, numerics);

    /* Multiply the fluxes with the integration weight of the corresponding
       integration point. */
    for(unsigned short i=0; i<nInt; ++i) {
      su2double *flux = fluxes + i*NPad;

      for(unsigned short j=0; j<NPad; ++j)
        flux[j] *= weights[i];
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the contribution to the residuals from the       ---*/
    /*---         integration over this chunk of internal matching faces.  ---*/
    /*------------------------------------------------------------------------*/

    /* Get the correct form of the basis functions needed for the matrix
       multiplication to compute the residual. Use solIntL as a temporary
       buffer to store this product. */
    const su2double *basisFaceTrans = standardMatchingFacesSol[ind].GetBasisFaceIntegrationTransposeSide0();
    su2double *resSide0 = solIntL;

    /* Call the general function to carry out the matrix product. */
    blasFunctions->gemm(nDOFsFace0, NPad, nInt, basisFaceTrans, fluxes, resSide0, config);

    /* Check if the number of DOFs on both sides of the face is different.
       In that case also the matrix product with the basis functions on side 1
       must be carried out. Use solIntR as a temporary buffer to store
       this product. */
    su2double *resSide1 = solIntR;
    if(nDOFsFace1 != nDOFsFace0) {
      basisFaceTrans = standardMatchingFacesSol[ind].GetBasisFaceIntegrationTransposeSide1();
      blasFunctions->gemm(nDOFsFace1, NPad, nInt, basisFaceTrans, fluxes, resSide1, config);
    }

    /* Loop over the number of faces in this chunk. */
    for(unsigned short ll=0; ll<llEnd; ++ll) {

      /* Set the pointer in the residual array where the residual of side 0
         and side 1 of this face must be stored. Update the counter. */
      su2double *resFace0 = VecResFaces.data() + indResFaces*nVar;
      indResFaces        += nDOFsFace0;
      su2double *resFace1 = VecResFaces.data() + indResFaces*nVar;
      indResFaces        += nDOFsFace1;

      /* Loop over the DOFs and copy the data for side 0. */
      for(unsigned short i=0; i<nDOFsFace0; ++i)
        for(unsigned short mm=0; mm<nVar; ++mm)
          resFace0[nVar*i+mm] = resSide0[NPad*i+nVar*ll+mm];

      /* If the number of DOFs on both sides is the same, then the residual
         of side 1 is obtained by simply negating the residual of side 0.
         Otherwise a copy must be carried out and the residual negated,
         because the normal is pointing into the adjacent element for side 1. */
      if(nDOFsFace1 == nDOFsFace0) {
        for(unsigned short i=0; i<(nVar*nDOFsFace1); ++i)
          resFace1[i] = -resFace0[i];
      }
      else {
        for(unsigned short i=0; i<nDOFsFace1; ++i)
          for(unsigned short mm=0; mm<nVar; ++mm)
            resFace1[nVar*i+mm] = -resSide1[NPad*i+nVar*ll+mm];
      }
    }

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::InviscidFluxesInternalMatchingFace(
                                              CConfig              *config,
                                              const unsigned long  lBeg,
                                              const unsigned long  lEnd,
                                              const unsigned short NPad,
                                              su2double            *solIntL,
                                              su2double            *solIntR,
                                              su2double            *fluxes,
                                              CNumerics            *numerics) {

  /* Set the pointer solFace to fluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the face and the fluxes. */
  su2double *solFace = fluxes;

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Interpolate the left state in the integration points of  ---*/
  /*---         the face.                                                ---*/
  /*------------------------------------------------------------------------*/

  /* Get the required information from the corresponding standard face,
     which is the same for the chunk of faces considered. */
  const unsigned short ind        = matchingInternalFaces[lBeg].indStandardElement;
  const unsigned short nInt       = standardMatchingFacesSol[ind].GetNIntegration();
  const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
  const su2double     *basisFace0 = standardMatchingFacesSol[ind].GetBasisFaceIntegrationSide0();

  /* Determine the number of faces treated simultaneously and define the arrays
     of pointers to be used later on in the call to ComputeInviscidFluxesFace. */
  const unsigned short llEnd = lEnd - lBeg;
  const su2double* arrNorm[SIZE_ARR_NORM];
  const su2double* arrGridVel[SIZE_ARR_NORM];

  if(llEnd > SIZE_ARR_NORM)
    SU2_MPI::Error("SIZE_ARR_NORM is too small. Increase it or decrease ALIGNED_BYTES_MATMUL",
                   CURRENT_FUNCTION);

  /*--- Loop over the chunk of faces to copy the states in the DOFs
        of the faces in the correct sequence. ---*/
  for(unsigned long l=lBeg; l<lEnd; ++l) {
    const unsigned short ll     = l-lBeg;
    const unsigned short llNVar = ll*nVar;

    /* Store the pointer to the normal and grid velocity information of this
       face. This is used later on in the call to ComputeInviscidFluxesFace. */
    arrNorm[ll]    = matchingInternalFaces[l].metricNormalsFace.data();
    arrGridVel[ll] = matchingInternalFaces[l].gridVelocities.data();

    /* Determine the time level of the face. */
    const unsigned long  elem0         = matchingInternalFaces[l].elemID0;
    const unsigned long  elem1         = matchingInternalFaces[l].elemID1;
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

    /*--- Store the solution of the DOFs of side 0 of the face in contiguous
          memory, such that the function blasFunctions->gemm can be used to
          compute the left states in the integration points. ---*/
    const unsigned long *DOFs = matchingInternalFaces[l].DOFsSolFaceSide0.data();
    for(unsigned short i=0; i<nDOFsFace0; ++i) {
      const su2double *solDOF = VecWorkSolDOFs[timeLevelFace].data()
                              + nVar*(DOFs[i] - offset);
      su2double       *sol    = solFace + NPad*i + llNVar;
      for(unsigned short mm=0; mm<nVar; ++mm)
        sol[mm] = solDOF[mm];
    }
  }

  /* Compute the left states. Call the general function to
     carry out the matrix product. */
  blasFunctions->gemm(nInt, NPad, nDOFsFace0, basisFace0, solFace, solIntL, config);

  /*------------------------------------------------------------------------*/
  /*--- Step 2: Interpolate the right state in the integration points of ---*/
  /*---         the face.                                                ---*/
  /*------------------------------------------------------------------------*/

  /* Get the required information from the corresponding standard face,
     which is the same for the chunk of faces considered. */
  const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();
  const su2double     *basisFace1 = standardMatchingFacesSol[ind].GetBasisFaceIntegrationSide1();

  /*--- Loop over the chunk of elements to copy the states in the DOFs
        of the faces in the correct sequence. ---*/
  for(unsigned long l=lBeg; l<lEnd; ++l) {
    const unsigned short ll     = l-lBeg;
    const unsigned short llNVar = ll*nVar;

    /* Determine the time level of the face. */
    const unsigned long  elem0         = matchingInternalFaces[l].elemID0;
    const unsigned long  elem1         = matchingInternalFaces[l].elemID1;
    const unsigned short timeLevelFace = min(volElem[elem0].timeLevel,
                                             volElem[elem1].timeLevel);

    /* Determine the offset between the numbering used in the DOFs of the face
       and the numbering used in the working solution vectors, see above for a
       more detailed explanation. */
    unsigned long offset;
    if(timeLevelFace == volElem[elem1].timeLevel)
      offset = volElem[elem1].offsetDOFsSolLocal - volElem[elem1].offsetDOFsSolThisTimeLevel;
    else
      offset = volElem[elem1].offsetDOFsSolLocal - volElem[elem1].offsetDOFsSolPrevTimeLevel;

    /*--- Store the solution of the DOFs of side 1 of the face in contiguous
          memory, such that the function blasFunctions->gemm can be used to
          compute the left states in the integration points. ---*/
    const unsigned long *DOFs = matchingInternalFaces[l].DOFsSolFaceSide1.data();
    for(unsigned short i=0; i<nDOFsFace1; ++i) {
      const su2double *solDOF = VecWorkSolDOFs[timeLevelFace].data()
                              + nVar*(DOFs[i] - offset);
      su2double       *sol    = solFace + NPad*i + llNVar;
      for(unsigned short mm=0; mm<nVar; ++mm)
        sol[mm] = solDOF[mm];
    }
  }

  /* Compute the right states. Call the general function to
     carry out the matrix product. */
  blasFunctions->gemm(nInt, NPad, nDOFsFace1, basisFace1, solFace, solIntR, config);

  /*------------------------------------------------------------------------*/
  /*--- Step 3: Compute the fluxes in the integration points using the   ---*/
  /*---         approximate Riemann solver.                              ---*/
  /*------------------------------------------------------------------------*/

  /* General function to compute the fluxes in the integration points
     of the faces. */
  ComputeInviscidFluxesFace(config, llEnd, NPad, nInt, arrNorm, arrGridVel,
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
                                              const unsigned long elemEnd,
                                              su2double           *workArray) {

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
    if( !volElem[l].lumpedMassMatrix.empty() ) {

      /* Multiply the residual with the inverse of the lumped mass matrix. */
      for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {

        su2double *resDOF = res + nVar*i;
        su2double lMInv   = 1.0/volElem[l].lumpedMassMatrix[i];
        for(unsigned short k=0; k<nVar; ++k) resDOF[k] *= lMInv;
      }
    }
    else {

      /* Multiply the residual with the inverse of the mass matrix.
         Use the array workArray as temporary storage. */
      for(unsigned short mm=0; mm<nVar*volElem[l].nDOFsSol; ++mm)
        workArray[mm] = res[mm];

      blasFunctions->gemm(volElem[l].nDOFsSol, nVar, volElem[l].nDOFsSol,
                          volElem[l].invMassMatrix.data(), workArray, res, config);
    }
  }
}

void CFEM_DG_EulerSolver::Pressure_Forces(const CGeometry* geometry, const CConfig* config) {

  /* Allocate the memory for the work array and initialize it to zero to avoid
     warnings in debug mode  about uninitialized memory when padding is applied. */
  vector<su2double> workArrayVec(sizeWorkArray, 0.0);
  su2double *workArray = workArrayVec.data();

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Get the information of the angle of attack, reference area, etc. ---*/
  const su2double Alpha        = config->GetAoA()*PI_NUMBER/180.0;
  const su2double Beta         = config->GetAoS()*PI_NUMBER/180.0;
  const su2double RefArea      = config->GetRefArea();
  const su2double RefLength    = config->GetRefLength();
  const su2double Gas_Constant = config->GetGas_ConstantND();
  auto Origin                  = config->GetRefOriginMoment(0);
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

  AeroCoeffForceRef = 0.5 * RefDensity * RefArea * RefVel2;
  const su2double factor = 1.0 / AeroCoeffForceRef;

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

        /*--- Loop over the faces of this boundary. Multiple faces are treated
              simultaneously to improve the performance of the matrix
              multiplications. As a consequence, the update of the counter l
              happens at the end of this loop section. ---*/
        for(unsigned long l=0; l<nSurfElem;) {

          /* Determine the end index for this chunk of faces and the padded
             N value in the gemm computations. */
          unsigned long lEnd;
          unsigned short ind, llEnd, NPad;

          MetaDataChunkOfElem(surfElem, l, nSurfElem, nFaceSimul,
                              nPadMin, lEnd, ind, llEnd, NPad);

          /* Get the required information from the corresponding standard face. */
          const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
          const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();
          const su2double *basisFace = standardBoundaryFacesSol[ind].GetBasisFaceIntegration();
          const su2double *weights   = standardBoundaryFacesSol[ind].GetWeightsIntegration();

          /* Set the pointers for the local work arrays. */
          su2double *solInt    = workArray;
          su2double *workArray = solInt + NPad*nInt;

          /* Loop over the faces that are treated simultaneously. */
          for(unsigned short ll=0; ll<llEnd; ++ll) {
            const unsigned short llNVar = ll*nVar;

            /* Easier storage of the DOFs of the face. */
            const unsigned long *DOFs = surfElem[l+ll].DOFsSolFace.data();

            /* Copy the solution of the DOFs of the face such that it is
               contiguous in memory. */
            for(unsigned short i=0; i<nDOFs; ++i) {
              const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
              su2double       *sol    = workArray + NPad*i + llNVar;
              for(unsigned short mm=0; mm<nVar; ++mm)
                sol[mm] = solDOF[mm];
            }
          }

          /* Call the general function to carry out the matrix product to determine
             the solution in the integration points. */
          blasFunctions->gemm(nInt, NPad, nDOFs, basisFace, workArray, solInt, config);

          /* Make a distinction between two and three space dimensions
             in order to have the most efficient code. */
          switch( nDim ) {

            case 2: {

              /* Two dimensional simulation. Loop over the number of faces treated
                 simultaneously. */
              for(unsigned short ll=0; ll<llEnd; ++ll) {
                const unsigned short llNVar = ll*nVar;
                const unsigned long  lll    = l + ll;

                /* Loop over the integration points of this surface element. */
                for(unsigned short i=0; i<nInt; ++i) {

                  /* Easier storage of the solution, the normals and the coordinates
                     for this integration point. */
                  const su2double *sol     = solInt + i*NPad + llNVar;
                  const su2double *normals = surfElem[lll].metricNormalsFace.data()
                                           + i*(nDim+1);
                  const su2double *Coord   = surfElem[lll].coorIntegrationPoints.data()
                                           + i*nDim;

                  /*--- Compute the pressure in this integration point. ---*/
                  const su2double DensityInv   = 1.0/sol[0];
                  const su2double u            = sol[1]*DensityInv;
                  const su2double v            = sol[2]*DensityInv;
                  const su2double StaticEnergy = sol[3]*DensityInv - 0.5*(u*u + v*v);

                  FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
                  const su2double Pressure = FluidModel->GetPressure();

                  /*-- Compute the vector from the reference point to the integration
                       point and update the inviscid force. Note that the normal points
                       into the geometry, hence no minus sign. ---*/
                  const su2double DistX = Coord[0] - Origin[0];
                  const su2double DistY = Coord[1] - Origin[1];

                  const su2double ForceMag = (Pressure - Pressure_Inf)*weights[i]
                                           * normals[nDim]*factor;
                  const su2double ForceX   = ForceMag*normals[0];
                  const su2double ForceY   = ForceMag*normals[1];

                  ForceInviscid[0] += ForceX;
                  ForceInviscid[1] += ForceY;

                  /*--- Update the inviscid moment. ---*/
                  MomentInviscid[2] += (ForceY*DistX - ForceX*DistY)/RefLength;
                }
              }

              break;
            }

            /*----------------------------------------------------------------*/

            case 3: {

              /* Three dimensional simulation. Loop over the number of faces treated
                 simultaneously. */
              for(unsigned short ll=0; ll<llEnd; ++ll) {
                const unsigned short llNVar = ll*nVar;
                const unsigned long  lll    = l + ll;

                /* Loop over the integration points of this surface element. */
                for(unsigned short i=0; i<nInt; ++i) {

                  /* Easier storage of the solution, the normals and the coordinates
                     for this integration point. */
                  const su2double *sol     = solInt + i*NPad + llNVar;
                  const su2double *normals = surfElem[lll].metricNormalsFace.data()
                                           + i*(nDim+1);
                  const su2double *Coord   = surfElem[lll].coorIntegrationPoints.data()
                                           + i*nDim;

                  /*--- Compute the pressure in this integration point. ---*/
                  const su2double DensityInv   = 1.0/sol[0];
                  const su2double u            = sol[1]*DensityInv;
                  const su2double v            = sol[2]*DensityInv;
                  const su2double w            = sol[3]*DensityInv;
                  const su2double StaticEnergy = sol[4]*DensityInv - 0.5*(u*u + v*v + w*w);

                  FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
                  const su2double Pressure = FluidModel->GetPressure();

                  /*-- Compute the vector from the reference point to the integration
                       point and update the inviscid force. Note that the normal points
                       into the geometry, hence no minus sign. ---*/
                  const su2double DistX = Coord[0] - Origin[0];
                  const su2double DistY = Coord[1] - Origin[1];
                  const su2double DistZ = Coord[2] - Origin[2];

                  const su2double ForceMag = (Pressure - Pressure_Inf)*weights[i]
                                           * normals[nDim]*factor;
                  const su2double ForceX   = ForceMag*normals[0];
                  const su2double ForceY   = ForceMag*normals[1];
                  const su2double ForceZ   = ForceMag*normals[2];

                  ForceInviscid[0] += ForceX;
                  ForceInviscid[1] += ForceY;
                  ForceInviscid[2] += ForceZ;

                  /*--- Update the inviscid moment. ---*/
                  MomentInviscid[0] += (ForceZ*DistY - ForceY*DistZ)/RefLength;
                  MomentInviscid[1] += (ForceX*DistZ - ForceZ*DistX)/RefLength;
                  MomentInviscid[2] += (ForceY*DistX - ForceX*DistY)/RefLength;
                }
              }

              break;
            }
          }

          /* Update the value of the counter l to the end index of the
             current chunk. */
          l = lEnd;
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
  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(locBuf.data(), globBuf.data(), nCommSize, MPI_DOUBLE,
                       MPI_SUM, SU2_MPI::GetComm());
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

  SetResToZero();

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

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar*volElem[l].offsetDOFsSolLocal;
    const su2double *res        = VecResDOFs.data() + offset;
    const su2double *solDOFsOld = VecSolDOFs.data() + offset;
    su2double *solDOFs          = solNew            + offset;

    /* Loop over the DOFs for this element and update the solution. */
    const su2double tmp = RK_AlphaCoeff*VecDeltaTime[l];
    const unsigned short nVarNDOFs = nVar*volElem[l].nDOFsSol;

    for(unsigned short j=0; j<nVarNDOFs; ++j)
      solDOFs[j] = solDOFsOld[j] - tmp*res[j];
  }

  /*--- Test for the last RK step. ---*/
  if(iRKStep == (nRKStages-1)) {

    /*--- Compute the root mean square residual. Note that the SetResidual_RMS
          function of CSolver cannot be used, because that is for the FV solver. ---*/
    SetResidual_RMS_FEM(geometry, config);

    /*--- For verification cases, compute the global error metrics. ---*/
    ComputeVerificationError(geometry, config);
  }
}

void CFEM_DG_EulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                                 CConfig *config, unsigned short iRKStep) {

  /*--- Hard-coded classical RK4 coefficients. Will be added to config. ---*/
  su2double RK_FuncCoeff[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
  su2double RK_TimeCoeff[4] = {0.5, 0.5, 1.0, 1.0};

  SetResToZero();

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

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

    const unsigned short nVarNDOFs = nVar*volElem[l].nDOFsSol;

    if (iRKStep < 3) {
      for(unsigned short j=0; j<nVarNDOFs; ++j) {
        solDOFsNew[j] += tmp_func*res[j];
        solDOFs[j]     = solDOFsOld[j] + tmp_time*res[j];
      }
    }
    else {
      for(unsigned short j=0; j<nVarNDOFs; ++j)
        solDOFs[j] = solDOFsNew[j] + tmp_func*res[j];
    }
  }

  /*--- Test for the last RK step. ---*/
  if(iRKStep == 3) {

    /*--- Compute the root mean square residual. Note that the SetResidual_RMS
          function of CSolver cannot be used, because that is for the FV solver. ---*/
    SetResidual_RMS_FEM(geometry, config);

    /*--- For verification cases, compute the global error metrics. ---*/
    ComputeVerificationError(geometry, config);
  }
}

void CFEM_DG_EulerSolver::SetResidual_RMS_FEM(CGeometry *geometry,
                                              CConfig *config) {

  /* Initialize the residuals to zero. */
  SetResToZero();

  /*--- Loop over the owned elements. It is not possible to loop directly over the owned
        DOFs, because the coordinates of the DOFs are only known in the volume class. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Set the pointer for the residual for this element. */
    const unsigned long offset  = nVar*volElem[l].offsetDOFsSolLocal;
    const su2double *res        = VecResDOFs.data() + offset;

    /* Loop over the DOFs for this element and update the norms. */
    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      const su2double *coor = volElem[l].coorSolDOFs.data() + j*nDim;

      for(unsigned short iVar=0; iVar<nVar; ++iVar, ++i) {

        Residual_RMS[iVar] += res[i]*res[i];
        AddRes_Max(iVar, fabs(res[i]), globalIndex, coor);
      }
    }
  }

#ifdef HAVE_MPI
  /* Parallel mode. Disable the reduce for the residual to avoid overhead if requested. */
  if (config->GetComm_Level() == COMM_FULL) {

    /*--- The local L2 norms must be added to obtain the
          global value. Also check for divergence. ---*/
    vector<su2double> rbufRes(nVar);
    SU2_MPI::Allreduce(Residual_RMS.data(), rbufRes.data(), nVar, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

    for(unsigned short iVar=0; iVar<nVar; ++iVar) {
      if (rbufRes[iVar] != rbufRes[iVar])
        SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

      Residual_RMS[iVar] = max(EPS*EPS, sqrt(rbufRes[iVar]/nDOFsGlobal));
    }

    /*--- The global maximum norms must be obtained. ---*/
    rbufRes.resize(nVar*size);
    SU2_MPI::Allgather(Residual_Max.data(), nVar, MPI_DOUBLE, rbufRes.data(),
                       nVar, MPI_DOUBLE, SU2_MPI::GetComm());

    vector<unsigned long> rbufPoint(nVar*size);
    SU2_MPI::Allgather(Point_Max.data(), nVar, MPI_UNSIGNED_LONG, rbufPoint.data(),
                       nVar, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

    vector<su2double> sbufCoor(nDim*nVar);
    for(unsigned short iVar=0; iVar<nVar; ++iVar) {
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        sbufCoor[iVar*nDim+iDim] = Point_Max_Coord[iVar][iDim];
    }

    vector<su2double> rbufCoor(nDim*nVar*size);
    SU2_MPI::Allgather(sbufCoor.data(), nVar*nDim, MPI_DOUBLE, rbufCoor.data(),
                       nVar*nDim, MPI_DOUBLE, SU2_MPI::GetComm());

    for(unsigned short iVar=0; iVar<nVar; ++iVar) {
      for(int proc=0; proc<size; ++proc)
        AddRes_Max(iVar, rbufRes[proc*nVar+iVar], rbufPoint[proc*nVar+iVar],
                   &rbufCoor[proc*nVar*nDim+iVar*nDim]);
    }
  }

#else
  /*--- Sequential mode. Check for a divergence of the solver and compute
   the L2-norm of the residuals. ---*/
  for(unsigned short iVar=0; iVar<nVar; ++iVar) {

    if(GetRes_RMS(iVar) != GetRes_RMS(iVar))
      SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

    Residual_RMS[iVar] = max(EPS*EPS, sqrt(GetRes_RMS(iVar)/nDOFsGlobal));
  }

#endif
}

void CFEM_DG_EulerSolver::ComputeVerificationError(CGeometry *geometry,
                                                   CConfig   *config) {

  /*--- The errors only need to be computed on the finest grid. ---*/
  if(MGLevel != MESH_0) return;

  /*--- If this is a verification case, we can compute the global
   error metrics by using the difference between the local error
   and the known solution at each DOF. This is then collected into
   RMS (L2) and maximum (Linf) global error norms. From these
   global measures, one can compute the order of accuracy. ---*/

  bool write_heads = ((((config->GetTimeIter() % (config->GetScreen_Wrt_Freq(2)*40)) == 0)
                       && (config->GetTimeIter()!= 0))
                      || (config->GetTimeIter() == 1));
  if( !write_heads ) return;

  /*--- Check if there actually is an exact solution for this
        verification case, if computed at all. ---*/
  if (VerificationSolution) {
    if (VerificationSolution->ExactSolutionKnown()) {

      /*--- Get the physical time if necessary. ---*/
      su2double time = 0.0;
      if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

      /*--- Reset the global error measures to zero. ---*/
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        VerificationSolution->SetError_RMS(iVar, 0.0);
        VerificationSolution->SetError_Max(iVar, 0.0, 0);
      }

      /*--- Loop over the owned elements. It is not possible to loop directly over the owned
            DOFs, because the coordinates of the DOFs are only known in the volume element class. ---*/
      for(unsigned long l=0; l<nVolElemOwned; ++l) {

        /* Set the pointer for the solution for this element. */
        const unsigned long offset = nVar*volElem[l].offsetDOFsSolLocal;
        const su2double *solDOFs      = VecSolDOFs.data() + offset;

        /* Loop over the DOFs for this element. */
        for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
          const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
          const su2double *solDOF = solDOFs + j*nVar;
          su2double *coor = volElem[l].coorSolDOFs.data() + j*nDim;

          /* Get local error from the verification solution class. */
          vector<su2double> error(nVar,0.0);
          VerificationSolution->GetLocalError(coor, time, solDOF, error.data());

          /* Increment the global error measures */
          for (unsigned short iVar = 0; iVar < nVar; iVar++) {
            VerificationSolution->AddError_RMS(iVar, error[iVar]*error[iVar]);
            VerificationSolution->AddError_Max(iVar, fabs(error[iVar]), globalIndex, coor);
          }
        }
      }

      /* Finalize the calculation of the global error measures. */
      VerificationSolution->SetVerificationError(nDOFsGlobal, config);

      /*--- Screen output of the error metrics. This can be improved
       once the new output classes are in place. ---*/

      if (rank == MASTER_NODE) {

        cout.precision(5);
        cout.setf(ios::scientific, ios::floatfield);

        if (!config->GetDiscrete_Adjoint()) {

          cout << endl   << "------------------------ Global Error Analysis --------------------------" << endl;

          cout << setw(20) << "RMS Error  [Rho]: " << setw(12) << VerificationSolution->GetError_RMS(0) << "     | ";
          cout << setw(20) << "Max Error  [Rho]: " << setw(12) << VerificationSolution->GetError_Max(0);
          cout << endl;

          cout << setw(20) << "RMS Error [RhoU]: " << setw(12) << VerificationSolution->GetError_RMS(1) << "     | ";
          cout << setw(20) << "Max Error [RhoU]: " << setw(12) << VerificationSolution->GetError_Max(1);
          cout << endl;

          cout << setw(20) << "RMS Error [RhoV]: " << setw(12) << VerificationSolution->GetError_RMS(2) << "     | ";
          cout << setw(20) << "Max Error [RhoV]: " << setw(12) << VerificationSolution->GetError_Max(2);
          cout << endl;

          if (nDim == 3) {
            cout << setw(20) << "RMS Error [RhoW]: " << setw(12) << VerificationSolution->GetError_RMS(3) << "     | ";
            cout << setw(20) << "Max Error [RhoW]: " << setw(12) << VerificationSolution->GetError_Max(3);
            cout << endl;
          }

          cout << setw(20) << "RMS Error [RhoE]: " << setw(12) << VerificationSolution->GetError_RMS(nDim+1) << "     | ";
          cout << setw(20) << "Max Error [RhoE]: " << setw(12) << VerificationSolution->GetError_Max(nDim+1);
          cout << endl;

          cout << "-------------------------------------------------------------------------" << endl << endl;
          cout.unsetf(ios_base::floatfield);

        }
      }
    }
  }
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
                                                    const unsigned short     nFaceSimul,
                                                    const unsigned short     NPad,
                                                    const CSurfaceElementFEM *surfElem,
                                                    const su2double          *solIntL,
                                                    su2double                *solIntR) {

  /*--- Apply the inviscid wall boundary conditions to compute the right
        state in the integration points. There are two options. Either the
        normal velocity is negated or the normal velocity is set to zero,
        relative to the grid motion. Some experiments are needed to see which
        formulation gives better results. ---*/

  /* Get the required information from the standard element, which is the
     same for all simultaneously treated faces. */
  const unsigned short ind  = surfElem[0].indStandardElement;
  const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

  /* Loop over the faces that are treated simultaneously. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {
    const unsigned short llNVar = l*nVar;

    /* Loop over the integration points of the face. */
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution, the normals and the
         grid velocity for this integration point. */
      const su2double *UL      = solIntL + i*NPad + llNVar;
            su2double *UR      = solIntR + i*NPad + llNVar;
      const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);
      const su2double *gridVel = surfElem[l].gridVelocities.data() + i*nDim;

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
}

void CFEM_DG_EulerSolver::BoundaryStates_Inlet(CConfig                  *config,
                                               const unsigned short     nFaceSimul,
                                               const unsigned short     NPad,
                                               const CSurfaceElementFEM *surfElem,
                                               unsigned short           val_marker,
                                               const su2double          *solIntL,
                                               su2double                *solIntR) {

  /*--- Retrieve the specified total conditions for this inlet. ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double P_Total   = config->GetInlet_Ptotal(Marker_Tag);
  su2double T_Total   = config->GetInlet_Ttotal(Marker_Tag);
  auto Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

  /*--- Non-dim. the inputs if necessary, and compute the total enthalpy. ---*/
  P_Total /= config->GetPressure_Ref();
  T_Total /= config->GetTemperature_Ref();

  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double H_Total      = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;

  /* Get the required information from the standard element, which is the
     same for all simultaneously treated faces. */
  const unsigned short ind  = surfElem[0].indStandardElement;
  const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

  /* Loop over the faces that are treated simultaneously. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {
    const unsigned short llNVar = l*nVar;

    /* Loop over the integration points of the face. */
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution and the normals
         for this integration point. */
      const su2double *UL      = solIntL + i*NPad + llNVar;
            su2double *UR      = solIntR + i*NPad + llNVar;
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
  }
}

void CFEM_DG_EulerSolver::BoundaryStates_Outlet(CConfig                  *config,
                                                const unsigned short     nFaceSimul,
                                                const unsigned short     NPad,
                                                const CSurfaceElementFEM *surfElem,
                                                unsigned short           val_marker,
                                                const su2double          *solIntL,
                                                su2double                *solIntR) {

  /*--- Retrieve the specified back pressure for this outlet.
        Nondimensionalize, if necessary. ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double P_Exit = config->GetOutlet_Pressure(Marker_Tag);
  P_Exit = P_Exit/config->GetPressure_Ref();

  /* Get the required information from the standard element, which is the
     same for all simultaneously treated faces. */
  const unsigned short ind  = surfElem[0].indStandardElement;
  const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

  /* Loop over the faces that are treated simultaneously. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {
    const unsigned short llNVar = l*nVar;

    /* Loop over the integration points of the face. */
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the left and right solution and the normals
         for this integration point. */
      const su2double *UL      = solIntL + i*NPad + llNVar;
            su2double *UR      = solIntR + i*NPad + llNVar;
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
  }
}

void CFEM_DG_EulerSolver::BoundaryStates_Riemann(CConfig                  *config,
                                                 const unsigned short     nFaceSimul,
                                                 const unsigned short     NPad,
                                                 const CSurfaceElementFEM *surfElem,
                                                 unsigned short           val_marker,
                                                 const su2double          *solIntL,
                                                 su2double                *solIntR) {

  /* Retrieve the corresponding string for this marker. */
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /* Determine the number of integration points for the face, which is
     the same for all simultaneously treated faces. */
  const unsigned short ind  = surfElem[0].indStandardElement;
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
      auto Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);

      P_Total /= config->GetPressure_Ref();
      T_Total /= config->GetTemperature_Ref();

      /* Compute the total enthalpy and entropy from these values. */
      FluidModel->SetTDState_PT(P_Total, T_Total);

      const su2double Enthalpy_e = FluidModel->GetStaticEnergy()
                                 + FluidModel->GetPressure()/FluidModel->GetDensity();
      const su2double Entropy_e  = FluidModel->GetEntropy();

      /* Loop over the faces that are treated simultaneously. */
      for(unsigned short l=0; l<nFaceSimul; ++l) {
        const unsigned short llNVar = l*nVar;

        /* Loop over the integration points of the face. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution
             for this integration point. */
          const su2double *UL = solIntL + i*NPad + llNVar;
          su2double       *UR = solIntR + i*NPad + llNVar;

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
      }

      break;
    }

    case STATIC_SUPERSONIC_INFLOW_PT: {

      /* Supersonic inlet with specified pressure, temperature and Mach
         number. Retrieve the non-dimensional static pressure and
         temperature as well as the three components of the Mach number. */
      su2double P_static = config->GetRiemann_Var1(Marker_Tag);
      su2double T_static = config->GetRiemann_Var2(Marker_Tag);
      auto Mach = config->GetRiemann_FlowDir(Marker_Tag);

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

      /* Loop over the number of faces that are treated simultaneously. */
      for(unsigned short l=0; l<nFaceSimul; ++l) {
        const unsigned short llNVar = l*nVar;

        /* Loop over the integration points and set the right state. */
        for(unsigned short i=0; i<nInt; ++i) {
          su2double *UR = solIntR + i*NPad + llNVar;
          for(unsigned short mm=0; mm<nVar; ++mm)
            UR[mm] = UCons[mm];
        }
      }

      break;
    }

    case STATIC_SUPERSONIC_INFLOW_PD: {

      /* Supersonic inlet with specified pressure, temperature and Mach
         number. Retrieve the non-dimensional static pressure and
         temperature as well as the three components of the Mach number. */
      su2double P_static   = config->GetRiemann_Var1(Marker_Tag);
      su2double Rho_static = config->GetRiemann_Var2(Marker_Tag);
      auto Mach            = config->GetRiemann_FlowDir(Marker_Tag);

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

      /* Loop over the number of faces that are treated simultaneously. */
      for(unsigned short l=0; l<nFaceSimul; ++l) {
        const unsigned short llNVar = l*nVar;

        /* Loop over the integration points and set the right state. */
        for(unsigned short i=0; i<nInt; ++i) {
          su2double *UR = solIntR + i*NPad + llNVar;
          for(unsigned short mm=0; mm<nVar; ++mm)
            UR[mm] = UCons[mm];
        }
      }

      break;
    }

    case DENSITY_VELOCITY: {

      /* Subsonic inlet with specified density, velocity magnitude and
         flow direction. Retrieve the non-dimensional data. */
      su2double Density_e = config->GetRiemann_Var1(Marker_Tag);
      su2double VelMag_e  = config->GetRiemann_Var2(Marker_Tag);
      auto Flow_Dir       = config->GetRiemann_FlowDir(Marker_Tag);

      Density_e /= config->GetDensity_Ref();
      VelMag_e  /= config->GetVelocity_Ref();

      /* Loop over the faces that are treated simultaneously. */
      for(unsigned short l=0; l<nFaceSimul; ++l) {
        const unsigned short llNVar = l*nVar;

        /* Loop over the integration points of the face. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution
             for this integration point. */
          const su2double *UL = solIntL + i*NPad + llNVar;
          su2double       *UR = solIntR + i*NPad + llNVar;

          /* Compute the total energy per unit mass, which is extrapolated. */
          const su2double Energy_e = UL[nDim+1]/UL[0];

          /* Set the conservative variables of the right state. */
          UR[0]      = Density_e;
          UR[nDim+1] = Density_e*Energy_e;

          const su2double momTot = Density_e*VelMag_e;
          for(unsigned short iDim=0; iDim<nDim; ++iDim)
            UR[iDim+1] = momTot*Flow_Dir[iDim];
        }
      }

      break;
    }

    case STATIC_PRESSURE: {

      /* Subsonic outlet with specified pressure.
         Retrieve the non-dimensional data. */
      const su2double Pressure_e = config->GetRiemann_Var1(Marker_Tag)
                                 / config->GetPressure_Ref();

      /* Loop over the faces that are treated simultaneously. */
      for(unsigned short l=0; l<nFaceSimul; ++l) {
        const unsigned short llNVar = l*nVar;

        /* Loop over the integration points of the face. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution
             for this integration point. */
          const su2double *UL = solIntL + i*NPad + llNVar;
          su2double       *UR = solIntR + i*NPad + llNVar;

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

      /* Two-dimensional simulation.
         Loop over the faces that are treated simultaneously. */
      for(unsigned short l=0; l<nFaceSimul; ++l) {
        const unsigned short llNVar = l*nVar;

        /* Loop over the integration points of the face. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution and the normals
             for this integration point. */
          const su2double *UL      = solIntL + i*NPad + llNVar;
                su2double *UR      = solIntR + i*NPad + llNVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

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
      }

      break;
    }

    case 3: {

      /* Three-dimensional simulation.
         Loop over the faces that are treated simultaneously. */
      for(unsigned short l=0; l<nFaceSimul; ++l) {
        const unsigned short llNVar = l*nVar;

        /* Loop over the integration points of the face. */
        for(unsigned short i=0; i<nInt; ++i) {

          /* Easier storage of the left and right solution and the normals
             for this integration point. */
          const su2double *UL      = solIntL + i*NPad + llNVar;
                su2double *UR      = solIntR + i*NPad + llNVar;
          const su2double *normals = surfElem[l].metricNormalsFace.data() + i*(nDim+1);

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
                                        CNumerics                *conv_numerics,
                                        su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array work is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Compute the right state by applying the inviscid wall BC's. */
    BoundaryStates_Euler_Wall(config, llEnd, NPad, &surfElem[l],
                              solIntL, solIntR);

    /* The remainder of the contribution of ths boundary faces to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, llEnd, NPad, conv_numerics, &surfElem[l],
                                 solIntL, solIntR, work, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::BC_Far_Field(CConfig                  *config,
                                       const unsigned long      surfElemBeg,
                                       const unsigned long      surfElemEnd,
                                       const CSurfaceElementFEM *surfElem,
                                       su2double                *resFaces,
                                       CNumerics                *conv_numerics,
                                       su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array work is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Set the right state in the integration points to the free stream value. */
    for(unsigned short ll=0; ll<llEnd; ++ll) {
      const unsigned short llNVar = ll*nVar;
      for(unsigned short i=0; i<nInt; ++i) {
        su2double *UR = solIntR + i*NPad + llNVar;
        for(unsigned short j=0; j<nVar; ++j)
          UR[j] = ConsVarFreeStream[j];
      }
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, llEnd, NPad, conv_numerics, &surfElem[l],
                                 solIntL, solIntR, work, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::BC_Sym_Plane(CConfig                  *config,
                                       const unsigned long      surfElemBeg,
                                       const unsigned long      surfElemEnd,
                                       const CSurfaceElementFEM *surfElem,
                                       su2double                *resFaces,
                                       CNumerics                *conv_numerics,
                                       su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array work is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Make a distinction between two and three space dimensions
       in order to have the most efficient code. */
    switch( nDim ) {

      case 2: {

        /* 2D simulation. Loop over the faces that are treated simultaneously
           and the integration points to apply the symmetry condition. */
        for(unsigned short ll=0; ll<llEnd; ++ll) {
          const unsigned short llNVar = ll*nVar;
          for(unsigned short i=0; i<nInt; ++i) {

            /* Easier storage of the left and right solution and the normals
               for this integration point. */
            const su2double *UL      = solIntL + NPad*i + llNVar;
                  su2double *UR      = solIntR + NPad*i + llNVar;
            const su2double *normals = surfElem[ll+l].metricNormalsFace.data() + i*(nDim+1);

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
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /* 3D simulation. Loop over the faces that are treated simultaneously
           and the integration points to apply the symmetry condition. */
        for(unsigned short ll=0; ll<llEnd; ++ll) {
          const unsigned short llNVar = ll*nVar;
          for(unsigned short i=0; i<nInt; ++i) {

            /* Easier storage of the left and right solution and the normals
               for this integration point. */
            const su2double *UL      = solIntL + NPad*i + llNVar;
                  su2double *UR      = solIntR + NPad*i + llNVar;
            const su2double *normals = surfElem[ll+l].metricNormalsFace.data() + i*(nDim+1);

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
        }

        break;
      }
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, llEnd, NPad, conv_numerics, &surfElem[l],
                                 solIntL, solIntR, work, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::BC_Supersonic_Outlet(CConfig                  *config,
                                               const unsigned long      surfElemBeg,
                                               const unsigned long      surfElemEnd,
                                               const CSurfaceElementFEM *surfElem,
                                               su2double                *resFaces,
                                               CNumerics                *conv_numerics,
                                               su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array work is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Set the right state in the integration points to the left state, i.e.
       no boundary condition is applied for a supersonic outlet. */
    for(unsigned short mm=0; mm<NPad*nInt; ++mm)
       solIntR[mm] = solIntL[mm];

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, llEnd, NPad, conv_numerics, &surfElem[l],
                                 solIntL, solIntR, work, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::BC_Inlet(CConfig                  *config,
                                   const unsigned long      surfElemBeg,
                                   const unsigned long      surfElemEnd,
                                   const CSurfaceElementFEM *surfElem,
                                   su2double                *resFaces,
                                   CNumerics                *conv_numerics,
                                   unsigned short           val_marker,
                                   su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array work is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Compute the right state by applying the subsonic inlet BC's. */
    BoundaryStates_Inlet(config, llEnd, NPad, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, llEnd, NPad, conv_numerics, &surfElem[l],
                                 solIntL, solIntR, work, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::BC_Outlet(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    unsigned short           val_marker,
                                    su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array work is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Compute the right state by applying the subsonic inlet BC's. */
    BoundaryStates_Outlet(config, llEnd, NPad, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, llEnd, NPad, conv_numerics, &surfElem[l],
                                 solIntL, solIntR, work, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::BC_Riemann(CConfig                  *config,
                                     const unsigned long      surfElemBeg,
                                     const unsigned long      surfElemEnd,
                                     const CSurfaceElementFEM *surfElem,
                                     su2double                *resFaces,
                                     CNumerics                *conv_numerics,
                                     unsigned short           val_marker,
                                     su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array work is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Compute the right state by applying the subsonic inlet BC's. */
    BoundaryStates_Riemann(config, llEnd, NPad, &surfElem[l], val_marker, solIntL, solIntR);

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, llEnd, NPad, conv_numerics, &surfElem[l],
                                 solIntL, solIntR, work, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::BC_Custom(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray) {

  /* Initialization of the counter in resFaces. */
  unsigned long indResFaces = 0;

  /* Determine the number of faces that are treated simultaneously
     in the matrix products to obtain good gemm performance. */
  const unsigned short nPadInput  = config->GetSizeMatMulPadding();
  const unsigned short nFaceSimul = nPadInput/nVar;

  /* Determine the minimum padded size in the matrix multiplications, which
     corresponds to 64 byte alignment. */
  const unsigned short nPadMin = 64/sizeof(passivedouble);

  /*--- Get the physical time if necessary. ---*/
  su2double time = 0.0;
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

  /*--- Loop over the requested range of surface faces. Multiple faces
        are treated simultaneously to improve the performance of the matrix
        multiplications. As a consequence, the update of the counter l
        happens at the end of this loop section. ---*/
  for(unsigned long l=surfElemBeg; l<surfElemEnd;) {

    /* Determine the end index for this chunk of faces and the padded
       N value in the gemm computations. */
    unsigned long lEnd;
    unsigned short ind, llEnd, NPad;

    MetaDataChunkOfElem(surfElem, l, surfElemEnd, nFaceSimul,
                        nPadMin, lEnd, ind, llEnd, NPad);

    /* Get the necessary data from the standard face. */
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Set the pointers for the local arrays. ---*/
    su2double *solIntL = workArray;
    su2double *solIntR = solIntL + NPad*nInt;
    su2double *work    = solIntR + NPad*nInt;

    /* Compute the left states in the integration points of the chunk of
       faces. The array work is used as temporary storage inside the
       function LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(config, llEnd, NPad, &surfElem[l],
                                            work, solIntL);

    /* Check for a verification solution. */
    if( VerificationSolution ) {

      /*--- Loop over the number of simultaneously treated faces and integration points
            to compute the right state for the boundary conditions. ---*/
      for(unsigned short ll=0; ll<llEnd; ++ll) {
        for(unsigned short i=0; i<nInt; ++i) {

          /* Determine the pointer to the coordinates of this integration
             point and the pointer to the solution and call the function
             GetBCState to determine the actual boundary state. */
          const su2double *coor = surfElem[ll+l].coorIntegrationPoints.data() + i*nDim;
          su2double *UR   = solIntR + NPad*i + ll*nVar;

          VerificationSolution->GetBCState(coor, time, UR);
        }
      }
    }
    else {

      /* The user must specify the custom BC's here. */
      SU2_MPI::Error("Implement customized boundary conditions here.", CURRENT_FUNCTION);
    }

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualInviscidBoundaryFace(config, llEnd, NPad, conv_numerics, &surfElem[l],
                                 solIntL, solIntR, work, resFaces, indResFaces);

    /* Update the value of the counter l to the end index of the
       current chunk. */
    l = lEnd;
  }
}

void CFEM_DG_EulerSolver::ResidualInviscidBoundaryFace(
                                      CConfig                  *config,
                                      const unsigned short     nFaceSimul,
                                      const unsigned short     NPad,
                                      CNumerics                *conv_numerics,
                                      const CSurfaceElementFEM *surfElem,
                                      su2double                *solInt0,
                                      su2double                *solInt1,
                                      su2double                *fluxes,
                                      su2double                *resFaces,
                                      unsigned long            &indResFaces) {

  /*--- Get the required information from the standard face, which is the
        same for all faces considered. ---*/
  const unsigned short ind     = surfElem[0].indStandardElement;
  const unsigned short nInt    = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFs   = standardBoundaryFacesSol[ind].GetNDOFsFace();
  const su2double     *weights = standardBoundaryFacesSol[ind].GetWeightsIntegration();

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Compute the fluxes in the integration points using the   ---*/
  /*---         approximate Riemann solver.                              ---*/
  /*------------------------------------------------------------------------*/

  /* Store the normals and grid velocities of the faces in an array of
     pointers to be used in the call to ComputeInviscidFluxesFace. */
  const su2double* arrNorm[SIZE_ARR_NORM];
  const su2double* arrGridVel[SIZE_ARR_NORM];

  if(nFaceSimul > SIZE_ARR_NORM)
    SU2_MPI::Error("SIZE_ARR_NORM is too small. Increase it or decrease ALIGNED_BYTES_MATMUL",
                   CURRENT_FUNCTION);

  for(unsigned short l=0; l<nFaceSimul; ++l) {
    arrNorm[l]    = surfElem[l].metricNormalsFace.data();
    arrGridVel[l] = surfElem[l].gridVelocities.data();
  }

  /* General function to compute the inviscid fluxes, using an approximate
     Riemann solver in the integration points. */
  ComputeInviscidFluxesFace(config, nFaceSimul, NPad, nInt, arrNorm, arrGridVel,
                            solInt0, solInt1, fluxes, conv_numerics);

  /* Multiply the fluxes with the integration weight of the corresponding
     integration point. */
  for(unsigned short i=0; i<nInt; ++i) {
    su2double *flux = fluxes + i*NPad;

    for(unsigned short j=0; j<NPad; ++j)
      flux[j] *= weights[i];
  }

  /*------------------------------------------------------------------------*/
  /*--- Step 2: Compute the contribution to the residuals from the       ---*/
  /*---         integration over the surface element.                    ---*/
  /*------------------------------------------------------------------------*/

  /* Get the correct form of the basis functions needed for the matrix
     multiplication to compute the residual. */
  const su2double *basisFaceTrans = standardBoundaryFacesSol[ind].GetBasisFaceIntegrationTranspose();

  /* Call the general function to carry out the matrix product. Use solint0
     as temporary storage for the result. */
  blasFunctions->gemm(nDOFs, NPad, nInt, basisFaceTrans, fluxes, solInt0, config);

  /* Loop over the number of simultaneously treated faces to store the
     residual in the correct locations in resFaces. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {
    const unsigned short llNVar = l*nVar;

    /* Easier storage of the position in the residual array for this face
       and update the corresponding counter. */
    su2double *resFace = resFaces + indResFaces*nVar;
    indResFaces       += nDOFs;

    /* Loop over the DOFs and copy the data. */
    for(unsigned short i=0; i<nDOFs; ++i)
      for(unsigned short mm=0; mm<nVar; ++mm)
        resFace[nVar*i+mm] = solInt0[NPad*i+llNVar+mm];
  }
}

void CFEM_DG_EulerSolver::LeftStatesIntegrationPointsBoundaryFace(
                                             CConfig                  *config,
                                             const unsigned short     nFaceSimul,
                                             const unsigned short     NPad,
                                             const CSurfaceElementFEM *surfElem,
                                             su2double                *solFace,
                                             su2double                *solIntL) {

  /* Get the required information from the corresponding standard face, which is the
     same for all simultaneously treated faces. */
  const unsigned short ind   = surfElem[0].indStandardElement;
  const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();
  const su2double *basisFace = standardBoundaryFacesSol[ind].GetBasisFaceIntegration();

  /* Loop over the faces that are treated simultaneously. */
  for(unsigned short l=0; l<nFaceSimul; ++l) {
    const unsigned short llNVar = l*nVar;

    /* Easier storage of the DOFs of the face. */
    const unsigned long *DOFs = surfElem[l].DOFsSolFace.data();

    /* The numbering of the DOFs is given for the entire grid. However, in the
       working vectors for the solution a numbering per time level is used.
       Therefore an offset must be applied, which can be obtained from the
       corresponding volume element. Note that this offset is non-negative. */
    const unsigned long  volID     = surfElem[l].volElemID;
    const unsigned short timeLevel = volElem[volID].timeLevel;
    const unsigned long  offset    = volElem[volID].offsetDOFsSolLocal
                                   - volElem[volID].offsetDOFsSolThisTimeLevel;

    /* Copy the solution of the DOFs of the face such that it is contiguous
       in memory. */
    for(unsigned short i=0; i<nDOFs; ++i) {

      const su2double *solDOF = VecWorkSolDOFs[timeLevel].data() + nVar*(DOFs[i]-offset);
      su2double       *sol    = solFace + NPad*i + llNVar;
      for(unsigned short mm=0; mm<nVar; ++mm)
        sol[mm] = solDOF[mm];
    }
  }

  /* Call the general function to carry out the matrix product. */
  blasFunctions->gemm(nInt, NPad, nDOFs, basisFace, solFace, solIntL, config);
}

void CFEM_DG_EulerSolver::ComputeInviscidFluxesFace(CConfig              *config,
                                                    const unsigned short nFaceSimul,
                                                    const unsigned short NPad,
                                                    const unsigned long  nPoints,
                                                    const su2double      *normalsFace[],
                                                    const su2double      *gridVelsFace[],
                                                    const su2double      *solL,
                                                    const su2double      *solR,
                                                    su2double            *fluxes,
                                                    CNumerics            *numerics) {

  /* Easier storage of the specific heat ratio. */
  const su2double gm1 = Gamma_Minus_One;

  /* Make a distinction between the several Riemann solvers. */
  switch( config->GetRiemann_Solver_FEM() ) {

    case UPWIND::ROE: {

      /* Roe's approximate Riemann solver. Easier storage of the cut off
         value for the entropy correction. */
      const su2double Delta = config->GetEntropyFix_Coeff();

      /* Make a distinction between two and three space dimensions
         in order to have the most efficient code. */
      switch( nDim ) {

        case 2: {

          /* Two dimensional simulation. Loop over the number of faces treated
             simultaneously. */
          for(unsigned short l=0; l<nFaceSimul; ++l) {

            /* Easier storage for some variables of this face. */
            const su2double *normals   = normalsFace[l];
            const su2double *gridVels  = gridVelsFace[l];
            const unsigned short lNVar = l*nVar;

            /* Loop over the number of points for this face. */
            for(unsigned long i=0; i<nPoints; ++i) {

              /* Easier storage of the left and right solution, the face normals,
                 the grid velocities and the flux vector for this point. */
              const unsigned long offPointer = i*NPad + lNVar;

              const su2double *UL      = solL + offPointer;
              const su2double *UR      = solR + offPointer;
              const su2double *norm    = normals + i*(nDim+1);
              const su2double *gridVel = gridVels + i*nDim;
                    su2double *flux    = fluxes + offPointer;

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
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /* Three dimensional simulation. Loop over the number of faces treated
             simultaneously. */
          for(unsigned short l=0; l<nFaceSimul; ++l) {

            /* Easier storage for some variables of this face. */
            const su2double *normals   = normalsFace[l];
            const su2double *gridVels  = gridVelsFace[l];
            const unsigned short lNVar = l*nVar;

            /* Loop over the number of points for this face. */
            for(unsigned long i=0; i<nPoints; ++i) {

              /* Easier storage of the left and right solution, the face normals,
                 the grid velocities and the flux vector for this point. */
              const unsigned long offPointer = i*NPad + lNVar;

              const su2double *UL      = solL + offPointer;
              const su2double *UR      = solR + offPointer;
              const su2double *norm    = normals + i*(nDim+1);
              const su2double *gridVel = gridVels + i*nDim;
                    su2double *flux    = fluxes + offPointer;

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
          }

          break;
        }
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case UPWIND::LAX_FRIEDRICH: {

      /* Local Lax-Friedrich (Rusanov) flux Make a distinction between two and
         three space dimensions in order to have the most efficient code. */
      switch( nDim ) {

        case 2: {

          /* Two dimensional simulation. Loop over the number of faces treated
             simultaneously. */
          for(unsigned short l=0; l<nFaceSimul; ++l) {

            /* Easier storage for some variables of this face. */
            const su2double *normals   = normalsFace[l];
            const su2double *gridVels  = gridVelsFace[l];
            const unsigned short lNVar = l*nVar;

            /* Loop over the number of points for this face. */
            for(unsigned long i=0; i<nPoints; ++i) {

              /* Easier storage of the left and right solution, the face normals,
                 the grid velocities and the flux vector for this point. */
              const unsigned long offPointer = i*NPad + lNVar;

              const su2double *UL      = solL + offPointer;
              const su2double *UR      = solR + offPointer;
              const su2double *norm    = normals + i*(nDim+1);
              const su2double *gridVel = gridVels + i*nDim;
                    su2double *flux    = fluxes + offPointer;

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
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /* Three dimensional simulation. Loop over the number of faces treated
             simultaneously. */
          for(unsigned short l=0; l<nFaceSimul; ++l) {

            /* Easier storage for some variables of this face. */
            const su2double *normals   = normalsFace[l];
            const su2double *gridVels  = gridVelsFace[l];
            const unsigned short lNVar = l*nVar;

            /* Loop over the number of points for this face. */
            for(unsigned long i=0; i<nPoints; ++i) {

              /* Easier storage of the left and right solution, the face normals,
                 the grid velocities and the flux vector for this point. */
              const unsigned long offPointer = i*NPad + lNVar;

              const su2double *UL      = solL + offPointer;
              const su2double *UR      = solR + offPointer;
              const su2double *norm    = normals + i*(nDim+1);
              const su2double *gridVel = gridVels + i*nDim;
                    su2double *flux    = fluxes + offPointer;

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

      /* Loop over the number of faces treated simultaneously. */
      for(unsigned short l=0; l<nFaceSimul; ++l) {

        /* Easier storage for some variables of this face. */
        const su2double *normals   = normalsFace[l];
        const su2double *gridVels  = gridVelsFace[l];
        const unsigned short lNVar = l*nVar;

        /* Loop over the number of points for this face. */
        for(unsigned long i=0; i<nPoints; ++i) {

          /* Easier storage of the left and right solution, the face normals,
             the grid velocities and the flux vector for this point. */
          const unsigned long offPointer = i*NPad + lNVar;

          const su2double *UL      = solL + offPointer;
          const su2double *UR      = solR + offPointer;
          const su2double *norm    = normals + i*(nDim+1);
          const su2double *gridVel = gridVels + i*nDim;
                su2double *flux    = fluxes + offPointer;

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
      }

      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        delete [] Jacobian_i[iVar];
        delete [] Jacobian_j[iVar];
      }
      delete [] Jacobian_i;
      delete [] Jacobian_j;

      Jacobian_i = nullptr;
      Jacobian_j = nullptr;
    }
  }
}

void CFEM_DG_EulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iVar;
  unsigned long index;

  const bool compressible = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  string restart_filename = config->GetSolution_FileName();

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned short rbuf_NotMatching = 0;
  unsigned long nDOF_Read = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  restart_filename = config->GetFilename(restart_filename, "", val_iter);

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
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

  if (rbuf_NotMatching != 0)
    SU2_MPI::Error(string("The solution file ") + restart_filename +
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
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    unsigned long nBadDOFsLoc = nBadDOFs;
    SU2_MPI::Reduce(&nBadDOFsLoc, &nBadDOFs, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
#endif

    if((rank == MASTER_NODE) && (nBadDOFs != 0))
      cout << "Warning. The initial solution contains "<< nBadDOFs << " DOFs that are not physical." << endl;
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars;
  delete [] Restart_Data;
  Restart_Vars = nullptr; Restart_Data = nullptr;

}
