/*!
 * \file solution_direct_mean_fem.cpp
 * \brief Main subroutines for solving finite element flow problems (Euler, Navier-Stokes, etc.).
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 4.0.2 "Cardinal"
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

/* LIBXSMM include files, if supported. */
#ifdef HAVE_LIBXSMM
#include "libxsmm.h"
#endif

/* MKL or BLAS include files, if supported. */
#ifdef HAVE_MKL
#include "mkl.h"
#elif HAVE_CBLAS
#include "cblas.h"
#endif

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(void) : CSolver() {
  
  /*--- Basic array initialization ---*/
  
  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  
  /*--- Surface-based array initialization ---*/
  
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL; Surface_CSideForce_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;
  
  Surface_CLift = NULL; Surface_CDrag = NULL; Surface_CSideForce = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Initialization of the boolean symmetrizingTermsPresent. ---*/
  symmetrizingTermsPresent = true;
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
  
  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL; CEff_Inv = NULL;
  CMx_Inv = NULL;   CMy_Inv = NULL;   CMz_Inv = NULL;
  CFx_Inv = NULL;   CFy_Inv = NULL;   CFz_Inv = NULL;
  
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL; Surface_CSideForce_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL;   Surface_CFy_Inv = NULL;   Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL;   Surface_CMy_Inv = NULL;   Surface_CMz_Inv = NULL;
  
  Surface_CLift = NULL; Surface_CDrag = NULL; Surface_CSideForce = NULL; Surface_CEff = NULL;
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

  nMeshPoints = DGGeometry->GetNMeshPoints();
  meshPoints  = DGGeometry->GetMeshPoints();

  nMatchingInternalFaces = DGGeometry->GetNMatchingFaces();
  matchingInternalFaces  = DGGeometry->GetMatchingFaces();

  boundaries = DGGeometry->GetBoundaries();

  nStandardBoundaryFacesSol = DGGeometry->GetNStandardBoundaryFacesSol();
  nStandardElementsSol      = DGGeometry->GetNStandardElementsSol();
  nStandardMatchingFacesSol = DGGeometry->GetNStandardMatchingFacesSol();

  standardBoundaryFacesSol = DGGeometry->GetStandardBoundaryFacesSol();
  standardElementsSol      = DGGeometry->GetStandardElementsSol();
  standardMatchingFacesSol = DGGeometry->GetStandardMatchingFacesSol();

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
  
  /*--- Perform the non-dimensionalization for the flow equations using the
        specified reference values. ---*/
  
  SetNondimensionalization(geometry, config, iMesh);

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
  
  CDrag_Inv         = new su2double[nMarker];
  CLift_Inv         = new su2double[nMarker];
  CSideForce_Inv    = new su2double[nMarker];
  CMx_Inv           = new su2double[nMarker];
  CMy_Inv           = new su2double[nMarker];
  CMz_Inv           = new su2double[nMarker];
  CEff_Inv          = new su2double[nMarker];
  CFx_Inv           = new su2double[nMarker];
  CFy_Inv           = new su2double[nMarker];
  CFz_Inv           = new su2double[nMarker];
  
  Surface_CLift_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce_Inv = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CLift          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Init total coefficients ---*/
  
  Total_CDrag   = 0.0;	Total_CLift        = 0.0;  Total_CSideForce   = 0.0;
  Total_CMx     = 0.0;	Total_CMy          = 0.0;  Total_CMz          = 0.0;
  Total_CFx     = 0.0;	Total_CFy          = 0.0;  Total_CFz          = 0.0;
  Total_CEff    = 0.0;
  
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

  /*--- Determine the global number of DOFs. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nDOFsLocOwned, &nDOFsGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nDOFsGlobal = nDOFsLocOwned;
#endif

  /*--- Allocate the memory to store the time steps, residuals, etc. ---*/

  VecDeltaTime.resize(nVolElemOwned);
  VecResDOFs.resize(nVar*nDOFsLocOwned);
  nEntriesResFaces.assign(nDOFsLocOwned+1, 0);
  startLocResFacesMarkers.resize(nMarker);

  /*--- Determine the size of the vector to store residuals that come from the
        integral over the faces and determine the number of entries in this
        vector for the owned DOFs. ---*/
  symmetrizingTermsPresent = false;
  if(config->GetViscous() && (fabs(config->GetTheta_Interior_Penalty_DGFEM()) > 1.e-8))
    symmetrizingTermsPresent = true;

  /*--- First the internal matching faces. ---*/
  unsigned long sizeVecResFaces = 0;
  for(unsigned long i=0; i<nMatchingInternalFaces; ++i) {

    /* The terms that only contribute to the DOFs located on the face. */
    const unsigned short ind = matchingInternalFaces[i].indStandardElement;
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    sizeVecResFaces += nDOFsFace0;
    for(unsigned short j=0; j<nDOFsFace0; ++j)
      ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolFaceSide0[j]+1];

    if(matchingInternalFaces[i].elemID1 < nVolElemOwned) {
      sizeVecResFaces += nDOFsFace1;
      for(unsigned short j=0; j<nDOFsFace1; ++j)
        ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolFaceSide1[j]+1];
    }

    /* The symmetrizing terms, if present, contribute to all
       the DOFs of the adjacent elements. */
    if( symmetrizingTermsPresent ) {
      const unsigned short nDOFsElem0 = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
      const unsigned short nDOFsElem1 = standardMatchingFacesSol[ind].GetNDOFsElemSide1();

      sizeVecResFaces += nDOFsElem0;
      for(unsigned short j=0; j<nDOFsElem0; ++j)
        ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolElementSide0[j]+1];

      if(matchingInternalFaces[i].elemID1 < nVolElemOwned) {
        sizeVecResFaces += nDOFsElem1;
        for(unsigned short j=0; j<nDOFsElem1; ++j)
          ++nEntriesResFaces[matchingInternalFaces[i].DOFsSolElementSide1[j]+1];
      }
    }
  }

  /* And the physical boundary faces. Exclude the periodic boundaries,
     because these are not physical boundaries. */
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    startLocResFacesMarkers[iMarker] = sizeVecResFaces;

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
      }
    }
  }

  /*--- Put nEntriesResFaces in cumulative storage format and allocate the
        memory for entriesResFaces and VecResFaces. ---*/
  for(unsigned long i=0; i<nDOFsLocOwned; ++i)
    nEntriesResFaces[i+1] += nEntriesResFaces[i];

  entriesResFaces.resize(nEntriesResFaces[nDOFsLocOwned]);
  VecResFaces.resize(nVar*sizeVecResFaces);

  /*--- Repeat the loops over the internal and boundary faces, but now store
        the enties in entriesResFaces. A counter variable is needed to keep
        track of the appropriate location in entriesResFaces. ---*/
  vector<unsigned long> counterEntries = nEntriesResFaces;

  /* First the loop over the internal matching faces. */
  sizeVecResFaces = 0;
  for(unsigned long i=0; i<nMatchingInternalFaces; ++i) {

    /* The terms that only contribute to the DOFs located on the face. */
    const unsigned short ind = matchingInternalFaces[i].indStandardElement;
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    for(unsigned short j=0; j<nDOFsFace0; ++j) {
      unsigned long jj    = counterEntries[matchingInternalFaces[i].DOFsSolFaceSide0[j]]++;
      entriesResFaces[jj] = sizeVecResFaces++;
    }

    if(matchingInternalFaces[i].elemID1 < nVolElemOwned) {
      for(unsigned short j=0; j<nDOFsFace1; ++j) {
        unsigned long jj    = counterEntries[matchingInternalFaces[i].DOFsSolFaceSide1[j]]++;
        entriesResFaces[jj] = sizeVecResFaces++;
      }
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

      if(matchingInternalFaces[i].elemID1 < nVolElemOwned) {
        for(unsigned short j=0; j<nDOFsElem1; ++j) {
          unsigned long jj    = counterEntries[matchingInternalFaces[i].DOFsSolElementSide1[j]]++;
          entriesResFaces[jj] = sizeVecResFaces++;
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

  /*--- Set up the persistent communication for the conservative variables. ---*/

  Prepare_MPI_Communication(DGGeometry, config);
  
  /*--- Perform the MPI communication of the solution including the possible
        self communication for the periodic data. Correct for rotational
        periodicity afterwards, if needed. ---*/
#ifdef HAVE_MPI
  if( nCommRequests ) {
    MPI_Startall(nCommRequests, commRequests.data());
    SU2_MPI::Waitall(nCommRequests, commRequests.data(), MPI_STATUSES_IGNORE);
  }
#endif

  SelfCommunication();
  CorrectForRotationalPeriodicity();
}

CFEM_DG_EulerSolver::~CFEM_DG_EulerSolver(void) {

  if( FluidModel) delete FluidModel;

  /*--- Array deallocation ---*/
  if (CDrag_Inv != NULL)         delete [] CDrag_Inv;
  if (CLift_Inv != NULL)         delete [] CLift_Inv;
  if (CSideForce_Inv != NULL)    delete [] CSideForce_Inv;
  if (CMx_Inv != NULL)           delete [] CMx_Inv;
  if (CMy_Inv != NULL)           delete [] CMy_Inv;
  if (CMz_Inv != NULL)           delete [] CMz_Inv;
  if (CFx_Inv != NULL)           delete [] CFx_Inv;
  if (CFy_Inv != NULL)           delete [] CFy_Inv;
  if (CFz_Inv != NULL)           delete [] CFz_Inv;
  if (Surface_CLift_Inv != NULL) delete[] Surface_CLift_Inv;
  if (Surface_CDrag_Inv != NULL) delete[] Surface_CDrag_Inv;
  if (Surface_CSideForce_Inv != NULL) delete[] Surface_CSideForce_Inv;
  if (Surface_CEff_Inv != NULL) delete[] Surface_CEff_Inv;
  if (Surface_CFx_Inv != NULL)  delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv != NULL)  delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv != NULL)  delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;
  if (Surface_CLift != NULL)    delete [] Surface_CLift;
  if (Surface_CDrag != NULL)    delete [] Surface_CDrag;
  if (Surface_CSideForce != NULL) delete [] Surface_CSideForce;
  if (Surface_CEff != NULL) delete [] Surface_CEff;
  if (Surface_CFx != NULL)      delete [] Surface_CFx;
  if (Surface_CFy != NULL)      delete [] Surface_CFy;
  if (Surface_CFz != NULL)      delete [] Surface_CFz;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)          delete [] CEff_Inv;
  
  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;

#ifdef HAVE_MPI

  /*--- Release the memory of the persistent communication and the derived
        data types. ---*/
  for(int i=0; i<nCommRequests; ++i) MPI_Request_free(&commRequests[i]);
  for(int i=0; i<nCommRequests; ++i) MPI_Type_free(&commTypes[i]);

#endif
}

void CFEM_DG_EulerSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) {
  
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
  
  /*--- Write output to the console if this is the master node and first domain ---*/
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {
    
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
  const vector<int>                    &ranksComm       = FEMGeometry->GetRanksComm();
  const vector<vector<unsigned long> > &elementsSend    = FEMGeometry->GetEntitiesSend();
  const vector<vector<unsigned long> > &elementsReceive = FEMGeometry->GetEntitiesReceive();

  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Find out whether or not self communicationn is present.    ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the rank inside MPI_COMM_WORLD. */
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* Loop over the entries of ranksComm and check for self communication. */
  for(unsigned long i=0; i<ranksComm.size(); ++i) {
    if(ranksComm[i] == rank) {

      /* The send and receive elements for this entry in elementsSend and
         elementsReceive correspond to self communication. Copy the data. */
      elementsReceiveSelfComm = elementsReceive[i];
      elementsSendSelfComm    = elementsSend[i];
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Set up the persistent MPI communication.                   ---*/
  /*--------------------------------------------------------------------------*/

#ifdef HAVE_MPI

  /*--- Determine the number of ranks with which communication takes place,
        multiply this number by two to obtain the number of communication
        requests (send and receive) and allocate the necessary memory. ---*/
  nCommRequests = ranksComm.size();
  if( elementsReceiveSelfComm.size() ) --nCommRequests;
  nCommRequests *= 2;

  commRequests.resize(nCommRequests);
  commTypes.resize(nCommRequests);

  /*--- Loop over the ranks for which communication takes place and exclude
        the self communication. ---*/
  unsigned int nn = 0;
  for(unsigned long i=0; i<ranksComm.size(); ++i) {
    if(ranksComm[i] != rank) {

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

      /* Create the communication request for this send operation. Afterwards
         update the counter nn for the next request. */
      MPI_Send_init(VecSolDOFs.data(), 1, commTypes[nn], ranksComm[i],
                    ranksComm[i], MPI_COMM_WORLD, &commRequests[nn]);
      ++nn;

      /*--- Determine the derived data type for receiving the data. ---*/
      const unsigned int nElemRecv = elementsReceive[i].size();
      blockLen.resize(nElemRecv);
      displ.resize(nElemRecv);

      for(unsigned int j=0; j<nElemRecv; ++j) {
        const unsigned long jj = elementsReceive[i][j];
        blockLen[j] = nVar*volElem[jj].nDOFsSol;
        displ[j]    = nVar*volElem[jj].offsetDOFsSolLocal;
      }

      MPI_Type_indexed(nElemRecv, blockLen.data(), displ.data(),
                       MPI_DOUBLE, &commTypes[nn]);
      MPI_Type_commit(&commTypes[nn]);

      /* Create the communication request for this receive operation. Afterwards
         update the counter nn for the next request. */
      MPI_Recv_init(VecSolDOFs.data(), 1, commTypes[nn], ranksComm[i],
                    rank, MPI_COMM_WORLD, &commRequests[nn]);
      ++nn;
    }
  }

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 3. Store the information for the rotational periodic          ---*/
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

void CFEM_DG_EulerSolver::CorrectForRotationalPeriodicity(void) {

  /*--- Loop over the markers for which a rotational periodic
        correction must be applied to the momentum variables. ---*/
  unsigned int ii;
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

void CFEM_DG_EulerSolver::SelfCommunication(void) {

  /* Easier storage of the number of variables times the size of su2double. */
  const unsigned long sizeEntity = nVar*sizeof(su2double);

  /* Loop over the number of elements involved and copy the data of the DOFs. */
  const unsigned long nSelfElements = elementsSendSelfComm.size();
  for(unsigned long i=0; i<nSelfElements; ++i) {
    const unsigned long indS   = nVar*volElem[elementsSendSelfComm[i]].offsetDOFsSolLocal;
    const unsigned long indR   = nVar*volElem[elementsReceiveSelfComm[i]].offsetDOFsSolLocal;
    const unsigned long nBytes = volElem[elementsSendSelfComm[i]].nDOFsSol * sizeEntity;

    memcpy(VecSolDOFs.data()+indR, VecSolDOFs.data()+indS, nBytes);
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
  /* Carry out the self communication. */
  SelfCommunication();

  /*--- Complete the MPI communication, if needed. ---*/

#ifdef HAVE_MPI
  if( nCommRequests )
    SU2_MPI::Waitall(nCommRequests, commRequests.data(), MPI_STATUSES_IGNORE);
#endif

  /* Correct the vector quantities in the rotational periodic halo's. */
  CorrectForRotationalPeriodicity();
}


void CFEM_DG_EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
}

void CFEM_DG_EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

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

void CFEM_DG_EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /* Easier storage whether or not a viscous computation is carried out. */
  const bool viscous = config->GetViscous();

  /* Initialize the minimum and maximum time step. */
  Min_Delta_Time = 1.e25; Max_Delta_Time = 0.0;

  /* Easier storage of the CFL number. */
  const su2double CFL = config->GetCFL(iMesh);

  /*--- Check for a compressible solver. ---*/
  if(config->GetKind_Regime() == COMPRESSIBLE) {

    /*--- Loop over the owned volume elements. ---*/
    for(unsigned long i=0; i<nVolElemOwned; ++i) {

      /*--- Loop over the DOFs of this element and determine the maximum wave speed
            and the maximum value of the kinematic viscosity (if needed). ---*/
      su2double charVel2Max = 0.0, nuMax = 0.0;
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

        /* Update the kinematic viscosity, if a viscous computation is carried out. */
        if( viscous ) {
          const su2double nu = DensityInv*FluidModel->GetLaminarViscosity();
          nuMax = max(nuMax, nu);
        }
      }

      /*--- Compute the time step for the element and update the minimum and
            maximum value. Note that in the length scale the polynomial degree
            must be taken into account for the high order element. ---*/
      const unsigned short ind = volElem[i].indStandardElement;
      unsigned short nPoly = standardElementsSol[ind].GetNPoly();
      if(nPoly == 0) nPoly = 1;

      const su2double lenScaleInv = nPoly/volElem[i].lenScale;
      const su2double dtInv       = lenScaleInv*(charVel2Max + nuMax*lenScaleInv);

      VecDeltaTime[i] = CFL/dtInv;

      Min_Delta_Time = min(Min_Delta_Time, VecDeltaTime[i]);
      Max_Delta_Time = max(Max_Delta_Time, VecDeltaTime[i]);
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

  /*--- Compute the max and the min dt (in parallel) ---*/
#ifdef HAVE_MPI
  su2double rbuf_time = Min_Delta_Time;
  SU2_MPI::Allreduce(&rbuf_time, &Min_Delta_Time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  rbuf_time = Max_Delta_Time;
  SU2_MPI::Allreduce(&rbuf_time, &Max_Delta_Time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

void CFEM_DG_EulerSolver::Internal_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                            CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  /*--- Allocate the memory for some temporary storage needed to compute the
        internal residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. ---*/
  su2double *solInt, *fluxes;

#ifdef HAVE_MKL
  solInt = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  fluxes = (su2double *) mkl_malloc(nIntegrationMax*nVar*nDim*sizeof(su2double), 64);
#else
  vector<su2double> helpSolInt(nIntegrationMax*nVar);
  vector<su2double> helpFluxes(nIntegrationMax*nVar*nDim);
  solInt = helpSolInt.data();
  fluxes = helpFluxes.data();
#endif

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

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
    su2double *solDOFs = VecSolDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Call the general function to carry out the matrix product. */
    MatrixProduct(nInt, nVar, nDOFs, matBasisInt, solDOFs, solInt);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the inviscid fluxes, multiplied by minus the     ---*/
    /*---         integration weight, in the integration points.           ---*/
    /*------------------------------------------------------------------------*/ 

    /* Loop over the integration points. */
    unsigned short ll = 0;
    for(unsigned short i=0; i<nInt; ++i) {

      /*--- Compute the velocities and pressure in this integration point. ---*/
      const su2double *sol = solInt + nVar*i;

      const su2double DensityInv = 1.0/sol[0];
      su2double vel[3], Velocity2 = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {
        vel[iDim] = sol[iDim+1]*DensityInv;
        Velocity2 += vel[iDim]*vel[iDim];
      }

      su2double StaticEnergy = sol[nDim+1]*DensityInv - 0.5*Velocity2;

      FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
      const su2double Pressure = FluidModel->GetPressure();

      /* Easier storage of the metric terms in this integration point. The +1
         is present, because the first element of the metric terms is the
         Jacobian in the integration point. */
      const su2double *metricTerms = volElem[l].metricTerms + i*nMetricPerPoint + 1;

      /*--- Loop over the number of dimensions to compute the fluxes in the
            direction of the parametric coordinates. ---*/
      for(unsigned short iDim=0; iDim<nDim; ++iDim) {

        /* Pointer to the metric terms for this direction. */
        const su2double *metric = metricTerms + iDim*nDim;

        /* Compute the velocity in the direction of the current parametric coordinate. */
        su2double vPar = 0.0;
        for(unsigned short jDim=0; jDim<nDim; ++jDim)
          vPar += vel[jDim]*metric[jDim];

        /* Compute the flux, multiplied by minus the integration weight. */
        fluxes[ll++] = -weights[i]*sol[0]*vPar;

        for(unsigned short jDim=0; jDim<nDim; ++jDim)
          fluxes[ll++] = -weights[i]*(sol[jDim+1]*vPar + Pressure*metric[jDim]);

        fluxes[ll++] = -weights[i]*(sol[nDim+1] + Pressure)*vPar;
      }
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
    /*---         integration over the volume element.                     ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the residuals for this volume element. */
    su2double *res = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Call the general function to carry out the matrix product. */
    MatrixProduct(nDOFs, nVar, nInt*nDim, matDerBasisIntTrans, fluxes, res);
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solInt);
  mkl_free(fluxes);
#endif
}

void CFEM_DG_EulerSolver::External_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                                  CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. ---*/
  su2double *solIntL, *solIntR, *fluxes;

#ifdef HAVE_MKL
  solIntL = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  fluxes  = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  solIntL = helpSolIntL.data();
  solIntR = helpSolIntR.data();
  fluxes  = helpFluxes.data();
#endif

  /*--------------------------------------------------------------------------*/
  /*--- Part 1: Compute the contribution to the contour integral in the    ---*/
  /*---         FEM formulation from the matching internal faces.          ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the internal matching faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nMatchingInternalFaces; ++l) {

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
    MatrixProduct(nDOFsFace0, nVar, nInt, basisFaceTrans, fluxes, resFace0);

    /* Check if the element to the right is an owned element. Only then
       the residual needs to be computed. */
    if(matchingInternalFaces[l].elemID1 < nVolElemOwned) {

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
        MatrixProduct(nDOFsFace1, nVar, nInt, basisFaceTrans, fluxes, resFace1);

        for(unsigned short i=0; i<(nVar*nDOFsFace1); ++i)
          resFace1[i] = -resFace1[i];
      }
    }
  }
  
  /*--------------------------------------------------------------------------*/
  /*--- Part 2: Accumulate the residuals for the owned DOFs and multiply   ---*/
  /*---         by the inverse of either the mass matrix or the lumped     ---*/
  /*---         mass matrix.                                               ---*/
  /*--------------------------------------------------------------------------*/

  /* Call the function CreateFinalResidual to carry out the task. The array
     fluxes is passed as temporary storage. */
  CreateFinalResidual(fluxes);

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(fluxes);
#endif
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
        such that the function MatrixProduct can be used to compute the left
        states in the integration points of the face, i.e. side 0. ---*/
  const unsigned long *DOFs = internalFace->DOFsSolFaceSide0;
  for(unsigned short i=0; i<nDOFsFace0; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
    su2double       *sol    = solFace + nVar*i;
    for(unsigned short j=0; j<nVar; ++j)
      sol[j] = solDOF[j];
  }

  /* Compute the left states. Call the general function to
     carry out the matrix product. */
  MatrixProduct(nInt, nVar, nDOFsFace0, basisFace0, solFace, solIntL);

  /*------------------------------------------------------------------------*/
  /*--- Step 2: Interpolate the right state in the integration points of ---*/
  /*---         the face.                                                ---*/
  /*------------------------------------------------------------------------*/

  /* Get the required information from the corresponding standard face. */
  const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();
  const su2double     *basisFace1 = standardMatchingFacesSol[ind].GetBasisFaceIntegrationSide1();

  /*--- Store the solution of the DOFs of side 1 of the face in contiguous memory
        such that the function MatrixProduct can be used to compute the right
        states in the integration points of the face, i.e. side 1. ---*/
  DOFs = internalFace->DOFsSolFaceSide1;
  for(unsigned short i=0; i<nDOFsFace1; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
    su2double       *sol    = solFace + nVar*i;
    for(unsigned short j=0; j<nVar; ++j)
      sol[j] = solDOF[j];
  }

  /* Compute the right states. Call the general function to
     carry out the matrix product. */
  MatrixProduct(nInt, nVar, nDOFsFace1, basisFace1, solFace, solIntR);

  /*------------------------------------------------------------------------*/
  /*--- Step 3: Compute the fluxes in the integration points using the   ---*/
  /*---         approximate Riemann solver.                              ---*/
  /*------------------------------------------------------------------------*/

  /* General function to compute the fluxes in the integration points. */
  ComputeInviscidFluxesFace(config, nInt, internalFace->metricNormalsFace,
                            solIntL, solIntR, fluxes, numerics);
}

void CFEM_DG_EulerSolver::CreateFinalResidual(su2double *tmpRes) {

  /* Loop over the owned volume elements. */
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Easier storage of the residuals for this volume element. */
    su2double *res = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /*--- Loop over the DOFs of the element and accumulate the residuals.
          This accumulation is carried inside the loop over the owned volume
          elements, such that an OpenMP parallelization of the loop over the
          volume elements is straightforward. ---*/
    for(unsigned short i=0; i<volElem[l].nDOFsSol; ++i) {

      su2double *resDOF      = res + nVar*i;
      const unsigned long ii = volElem[l].offsetDOFsSolLocal + i;
      for(unsigned long j=nEntriesResFaces[ii]; j<nEntriesResFaces[ii+1]; ++j) {
        const su2double *resFace = VecResFaces.data() + nVar*entriesResFaces[j];
        for(unsigned short k=0; k<nVar; ++k)
          resDOF[k] += resFace[k];
      }
    }

    /* Check whether a multiplication must be carried out with the inverse of
       the lumped mass matrix or the full mass matrix. Note that it is crucial
       that the test is performed with the pointer lumpedMassMatrix and not
       with massMatrix. The reason is that for implicit time stepping schemes
       both arrays are in use. */
    if( volElem[l].lumpedMassMatrix ) {

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
      MatrixProduct(volElem[l].nDOFsSol, nVar, volElem[l].nDOFsSol,
                    volElem[l].massMatrix, tmpRes, res);
    }
  }
}

void CFEM_DG_EulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {

  /*--- Allocate the memory for the storage of the solution in the DOFs
        and in the integration points. Note that when the MKL library is used
        a special allocation is used to optimize performance. ---*/
  su2double *solIntL, *solDOFs;

#ifdef HAVE_MKL
  solIntL = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solDOFs = (su2double *) mkl_malloc(nDOFsMax*nVar*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolDOFs(nDOFsMax*nVar);
  solIntL = helpSolIntL.data();
  solDOFs = helpSolDOFs.data();
#endif

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
  Total_CDrag = 0.0; Total_CLift = 0.0; Total_CSideForce = 0.0; Total_CEff = 0.0;
  Total_CMx = 0.0;   Total_CMy = 0.0;   Total_CMz = 0.0;
  Total_CFx = 0.0;   Total_CFy = 0.0;   Total_CFz = 0.0;

  AllBound_CDrag_Inv = 0.0;  AllBound_CLift_Inv = 0.0; AllBound_CSideForce_Inv = 0.0;
  AllBound_CMx_Inv = 0.0;    AllBound_CMy_Inv = 0.0;   AllBound_CMz_Inv = 0.0;
  AllBound_CFx_Inv = 0.0;    AllBound_CFy_Inv = 0.0;   AllBound_CFz_Inv = 0.0;
  AllBound_CEff_Inv = 0.0;

  for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring(); ++iMarker_Monitoring) {
    Surface_CLift_Inv[iMarker_Monitoring]      = 0.0; Surface_CDrag_Inv[iMarker_Monitoring] = 0.0;
    Surface_CSideForce_Inv[iMarker_Monitoring] = 0.0; Surface_CEff_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0; Surface_CFy_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0; Surface_CMx_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0; Surface_CMz_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CLift[iMarker_Monitoring]          = 0.0; Surface_CDrag[iMarker_Monitoring]     = 0.0;
    Surface_CSideForce[iMarker_Monitoring]     = 0.0; Surface_CEff[iMarker_Monitoring]      = 0.0;
    Surface_CFx[iMarker_Monitoring]            = 0.0; Surface_CFy[iMarker_Monitoring]       = 0.0;
    Surface_CFz[iMarker_Monitoring]            = 0.0; Surface_CMx[iMarker_Monitoring]       = 0.0;
    Surface_CMy[iMarker_Monitoring]            = 0.0; Surface_CMz[iMarker_Monitoring]       = 0.0;
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
        string Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
        string Marker_Tag     = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }

      /* Check for a boundary for which the forces must be computed. */
      if((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
         (Boundary == ISOTHERMAL)) {

        /*--- Initialization for this marker ---*/
        CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
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

          /* Compute the left states in the integration points of the face. */
          LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], solDOFs, solIntL);

          /*--- Get the number of integration points and the integration
                weights from the corresponding standard element. ---*/
          const unsigned short ind   = surfElem[l].indStandardElement;
          const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
          const su2double *weights   = standardBoundaryFacesSol[ind].GetWeightsIntegration();

          /* Loop over the integration points of this surface element. */
          for(unsigned short i=0; i<nInt; ++i) {

            /* Easier storage of the solution, the normals and the coordinates
               for this integration point. */
            const su2double *sol     = solIntL + i*nVar;
            const su2double *normals = surfElem[l].metricNormalsFace + i*(nDim+1);
            const su2double *Coord   = surfElem[l].coorIntegrationPoints + i*nDim;

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
            su2double MomentDist[3], Force[3];
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
          CDrag_Inv[iMarker]  =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
          CLift_Inv[iMarker]  = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
          CEff_Inv[iMarker]   =  CLift_Inv[iMarker] / (CDrag_Inv[iMarker]+EPS);
          CMz_Inv[iMarker]    =  MomentInviscid[2];
          CFx_Inv[iMarker]    =  ForceInviscid[0];
          CFy_Inv[iMarker]    =  ForceInviscid[1];
        }
        if(nDim == 3) {
          CDrag_Inv[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta)
                                  +  ForceInviscid[1]*sin(Beta)
                                  +  ForceInviscid[2]*sin(Alpha)*cos(Beta);
          CLift_Inv[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
          CSideForce_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha)
                                  +  ForceInviscid[1]*cos(Beta)
                                  -  ForceInviscid[2]*sin(Beta)*sin(Alpha);
          CEff_Inv[iMarker]       =  CLift_Inv[iMarker] / (CDrag_Inv[iMarker] + EPS);
          CMx_Inv[iMarker]        =  MomentInviscid[0];
          CMy_Inv[iMarker]        =  MomentInviscid[1];
          CMz_Inv[iMarker]        =  MomentInviscid[2];
          CFx_Inv[iMarker]        =  ForceInviscid[0];
          CFy_Inv[iMarker]        =  ForceInviscid[1];
          CFz_Inv[iMarker]        =  ForceInviscid[2];
        }

        AllBound_CDrag_Inv        += CDrag_Inv[iMarker];
        AllBound_CLift_Inv        += CLift_Inv[iMarker];
        AllBound_CSideForce_Inv   += CSideForce_Inv[iMarker];
        AllBound_CEff_Inv          = AllBound_CLift_Inv / (AllBound_CDrag_Inv + EPS);
        AllBound_CMx_Inv          += CMx_Inv[iMarker];
        AllBound_CMy_Inv          += CMy_Inv[iMarker];
        AllBound_CMz_Inv          += CMz_Inv[iMarker];
        AllBound_CFx_Inv          += CFx_Inv[iMarker];
        AllBound_CFy_Inv          += CFy_Inv[iMarker];
        AllBound_CFz_Inv          += CFz_Inv[iMarker];

        /*--- Compute the coefficients per surface ---*/
        for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                         ++iMarker_Monitoring) {
          string Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
          string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CLift_Inv[iMarker_Monitoring]      += CLift_Inv[iMarker];
            Surface_CDrag_Inv[iMarker_Monitoring]      += CDrag_Inv[iMarker];
            Surface_CSideForce_Inv[iMarker_Monitoring] += CSideForce_Inv[iMarker];
            Surface_CEff_Inv[iMarker_Monitoring]        = Surface_CLift_Inv[iMarker_Monitoring]
                                                        / (Surface_CDrag_Inv[iMarker_Monitoring] + EPS);
            Surface_CFx_Inv[iMarker_Monitoring]        += CFx_Inv[iMarker];
            Surface_CFy_Inv[iMarker_Monitoring]        += CFy_Inv[iMarker];
            Surface_CFz_Inv[iMarker_Monitoring]        += CFz_Inv[iMarker];
            Surface_CMx_Inv[iMarker_Monitoring]        += CMx_Inv[iMarker];
            Surface_CMy_Inv[iMarker_Monitoring]        += CMy_Inv[iMarker];
            Surface_CMz_Inv[iMarker_Monitoring]        += CMz_Inv[iMarker];
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
  locBuf[ii++] = AllBound_CDrag_Inv;      locBuf[ii++] = AllBound_CLift_Inv;
  locBuf[ii++] = AllBound_CSideForce_Inv; locBuf[ii++] = AllBound_CMx_Inv;
  locBuf[ii++] = AllBound_CMy_Inv;        locBuf[ii++] = AllBound_CMz_Inv;
  locBuf[ii++] = AllBound_CFx_Inv;        locBuf[ii++] = AllBound_CFy_Inv;
  locBuf[ii++] = AllBound_CFz_Inv;

  for(unsigned short i=0; i<config->GetnMarker_Monitoring(); ++i) {
    locBuf[ii++] = Surface_CLift_Inv[i];      locBuf[ii++] = Surface_CDrag_Inv[i];
    locBuf[ii++] = Surface_CSideForce_Inv[i]; locBuf[ii++] = Surface_CFx_Inv[i];
    locBuf[ii++] = Surface_CFy_Inv[i];        locBuf[ii++] = Surface_CFz_Inv[i];
    locBuf[ii++] = Surface_CMx_Inv[i];        locBuf[ii++] = Surface_CMy_Inv[i];
    locBuf[ii++] = Surface_CMz_Inv[i];
  }

  /* Sum up all the data from all ranks. The result will be available on all ranks. */
  SU2_MPI::Allreduce(locBuf.data(), globBuf.data(), nCommSize, MPI_DOUBLE,
                     MPI_SUM, MPI_COMM_WORLD);

  /*--- Copy the data back from globBuf into the required variables. ---*/
  ii = 0;
  AllBound_CDrag_Inv      = globBuf[ii++]; AllBound_CLift_Inv = globBuf[ii++];
  AllBound_CSideForce_Inv = globBuf[ii++]; AllBound_CMx_Inv   = globBuf[ii++];
  AllBound_CMy_Inv        = globBuf[ii++]; AllBound_CMz_Inv   = globBuf[ii++]; 
  AllBound_CFx_Inv        = globBuf[ii++]; AllBound_CFy_Inv   = globBuf[ii++];
  AllBound_CFz_Inv        = globBuf[ii++];

  AllBound_CEff_Inv = AllBound_CLift_Inv/(AllBound_CDrag_Inv + EPS);

  for(unsigned short i=0; i<config->GetnMarker_Monitoring(); ++i) {
    Surface_CLift_Inv[i]      = globBuf[ii++]; Surface_CDrag_Inv[i] = globBuf[ii++];
    Surface_CSideForce_Inv[i] = globBuf[ii++]; Surface_CFx_Inv[i]   = globBuf[ii++];
    Surface_CFy_Inv[i]        = globBuf[ii++]; Surface_CFz_Inv[i]   = globBuf[ii++];
    Surface_CMx_Inv[i]        = globBuf[ii++]; Surface_CMy_Inv[i]   = globBuf[ii++];
    Surface_CMz_Inv[i]        = globBuf[ii++];

    Surface_CEff_Inv[i] = Surface_CLift_Inv[i]/(Surface_CDrag_Inv[i] + EPS);
  }

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
  Total_CDrag      = AllBound_CDrag_Inv;
  Total_CLift      = AllBound_CLift_Inv;
  Total_CSideForce = AllBound_CSideForce_Inv;
  Total_CEff       = Total_CLift / (Total_CDrag + EPS);
  Total_CMx        = AllBound_CMx_Inv;
  Total_CMy        = AllBound_CMy_Inv;
  Total_CMz        = AllBound_CMz_Inv;
  Total_CFx        = AllBound_CFx_Inv;
  Total_CFy        = AllBound_CFy_Inv;
  Total_CFz        = AllBound_CFz_Inv;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                   ++iMarker_Monitoring) {
    Surface_CLift[iMarker_Monitoring]      = Surface_CLift_Inv[iMarker_Monitoring];
    Surface_CDrag[iMarker_Monitoring]      = Surface_CDrag_Inv[iMarker_Monitoring];
    Surface_CSideForce[iMarker_Monitoring] = Surface_CSideForce_Inv[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CLift_Inv[iMarker_Monitoring] 
                                           / (Surface_CDrag_Inv[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solDOFs);
#endif
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
    const su2double *solDOFSOld = VecSolDOFsOld.data() + offset;
    su2double *solDOFs          = VecSolDOFs.data()    + offset;

    /* Loop over the DOFs for this element and update the solution and the L2 norm. */
    const su2double tmp = RK_AlphaCoeff*VecDeltaTime[l];

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      for(unsigned short iVar=0; iVar<nVar; ++iVar, ++i) {
        solDOFs[i] = solDOFSOld[i] - tmp*res[i];

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
  SU2_MPI::Allreduce(Residual_RMS, rbuf.data(), nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(unsigned short iVar=0; iVar<nVar; ++iVar) {

    if (rbuf[iVar] != rbuf[iVar]) {
      if (rank == MASTER_NODE)
        cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    SetRes_RMS(iVar, max(EPS*EPS, sqrt(rbuf[iVar]/nDOFsGlobal)));
  }

  /*--- Communicate the information about the maximum residual. ---*/

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

void CFEM_DG_EulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                        CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. ---*/
  su2double *solIntL, *solIntR, *fluxes;

#ifdef HAVE_MKL
  solIntL = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  fluxes  = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  solIntL = helpSolIntL.data();
  solIntR = helpSolIntR.data();
  fluxes  = helpFluxes.data();
#endif

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], fluxes, solIntL);

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
      const su2double *normals = surfElem[l].metricNormalsFace + i*(nDim+1);

      /* Compute the normal component of the momentum variables. */
      su2double rVn = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        rVn += UL[iDim+1]*normals[iDim];

      /* If the normal velocity must be mirrored instead of set to zero,
         the normal component that must be subtracted must be doubled. If the
         normal velocity must be set to zero, simply comment this line. */
      rVn *= 2.0;

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

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(fluxes);
#endif
}

void CFEM_DG_EulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. ---*/
  su2double *solIntL, *solIntR, *fluxes;

#ifdef HAVE_MKL
  solIntL = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  fluxes  = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  solIntL = helpSolIntL.data();
  solIntR = helpSolIntR.data();
  fluxes  = helpFluxes.data();
#endif

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       Use fluxes as a temporary storage array. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], fluxes, solIntL);

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

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(fluxes);
#endif
}

void CFEM_DG_EulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler wall boundary condition.  ---*/
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
}

void CFEM_DG_EulerSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  cout << "CFEM_DG_EulerSolver::BC_Custom: Not implemented yet" << endl;
  exit(1);
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

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Compute the fluxes in the integration points using the   ---*/
  /*---         approximate Riemann solver.                              ---*/
  /*------------------------------------------------------------------------*/

  /* General function to compute the fluxes in the integration points. */
  ComputeInviscidFluxesFace(config, nInt, surfElem->metricNormalsFace,
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
  MatrixProduct(nDOFs, nVar, nInt, basisFaceTrans, fluxes, resFace);
}

void CFEM_DG_EulerSolver::LeftStatesIntegrationPointsBoundaryFace(const CSurfaceElementFEM *surfElem,
                                                                  su2double *solFace, su2double *solIntL) {

  /* Get the required information from the corresponding standard face. */
  const unsigned short ind   = surfElem->indStandardElement;
  const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();
  const su2double *basisFace = standardBoundaryFacesSol[ind].GetBasisFaceIntegration();

  /* Easier storage of the DOFs of the face. */
  const unsigned long *DOFs = surfElem->DOFsSolFace;

  /* Copy the solution of the DOFs of the face such that it is contigious
     in memory. */
  for(unsigned short i=0; i<nDOFs; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
    su2double       *sol    = solFace + nVar*i;
    for(unsigned short j=0; j<nVar; ++j)
      sol[j] = solDOF[j];
  }

  /* Call the general function to carry out the matrix product. */
  MatrixProduct(nInt, nVar, nDOFs, basisFace, solFace, solIntL);
}

void CFEM_DG_EulerSolver::ComputeInviscidFluxesFace(CConfig             *config,
                                                    const unsigned long nPoints,
                                                    const su2double     *normalsFace,
                                                    const su2double     *solL,
                                                    const su2double     *solR,
                                                    su2double           *fluxes,
                                                    CNumerics           *numerics) {

  /* Easier storage of the specific heat ratio. */
  const su2double gm1   = Gamma - 1.0;
  
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

void CFEM_DG_EulerSolver::MatrixProduct(const int M,        const int N,        const int K,
                                        const su2double *A, const su2double *B, su2double *C) {

#ifdef HAVE_LIBXSMM

  /* The gemm function of libxsmm is used to carry out the multiplication.
     Note that libxsmm_gemm expects the matrices in column major order. That's
     why the calling sequence is different from cblas_dgemm. */
  libxsmm_gemm(NULL, NULL, N, M, K, NULL, B, NULL, A, NULL, NULL, C, NULL);

#elif defined (HAVE_CBLAS) || defined(HAVE_MKL)

  /* The standard blas routine dgemm is used for the multiplication. */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K,
              1.0, A, K, B, N, 0.0, C, N);
#else

  /* Standard internal implementation of the matrix matrix multiplication. */
  for(int i=0; i<M; ++i) {
    const int jj = i*K;
    const int kk = i*N;
    for(int j=0; j<N; ++j) {
      const int ii = kk + j;
      C[ii] = 0.0;
      for(int k=0; k<K; ++k)
        C[ii] += A[jj+k]*B[k*nVar+j];
    }
  }

#endif
}

CFEM_DG_NSSolver::CFEM_DG_NSSolver(void) : CFEM_DG_EulerSolver() {
  
  /*--- Basic array initialization ---*/
  
  CDrag_Visc = NULL; CLift_Visc = NULL; CSideForce_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  
  ForceViscous = NULL; MomentViscous = NULL; CSkinFriction = NULL;
  
  /*--- Surface-based array initialization ---*/
  
  Surface_CLift_Visc = NULL; Surface_CDrag_Visc = NULL; Surface_CSideForce_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  
}

CFEM_DG_NSSolver::CFEM_DG_NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
 : CFEM_DG_EulerSolver(geometry, config, iMesh) {
  
  /*--- Array initialization ---*/
  
  CDrag_Visc = NULL; CLift_Visc = NULL; CSideForce_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  
  Surface_CLift_Visc = NULL; Surface_CDrag_Visc = NULL; Surface_CSideForce_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  
  ForceViscous = NULL; MomentViscous = NULL;
  CSkinFriction = NULL;    Cauchy_Serie = NULL;

  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  //LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  //LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Non dimensional coefficients ---*/

  ForceViscous     = new su2double[3];
  MomentViscous    = new su2double[3];
  CDrag_Visc       = new su2double[nMarker];
  CLift_Visc       = new su2double[nMarker];
  CSideForce_Visc  = new su2double[nMarker];
  CMx_Visc         = new su2double[nMarker];
  CMy_Visc         = new su2double[nMarker];
  CMz_Visc         = new su2double[nMarker];
  CEff_Visc        = new su2double[nMarker];
  CFx_Visc         = new su2double[nMarker];
  CFy_Visc         = new su2double[nMarker];
  CFz_Visc         = new su2double[nMarker];
  
  Surface_CLift_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce_Visc = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Init total coefficients ---*/
  
  Total_CDrag   = 0.0;	Total_CLift        = 0.0;  Total_CSideForce   = 0.0;
  Total_CMx     = 0.0;	Total_CMy          = 0.0;  Total_CMz          = 0.0;
  Total_CEff    = 0.0;
  Total_CFx     = 0.0;	Total_CFy          = 0.0;  Total_CFz          = 0.0;
  
  /*--- Read farfield conditions from config ---*/
  
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  Prandtl_Lam   = config->GetPrandtl_Lam();
  Prandtl_Turb  = config->GetPrandtl_Turb();
  Tke_Inf       = config->GetTke_FreeStreamND();
}

CFEM_DG_NSSolver::~CFEM_DG_NSSolver(void) {
  unsigned short iMarker;
  
  if (CDrag_Visc != NULL)       delete [] CDrag_Visc;
  if (CLift_Visc != NULL)       delete [] CLift_Visc;
  if (CSideForce_Visc != NULL)  delete [] CSideForce_Visc;
  if (CMx_Visc != NULL)         delete [] CMx_Visc;
  if (CMy_Visc != NULL)         delete [] CMy_Visc;
  if (CMz_Visc != NULL)         delete [] CMz_Visc;
  if (CFx_Visc != NULL)         delete [] CFx_Visc;
  if (CFy_Visc != NULL)         delete [] CFy_Visc;
  if (CFz_Visc != NULL)         delete [] CFz_Visc;
  if (CEff_Visc != NULL)        delete [] CEff_Visc;
  if (ForceViscous != NULL)     delete [] ForceViscous;
  if (MomentViscous != NULL)    delete [] MomentViscous;
  
  
  if (Surface_CLift_Visc != NULL)      delete [] Surface_CLift_Visc;
  if (Surface_CDrag_Visc != NULL)      delete [] Surface_CDrag_Visc;
  if (Surface_CSideForce_Visc != NULL) delete [] Surface_CSideForce_Visc;
  if (Surface_CEff_Visc != NULL)       delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != NULL)        delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != NULL)        delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != NULL)        delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != NULL)        delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)        delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)        delete [] Surface_CMz_Visc;
  
  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;
  
  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }
  
}

void CFEM_DG_NSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  /*--- Collect the number of non-physical points for this iteration. ---*/
  CFEM_DG_EulerSolver::Preprocessing(geometry, solver_container, config, iMesh,
                                     iRKStep, RunTime_EqSystem, Output);
}

void CFEM_DG_NSSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

//if (rank == MASTER_NODE) cout << " Computing viscous forces." << endl;
}

void CFEM_DG_NSSolver::Internal_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  /* Constant factor present in the heat flux vector. */
  const su2double factHeatFlux = Gamma/Prandtl_Lam;

  /*--- Allocate the memory for some temporary storage needed to compute the
        internal residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. ---*/
  su2double *solAndGradInt, *fluxes;

#ifdef HAVE_MKL
  solAndGradInt = (su2double *) mkl_malloc(nIntegrationMax*nVar*(nDim+1)*sizeof(su2double), 64);
  fluxes        = (su2double *) mkl_malloc(nIntegrationMax*nVar*nDim    *sizeof(su2double), 64);
#else
  vector<su2double> helpSolInt(nIntegrationMax*nVar*(nDim+1));
  vector<su2double> helpFluxes(nIntegrationMax*nVar*nDim);
  solAndGradInt = helpSolInt.data();
  fluxes        = helpFluxes.data();
#endif

  /* Store the number of metric points per integration point, which depends
     on the number of dimensions. */
  const unsigned short nMetricPerPoint = nDim*nDim + 1;

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
    /*--- Step 1: Determine the solution variables and their gradients     ---*/
    /*---         w.r.t. the parametric coordinates in the integration     ---*/
    /*---         points of the element.                                   ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the solution variables for this element. */
    su2double *solDOFs = VecSolDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Call the general function to carry out the matrix product. */
    MatrixProduct(nInt*(nDim+1), nVar, nDOFs, matBasisInt, solDOFs, solAndGradInt);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the total fluxes (inviscid fluxes minus the      ---*/
    /*---         viscous fluxes), multiplied by minus the integration     ---*/
    /*---         weight, in the integration points.                       ---*/
    /*------------------------------------------------------------------------*/

    /* Determine the offset between the solution variables and the r-derivatives,
       which is also the offset between the r- and s-derivatives and the offset
       between s- and t-derivatives. */
    const unsigned short offDeriv = nVar*nInt;

    /* Loop over the integration points, ll is the counter for the fluxes
       in the integration points. */
    unsigned short ll = 0;
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the metric terms in this integration point. First
         compute the inverse of the Jacobian, the Jacobian is the first entry
         in the metric terms, and afterwards update the metric terms by 1. */
      const su2double *metricTerms = volElem[l].metricTerms + i*nMetricPerPoint;
      const su2double Jac          = metricTerms[0];
      const su2double JacInv       = 1.0/Jac;
      metricTerms                 += 1;

      /* Easier storage of the location where the solution data of this
         integration point starts. */
      const su2double *sol = solAndGradInt + nVar*i;

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

      /*--- Compute the pressure and the laminar viscosity. ---*/
      FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
      const su2double Pressure  = FluidModel->GetPressure();
      const su2double Viscosity = FluidModel->GetLaminarViscosity();

      /*--- Set the value of the second viscosity and compute the divergence
            term in the viscous normal stresses. ---*/
      const su2double lambda     = -2.0*Viscosity/3.0;
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

      /*--- Compute the Cartesian flux vectors in this integration point. This
            is the sum of the inviscid flux and the negative of the viscous
            flux. ---*/
      su2double fluxCart[5][3];
      const su2double rH = sol[nDim+1] + Pressure;
      for(unsigned short k=0; k<nDim; ++k) {

        /* Density flux vector. */
        fluxCart[0][k] = sol[k+1];

        /* Momentum flux vector. */
        for(unsigned short j=0; j<nDim; ++j)
          fluxCart[j+1][k] = sol[j+1]*vel[k] - tauVis[j][k];
        fluxCart[k+1][k]  += Pressure;

        /* Energy flux vector. */
        fluxCart[nDim+1][k] = rH*vel[k]                                   // Inviscid part
                            - Viscosity*factHeatFlux*StaticEnergyGrad[k]; // Heat flux part
        for(unsigned short j=0; j<nDim; ++j)
          fluxCart[nDim+1][k] -= tauVis[j][k]*vel[j];    // Work of the viscous forces part
      }

      /*--- Loop over the number of dimensions to compute the fluxes in the
            direction of the parametric coordinates. ---*/
      for(unsigned short k=0; k<nDim; ++k) {

        /* Pointer to the metric terms for this direction. */
        const su2double *metric = metricTerms + k*nDim;

        /*--- Loop over the number of variables in the flux vector.
              Note that also the counter ll must be updated here. ---*/
        for(unsigned short j=0; j<nVar; ++j, ++ll) {

          /* Carry out the dot product of the Cartesian fluxes and the
             metric terms to obtain the correct flux vector in this
             parametric direction. */
          fluxes[ll] = 0.0;
          for(unsigned short iDim=0; iDim<nDim; ++iDim)
            fluxes[ll] += fluxCart[j][iDim]*metric[iDim];

          /* Multiply the flux by minus the integration weight to obtain the
             correct expression in the weak formulation. */
          fluxes[ll] *= -weights[i];
        }
      }
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
    /*---         integration over the volume element.                     ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the residuals for this volume element. */
    su2double *res = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Call the general function to carry out the matrix product. */
    MatrixProduct(nDOFs, nVar, nInt*nDim, matDerBasisIntTrans, fluxes, res);
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solAndGradInt);
  mkl_free(fluxes);
#endif
}

void CFEM_DG_NSSolver::External_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. Furthermore, note
        the size for fluxes and gradSolInt and the max function for the viscFluxes.
        This is because these arrays are also used as temporary storage for the other
        purposes. ---*/
  su2double *solIntL, *solIntR, *gradSolInt, *fluxes, *viscFluxes;
  su2double *viscosityIntL, *viscosityIntR;

  unsigned short sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, nDOFsMax);

  const unsigned short sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

#ifdef HAVE_MKL
  solIntL       = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR       = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  gradSolInt    = (su2double *) mkl_malloc(sizeGradSolInt*sizeof(su2double), 64);
  fluxes        = (su2double *) mkl_malloc(sizeFluxes*sizeof(su2double), 64);
  viscFluxes    = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
  viscosityIntL = (su2double *) mkl_malloc(nIntegrationMax*sizeof(su2double), 64);
  viscosityIntR = (su2double *) mkl_malloc(nIntegrationMax*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpGradSolInt(sizeGradSolInt);
  vector<su2double> helpFluxes(sizeFluxes);
  vector<su2double> helpViscFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  vector<su2double> helpViscosityIntL(nIntegrationMax);
  vector<su2double> helpViscosityIntR(nIntegrationMax);
  solIntL       = helpSolIntL.data();
  solIntR       = helpSolIntR.data();
  gradSolInt    = helpGradSolInt.data();
  fluxes        = helpFluxes.data();
  viscFluxes    = helpViscFluxes.data();
  viscosityIntL = helpViscosityIntL.data();
  viscosityIntR = helpViscosityIntR.data();
#endif

  /*--------------------------------------------------------------------------*/
  /*--- Part 1: Compute the contribution to the contour integral in the    ---*/
  /*---         FEM formulation from the matching internal faces.          ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the internal matching faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nMatchingInternalFaces; ++l) {

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Compute the inviscid fluxes in the integration points of ---*/
    /*---         this matching face.                                      ---*/
    /*------------------------------------------------------------------------*/

    /* Compute the inviscid fluxes in the integration points. */
    InviscidFluxesInternalMatchingFace(config, &matchingInternalFaces[l],
                                       solIntL, solIntR, fluxes, numerics);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the viscous fluxes in the integration points of  ---*/
    /*---         this matching face and subtract them from the already    ---*/
    /*---         computed inviscid fluxes.                                ---*/
    /*------------------------------------------------------------------------*/

    /*--- Get the information from the standard face to compute the viscous
          fluxes in the integration points on the left side, i.e. side 0. ---*/
    const unsigned short ind          = matchingInternalFaces[l].indStandardElement;
    const unsigned short nInt         = standardMatchingFacesSol[ind].GetNIntegration();
          unsigned short nDOFsElem    = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
    const su2double     *derBasisElem = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationSide0();

    /* Call the general function to compute the viscous flux in normal
       direction for side 0. */
    ViscousNormalFluxFace(nInt, nDOFsElem, 0.0, false, derBasisElem, solIntL,
                          matchingInternalFaces[l].DOFsSolElementSide0,
                          matchingInternalFaces[l].metricCoorDerivFace0,
                          matchingInternalFaces[l].metricNormalsFace,
                          gradSolInt, viscFluxes, viscosityIntL);

    /*--- Subtract half of the viscous fluxes from the inviscid fluxes. The
          factor 0.5 comes from the fact that the average of the viscous fluxes
          of side 0 and side 1 must be taken in the DG-FEM formulation. ---*/
    for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] -= 0.5*viscFluxes[j];

    /*--- Get the information from the standard face to compute the viscous
          fluxes in the integration points on the right side, i.e. side 1. ---*/
    nDOFsElem    = standardMatchingFacesSol[ind].GetNDOFsElemSide1();
    derBasisElem = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationSide1();

    /* Call the general function to compute the viscous flux in normal
       direction for side 1. */
    ViscousNormalFluxFace(nInt, nDOFsElem, 0.0, false, derBasisElem, solIntR,
                          matchingInternalFaces[l].DOFsSolElementSide1,
                          matchingInternalFaces[l].metricCoorDerivFace1,
                          matchingInternalFaces[l].metricNormalsFace,
                          gradSolInt, viscFluxes, viscosityIntR);

    /*--- Subtract half of the viscous fluxes from the inviscid fluxes. ---*/
    for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] -= 0.5*viscFluxes[j];

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the penalty terms in the integration points of   ---*/
    /*---         this matching face and them to the already stored        ---*/
    /*---         inviscid and viscous fluxes.                             ---*/
    /*------------------------------------------------------------------------*/

    /* Get the required constants needed for the penalty terms. */
    const su2double ConstPenFace = standardMatchingFacesSol[ind].GetPenaltyConstant();
    const su2double lenScale0    = volElem[matchingInternalFaces[l].elemID0].lenScale;
    const su2double lenScale1    = volElem[matchingInternalFaces[l].elemID1].lenScale;

    /* Call the function PenaltyTermsFluxFace to compute the actual penalty
       terms. Use the array viscFluxes as storage. */
    PenaltyTermsFluxFace(nInt, solIntL, solIntR, viscosityIntL, viscosityIntR,
                         ConstPenFace, lenScale0, lenScale1,
                         matchingInternalFaces[l].metricNormalsFace,
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
    MatrixProduct(nDOFsFace0, nVar, nInt, basisFaceTrans, fluxes, resFace0);

    /* Check if the element to the right is an owned element. Only then
       the residual needs to be computed. */
    if(matchingInternalFaces[l].elemID1 < nVolElemOwned) {

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
        MatrixProduct(nDOFsFace1, nVar, nInt, basisFaceTrans, fluxes, resFace1);

        for(unsigned short i=0; i<(nVar*nDOFsFace1); ++i)
          resFace1[i] = -resFace1[i];
      }
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 5: Compute the symmetrizing terms, if present, in the       ---*/
    /*---         integration points of this matching face.                ---*/
    /*------------------------------------------------------------------------*/

    if( symmetrizingTermsPresent ) {

      /* Compute the symmetrizing fluxes in the nDim directions. */
      SymmetrizingFluxesFace(nInt, solIntL, solIntR, viscosityIntL, viscosityIntR,
                             matchingInternalFaces[l].metricNormalsFace, fluxes);

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
      /*--- Step 6: Distribute the symmetrizing terms to the DOFs. Note that ---*/
      /*---         these terms must be distributed to all the DOFs of the   ---*/
      /*---         adjacent elements, not only to the DOFs of the face.     ---*/
      /*------------------------------------------------------------------------*/

      /* Get the element information of side 0 of the face. */
      const unsigned short nDOFsElem0    = standardMatchingFacesSol[ind].GetNDOFsElemSide0();
      const su2double *derBasisElemTrans = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationTransposeSide0();

      /*--- Create the Cartesian derivatives of the basis functions in the integration
            points. The array gradSolInt is used to store these derivatives. ---*/
      unsigned int ii = 0;
      for(unsigned short j=0; j<nDOFsElem0; ++j) {
        for(unsigned short i=0; i<nInt; ++i, ii+=nDim) {

          /* Easier storage of the derivatives of the basis function w.r.t. the
             parametric coordinates, the location where to store the Cartesian
             derivatives of the basis functions, and the metric terms in this
             integration point. */
          const su2double *derParam    = derBasisElemTrans + ii;
          const su2double *metricTerms = matchingInternalFaces[l].metricCoorDerivFace0 + i*nDim*nDim;
                su2double *derCar      = gradSolInt + ii;

          /*--- Loop over the dimensions to compute the Cartesian derivatives
                of the basis functions. ---*/
          for(unsigned short k=0; k<nDim; ++k) {
            derCar[k] = 0.0;
            for(unsigned short l=0; l<nDim; ++l)
              derCar[k] += derParam[l]*metricTerms[k+l*nDim];
          }
        }
      }

      /* Set the pointer where to store the current residual and update the
         counter indResFaces. */
      su2double *resElem0 = VecResFaces.data() + indResFaces*nVar;
      indResFaces        += nDOFsElem0;

      /* Call the general function to carry out the matrix product to compute
         the residual for side 0. */
      MatrixProduct(nDOFsElem0, nVar, nInt*nDim, gradSolInt, fluxes, resElem0);
    
      /* Check if the element to the right is an owned element. Only then
         the residual needs to be computed. */
      if(matchingInternalFaces[l].elemID1 < nVolElemOwned) {

        /* Get the element information of side 1 of the face. */
        const unsigned short nDOFsElem1 = standardMatchingFacesSol[ind].GetNDOFsElemSide1();
        derBasisElemTrans = standardMatchingFacesSol[ind].GetMatDerBasisElemIntegrationTransposeSide1();

        /*--- Create the Cartesian derivatives of the basis functions in the integration
              points. The array gradSolInt is used to store these derivatives. ---*/
        ii = 0;
        for(unsigned short j=0; j<nDOFsElem1; ++j) {
          for(unsigned short i=0; i<nInt; ++i, ii+=nDim) {

            /* Easier storage of the derivatives of the basis function w.r.t. the
               parametric coordinates, the location where to store the Cartesian
               derivatives of the basis functions, and the metric terms in this
               integration point. */
            const su2double *derParam    = derBasisElemTrans + ii;
            const su2double *metricTerms = matchingInternalFaces[l].metricCoorDerivFace1 + i*nDim*nDim;
                  su2double *derCar      = gradSolInt + ii;

            /*--- Loop over the dimensions to compute the Cartesian derivatives
                  of the basis functions. ---*/
            for(unsigned short k=0; k<nDim; ++k) {
              derCar[k] = 0.0;
              for(unsigned short l=0; l<nDim; ++l)
                derCar[k] += derParam[l]*metricTerms[k+l*nDim];
            }
          }
        }

        /* Set the pointer where to store the current residual and update the
           counter indResFaces. */
        su2double *resElem1 = VecResFaces.data() + indResFaces*nVar;
        indResFaces        += nDOFsElem1;

        /* Call the general function to carry out the matrix product to compute
           the residual for side 1. Note that the symmetrizing residual should not
           be negated, because two minus signs enter the formulation for side 1,
           which cancel each other. */
        MatrixProduct(nDOFsElem1, nVar, nInt*nDim, gradSolInt, fluxes, resElem1);
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Part 2: Accumulate the residuals for the owned DOFs and multiply   ---*/
  /*---         by the inverse of either the mass matrix or the lumped     ---*/
  /*---         mass matrix.                                               ---*/
  /*--------------------------------------------------------------------------*/

  /* Call the function CreateFinalResidual to carry out the task. The array
     fluxes is passed as temporary storage. */
  CreateFinalResidual(fluxes);

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(gradSolInt);
  mkl_free(fluxes);
  mkl_free(viscFluxes);
  mkl_free(viscosityIntL);
  mkl_free(viscosityIntR);
#endif
}

void CFEM_DG_NSSolver::ViscousNormalFluxFace(const unsigned short nInt,
                                             const unsigned short nDOFsElem,
                                             const su2double      Wall_HeatFlux,
                                             const bool           HeatFlux_Prescribed,
                                             const su2double      *derBasisElem,
                                             const su2double      *solInt,
                                             const unsigned long  *DOFsElem,
                                             const su2double      *metricCoorDerivFace,
                                             const su2double      *metricNormalsFace,
                                             su2double            *gradSolInt,
                                             su2double            *viscNormFluxes,
                                             su2double            *viscosityInt) {

  /* Constant factor present in the heat flux vector. Set it to zero if the heat
     flux is prescribed, such that no if statements are needed in the loop. */
  const su2double factHeatFlux = HeatFlux_Prescribed ? 0.0: Gamma/Prandtl_Lam;

  /* Set the value of the prescribed heat flux for the same reason. */
  const su2double HeatFlux = HeatFlux_Prescribed ? Wall_HeatFlux : 0.0;

  /* Set the pointer solElem to viscNormFluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the element and the fluxes to be computed. */
  su2double *solElem = viscNormFluxes;

  /*--- Store the solution of the DOFs of the adjacent element in contiguous
        memory such that the function MatrixProduct can be used to compute the
        gradients solution variables in the integration points of the face. ---*/
  for(unsigned short i=0; i<nDOFsElem; ++i) {
    const su2double *solDOF = VecSolDOFs.data() + nVar*DOFsElem[i];
    su2double       *sol    = solElem + nVar*i;
    for(unsigned short j=0; j<nVar; ++j)
      sol[j] = solDOF[j];
  }

  /* Compute the gradients in the integration points. Call the general function to
     carry out the matrix product. */
  MatrixProduct(nInt*nDim, nVar, nDOFsElem, derBasisElem, solElem, gradSolInt);

  /* Determine the offset between r- and -s-derivatives, which is also the
     offset between s- and t-derivatives. */
  const unsigned short offDeriv = nVar*nInt;

  /*--- Loop over the integration points to compute the viscous fluxes. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the metric terms needed to compute the Cartesian
       gradients in this integration point and the starting locations of
       the solution and the gradients, w.r.t. the parametric coordinates
       of this solution. */
    const su2double *metricTerms = metricCoorDerivFace + i*nDim*nDim;
    const su2double *sol         = solInt + nVar*i;
    const su2double *gradSol     = gradSolInt + nVar*i;

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

    /*--- Compute the laminar viscosity, which is stored for later use. ---*/
    FluidModel->SetTDState_rhoe(sol[0], StaticEnergy);
    const su2double Viscosity = FluidModel->GetLaminarViscosity();
    viscosityInt[i] = Viscosity;

    /*--- Set the value of the second viscosity and compute the divergence
          term in the viscous normal stresses. ---*/
    const su2double lambda     = -2.0*Viscosity/3.0;
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

    /*--- Compute the Cartesian viscous flux vectors in this
          integration point. ---*/
    su2double fluxCart[5][3];
    for(unsigned short k=0; k<nDim; ++k) {

      /* Momentum flux vector. */
      for(unsigned short j=0; j<nDim; ++j)
        fluxCart[j+1][k] = tauVis[j][k];

      /* Energy flux vector. */
      fluxCart[nDim+1][k] = Viscosity*factHeatFlux*StaticEnergyGrad[k]; // Heat flux part
      for(unsigned short j=0; j<nDim; ++j)
        fluxCart[nDim+1][k] += tauVis[j][k]*vel[j];    // Work of the viscous forces part
    }

    /*--- Set the pointer to the normals and the pointer to the
          viscous normal fluxes of this integration point.
          Initialize the fluxes to zero. ---*/
    const su2double *normal = metricNormalsFace + i*(nDim+1);
    su2double *normalFlux   = viscNormFluxes + i*nVar;
    for(unsigned short j=0; j<nVar; ++j) normalFlux[j] = 0.0;

    /* Set the energy flux to the prescribed heat flux. If the heat flux is
       not prescribed HeatFlux is equal to zero. */
    normalFlux[nVar-1] = normal[nDim]*HeatFlux;
      
    /*--- Loop over the number of dimensions to compute the fluxes in the
          direction of the parametric coordinates. ---*/
    for(unsigned short k=0; k<nDim; ++k) {
      const su2double unscaledNorm = normal[k]*normal[nDim];

      /*--- Loop over the variables to update the normal flux. Note
            that the loop starts at 1, because the mass conservation
            does not have a viscous contribution. ---*/
      for(unsigned short j=1; j<nVar; ++j) 
        normalFlux[j] += fluxCart[j][k]*unscaledNorm;
    }
  }
}

void CFEM_DG_NSSolver::PenaltyTermsFluxFace(const unsigned short nInt,
                                            const su2double      *solInt0,
                                            const su2double      *solInt1,
                                            const su2double      *viscosityInt0,
                                            const su2double      *viscosityInt1,
                                            const su2double      ConstPenFace,
                                            const su2double      lenScale0,
                                            const su2double      lenScale1,
                                            const su2double      *metricNormalsFace,
                                            su2double            *penaltyFluxes) {

  /* Constant factor present in the heat flux vector and the ratio of the
     second viscosity and the viscosity itself. */
  const su2double factHeatFlux =  Gamma/Prandtl_Lam;
  const su2double lambdaOverMu = -2.0/3.0;

  /* Compute the maximum value of the scaled spectral radius of the viscous
     operator. Scaled means divided by nu. Multiply it by ConstPenFace. */
  const su2double radOverNu = ConstPenFace*max(2.0+lambdaOverMu, max(1.0,factHeatFlux));

  /*--- Loop over the integration points to compute the penalty fluxes. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the variables for this integration point. */
    const su2double *sol0   = solInt0 + nVar*i;
    const su2double *sol1   = solInt1 + nVar*i;
    const su2double *normal = metricNormalsFace + i*(nDim+1);
    su2double       *flux   = penaltyFluxes + nVar*i;

    /* Compute the penalty parameter of this face as the maximum of the
       penalty parameter from both sides. Multiply by the area to obtain the
       correct expression. */
    const su2double pen0 = radOverNu*viscosityInt0[i]/(lenScale0*sol0[0]);
    const su2double pen1 = radOverNu*viscosityInt1[i]/(lenScale1*sol1[0]);

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
                                              const su2double      *metricNormalsFace,
                                              su2double            *symmFluxes) {

  /* Constant factor present in the heat flux vector and the ratio of the
     second viscosity and the viscosity itself. */
  const su2double factHeatFlux =  Gamma/Prandtl_Lam;
  const su2double lambdaOverMu = -2.0/3.0;

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

  /*--- Loop over the number of integration points of the face. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    /* Easier storage of the variables for this integration point. */
    const su2double *sol0   = solInt0 + nVar*i;
    const su2double *sol1   = solInt1 + nVar*i;
    const su2double *normal = metricNormalsFace + i*(nDim+1);
    su2double       *flux   = symmFluxes + nVar*i*nDim;

    /* Determine the difference in conservative variables. Multiply these
       differences by the length of the normal vector to obtain the correct
       dimensions for the symmetrizing fluxes. */
    su2double dSol[5];
    for(unsigned short j=0; j<nVar; ++j)
      dSol[j] = normal[nDim]*(sol0[j] - sol1[j]);

    /*--- Compute the terms that occur in the symmetrizing fluxes
          for state 0 and state 1. ---*/
    const su2double DensityInv0 = 1.0/sol0[0],          DensityInv1 = 1.0/sol1[0];
    const su2double Etot0 = DensityInv0*sol0[nDim+1],   Etot1 = DensityInv1*sol1[nDim+1];
    const su2double nu0 = DensityInv0*viscosityInt0[i], nu1 = DensityInv1*viscosityInt1[i];

    su2double velNorm0 = 0.0, velNorm1 = 0.0, velSquared0 = 0.0, velSquared1 = 0.0;
    su2double vel0[3], vel1[3];
    for(unsigned short j=0; j<nDim; ++j) {
      vel0[j] = DensityInv0*sol0[j+1];
      vel1[j] = DensityInv1*sol1[j+1];

      velNorm0 += vel0[j]*normal[j];
      velNorm1 += vel1[j]*normal[j];

      velSquared0 += vel0[j]*vel0[j];
      velSquared1 += vel1[j]*vel1[j];
    }

    /*--- Compute the average of the terms that occur in the symmetrizing
          fluxes. The average of the left and right terms is taken, rather
          than the terms evaluated at the average state, because the viscous
          fluxes are also computed as the average of the fluxes and not the
          fluxes of the averaged state. ---*/
    const su2double nuAvg = 0.5*(nu0 + nu1);
    const su2double nuVelSquaredAvg     = 0.5*(nu0*velSquared0 + nu1*velSquared1);
    const su2double nuVelNormAve        = 0.5*(nu0*velNorm0    + nu1*velNorm1);
    const su2double nuEminVelSquaredAve = 0.5*(nu0*(Etot0-velSquared0)
                                        +      nu1*(Etot1-velSquared1));

    su2double nuVelAvg[3], nuVelVelAvg[3];
    for(unsigned short j=0; j<nDim; ++j) {
      nuVelAvg[j]    = 0.5*(nu0*vel0[j]          + nu1*vel1[j]);
      nuVelVelAvg[j] = 0.5*(nu0*vel0[j]*velNorm0 + nu1*vel1[j]*velNorm1);
    }

    /*--- Abbreviations to make the flux computations a bit more efficient. ---*/
    su2double abv1 = 0.0, abv2 = 0.0;
    for(unsigned short j=0; j<nDim; ++j) {
      abv1 += normal[j]  *dSol[j+1];
      abv2 += nuVelAvg[j]*dSol[j+1];
    }

    const su2double abv3 = beta*(nuAvg*abv1 - nuVelNormAve*dSol[0]);
    const su2double abv4 = factHeatFlux*(nuAvg*dSol[nDim+1] - abv2
                         -               nuEminVelSquaredAve*dSol[0]) + abv2;

    /*--- Loop over the dimensions to compute the symmetrizing fluxes.
          ll is the counter for flux. ---*/
    unsigned short ll = 0;
    for(unsigned short k=0; k<nDim; ++k) {

      /* The symmetrizing density flux, which is zero. */
      flux[ll++] = 0.0;

      /* Loop over the dimensions to compute the symmetrizing momentum fluxes. */
      for(unsigned short j=0; j<nDim; ++j, ++ll) {

        /* Make a distinction between k == j and k != j and compute
           the momentum fluxes accordingly. */
        if(k == j) flux[ll] = abv3 + alphaP1*normal[j]*(nuAvg*dSol[j+1]
                                   -                    nuVelAvg[j]*dSol[0]);
        else       flux[ll] = nuAvg*(normal[k]*dSol[j+1] + alpha*normal[j]*dSol[k+1])
                            - (normal[k]*nuVelAvg[j] + alpha*normal[j]*nuVelAvg[k])*dSol[0];
      }

      /* The symmetrizing energy flux. */
      flux[ll++] = normal[k]*abv4
                 - (lambdaP1*nuVelVelAvg[k] + nuVelSquaredAvg*normal[k])*dSol[0]
                 + alpha*nuVelNormAve*dSol[k+1] + beta*nuVelAvg[k]*abv1;
    }
  }
}

void CFEM_DG_NSSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. Furthermore, note
        the size for fluxes and gradSolInt and the max function for the viscFluxes.
        This is because these arrays are also used as temporary storage for the other
        purposes. ---*/
  su2double *solIntL, *solIntR, *gradSolInt, *fluxes, *viscFluxes, *viscosityInt;

  unsigned short sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, nDOFsMax);

  const unsigned short sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

#ifdef HAVE_MKL
  solIntL      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  gradSolInt   = (su2double *) mkl_malloc(sizeGradSolInt*sizeof(su2double), 64);
  fluxes       = (su2double *) mkl_malloc(sizeFluxes*sizeof(su2double), 64);
  viscFluxes   = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
  viscosityInt = (su2double *) mkl_malloc(nIntegrationMax*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpGradSolInt(sizeGradSolInt);
  vector<su2double> helpFluxes(sizeFluxes);
  vector<su2double> helpViscFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  vector<su2double> helpViscosityInt(nIntegrationMax);
  solIntL      = helpSolIntL.data();
  solIntR      = helpSolIntR.data();
  gradSolInt   = helpGradSolInt.data();
  fluxes       = helpFluxes.data();
  viscFluxes   = helpViscFluxes.data();
  viscosityInt = helpViscosityInt.data();
#endif

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], fluxes, solIntL);

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
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], 0.0,
                                false, solIntL, solIntR, gradSolInt, fluxes,
                                viscFluxes, viscosityInt, resFaces, indResFaces);
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(gradSolInt);
  mkl_free(fluxes);
  mkl_free(viscFluxes);
  mkl_free(viscosityInt);
#endif
}

void CFEM_DG_NSSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler wall boundary condition.  ---*/
  CFEM_DG_NSSolver::BC_Euler_Wall(geometry, solver_container, conv_numerics,
                                  config, val_marker);
}

void CFEM_DG_NSSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. Furthermore, note
        the size for fluxes and gradSolInt and the max function for the viscFluxes.
        This is because these arrays are also used as temporary storage for the other
        purposes. ---*/
  su2double *solIntL, *solIntR, *gradSolInt, *fluxes, *viscFluxes, *viscosityInt;

  unsigned short sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, nDOFsMax);

  const unsigned short sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

#ifdef HAVE_MKL
  solIntL      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  gradSolInt   = (su2double *) mkl_malloc(sizeGradSolInt*sizeof(su2double), 64);
  fluxes       = (su2double *) mkl_malloc(sizeFluxes*sizeof(su2double), 64);
  viscFluxes   = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
  viscosityInt = (su2double *) mkl_malloc(nIntegrationMax*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpGradSolInt(sizeGradSolInt);
  vector<su2double> helpFluxes(sizeFluxes);
  vector<su2double> helpViscFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  vector<su2double> helpViscosityInt(nIntegrationMax);
  solIntL      = helpSolIntL.data();
  solIntR      = helpSolIntR.data();
  gradSolInt   = helpGradSolInt.data();
  fluxes       = helpFluxes.data();
  viscFluxes   = helpViscFluxes.data();
  viscosityInt = helpViscosityInt.data();
#endif

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], fluxes, solIntL);

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
      const su2double *normals = surfElem[l].metricNormalsFace + i*(nDim+1);

      /* Compute the normal component of the momentum variables. */
      su2double rVn = 0.0;
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        rVn += UL[iDim+1]*normals[iDim];

      /* If the normal velocity must be mirrored instead of set to zero,
         the normal component that must be subtracted must be doubled. If the
         normal velocity must be set to zero, simply comment this line. */
      rVn *= 2.0;

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
    ResidualViscousBoundaryFace(config, numerics, &surfElem[l], 0.0,
                                false, solIntL, solIntR, gradSolInt, fluxes,
                                viscFluxes, viscosityInt, resFaces, indResFaces);
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(gradSolInt);
  mkl_free(fluxes);
  mkl_free(viscFluxes);
  mkl_free(viscosityInt);
#endif
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

  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. Furthermore, note
        the size for fluxes and gradSolInt and the max function for the viscFluxes.
        This is because these arrays are also used as temporary storage for the other
        purposes. ---*/
  su2double *solIntL, *solIntR, *gradSolInt, *fluxes, *viscFluxes, *viscosityInt;

  unsigned short sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, nDOFsMax);

  const unsigned short sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

#ifdef HAVE_MKL
  solIntL      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  gradSolInt   = (su2double *) mkl_malloc(sizeGradSolInt*sizeof(su2double), 64);
  fluxes       = (su2double *) mkl_malloc(sizeFluxes*sizeof(su2double), 64);
  viscFluxes   = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
  viscosityInt = (su2double *) mkl_malloc(nIntegrationMax*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpGradSolInt(sizeGradSolInt);
  vector<su2double> helpFluxes(sizeFluxes);
  vector<su2double> helpViscFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  vector<su2double> helpViscosityInt(nIntegrationMax);
  solIntL      = helpSolIntL.data();
  solIntR      = helpSolIntR.data();
  gradSolInt   = helpGradSolInt.data();
  fluxes       = helpFluxes.data();
  viscFluxes   = helpViscFluxes.data();
  viscosityInt = helpViscosityInt.data();
#endif

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], fluxes, solIntL);

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

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], Wall_HeatFlux,
                                true, solIntL, solIntR, gradSolInt, fluxes,
                                viscFluxes, viscosityInt, resFaces, indResFaces);
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(gradSolInt);
  mkl_free(fluxes);
  mkl_free(viscFluxes);
  mkl_free(viscosityInt);
#endif
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
  
  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. Furthermore, note
        the size for fluxes and gradSolInt and the max function for the viscFluxes.
        This is because these arrays are also used as temporary storage for the other
        purposes. ---*/
  su2double *solIntL, *solIntR, *gradSolInt, *fluxes, *viscFluxes, *viscosityInt;

  unsigned short sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, nDOFsMax);

  const unsigned short sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

#ifdef HAVE_MKL
  solIntL      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  gradSolInt   = (su2double *) mkl_malloc(sizeGradSolInt*sizeof(su2double), 64);
  fluxes       = (su2double *) mkl_malloc(sizeFluxes*sizeof(su2double), 64);
  viscFluxes   = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
  viscosityInt = (su2double *) mkl_malloc(nIntegrationMax*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpGradSolInt(sizeGradSolInt);
  vector<su2double> helpFluxes(sizeFluxes);
  vector<su2double> helpViscFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  vector<su2double> helpViscosityInt(nIntegrationMax);
  solIntL      = helpSolIntL.data();
  solIntR      = helpSolIntR.data();
  gradSolInt   = helpGradSolInt.data();
  fluxes       = helpFluxes.data();
  viscFluxes   = helpViscFluxes.data();
  viscosityInt = helpViscosityInt.data();
#endif

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], fluxes, solIntL);

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

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, conv_numerics, &surfElem[l], 0.0,
                                false, solIntL, solIntR, gradSolInt, fluxes,
                                viscFluxes, viscosityInt, resFaces, indResFaces);
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(gradSolInt);
  mkl_free(fluxes);
  mkl_free(viscFluxes);
  mkl_free(viscosityInt);
#endif
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

  /*--- Allocate the memory for some temporary storage needed to compute the
        residual efficiently. Note that when the MKL library is used
        a special allocation is used to optimize performance. Furthermore, note
        the size for fluxes and gradSolInt and the max function for the viscFluxes.
        This is because these arrays are also used as temporary storage for the other
        purposes. ---*/
  su2double *solIntL, *solIntR, *gradSolInt, *fluxes, *viscFluxes, *viscosityInt;

  unsigned short sizeFluxes = nIntegrationMax*nDim;
  sizeFluxes = nVar*max(sizeFluxes, nDOFsMax);

  const unsigned short sizeGradSolInt = nIntegrationMax*nDim*max(nVar,nDOFsMax);

#ifdef HAVE_MKL
  solIntL      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  solIntR      = (su2double *) mkl_malloc(nIntegrationMax*nVar*sizeof(su2double), 64);
  gradSolInt   = (su2double *) mkl_malloc(sizeGradSolInt*sizeof(su2double), 64);
  fluxes       = (su2double *) mkl_malloc(sizeFluxes*sizeof(su2double), 64);
  viscFluxes   = (su2double *) mkl_malloc(max(nIntegrationMax,nDOFsMax)*nVar*sizeof(su2double), 64);
  viscosityInt = (su2double *) mkl_malloc(nIntegrationMax*sizeof(su2double), 64);
#else
  vector<su2double> helpSolIntL(nIntegrationMax*nVar);
  vector<su2double> helpSolIntR(nIntegrationMax*nVar);
  vector<su2double> helpGradSolInt(sizeGradSolInt);
  vector<su2double> helpFluxes(sizeFluxes);
  vector<su2double> helpViscFluxes(max(nIntegrationMax,nDOFsMax)*nVar);
  vector<su2double> helpViscosityInt(nIntegrationMax);
  solIntL      = helpSolIntL.data();
  solIntR      = helpSolIntR.data();
  gradSolInt   = helpGradSolInt.data();
  fluxes       = helpFluxes.data();
  viscFluxes   = helpViscFluxes.data();
  viscosityInt = helpViscosityInt.data();
#endif

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /* Compute the left states in the integration points of the face.
       The array fluxes is used as temporary storage inside the function
       LeftStatesIntegrationPointsBoundaryFace. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], fluxes, solIntL);

    /* Determine the number of integration points. */
    const unsigned short ind  = surfElem[l].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Loop over the integration points to compute the right state via the
          customized boundary conditions. ---*/
    for(unsigned short i=0; i<nInt; ++i) {

      /* Easier storage of the right solution for this integration point and
         the coordinates of this integration point. */
      const su2double *coor = surfElem[l].coorIntegrationPoints + i*nDim;
            su2double *UR   = solIntR + i*nVar;

#ifdef CUSTOM_BC_NSUNITQUAD

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

    /* The remainder of the contribution of this boundary face to the residual
       is the same for all boundary conditions. Hence a generic function can
       be used to carry out this task. */
    ResidualViscousBoundaryFace(config, numerics, &surfElem[l], 0.0,
                                false, solIntL, solIntR, gradSolInt, fluxes,
                                viscFluxes, viscosityInt, resFaces, indResFaces);
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(gradSolInt);
  mkl_free(fluxes);
  mkl_free(viscFluxes);
  mkl_free(viscosityInt);
#endif
}

void CFEM_DG_NSSolver::ResidualViscousBoundaryFace(
                                      CConfig                  *config,
                                      CNumerics                *conv_numerics,
                                      const CSurfaceElementFEM *surfElem,
                                      const su2double          Wall_HeatFlux,
                                      const bool               HeatFlux_Prescribed,
                                      const su2double          *solInt0,
                                      const su2double          *solInt1,
                                      su2double                *gradSolInt,
                                      su2double                *fluxes,
                                      su2double                *viscFluxes,
                                      su2double                *viscosityInt,
                                      su2double                *resFaces,
                                      unsigned long            &indResFaces) {

  /*--- Get the required information from the standard element. ---*/
  const unsigned short ind          = surfElem->indStandardElement;
  const unsigned short nInt         = standardBoundaryFacesSol[ind].GetNIntegration();
  const unsigned short nDOFs        = standardBoundaryFacesSol[ind].GetNDOFsFace();
  const unsigned short nDOFsElem    = standardBoundaryFacesSol[ind].GetNDOFsElem();
  const su2double      ConstPenFace = standardBoundaryFacesSol[ind].GetPenaltyConstant();
  const su2double     *derBasisElem = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegration();
  const su2double     *weights      = standardBoundaryFacesSol[ind].GetWeightsIntegration();

  /*------------------------------------------------------------------------*/
  /*--- Step 1: Compute the fluxes in the integration points.            ---*/
  /*------------------------------------------------------------------------*/

  /* General function to compute the inviscid fluxes, using an approximate
     Riemann solver in the integration points. */
  ComputeInviscidFluxesFace(config, nInt, surfElem->metricNormalsFace,
                            solInt0, solInt1, fluxes, conv_numerics);

  /* Call the general function to compute the viscous flux in normal
     direction for the face. */
  ViscousNormalFluxFace(nInt, nDOFsElem, Wall_HeatFlux, HeatFlux_Prescribed,
                        derBasisElem, solInt0, surfElem->DOFsSolElement,
                        surfElem->metricCoorDerivFace, surfElem->metricNormalsFace,
                        gradSolInt, viscFluxes, viscosityInt);

  /* Subtract the viscous fluxes from the inviscid fluxes. */
  for(unsigned short j=0; j<(nVar*nInt); ++j) fluxes[j] -= viscFluxes[j];

  /* Get the length scale for the adjacent element. */
  const su2double lenScale = volElem[surfElem->volElemID].lenScale;

  /* Call the function PenaltyTermsFluxFace to compute the actual penalty
     terms. Use the array viscFluxes as storage. */
  PenaltyTermsFluxFace(nInt, solInt0, solInt1, viscosityInt, viscosityInt,
                       ConstPenFace, lenScale, lenScale,
                       surfElem->metricNormalsFace, viscFluxes);

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
  MatrixProduct(nDOFs, nVar, nInt, basisFaceTrans, fluxes, resFace);

  /*------------------------------------------------------------------------*/
  /*--- Step 3: Compute the symmetrizing terms, if present, in the       ---*/
  /*---         integration points of this boundary face.                ---*/
  /*------------------------------------------------------------------------*/

  if( symmetrizingTermsPresent ) {

    /* Compute the symmetrizing fluxes in the nDim directions. */
    SymmetrizingFluxesFace(nInt, solInt0, solInt1, viscosityInt, viscosityInt,
                           surfElem->metricNormalsFace, fluxes);

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

    /* Get the correct form of the basis functions needed for the matrix
       multiplication to compute the residual. */
    const su2double *derBasisElemTrans = standardBoundaryFacesSol[ind].GetMatDerBasisElemIntegrationTranspose();

    /*--- Create the Cartesian derivatives of the basis functions in the integration
          points. The array gradSolInt is used to store these derivatives. ---*/
    unsigned int ii = 0;
    for(unsigned short j=0; j<nDOFsElem; ++j) {
      for(unsigned short i=0; i<nInt; ++i, ii+=nDim) {

        /* Easier storage of the derivatives of the basis function w.r.t. the
           parametric coordinates, the location where to store the Cartesian
           derivatives of the basis functions, and the metric terms in this
           integration point. */
        const su2double *derParam    = derBasisElemTrans + ii;
        const su2double *metricTerms = surfElem->metricCoorDerivFace + i*nDim*nDim;
              su2double *derCar      = gradSolInt + ii;

        /*--- Loop over the dimensions to compute the Cartesian derivatives
              of the basis functions. ---*/
        for(unsigned short k=0; k<nDim; ++k) {
          derCar[k] = 0.0;
          for(unsigned short l=0; l<nDim; ++l)
            derCar[k] += derParam[l]*metricTerms[k+l*nDim];
        }
      }
    }

    /* Call the general function to carry out the matrix product to compute
       the residual. */
    MatrixProduct(nDOFsElem, nVar, nInt*nDim, gradSolInt, fluxes, resElem);
  }
}
