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

  /*--- Determine the maximum number of DOFs used. This is for the volume elements. */
  nDOFsMax = 0;
  for(unsigned short i=0; i<nStandardElementsSol; ++i) {
    const unsigned short nDOFs = standardElementsSol[i].GetNDOFs();
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

  /* First the internal matching faces. */
  unsigned long sizeVecResFaces = 0;
  for(unsigned long i=0; i<nMatchingInternalFaces; ++i) {
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

        sizeVecResFaces += nDOFsFace;
        for(unsigned short j=0; j<nDOFsFace; ++j)
          ++nEntriesResFaces[surfElem[i].DOFsSolFace[j]+1];
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
        const unsigned short ind       = surfElem[i].indStandardElement;
        const unsigned short nDOFsFace = standardBoundaryFacesSol[ind].GetNDOFsFace();

        for(unsigned short j=0; j<nDOFsFace; ++j) {
          unsigned long jj    = counterEntries[surfElem[i].DOFsSolFace[j]]++;
          entriesResFaces[jj] = sizeVecResFaces++;
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
    su2double nBadDOFsLoc = nBadDOFs;
    SU2_MPI::Reduce(&nBadDOFsLoc, &nBadDOFs, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#endif

    if((rank == MASTER_NODE) && (nBadDOFs != 0))
      cout << "Warning. The original solution contains "<< nBadDOFs << " DOFs that are not physical." << endl;
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
  bool turbulent          = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
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

void CFEM_DG_EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /* Initialize the minimum and maximum time step. */
  Min_Delta_Time = 1.e25; Max_Delta_Time = 0.0;

  /* Easier storage of the CFL number. */
  const su2double CFL = config->GetCFL(iMesh);

  /*--- Check for a compressible solver. ---*/
  if(config->GetKind_Regime() == COMPRESSIBLE) {

    /*--- Loop over the owned volume elements. ---*/
    for(unsigned long i=0; i<nVolElemOwned; ++i) {

      /*--- Loop over the DOFs of this element and determine the maximum wave speed. ---*/
      su2double charVel2Max = 0.0;
      for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
        const su2double *solDOF = VecSolDOFs.data() + nVar*(volElem[i].offsetDOFsSolLocal + j);

        su2double velAbs[3];
        const su2double DensityInv = 1.0/solDOF[0];
        su2double Velocity2 = 0.0;
        for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
          const su2double vel = solDOF[iDim]*DensityInv;
          velAbs[iDim-1] = fabs(vel);
          Velocity2 += vel*vel;
        }

        su2double StaticEnergy = solDOF[nDim+1]*DensityInv - 0.5*Velocity2;
        FluidModel->SetTDState_rhoe(solDOF[0], StaticEnergy);
        su2double SoundSpeed2 = FluidModel->GetSoundSpeed2();

        su2double charVel2 = 0.0;
        for(unsigned short iDim=0; iDim<nDim; ++iDim)
          charVel2 += velAbs[iDim]*velAbs[iDim] + SoundSpeed2;

        charVel2Max = max(charVel2Max, charVel2);
      }

      /* Compute the time step for the element and update the minimum and
         maximum value. */
      VecDeltaTime[i] = CFL*volElem[i].lenScale/charVel2Max;

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
    /*---         the elements. This is a matrix matrix multiplication,    ---*/
    /*---         which can be carried out with either LIBXSMM, BLAS or    ---*/
    /*---         a standard internal implementation.                      ---*/ 
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

  /* Set the pointer solFace to fluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the face and the fluxes. */
  su2double *solFace = fluxes;

  /*--- Loop over the internal matching faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nMatchingInternalFaces; ++l) {
  
    /*------------------------------------------------------------------------*/
    /*--- Step 1: Interpolate the left and right states in the integration ---*/
    /*---         points of the face.                                      ---*/
    /*------------------------------------------------------------------------*/

    /* Get the required information from the corresponding standard face. */
    const unsigned short ind        = matchingInternalFaces[l].indStandardElement;
    const unsigned short nInt       = standardMatchingFacesSol[ind].GetNIntegration();
    const unsigned short nDOFsFace0 = standardMatchingFacesSol[ind].GetNDOFsFaceSide0();
    const unsigned short nDOFsFace1 = standardMatchingFacesSol[ind].GetNDOFsFaceSide1();

    const su2double *basisFace0 = standardMatchingFacesSol[ind].GetBasisFaceIntegrationSide0();
    const su2double *basisFace1 = standardMatchingFacesSol[ind].GetBasisFaceIntegrationSide1();
    const su2double *weights    = standardMatchingFacesSol[ind].GetWeightsIntegration();

    /* Compute the left states in the integration points of the face, i.e. side 0. */
    const unsigned long *DOFs = matchingInternalFaces[l].DOFsSolFaceSide0;
    for(unsigned short i=0; i<nDOFsFace0; ++i) {
      const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
      su2double       *sol    = solFace + nVar*i;
      for(unsigned short j=0; j<nVar; ++j)
        sol[j] = solDOF[j];
    }

    /* Call the general function to carry out the matrix product. */
    MatrixProduct(nInt, nVar, nDOFsFace0, basisFace0, solFace, solIntL);

    /* Compute the right states in the integration points of the face, i.e. side 1. */
    DOFs = matchingInternalFaces[l].DOFsSolFaceSide1;
    for(unsigned short i=0; i<nDOFsFace1; ++i) {
      const su2double *solDOF = VecSolDOFs.data() + nVar*DOFs[i];
      su2double       *sol    = solFace + nVar*i;
      for(unsigned short j=0; j<nVar; ++j)
        sol[j] = solDOF[j];
    }

    /* Call the general function to carry out the matrix product. */
    MatrixProduct(nInt, nVar, nDOFsFace1, basisFace1, solFace, solIntR);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the fluxes in the integration points using the   ---*/
    /*---         approximate Riemann solver.                              ---*/
    /*------------------------------------------------------------------------*/

    /* General function to compute the fluxes in the integration points. */
    ComputeInviscidFluxesFace(config, nInt, matchingInternalFaces[l].metricNormalsFace,
                              solIntL, solIntR, fluxes);

    /* Multiply the fluxes with the integration weight of the corresponding
       integration point. */
    for(unsigned short i=0; i<nInt; ++i) {
      su2double *flux = fluxes + i*nVar;

      for(unsigned short j=0; j<nVar; ++j)
        flux[j] *= weights[i];
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
    /*---         integration over this internal matching face.            ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the position in the residual array for side 0 of
       this face and update the corresponding counter. */
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
  /*---         mass matrix. The accumulation of the residuals is carried  ---*/
  /*---         inside the loop over the owned volume elements, such that  ---*/
  /*---         an OpenMP parallelization of this loop can be done.        ---*/
  /*--------------------------------------------------------------------------*/

  /* Loop over the owned volume elements. */
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Easier storage of the residuals for this volume element. */
    su2double *res = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /*--- Loop over the DOFs of the element and accumulate the residuals. ---*/
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
         Use the array fluxes as temporary storage. */
      memcpy(fluxes, res, nVar*volElem[l].nDOFsSol*sizeof(su2double));
      MatrixProduct(volElem[l].nDOFsSol, nVar, volElem[l].nDOFsSol,
                    volElem[l].massMatrix, fluxes, res);
    }
  }

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(fluxes);
#endif

}

void CFEM_DG_EulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE) cout << " Computing inviscid forces." << endl;
  
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
    const su2double *res = VecResDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;
    su2double *solDOFs   = VecSolDOFs.data() + nVar*volElem[l].offsetDOFsSolLocal;

    /* Loop over the DOFs for this element and update the solution and the L2 norm. */
    const su2double tmp = RK_AlphaCoeff*VecDeltaTime[l];

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      for(unsigned short iVar=0; iVar<nVar; ++iVar, ++i) {
        solDOFs[i] -= tmp*res[i];

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

  /* Set the pointer solFace to fluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the face and the fluxes. */
  su2double *solFace = fluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Interpolate the left states in the integration points of ---*/
    /*---         the face. From these left states, construct the right    ---*/
    /*---         states using the boundary conditions.                    ---*/
    /*------------------------------------------------------------------------*/

    /* Compute the left states in the integration points of the face. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], solFace, solIntL);

    /*--- Apply the inviscid wall boundary conditions to compute the right
          state. There are two options. Either the normal velocity is negated
          or the normal velocity is set to zero. Some experiments are needed
          to see which formulation gives better results. ---*/
    const unsigned short ind   = surfElem[l].indStandardElement;
    const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
    const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();
    const su2double *weights   = standardBoundaryFacesSol[ind].GetWeightsIntegration();

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

      /*--- Actually, the internal energy of UR is equal to UL. If the kinetic
            energy differs for UL and UR, the difference must be subtracted
            from the total energy of UR to obtain the correct value. ---*/
      su2double DensityInv = 1.0/UL[0];
      su2double diffKin    = 0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double velL = DensityInv*UL[iDim];
        const su2double velR = DensityInv*UR[iDim];
        diffKin += velL*velL - velR*velR;
      }

      UR[nDim+1] -= 0.5*UL[0]*diffKin;
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the fluxes in the integration points using the   ---*/
    /*---         approximate Riemann solver.                              ---*/
    /*------------------------------------------------------------------------*/

    /* General function to compute the fluxes in the integration points. */
    ComputeInviscidFluxesFace(config, nInt, surfElem[l].metricNormalsFace,
                              solIntL, solIntR, fluxes);

    /* Multiply the fluxes with the integration weight of the corresponding
       integration point. */
    for(unsigned short i=0; i<nInt; ++i) {
      su2double *flux = fluxes + i*nVar;

      for(unsigned short j=0; j<nVar; ++j)
        flux[j] *= weights[i];
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
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

  /* Set the pointer solFace to fluxes. This is just for readability, as the
     same memory can be used for the storage of the solution of the DOFs of
     the face and the fluxes. */
  su2double *solFace = fluxes;

  /* Set the starting position in the vector for the face residuals for
     this boundary marker. */
  su2double *resFaces = VecResFaces.data() + nVar*startLocResFacesMarkers[val_marker];

  /* Easier storage of the boundary faces for this boundary marker. */
  const unsigned long      nSurfElem = boundaries[val_marker].surfElem.size();
  const CSurfaceElementFEM *surfElem = boundaries[val_marker].surfElem.data();

  /*--- Loop over the boundary faces. ---*/
  unsigned long indResFaces = 0;
  for(unsigned long l=0; l<nSurfElem; ++l) {

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Interpolate the left states in the integration points of ---*/
    /*---         the face. From these left states, construct the right    ---*/
    /*---         states using the boundary conditions.                    ---*/
    /*------------------------------------------------------------------------*/

    /* Compute the left states in the integration points of the face. */
    LeftStatesIntegrationPointsBoundaryFace(&surfElem[l], solFace, solIntL);

    /*--- Apply the farfield boundary conditions to compute the right state. ---*/
    const unsigned short ind   = surfElem[l].indStandardElement;
    const unsigned short nInt  = standardBoundaryFacesSol[ind].GetNIntegration();
    const unsigned short nDOFs = standardBoundaryFacesSol[ind].GetNDOFsFace();
    const su2double *weights   = standardBoundaryFacesSol[ind].GetWeightsIntegration();

    for(unsigned short i=0; i<nInt; ++i) {
      su2double *UR = solIntR + i*nVar;
      for(unsigned short j=0; j<nVar; ++j)
        UR[j] = ConsVarFreeStream[j];
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the fluxes in the integration points using the   ---*/
    /*---         approximate Riemann solver.                              ---*/
    /*------------------------------------------------------------------------*/

    /* General function to compute the fluxes in the integration points. */
    ComputeInviscidFluxesFace(config, nInt, surfElem[l].metricNormalsFace,
                              solIntL, solIntR, fluxes);

    /* Multiply the fluxes with the integration weight of the corresponding
       integration point. */
    for(unsigned short i=0; i<nInt; ++i) {
      su2double *flux = fluxes + i*nVar;

      for(unsigned short j=0; j<nVar; ++j)
        flux[j] *= weights[i];
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
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

  /*--- If the MKL is used the temporary storage must be released again. ---*/
#ifdef HAVE_MKL
  mkl_free(solIntL);
  mkl_free(solIntR);
  mkl_free(fluxes);
#endif
}

void CFEM_DG_EulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler residual ---*/
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
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
                                                    su2double           *fluxes) {

  /*****************************************************************************/
  /* THIS IS A TEMPORARY IMPLEMENTATION, WHICH USES ROE'S APPROXIMATE RIEMANN  */
  /* SOLVER. FOR THE ACTUAL IMPLEMENTATION THE RIEMANN SOLVERS OF THE FINITE   */
  /* VOLUME SOLVER MUST BE USED. IN ORDER TO ACCOMPLISH THIS THE WAY OF        */
  /* CALLING THESE ROUTINES MUST LIKELY BE CHANGED.                            */
  /*****************************************************************************/

  /* Easier storage of the specific heat ratio. */
  const su2double gamma = config->GetGamma();
  const su2double gm1   = gamma - 1.0;

  /* Make a distinction between two and three space dimensions. */
  switch( nDim ) {

    case 2: {

      /* Two dimensional simulation. Loop over the number of points. */
      for(unsigned long i=0; i<nPoints; ++i) {
        
        /* Easier storage of the left and right solution, the face normals and
           the flux vector for this point. */
        const su2double *UL   = solL + i*nVar;
        const su2double *UR   = solR + i*nVar;
        const su2double *norm = normalsFace + i*(nDim+1);
              su2double *flux = fluxes + i*nVar;

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

        /*--- Compute the coefficient eta for the entropy correction. At the moment a
              1D entropy correction is used, which removes expansion shocks. ---*/
        const su2double aL  = sqrt(gamma*pL/UL[0]);
        const su2double aR  = sqrt(gamma*pR/UR[0]);
        const su2double eta = 0.5*(fabs((vxL - vxR)*nx + (vyL - vyR)*ny)
                            +      fabs(aL - aR));        

        /*--- Compute the absolute values of the three eigenvalues, apply the
              entropy correction and multiply by the area to obtain the correct
              values for the dissipation term. ---*/
        su2double lam1 = fabs(vnAvg + aAvg);
        su2double lam2 = fabs(vnAvg - aAvg);
        su2double lam3 = fabs(vnAvg);

        tmp = 2.0*eta;
        if(lam1 < tmp) lam1 = eta + 0.25*lam1*lam1/eta;
        if(lam2 < tmp) lam2 = eta + 0.25*lam2*lam2/eta;
        if(lam3 < tmp) lam3 = eta + 0.25*lam3*lam3/eta;

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

    case 3: {

      /* Three dimensional simulation. Loop over the number of points. */
      for(unsigned long i=0; i<nPoints; ++i) {

        /* Easier storage of the left and right solution, the face normals and
           the flux vector for this point. */
        const su2double *UL   = solL + i*nVar;
        const su2double *UR   = solR + i*nVar;
        const su2double *norm = normalsFace + i*(nDim+1);
              su2double *flux = fluxes + i*nVar;

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

        /*--- Compute the coefficient eta for the entropy correction. At the moment a
              1D entropy correction is used, which removes expansion shocks. ---*/
        const su2double aL  = sqrt(gamma*pL/UL[0]);
        const su2double aR  = sqrt(gamma*pR/UR[0]);
        const su2double eta = 0.5*(fabs((vxL-vxR)*nx + (vyL-vyR)*ny + (vzL-vzR)*nz)
                            +      fabs(aL-aR));

        /*--- Compute the absolute values of the three eigenvalues, apply the
              entropy correction and multiply by the area to obtain the correct
              values for the dissipation term. ---*/
        su2double lam1 = fabs(vnAvg + aAvg);
        su2double lam2 = fabs(vnAvg - aAvg);
        su2double lam3 = fabs(vnAvg);

        tmp = 2.0*eta;
        if(lam1 < tmp) lam1 = eta + 0.25*lam1*lam1/eta;
        if(lam2 < tmp) lam2 = eta + 0.25*lam2*lam2/eta;
        if(lam3 < tmp) lam3 = eta + 0.25*lam3*lam3/eta;

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
}

void CFEM_DG_EulerSolver::MatrixProduct(const int M,        const int N,        const int K,
                                        const su2double *A, const su2double *B, su2double *C) {
#ifdef HAVE_LIBXSMM

  /* The gemm function of libxsmm is used to carry out the multiplication. */
  libxsmm_gemm(NULL, NULL, M, N, K, NULL, A, NULL, B, NULL, NULL, C, NULL);

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
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double StrainMag = 0.0, Omega = 0.0, *Vorticity;
  
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  bool limiter_visc = config->GetViscous_Limiter_Flow();
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  
  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  
  StrainMag_Max = 0.0, Omega_Max = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    solver_container[FLOW_SOL]->node[iPoint]->SetVorticity(limiter_visc);
    solver_container[FLOW_SOL]->node[iPoint]->SetStrainMag(limiter_visc);
    
    StrainMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
    Vorticity = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
    Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);
    
    StrainMag_Max = max(StrainMag_Max, StrainMag);
    Omega_Max = max(Omega_Max, Omega);
    
  }
  
  /*--- Collect the number of non-physical points for this iteration. ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    su2double MyOmega_Max = Omega_Max; Omega_Max = 0.0;
    su2double MyStrainMag_Max = StrainMag_Max; StrainMag_Max = 0.0;
    
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    if (iMesh == MESH_0) {
      config->SetNonphysical_Points(ErrorCounter);
      solver_container[FLOW_SOL]->SetStrainMag_Max(StrainMag_Max);
      solver_container[FLOW_SOL]->SetOmega_Max(Omega_Max);
    }
    
  }
  
}

void CFEM_DG_NSSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE) cout << " Computing viscous forces." << endl;
}

void CFEM_DG_NSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
  
}

void CFEM_DG_NSSolver::Internal_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                                  CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Compute finite element convective residual and store here. ---*/
  
  if (rank == MASTER_NODE) cout << " Computing internal contributions to residual." << endl;
  
  // numerics holds the CDiscGalerkin class for compute residual here
  
  if ((config->GetRiemann_Solver_FEM() == ROE) && (rank == MASTER_NODE))
    cout << " Using a Roe Riemann solver for the DG method." << endl;
  
  su2double quad_fact_straight = config->GetQuadrature_Factor_Straight();
  if (rank == MASTER_NODE) cout << " Quad factor straight: " << quad_fact_straight << endl;
  
  su2double quad_fact_curved = config->GetQuadrature_Factor_Curved();
  if (rank == MASTER_NODE) cout << " Quad factor curved: " << quad_fact_curved << endl;
  
}

void CFEM_DG_NSSolver::External_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Compute finite element convective residual and store here. ---*/
  
  if (rank == MASTER_NODE) cout << " Computing contour contributions to residual." << endl;
  
}

void CFEM_DG_NSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Apply adiabatic / heat flux wall boundary condition here. ---*/
  
  if (rank == MASTER_NODE) cout << " Heat flux wall BC." << endl;
  
}

void CFEM_DG_NSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Apply isothermal wall boundary condition here. ---*/
  
  if (rank == MASTER_NODE) cout << " Isothermal wall BC." << endl;
  
}
