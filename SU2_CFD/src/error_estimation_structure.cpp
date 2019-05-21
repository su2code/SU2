/*!
 * \file error_estimation_structure.cpp
 * \brief Functions for error estimation.
 * \author B. Munguía
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../include/error_estimation_structure.hpp"

CErrorEstimationDriver::CErrorEstimationDriver(char* confFile,
                                               unsigned short val_nZone,
                                               unsigned short val_nDim, 
                                               bool val_periodic,
                                               SU2_Comm MPICommunicator):config_file_name(confFile), nZone(val_nZone), nDim(val_nDim), fsi(false), fem_solver(false) {

  char zone_file_name[MAX_STRING_SIZE];

  SU2_MPI::SetComm(MPICommunicator);

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  geometry_container            = NULL;
  solver_container              = NULL;
  config_container              = NULL;
  output                        = NULL;
  iteration_container           = NULL;
  direct_iteration              = NULL;
  integration_container         = NULL;
  numerics_container            = NULL;
  surface_movement              = NULL;
  grid_movement                 = NULL;
  FFDBox                        = NULL;
  nInst                         = NULL;

  /*--- Definition of the containers for all possible zones. ---*/

  geometry_container            = new CGeometry***[nZone];
  solver_container              = new CSolver****[nZone];
  config_container              = new CConfig*[nZone];
  iteration_container           = new CIteration**[nZone];
  direct_iteration              = new CIteration*[nZone];
  integration_container         = new CIntegration***[nZone];
  numerics_container            = new CNumerics*****[nZone];
  nInst                         = new unsigned short[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone]                   = NULL;
    geometry_container[iZone]                 = NULL;
    config_container[iZone]                   = NULL;
    integration_container[iZone]              = NULL;
    numerics_container[iZone]                 = NULL;
    nInst[iZone]                              = 1;
  }

  /*--- Initialize the configuration of the driver ---*/

  driver_config = new CConfig(config_file_name, SU2_MET, ZONE_0, nZone, nDim, VERB_NONE);

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    if (driver_config->GetKind_Solver() == MULTIZONE){
      strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
      config_container[iZone] = new CConfig(zone_file_name, SU2_MET, iZone, nZone, nDim, VERB_HIGH);
    }
    else{
      config_container[iZone] = new CConfig(config_file_name, SU2_MET, iZone, nZone, nDim, VERB_HIGH);
    }

    /*--- Set the MPI communicator ---*/

    config_container[iZone]->SetMPICommunicator(MPICommunicator);

    /* --- For the output config, enable restart reading (to export sensor fields) and change grid file to target mesh ---*/

    config_container[iZone]->SetRestart(true);

    config_container[iZone]->SetMGLevels(0);

  }

  /*--- Set the multizone part of the problem. ---*/
  if (driver_config->GetKind_Solver() == MULTIZONE){
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config_container[iZone]->SetMultizone(driver_config, config_container);
    }
  }

  /*--- Preprocessing of the config and mesh files. In this routine, the config file is read
   and it is determined whether a problem is single physics or multiphysics. . ---*/

  Input_Preprocessing(config_container, geometry_container, val_periodic);

  /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
   based data structure is constructed, i.e. node and cell neighbors are
   identified and linked, and face areas are computed ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Geometry Preprocessing ------------------------" << endl;

  /*--- Determine whether or not the FEM solver is used, which decides the
   type of geometry classes that are instantiated. Only adapted for single-zone problems ---*/
  fem_solver = ((config_container[ZONE_0]->GetKind_Solver() == FEM_EULER)          ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_NAVIER_STOKES)  ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_RANS)           ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_LES)            ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_EULER) ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_NS)    ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_RANS));

  if( fem_solver ) {
    switch( config_container[ZONE_0]->GetKind_FEM_Flow() ) {
      case DG: {
        Geometrical_Preprocessing_DGFEM(config_container, geometry_container);
        break;
      }
    }
  }
  else {
    Geometrical_Preprocessing(config_container, geometry_container);
  }

  for (iZone = 0; iZone < nZone; iZone++) {

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      /*--- Computation of positive surface area in the z-plane which is used for
        the calculation of force coefficient (non-dimensionalization). ---*/

      geometry_container[iZone][iInst][MESH_0]->SetPositive_ZArea(config_container[iZone]);

    }

  }

  for (iZone = 0; iZone < nZone; iZone++) {

  	/*--- Instantiate the type of physics iteration to be executed within each zone. For
     example, one can execute the same physics across multiple zones (mixing plane),
     different physics in different zones (fluid-structure interaction), or couple multiple
     systems tightly within a single zone by creating a new iteration class (e.g., RANS). ---*/
    
    if (rank == MASTER_NODE) {
      cout << endl <<"------------------------ Iteration Preprocessing ------------------------" << endl;
    }

    iteration_container[iZone] = new CIteration* [nInst[iZone]];
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      iteration_container[iZone][iInst] = NULL;
    }

    Iteration_Preprocessing();

    RecordingState = NONE;
    ExtIter        = 0;

    if(config_container[iZone]->GetBoolTurbomachinery()){
      direct_iteration[iZone] = new CTurboIteration(config_container[iZone]);
    }
    else{
      direct_iteration[iZone] = new CFluidIteration(config_container[iZone]);
    }

    /*--- Definition of the solver class: solver_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS].
    The solver classes are specific to a particular set of governing equations,
    and they contain the subroutines with instructions for computing each spatial
    term of the PDE, i.e. loops over the edges to compute convective and viscous
    fluxes, loops over the nodes to compute source terms, and routines for
    imposing various boundary condition type for the PDE. ---*/

    if (rank == MASTER_NODE)
      cout << endl <<"------------------------- Solver Preprocessing --------------------------" << endl;

    solver_container[iZone] = new CSolver*** [nInst[iZone]];


    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      solver_container[iZone][iInst] = NULL;
      solver_container[iZone][iInst] = new CSolver** [MESH_0+1];
      solver_container[iZone][iInst][MESH_0] = NULL;
      solver_container[iZone][iInst][MESH_0] = new CSolver* [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        solver_container[iZone][iInst][MESH_0][iSol] = NULL;
      
      Solver_Preprocessing(solver_container[iZone], geometry_container[iZone],
                           config_container[iZone], iInst);

    } // End of loop over iInst

  }

  for (iZone = 0; iZone < nZone; iZone++) {

    if (rank == MASTER_NODE)
      cout << endl <<"----------------- Integration and Numerics Preprocessing ----------------" << endl;

    /*--- Definition of the integration class: integration_container[#ZONES][#INSTANCES][#EQ_SYSTEMS].
     The integration class orchestrates the execution of the spatial integration
     subroutines contained in the solver class (including multigrid) for computing
     the residual at each node, R(U) and then integrates the equations to a
     steady state or time-accurately. ---*/

    integration_container[iZone] = new CIntegration** [nInst[iZone]];
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      integration_container[iZone][iInst] = NULL;

      integration_container[iZone][iInst] = new CIntegration*[MAX_SOLS];
      Integration_Preprocessing(integration_container[iZone], geometry_container[iZone],
                                config_container[iZone], iInst);
    }
    
    if (rank == MASTER_NODE) cout << "Integration Preprocessing." << endl;

    /*--- Definition of the numerical method class:
     numerics_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
     The numerics class contains the implementation of the numerical methods for
     evaluating convective or viscous fluxes between any two nodes in the edge-based
     data structure (centered, upwind, galerkin), as well as any source terms
     (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/

    numerics_container[iZone] = new CNumerics****[nInst[iZone]];
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      numerics_container[iZone][iInst] = NULL;

      numerics_container[iZone][iInst] = new CNumerics***[config_container[iZone]->GetnMGLevels()+1];

      Numerics_Preprocessing(numerics_container[iZone], solver_container[iZone],
                             geometry_container[iZone], config_container[iZone], iInst);
    }

    if (rank == MASTER_NODE) cout << "Numerics Preprocessing." << endl;

  }

  /*--- Definition of the output class (one for all zones). The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/

  output = new COutput(config_container[ZONE_0]);

}

void CErrorEstimationDriver::Input_Preprocessing(CConfig **config_container, CGeometry ****geometry_container, bool val_periodic) {

  bool fem_solver = false; 

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Determine whether or not the FEM solver is used, which decides the
     type of geometry classes that are instantiated. ---*/
    fem_solver = ((config_container[iZone]->GetKind_Solver() == FEM_EULER)         ||
                  (config_container[iZone]->GetKind_Solver() == FEM_NAVIER_STOKES) ||
                  (config_container[iZone]->GetKind_Solver() == FEM_RANS)          ||
                  (config_container[iZone]->GetKind_Solver() == FEM_LES)           ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_EULER) ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_NS)    ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_RANS));

    /*--- Read the number of instances for each zone ---*/

    nInst[iZone] = config_container[iZone]->GetnTimeInstances();

    geometry_container[iZone] = new CGeometry** [nInst[iZone]];

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      config_container[iZone]->SetiInst(iInst);

      /*--- Definition of the geometry class to store the primal grid in the
     partitioning process. ---*/

      CGeometry *geometry_aux = NULL;

      /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

      geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

      /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

      if ( fem_solver ) geometry_aux->SetColorFEMGrid_Parallel(config_container[iZone]);
      else              geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

      /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

      geometry_container[iZone][iInst]            = NULL;
      geometry_container[iZone][iInst]            = new CGeometry *[MESH_0+1];
      geometry_container[iZone][iInst][MESH_0]    = NULL;


      if( fem_solver ) {
        switch( config_container[iZone]->GetKind_FEM_Flow() ) {
          case DG: {
            geometry_container[iZone][iInst][MESH_0] = new CMeshFEM_DG(geometry_aux, config_container[iZone]);
            break;
          }

          default: {
            SU2_MPI::Error("Unknown FEM flow solver.", CURRENT_FUNCTION);
            break;
          }
        }
      }
      else {

        /*--- Build the grid data structures using the ParMETIS coloring. ---*/
        
        geometry_container[iZone][iInst][MESH_0] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);

      }

      /*--- Deallocate the memory of geometry_aux and solver_aux ---*/

      delete geometry_aux;

      /*--- Add the Send/Receive boundaries ---*/
      geometry_container[iZone][iInst][MESH_0]->SetSendReceive(config_container[iZone]);

      /*--- Add the Send/Receive boundaries ---*/
      geometry_container[iZone][iInst][MESH_0]->SetBoundaries(config_container[iZone]);

    }

  }

}

void CErrorEstimationDriver::Geometrical_Preprocessing(CConfig **config_container, CGeometry ****geometry_container) {

  bool fea = false;

  for (iZone = 0; iZone < nZone; iZone++) {

    fea = ((config_container[iZone]->GetKind_Solver() == FEM_ELASTICITY) ||
        (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM));

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      /*--- Compute elements surrounding points, points surrounding points ---*/

      if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetPoint_Connectivity();

      /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/

      // if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
      // geometry_container[iZone][iInst][MESH_0]->SetRCM_Ordering(config_container[iZone]);

      /*--- recompute elements surrounding points, points surrounding points ---*/

      if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetPoint_Connectivity();

      /*--- Compute elements surrounding elements ---*/

      if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetElement_Connectivity();

      /*--- Check the orientation before computing geometrical quantities ---*/

      geometry_container[iZone][iInst][MESH_0]->SetBoundVolume();
      if (config_container[iZone]->GetReorientElements()) {
        if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
        geometry_container[iZone][iInst][MESH_0]->Check_IntElem_Orientation(config_container[iZone]);
        geometry_container[iZone][iInst][MESH_0]->Check_BoundElem_Orientation(config_container[iZone]);
      }

      /*--- Create the edge structure ---*/

      if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetEdges();
      geometry_container[iZone][iInst][MESH_0]->SetVertex(config_container[iZone]);

      /*--- Compute cell center of gravity ---*/

      if ((rank == MASTER_NODE) && (!fea)) cout << "Computing centers of gravity." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetCoord_CG();

      /*--- Create the control volume structures ---*/

      if ((rank == MASTER_NODE) && (!fea)) cout << "Setting the control volume structure." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetControlVolume(config_container[iZone], ALLOCATE);
      geometry_container[iZone][iInst][MESH_0]->SetBoundControlVolume(config_container[iZone], ALLOCATE);

      /*--- Visualize a dual control volume if requested ---*/

      if ((config_container[iZone]->GetVisualize_CV() >= 0) &&
          (config_container[iZone]->GetVisualize_CV() < (long)geometry_container[iZone][iInst][MESH_0]->GetnPointDomain()))
        geometry_container[iZone][iInst][MESH_0]->VisualizeControlVolume(config_container[iZone], UPDATE);

      /*--- Identify closest normal neighbor ---*/

      if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
      geometry_container[iZone][iInst][MESH_0]->FindNormal_Neighbor(config_container[iZone]);

      /*--- Store the global to local mapping. ---*/

      if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetGlobal_to_Local_Point();

      /*--- Compute the surface curvature ---*/

      if ((rank == MASTER_NODE) && (!fea)) cout << "Compute the surface curvature." << endl;
      geometry_container[iZone][iInst][MESH_0]->ComputeSurf_Curvature(config_container[iZone]);

      /*--- Check for periodicity and disable MG if necessary. ---*/

      if (rank == MASTER_NODE) cout << "Checking for periodicity." << endl;
      geometry_container[iZone][iInst][MESH_0]->Check_Periodicity(config_container[iZone]);

    }

  }

  /*--- Create the data structure for MPI point-to-point communications. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        geometry_container[iZone][iInst][MESH_0]->PreprocessP2PComms(geometry_container[iZone][iInst][MESH_0], config_container[iZone]);
    }
  }
  
  /*--- Perform a few preprocessing routines and communications. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        
      /*--- Compute the max length. ---*/
        
      if ((rank == MASTER_NODE) && (!fea)) cout << "Finding max control volume width." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetMaxLength(config_container[iZone]);
        
      /*--- Communicate the number of neighbors. This is needed for
       some centered schemes and for multigrid in parallel. ---*/
        
      if ((rank == MASTER_NODE) && (size > SINGLE_NODE) && (!fea)) cout << "Communicating number of neighbors." << endl;
      geometry_container[iZone][iInst][MESH_0]->InitiateComms(geometry_container[iZone][iInst][MESH_0], config_container[iZone], NEIGHBORS);
      geometry_container[iZone][iInst][MESH_0]->CompleteComms(geometry_container[iZone][iInst][MESH_0], config_container[iZone], NEIGHBORS);
    }
  }

}

void CErrorEstimationDriver::Geometrical_Preprocessing_DGFEM(CConfig **config_container, CGeometry ****geometry_container) {

  //*--- Loop over the number of zones of the fine grid. ---*/

  for(unsigned short iZone = 0; iZone < nZone; iZone++) {

    /*--- Loop over the time instances of this zone. ---*/
    for(unsigned short iInst = 0; iInst < nInst[iZone]; iInst++) {

      /*--- Carry out a dynamic cast to CMeshFEM_DG, such that it is not needed to
       define all virtual functions in the base class CGeometry. ---*/
      CMeshFEM_DG *DGMesh = dynamic_cast<CMeshFEM_DG *>(geometry_container[iZone][iInst][MESH_0]);

      /*--- Determine the standard elements for the volume elements. ---*/
      if (rank == MASTER_NODE) cout << "Creating standard volume elements." << endl;
      DGMesh->CreateStandardVolumeElements(config_container[iZone]);

      /*--- Create the face information needed to compute the contour integral
       for the elements in the Discontinuous Galerkin formulation. ---*/
      if (rank == MASTER_NODE) cout << "Creating face information." << endl;
      DGMesh->CreateFaces(config_container[iZone]);

      /*--- Compute the metric terms of the volume elements. ---*/
      if (rank == MASTER_NODE) cout << "Computing metric terms volume elements." << endl;
      DGMesh->MetricTermsVolumeElements(config_container[iZone]);

      /*--- Compute the metric terms of the surface elements. ---*/
      if (rank == MASTER_NODE) cout << "Computing metric terms surface elements." << endl;
      DGMesh->MetricTermsSurfaceElements(config_container[iZone]);

      /*--- Compute a length scale of the volume elements. ---*/
      if (rank == MASTER_NODE) cout << "Computing length scale volume elements." << endl;
      DGMesh->LengthScaleVolumeElements();

      /*--- Compute the coordinates of the integration points. ---*/
      if (rank == MASTER_NODE) cout << "Computing coordinates of the integration points." << endl;
      DGMesh->CoordinatesIntegrationPoints();

      /*--- Compute the coordinates of the location of the solution DOFs. This is different
            from the grid points when a different polynomial degree is used to represent the
            geometry and solution. ---*/
      if (rank == MASTER_NODE) cout << "Computing coordinates of the solution DOFs." << endl;
      DGMesh->CoordinatesSolDOFs();

      /*--- Initialize the static mesh movement, if necessary. ---*/
      const unsigned short Kind_Grid_Movement = config_container[iZone]->GetKind_GridMovement(iZone);
      const bool initStaticMovement = (config_container[iZone]->GetGrid_Movement() &&
                                      (Kind_Grid_Movement == MOVING_WALL    ||
                                       Kind_Grid_Movement == ROTATING_FRAME ||
                                       Kind_Grid_Movement == STEADY_TRANSLATION));

      if(initStaticMovement){
        if (rank == MASTER_NODE) cout << "Initialize Static Mesh Movement" << endl;
        DGMesh->InitStaticMeshMovement(config_container[iZone], Kind_Grid_Movement, iZone);
      }

      /*--- Perform the preprocessing tasks when wall functions are used. ---*/
      if (rank == MASTER_NODE) cout << "Preprocessing for the wall functions. " << endl;
      DGMesh->WallFunctionPreprocessing(config_container[iZone]);

      /*--- Store the global to local mapping. ---*/
      if (rank == MASTER_NODE) cout << "Storing a mapping from global to local DOF index." << endl;
      geometry_container[iZone][iInst][MESH_0]->SetGlobal_to_Local_Point();
    }

    /*--- Loop to create the coarser grid levels. ---*/

    for(unsigned short iMGlevel=1; iMGlevel<=config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

      SU2_MPI::Error("Geometrical_Preprocessing_DGFEM: Coarse grid levels not implemented yet.",
                     CURRENT_FUNCTION);
    }
  }
}

void CErrorEstimationDriver::Iteration_Preprocessing() {

  for (iInst = 0; iInst < nInst[iZone]; iInst++)  {

    /*--- Initial print to console for this zone. ---*/

    if (rank == MASTER_NODE) cout << "Zone " << iZone+1;
    if ((rank == MASTER_NODE) && (nInst[iZone] > 1)) cout << ", instance: " << iInst+1;

    /*--- Loop over all zones and instantiate the physics iteration. ---*/

    switch (config_container[iZone]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:

      if(config_container[iZone]->GetBoolTurbomachinery()){
        if (rank == MASTER_NODE)
          cout << ": Euler/Navier-Stokes/RANS turbomachinery fluid iteration." << endl;
        iteration_container[iZone][iInst] = new CTurboIteration(config_container[iZone]);

      }
      else{
        if (rank == MASTER_NODE)
          cout << ": Euler/Navier-Stokes/RANS fluid iteration." << endl;
        // iteration_container[iZone][iInst] = new CFluidIteration(config_container[iZone]);
        iteration_container[iZone][iInst] = new CDiscAdjFluidIteration(config_container[iZone]);
      }
      break;

    case FEM_EULER: case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES:
      if (rank == MASTER_NODE)
        cout << ": finite element Euler/Navier-Stokes/RANS/LES flow iteration." << endl;
      iteration_container[iZone][iInst] = new CFEMFluidIteration(config_container[iZone]);
      break;

    case HEAT_EQUATION_FVM:
      if (rank == MASTER_NODE)
        cout << ": heat iteration (finite volume method)." << endl;
      iteration_container[iZone][iInst] = new CHeatIteration(config_container[iZone]);
      break;

    case FEM_ELASTICITY:
      if (rank == MASTER_NODE)
        cout << ": FEM iteration." << endl;
      iteration_container[iZone][iInst] = new CFEAIteration(config_container[iZone]);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << ": adjoint Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration_container[iZone][iInst] = new CAdjFluidIteration(config_container[iZone]);
      break;

    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << ": discrete adjoint Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration_container[iZone][iInst] = new CDiscAdjFluidIteration(config_container[iZone]);
      break;

    case DISC_ADJ_FEM_EULER : case DISC_ADJ_FEM_NS : case DISC_ADJ_FEM_RANS :
      if (rank == MASTER_NODE)
        cout << ": discrete adjoint finite element Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration_container[iZone][iInst] = new CDiscAdjFluidIteration(config_container[iZone]);
      break;

    case DISC_ADJ_FEM:
      if (rank == MASTER_NODE)
        cout << ": discrete adjoint FEM structural iteration." << endl;
      iteration_container[iZone][iInst] = new CDiscAdjFEAIteration(config_container[iZone]);
      break;

    case DISC_ADJ_HEAT:
      if (rank == MASTER_NODE)
        cout << ": discrete adjoint heat iteration." << endl;
      iteration_container[iZone][iInst] = new CDiscAdjHeatIteration(config_container[iZone]);
      break;
    }

  }

}

void CErrorEstimationDriver::Solver_Preprocessing(CSolver ****solver_container, CGeometry ***geometry,
                                                   CConfig *config, unsigned short val_iInst) {
  
  bool euler, ns, turbulent,
  fem_euler, fem_ns, fem_turbulent, fem_transition,
  adj_euler, adj_ns, adj_turb,
  heat_fvm,
  fem, disc_adj_fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb, disc_adj_heat,
  fem_dg_flow, fem_dg_shock_persson,
  e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;
  
  /*--- Count the number of DOFs per solution point. ---*/
  
  DOFsPerPoint = 0;
  
  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns              = false;  turbulent     = false;
  fem_euler        = false;  fem_ns          = false;  fem_turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb      = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  heat_fvm         = false;  disc_adj_heat    = false;
  transition       = false;  fem_transition   = false;
  template_solver  = false;
  fem_dg_flow      = false;  fem_dg_shock_persson = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;
  
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; fem_turbulent = true; if(config->GetKind_Trans_Model() == LM) fem_transition = true; break;
    case FEM_LES : fem_ns = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; fem_turbulent = true; disc_adj = true; if(config->GetKind_Trans_Model() == LM) fem_transition = true; break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }
  
  /*--- Determine the kind of FEM solver used for the flow. ---*/

  switch( config->GetKind_FEM_Flow() ) {
    case DG: fem_dg_flow = true; break;
  }

  /*--- Determine the kind of shock capturing method for FEM DG solver. ---*/

  switch( config->GetKind_FEM_DG_Shock() ) {
    case PERSSON: fem_dg_shock_persson = true; break;
  }

  /*--- Assign turbulence model booleans ---*/

  if (turbulent || fem_turbulent){
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;           break;
      case SA_E:   e_spalart_allmaras = true;   break;
      case SA_COMP: comp_spalart_allmaras = true; break;
      case SA_E_COMP: e_comp_spalart_allmaras = true; break;
      default: SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION); break;
    }
  }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][INSTANCE][MESH_0][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/


    /*--- Allocate solution for a template problem ---*/

  if (template_solver) {
    solver_container[val_iInst][MESH_0][TEMPLATE_SOL] = new CTemplateSolver(geometry[val_iInst][MESH_0], config);
    DOFsPerPoint += solver_container[val_iInst][MESH_0][TEMPLATE_SOL]->GetnVar();
  }

    /*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/
    /*--- Note that we need both direct and adjoint for the error estimate ---*/

  if (euler || disc_adj) {
    if (compressible) {
      solver_container[val_iInst][MESH_0][FLOW_SOL] = new CEulerSolver(geometry[val_iInst][MESH_0], config, MESH_0);
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
    if (incompressible) {
      solver_container[val_iInst][MESH_0][FLOW_SOL] = new CIncEulerSolver(geometry[val_iInst][MESH_0], config, MESH_0);
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
    DOFsPerPoint += solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  }
  if (ns) {
    if (compressible) {
      solver_container[val_iInst][MESH_0][FLOW_SOL] = new CNSSolver(geometry[val_iInst][MESH_0], config, MESH_0);
    }
    if (incompressible) {
      solver_container[val_iInst][MESH_0][FLOW_SOL] = new CIncNSSolver(geometry[val_iInst][MESH_0], config, MESH_0);
    }
    DOFsPerPoint += solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  }
  if (turbulent) {
    if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras || neg_spalart_allmaras) {
      solver_container[val_iInst][MESH_0][TURB_SOL] = new CTurbSASolver(geometry[val_iInst][MESH_0], config, MESH_0, solver_container[val_iInst][MESH_0][FLOW_SOL]->GetFluidModel() );
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      solver_container[val_iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0);
    }
    else if (menter_sst) {
      solver_container[val_iInst][MESH_0][TURB_SOL] = new CTurbSSTSolver(geometry[val_iInst][MESH_0], config, MESH_0);
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      solver_container[val_iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0);
      solver_container[val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iInst][MESH_0], solver_container[val_iInst][MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
    DOFsPerPoint += solver_container[val_iInst][MESH_0][TURB_SOL]->GetnVar();
    if (transition) {
      solver_container[val_iInst][MESH_0][TRANS_SOL] = new CTransLMSolver(geometry[val_iInst][MESH_0], config, MESH_0);
      DOFsPerPoint += solver_container[val_iInst][MESH_0][TRANS_SOL]->GetnVar();
    }
  }

    /*--- Allocate solution for adjoint problem ---*/

  if (disc_adj || euler || ns || turbulent) {
    solver_container[val_iInst][MESH_0][ADJFLOW_SOL] = new CDiscAdjSolver(geometry[val_iInst][MESH_0], config, solver_container[val_iInst][MESH_0][FLOW_SOL], RUNTIME_FLOW_SYS, MESH_0);
    DOFsPerPoint += solver_container[val_iInst][MESH_0][ADJFLOW_SOL]->GetnVar();
    if (disc_adj_turb || turbulent) {
      solver_container[val_iInst][MESH_0][ADJTURB_SOL] = new CDiscAdjSolver(geometry[val_iInst][MESH_0], config, solver_container[val_iInst][MESH_0][TURB_SOL], RUNTIME_TURB_SYS, MESH_0);
      DOFsPerPoint += solver_container[val_iInst][MESH_0][ADJTURB_SOL]->GetnVar();
    }
  }  

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  bool update_geo = true;
  if (config->GetFSI_Simulation()) update_geo = false;

  Solver_Restart(solver_container, geometry, config, update_geo, val_iInst);

}

void CErrorEstimationDriver::Integration_Preprocessing(CIntegration ***integration_container,
    CGeometry ***geometry, CConfig *config, unsigned short val_iInst) {

  bool euler, adj_euler, ns, adj_ns, turbulent, adj_turb, fem,
      fem_euler, fem_ns, fem_turbulent,
      heat_fvm, template_solver, transition, disc_adj, disc_adj_fem, disc_adj_heat;

  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false;
  ns               = false; adj_ns           = false;
  turbulent        = false; adj_turb         = false;
  disc_adj         = false;
  fem_euler        = false;
  fem_ns           = false;
  fem_turbulent    = false;
  heat_fvm         = false; disc_adj_heat    = false;
  fem 			       = false; disc_adj_fem     = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true;  heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; fem_turbulent = true; break;
    case FEM_LES :  fem_ns = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER : euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; fem_turbulent = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS : ns = true; turbulent = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }

  /*--- Allocate solution for a template problem ---*/
  if (template_solver) integration_container[val_iInst][TEMPLATE_SOL] = new CSingleGridIntegration(config);

  /*--- Allocate solution for direct problem ---*/
  if (euler) integration_container[val_iInst][FLOW_SOL] = new CMultiGridIntegration(config);
  if (ns) integration_container[val_iInst][FLOW_SOL] = new CMultiGridIntegration(config);
  if (turbulent) integration_container[val_iInst][TURB_SOL] = new CSingleGridIntegration(config);
  if (transition) integration_container[val_iInst][TRANS_SOL] = new CSingleGridIntegration(config);
  if (heat_fvm) integration_container[val_iInst][HEAT_SOL] = new CSingleGridIntegration(config);
  if (fem) integration_container[val_iInst][FEA_SOL] = new CStructuralIntegration(config);

  /*--- Allocate integration container for finite element flow solver. ---*/

  if (fem_euler) integration_container[val_iInst][FLOW_SOL] = new CFEM_DG_Integration(config);
  if (fem_ns)    integration_container[val_iInst][FLOW_SOL] = new CFEM_DG_Integration(config);
  //if (fem_turbulent) integration_container[val_iInst][FEM_TURB_SOL] = new CSingleGridIntegration(config);

  if (fem_turbulent)
    SU2_MPI::Error("No turbulent FEM solver yet", CURRENT_FUNCTION);

  /*--- Allocate solution for adjoint problem ---*/
  if (adj_euler) integration_container[val_iInst][ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_ns) integration_container[val_iInst][ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_turb) integration_container[val_iInst][ADJTURB_SOL] = new CSingleGridIntegration(config);

  if (disc_adj) integration_container[val_iInst][ADJFLOW_SOL] = new CIntegration(config);
  if (disc_adj_fem) integration_container[val_iInst][ADJFEA_SOL] = new CIntegration(config);
  if (disc_adj_heat) integration_container[val_iInst][ADJHEAT_SOL] = new CIntegration(config);

}

void CErrorEstimationDriver::Numerics_Preprocessing(CNumerics *****numerics_container,
                                     CSolver ****solver_container, CGeometry ***geometry,
                                     CConfig *config, unsigned short val_iInst) {

  unsigned short iMGlevel, iSol, nDim,
  
  nVar_Template         = 0,
  nVar_Flow             = 0,
  nVar_Trans            = 0,
  nVar_Turb             = 0,
  nVar_Adj_Flow         = 0,
  nVar_Adj_Turb         = 0,
  nVar_FEM              = 0,
  nVar_Heat             = 0;
  
  su2double *constants = NULL;
  
  bool
  euler, adj_euler,
  ns, adj_ns,
  turbulent, adj_turb,
  fem_euler, fem_ns, fem_turbulent,
  spalart_allmaras, neg_spalart_allmaras, menter_sst,
  fem,
  heat_fvm,
  transition,
  template_solver;
  bool e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );
  bool roe_low_dissipation = config->GetKind_RoeLowDiss() != NO_ROELOWDISS;
  
  /*--- Initialize some useful booleans ---*/
  euler            = false; ns     = false; turbulent     = false;
  fem_euler        = false; fem_ns = false; fem_turbulent = false;
  adj_euler        = false;   adj_ns           = false;   adj_turb         = false;
  heat_fvm         = false;
  fem              = false;
  spalart_allmaras = false; neg_spalart_allmaras = false;	menter_sst       = false;
  transition       = false;
  template_solver  = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case DISC_ADJ_EULER: euler = true; break;
    case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : case DISC_ADJ_RANS:  ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : case DISC_ADJ_FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: case DISC_ADJ_FEM_NS : fem_ns = true; break;
    case FEM_RANS : case DISC_ADJ_FEM_RANS : fem_ns = true; fem_turbulent = true; break;
    case FEM_LES :  fem_ns = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: case DISC_ADJ_FEM: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
  }
  
  /*--- Assign turbulence model booleans ---*/

  if (turbulent || fem_turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SA_E:   e_spalart_allmaras = true; break;
      case SA_COMP:   comp_spalart_allmaras = true; break;
      case SA_E_COMP:   e_comp_spalart_allmaras = true; break;
      case SST:    menter_sst = true; constants = solver_container[val_iInst][MESH_0][TURB_SOL]->GetConstants(); break;
      default: SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION); break;
    }
  
  /*--- Number of variables for the template ---*/
  
  if (template_solver) nVar_Flow = solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  
  /*--- Number of variables for direct problem ---*/

  if (euler)        nVar_Flow = solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  if (ns)           nVar_Flow = solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  if (turbulent)    nVar_Turb = solver_container[val_iInst][MESH_0][TURB_SOL]->GetnVar();
  if (transition)   nVar_Trans = solver_container[val_iInst][MESH_0][TRANS_SOL]->GetnVar();

  if (fem_euler)        nVar_Flow = solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  if (fem_ns)           nVar_Flow = solver_container[val_iInst][MESH_0][FLOW_SOL]->GetnVar();
  //if (fem_turbulent)    nVar_Turb = solver_container[val_iInst][MESH_0][FEM_TURB_SOL]->GetnVar();
  
  if (fem)          nVar_FEM = solver_container[val_iInst][MESH_0][FEA_SOL]->GetnVar();
  if (heat_fvm)     nVar_Heat = solver_container[val_iInst][MESH_0][HEAT_SOL]->GetnVar();

  /*--- Number of variables for adjoint problem ---*/
  
  if (adj_euler)        nVar_Adj_Flow = solver_container[val_iInst][MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_ns)           nVar_Adj_Flow = solver_container[val_iInst][MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_turb)         nVar_Adj_Turb = solver_container[val_iInst][MESH_0][ADJTURB_SOL]->GetnVar();
  
  /*--- Number of dimensions ---*/
  
  nDim = geometry[val_iInst][MESH_0]->GetnDim();
  
  /*--- Definition of the Class for the numerical method: numerics_container[INSTANCE_LEVEL][MESH_LEVEL][EQUATION][EQ_TERM] ---*/
  if (fem){
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[val_iInst][iMGlevel] = new CNumerics** [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        numerics_container[val_iInst][iMGlevel][iSol] = new CNumerics* [MAX_TERMS_FEA];
    }
  }
  else{
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[val_iInst][iMGlevel] = new CNumerics** [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        numerics_container[val_iInst][iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
    }
  }
  
  /*--- Solver definition for the template problem ---*/
  if (template_solver) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[val_iInst][iMGlevel][TEMPLATE_SOL][CONV_TERM] = new CConvective_Template(nDim, nVar_Template, config);
        break;
      default : SU2_MPI::Error("Convective scheme not implemented (template_solver).", CURRENT_FUNCTION); break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics_container[val_iInst][iMGlevel][TEMPLATE_SOL][VISC_TERM] = new CViscous_Template(nDim, nVar_Template, config);
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics_container[val_iInst][iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM] = new CSource_Template(nDim, nVar_Template, config);
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[val_iInst][iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM] = new CConvective_Template(nDim, nVar_Template, config);
    }
    
  }
  
  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        SU2_MPI::Error("No convective scheme.", CURRENT_FUNCTION);
        break;
        
      case SPACE_CENTERED :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[val_iInst][MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics_container[val_iInst][MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim, nVar_Flow, config); break;
            case JST_KE : numerics_container[val_iInst][MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_KE_Flow(nDim, nVar_Flow, config); break;
            default : SU2_MPI::Error("Centered scheme not implemented.", CURRENT_FUNCTION); break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config, false);
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use preconditioning method ---*/
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[val_iInst][MESH_0][FLOW_SOL][CONV_TERM] = new CCentLaxInc_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics_container[val_iInst][MESH_0][FLOW_SOL][CONV_TERM] = new CCentJSTInc_Flow(nDim, nVar_Flow, config); break;
            default : SU2_MPI::Error("Centered scheme not implemented.\n Currently, only JST and LAX-FRIEDRICH are available for incompressible flows.", CURRENT_FUNCTION); break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLaxInc_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);
          
        }
        break;
      case SPACE_UPWIND :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              if (ideal_gas) {
                
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                  numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config, false);
                }
              } else {
                
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                  numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                }
              }
              break;
              
            case AUSM:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              }
              break;

	    case AUSMPLUSUP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSMPLUSUP_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSMPLUSUP_Flow(nDim, nVar_Flow, config);
              }
              break;

	    case AUSMPLUSUP2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSMPLUSUP2_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSMPLUSUP2_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case TURKEL:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
              }
              break;
                  
            case L2ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwL2Roe_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwL2Roe_Flow(nDim, nVar_Flow, config);
              }
              break;
            case LMROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwLMRoe_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwLMRoe_Flow(nDim, nVar_Flow, config);
              }
              break;

            case SLAU:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwSLAU_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwSLAU_Flow(nDim, nVar_Flow, config, false);
              }
              break;
              
            case SLAU2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwSLAU2_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwSLAU2_Flow(nDim, nVar_Flow, config, false);
              }
              break;
              
            case HLLC:
              if (ideal_gas) {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                  numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              else {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                  numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              break;
              
            case MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            default : SU2_MPI::Error("Upwind scheme not implemented.", CURRENT_FUNCTION); break;
          }
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case FDS:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);
                numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);
              }
              break;
            default : SU2_MPI::Error("Upwind scheme not implemented.\n Currently, only FDS is available for incompressible flows.", CURRENT_FUNCTION); break;
          }
        }
        break;
        
      default :
        SU2_MPI::Error("Convective scheme not implemented (Euler and Navier-Stokes).", CURRENT_FUNCTION);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible) {
      if (ideal_gas) {
        
        /*--- Compressible flow Ideal gas ---*/
        numerics_container[val_iInst][MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, true, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, false, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, false, config);
        
      } else {
        
        /*--- Compressible flow Realgas ---*/
        numerics_container[val_iInst][MESH_0][FLOW_SOL][VISC_TERM] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, true, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][VISC_TERM] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, false, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, false, config);
        
      }
    }
    if (incompressible) {
      /*--- Incompressible flow, use preconditioning method ---*/
      numerics_container[val_iInst][MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradInc_Flow(nDim, nVar_Flow, true, config);
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics_container[val_iInst][iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradInc_Flow(nDim, nVar_Flow, false, config);
      
      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics_container[val_iInst][iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradInc_Flow(nDim, nVar_Flow, false, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      
      if (config->GetBody_Force() == YES)
        if (incompressible) numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceIncBodyForce(nDim, nVar_Flow, config);
        else numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceBodyForce(nDim, nVar_Flow, config);
      else if (incompressible && (config->GetKind_DensityModel() == BOUSSINESQ))
        numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceBoussinesq(nDim, nVar_Flow, config);
      else if (config->GetRotating_Frame() == YES)
        numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_Flow(nDim, nVar_Flow, config);
      else if (config->GetAxisymmetric() == YES)
        if (incompressible) numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceIncAxisymmetric_Flow(nDim, nVar_Flow, config);
      else numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_Flow(nDim, nVar_Flow, config);
      else if (config->GetGravityForce() == YES)
        numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceGravity(nDim, nVar_Flow, config);
      else if (config->GetWind_Gust() == YES)
        numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceWindGust(nDim, nVar_Flow, config);
      else
        numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
      
      numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
    }
    
  }

  /*--- Riemann solver definition for the Euler, Navier-Stokes problems for the FEM discretization. ---*/
  if ((fem_euler) || (fem_ns)) {

    switch (config->GetRiemann_Solver_FEM()) {
      case NO_UPWIND : cout << "Riemann solver disabled." << endl; break;
      case ROE:
      case LAX_FRIEDRICH:
        /* Hard coded optimized implementation is used in the DG solver. No need to allocate the
           corresponding entry in numerics. */
        break;

      case AUSM:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
        }
        break;

      case TURKEL:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
        }
        break;

      case HLLC:
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
            numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
            numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
          }
        break;

      case MSW:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
        }
        break;

      case CUSP:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
          numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
        }
        break;

      default :
        SU2_MPI::Error("Riemann solver not implemented.", CURRENT_FUNCTION);
    }

  }

  /*--- Solver definition for the turbulent model problem ---*/
  
  if (turbulent) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          if (spalart_allmaras || neg_spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras ) {
            numerics_container[val_iInst][iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
          }
          else if (menter_sst) numerics_container[val_iInst][iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        }
        break;
      default :
        SU2_MPI::Error("Convective scheme not implemented (turbulent).", CURRENT_FUNCTION);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras){
        numerics_container[val_iInst][iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, true, config);
      }
      else if (neg_spalart_allmaras) numerics_container[val_iInst][iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSA_Neg(nDim, nVar_Turb, true, config);
      else if (menter_sst) numerics_container[val_iInst][iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, true, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
      else if (e_spalart_allmaras) numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA_E(nDim, nVar_Turb, config);
      else if (comp_spalart_allmaras) numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA_COMP(nDim, nVar_Turb, config);
      else if (e_comp_spalart_allmaras) numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA_E_COMP(nDim, nVar_Turb, config);
      else if (neg_spalart_allmaras) numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA_Neg(nDim, nVar_Turb, config);
      else if (menter_sst) numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSST(nDim, nVar_Turb, constants, config);
      numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Turb, config);
    }
    
    /*--- Definition of the boundary condition method ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
        numerics_container[val_iInst][iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[val_iInst][iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, false, config);
      }
      else if (neg_spalart_allmaras) {
        numerics_container[val_iInst][iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[val_iInst][iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA_Neg(nDim, nVar_Turb, false, config);
      }
      else if (menter_sst) {
        numerics_container[val_iInst][iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        numerics_container[val_iInst][iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, false, config);
      }
    }
  }
  
  /*--- Solver definition for the transition model problem ---*/
  if (transition) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics_container[val_iInst][iMGlevel][TRANS_SOL][CONV_TERM] = new CUpwSca_TransLM(nDim, nVar_Trans, config);
        }
        break;
      default :
        SU2_MPI::Error("Convective scheme not implemented (transition).", CURRENT_FUNCTION);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[val_iInst][iMGlevel][TRANS_SOL][VISC_TERM] = new CAvgGradCorrected_TransLM(nDim, nVar_Trans, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[val_iInst][iMGlevel][TRANS_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TransLM(nDim, nVar_Trans, config);
      numerics_container[val_iInst][iMGlevel][TRANS_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Trans, config);
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[val_iInst][iMGlevel][TRANS_SOL][CONV_BOUND_TERM] = new CUpwLin_TransLM(nDim, nVar_Trans, config);
    }
  }
  
  /*--- Solver definition of the finite volume heat solver  ---*/
  if (heat_fvm) {

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      numerics_container[val_iInst][iMGlevel][HEAT_SOL][VISC_TERM] = new CAvgGradCorrected_Heat(nDim, nVar_Heat, config);
      numerics_container[val_iInst][iMGlevel][HEAT_SOL][VISC_BOUND_TERM] = new CAvgGrad_Heat(nDim, nVar_Heat, config);

      switch (config->GetKind_ConvNumScheme_Heat()) {

        case SPACE_UPWIND :
          numerics_container[val_iInst][iMGlevel][HEAT_SOL][CONV_TERM] = new CUpwSca_Heat(nDim, nVar_Heat, config);
          numerics_container[val_iInst][iMGlevel][HEAT_SOL][CONV_BOUND_TERM] = new CUpwSca_Heat(nDim, nVar_Heat, config);
          break;

        case SPACE_CENTERED :
          numerics_container[val_iInst][iMGlevel][HEAT_SOL][CONV_TERM] = new CCentSca_Heat(nDim, nVar_Heat, config);
          numerics_container[val_iInst][iMGlevel][HEAT_SOL][CONV_BOUND_TERM] = new CUpwSca_Heat(nDim, nVar_Heat, config);
        break;

        default :
          cout << "Convective scheme not implemented (heat)." << endl; exit(EXIT_FAILURE);
        break;
      }
    }
  }
  
  /*--- Solver definition for the flow adjoint problem ---*/
  
  if (adj_euler || adj_ns) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    
    switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case NO_CONVECTIVE :
        SU2_MPI::Error("No convective scheme.", CURRENT_FUNCTION);
        break;
        
      case SPACE_CENTERED :
        
        if (compressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[val_iInst][MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics_container[val_iInst][MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default : SU2_MPI::Error("Centered scheme not implemented.", CURRENT_FUNCTION); break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        
        if (incompressible) {

          SU2_MPI::Error("Schemes not implemented for incompressible continuous adjoint.", CURRENT_FUNCTION);

        }
        
        break;
        
      case SPACE_UPWIND :
        
        if (compressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : SU2_MPI::Error("Upwind scheme not implemented.", CURRENT_FUNCTION); break;
          }
        }
        
        if (incompressible) {
          
          SU2_MPI::Error("Schemes not implemented for incompressible continuous adjoint.", CURRENT_FUNCTION);

        }
        
        break;
        
      default :
        SU2_MPI::Error("Convective scheme not implemented (adj_euler and adj_ns).", CURRENT_FUNCTION);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    
    if (compressible) {
      
      /*--- Compressible flow ---*/
      
      numerics_container[val_iInst][MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
      numerics_container[val_iInst][MESH_0][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
      
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
        numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
      }
      
    }
    
    if (incompressible) {
      
      SU2_MPI::Error("Schemes not implemented for incompressible continuous adjoint.", CURRENT_FUNCTION);

    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      
      /*--- Note that RANS is incompatible with Axisymmetric or Rotational (Fix it!) ---*/
      
      if (compressible) {
        
        if (adj_ns) {
          
          numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceViscous_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          if (config->GetRotating_Frame() == YES)
            numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        
        else {
          
          if (config->GetRotating_Frame() == YES)
            numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else if (config->GetAxisymmetric() == YES)
            numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
          
          numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
          
        }
        
      }
      
      if (incompressible) {
        
        SU2_MPI::Error("Schemes not implemented for incompressible continuous adjoint.", CURRENT_FUNCTION);

      }
      
    }
    
  }
  
  /*--- Solver definition for the turbulent adjoint problem ---*/
  if (adj_turb) {
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTurb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          if (spalart_allmaras) {
            numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][CONV_TERM] = new CUpwSca_AdjTurb(nDim, nVar_Adj_Turb, config);
          }
          else if (neg_spalart_allmaras) {SU2_MPI::Error("Adjoint Neg SA turbulence model not implemented.", CURRENT_FUNCTION);}
          else if (menter_sst) {SU2_MPI::Error("Adjoint SST turbulence model not implemented.", CURRENT_FUNCTION);}
          else if (e_spalart_allmaras) {SU2_MPI::Error("Adjoint Edward's SA turbulence model not implemented.", CURRENT_FUNCTION);}
          else if (comp_spalart_allmaras) {SU2_MPI::Error("Adjoint CC SA turbulence model not implemented.", CURRENT_FUNCTION);}
          else if (e_comp_spalart_allmaras) {SU2_MPI::Error("Adjoint CC Edward's SA turbulence model not implemented.", CURRENT_FUNCTION);}
        break;
      default :
        SU2_MPI::Error("Convective scheme not implemented (adj_turb).", CURRENT_FUNCTION);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][VISC_TERM] = new CAvgGradCorrected_AdjTurb(nDim, nVar_Adj_Turb, config);
      }

      else if (neg_spalart_allmaras) {SU2_MPI::Error("Adjoint Neg SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (menter_sst) {SU2_MPI::Error("Adjoint SST turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (e_spalart_allmaras) {SU2_MPI::Error("Adjoint Edward's SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (comp_spalart_allmaras) {SU2_MPI::Error("Adjoint CC SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (e_comp_spalart_allmaras) {SU2_MPI::Error("Adjoint CC Edward's SA turbulence model not implemented.", CURRENT_FUNCTION);}
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_AdjTurb(nDim, nVar_Adj_Turb, config);
        numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjTurb(nDim, nVar_Adj_Turb, config);
      }
      else if (neg_spalart_allmaras) {SU2_MPI::Error("Adjoint Neg SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (menter_sst) {SU2_MPI::Error("Adjoint SST turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (e_spalart_allmaras) {SU2_MPI::Error("Adjoint Edward's SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (comp_spalart_allmaras) {SU2_MPI::Error("Adjoint CC SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (e_comp_spalart_allmaras) {SU2_MPI::Error("Adjoint CC Edward's SA turbulence model not implemented.", CURRENT_FUNCTION);}
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][CONV_BOUND_TERM] = new CUpwLin_AdjTurb(nDim, nVar_Adj_Turb, config);
      else if (neg_spalart_allmaras) {SU2_MPI::Error("Adjoint Neg SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (menter_sst) {SU2_MPI::Error("Adjoint SST turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (e_spalart_allmaras) {SU2_MPI::Error("Adjoint Edward's SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (comp_spalart_allmaras) {SU2_MPI::Error("Adjoint CC SA turbulence model not implemented.", CURRENT_FUNCTION);}
      else if (e_comp_spalart_allmaras) {SU2_MPI::Error("Adjoint CC Edward's SA turbulence model not implemented.", CURRENT_FUNCTION);}
    }
    
  }

  /*--- Solver definition for the FEM problem ---*/
  if (fem) {

  /*--- Initialize the container for FEA_TERM. This will be the only one for most of the cases ---*/
  switch (config->GetGeometricConditions()) {
      case SMALL_DEFORMATIONS :
        switch (config->GetMaterialModel()) {
          case LINEAR_ELASTIC: numerics_container[val_iInst][MESH_0][FEA_SOL][FEA_TERM] = new CFEALinearElasticity(nDim, nVar_FEM, config); break;
          case NEO_HOOKEAN : SU2_MPI::Error("Material model does not correspond to geometric conditions.", CURRENT_FUNCTION); break;
          default: SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION); break;
        }
        break;
      case LARGE_DEFORMATIONS :
        switch (config->GetMaterialModel()) {
          case LINEAR_ELASTIC: SU2_MPI::Error("Material model does not correspond to geometric conditions.", CURRENT_FUNCTION); break;
          case NEO_HOOKEAN :
            switch (config->GetMaterialCompressibility()) {
              case COMPRESSIBLE_MAT : numerics_container[val_iInst][MESH_0][FEA_SOL][FEA_TERM] = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config); break;
              case INCOMPRESSIBLE_MAT : numerics_container[val_iInst][MESH_0][FEA_SOL][FEA_TERM] = new CFEM_NeoHookean_Incomp(nDim, nVar_FEM, config); break;
              default: SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION); break;
            }
            break;
          case KNOWLES:
            switch (config->GetMaterialCompressibility()) {
              case NEARLY_INCOMPRESSIBLE_MAT : numerics_container[val_iInst][MESH_0][FEA_SOL][FEA_TERM] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config); break;
              case INCOMPRESSIBLE_MAT : numerics_container[val_iInst][MESH_0][FEA_SOL][FEA_TERM] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config); break;
              default:  SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION); break;
            }
            break;
          case IDEAL_DE:
            switch (config->GetMaterialCompressibility()) {
              case NEARLY_INCOMPRESSIBLE_MAT : numerics_container[val_iInst][MESH_0][FEA_SOL][FEA_TERM] = new CFEM_IdealDE(nDim, nVar_FEM, config); break;
              default:  SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION); break;
            }
            break;
          default:  SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION); break;
        }
        break;
      default:  SU2_MPI::Error("Solver not implemented.", CURRENT_FUNCTION);  break;
    }

  /*--- The following definitions only make sense if we have a non-linear solution ---*/
  if (config->GetGeometricConditions() == LARGE_DEFORMATIONS){

      /*--- This allocates a container for electromechanical effects ---*/

      bool de_effects = config->GetDE_Effects();
      if (de_effects) numerics_container[val_iInst][MESH_0][FEA_SOL][DE_TERM] = new CFEM_DielectricElastomer(nDim, nVar_FEM, config);

      string filename;
      ifstream properties_file;

      filename = config->GetFEA_FileName();
      if (nZone > 1)
        filename = config->GetMultizone_FileName(filename, iZone);

      properties_file.open(filename.data(), ios::in);

      /*--- In case there is a properties file, containers are allocated for a number of material models ---*/

      if (!(properties_file.fail())) {

          numerics_container[val_iInst][MESH_0][FEA_SOL][MAT_NHCOMP]  = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config);
          numerics_container[val_iInst][MESH_0][FEA_SOL][MAT_NHINC]   = new CFEM_NeoHookean_Incomp(nDim, nVar_FEM, config);
          numerics_container[val_iInst][MESH_0][FEA_SOL][MAT_IDEALDE] = new CFEM_IdealDE(nDim, nVar_FEM, config);
          numerics_container[val_iInst][MESH_0][FEA_SOL][MAT_KNOWLES] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config);

          properties_file.close();
      }
  }

  }

}

void CErrorEstimationDriver::Solver_Restart(CSolver ****solver_container, CGeometry ***geometry,
                             CConfig *config, bool update_geo, unsigned short val_iInst) {

  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat_fvm, fem, fem_euler, fem_ns, fem_dg_flow,
  template_solver, disc_adj, disc_adj_fem, disc_adj_turb, disc_adj_heat;
  int val_iter = 0;

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns           = false;  turbulent   = false;
  adj_euler        = false;  adj_ns       = false;  adj_turb    = false;
  fem_euler        = false;  fem_ns       = false;  fem_dg_flow = false;
  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  disc_adj_turb    = false;
  heat_fvm         = false;  disc_adj_heat    = false;
  template_solver  = false;

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  bool restart      = config->GetRestart();
  bool restart_flow = config->GetRestart_Flow();
  bool no_restart   = false;

  /*--- Adjust iteration number for unsteady restarts. ---*/

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool adjoint = (config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint());
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC); // Dynamic simulation (FSI).

  if (dual_time) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
    else val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
  }

  if (time_stepping) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
  }

  /*--- Assign booleans ---*/

  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; break;
    case FEM_LES : fem_ns = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;

  }

  /*--- Determine the kind of FEM solver used for the flow. ---*/

  switch( config->GetKind_FEM_Flow() ) {
    case DG: fem_dg_flow = true; break;
  }

  /*--- Load restarts for any of the active solver containers. Note that
   these restart routines fill the fine grid and interpolate to all MG levels. ---*/

  if (restart || restart_flow) {
    if (euler || ns) {
      solver_container[val_iInst][MESH_0][FLOW_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (turbulent) {
      solver_container[val_iInst][MESH_0][TURB_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (fem) {
      if (dynamic) val_iter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
      solver_container[val_iInst][MESH_0][FEA_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (fem_euler || fem_ns) {
      if (fem_dg_flow)
        solver_container[val_iInst][MESH_0][FLOW_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (heat_fvm) {
      solver_container[val_iInst][MESH_0][HEAT_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
  }

  if (restart) {
    if (template_solver) {
      no_restart = true;
    }
    if (heat_fvm) {
      solver_container[val_iInst][MESH_0][HEAT_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (adj_euler || adj_ns) {
      solver_container[val_iInst][MESH_0][ADJFLOW_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (adj_turb) {
      no_restart = true;
    }
    if (disc_adj) {
      solver_container[val_iInst][MESH_0][ADJFLOW_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
      if (disc_adj_turb)
        solver_container[val_iInst][MESH_0][ADJTURB_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
      if (disc_adj_heat)
        solver_container[val_iInst][MESH_0][ADJHEAT_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (disc_adj_fem) {
        if (dynamic) val_iter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
        solver_container[val_iInst][MESH_0][ADJFEA_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
    if (disc_adj_heat) {
      solver_container[val_iInst][MESH_0][ADJHEAT_SOL]->LoadRestart(geometry[val_iInst], solver_container[val_iInst], config, val_iter, update_geo);
    }
  }

  /*--- Exit if a restart was requested for a solver that is not available. ---*/

  if (no_restart) {
    SU2_MPI::Error(string("A restart capability has not been implemented yet for this solver.\n") +
                   string("Please set RESTART_SOL= NO and try again."), CURRENT_FUNCTION);
  }

  /*--- Think about calls to pre / post-processing here, plus realizability checks. ---*/
  
}

void CErrorEstimationDriver::ComputeMetric() {

  unsigned short iZone, jZone, checkConvergence;
  unsigned long IntIter = 0, nIntIter = 1;
  bool unsteady;

  CSolver *solver_flow = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL],
          *solver_adj  = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];

  if (rank == MASTER_NODE)
    cout << endl <<"----------------------------- Compute Metric ----------------------------" << endl;

  //--- Volume flow grad
  if(rank == MASTER_NODE) cout << "Computing flow volume gradient via L2 Projection." << endl;
  if(nDim == 2)
    solver_flow->SetGradient_L2Proj2(geometry_container[ZONE_0][INST_0][MESH_0], 
                                     config_container[ZONE_0]);
  else
    solver_flow->SetGradient_L2Proj3(geometry_container[ZONE_0][INST_0][MESH_0], 
                                     config_container[ZONE_0]);

  //--- Volume flow Hess
  if(rank == MASTER_NODE) cout << "Computing flow volume Hessian via L2 Projection." << endl;
  if(nDim == 2)
    solver_flow->SetHessian_L2Proj2(geometry_container[ZONE_0][INST_0][MESH_0], 
                                    config_container[ZONE_0]);
  else
    solver_flow->SetHessian_L2Proj3(geometry_container[ZONE_0][INST_0][MESH_0], 
                                    config_container[ZONE_0]);

  //--- Volume adj grad
  if(rank == MASTER_NODE) cout << "Computing adjoint volume gradient via L2 Projection." << endl;
  if(nDim == 2)
    solver_adj->SetGradient_L2Proj2(geometry_container[ZONE_0][INST_0][MESH_0], 
                                    config_container[ZONE_0]);
  else
    solver_adj->SetGradient_L2Proj3(geometry_container[ZONE_0][INST_0][MESH_0], 
                                    config_container[ZONE_0]);

  if(rank == MASTER_NODE) cout << "Computing goal-oriented metric tensor." << endl;
  SumWeightedHessian(solver_flow, solver_adj);

}

void CErrorEstimationDriver::SumWeightedHessian(CSolver* solver_flow,
                                                CSolver* solver_adj) {
  unsigned long iPoint, iElem;
  unsigned long nPointDomain = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain(),
                nElem = geometry_container[ZONE_0][INST_0][MESH_0]->GetnElem();
  unsigned short iVar, iDim, iFlux;
  unsigned short nVarMetr = 4, nFluxMetr = 2;  //--- TODO: adjust size for goal (for different solvers, currently Euler) vs. feature
  unsigned short nMetr = 3*(nDim-1);
  unsigned short nDim = geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    //--- initialize metric tensor
    for(unsigned short im = 0; im < 3*(nDim-1); ++im)
      solver_flow->node[iPoint]->SetAnisoMetr(im, 0.0);

    //--- perform summation of weighted Hessians
    for (unsigned short iVar = 0; iVar < nVarMetr; ++iVar) {

      for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
        const unsigned short ig = iVar*nDim + iFlux;
        const su2double grad = solver_adj->node[iPoint]->GetAnisoGrad(ig);

        for (unsigned short im = 0; im < 3*(nDim-1); ++im) {

          const unsigned short ih = iFlux*nVarMetr*nMetr + iVar*nMetr + im;  
          const su2double hess = solver_flow->node[iPoint]->GetAnisoHess(ih);
          const su2double part = abs(grad)*hess;
          solver_flow->node[iPoint]->AddAnisoMetr(im,part);

        }

      }

    }

  }

  //--- normalize to obtain Lp metric
  su2double p = 1.0;                                                                 // For now, hardcode L1 metric
  su2double Complexity = su2double(config_container[ZONE_0]->GetMesh_Complexity());  // Constraint mesh complexity

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    CVariable *var = solver_flow->node[iPoint];

    const su2double Metr[3] = {var->GetAnisoMetr(0), 
                               var->GetAnisoMetr(1), 
                               var->GetAnisoMetr(2)};

    const su2double DetH = Metr[0]*Metr[2]
                         - Metr[1]*Metr[1];
    const su2double TraH = Metr[0]+Metr[2];

    const su2double Lam1 = TraH/2.0 + sqrt(TraH*TraH/4.-DetH);
    const su2double Lam2 = TraH/2.0 - sqrt(TraH*TraH/4.-DetH);

    const su2double factor = Complexity
                           * pow(abs(Lam1*Lam2), -1./5.);

    var->SetAnisoMetr(0, factor*Metr[0]);
    var->SetAnisoMetr(1, factor*Metr[1]);
    var->SetAnisoMetr(2, factor*Metr[2]);

  }
}

void CErrorEstimationDriver::Run() {
  
  unsigned short iZone = 0, checkConvergence;
  unsigned long IntIter, nIntIter;

  bool unsteady;

  unsteady = (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  /*--- Begin Unsteady pseudo-time stepping internal loop, if not unsteady it does only one step --*/

  if (unsteady)
    nIntIter = config_container[MESH_0]->GetUnst_nIntIter();
  else
    nIntIter = 1;

  for (iZone = 0; iZone < nZone; iZone++) {

    iteration_container[iZone][INST_0]->Preprocess(output, integration_container, geometry_container,
                                                     solver_container, numerics_container, config_container,
                                                     surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }


  /*--- For the adjoint iteration we need the derivatives of the iteration function with
   *    respect to the conservative flow variables. Since these derivatives do not change in the steady state case
   *    we only have to record if the current recording is different from cons. variables. ---*/

  if (RecordingState != FLOW_CONS_VARS || unsteady){

    /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
     *    as argument ensures that all information from a previous recording is removed. ---*/

    SetRecording(NONE);

    /*--- Store the computational graph of one direct iteration with the conservative variables as input. ---*/

    SetRecording(FLOW_CONS_VARS);

  }

  for (IntIter = 0; IntIter < nIntIter; IntIter++) {


    /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {

      config_container[iZone]->SetIntIter(IntIter);

      iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container, config_container, iZone, INST_0);

    }

    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdj_ObjFunction();

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      iteration_container[iZone][INST_0]->Iterate(output, integration_container, geometry_container,
                                          solver_container, numerics_container, config_container,
                                          surface_movement, grid_movement, FFDBox, iZone, INST_0);
    }

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

  }

}

void CErrorEstimationDriver::SetRecording(unsigned short kind_recording){
  unsigned short iZone, iMesh;

  AD::Reset();

  /*--- Prepare for recording by resetting the flow solution to the initial converged solution---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
      solver_container[iZone][INST_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[iZone][INST_0][iMesh], config_container[iZone]);
    }
    if (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS && !config_container[iZone]->GetFrozen_Visc_Disc()) {
      solver_container[iZone][INST_0][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
    }
  }


  /*---Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if (kind_recording != NONE){

    AD::StartRecording();

    if (rank == MASTER_NODE && ((ExtIter == 0)) && kind_recording == FLOW_CONS_VARS) {
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << "Direct iteration to store computational graph." << endl;
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
      cout << "-------------------------------------------------------------------------" << endl << endl;
    }
    for (iZone = 0; iZone < nZone; iZone++) {
      iteration_container[iZone][INST_0]->RegisterInput(solver_container, geometry_container, config_container, iZone, INST_0, kind_recording);
    }

  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->SetDependencies(solver_container, geometry_container, config_container, iZone, INST_0, kind_recording);
  }

  /*--- Do one iteration of the direct flow solver ---*/

  DirectRun();

  /*--- Read the target pressure ---*/

  if (config_container[ZONE_0]->GetInvDesign_Cp() == YES)
    output->SetCp_InverseDesign(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL],
        geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], ExtIter);

  /*--- Read the target heat flux ---*/

  if (config_container[ZONE_0]->GetInvDesign_HeatFlux() == YES)
    output->SetHeatFlux_InverseDesign(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL],
        geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], ExtIter);

  /*--- Print residuals in the first iteration ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    if (rank == MASTER_NODE && ((ExtIter == 0) || (config_container[iZone]->GetUnsteady_Simulation() != STEADY)) && (kind_recording == FLOW_CONS_VARS)) {
      cout << " Zone " << iZone << ": log10[Conservative 0]: "<< log10(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetRes_RMS(0)) << endl;
      if ( config_container[iZone]->GetKind_Turb_Model() != NONE && !config_container[iZone]->GetFrozen_Visc_Disc()) {
        cout <<"       log10[RMS k]: " << log10(solver_container[iZone][INST_0][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
      }
    }
  }

  RecordingState = kind_recording;

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->RegisterOutput(solver_container, geometry_container, config_container, output, iZone, INST_0);
  }

  /*--- Extract the objective function and store it --- */

  SetObjFunction();

  AD::StopRecording();

}

void CErrorEstimationDriver::SetAdj_ObjFunction(){

  bool time_stepping = config_container[ZONE_0]->GetUnsteady_Simulation() != STEADY;
  unsigned long IterAvg_Obj = config_container[ZONE_0]->GetIter_Avg_Objective();
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  su2double seeding = 1.0;

  if (time_stepping){
    if (ExtIter < IterAvg_Obj){
      seeding = 1.0/((su2double)IterAvg_Obj);
    }
    else{
      seeding = 0.0;
    }
  }

  if (rank == MASTER_NODE){
    SU2_TYPE::SetDerivative(ObjFunc, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc, 0.0);
  }

}

void CErrorEstimationDriver::SetObjFunction(){

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
  bool heat         = (config_container[ZONE_0]->GetWeakly_Coupled_Heat());

  ObjFunc = 0.0;

  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->SetTotal_ComboObj(0.0);
  }

  /*--- Specific scalar objective functions ---*/

  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_Solver()) {
      case EULER:                    case NAVIER_STOKES:                   case RANS:
      case DISC_ADJ_EULER:           case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
      case DISC_ADJ_FEM_EULER:       case DISC_ADJ_FEM_NS:                 case DISC_ADJ_FEM_RANS:
        
        if (config_container[ZONE_0]->GetnMarker_Analyze() != 0)
          output->SpecialOutput_AnalyzeSurface(solver_container[iZone][INST_0][MESH_0][FLOW_SOL], geometry_container[iZone][INST_0][MESH_0], config_container[iZone], false);
        
        if ((config_container[ZONE_0]->GetnMarker_Analyze() != 0) && compressible)
          output->SpecialOutput_Distortion(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], false);
        
        if (config_container[ZONE_0]->GetnMarker_NearFieldBound() != 0)
          output->SpecialOutput_SonicBoom(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], false);
          
        if (config_container[ZONE_0]->GetPlot_Section_Forces())
          output->SpecialOutput_SpanLoad(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], false);
        
        break;
    }
  }

  /*--- Surface based obj. function ---*/

  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Evaluate_ObjFunc(config_container[iZone]);
    ObjFunc += solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_ComboObj();
    if (heat){
      if (config_container[iZone]->GetKind_ObjFunc() == TOTAL_HEATFLUX) {
        ObjFunc += solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetTotal_HeatFlux();
      }
      else if (config_container[iZone]->GetKind_ObjFunc() == TOTAL_AVG_TEMPERATURE) {
        ObjFunc += solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetTotal_AvgTemperature();
      }
    }
  }

  if (rank == MASTER_NODE){
    AD::RegisterOutput(ObjFunc);
  }
  
}

void CErrorEstimationDriver::DirectRun(){


  unsigned short iZone, jZone;
  bool unsteady = config_container[ZONE_0]->GetUnsteady_Simulation() != STEADY;

  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  unsteady = (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  /*--- Zone preprocessing ---*/

  for (iZone = 0; iZone < nZone; iZone++)
    direct_iteration[iZone]->Preprocess(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

  /*--- For each zone runs one single iteration ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]->SetIntIter(1);
    direct_iteration[iZone]->Iterate(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

}

void CErrorEstimationDriver::Output() {

  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter_OffSet();

  /*--- Set Kind_Solver to primal since adaptation requires Mach (or pressure) ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
    case DISC_ADJ_EULER: config_container[ZONE_0]->SetKind_Solver(EULER); break;
    case DISC_ADJ_NAVIER_STOKES: config_container[ZONE_0]->SetKind_Solver(NAVIER_STOKES); break;
    case DISC_ADJ_RANS: config_container[ZONE_0]->SetKind_Solver(RANS); break;
  }

  /*--- Set DiscreteAdjoint flag to false so we write to correct files ---*/

  config_container[ZONE_0]->SetDiscrete_Adjoint(false);

  if (rank == MASTER_NODE) cout << endl << "-------------------------- File Output Summary --------------------------";

    /*--- Execute the routine for writing restart, volume solution,
     surface solution, and surface comma-separated value files. ---*/

  output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, ExtIter, nZone);

  if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl << endl;

  
}

void CErrorEstimationDriver::Postprocessing() {

  if (rank == MASTER_NODE)
    cout << endl <<"------------------- Error Estimation Postprocessing ---------------------" << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Numerics_Postprocessing(numerics_container[iZone], solver_container[iZone][iInst],
          geometry_container[iZone][iInst], config_container[iZone], iInst);
    } // End of loop over iInst
    delete [] numerics_container[iZone];
  }
  delete [] numerics_container;
  if (rank == MASTER_NODE) cout << "Deleted CNumerics container." << endl;
  
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Integration_Postprocessing(integration_container[iZone],
          geometry_container[iZone][iInst],
          config_container[iZone],
          iInst);
    } // End of loop over iInst
    delete [] integration_container[iZone];
  }
  delete [] integration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIntegration container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      Solver_Postprocessing(solver_container[iZone], geometry_container[iZone][iInst],
                            config_container[iZone], iInst);
    } // End of loop over iInst
  }
  if (rank == MASTER_NODE) cout << "Deleted CSolver containers." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
	for (iInst = 0; iInst < nInst[iZone]; iInst++)
    delete iteration_container[iZone][iInst];
    delete [] iteration_container[iZone];
  }
  delete [] iteration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIteration container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    if (geometry_container[iZone] != NULL) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        for (unsigned short iMGlevel = 0; iMGlevel < config_container[iZone]->GetnMGLevels()+1; iMGlevel++) {
          if (geometry_container[iZone][iInst][iMGlevel] != NULL) delete geometry_container[iZone][iInst][iMGlevel];
        }
        if (geometry_container[iZone][iInst] != NULL) delete [] geometry_container[iZone][iInst];
      } // End of loop over iInst
      delete [] geometry_container[iZone];
    }
  }
  delete [] geometry_container;
  if (rank == MASTER_NODE) cout << "Deleted CGeometry containers." << endl;

  /*--- Deallocate config container ---*/
  if (config_container!= NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != NULL) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig containers." << endl;

}

void CErrorEstimationDriver::Solver_Postprocessing(CSolver ****solver_container, CGeometry **geometry,
                                                   CConfig *config, unsigned short val_iInst) {
  
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat_fvm, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb, disc_adj_fem, disc_adj_heat,
  e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;

  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj        = false;
  fem              = false;  disc_adj_fem    = false;
  heat_fvm        = false;   disc_adj_heat   = false;
  transition       = false;
  template_solver  = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;

  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : euler = true; break;
    case FEM_NAVIER_STOKES:
    case FEM_LES: ns = true; break;
    case FEM_RANS: ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
    case SA:     spalart_allmaras = true;     break;
    case SA_NEG: neg_spalart_allmaras = true; break;
    case SST:    menter_sst = true;           break;
    case SA_E: e_spalart_allmaras = true; break;
    case SA_COMP: comp_spalart_allmaras = true; break;
    case SA_E_COMP: e_comp_spalart_allmaras = true; break;
    }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    
    /*--- DeAllocate solution for a template problem ---*/
    
    if (template_solver) {
      delete solver_container[val_iInst][iMGlevel][TEMPLATE_SOL];
    }

    /*--- DeAllocate solution for adjoint problem ---*/
    
    if (adj_euler || adj_ns || disc_adj) {
      delete solver_container[val_iInst][iMGlevel][ADJFLOW_SOL];
      if (disc_adj_turb || adj_turb) {
        delete solver_container[val_iInst][iMGlevel][ADJTURB_SOL];
      }
      if (heat_fvm) {
        delete solver_container[val_iInst][iMGlevel][ADJHEAT_SOL];
      }
    }

    if (disc_adj_heat) {
      delete solver_container[val_iInst][iMGlevel][ADJHEAT_SOL];
    }

    /*--- DeAllocate solution for direct problem ---*/
    
    if (euler || ns) {
      delete solver_container[val_iInst][iMGlevel][FLOW_SOL];
    }

    if (turbulent) {
      if (spalart_allmaras || neg_spalart_allmaras || menter_sst || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
        delete solver_container[val_iInst][iMGlevel][TURB_SOL];
      }
      if (transition) {
        delete solver_container[val_iInst][iMGlevel][TRANS_SOL];
      }
    }
    if (heat_fvm) {
      delete solver_container[val_iInst][iMGlevel][HEAT_SOL];
    }
    if (fem) {
      delete solver_container[val_iInst][iMGlevel][FEA_SOL];
    }
    if (disc_adj_fem) {
      delete solver_container[val_iInst][iMGlevel][ADJFEA_SOL];
    }
    
    delete [] solver_container[val_iInst][iMGlevel];
  }
  
  delete [] solver_container[val_iInst];

}

void CErrorEstimationDriver::Integration_Postprocessing(CIntegration ***integration_container,
    CGeometry **geometry, CConfig *config, unsigned short val_iInst) {
  bool euler, adj_euler, ns, adj_ns, turbulent, adj_turb, fem,
      fem_euler, fem_ns, fem_turbulent,
      heat_fvm, template_solver, transition, disc_adj, disc_adj_fem, disc_adj_heat;

  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false;
  ns               = false; adj_ns           = false;
  turbulent        = false; adj_turb         = false;
  disc_adj         = false;
  fem_euler        = false;
  fem_ns           = false;
  fem_turbulent    = false;
  heat_fvm         = false; disc_adj_heat    = false;
  fem              = false; disc_adj_fem     = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; fem_turbulent = true; break;
    case FEM_LES :  fem_ns = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER : euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS : ns = true; turbulent = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; fem_turbulent = true; disc_adj = true; break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }

  /*--- DeAllocate solution for a template problem ---*/
  if (template_solver) integration_container[val_iInst][TEMPLATE_SOL] = new CSingleGridIntegration(config);

  /*--- DeAllocate solution for direct problem ---*/
  if (euler || ns) delete integration_container[val_iInst][FLOW_SOL];
  if (turbulent) delete integration_container[val_iInst][TURB_SOL];
  if (transition) delete integration_container[val_iInst][TRANS_SOL];
  if (heat_fvm) delete integration_container[val_iInst][HEAT_SOL];
  if (fem) delete integration_container[val_iInst][FEA_SOL];
  if (disc_adj_fem) delete integration_container[val_iInst][ADJFEA_SOL];
  if (disc_adj_heat) delete integration_container[val_iInst][ADJHEAT_SOL];

  /*--- DeAllocate solution for adjoint problem ---*/
  if (adj_euler || adj_ns || disc_adj) delete integration_container[val_iInst][ADJFLOW_SOL];
  if (adj_turb) delete integration_container[val_iInst][ADJTURB_SOL];

  /*--- DeAllocate integration container for finite element flow solver. ---*/
  if (fem_euler || fem_ns) delete integration_container[val_iInst][FLOW_SOL];
  //if (fem_turbulent)     delete integration_container[val_iInst][FEM_TURB_SOL];

  if (fem_turbulent)
    SU2_MPI::Error("No turbulent FEM solver yet", CURRENT_FUNCTION);

  delete [] integration_container[val_iInst];
}

void CErrorEstimationDriver::Numerics_Postprocessing(CNumerics *****numerics_container,
                                      CSolver ***solver_container, CGeometry **geometry,
                                      CConfig *config, unsigned short val_iInst) {
  
  unsigned short iMGlevel, iSol;
  
  
  bool
  euler, adj_euler,
  ns, adj_ns,
  fem_euler, fem_ns, fem_turbulent,
  turbulent, adj_turb,
  spalart_allmaras, neg_spalart_allmaras, menter_sst,
  fem,
  heat_fvm,
  transition,
  template_solver;

  bool e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;

  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Initialize some useful booleans ---*/
  euler            = false; ns     = false; turbulent     = false;
  fem_euler        = false; fem_ns = false; fem_turbulent = false;
  adj_euler        = false;   adj_ns           = false;   adj_turb         = false;
  fem        = false;
  spalart_allmaras = false;   neg_spalart_allmaras = false; menter_sst       = false;
  transition       = false;   heat_fvm         = false;
  template_solver  = false;
    
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case DISC_ADJ_EULER: euler = true;  heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES: ns = true;  heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : case DISC_ADJ_RANS:  ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case FEM_EULER : case DISC_ADJ_FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: case DISC_ADJ_FEM_NS : fem_ns = true; break;
    case FEM_RANS : case DISC_ADJ_FEM_RANS : fem_ns = true; fem_turbulent = true; break;
    case FEM_LES :  fem_ns = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: case DISC_ADJ_FEM: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
  }
  
  /*--- Assign turbulence model booleans ---*/

  if (turbulent || fem_turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;  break;
      case SA_COMP: comp_spalart_allmaras = true; break;
      case SA_E: e_spalart_allmaras = true; break;
      case SA_E_COMP: e_comp_spalart_allmaras = true; break;

    }
  
  /*--- Solver definition for the template problem ---*/
  if (template_solver) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          delete numerics_container[val_iInst][iMGlevel][TEMPLATE_SOL][CONV_TERM];
        break;
    }
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      /*--- Definition of the viscous scheme for each equation and mesh level ---*/
      delete numerics_container[val_iInst][iMGlevel][TEMPLATE_SOL][VISC_TERM];
      /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
      delete numerics_container[val_iInst][iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM];
      /*--- Definition of the boundary condition method ---*/
      delete numerics_container[val_iInst][iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM];
    }
    
  }
  
  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
        
      case SPACE_CENTERED :
        if (compressible) {
          
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_Flow()) {
            case LAX : case JST :  case JST_KE : delete numerics_container[val_iInst][MESH_0][FLOW_SOL][CONV_TERM]; break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM];
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use preconditioning method ---*/
          switch (config->GetKind_Centered_Flow()) {
            case LAX : case JST : delete numerics_container[val_iInst][MESH_0][FLOW_SOL][CONV_TERM]; break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM];
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
          
        }
        break;
      case SPACE_UPWIND :
        
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE: case AUSM : case TURKEL: case HLLC: case MSW:  case CUSP: case L2ROE: case LMROE: case SLAU: case SLAU2: case AUSMPLUSUP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM];
                delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
              }
              
              break;
          }
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use preconditioning method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case FDS:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM];
                delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
              }
              break;
          }
        }
        
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible||incompressible) {
      /*--- Compressible flow Ideal gas ---*/
      delete numerics_container[val_iInst][MESH_0][FLOW_SOL][VISC_TERM];
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][VISC_TERM];
      
      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][VISC_BOUND_TERM];
      
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM];
      delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM];
    }
    
  }

  /*--- DG-FEM solver definition for Euler, Navier-Stokes problems ---*/

  if ((fem_euler) || (fem_ns)) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetRiemann_Solver_FEM()) {
      case AUSM: case TURKEL: case HLLC: case MSW: /* Note that not all need to be deleted. */

        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_TERM];
          delete numerics_container[val_iInst][iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
        }
        break;
    }
  }

  /*--- Solver definition for the turbulent model problem ---*/
  
  if (turbulent) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          if (spalart_allmaras || neg_spalart_allmaras ||menter_sst|| comp_spalart_allmaras || e_spalart_allmaras || e_comp_spalart_allmaras)
            delete numerics_container[val_iInst][iMGlevel][TURB_SOL][CONV_TERM];
        }
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    
      if (spalart_allmaras || neg_spalart_allmaras ||menter_sst|| comp_spalart_allmaras || e_spalart_allmaras || e_comp_spalart_allmaras){
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          delete numerics_container[val_iInst][iMGlevel][TURB_SOL][VISC_TERM];
          delete numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_FIRST_TERM];
          delete numerics_container[val_iInst][iMGlevel][TURB_SOL][SOURCE_SECOND_TERM];
          /*--- Definition of the boundary condition method ---*/
          delete numerics_container[val_iInst][iMGlevel][TURB_SOL][CONV_BOUND_TERM];
          delete numerics_container[val_iInst][iMGlevel][TURB_SOL][VISC_BOUND_TERM];

      }
    }
    
  }
  
  /*--- Solver definition for the transition model problem ---*/
  if (transition) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          delete numerics_container[val_iInst][iMGlevel][TRANS_SOL][CONV_TERM];
        }
        break;
    }
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      /*--- Definition of the viscous scheme for each equation and mesh level ---*/
      delete numerics_container[val_iInst][iMGlevel][TRANS_SOL][VISC_TERM];
      /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
      delete numerics_container[val_iInst][iMGlevel][TRANS_SOL][SOURCE_FIRST_TERM];
      delete numerics_container[val_iInst][iMGlevel][TRANS_SOL][SOURCE_SECOND_TERM];
      /*--- Definition of the boundary condition method ---*/
      delete numerics_container[val_iInst][iMGlevel][TRANS_SOL][CONV_BOUND_TERM];
    }
  }

  if (heat_fvm) {

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      delete numerics_container[val_iInst][iMGlevel][HEAT_SOL][VISC_TERM];
      delete numerics_container[val_iInst][iMGlevel][HEAT_SOL][VISC_BOUND_TERM];

      switch (config->GetKind_ConvNumScheme_Heat()) {
        case SPACE_UPWIND :

          delete numerics_container[val_iInst][iMGlevel][HEAT_SOL][CONV_TERM];
          delete numerics_container[val_iInst][iMGlevel][HEAT_SOL][CONV_BOUND_TERM];
          break;

        case SPACE_CENTERED :

          delete numerics_container[val_iInst][iMGlevel][HEAT_SOL][CONV_TERM];
          delete numerics_container[val_iInst][iMGlevel][HEAT_SOL][CONV_BOUND_TERM];
        break;
      }
    }
  }
  
  /*--- Solver definition for the flow adjoint problem ---*/
  
  if (adj_euler || adj_ns ) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    
    switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case SPACE_CENTERED :
        
        if (compressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Centered_AdjFlow()) {
            case LAX : case JST:
              delete numerics_container[val_iInst][MESH_0][ADJFLOW_SOL][CONV_TERM];
              break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_TERM];
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];
          
        }
        
        if (incompressible) {
          
          /*--- Incompressible flow, use artificial compressibility method ---*/
          
          switch (config->GetKind_Centered_AdjFlow()) {
            case LAX : case JST:
              delete numerics_container[val_iInst][MESH_0][ADJFLOW_SOL][CONV_TERM]; break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_TERM];
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];
          
        }
        
        break;
        
      case SPACE_UPWIND :
        
        if (compressible || incompressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Upwind_AdjFlow()) {
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_TERM];
                delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];
              }
              break;
          }
        }
        
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    
    if (compressible || incompressible) {
      
      /*--- Compressible flow ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][VISC_TERM];
        delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM];
      }
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      
      
      if (compressible || incompressible) {
        
        delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM];
        delete numerics_container[val_iInst][iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM];
        
      }
    }
    
  }
  
  
  /*--- Solver definition for the turbulent adjoint problem ---*/
  if (adj_turb) {
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTurb()) {
        
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          if (spalart_allmaras) {
            delete numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][CONV_TERM];
          }
        break;
    }
    
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        /*--- Definition of the viscous scheme for each equation and mesh level ---*/
        delete numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][VISC_TERM];
        /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
        delete numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][SOURCE_FIRST_TERM];
        delete numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][SOURCE_SECOND_TERM];
        /*--- Definition of the boundary condition method ---*/
        delete numerics_container[val_iInst][iMGlevel][ADJTURB_SOL][CONV_BOUND_TERM];
      }
    }
  }
  
  /*--- Solver definition for the FEA problem ---*/
  if (fem) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    delete numerics_container[val_iInst][MESH_0][FEA_SOL][FEA_TERM];
    
  }
  
  /*--- Definition of the Class for the numerical method: numerics_container[INST_LEVEL][MESH_LEVEL][EQUATION][EQ_TERM] ---*/
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    for (iSol = 0; iSol < MAX_SOLS; iSol++) {
      delete [] numerics_container[val_iInst][iMGlevel][iSol];
    }
    delete[] numerics_container[val_iInst][iMGlevel];
  }
  
  delete[] numerics_container[val_iInst];

}

void CErrorEstimationDriver::Solver_Deletion(CSolver ****solver_container,
                                    CConfig *config, unsigned short val_iInst) {
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat_fvm, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb, disc_adj_fem, disc_adj_heat,
  e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;

  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj        = false;
  fem              = false;  disc_adj_fem    = false;
  heat_fvm        = false;   disc_adj_heat   = false;
  transition       = false;
  template_solver  = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;

  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : euler = true; break;
    case FEM_NAVIER_STOKES:
    case FEM_LES: ns = true; break;
    case FEM_RANS: ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent){
    switch (config->GetKind_Turb_Model()) {
    case SA:     spalart_allmaras = true;     break;
    case SA_NEG: neg_spalart_allmaras = true; break;
    case SST:    menter_sst = true;           break;
    case SA_E: e_spalart_allmaras = true; break;
    case SA_COMP: comp_spalart_allmaras = true; break;
    case SA_E_COMP: e_comp_spalart_allmaras = true; break;
    }
  }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  
    
    /*--- DeAllocate solution for a template problem ---*/
    
  if (template_solver) {
    delete solver_container[val_iInst][MESH_0][TEMPLATE_SOL];
  }

    /*--- DeAllocate solution for adjoint problem ---*/

  if (adj_euler || adj_ns || disc_adj) {
    delete solver_container[val_iInst][MESH_0][ADJFLOW_SOL];
    if (disc_adj_turb || adj_turb) {
      delete solver_container[val_iInst][MESH_0][ADJTURB_SOL];
    }
    if (heat_fvm) {
      delete solver_container[val_iInst][MESH_0][ADJHEAT_SOL];
    }
  }

  if (disc_adj_heat) {
    delete solver_container[val_iInst][MESH_0][ADJHEAT_SOL];
  }

    /*--- DeAllocate solution for direct problem ---*/

  if (euler || ns) {
    delete solver_container[val_iInst][MESH_0][FLOW_SOL];
  }

  if (turbulent) {
    if (spalart_allmaras || neg_spalart_allmaras || menter_sst || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
      delete solver_container[val_iInst][MESH_0][TURB_SOL];
    }
    if (transition) {
      delete solver_container[val_iInst][MESH_0][TRANS_SOL];
    }
  }
  if (heat_fvm) {
    delete solver_container[val_iInst][MESH_0][HEAT_SOL];
  }
  if (fem) {
    delete solver_container[val_iInst][MESH_0][FEA_SOL];
  }
  if (disc_adj_fem) {
    delete solver_container[val_iInst][MESH_0][ADJFEA_SOL];
  }

  delete [] solver_container[val_iInst][MESH_0];

  delete [] solver_container[val_iInst];

}