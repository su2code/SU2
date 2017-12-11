/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Economon, H. Kline, R. Sanchez, F. Palacios
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/driver_structure.hpp"
#include "../include/definition_structure.hpp"

CDriver::CDriver(char* confFile,
                 unsigned short val_nZone,
                 unsigned short val_nDim,
                 SU2_Comm MPICommunicator):config_file_name(confFile), StartTime(0.0), StopTime(0.0), UsedTime(0.0), ExtIter(0), nZone(val_nZone), nDim(val_nDim), StopCalc(false), fsi(false) {


  unsigned short jZone, iSol;
  unsigned short Kind_Grid_Movement;
  bool initStaticMovement;



  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPICommunicator, &rank);
#endif

  /*--- Create pointers to all of the classes that may be used throughout
   the SU2_CFD code. In general, the pointers are instantiated down a
   hierarchy over all zones, multigrid levels, equation sets, and equation
   terms as described in the comments below. ---*/

  iteration_container            = NULL;
  output                         = NULL;
  integration_container          = NULL;
  geometry_container             = NULL;
  solver_container               = NULL;
  numerics_container             = NULL;
  config_container               = NULL;
  surface_movement               = NULL;
  grid_movement                  = NULL;
  FFDBox                         = NULL;
  interpolator_container         = NULL;
  transfer_container             = NULL;

  /*--- Definition and of the containers for all possible zones. ---*/

  iteration_container            = new CIteration*[nZone];
  solver_container               = new CSolver***[nZone];
  integration_container          = new CIntegration**[nZone];
  numerics_container             = new CNumerics****[nZone];
  config_container               = new CConfig*[nZone];
  geometry_container             = new CGeometry**[nZone];
  surface_movement               = new CSurfaceMovement*[nZone];
  grid_movement                  = new CVolumetricMovement*[nZone];
  FFDBox                         = new CFreeFormDefBox**[nZone];
  interpolator_container         = new CInterpolator**[nZone];
  transfer_container             = new CTransfer**[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone]               = NULL;
    integration_container[iZone]          = NULL;
    numerics_container[iZone]             = NULL;
    config_container[iZone]               = NULL;
    geometry_container[iZone]             = NULL;
    surface_movement[iZone]               = NULL;
    grid_movement[iZone]                  = NULL;
    FFDBox[iZone]                         = NULL;
    interpolator_container[iZone]         = NULL;
    transfer_container[iZone]             = NULL;
  }

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    config_container[iZone] = new CConfig(config_file_name, SU2_CFD, iZone, nZone, nDim, VERB_HIGH);

    /*--- Set the MPI communicator ---*/

    config_container[iZone]->SetMPICommunicator(MPICommunicator);

    /*--- Definition of the geometry class to store the primal grid in the
     partitioning process. ---*/

    CGeometry *geometry_aux = NULL;

    /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

    geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

    /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

    geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

    /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

    geometry_container[iZone] = new CGeometry *[config_container[iZone]->GetnMGLevels()+1];
    geometry_container[iZone][MESH_0] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);

    /*--- Deallocate the memory of geometry_aux ---*/

    delete geometry_aux;

    /*--- Add the Send/Receive boundaries ---*/

    geometry_container[iZone][MESH_0]->SetSendReceive(config_container[iZone]);

    /*--- Add the Send/Receive boundaries ---*/

    geometry_container[iZone][MESH_0]->SetBoundaries(config_container[iZone]);

  }

  /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
   based data structure is constructed, i.e. node and cell neighbors are
   identified and linked, face areas and volumes of the dual mesh cells are
   computed, and the multigrid levels are created using an agglomeration procedure. ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Geometry Preprocessing ------------------------" << endl;

  Geometrical_Preprocessing();

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Computation of wall distances for turbulence modeling ---*/

    if (rank == MASTER_NODE)
      cout << "Computing wall distances." << endl;

    if ((config_container[iZone]->GetKind_Solver() == RANS) ||
        (config_container[iZone]->GetKind_Solver() == ADJ_RANS) ||
        (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS))
      geometry_container[iZone][MESH_0]->ComputeWall_Distance(config_container[iZone]);

    /*--- Computation of positive surface area in the z-plane which is used for
     the calculation of force coefficient (non-dimensionalization). ---*/

    geometry_container[iZone][MESH_0]->SetPositive_ZArea(config_container[iZone]);

    /*--- Set the near-field, interface and actuator disk boundary conditions, if necessary. ---*/

    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
      geometry_container[iZone][iMesh]->MatchNearField(config_container[iZone]);
      geometry_container[iZone][iMesh]->MatchInterface(config_container[iZone]);
      geometry_container[iZone][iMesh]->MatchActuator_Disk(config_container[iZone]);
    }

  }

  /*--- If activated by the compile directive, perform a partition analysis. ---*/
#if PARTITION
  Partition_Analysis(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
#endif

  /*--- Output some information about the driver that has been instantiated for the problem. ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Driver information --------------------------" << endl;

  fsi = config_container[ZONE_0]->GetFSI_Simulation();
  bool stat_fsi = ((config_container[ZONE_0]->GetDynamic_Analysis() == STATIC) && (config_container[ZONE_0]->GetUnsteady_Simulation() == STEADY));
  bool disc_adj_fsi = (config_container[ZONE_0]->GetDiscrete_Adjoint());

  if ( (config_container[ZONE_0]->GetKind_Solver() == FEM_ELASTICITY ||
        config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM ||
        config_container[ZONE_0]->GetKind_Solver() == POISSON_EQUATION ||
        config_container[ZONE_0]->GetKind_Solver() == WAVE_EQUATION ||
        config_container[ZONE_0]->GetKind_Solver() == HEAT_EQUATION) ) {
    if (rank == MASTER_NODE) cout << "A General driver has been instantiated." << endl;
  }
  else if (config_container[ZONE_0]->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    if (rank == MASTER_NODE) cout << "A Harmonic Balance driver has been instantiated." << endl;
  }
  else if (nZone == 2 && fsi) {
    if (disc_adj_fsi) {
      if (stat_fsi)
        if (rank == MASTER_NODE) cout << "A Discrete-Adjoint driver for Fluid-Structure Interaction has been instantiated." << endl;
    }
    else{
      if (stat_fsi){if (rank == MASTER_NODE) cout << "A Static Fluid-Structure Interaction driver has been instantiated." << endl;}
      else{if (rank == MASTER_NODE) cout << "A Dynamic Fluid-Structure Interaction driver has been instantiated." << endl;}
    }

  }
  else {
    if (rank == MASTER_NODE) cout << "A Fluid driver has been instantiated." << endl;
  }

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Instantiate the type of physics iteration to be executed within each zone. For
     example, one can execute the same physics across multiple zones (mixing plane),
     different physics in different zones (fluid-structure interaction), or couple multiple
     systems tightly within a single zone by creating a new iteration class (e.g., RANS). ---*/
    
    if (rank == MASTER_NODE) {
      cout << endl <<"------------------------ Iteration Preprocessing ------------------------" << endl;
    }
    Iteration_Preprocessing();

    /*--- Definition of the solver class: solver_container[#ZONES][#MG_GRIDS][#EQ_SYSTEMS].
     The solver classes are specific to a particular set of governing equations,
     and they contain the subroutines with instructions for computing each spatial
     term of the PDE, i.e. loops over the edges to compute convective and viscous
     fluxes, loops over the nodes to compute source terms, and routines for
     imposing various boundary condition type for the PDE. ---*/

    if (rank == MASTER_NODE)
      cout << endl <<"------------------------- Solver Preprocessing --------------------------" << endl;

    solver_container[iZone] = new CSolver** [config_container[iZone]->GetnMGLevels()+1];
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++)
      solver_container[iZone][iMesh] = NULL;

    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
      solver_container[iZone][iMesh] = new CSolver* [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        solver_container[iZone][iMesh][iSol] = NULL;
    }
    Solver_Preprocessing(solver_container[iZone], geometry_container[iZone],
        config_container[iZone]);

    if (rank == MASTER_NODE)
      cout << endl <<"----------------- Integration and Numerics Preprocessing ----------------" << endl;

    /*--- Definition of the integration class: integration_container[#ZONES][#EQ_SYSTEMS].
     The integration class orchestrates the execution of the spatial integration
     subroutines contained in the solver class (including multigrid) for computing
     the residual at each node, R(U) and then integrates the equations to a
     steady state or time-accurately. ---*/

    integration_container[iZone] = new CIntegration*[MAX_SOLS];
    Integration_Preprocessing(integration_container[iZone], geometry_container[iZone],
                              config_container[iZone]);

    
    if (rank == MASTER_NODE) cout << "Integration Preprocessing." << endl;

    /*--- Definition of the numerical method class:
     numerics_container[#ZONES][#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
     The numerics class contains the implementation of the numerical methods for
     evaluating convective or viscous fluxes between any two nodes in the edge-based
     data structure (centered, upwind, galerkin), as well as any source terms
     (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/

    numerics_container[iZone] = new CNumerics***[config_container[iZone]->GetnMGLevels()+1];
    Numerics_Preprocessing(numerics_container[iZone], solver_container[iZone],
        geometry_container[iZone], config_container[iZone]);

    if (rank == MASTER_NODE) cout << "Numerics Preprocessing." << endl;

  }

  /*--- Definition of the interface and transfer conditions between different zones.
   *--- The transfer container is defined for zones paired one to one.
   *--- This only works for a multizone FSI problem (nZone > 1).
   *--- Also, at the moment this capability is limited to two zones (nZone < 3).
   *--- This will change in the future. ---*/

  if ((rank == MASTER_NODE) && nZone > 1)
    cout << endl <<"------------------- Multizone Interface Preprocessing -------------------" << endl;

  if ( nZone > 1 ) {
    for (iZone = 0; iZone < nZone; iZone++){
      transfer_container[iZone] = new CTransfer*[nZone];
      interpolator_container[iZone] = new CInterpolator*[nZone];
      for (jZone = 0; jZone < nZone; jZone++){
        transfer_container[iZone][jZone]             = NULL;
        interpolator_container[iZone][jZone]         = NULL;
      }
    }

    Interface_Preprocessing();
  }

  /*--- Instantiate the geometry movement classes for the solution of unsteady
   flows on dynamic meshes, including rigid mesh transformations, dynamically
   deforming meshes, and preprocessing of harmonic balance. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    if (config_container[iZone]->GetGrid_Movement() ||
        (config_container[iZone]->GetDirectDiff() == D_DESIGN)) {
      if (rank == MASTER_NODE)
        cout << "Setting dynamic mesh structure for zone "<< iZone<<"." << endl;
      grid_movement[iZone] = new CVolumetricMovement(geometry_container[iZone][MESH_0], config_container[iZone]);
      FFDBox[iZone] = new CFreeFormDefBox*[MAX_NUMBER_FFD];
      surface_movement[iZone] = new CSurfaceMovement();
      surface_movement[iZone]->CopyBoundary(geometry_container[iZone][MESH_0], config_container[iZone]);
      if (config_container[iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE)
        iteration_container[iZone]->SetGrid_Movement(geometry_container, surface_movement, grid_movement, FFDBox, solver_container, config_container, iZone, 0, 0);
    }

    if (config_container[iZone]->GetDirectDiff() == D_DESIGN) {
      if (rank == MASTER_NODE)
        cout << "Setting surface/volume derivatives." << endl;

      /*--- Set the surface derivatives, i.e. the derivative of the surface mesh nodes with respect to the design variables ---*/

      surface_movement[iZone]->SetSurface_Derivative(geometry_container[iZone][MESH_0],config_container[iZone]);

      /*--- Call the volume deformation routine with derivative mode enabled.
       This computes the derivative of the volume mesh with respect to the surface nodes ---*/

      grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone][MESH_0],config_container[iZone], true, true);

      /*--- Update the multi-grid structure to propagate the derivative information to the coarser levels ---*/

      geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);

      /*--- Set the derivative of the wall-distance with respect to the surface nodes ---*/

      if ( (config_container[iZone]->GetKind_Solver() == RANS) ||
          (config_container[iZone]->GetKind_Solver() == ADJ_RANS) ||
          (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS))
        geometry_container[iZone][MESH_0]->ComputeWall_Distance(config_container[iZone]);
    }

    if (config_container[iZone]->GetKind_GridMovement(iZone) == FLUID_STRUCTURE_STATIC){
      if (rank == MASTER_NODE)
        cout << "Setting moving mesh structure for static FSI problems." << endl;
        /*--- Instantiate the container for the grid movement structure ---*/
        grid_movement[iZone]    = new CElasticityMovement(geometry_container[iZone][MESH_0], config_container[iZone]);
    }

  }

  if(fsi && (config_container[ZONE_0]->GetRestart() || config_container[ZONE_0]->GetDiscrete_Adjoint())){
    if (rank == MASTER_NODE)cout << endl <<"Restarting Fluid and Structural Solvers." << endl;

    for (iZone = 0; iZone < nZone; iZone++) {
        Solver_Restart(solver_container[iZone], geometry_container[iZone],
                       config_container[iZone], true);
    }

  }

  /*---If the Grid Movement is static initialize the static mesh movment ---*/
  Kind_Grid_Movement = config_container[ZONE_0]->GetKind_GridMovement(ZONE_0);
  initStaticMovement = (config_container[ZONE_0]->GetGrid_Movement() && (Kind_Grid_Movement == MOVING_WALL
                        || Kind_Grid_Movement == ROTATING_FRAME || Kind_Grid_Movement == STEADY_TRANSLATION));


  if(initStaticMovement){
    if (rank == MASTER_NODE)cout << endl <<"--------------------- Initialize Static Mesh Movement --------------------" << endl;

      InitStaticMeshMovement();
  }

 if (config_container[ZONE_0]->GetBoolTurbomachinery()){
   if (rank == MASTER_NODE)cout << endl <<"---------------------- Turbomachinery Preprocessing ---------------------" << endl;
      TurbomachineryPreprocessing();
  }


  /*--- Definition of the output class (one for all zones). The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/

  output = new COutput(config_container[ZONE_0]);

  /*--- Open the convergence history file ---*/

  if (rank == MASTER_NODE){
    ConvHist_file = new ofstream[nZone];
    for (iZone = 0; iZone < nZone; iZone++) {
      output->SetConvHistory_Header(&ConvHist_file[iZone], config_container[ZONE_0], iZone);
    }
  }
  /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
  if (config_container[ZONE_0]->GetWrt_Unsteady() && config_container[ZONE_0]->GetRestart())
    ExtIter = config_container[ZONE_0]->GetUnst_RestartIter();

  /*--- Check for a dynamic restart (structural analysis). Update ExtIter if necessary. ---*/
  if (config_container[ZONE_0]->GetKind_Solver() == FEM_ELASTICITY
      && config_container[ZONE_0]->GetWrt_Dynamic() && config_container[ZONE_0]->GetRestart())
    ExtIter = config_container[ZONE_0]->GetDyn_RestartIter();

  /*--- Initialize some variables used for external communications trough the Py wrapper. ---*/
  PyWrapVarCoord[0] = 0.0;
  PyWrapVarCoord[1] = 0.0;
  PyWrapVarCoord[2] = 0.0;
  PyWrapNodalForce[0] = 0.0;
  PyWrapNodalForce[1] = 0.0;
  PyWrapNodalForce[2] = 0.0;
  PyWrapNodalForceDensity[0] = 0.0;
  PyWrapNodalForceDensity[1] = 0.0;
  PyWrapNodalForceDensity[2] = 0.0;

  /*--- Set up a timer for performance benchmarking (preprocessing time is not included) ---*/

#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif

}

void CDriver::Postprocessing() {

  int rank = MASTER_NODE;
  int size = SINGLE_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    /*--- Output some information to the console. ---*/

  if (rank == MASTER_NODE) {

    /*--- Print out the number of non-physical points and reconstructions ---*/

    if (config_container[ZONE_0]->GetNonphysical_Points() > 0)
      cout << "Warning: there are " << config_container[ZONE_0]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
    if (config_container[ZONE_0]->GetNonphysical_Reconstr() > 0)
      cout << "Warning: " << config_container[ZONE_0]->GetNonphysical_Reconstr() << " reconstructed states for upwinding are non-physical." << endl;

    /*--- Close the convergence history file. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      ConvHist_file[iZone].close();
    }

  }

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
     Numerics_Postprocessing(numerics_container[iZone], solver_container[iZone],
     geometry_container[iZone], config_container[iZone]);
    delete [] numerics_container[iZone];
  }
  delete [] numerics_container;
  if (rank == MASTER_NODE) cout << "Deleted CNumerics container." << endl;
  
  for (iZone = 0; iZone < nZone; iZone++) {
    Integration_Postprocessing(integration_container[iZone],
                               geometry_container[iZone],
                               config_container[iZone]);
    delete [] integration_container[iZone];
  }
  delete [] integration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIntegration container." << endl;
  
  for (iZone = 0; iZone < nZone; iZone++) {
    Solver_Postprocessing(solver_container[iZone],
                          geometry_container[iZone],
                          config_container[iZone]);
    delete [] solver_container[iZone];
  }
  delete [] solver_container;
  if (rank == MASTER_NODE) cout << "Deleted CSolver container." << endl;
  
  for (iZone = 0; iZone < nZone; iZone++) {
    delete iteration_container[iZone];
  }
  delete [] iteration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIteration container." << endl;
  
  if (interpolator_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
    if (interpolator_container[iZone] != NULL){
            delete [] interpolator_container[iZone];
      }
    }
    delete [] interpolator_container;
    if (rank == MASTER_NODE) cout << "Deleted CInterpolator container." << endl;
  }
  
  if (transfer_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
        if (transfer_container[iZone] != NULL)
          delete [] transfer_container[iZone];
    }
    delete [] transfer_container;
    if (rank == MASTER_NODE) cout << "Deleted CTransfer container." << endl;
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    if (geometry_container[iZone] != NULL) {
      for (unsigned short iMGlevel = 0; iMGlevel < config_container[iZone]->GetnMGLevels()+1; iMGlevel++) {
        if (geometry_container[iZone][iMGlevel] != NULL) delete geometry_container[iZone][iMGlevel];
      }
      delete [] geometry_container[iZone];
    }
  }
  delete [] geometry_container;
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    delete [] FFDBox[iZone];
  }
  delete [] FFDBox;
  if (rank == MASTER_NODE) cout << "Deleted CFreeFormDefBox class." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    delete surface_movement[iZone];
  }
  delete [] surface_movement;
  if (rank == MASTER_NODE) cout << "Deleted CSurfaceMovement class." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    delete grid_movement[iZone];
  }
  delete [] grid_movement;
  if (rank == MASTER_NODE) cout << "Deleted CVolumetricMovement class." << endl;

  if (config_container!= NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != NULL) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

  /*--- Deallocate output container ---*/
  if (output!= NULL) delete output;
  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

  if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl;


  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_CFD) ------------------------" << endl << endl;

}

void CDriver::Geometrical_Preprocessing() {

  unsigned short iMGlevel;
  unsigned short requestedMGlevels = config_container[ZONE_0]->GetnMGLevels();
  unsigned long iPoint;
  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Compute elements surrounding points, points surrounding points ---*/

    if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
    geometry_container[iZone][MESH_0]->SetPoint_Connectivity();

    /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/

    if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
    geometry_container[iZone][MESH_0]->SetRCM_Ordering(config_container[iZone]);

    /*--- recompute elements surrounding points, points surrounding points ---*/

    if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
    geometry_container[iZone][MESH_0]->SetPoint_Connectivity();

    /*--- Compute elements surrounding elements ---*/

    if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
    geometry_container[iZone][MESH_0]->SetElement_Connectivity();

    /*--- Check the orientation before computing geometrical quantities ---*/

    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry_container[iZone][MESH_0]->SetBoundVolume();
    geometry_container[iZone][MESH_0]->Check_IntElem_Orientation(config_container[iZone]);
    geometry_container[iZone][MESH_0]->Check_BoundElem_Orientation(config_container[iZone]);

    /*--- Create the edge structure ---*/

    if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
    geometry_container[iZone][MESH_0]->SetEdges();
    geometry_container[iZone][MESH_0]->SetVertex(config_container[iZone]);

    /*--- Compute cell center of gravity ---*/

    if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
    geometry_container[iZone][MESH_0]->SetCoord_CG();

    /*--- Create the control volume structures ---*/

    if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
    geometry_container[iZone][MESH_0]->SetControlVolume(config_container[iZone], ALLOCATE);
    geometry_container[iZone][MESH_0]->SetBoundControlVolume(config_container[iZone], ALLOCATE);

    /*--- Visualize a dual control volume if requested ---*/

    if ((config_container[iZone]->GetVisualize_CV() >= 0) &&
        (config_container[iZone]->GetVisualize_CV() < (long)geometry_container[iZone][MESH_0]->GetnPointDomain()))
      geometry_container[iZone][MESH_0]->VisualizeControlVolume(config_container[iZone], UPDATE);

    /*--- Identify closest normal neighbor ---*/

    if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
    geometry_container[iZone][MESH_0]->FindNormal_Neighbor(config_container[iZone]);

    /*--- Store the global to local mapping. ---*/

    if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
    geometry_container[iZone][MESH_0]->SetGlobal_to_Local_Point();

    /*--- Compute the surface curvature ---*/

    if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
    geometry_container[iZone][MESH_0]->ComputeSurf_Curvature(config_container[iZone]);

    /*--- Check for periodicity and disable MG if necessary. ---*/

    if (rank == MASTER_NODE) cout << "Checking for periodicity." << endl;
    geometry_container[iZone][MESH_0]->Check_Periodicity(config_container[iZone]);

    if ((config_container[iZone]->GetnMGLevels() != 0) && (rank == MASTER_NODE))
      cout << "Setting the multigrid structure." << endl;

  }

  /*--- Loop over all zones at each grid level. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Loop over all the new grid ---*/

    for (iMGlevel = 1; iMGlevel <= config_container[iZone]->GetnMGLevels(); iMGlevel++) {

      /*--- Create main agglomeration structure ---*/

      geometry_container[iZone][iMGlevel] = new CMultiGridGeometry(geometry_container, config_container, iMGlevel, iZone);

      /*--- Compute points surrounding points. ---*/

      geometry_container[iZone][iMGlevel]->SetPoint_Connectivity(geometry_container[iZone][iMGlevel-1]);

      /*--- Create the edge structure ---*/

      geometry_container[iZone][iMGlevel]->SetEdges();
      geometry_container[iZone][iMGlevel]->SetVertex(geometry_container[iZone][iMGlevel-1], config_container[iZone]);

      /*--- Create the control volume structures ---*/

      geometry_container[iZone][iMGlevel]->SetControlVolume(config_container[iZone], geometry_container[iZone][iMGlevel-1], ALLOCATE);
      geometry_container[iZone][iMGlevel]->SetBoundControlVolume(config_container[iZone], geometry_container[iZone][iMGlevel-1], ALLOCATE);
      geometry_container[iZone][iMGlevel]->SetCoord(geometry_container[iZone][iMGlevel-1]);

      /*--- Find closest neighbor to a surface point ---*/

      geometry_container[iZone][iMGlevel]->FindNormal_Neighbor(config_container[iZone]);

      /*--- Protect against the situation that we were not able to complete
       the agglomeration for this level, i.e., there weren't enough points.
       We need to check if we changed the total number of levels and delete
       the incomplete CMultiGridGeometry object. ---*/

      if (config_container[iZone]->GetnMGLevels() != requestedMGlevels) {
        delete geometry_container[iZone][iMGlevel];
        break;
      }

    }

  }

  /*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    if (config_container[iZone]->GetUnsteady_Simulation() && config_container[iZone]->GetGrid_Movement()) {
      for (iMGlevel = 0; iMGlevel <= config_container[iZone]->GetnMGLevels(); iMGlevel++) {
        for (iPoint = 0; iPoint < geometry_container[iZone][iMGlevel]->GetnPoint(); iPoint++) {

          /*--- Update cell volume ---*/

          geometry_container[iZone][iMGlevel]->node[iPoint]->SetVolume_n();
          geometry_container[iZone][iMGlevel]->node[iPoint]->SetVolume_nM1();

          /*--- Update point coordinates ---*/
          geometry_container[iZone][iMGlevel]->node[iPoint]->SetCoord_n();
          geometry_container[iZone][iMGlevel]->node[iPoint]->SetCoord_n1();

        }
      }
    }
  }

}

void CDriver::Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry,
                                   CConfig *config) {
  
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  poisson, wave, heat,
  fem, disc_adj_fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb;

  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  poisson          = false;  neg_spalart_allmaras = false;
  wave             = false;	 disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  heat             = false;
  transition       = false;
  template_solver  = false;
  
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;           break;
        
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(EXIT_FAILURE); break;
    }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    
    /*--- Allocate solution for a template problem ---*/
    
    if (template_solver) {
      solver_container[iMGlevel][TEMPLATE_SOL] = new CTemplateSolver(geometry[iMGlevel], config);
    }
    
    /*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/
    
    if (euler) {
      if (compressible) {
        solver_container[iMGlevel][FLOW_SOL] = new CEulerSolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
      if (incompressible) {
        solver_container[iMGlevel][FLOW_SOL] = new CIncEulerSolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
    }
    if (ns) {
      if (compressible) {
        solver_container[iMGlevel][FLOW_SOL] = new CNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        solver_container[iMGlevel][FLOW_SOL] = new CIncNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
    }
    if (turbulent) {
      if (spalart_allmaras) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel, solver_container[iMGlevel][FLOW_SOL]->GetFluidModel() );
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
      else if (neg_spalart_allmaras) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel, solver_container[iMGlevel][FLOW_SOL]->GetFluidModel() );
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
      else if (menter_sst) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSSTSolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
      if (transition) {
        solver_container[iMGlevel][TRANS_SOL] = new CTransLMSolver(geometry[iMGlevel], config, iMGlevel);
      }
    }
    if (poisson) {
      solver_container[iMGlevel][POISSON_SOL] = new CPoissonSolver(geometry[iMGlevel], config);
    }
    if (wave) {
      solver_container[iMGlevel][WAVE_SOL] = new CWaveSolver(geometry[iMGlevel], config);
    }
    if (heat) {
      solver_container[iMGlevel][HEAT_SOL] = new CHeatSolver(geometry[iMGlevel], config);
    }
    if (fem) {
      solver_container[iMGlevel][FEA_SOL] = new CFEM_ElasticitySolver(geometry[iMGlevel], config);
    }
    
    /*--- Allocate solution for adjoint problem ---*/
    
    if (adj_euler) {
      if (compressible) {
        solver_container[iMGlevel][ADJFLOW_SOL] = new CAdjEulerSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        solver_container[iMGlevel][ADJFLOW_SOL] = new CAdjIncEulerSolver(geometry[iMGlevel], config, iMGlevel);
      }
    }
    if (adj_ns) {
      if (compressible) {
        solver_container[iMGlevel][ADJFLOW_SOL] = new CAdjNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        solver_container[iMGlevel][ADJFLOW_SOL] = new CAdjIncNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
    }
    if (adj_turb) {
      solver_container[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolver(geometry[iMGlevel], config, iMGlevel);
    }
    
    if (disc_adj) {
      solver_container[iMGlevel][ADJFLOW_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver_container[iMGlevel][FLOW_SOL], RUNTIME_FLOW_SYS, iMGlevel);
      if (disc_adj_turb)
        solver_container[iMGlevel][ADJTURB_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver_container[iMGlevel][TURB_SOL], RUNTIME_TURB_SYS, iMGlevel);
    }

    if (disc_adj_fem) {
      solver_container[iMGlevel][ADJFEA_SOL] = new CDiscAdjFEASolver(geometry[iMGlevel], config, solver_container[iMGlevel][FEA_SOL], RUNTIME_FEA_SYS, iMGlevel);
    }
  }
  

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  bool update_geo = true;
  if (config->GetFSI_Simulation()) update_geo = false;

  Solver_Restart(solver_container, geometry, config, update_geo);

}

void CDriver::Solver_Restart(CSolver ***solver_container, CGeometry **geometry,
                             CConfig *config, bool update_geo) {

  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  poisson, wave, heat,
  fem,
  template_solver, disc_adj, disc_adj_fem, disc_adj_turb;
  int val_iter = 0;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  poisson          = false;
  wave             = false;  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  heat             = false;  disc_adj_turb    = false;
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
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
  }


  /*--- Load restarts for any of the active solver containers. Note that
   these restart routines fill the fine grid and interpolate to all MG levels. ---*/

  if (restart || restart_flow) {
    if (euler || ns) {
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, val_iter, update_geo);
    }
    if (turbulent) {
      solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, val_iter, update_geo);
    }
    if (fem) {
      if (dynamic) val_iter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
      solver_container[MESH_0][FEA_SOL]->LoadRestart(geometry, solver_container, config, val_iter, update_geo);
    }
  }

  if (restart) {
    if (template_solver) {
      no_restart = true;
    }
    if (poisson) {
      no_restart = true;
    }
    if (wave) {
      no_restart = true;
    }
    if (heat) {
      no_restart = true;
    }
    if (adj_euler || adj_ns) {
      solver_container[MESH_0][ADJFLOW_SOL]->LoadRestart(geometry, solver_container, config, val_iter, update_geo);
    }
    if (adj_turb) {
      no_restart = true;
    }
    if (disc_adj) {
      solver_container[MESH_0][ADJFLOW_SOL]->LoadRestart(geometry, solver_container, config, val_iter, update_geo);
      if (disc_adj_turb)
        solver_container[MESH_0][ADJTURB_SOL]->LoadRestart(geometry, solver_container, config, val_iter, update_geo);
    }
    if (disc_adj_fem) {
        if (dynamic) val_iter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
        solver_container[MESH_0][ADJFEA_SOL]->LoadRestart(geometry, solver_container, config, val_iter, update_geo);
    }
  }

  /*--- Exit if a restart was requested for a solver that is not available. ---*/

  if (no_restart) {
    if (rank == MASTER_NODE) {
      cout << endl << "A restart capability has not been implemented yet for this solver. " << endl;
      cout << "Please set RESTART_SOL= NO and try again." << endl << endl;
    }
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Think about calls to pre / post-processing here, plus realizability checks. ---*/
  

}

void CDriver::Solver_Postprocessing(CSolver ***solver_container, CGeometry **geometry,
                                    CConfig *config) {
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  poisson, wave, heat, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb, disc_adj_fem;
  
  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  poisson          = false;  neg_spalart_allmaras = false;
  wave             = false;  disc_adj        = false;
  fem              = false;  disc_adj_fem    = false;
  heat             = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;           break;
    }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    
    /*--- DeAllocate solution for a template problem ---*/
    
    if (template_solver) {
      delete solver_container[iMGlevel][TEMPLATE_SOL];
    }
    
    /*--- DeAllocate solution for adjoint problem ---*/
    
    if (adj_euler || adj_ns || disc_adj) {
      delete solver_container[iMGlevel][ADJFLOW_SOL];
      if (disc_adj_turb || adj_turb) {
        delete solver_container[iMGlevel][ADJTURB_SOL];
      }
    }
    
    /*--- DeAllocate solution for direct problem ---*/
    
    if (euler || ns) {
      delete solver_container[iMGlevel][FLOW_SOL];
    }
    
    if (turbulent) {
      if (spalart_allmaras || neg_spalart_allmaras || menter_sst ) {
        delete solver_container[iMGlevel][TURB_SOL];
      }
      if (transition) {
        delete solver_container[iMGlevel][TRANS_SOL];
      }
    }
    if (poisson) {
      delete solver_container[iMGlevel][POISSON_SOL];
    }
    if (wave) {
      delete solver_container[iMGlevel][WAVE_SOL];
    }
    if (heat) {
      delete solver_container[iMGlevel][HEAT_SOL];
    }
    if (fem) {
      delete solver_container[iMGlevel][FEA_SOL];
    }
    if (disc_adj_fem) {
      delete solver_container[iMGlevel][ADJFEA_SOL];
    }
    
    delete [] solver_container[iMGlevel];
  }
  
}

void CDriver::Integration_Preprocessing(CIntegration **integration_container,
    CGeometry **geometry, CConfig *config) {

  bool euler, adj_euler, ns, adj_ns, turbulent, adj_turb, poisson, wave, fem,
      heat, template_solver, transition, disc_adj, disc_adj_fem;
  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false;
  ns               = false; adj_ns           = false;
  turbulent        = false; adj_turb         = false;
  poisson          = false; disc_adj         = false;
  wave             = false;
  heat             = false;
  fem 			       = false; disc_adj_fem     = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER : euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS : ns = true; turbulent = true; disc_adj = true; break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
  }

  /*--- Allocate solution for a template problem ---*/
  if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);

  /*--- Allocate solution for direct problem ---*/
  if (euler) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (ns) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (turbulent) integration_container[TURB_SOL] = new CSingleGridIntegration(config);
  if (transition) integration_container[TRANS_SOL] = new CSingleGridIntegration(config);
  if (poisson) integration_container[POISSON_SOL] = new CSingleGridIntegration(config);
  if (wave) integration_container[WAVE_SOL] = new CSingleGridIntegration(config);
  if (heat) integration_container[HEAT_SOL] = new CSingleGridIntegration(config);
  if (fem) integration_container[FEA_SOL] = new CStructuralIntegration(config);

  /*--- Allocate solution for adjoint problem ---*/
  if (adj_euler) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_ns) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_turb) integration_container[ADJTURB_SOL] = new CSingleGridIntegration(config);

  if (disc_adj) integration_container[ADJFLOW_SOL] = new CIntegration(config);
  if (disc_adj_fem) integration_container[ADJFEA_SOL] = new CIntegration(config);

}

void CDriver::Integration_Postprocessing(CIntegration **integration_container,
    CGeometry **geometry, CConfig *config) {
  bool euler, adj_euler, ns, adj_ns, turbulent, adj_turb, poisson, wave, fem,
      heat, template_solver, transition, disc_adj, disc_adj_fem;

  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false;
  ns               = false; adj_ns           = false;
  turbulent        = false; adj_turb         = false;
  poisson          = false; disc_adj         = false;
  wave             = false;
  heat             = false;
  fem              = false; disc_adj_fem     = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER : euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS : ns = true; turbulent = true; disc_adj = true; break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
  }

  /*--- DeAllocate solution for a template problem ---*/
  if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);

  /*--- DeAllocate solution for direct problem ---*/
  if (euler || ns) delete integration_container[FLOW_SOL];
  if (turbulent) delete integration_container[TURB_SOL];
  if (transition) delete integration_container[TRANS_SOL];
  if (poisson) delete integration_container[POISSON_SOL];
  if (wave) delete integration_container[WAVE_SOL];
  if (heat) delete integration_container[HEAT_SOL];
  if (fem) delete integration_container[FEA_SOL];
  if (disc_adj_fem) delete integration_container[ADJFEA_SOL];

  /*--- DeAllocate solution for adjoint problem ---*/
  if (adj_euler || adj_ns || disc_adj) delete integration_container[ADJFLOW_SOL];
  if (adj_turb) delete integration_container[ADJTURB_SOL];
  

}

void CDriver::Numerics_Preprocessing(CNumerics ****numerics_container,
                                     CSolver ***solver_container, CGeometry **geometry,
                                     CConfig *config) {
  
  unsigned short iMGlevel, iSol, nDim,
  
  nVar_Template         = 0,
  nVar_Flow             = 0,
  nVar_Trans            = 0,
  nVar_Turb             = 0,
  nVar_Adj_Flow         = 0,
  nVar_Adj_Turb         = 0,
  nVar_Poisson          = 0,
  nVar_FEM              = 0,
  nVar_Wave             = 0,
  nVar_Heat             = 0;
  
  su2double *constants = NULL;
  
  bool
  euler, adj_euler,
  ns, adj_ns,
  turbulent, adj_turb,
  spalart_allmaras, neg_spalart_allmaras, menter_sst,
  poisson,
  wave,
  fem,
  heat,
  transition,
  template_solver;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;   ns               = false;   turbulent        = false;
  poisson          = false;
  adj_euler        = false;   adj_ns           = false;   adj_turb         = false;
  wave             = false;   heat             = false;
  fem              = false;
  spalart_allmaras = false; neg_spalart_allmaras = false;	menter_sst       = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case DISC_ADJ_EULER: euler = true; break;
    case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES: ns = true; break;
    case RANS : case DISC_ADJ_RANS:  ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: case DISC_ADJ_FEM: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true; constants = solver_container[MESH_0][TURB_SOL]->GetConstants(); break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(EXIT_FAILURE); break;
    }
  
  /*--- Number of variables for the template ---*/
  
  if (template_solver) nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  
  /*--- Number of variables for direct problem ---*/
  
  if (euler)        nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  if (ns)           nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  if (turbulent)    nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
  if (transition)   nVar_Trans = solver_container[MESH_0][TRANS_SOL]->GetnVar();
  if (poisson)      nVar_Poisson = solver_container[MESH_0][POISSON_SOL]->GetnVar();
  
  if (wave)         nVar_Wave = solver_container[MESH_0][WAVE_SOL]->GetnVar();
  if (fem)          nVar_FEM = solver_container[MESH_0][FEA_SOL]->GetnVar();
  if (heat)         nVar_Heat = solver_container[MESH_0][HEAT_SOL]->GetnVar();
  
  /*--- Number of variables for adjoint problem ---*/
  
  if (adj_euler)        nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_ns)           nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_turb)         nVar_Adj_Turb = solver_container[MESH_0][ADJTURB_SOL]->GetnVar();
  
  /*--- Number of dimensions ---*/
  
  nDim = geometry[MESH_0]->GetnDim();
  
  /*--- Definition of the Class for the numerical method: numerics_container[MESH_LEVEL][EQUATION][EQ_TERM] ---*/
  if (fem){
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel] = new CNumerics** [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        numerics_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS_FEA];
    }
  }
  else{
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel] = new CNumerics** [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        numerics_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
    }
  }
  
  /*--- Solver definition for the template problem ---*/
  if (template_solver) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][TEMPLATE_SOL][CONV_TERM] = new CConvective_Template(nDim, nVar_Template, config);
        break;
      default : cout << "Convective scheme not implemented (template_solver)." << endl; exit(EXIT_FAILURE); break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics_container[iMGlevel][TEMPLATE_SOL][VISC_TERM] = new CViscous_Template(nDim, nVar_Template, config);
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics_container[iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM] = new CSource_Template(nDim, nVar_Template, config);
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM] = new CConvective_Template(nDim, nVar_Template, config);
    }
    
  }
  
  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(EXIT_FAILURE);
        break;
        
      case SPACE_CENTERED :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim, nVar_Flow, config); break;
            case JST_KE : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_KE_Flow(nDim, nVar_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJSTArtComp_Flow(nDim, nVar_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
          
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
                  numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                  numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                }
              } else {
                
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                  numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                }
              }
              break;
              
            case AUSM:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case TURKEL:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case HLLC:
              if (ideal_gas) {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                  numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              else {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                  numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              break;
              
            case MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
        }
        break;
        
      default :
        cout << "Convective scheme not implemented (euler and ns)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible) {
      if (ideal_gas) {
        
        /*--- Compressible flow Ideal gas ---*/
        numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrected_Flow(nDim, nVar_Flow, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        
      } else {
        
        /*--- Compressible flow Realgas ---*/
        numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CGeneralAvgGradCorrected_Flow(nDim, nVar_Flow, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, config);
        
      }
    }
    if (incompressible) {
      /*--- Incompressible flow, use artificial compressibility method ---*/
      numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
      
      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      
      if (config->GetBody_Force() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceBodyForce(nDim, nVar_Flow, config);
      else if (config->GetRotating_Frame() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_Flow(nDim, nVar_Flow, config);
      else if (config->GetAxisymmetric() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_Flow(nDim, nVar_Flow, config);
      else if (config->GetGravityForce() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceGravity(nDim, nVar_Flow, config);
      else if (config->GetWind_Gust() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceWindGust(nDim, nVar_Flow, config);
      else
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
      
      numerics_container[iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
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
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
          else if (neg_spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        }
        break;
      default :
        cout << "Convective scheme not implemented (turbulent)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, true, config);
      else if (neg_spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSA_Neg(nDim, nVar_Turb, true, config);
      else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, true, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
      else if (neg_spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA_Neg(nDim, nVar_Turb, config);
      else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSST(nDim, nVar_Turb, constants, config);
      numerics_container[iMGlevel][TURB_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Turb, config);
    }
    
    /*--- Definition of the boundary condition method ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, false, config);
      }
      else if (neg_spalart_allmaras) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA_Neg(nDim, nVar_Turb, false, config);
      }
      else if (menter_sst) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, false, config);
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
          numerics_container[iMGlevel][TRANS_SOL][CONV_TERM] = new CUpwSca_TransLM(nDim, nVar_Trans, config);
        }
        break;
      default :
        cout << "Convective scheme not implemented (transition)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TRANS_SOL][VISC_TERM] = new CAvgGradCorrected_TransLM(nDim, nVar_Trans, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TRANS_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TransLM(nDim, nVar_Trans, config);
      numerics_container[iMGlevel][TRANS_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Trans, config);
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TRANS_SOL][CONV_BOUND_TERM] = new CUpwLin_TransLM(nDim, nVar_Trans, config);
    }
  }
  
  /*--- Solver definition for the poisson potential problem ---*/
  if (poisson) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][POISSON_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Poisson, config);
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][POISSON_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Poisson, config);
    numerics_container[MESH_0][POISSON_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Poisson, config);
    
  }
  
  /*--- Solver definition for the poisson potential problem ---*/
  if (heat) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][HEAT_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Heat, config);
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][HEAT_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Heat, config);
    numerics_container[MESH_0][HEAT_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Heat, config);
    
  }
  
  /*--- Solver definition for the flow adjoint problem ---*/
  
  if (adj_euler || adj_ns) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    
    switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(EXIT_FAILURE);
        break;
        
      case SPACE_CENTERED :
        
        if (compressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        
        if (incompressible) {
          
          /*--- Incompressible flow, use artificial compressibility method ---*/
          
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJSTArtComp_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        
        break;
        
      case SPACE_UPWIND :
        
        if (compressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
        }
        
        if (incompressible) {
          
          /*--- Incompressible flow, use artificial compressibility method ---*/
          
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
        }
        
        break;
        
      default :
        cout << "Convective scheme not implemented (adj_euler and adj_ns)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    
    if (compressible) {
      
      /*--- Compressible flow ---*/
      
      numerics_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
      numerics_container[MESH_0][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
      
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
        numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
      }
      
    }
    
    if (incompressible) {
      
      /*--- Incompressible flow, use artificial compressibility method ---*/
      
      numerics_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
      numerics_container[MESH_0][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
      
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
        numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
      }
      
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      
      /*--- Note that RANS is incompatible with Axisymmetric or Rotational (Fix it!) ---*/
      
      if (compressible) {
        
        if (adj_ns) {
          
          numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceViscous_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          if (config->GetRotating_Frame() == YES)
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        
        else {
          
          if (config->GetRotating_Frame() == YES)
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else if (config->GetAxisymmetric() == YES)
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
          
          numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
          
        }
        
      }
      
      if (incompressible) {
        
        numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
        numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
        
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
            numerics_container[iMGlevel][ADJTURB_SOL][CONV_TERM] = new CUpwSca_AdjTurb(nDim, nVar_Adj_Turb, config);
          }
          else if (neg_spalart_allmaras) {cout << "Adjoint Neg SA turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
          else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
        break;
      default :
        cout << "Convective scheme not implemented (adj_turb)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        numerics_container[iMGlevel][ADJTURB_SOL][VISC_TERM] = new CAvgGradCorrected_AdjTurb(nDim, nVar_Adj_Turb, config);
      }
      else if (neg_spalart_allmaras) {cout << "Adjoint Neg SA turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
      else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_AdjTurb(nDim, nVar_Adj_Turb, config);
        numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjTurb(nDim, nVar_Adj_Turb, config);
      }
      else if (neg_spalart_allmaras) {cout << "Adjoint Neg SA turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
      else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[iMGlevel][ADJTURB_SOL][CONV_BOUND_TERM] = new CUpwLin_AdjTurb(nDim, nVar_Adj_Turb, config);
      else if (neg_spalart_allmaras) {cout << "Adjoint Neg SA turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
      else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
    }
    
  }

  /*--- Solver definition for the wave problem ---*/
  if (wave) {

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][WAVE_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Wave, config);

  }

  /*--- Solver definition for the FEM problem ---*/
  if (fem) {

  /*--- Initialize the container for FEA_TERM. This will be the only one for most of the cases ---*/
  switch (config->GetGeometricConditions()) {
      case SMALL_DEFORMATIONS :
        switch (config->GetMaterialModel()) {
          case LINEAR_ELASTIC: numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_LinearElasticity(nDim, nVar_FEM, config); break;
          case NEO_HOOKEAN : cout << "Material model does not correspond to geometric conditions." << endl; exit(EXIT_FAILURE); break;
          default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
        }
        break;
      case LARGE_DEFORMATIONS :
        switch (config->GetMaterialModel()) {
        case LINEAR_ELASTIC: cout << "Material model does not correspond to geometric conditions." << endl; exit(EXIT_FAILURE); break;
          case NEO_HOOKEAN :
            switch (config->GetMaterialCompressibility()) {
              case COMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config); break;
              case INCOMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_NeoHookean_Incomp(nDim, nVar_FEM, config); break;
              default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
            }
            break;
          case KNOWLES:
            switch (config->GetMaterialCompressibility()) {
              case NEARLY_INCOMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config); break;
              case INCOMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config); break;
              default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
            }
            break;
          case IDEAL_DE:
            switch (config->GetMaterialCompressibility()) {
              case NEARLY_INCOMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_IdealDE(nDim, nVar_FEM, config); break;
              default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
            }
            break;
          default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
        }
        break;
      default: cout << " Solver not implemented." << endl; exit(EXIT_FAILURE); break;
    }

  /*--- The following definitions only make sense if we have a non-linear solution ---*/
  if (config->GetGeometricConditions() == LARGE_DEFORMATIONS){

      /*--- This allocates a container for electromechanical effects ---*/

      bool de_effects = config->GetDE_Effects();
      if (de_effects) numerics_container[MESH_0][FEA_SOL][DE_TERM] = new CFEM_DielectricElastomer(nDim, nVar_FEM, config);

      string filename;
      ifstream properties_file;

      filename = config->GetFEA_FileName();
      if (nZone > 1)
        filename = config->GetMultizone_FileName(filename, iZone);

      properties_file.open(filename.data(), ios::in);

      /*--- In case there is a properties file, containers are allocated for a number of material models ---*/

      if (!(properties_file.fail())) {

          numerics_container[MESH_0][FEA_SOL][MAT_NHCOMP]  = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config);
          numerics_container[MESH_0][FEA_SOL][MAT_NHINC]   = new CFEM_NeoHookean_Incomp(nDim, nVar_FEM, config);
          numerics_container[MESH_0][FEA_SOL][MAT_IDEALDE] = new CFEM_IdealDE(nDim, nVar_FEM, config);
          numerics_container[MESH_0][FEA_SOL][MAT_KNOWLES] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config);

          properties_file.close();
      }
  }

  }

}

void CDriver::Numerics_Postprocessing(CNumerics ****numerics_container,
                                      CSolver ***solver_container, CGeometry **geometry,
                                      CConfig *config) {
  
  unsigned short iMGlevel, iSol;
  
  
  bool
  euler, adj_euler,
  ns, adj_ns,
  turbulent, adj_turb,
  spalart_allmaras, neg_spalart_allmaras, menter_sst,
  poisson,
  wave,
  fem,
  heat,
  transition,
  template_solver;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;   ns               = false;   turbulent        = false;
  poisson          = false;
  adj_euler        = false;   adj_ns           = false;   adj_turb         = false;
  wave             = false;   heat             = false;   fem        = false;
  spalart_allmaras = false; neg_spalart_allmaras = false; menter_sst       = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case DISC_ADJ_EULER: euler = true; break;
    case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES: ns = true; break;
    case RANS : case DISC_ADJ_RANS:  ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: case DISC_ADJ_FEM: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;  break;
        
    }
  
  /*--- Solver definition for the template problem ---*/
  if (template_solver) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          delete numerics_container[iMGlevel][TEMPLATE_SOL][CONV_TERM];
        break;
    }
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      /*--- Definition of the viscous scheme for each equation and mesh level ---*/
      delete numerics_container[iMGlevel][TEMPLATE_SOL][VISC_TERM];
      /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
      delete numerics_container[iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM];
      /*--- Definition of the boundary condition method ---*/
      delete numerics_container[iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM];
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
            case LAX : case JST :  case JST_KE : delete numerics_container[MESH_0][FLOW_SOL][CONV_TERM]; break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][FLOW_SOL][CONV_TERM];
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Centered_Flow()) {
              
            case LAX : case JST : delete numerics_container[MESH_0][FLOW_SOL][CONV_TERM]; break;
              
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][FLOW_SOL][CONV_TERM];
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
          
        }
        break;
      case SPACE_UPWIND :
        
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE: case AUSM : case TURKEL: case HLLC: case MSW:  case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[iMGlevel][FLOW_SOL][CONV_TERM];
                delete numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
              }
              
              break;
          }
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[iMGlevel][FLOW_SOL][CONV_TERM];
                delete numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
              }
              break;
          }
        }
        
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible||incompressible) {
      /*--- Compressible flow Ideal gas ---*/
      delete numerics_container[MESH_0][FLOW_SOL][VISC_TERM];
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        delete numerics_container[iMGlevel][FLOW_SOL][VISC_TERM];
      
      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        delete numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM];
      
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      delete numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM];
      delete numerics_container[iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM];
    }
    
  }
  
  
  /*--- Solver definition for the turbulent model problem ---*/
  
  if (turbulent) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          if (spalart_allmaras || neg_spalart_allmaras ||menter_sst)
            delete numerics_container[iMGlevel][TURB_SOL][CONV_TERM];
        }
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (spalart_allmaras || neg_spalart_allmaras || menter_sst) {
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        delete numerics_container[iMGlevel][TURB_SOL][VISC_TERM];
        delete numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM];
        delete numerics_container[iMGlevel][TURB_SOL][SOURCE_SECOND_TERM];
        /*--- Definition of the boundary condition method ---*/
        delete numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM];
        delete numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM];
        
      }
    }
    
  }
  
  /*--- Solver definition for the transition model problem ---*/
  if (transition) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          delete numerics_container[iMGlevel][TRANS_SOL][CONV_TERM];
        }
        break;
    }
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      /*--- Definition of the viscous scheme for each equation and mesh level ---*/
      delete numerics_container[iMGlevel][TRANS_SOL][VISC_TERM];
      /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
      delete numerics_container[iMGlevel][TRANS_SOL][SOURCE_FIRST_TERM];
      delete numerics_container[iMGlevel][TRANS_SOL][SOURCE_SECOND_TERM];
      /*--- Definition of the boundary condition method ---*/
      delete numerics_container[iMGlevel][TRANS_SOL][CONV_BOUND_TERM];
    }
  }
  
  /*--- Solver definition for the poisson potential problem ---*/
  if (poisson || heat) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    delete numerics_container[MESH_0][POISSON_SOL][VISC_TERM];
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    delete numerics_container[MESH_0][POISSON_SOL][SOURCE_FIRST_TERM];
    delete numerics_container[MESH_0][POISSON_SOL][SOURCE_SECOND_TERM];
    
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
              delete numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM];
              break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM];
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];
          
        }
        
        if (incompressible) {
          
          /*--- Incompressible flow, use artificial compressibility method ---*/
          
          switch (config->GetKind_Centered_AdjFlow()) {
            case LAX : case JST:
              delete numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM]; break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM];
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];
          
        }
        
        break;
        
      case SPACE_UPWIND :
        
        if (compressible || incompressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Upwind_AdjFlow()) {
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM];
                delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];
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
        delete numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM];
        delete numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM];
      }
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      
      
      if (compressible || incompressible) {
        
        delete numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM];
        delete numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM];
        
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
            delete numerics_container[iMGlevel][ADJTURB_SOL][CONV_TERM];
          }
        break;
    }
    
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        /*--- Definition of the viscous scheme for each equation and mesh level ---*/
        delete numerics_container[iMGlevel][ADJTURB_SOL][VISC_TERM];
        /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
        delete numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_FIRST_TERM];
        delete numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_SECOND_TERM];
        /*--- Definition of the boundary condition method ---*/
        delete numerics_container[iMGlevel][ADJTURB_SOL][CONV_BOUND_TERM];
      }
    }
  }
  
  /*--- Solver definition for the wave problem ---*/
  if (wave) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    delete numerics_container[MESH_0][WAVE_SOL][VISC_TERM];
    
  }
  
  /*--- Solver definition for the FEA problem ---*/
  if (fem) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    delete numerics_container[MESH_0][FEA_SOL][FEA_TERM];
    
  }
  
  /*--- Definition of the Class for the numerical method: numerics_container[MESH_LEVEL][EQUATION][EQ_TERM] ---*/
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    for (iSol = 0; iSol < MAX_SOLS; iSol++) {
      delete [] numerics_container[iMGlevel][iSol];
    }
    delete[] numerics_container[iMGlevel];
  }
  
}

void CDriver::Iteration_Preprocessing() {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Initial print to console for this zone. ---*/

  if (rank == MASTER_NODE) cout << "Zone " << iZone+1;

  /*--- Loop over all zones and instantiate the physics iteration. ---*/

  switch (config_container[iZone]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:

      if(config_container[iZone]->GetBoolTurbomachinery()){
        if (rank == MASTER_NODE)
          cout << ": Euler/Navier-Stokes/RANS turbomachinery fluid iteration." << endl;
      iteration_container[iZone] = new CTurboIteration(config_container[iZone]);

      }
      else{
        if (rank == MASTER_NODE)
          cout << ": Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration_container[iZone] = new CFluidIteration(config_container[iZone]);
      }
      break;

    case WAVE_EQUATION:
      if (rank == MASTER_NODE)
        cout << ": wave iteration." << endl;
      iteration_container[iZone] = new CWaveIteration(config_container[iZone]);
      break;

    case HEAT_EQUATION:
      if (rank == MASTER_NODE)
        cout << ": heat iteration." << endl;
      iteration_container[iZone] = new CHeatIteration(config_container[iZone]);
      break;

    case POISSON_EQUATION:
      if (rank == MASTER_NODE)
        cout << ": poisson iteration." << endl;
      iteration_container[iZone] = new CPoissonIteration(config_container[iZone]);
      break;

    case FEM_ELASTICITY:
      if (rank == MASTER_NODE)
        cout << ": FEM iteration." << endl;
      iteration_container[iZone] = new CFEM_StructuralAnalysis(config_container[iZone]);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << ": adjoint Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration_container[iZone] = new CAdjFluidIteration(config_container[iZone]);
      break;

    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << ": discrete adjoint Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration_container[iZone] = new CDiscAdjFluidIteration(config_container[iZone]);
      break;

    case DISC_ADJ_FEM:
      if (rank == MASTER_NODE)
        cout << ": discrete adjoint FEM structural iteration." << endl;
      iteration_container[iZone] = new CDiscAdjFEAIteration(config_container[iZone]);
      break;
  }
  
}

void CDriver::Interface_Preprocessing() {

  int rank = MASTER_NODE;
  unsigned short donorZone, targetZone;
  unsigned short nVar, nVarTransfer;

  unsigned short nMarkerTarget, iMarkerTarget, nMarkerDonor, iMarkerDonor;

  /*--- Initialize some useful booleans ---*/
  bool fluid_donor, structural_donor;
  bool fluid_target, structural_target;

  bool discrete_adjoint = config_container[ZONE_0]->GetDiscrete_Adjoint();

  int markDonor, markTarget, Donor_check, Target_check, iMarkerInt, nMarkerInt;

#ifdef HAVE_MPI
  int *Buffer_Recv_mark = NULL, iRank, nProcessor = 1;;

  MPI_Comm_rank(config_container[ZONE_0]->GetMPICommunicator(), &rank);
  MPI_Comm_size(config_container[ZONE_0]->GetMPICommunicator(), &nProcessor);

  if (rank == MASTER_NODE)
  Buffer_Recv_mark = new int[nProcessor];
#endif

  if (config_container[ZONE_0]->GetFSI_Simulation() && nZone != 2 && rank == MASTER_NODE) {
    cout << "Error, cannot run the FSI solver on more than 2 zones!" << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Coupling between zones ---*/
  // There's a limit here, the interface boundary must connect only 2 zones

  /*--- Loops over all target and donor zones to find which ones are connected through an interface boundary (fsi or sliding mesh) ---*/
  for (targetZone = 0; targetZone < nZone; targetZone++) {

    for (donorZone = 0; donorZone < nZone; donorZone++) {

      if ( donorZone == targetZone ) // We're processing the same zone, so skip the following
        continue;

      nMarkerInt = (int) ( config_container[donorZone]->GetMarker_n_ZoneInterface() / 2 );

      /*--- Loops on Interface markers to find if the 2 zones are sharing the boundary and to determine donor and target marker tag ---*/
      for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

        markDonor  = -1;
        markTarget = -1;

        /*--- On the donor side ---*/
        nMarkerDonor = config_container[donorZone]->GetnMarker_All();

        for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++) {

          /*--- If the tag GetMarker_All_ZoneInterface(iMarker) equals the index we are looping at ---*/
          if ( config_container[donorZone]->GetMarker_All_ZoneInterface(iMarkerDonor) == iMarkerInt ) {
            /*--- We have identified the identifier for the interface marker ---*/
            markDonor = iMarkerDonor;

            break;
          }
        }

        /*--- On the target side ---*/
        nMarkerTarget = config_container[targetZone]->GetnMarker_All();

      for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++) {

          /*--- If the tag GetMarker_All_ZoneInterface(iMarker) equals the index we are looping at ---*/
        if ( config_container[targetZone]->GetMarker_All_ZoneInterface(iMarkerTarget) == iMarkerInt ) {
            /*--- We have identified the identifier for the interface marker ---*/
            markTarget = iMarkerTarget;

            break;
        } 
        }

#ifdef HAVE_MPI

      Donor_check  = -1;
      Target_check = -1;

        /*--- We gather a vector in MASTER_NODE that determines if the boundary is not on the processor because of the partition or because the zone does not include it ---*/

        SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

      if (rank == MASTER_NODE) {
        for (iRank = 0; iRank < nProcessor; iRank++) {
          if( Buffer_Recv_mark[iRank] != -1 ) {
              Donor_check = Buffer_Recv_mark[iRank];

              break;
            }
          }
        }

        SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

        SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

      if (rank == MASTER_NODE){
        for (iRank = 0; iRank < nProcessor; iRank++){
          if( Buffer_Recv_mark[iRank] != -1 ){
              Target_check = Buffer_Recv_mark[iRank];

              break;
            }
          }
        }

        SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

#else
      Donor_check  = markDonor;
      Target_check = markTarget;  
#endif

      /* --- Check ifzones are actually sharing the interface boundary, if not skip ---*/        
      if(Target_check == -1 || Donor_check == -1)
          continue;

        /*--- Set some boolean to properly allocate data structure later ---*/
      fluid_target      = false; 
      structural_target = false;

      fluid_donor       = false; 
      structural_donor  = false;

      switch ( config_container[targetZone]->GetKind_Solver() ) {

        case EULER : case NAVIER_STOKES: case RANS: 
        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          fluid_target  = true;   

          break;

        case FEM_ELASTICITY: case DISC_ADJ_FEM:
          structural_target = true;   

          break;
        }


      switch ( config_container[donorZone]->GetKind_Solver() ) {

      case EULER : case NAVIER_STOKES: case RANS: 
      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
        fluid_donor  = true;   

          break;

      case FEM_ELASTICITY: case DISC_ADJ_FEM:
        structural_donor = true;  

          break;
        }

        /*--- Begin the creation of the communication pattern among zones ---*/

        /*--- Retrieve the number of conservative variables (for problems not involving structural analysis ---*/
        if (!structural_donor && !structural_target)
          nVar = solver_container[donorZone][MESH_0][FLOW_SOL]->GetnVar();
        else
          /*--- If at least one of the components is structural ---*/
          nVar = nDim;

      if (rank == MASTER_NODE) cout << "From zone " << donorZone << " to zone " << targetZone << ": ";

        /*--- Match Zones ---*/
      if (rank == MASTER_NODE) cout << "Setting coupling "<<endl;

        /*--- If the mesh is matching: match points ---*/
      if ( config_container[donorZone]->GetMatchingMesh() ) {
        if (rank == MASTER_NODE) 
            cout << "between matching meshes. " << endl;
        geometry_container[donorZone][MESH_0]->MatchZone(config_container[donorZone], geometry_container[targetZone][MESH_0], config_container[targetZone], donorZone, nZone);
        }
        /*--- Else: interpolate ---*/
        else {
        switch (config_container[donorZone]->GetKindInterpolation()) {

          case NEAREST_NEIGHBOR:
            interpolator_container[donorZone][targetZone] = new CNearestNeighbor(geometry_container, config_container, donorZone, targetZone);
            if (rank == MASTER_NODE) cout << "using a nearest-neighbor approach." << endl;

            break;

          case ISOPARAMETRIC:
            interpolator_container[donorZone][targetZone] = new CIsoparametric(geometry_container, config_container, donorZone, targetZone);
            if (rank == MASTER_NODE) cout << "using an isoparametric approach." << endl;

            break;

          case WEIGHTED_AVERAGE:
            interpolator_container[donorZone][targetZone] = new CSlidingMesh(geometry_container, config_container, donorZone, targetZone);
            if (rank == MASTER_NODE) cout << "using an sliding mesh approach." << endl;

            break;

          case CONSISTCONSERVE:
            if ( targetZone > 0 && structural_target ) {
              interpolator_container[donorZone][targetZone] = new CMirror(geometry_container, config_container, donorZone, targetZone);
              if (rank == MASTER_NODE) cout << "using a mirror approach: matching coefficients from opposite mesh." << endl;
            }
            else{
              interpolator_container[donorZone][targetZone] = new CIsoparametric(geometry_container, config_container, donorZone, targetZone);
              if (rank == MASTER_NODE) cout << "using an isoparametric approach." << endl;
            }
            if ( targetZone == 0 && structural_target ) {
              if (rank == MASTER_NODE) cout << "Consistent and conservative interpolation assumes the structure model mesh is evaluated second. Somehow this has not happened. The isoparametric coefficients will be calculated for both meshes, and are not guaranteed to be consistent." << endl;
            }


            break;

          }
        }

        /*--- Initialize the appropriate transfer strategy ---*/
      if (rank == MASTER_NODE) cout << "Transferring ";

        if (fluid_donor && structural_target && (!discrete_adjoint)) {
          nVarTransfer = 2;
        transfer_container[donorZone][targetZone] = new CTransfer_FlowTraction(nVar, nVarTransfer, config_container[donorZone]);
        if (rank == MASTER_NODE) cout << "flow tractions. "<< endl;
      }
      else if (structural_donor && fluid_target && (!discrete_adjoint)) {
        nVarTransfer = 0;
        transfer_container[donorZone][targetZone] = new CTransfer_StructuralDisplacements(nVar, nVarTransfer, config_container[donorZone]);
        if (rank == MASTER_NODE) cout << "structural displacements. "<< endl;
      }
      else if (fluid_donor && structural_target && discrete_adjoint) {
        nVarTransfer = 2;
        transfer_container[donorZone][targetZone] = new CTransfer_FlowTraction_DiscAdj(nVar, nVarTransfer, config_container[donorZone]);
        if (rank == MASTER_NODE) cout << "flow tractions. "<< endl;
      }
      else if (structural_donor && fluid_target && discrete_adjoint){
        nVarTransfer = 0;
        transfer_container[donorZone][targetZone] = new CTransfer_StructuralDisplacements_DiscAdj(nVar, nVarTransfer, config_container[donorZone]);
        if (rank == MASTER_NODE) cout << "structural displacements. "<< endl;
      }
      else if (!structural_donor && !structural_target) {
        nVarTransfer = 0;
        nVar = solver_container[donorZone][MESH_0][FLOW_SOL]->GetnPrimVar();
        transfer_container[donorZone][targetZone] = new CTransfer_SlidingInterface(nVar, nVarTransfer, config_container[donorZone]);
        if (rank == MASTER_NODE) cout << "sliding interface. " << endl;
      }
      else {
        nVarTransfer = 0;
        transfer_container[donorZone][targetZone] = new CTransfer_ConservativeVars(nVar, nVarTransfer, config_container[donorZone]);
        if (rank == MASTER_NODE) cout << "generic conservative variables. " << endl;  
        }

      break;

      }

      if (config_container[donorZone]->GetBoolMixingPlaneInterface()){
      	nVarTransfer = 0;
      	nVar = solver_container[donorZone][MESH_0][FLOW_SOL]->GetnVar();
      	transfer_container[donorZone][targetZone] = new CTransfer_MixingPlaneInterface(nVar, nVarTransfer, config_container[donorZone], config_container[targetZone]);
        if (rank == MASTER_NODE) cout << "Set mixing-plane interface from donor zone "<< donorZone << " to target zone " << targetZone <<"."<<endl;
      }

    }

  }

#ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
  delete [] Buffer_Recv_mark;
#endif

}

void CDriver::InitStaticMeshMovement(){

  unsigned short iMGlevel;
  unsigned short Kind_Grid_Movement;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for (iZone = 0; iZone < nZone; iZone++) {
    Kind_Grid_Movement = config_container[iZone]->GetKind_GridMovement(iZone);

    switch (Kind_Grid_Movement) {

    case MOVING_WALL:

      /*--- Fixed wall velocities: set the grid velocities only one time
         before the first iteration flow solver. ---*/
      if (rank == MASTER_NODE)
        cout << endl << " Setting the moving wall velocities." << endl;

      surface_movement[iZone]->Moving_Walls(geometry_container[iZone][MESH_0],
          config_container[iZone], iZone, 0);

      /*--- Update the grid velocities on the coarser multigrid levels after
           setting the moving wall velocities for the finest mesh. ---*/

      grid_movement[iZone]->UpdateMultiGrid(geometry_container[iZone], config_container[iZone]);
      break;


    case ROTATING_FRAME:

      /*--- Steadily rotating frame: set the grid velocities just once
         before the first iteration flow solver. ---*/

      if (rank == MASTER_NODE) {
        cout << endl << " Setting rotating frame grid velocities";
        cout << " for zone " << iZone << "." << endl;
      }

      /*--- Set the grid velocities on all multigrid levels for a steadily
           rotating reference frame. ---*/

      for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++){
        geometry_container[iZone][iMGlevel]->SetRotationalVelocity(config_container[iZone], iZone, true);
        geometry_container[iZone][iMGlevel]->SetShroudVelocity(config_container[iZone]);
      }

      break;

    case STEADY_TRANSLATION:

      /*--- Set the translational velocity and hold the grid fixed during
         the calculation (similar to rotating frame, but there is no extra
         source term for translation). ---*/

      if (rank == MASTER_NODE)
        cout << endl << " Setting translational grid velocities." << endl;

      /*--- Set the translational velocity on all grid levels. ---*/

      for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++)
        geometry_container[iZone][iMGlevel]->SetTranslationalVelocity(config_container[iZone], iZone, true);



      break;
    }
  }
}

void CDriver::TurbomachineryPreprocessing(){

  int rank = MASTER_NODE;
  unsigned short donorZone,targetZone, nMarkerInt, iMarkerInt;
  unsigned short nSpanMax = 0;
  bool restart   = (config_container[ZONE_0]->GetRestart() || config_container[ZONE_0]->GetRestart_Flow());
  mixingplane = config_container[ZONE_0]->GetBoolMixingPlaneInterface();
  bool discrete_adjoint = config_container[ZONE_0]->GetDiscrete_Adjoint();
  su2double areaIn, areaOut, nBlades, flowAngleIn, flowAngleOut;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Create turbovertex structure ---*/
  if (rank == MASTER_NODE) cout<<endl<<"Initialize Turbo Vertex Structure." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config_container[iZone]->GetBoolTurbomachinery()){
      geometry_container[iZone][MESH_0]->ComputeNSpan(config_container[iZone], iZone, INFLOW, true);
      geometry_container[iZone][MESH_0]->ComputeNSpan(config_container[iZone], iZone, OUTFLOW, true);
      if (rank == MASTER_NODE) cout <<"Number of span-wise sections in Zone "<< iZone<<": "<< config_container[iZone]->GetnSpanWiseSections() <<"."<< endl;
      if (config_container[iZone]->GetnSpanWiseSections() > nSpanMax){
        nSpanMax = config_container[iZone]->GetnSpanWiseSections();
      }

      config_container[ZONE_0]->SetnSpan_iZones(config_container[iZone]->GetnSpanWiseSections(), iZone);

      geometry_container[iZone][MESH_0]->SetTurboVertex(config_container[iZone], iZone, INFLOW, true);
      geometry_container[iZone][MESH_0]->SetTurboVertex(config_container[iZone], iZone, OUTFLOW, true);
    }
  }

  /*--- Set maximum number of Span among all zones ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config_container[iZone]->GetBoolTurbomachinery()){
      config_container[iZone]->SetnSpanMaxAllZones(nSpanMax);
    }
  }
  if (rank == MASTER_NODE) cout<<"Max number of span-wise sections among all zones: "<< nSpanMax<<"."<< endl;


  if (rank == MASTER_NODE) cout<<"Initialize solver containers for average and performance quantities." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone][MESH_0][FLOW_SOL]->InitTurboContainers(geometry_container[iZone][MESH_0],config_container[iZone]);
  }

//TODO(turbo) make it general for turbo HB
  if (rank == MASTER_NODE) cout<<"Compute inflow and outflow average geometric quantities." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    geometry_container[iZone][MESH_0]->SetAvgTurboValue(config_container[iZone], iZone, INFLOW, true);
    geometry_container[iZone][MESH_0]->SetAvgTurboValue(config_container[iZone],iZone, OUTFLOW, true);
    geometry_container[iZone][MESH_0]->GatherInOutAverageValues(config_container[iZone], true);
  }


  if(mixingplane){
    if (rank == MASTER_NODE) cout << "Set span-wise sections between zones on Mixing-Plane interface." << endl;
    for (donorZone = 0; donorZone < nZone; donorZone++) {
      for (targetZone = 0; targetZone < nZone; targetZone++) {
        if (targetZone != donorZone){
          transfer_container[donorZone][targetZone]->SetSpanWiseLevels(config_container[donorZone], config_container[targetZone]);
        }
      }
    }
  }

  if (rank == MASTER_NODE) cout << "Transfer average geometric quantities to zone 0." << endl;
  for (iZone = 1; iZone < nZone; iZone++) {
    transfer_container[iZone][ZONE_0]->GatherAverageTurboGeoValues(geometry_container[iZone][MESH_0],geometry_container[ZONE_0][MESH_0], iZone);
  }

  /*--- Transfer number of blade to ZONE_0 to correctly compute turbo performance---*/
  for (iZone = 1; iZone < nZone; iZone++) {
    nBlades = config_container[iZone]->GetnBlades(iZone);
    config_container[ZONE_0]->SetnBlades(iZone, nBlades);
  }

  if (rank == MASTER_NODE){
    for (iZone = 0; iZone < nZone; iZone++) {
    areaIn  = geometry_container[iZone][MESH_0]->GetSpanAreaIn(iZone, config_container[iZone]->GetnSpanWiseSections());
    areaOut = geometry_container[iZone][MESH_0]->GetSpanAreaOut(iZone, config_container[iZone]->GetnSpanWiseSections());
    nBlades = config_container[iZone]->GetnBlades(iZone);
    cout << "Inlet area for Row "<< iZone + 1<< ": " << areaIn*10000.0 <<" cm^2."  <<endl;
    cout << "Oulet area for Row "<< iZone + 1<< ": " << areaOut*10000.0 <<" cm^2."  <<endl;
    cout << "Recomputed number of blades for Row "<< iZone + 1 << ": " << nBlades<<"."  <<endl;
    }
  }


  if(mixingplane){
    if (rank == MASTER_NODE) cout<<"Preprocessing of the Mixing-Plane Interface." << endl;
    for (donorZone = 0; donorZone < nZone; donorZone++) {
      nMarkerInt     = config_container[donorZone]->GetnMarker_MixingPlaneInterface()/2;
      for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){
        for (targetZone = 0; targetZone < nZone; targetZone++) {
          if (targetZone != donorZone){
            transfer_container[donorZone][targetZone]->Preprocessing_InterfaceAverage(geometry_container[donorZone][MESH_0], geometry_container[targetZone][MESH_0],
                config_container[donorZone], config_container[targetZone],
                iMarkerInt);
          }
        }
      }
    }
  }

  if(!restart && !discrete_adjoint){
    if (rank == MASTER_NODE) cout<<"Initialize turbomachinery solution quantities." << endl;
    for(iZone = 0; iZone < nZone; iZone++) {
      solver_container[iZone][MESH_0][FLOW_SOL]->SetFreeStream_TurboSolution(config_container[iZone]);
    }
  }

  if (rank == MASTER_NODE) cout<<"Initialize inflow and outflow average solution quantities." << endl;
  for(iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone][MESH_0][FLOW_SOL]->PreprocessAverage(solver_container[iZone][MESH_0], geometry_container[iZone][MESH_0],config_container[iZone],INFLOW);
    solver_container[iZone][MESH_0][FLOW_SOL]->PreprocessAverage(solver_container[iZone][MESH_0], geometry_container[iZone][MESH_0],config_container[iZone],OUTFLOW);
    solver_container[iZone][MESH_0][FLOW_SOL]->TurboAverageProcess(solver_container[iZone][MESH_0], geometry_container[iZone][MESH_0],config_container[iZone],INFLOW);
    solver_container[iZone][MESH_0][FLOW_SOL]->TurboAverageProcess(solver_container[iZone][MESH_0], geometry_container[iZone][MESH_0],config_container[iZone],OUTFLOW);
    solver_container[iZone][MESH_0][FLOW_SOL]->GatherInOutAverageValues(config_container[iZone], geometry_container[iZone][MESH_0]);
    if (rank == MASTER_NODE){
      flowAngleIn = solver_container[iZone][MESH_0][FLOW_SOL]->GetTurboVelocityIn(iZone, config_container[iZone]->GetnSpanWiseSections())[1];
      flowAngleIn /= solver_container[iZone][MESH_0][FLOW_SOL]->GetTurboVelocityIn(iZone, config_container[iZone]->GetnSpanWiseSections())[0];
      flowAngleIn = atan(flowAngleIn)*180.0/PI_NUMBER;
      cout << "Inlet flow angle for Row "<< iZone + 1<< ": "<< flowAngleIn <<"°."  <<endl;
      flowAngleOut = solver_container[iZone][MESH_0][FLOW_SOL]->GetTurboVelocityOut(iZone, config_container[iZone]->GetnSpanWiseSections())[1];
      flowAngleOut /= solver_container[iZone][MESH_0][FLOW_SOL]->GetTurboVelocityOut(iZone, config_container[iZone]->GetnSpanWiseSections())[0];
      flowAngleOut = atan(flowAngleOut)*180.0/PI_NUMBER;
      cout << "Outlet flow angle for Row "<< iZone + 1<< ": "<< flowAngleOut <<"°."  <<endl;

    }
  }

}

void CDriver::StartSolver() {

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Main external loop of the solver. Within this loop, each iteration ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;

  while ( ExtIter < config_container[ZONE_0]->GetnExtIter() ) {

    /*--- Perform some external iteration preprocessing. ---*/

    PreprocessExtIter(ExtIter);

    /*--- Perform a single iteration of the chosen PDE solver. ---*/

    if (!fsi) {

      /*--- Perform a dynamic mesh update if required. ---*/

      DynamicMeshUpdate(ExtIter);

      /*--- Run a single iteration of the problem (fluid, elasticity, wave, heat, ...). ---*/

      Run();

      /*--- Update the solution for dual time stepping strategy ---*/

      Update();

    }
    
    /*--- In the FSIDriver case, mesh and solution updates are already included into the Run function ---*/
    
    else {
      
      Run();
      
    }

    /*--- Monitor the computations after each iteration. ---*/

    Monitor(ExtIter);

    /*--- Output the solution in files. ---*/

    Output(ExtIter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    ExtIter++;

  }

}

void CDriver::PreprocessExtIter(unsigned long ExtIter) {

  /*--- Set the value of the external iteration. ---*/

  for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetExtIter(ExtIter);
  

  /*--- Read the target pressure ---*/

  if (config_container[ZONE_0]->GetInvDesign_Cp() == YES)
    output->SetCp_InverseDesign(solver_container[ZONE_0][MESH_0][FLOW_SOL],
        geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);

  /*--- Read the target heat flux ---*/

  if (config_container[ZONE_0]->GetInvDesign_HeatFlux() == YES)
    output->SetHeatFlux_InverseDesign(solver_container[ZONE_0][MESH_0][FLOW_SOL],
        geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);

  /*--- Set the initial condition for EULER/N-S/RANS and for a non FSI simulation ---*/

  if ( (!fsi) &&
      ( (config_container[ZONE_0]->GetKind_Solver() ==  EULER) ||
       (config_container[ZONE_0]->GetKind_Solver() ==  NAVIER_STOKES) ||
       (config_container[ZONE_0]->GetKind_Solver() ==  RANS) ) ) {
        for(iZone = 0; iZone < nZone; iZone++) {
          solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
    }
  }

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

bool CDriver::Monitor(unsigned long ExtIter) {

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  
  UsedTime = (StopTime - StartTime);
  
  
  /*--- Check if there is any change in the runtime parameters ---*/
  
  CConfig *runtime = NULL;
  strcpy(runtime_file_name, "runtime.dat");
  runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
  runtime->SetExtIter(ExtIter);
  delete runtime;
  
  /*--- Update the convergence history file (serial and parallel computations). ---*/
  
  if (!fsi) {
    for (iZone = 0; iZone < nZone; iZone++) {
      output->SetConvHistory_Body(&ConvHist_file[iZone], geometry_container, solver_container,
          config_container, integration_container, false, UsedTime, iZone);
    }
  }

  /*--- Evaluate the new CFL number (adaptive). ---*/

  if (config_container[ZONE_0]->GetCFL_Adapt() == YES){
  	output->SetCFL_Number(solver_container, config_container, ZONE_0);
  }

  /*--- Check whether the current simulation has reached the specified
   convergence criteria, and set StopCalc to true, if so. ---*/
  
  switch (config_container[ZONE_0]->GetKind_Solver()) {
    case EULER: case NAVIER_STOKES: case RANS:
      StopCalc = integration_container[ZONE_0][FLOW_SOL]->GetConvergence(); break;
    case WAVE_EQUATION:
      StopCalc = integration_container[ZONE_0][WAVE_SOL]->GetConvergence(); break;
    case HEAT_EQUATION:
      StopCalc = integration_container[ZONE_0][HEAT_SOL]->GetConvergence(); break;
    case FEM_ELASTICITY:
      StopCalc = integration_container[ZONE_0][FEA_SOL]->GetConvergence(); break;
    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      StopCalc = integration_container[ZONE_0][ADJFLOW_SOL]->GetConvergence(); break;
  }
  
  return StopCalc;
  
}


void CDriver::Output(unsigned long ExtIter) {
  
  int rank = MASTER_NODE;
  unsigned long nExtIter = config_container[ZONE_0]->GetnExtIter();
  bool output_files = false;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Determine whether a solution needs to be written
   after the current iteration ---*/
  
  if (
      
      /*--- General if statements to print output statements ---*/
      
      (ExtIter+1 >= nExtIter) || (StopCalc) ||
      
      /*--- Fixed CL problem ---*/
      
      ((config_container[ZONE_0]->GetFixed_CL_Mode()) &&
       (config_container[ZONE_0]->GetnExtIter()-config_container[ZONE_0]->GetIter_dCL_dAlpha() - 1 == ExtIter)) ||
      
      /*--- Steady problems ---*/
      
      ((ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (ExtIter != 0) &&
       ((config_container[ZONE_0]->GetUnsteady_Simulation() == STEADY) ||
        (config_container[ZONE_0]->GetUnsteady_Simulation() == HARMONIC_BALANCE) ||
        (config_container[ZONE_0]->GetUnsteady_Simulation() == ROTATIONAL_FRAME))) ||
      
      /*--- Unsteady problems ---*/
      
      (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
        (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_STEPPING)) &&
       ((ExtIter == 0) || (ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0))) ||
      
      ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND) && (!fsi) &&
       ((ExtIter == 0) || ((ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0) ||
                           ((ExtIter-1) % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) ||
      
      ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND) && (fsi) &&
       ((ExtIter == 0) || ((ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) ||
      
      ((config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC) &&
       ((ExtIter == 0) || (ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))
      
      ) {
    
    output_files = true;
    
  }
  
  /*--- Determine whether a solution doesn't need to be written
   after the current iteration ---*/
  
  if (config_container[ZONE_0]->GetFixed_CL_Mode()) {
    if (config_container[ZONE_0]->GetnExtIter()-config_container[ZONE_0]->GetIter_dCL_dAlpha() - 1 < ExtIter) output_files = false;
    if (config_container[ZONE_0]->GetnExtIter() - 1 == ExtIter) output_files = true;
  }
  
  /*--- write the solution ---*/
  
  if (output_files) {
    
    if (rank == MASTER_NODE) cout << endl << "-------------------------- File Output Summary --------------------------";
    
    /*--- Execute the routine for writing restart, volume solution,
     surface solution, and surface comma-separated value files. ---*/
    
    output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, ExtIter, nZone);
    
    
    if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl << endl;
    
  }
  
}
CDriver::~CDriver(void) {}

su2double CDriver::Get_Drag() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag, RefDensity, RefArea, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CDrag = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CD();

  return CDrag*factor;
}

su2double CDriver::Get_Lift() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift, RefDensity, RefArea, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CLift = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CL();

  return CLift*factor;
}

su2double CDriver::Get_Mx(){

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMx, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();
  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CMx = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();

  return CMx*factor*RefLengthCoeff;

}

su2double CDriver::Get_My(){

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMy, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();
  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CMy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();

  return CMy*factor*RefLengthCoeff;

}

su2double CDriver::Get_Mz() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMz, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();
  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate moment around z-axis based on coefficients ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CMz = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();

  return CMz*factor*RefLengthCoeff;

}

su2double CDriver::Get_DragCoeff() {

    unsigned short val_iZone = ZONE_0;
    unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
    su2double CDrag;

    CDrag = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CD();

    return CDrag;
}

su2double CDriver::Get_LiftCoeff() {

    unsigned short val_iZone = ZONE_0;
    unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
    su2double CLift;

    CLift = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CL();

    return CLift;
}

unsigned short CDriver::GetMovingMarker() {

  unsigned short IDtoSend,iMarker, jMarker, Moving;
  string Marker_Tag, Moving_Tag;

  IDtoSend = 0;
  for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
    Moving = config_container[ZONE_0]->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config_container[ZONE_0]->GetnMarker_Moving(); jMarker++) {
        Moving_Tag = config_container[ZONE_0]->GetMarker_Moving_TagBound(jMarker);
        Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Moving_Tag) {
          IDtoSend = iMarker;
          break;
        }
      }
    }
  }

  return IDtoSend;

}

unsigned long CDriver::GetNumberVertices(unsigned short iMarker){

  return geometry_container[ZONE_0][MESH_0]->nVertex[iMarker];

}

unsigned long CDriver::GetNumberHaloVertices(unsigned short iMarker){

  unsigned long nHaloVertices, iVertex, iPoint;

  nHaloVertices = 0;
  for(iVertex = 0; iVertex < geometry_container[ZONE_0][MESH_0]->nVertex[iMarker]; iVertex++){
    iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
    if(!(geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain())) nHaloVertices += 1;
  }

  return nHaloVertices;

}

unsigned long CDriver::GetVertexGlobalIndex(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint, GlobalIndex;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  GlobalIndex = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetGlobalIndex();

  return GlobalIndex;

}

bool CDriver::IsAHaloNode(unsigned short iMarker, unsigned short iVertex) {
 
  unsigned long iPoint; 
  
  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  if(geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain()) return false;
  else return true;

}

unsigned long CDriver::GetnExtIter() {

    return config_container[ZONE_0]->GetnExtIter();
}

unsigned long CDriver::GetExtIter(){

  return ExtIter;
}

su2double CDriver::GetUnsteady_TimeStep(){

  return config_container[ZONE_0]->GetDelta_UnstTime();
}

su2double CDriver::GetVertexCoordX(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord();
  return Coord[0];

}

su2double CDriver::GetVertexCoordY(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord();
  return Coord[1];
}

su2double CDriver::GetVertexCoordZ(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  if(nDim == 3) {
    iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
    Coord = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord();
    return Coord[2];
  }
  else {
    return 0.0;
  }

}

bool CDriver::ComputeVertexForces(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;
  unsigned short iDim, jDim;
  su2double *Normal, AreaSquare, Area;
  bool halo;

  unsigned short FinestMesh = config_container[ZONE_0]->GetFinestMesh();

  /*--- Check the kind of fluid problem ---*/
  bool compressible       = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);
  bool viscous_flow       = ((config_container[ZONE_0]->GetKind_Solver() == NAVIER_STOKES) ||
                 (config_container[ZONE_0]->GetKind_Solver() == RANS) );

  /*--- Parameters for the calculations ---*/
  // Pn: Pressure
  // Pinf: Pressure_infinite
  // div_vel: Velocity divergence
  // Dij: Dirac delta
  su2double Pn = 0.0, div_vel = 0.0, Dij = 0.0;
  su2double Viscosity = 0.0;
  su2double Grad_Vel[3][3] = { {0.0, 0.0, 0.0} ,
              {0.0, 0.0, 0.0} ,
              {0.0, 0.0, 0.0} } ;
  su2double Tau[3][3] = { {0.0, 0.0, 0.0} ,
              {0.0, 0.0, 0.0} ,
              {0.0, 0.0, 0.0} } ;

  su2double Pinf = solver_container[ZONE_0][FinestMesh][FLOW_SOL]->GetPressure_Inf();

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  /*--- It is necessary to distinguish the halo nodes from the others, since they introduce non physical forces. ---*/
  if(geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain()) {
    /*--- Get the normal at the vertex: this normal goes inside the fluid domain. ---*/
    Normal = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
    AreaSquare = 0.0;
    for(iDim = 0; iDim < nDim; iDim++) {
      AreaSquare += Normal[iDim]*Normal[iDim];
    }
    Area = sqrt(AreaSquare);

    /*--- Get the values of pressure and viscosity ---*/
    Pn = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetPressure();
    if (viscous_flow) {
      for(iDim=0; iDim<nDim; iDim++) {
        for(jDim=0; jDim<nDim; jDim++) {
          Grad_Vel[iDim][jDim] = solver_container[ZONE_0][FinestMesh][FLOW_SOL]->node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
        }
      }
      Viscosity = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }

    /*--- Calculate the inviscid (pressure) part of tn in the fluid nodes (force units) ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
     PyWrapNodalForce[iDim] = -(Pn-Pinf)*Normal[iDim];     //NB : norm(Normal) = Area
    }

    /*--- Calculate the viscous (shear stress) part of tn in the fluid nodes (force units) ---*/
    if ((incompressible || compressible) && viscous_flow) {
      div_vel = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        div_vel += Grad_Vel[iDim][iDim];
     if (incompressible) div_vel = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) {
       for (jDim = 0 ; jDim < nDim; jDim++) {
         Dij = 0.0; if (iDim == jDim) Dij = 1.0;
         Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*Dij;
         PyWrapNodalForce[iDim] += Tau[iDim][jDim]*Normal[jDim];
        }
      }
    }

    //Divide by local are in case of force density communication.
   for(iDim = 0; iDim < nDim; iDim++) {
     PyWrapNodalForceDensity[iDim] = PyWrapNodalForce[iDim]/Area;
    }

    halo = false;
  }
  else {
    halo = true;
  }

  return halo;

}

su2double CDriver::GetVertexForceX(unsigned short iMarker, unsigned short iVertex) {

  return PyWrapNodalForce[0];

}

su2double CDriver::GetVertexForceY(unsigned short iMarker, unsigned short iVertex) {

  return PyWrapNodalForce[1];

}

su2double CDriver::GetVertexForceZ(unsigned short iMarker, unsigned short iVertex) {

  return PyWrapNodalForce[2];

}

su2double CDriver::GetVertexForceDensityX(unsigned short iMarker, unsigned short iVertex) {
  return PyWrapNodalForceDensity[0];
}

su2double CDriver::GetVertexForceDensityY(unsigned short iMarker, unsigned short iVertex) {
  return PyWrapNodalForceDensity[1];
}

su2double CDriver::GetVertexForceDensityZ(unsigned short iMarker, unsigned short iVertex) {
  return PyWrapNodalForceDensity[2];
}

void CDriver::SetVertexCoordX(unsigned short iMarker, unsigned short iVertex, su2double newPosX) {

  unsigned long iPoint;
  su2double *Coord;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord();

  PyWrapVarCoord[0] = newPosX - Coord[0];

}

void CDriver::SetVertexCoordY(unsigned short iMarker, unsigned short iVertex, su2double newPosY) {

  unsigned long iPoint;
  su2double *Coord;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord();

  PyWrapVarCoord[1] = newPosY - Coord[1];
}

void CDriver::SetVertexCoordZ(unsigned short iMarker, unsigned short iVertex, su2double newPosZ) {

  unsigned long iPoint;
  su2double *Coord;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord();

  if(nDim > 2) {
    PyWrapVarCoord[2] = newPosZ - Coord[2];
  }
  else {
    PyWrapVarCoord[2] = 0.0;
  }
}

su2double CDriver::SetVertexVarCoord(unsigned short iMarker, unsigned short iVertex) {

  su2double nodalVarCoordNorm;

    geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(PyWrapVarCoord);
    nodalVarCoordNorm = sqrt((PyWrapVarCoord[0])*(PyWrapVarCoord[0]) + (PyWrapVarCoord[1])*(PyWrapVarCoord[1]) + (PyWrapVarCoord[2])*(PyWrapVarCoord[2]));

  return nodalVarCoordNorm;

}

su2double CDriver::GetVertexTemperature(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;
  su2double vertexWallTemp(0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain() && compressible){
    vertexWallTemp = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetTemperature();
  }

  return vertexWallTemp;

}

void CDriver::SetVertexTemperature(unsigned short iMarker, unsigned short iVertex, su2double val_WallTemp){

  unsigned long iPoint;
  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if( geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain() && compressible){
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetImposedTemperature(val_WallTemp);
  }
}

bool CDriver::ComputeVertexHeatFluxes(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;
  unsigned short iDim;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;
  su2double GradT[3] = {0.0,0.0,0.0};

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
  bool halo;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain()){
    halo = false;
  }
  else{
    halo = true;
  }

  if(!halo && compressible){
    laminar_viscosity    = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetGradient_Primitive(0, iDim);
      PyWrapNodalHeatFlux[iDim] = -thermal_conductivity*GradT[iDim];
    }
  }

  return halo;
}

su2double CDriver::GetVertexHeatFluxX(unsigned short iMarker, unsigned short iVertex){

  return PyWrapNodalHeatFlux[0];
}

su2double CDriver::GetVertexHeatFluxY(unsigned short iMarker, unsigned short iVertex){

  return PyWrapNodalHeatFlux[1];
}

su2double CDriver::GetVertexHeatFluxZ(unsigned short iMarker, unsigned short iVertex){

  return PyWrapNodalHeatFlux[2];
}

su2double CDriver::GetVertexNormalHeatFlux(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;
  unsigned short iDim;
  su2double vertexWallHeatFlux;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double Area;
  su2double laminar_viscosity, thermal_conductivity, dTdn;
  su2double *Normal, GradT[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);

  vertexWallHeatFlux = 0.0;
  dTdn = 0.0;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain() && compressible){
    Normal = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
    Area = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);

    for (iDim = 0; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    laminar_viscosity    = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    /*Compute wall heat flux (normal to the wall) based on computed temperature gradient*/
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetGradient_Primitive(0, iDim);
      dTdn += GradT[iDim]*UnitNormal[iDim];
    }

    vertexWallHeatFlux = -thermal_conductivity*dTdn;
  }

  return vertexWallHeatFlux;
}

void CDriver::SetVertexNormalHeatFlux(unsigned short iMarker, unsigned short iVertex, su2double val_WallHeatFlux){

  unsigned long iPoint;
  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if( geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain() && compressible){
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetImposedHeatFlux(val_WallHeatFlux);
  }
}

su2double CDriver::GetThermalConductivity(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;

  iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  laminar_viscosity    = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
  thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);

  return thermal_conductivity;

}

vector<su2double> CDriver::GetVertexUnitNormal(unsigned short iMarker, unsigned short iVertex){

  unsigned short iDim;
  su2double *Normal;
  su2double Area;
  vector<su2double> ret_Normal(3, 0.0);

  Normal = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  ret_Normal[0] = Normal[0]/Area;
  ret_Normal[1] = Normal[1]/Area;
  if(nDim>2) ret_Normal[2] = Normal[2]/Area;

  return ret_Normal;


}

vector<string> CDriver::GetAllBoundaryMarkersTag(){

  vector<string> boundariesTagList;
  unsigned short iMarker,nBoundariesMarkers;
  string Marker_Tag;

  nBoundariesMarkers = config_container[ZONE_0]->GetnMarker_All();
  boundariesTagList.resize(nBoundariesMarkers);

  for(iMarker=0; iMarker < nBoundariesMarkers; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    boundariesTagList[iMarker] = Marker_Tag;
  }

  return boundariesTagList;
}

vector<string> CDriver::GetAllMovingMarkersTag(){

  vector<string> movingBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_Moving();
  movingBoundariesTagList.resize(nBoundariesMarker);

  for(iMarker=0; iMarker < nBoundariesMarker; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_Moving_TagBound(iMarker);
    movingBoundariesTagList[iMarker] = Marker_Tag;
  }

  return movingBoundariesTagList;
}

vector<string> CDriver::GetAllCHTMarkersTag(){

  vector<string> CHTBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_CHT();
  CHTBoundariesTagList.resize(nBoundariesMarker);

  for(iMarker=0; iMarker<nBoundariesMarker; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_CHT_TagBound(iMarker);
    CHTBoundariesTagList[iMarker] = Marker_Tag;
  }

  return CHTBoundariesTagList;
}

map<string, int> CDriver::GetAllBoundaryMarkers(){

  map<string, int>  allBoundariesMap;
  unsigned short iMarker, nBoundaryMarkers;
  string Marker_Tag;

  nBoundaryMarkers = config_container[ZONE_0]->GetnMarker_All();

  for(iMarker=0; iMarker < nBoundaryMarkers; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    allBoundariesMap[Marker_Tag] = iMarker;
  }

  return allBoundariesMap;
}

map<string, string> CDriver::GetAllBoundaryMarkersType(){

  map<string, string> allBoundariesTypeMap;
  unsigned short iMarker, KindBC;
  string Marker_Tag, Marker_Type;

  for(iMarker=0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    KindBC = config_container[ZONE_0]->GetMarker_All_KindBC(iMarker);
    switch(KindBC){
      case EULER_WALL:
        Marker_Type = "EULER_WALL";
        break;
      case FAR_FIELD:
        Marker_Type = "FARFIELD";
        break;
      case ISOTHERMAL:
        Marker_Type = "ISOTHERMAL";
        break;
      case HEAT_FLUX:
        Marker_Type = "HEATFLUX";
        break;
      case INLET_FLOW:
        Marker_Type = "INLET_FLOW";
        break;
      case OUTLET_FLOW:
        Marker_Type = "OUTLET_FLOW";
        break;
      case SYMMETRY_PLANE:
        Marker_Type = "SYMMETRY";
        break;
      case SEND_RECEIVE:
        Marker_Type = "SEND_RECEIVE";
        break;
      default:
        Marker_Type = "UNKNOWN_TYPE";
    }
    allBoundariesTypeMap[Marker_Tag] = Marker_Type;
  }

  return allBoundariesTypeMap;
}

CGeneralDriver::CGeneralDriver(char* confFile, unsigned short val_nZone,
                               unsigned short val_nDim, 
                               SU2_Comm MPICommunicator) : CDriver(confFile,
                                                                   val_nZone,
                                                                   val_nDim,
                                                                   MPICommunicator) { }

CGeneralDriver::~CGeneralDriver(void) { }

void CGeneralDriver::Run() {

  unsigned short iZone;

  /*--- Run a single iteration of a fem problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    iteration_container[iZone]->Preprocess(output, integration_container, geometry_container,
                                           solver_container, numerics_container, config_container,
                                           surface_movement, grid_movement, FFDBox, iZone);

    iteration_container[iZone]->Iterate(output, integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox, iZone);
  }

}

void CGeneralDriver::Update() {

  for (iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone]->Update(output, integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox, iZone);

  if (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM){
      iteration_container[ZONE_0]->Postprocess(output, integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox, ZONE_0);
  }

}

void CGeneralDriver::ResetConvergence() {

  switch (config_container[ZONE_0]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
    integration_container[ZONE_0][FLOW_SOL]->SetConvergence(false);
      if (config_container[ZONE_0]->GetKind_Solver() == RANS) integration_container[ZONE_0][TURB_SOL]->SetConvergence(false);
      if(config_container[ZONE_0]->GetKind_Trans_Model() == LM) integration_container[ZONE_0][TRANS_SOL]->SetConvergence(false);
    break;

  case WAVE_EQUATION:
    integration_container[ZONE_0][WAVE_SOL]->SetConvergence(false);
    break;

  case HEAT_EQUATION:
    integration_container[ZONE_0][HEAT_SOL]->SetConvergence(false);
    break;

  case POISSON_EQUATION:
    break;

  case FEM_ELASTICITY:
    integration_container[ZONE_0][FEA_SOL]->SetConvergence(false);
    break;

  case DISC_ADJ_FEM:
    integration_container[ZONE_0][ADJFEA_SOL]->SetConvergence(false);
    break;

  case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    integration_container[ZONE_0][ADJFLOW_SOL]->SetConvergence(false);
      if( (config_container[ZONE_0]->GetKind_Solver() == ADJ_RANS) || (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS) )
      integration_container[ZONE_0][ADJTURB_SOL]->SetConvergence(false);
    break;

  }

}

void CGeneralDriver::DynamicMeshUpdate(unsigned long ExtIter) {

  bool harmonic_balance;

  for (iZone = 0; iZone < nZone; iZone++) {
   harmonic_balance = (config_container[iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE);
    /*--- Dynamic mesh update ---*/
    if ((config_container[iZone]->GetGrid_Movement()) && (!harmonic_balance)) {
      iteration_container[iZone]->SetGrid_Movement(geometry_container, surface_movement, grid_movement, FFDBox, solver_container, config_container, iZone, 0, ExtIter );
    }
  }
}

void CGeneralDriver::StaticMeshUpdate() {

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if(rank == MASTER_NODE) cout << " Deforming the volume grid." << endl;
  grid_movement[ZONE_0]->SetVolume_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], true);

  if(rank == MASTER_NODE) cout << "No grid velocity to be computed : static grid deformation." << endl;

  if(rank == MASTER_NODE) cout << " Updating multigrid structure." << endl;
  grid_movement[ZONE_0]->UpdateMultiGrid(geometry_container[ZONE_0], config_container[ZONE_0]);

}

void CGeneralDriver::SetInitialMesh() {

  unsigned long iPoint;

  StaticMeshUpdate();

  /*--- Propagate the initial deformation to the past ---*/
  //if (!restart) {
  for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
      for(iPoint = 0; iPoint < geometry_container[ZONE_0][iMesh]->GetnPoint(); iPoint++) {
      //solver_container[ZONE_0][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
      //solver_container[ZONE_0][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
      geometry_container[ZONE_0][iMesh]->node[iPoint]->SetVolume_n();
      geometry_container[ZONE_0][iMesh]->node[iPoint]->SetVolume_nM1();
      geometry_container[ZONE_0][iMesh]->node[iPoint]->SetCoord_n();
      geometry_container[ZONE_0][iMesh]->node[iPoint]->SetCoord_n1();
    }
  }
  //}

}

CFluidDriver::CFluidDriver(char* confFile, unsigned short val_nZone, unsigned short val_nDim, SU2_Comm MPICommunicator) : CDriver(confFile, val_nZone, val_nDim, MPICommunicator) { }

CFluidDriver::~CFluidDriver(void) { }

void CFluidDriver::Run() {

  unsigned short iZone, jZone, checkConvergence;
  unsigned long IntIter, nIntIter;
  bool unsteady;

  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  unsteady = (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  /*--- Zone preprocessing ---*/

  for (iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone]->Preprocess(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone);

  /*--- Updating zone interface communication patterns,
   needed only for unsteady simulation since for steady problems
   this is done once in the interpolator_container constructor 
   at the beginning of the computation ---*/

  if ( unsteady ) {
    for (iZone = 0; iZone < nZone; iZone++) {   
      for (jZone = 0; jZone < nZone; jZone++)
        if(jZone != iZone && interpolator_container[iZone][jZone] != NULL)
        interpolator_container[iZone][jZone]->Set_TransferCoeff(config_container);
    }
  }

  /*--- Begin Unsteady pseudo-time stepping internal loop, if not unsteady it does only one step --*/

  if (unsteady) 
    nIntIter = config_container[MESH_0]->GetUnst_nIntIter();
  else
    nIntIter = 1;

  for (IntIter = 0; IntIter < nIntIter; IntIter++) {

    /*--- At each pseudo time-step updates transfer data ---*/
    for (iZone = 0; iZone < nZone; iZone++)   
      for (jZone = 0; jZone < nZone; jZone++)
        if(jZone != iZone && transfer_container[iZone][jZone] != NULL)
          Transfer_Data(iZone, jZone);

    /*--- For each zone runs one single iteration ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      config_container[iZone]->SetIntIter(IntIter);
      iteration_container[iZone]->Iterate(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone);
    }

    /*--- Check convergence in each zone --*/

    checkConvergence = 0;
    for (iZone = 0; iZone < nZone; iZone++)
    checkConvergence += (int) integration_container[iZone][FLOW_SOL]->GetConvergence();

    /*--- If convergence was reached in every zone --*/

  if (checkConvergence == nZone) break;
  }

}

void CFluidDriver::Transfer_Data(unsigned short donorZone, unsigned short targetZone) {

#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  bool MatchingMesh = config_container[targetZone]->GetMatchingMesh();

  /*--- Select the transfer method and the appropriate mesh properties (matching or nonmatching mesh) ---*/

  switch (config_container[targetZone]->GetKind_TransferMethod()) {

  case BROADCAST_DATA:
      if (MatchingMesh) {
        transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
            geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
            config_container[donorZone], config_container[targetZone]);
        if (config_container[targetZone]->GetKind_Solver() == RANS)
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][MESH_0][TURB_SOL],solver_container[targetZone][MESH_0][TURB_SOL],
              geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
              config_container[donorZone], config_container[targetZone]);
      }
      else {
        transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
            geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
            config_container[donorZone], config_container[targetZone]);
        if (config_container[targetZone]->GetKind_Solver() == RANS)
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][MESH_0][TURB_SOL],solver_container[targetZone][MESH_0][TURB_SOL],
              geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
              config_container[donorZone], config_container[targetZone]);
    }
    break;

  case SCATTER_DATA:
    if (MatchingMesh) {
      transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
          geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
          config_container[donorZone], config_container[targetZone]);
      if (config_container[targetZone]->GetKind_Solver() == RANS)
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][MESH_0][TURB_SOL],solver_container[targetZone][MESH_0][TURB_SOL],
            geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
            config_container[donorZone], config_container[targetZone]);
    }
    else {
      cout << "Scatter method not implemented for non-matching meshes. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  case ALLGATHER_DATA:
    if (MatchingMesh) {
      cout << "Allgather method not yet implemented for matching meshes. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    else {
      transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
          geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
          config_container[donorZone], config_container[targetZone]);
      if (config_container[targetZone]->GetKind_Solver() == RANS)
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][MESH_0][TURB_SOL],solver_container[targetZone][MESH_0][TURB_SOL],
            geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
            config_container[donorZone], config_container[targetZone]);
    }
    break;
  }

}

void CFluidDriver::Update() {

  for(iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone]->Update(output, integration_container, geometry_container,
         solver_container, numerics_container, config_container,
         surface_movement, grid_movement, FFDBox, iZone);
}

void CFluidDriver::ResetConvergence() {

  for(iZone = 0; iZone < nZone; iZone++) {
    switch (config_container[iZone]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
      integration_container[iZone][FLOW_SOL]->SetConvergence(false);
      if (config_container[iZone]->GetKind_Solver() == RANS) integration_container[iZone][TURB_SOL]->SetConvergence(false);
      if(config_container[iZone]->GetKind_Trans_Model() == LM) integration_container[iZone][TRANS_SOL]->SetConvergence(false);
      break;

    case WAVE_EQUATION:
      integration_container[iZone][WAVE_SOL]->SetConvergence(false);
      break;

    case HEAT_EQUATION:
      integration_container[iZone][HEAT_SOL]->SetConvergence(false);
      break;

    case POISSON_EQUATION:
      break;

    case FEM_ELASTICITY:
      integration_container[iZone][FEA_SOL]->SetConvergence(false);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      integration_container[iZone][ADJFLOW_SOL]->SetConvergence(false);
      if( (config_container[iZone]->GetKind_Solver() == ADJ_RANS) || (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS) )
        integration_container[iZone][ADJTURB_SOL]->SetConvergence(false);
      break;
    }
  }

}

void CFluidDriver::DynamicMeshUpdate(unsigned long ExtIter) {

  bool harmonic_balance;

  for (iZone = 0; iZone < nZone; iZone++) {
   harmonic_balance = (config_container[iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE);
    /*--- Dynamic mesh update ---*/
    if ((config_container[iZone]->GetGrid_Movement()) && (!harmonic_balance)) {
      iteration_container[iZone]->SetGrid_Movement(geometry_container, surface_movement, grid_movement, FFDBox, solver_container, config_container, iZone, 0, ExtIter );
    }
  }

}

void CFluidDriver::StaticMeshUpdate() {

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for(iZone = 0; iZone < nZone; iZone++) {
    if(rank == MASTER_NODE) cout << " Deforming the volume grid." << endl;
    grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone][MESH_0], config_container[iZone], true);

    if(rank == MASTER_NODE) cout << "No grid velocity to be computde : static grid deformation." << endl;

    if(rank == MASTER_NODE) cout << " Updating multigrid structure." << endl;
    grid_movement[iZone]->UpdateMultiGrid(geometry_container[iZone], config_container[iZone]);
  }
}

void CFluidDriver::SetInitialMesh() {

  unsigned long iPoint;

  StaticMeshUpdate();

  /*--- Propagate the initial deformation to the past ---*/
  //if (!restart) {
    for(iZone = 0; iZone < nZone; iZone++) {
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
        for(iPoint = 0; iPoint < geometry_container[iZone][iMesh]->GetnPoint(); iPoint++) {
        //solver_container[iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        //solver_container[iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
        geometry_container[iZone][iMesh]->node[iPoint]->SetVolume_n();
        geometry_container[iZone][iMesh]->node[iPoint]->SetVolume_nM1();
        geometry_container[iZone][iMesh]->node[iPoint]->SetCoord_n();
        geometry_container[iZone][iMesh]->node[iPoint]->SetCoord_n1();
      }
    }
  }
  //}
}


void CFluidDriver::MGUpdateBoundaryConditions_HeatFlux(unsigned short val_marker){

  unsigned short iMGfine, iMGlevel, nMGlevel;

  for(iZone = 0; iZone < nZone; iZone++){
    nMGlevel = config_container[iZone]->GetnMGLevels();
    for (iMGlevel=1; iMGlevel <= nMGlevel; iMGlevel++){
      iMGfine = iMGlevel-1;
      geometry_container[iZone][iMGlevel]->SetWallHeatFlux(geometry_container[iZone][iMGfine], val_marker);
    }
  }
}

void CFluidDriver::MGUpdateBoundaryConditions_Temperature(unsigned short val_marker){

  unsigned short iMGfine, iMGlevel, nMGlevel;

  for(iZone = 0; iZone < nZone; iZone++){
    nMGlevel = config_container[iZone]->GetnMGLevels();
    for (iMGlevel=1; iMGlevel <= nMGlevel; iMGlevel++){
      iMGfine = iMGlevel-1;
      geometry_container[iZone][iMGlevel]->SetWallTemperature(geometry_container[iZone][iMGfine], val_marker);
    }
  }
}

CTurbomachineryDriver::CTurbomachineryDriver(char* confFile,
    unsigned short val_nZone,
    unsigned short val_nDim, SU2_Comm MPICommunicator) : CFluidDriver(confFile,
        val_nZone,
        val_nDim,
        MPICommunicator) { }

CTurbomachineryDriver::~CTurbomachineryDriver(void) { }

void CTurbomachineryDriver::Run() {


//  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone]->Preprocess(output, integration_container, geometry_container,
                                           solver_container, numerics_container, config_container,
                                           surface_movement, grid_movement, FFDBox, iZone);
  }

  /* --- Update the mixing-plane interface ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if(mixingplane)SetMixingPlane(iZone);
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone]->Iterate(output, integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox, iZone);
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone]->Postprocess(output, integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox, iZone);
  }

  if (rank == MASTER_NODE){
    SetTurboPerformance(ZONE_0);
  }


}

void CTurbomachineryDriver::SetMixingPlane(unsigned short donorZone){

  unsigned short targetZone, nMarkerInt, iMarkerInt ;
  nMarkerInt     = config_container[donorZone]->GetnMarker_MixingPlaneInterface()/2;

  /* --- transfer the average value from the donorZone to the targetZone*/
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){
    for (targetZone = 0; targetZone < nZone; targetZone++) {
      if (targetZone != donorZone){
        transfer_container[donorZone][targetZone]->Allgather_InterfaceAverage(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
            geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
            config_container[donorZone], config_container[targetZone], iMarkerInt );
      }
    }
  }
}

void CTurbomachineryDriver::SetTurboPerformance(unsigned short targetZone){

  unsigned short donorZone;
  //IMPORTANT this approach of multi-zone performances rely upon the fact that turbomachinery markers follow the natural (stator-rotor) development of the real machine.
  /* --- transfer the local turboperfomance quantities (for each blade)  from all the donorZones to the targetZone (ZONE_0) ---*/
  for (donorZone = 1; donorZone < nZone; donorZone++) {
    transfer_container[donorZone][targetZone]->GatherAverageValues(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FLOW_SOL], donorZone);
  }

  /* --- compute turboperformance for each stage and the global machine ---*/

  output->ComputeTurboPerformance(solver_container[targetZone][MESH_0][FLOW_SOL], geometry_container[targetZone][MESH_0], config_container[targetZone]);

}


bool CTurbomachineryDriver::Monitor(unsigned long ExtIter) {

  su2double CFL;
  su2double rot_z_ini, rot_z_final ,rot_z;
  su2double outPres_ini, outPres_final, outPres;
  unsigned long rampFreq, finalRamp_Iter;
  unsigned short iMarker, KindBC, KindBCOption;
  string Marker_Tag;

  int rank = MASTER_NODE;
  bool print;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif

  UsedTime = (StopTime - StartTime);


  /*--- Check if there is any change in the runtime parameters ---*/
  CConfig *runtime = NULL;
  strcpy(runtime_file_name, "runtime.dat");
  runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
  runtime->SetExtIter(ExtIter);
  delete runtime;

  /*--- Update the convergence history file (serial and parallel computations). ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    output->SetConvHistory_Body(&ConvHist_file[iZone], geometry_container, solver_container,
        config_container, integration_container, false, UsedTime, iZone);
  }


  /*--- Evaluate the new CFL number (adaptive). ---*/
  if (config_container[ZONE_0]->GetCFL_Adapt() == YES) {
    if(mixingplane){
      CFL = 0;
      for (iZone = 0; iZone < nZone; iZone++){
        output->SetCFL_Number(solver_container, config_container, iZone);
        CFL += config_container[iZone]->GetCFL(MESH_0);
      }
      /*--- For fluid-multizone the new CFL number is the same for all the zones and it is equal to the zones' minimum value. ---*/
      for (iZone = 0; iZone < nZone; iZone++){
        config_container[iZone]->SetCFL(MESH_0, CFL/nZone);
      }
    }
    else{
      output->SetCFL_Number(solver_container, config_container, ZONE_0);
    }
  }


  /*--- ROTATING FRAME Ramp: Compute the updated rotational velocity. ---*/
  if (config_container[ZONE_0]->GetGrid_Movement() && config_container[ZONE_0]->GetRampRotatingFrame()) {
    rampFreq       = SU2_TYPE::Int(config_container[ZONE_0]->GetRampRotatingFrame_Coeff(1));
    finalRamp_Iter = SU2_TYPE::Int(config_container[ZONE_0]->GetRampRotatingFrame_Coeff(2));
    rot_z_ini = config_container[ZONE_0]->GetRampRotatingFrame_Coeff(0);
    print = false;
    if(ExtIter % rampFreq == 0 &&  ExtIter <= finalRamp_Iter){

      for (iZone = 0; iZone < nZone; iZone++) {
        rot_z_final = config_container[iZone]->GetFinalRotation_Rate_Z(iZone);
        if(abs(rot_z_final) > 0.0){
          rot_z = rot_z_ini + ExtIter*( rot_z_final - rot_z_ini)/finalRamp_Iter;
          config_container[iZone]->SetRotation_Rate_Z(rot_z, iZone);
          if(rank == MASTER_NODE && print && ExtIter > 0) {
            cout << endl << " Updated rotating frame grid velocities";
            cout << " for zone " << iZone << "." << endl;
          }
          geometry_container[iZone][MESH_0]->SetRotationalVelocity(config_container[iZone], iZone, print);
          geometry_container[iZone][MESH_0]->SetShroudVelocity(config_container[iZone]);
        }
      }

      for (iZone = 0; iZone < nZone; iZone++) {
        geometry_container[iZone][MESH_0]->SetAvgTurboValue(config_container[iZone], iZone, INFLOW, false);
        geometry_container[iZone][MESH_0]->SetAvgTurboValue(config_container[iZone],iZone, OUTFLOW, false);
        geometry_container[iZone][MESH_0]->GatherInOutAverageValues(config_container[iZone], false);

      }

      for (iZone = 1; iZone < nZone; iZone++) {
        transfer_container[iZone][ZONE_0]->GatherAverageTurboGeoValues(geometry_container[iZone][MESH_0],geometry_container[ZONE_0][MESH_0], iZone);
      }

    }
  }


  /*--- Outlet Pressure Ramp: Compute the updated rotational velocity. ---*/
  if (config_container[ZONE_0]->GetRampOutletPressure()) {
    rampFreq       = SU2_TYPE::Int(config_container[ZONE_0]->GetRampOutletPressure_Coeff(1));
    finalRamp_Iter = SU2_TYPE::Int(config_container[ZONE_0]->GetRampOutletPressure_Coeff(2));
    outPres_ini    = config_container[ZONE_0]->GetRampOutletPressure_Coeff(0);
    outPres_final  = config_container[ZONE_0]->GetFinalOutletPressure();

    if(ExtIter % rampFreq == 0 &&  ExtIter <= finalRamp_Iter){
      outPres = outPres_ini + ExtIter*(outPres_final - outPres_ini)/finalRamp_Iter;
      if(rank == MASTER_NODE) config_container[ZONE_0]->SetMonitotOutletPressure(outPres);

      for (iZone = 0; iZone < nZone; iZone++) {
        for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++) {
          KindBC = config_container[iZone]->GetMarker_All_KindBC(iMarker);
          switch (KindBC) {
          case RIEMANN_BOUNDARY:
            Marker_Tag         = config_container[iZone]->GetMarker_All_TagBound(iMarker);
            KindBCOption       = config_container[iZone]->GetKind_Data_Riemann(Marker_Tag);
            if(KindBCOption == STATIC_PRESSURE || KindBCOption == RADIAL_EQUILIBRIUM ){
              cout << "Outlet pressure ramp only implemented for NRBC" <<endl;
              exit(EXIT_FAILURE);
            }
            break;
          case GILES_BOUNDARY:
            Marker_Tag         = config_container[iZone]->GetMarker_All_TagBound(iMarker);
            KindBCOption       = config_container[iZone]->GetKind_Data_Giles(Marker_Tag);
            if(KindBCOption == STATIC_PRESSURE || KindBCOption == STATIC_PRESSURE_1D || KindBCOption == RADIAL_EQUILIBRIUM ){
              config_container[iZone]->SetGiles_Var1(outPres, Marker_Tag);
            }
            break;
          }
        }
      }
    }
  }


  /*--- Check whether the current simulation has reached the specified
   convergence criteria, and set StopCalc to true, if so. ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
  case EULER: case NAVIER_STOKES: case RANS:
    StopCalc = integration_container[ZONE_0][FLOW_SOL]->GetConvergence(); break;
  case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    StopCalc = integration_container[ZONE_0][ADJFLOW_SOL]->GetConvergence(); break;
  }

  return StopCalc;

}

CDiscAdjFluidDriver::CDiscAdjFluidDriver(char* confFile,
                                                 unsigned short val_nZone,
                                                 unsigned short val_nDim, SU2_Comm MPICommunicator) : CFluidDriver(confFile,
                                                                                    val_nZone,
                                                                                    val_nDim, MPICommunicator) {
  RecordingState = NONE;
  unsigned short iZone;

  direct_iteration = new CIteration*[nZone];

  for (iZone = 0; iZone < nZone; iZone++){
    if(config_container[iZone]->GetBoolTurbomachinery()){
      direct_iteration[iZone] = new CTurboIteration(config_container[iZone]);
    }
    else{
      direct_iteration[iZone] = new CFluidIteration(config_container[iZone]);
    }
  }

}

CDiscAdjFluidDriver::~CDiscAdjFluidDriver(){

  for (iZone = 0; iZone < nZone; iZone++){
    delete direct_iteration[iZone];
  }

  delete [] direct_iteration;

}

void CDiscAdjFluidDriver::Run() {

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

    iteration_container[iZone]->Preprocess(output, integration_container, geometry_container,
                                                     solver_container, numerics_container, config_container,
                                                     surface_movement, grid_movement, FFDBox, iZone);
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

      iteration_container[iZone]->InitializeAdjoint(solver_container, geometry_container, config_container, iZone);

    }

    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdj_ObjFunction();

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      iteration_container[iZone]->Iterate(output, integration_container, geometry_container,
                                          solver_container, numerics_container, config_container,
                                          surface_movement, grid_movement, FFDBox, iZone);
    }

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    /*--- Check convergence in each zone --*/

    checkConvergence = 0;
    for (iZone = 0; iZone < nZone; iZone++)
      checkConvergence += (int) integration_container[iZone][ADJFLOW_SOL]->GetConvergence();

    /*--- If convergence was reached in every zone --*/

    if (checkConvergence == nZone) break;

    /*--- Write the convergence history (only screen output) ---*/

    if (unsteady)
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

  }

  /*--- Compute the geometrical sensitivities ---*/

  if ((ExtIter+1 >= config_container[ZONE_0]->GetnExtIter()) ||
      integration_container[ZONE_0][ADJFLOW_SOL]->GetConvergence() ||
      (ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) || unsteady){

    /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
     * as argument ensures that all information from a previous recording is removed. ---*/

    SetRecording(NONE);

    /*--- Store the computational graph of one direct iteration with the mesh coordinates as input. ---*/

    SetRecording(MESH_COORDS);

    /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
     *    of the current iteration. The values are passed to the AD tool. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {

      iteration_container[iZone]->InitializeAdjoint(solver_container, geometry_container, config_container, iZone);

    }

    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdj_ObjFunction();

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed sensitivity values. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]);
    }

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();
  }
}

void CDiscAdjFluidDriver::SetRecording(unsigned short kind_recording){
  unsigned short iZone, iMesh;
  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  AD::Reset();

  /*--- Prepare for recording by resetting the flow solution to the initial converged solution---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
      solver_container[iZone][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[iZone][iMesh], config_container[iZone]);
    }
    if (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS && !config_container[iZone]->GetFrozen_Visc_Disc()) {
      solver_container[iZone][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[iZone][MESH_0], config_container[iZone]);
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
      iteration_container[iZone]->RegisterInput(solver_container, geometry_container, config_container, iZone, kind_recording);
    }

  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone]->SetDependencies(solver_container, geometry_container, config_container, iZone, kind_recording);
  }

  /*--- Do one iteration of the direct flow solver ---*/

  DirectRun();

  /*--- Print residuals in the first iteration ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    if (rank == MASTER_NODE && ((ExtIter == 0) || (config_container[iZone]->GetUnsteady_Simulation() != STEADY)) && (kind_recording == FLOW_CONS_VARS)) {
      cout << " Zone " << iZone << ": log10[Conservative 0]: "<< log10(solver_container[iZone][MESH_0][FLOW_SOL]->GetRes_RMS(0)) << endl;
      if ( config_container[iZone]->GetKind_Turb_Model() != NONE && !config_container[iZone]->GetFrozen_Visc_Disc()) {
        cout <<"       log10[RMS k]: " << log10(solver_container[iZone][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
      }
    }
  }

  RecordingState = kind_recording;

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone]->RegisterOutput(solver_container, geometry_container, config_container, output, iZone);
  }

  /*--- Extract the objective function and store it --- */

  SetObjFunction();

  AD::StopRecording();

}

void CDiscAdjFluidDriver::SetAdj_ObjFunction(){

  int rank = MASTER_NODE;

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

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE){
    SU2_TYPE::SetDerivative(ObjFunc, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc, 0.0);
  }

}

void CDiscAdjFluidDriver::SetObjFunction(){

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  ObjFunc = 0.0;

  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][FLOW_SOL]->SetTotal_ComboObj(0.0);
  }

  /*--- Specific scalar objective functions ---*/

  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
        
        if (config_container[ZONE_0]->GetnMarker_Analyze() != 0)
          output->SpecialOutput_AnalyzeSurface(solver_container[iZone][MESH_0][FLOW_SOL], geometry_container[iZone][MESH_0], config_container[iZone], false);
        
        if (config_container[ZONE_0]->GetnMarker_Analyze() != 0)
          output->SpecialOutput_Distortion(solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], false);
        
        if (config_container[ZONE_0]->GetnMarker_NearFieldBound() != 0)
          output->SpecialOutput_SonicBoom(solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], false);
          
        if (config_container[ZONE_0]->GetPlot_Section_Forces())
          output->SpecialOutput_SpanLoad(solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], false);
        
        break;
    }
  }

  /*--- Surface based obj. function ---*/

  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][FLOW_SOL]->Evaluate_ObjFunc(config_container[iZone]);
    ObjFunc += solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_ComboObj();
  }

  if (rank == MASTER_NODE){
    AD::RegisterOutput(ObjFunc);
  }
  
}

void CDiscAdjFluidDriver::DirectRun(){


  unsigned short iZone, jZone;
  bool unsteady = config_container[ZONE_0]->GetUnsteady_Simulation() != STEADY;

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif



  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  unsteady = (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  /*--- Zone preprocessing ---*/

  for (iZone = 0; iZone < nZone; iZone++)
    direct_iteration[iZone]->Preprocess(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone);

  /*--- Updating zone interface communication patterns,
   needed only for unsteady simulation since for steady problems
  this is done once in the interpolator_container constructor
   at the beginning of the computation ---*/

  if ( unsteady ) {
  for (iZone = 0; iZone < nZone; iZone++) {
      for (jZone = 0; jZone < nZone; jZone++)
        if(jZone != iZone && interpolator_container[iZone][jZone] != NULL)
        interpolator_container[iZone][jZone]->Set_TransferCoeff(config_container);
    }
  }

  /*--- Do one iteration of the direct solver  --*/

  /*--- At each pseudo time-step updates transfer data ---*/
  for (iZone = 0; iZone < nZone; iZone++)
    for (jZone = 0; jZone < nZone; jZone++)
      if(jZone != iZone && transfer_container[iZone][jZone] != NULL)
        Transfer_Data(iZone, jZone);

  /*--- For each zone runs one single iteration ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]->SetIntIter(1);
    direct_iteration[iZone]->Iterate(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone);
  }

}

CDiscAdjTurbomachineryDriver::CDiscAdjTurbomachineryDriver(char* confFile,
                                                           unsigned short val_nZone,
                                                           unsigned short val_nDim,
                                                           SU2_Comm MPICommunicator): CDiscAdjFluidDriver(confFile, val_nZone, val_nDim, MPICommunicator){ }
CDiscAdjTurbomachineryDriver::~CDiscAdjTurbomachineryDriver(){

}


void CDiscAdjTurbomachineryDriver::DirectRun(){

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    direct_iteration[iZone]->Preprocess(output, integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, iZone);

  }


  /* --- Update the mixing-plane interface ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if(mixingplane)SetMixingPlane(iZone);
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    direct_iteration[iZone]->Iterate(output, integration_container, geometry_container,
                                     solver_container, numerics_container, config_container,
                                     surface_movement, grid_movement, FFDBox, iZone);
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    direct_iteration[iZone]->Postprocess(output, integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox, iZone);
  }


  if (rank == MASTER_NODE){
    SetTurboPerformance(ZONE_0);
  }

}

void CDiscAdjTurbomachineryDriver::SetObjFunction(){


  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetTotal_ComboObj(0.0);

  switch (config_container[ZONE_0]->GetKind_ObjFunc()){
  case ENTROPY_GENERATION:
    solver_container[ZONE_0][MESH_0][FLOW_SOL]->AddTotal_ComboObj(output->GetEntropyGen(config_container[ZONE_0]->GetnMarker_TurboPerformance() - 1, config_container[ZONE_0]->GetnSpanWiseSections()));
    break;
  case FLOW_ANGLE_OUT:
      solver_container[ZONE_0][MESH_0][FLOW_SOL]->AddTotal_ComboObj(output->GetFlowAngleOut(config_container[ZONE_0]->GetnMarker_TurboPerformance() - 1, config_container[ZONE_0]->GetnSpanWiseSections()));
      break;
  case MASS_FLOW_IN:
    solver_container[ZONE_0][MESH_0][FLOW_SOL]->AddTotal_ComboObj(output->GetMassFlowIn(config_container[ZONE_0]->GetnMarker_TurboPerformance() - 1, config_container[ZONE_0]->GetnSpanWiseSections()));
    break;
  default:
    break;
  }

  ObjFunc = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_ComboObj();

  if (rank == MASTER_NODE){
    AD::RegisterOutput(ObjFunc);
  }
}

void CDiscAdjTurbomachineryDriver::SetMixingPlane(unsigned short donorZone){

  unsigned short targetZone, nMarkerInt, iMarkerInt ;
  nMarkerInt     = config_container[donorZone]->GetnMarker_MixingPlaneInterface()/2;

  /* --- transfer the average value from the donorZone to the targetZone*/
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){
    for (targetZone = 0; targetZone < nZone; targetZone++) {
      if (targetZone != donorZone){
        transfer_container[donorZone][targetZone]->Allgather_InterfaceAverage(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
            geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
            config_container[donorZone], config_container[targetZone], iMarkerInt );
      }
    }
  }
}


void CDiscAdjTurbomachineryDriver::SetTurboPerformance(unsigned short targetZone){

  unsigned short donorZone;
  //IMPORTANT this approach of multi-zone performances rely upon the fact that turbomachinery markers follow the natural (stator-rotor) development of the real machine.
  /* --- transfer the local turboperfomance quantities (for each blade)  from all the donorZones to the targetZone (ZONE_0) ---*/
  for (donorZone = 1; donorZone < nZone; donorZone++) {
    transfer_container[donorZone][targetZone]->GatherAverageValues(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FLOW_SOL], donorZone);
  }

  /* --- compute turboperformance for each stage and the global machine ---*/

  output->ComputeTurboPerformance(solver_container[targetZone][MESH_0][FLOW_SOL], geometry_container[targetZone][MESH_0], config_container[targetZone]);

}

CHBDriver::CHBDriver(char* confFile,
    unsigned short val_nZone,
    unsigned short val_nDim,
    SU2_Comm MPICommunicator) : CDriver(confFile,
        val_nZone,
        val_nDim,
        MPICommunicator) {
  unsigned short kZone;

  D = NULL;
  /*--- allocate dynamic memory for the Harmonic Balance operator ---*/
  D = new su2double*[nZone]; for (kZone = 0; kZone < nZone; kZone++) D[kZone] = new su2double[nZone];

}

CHBDriver::~CHBDriver(void) {

  unsigned short kZone;

  /*--- delete dynamic memory for the Harmonic Balance operator ---*/
  for (kZone = 0; kZone < nZone; kZone++) if (D[kZone] != NULL) delete [] D[kZone];
  if (D[kZone] != NULL) delete [] D;

}

void CHBDriver::Run() {

  /*--- Run a single iteration of a Harmonic Balance problem. Preprocess all
   all zones before beginning the iteration. ---*/

  for (iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone]->Preprocess(output, integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, iZone);

  for (iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone]->Iterate(output, integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, iZone);

}

void CHBDriver::Update() {

  for (iZone = 0; iZone < nZone; iZone++) {
    /*--- Compute the harmonic balance terms across all zones ---*/
    SetHarmonicBalance(iZone);

  }

  /*--- Precondition the harmonic balance source terms ---*/
  if (config_container[ZONE_0]->GetHB_Precondition() == YES) {
    StabilizeHarmonicBalance();

  }

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Update the harmonic balance terms across all zones ---*/
    iteration_container[iZone]->Update(output, integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, iZone);

  }

}

void CHBDriver::ResetConvergence() {

  for(iZone = 0; iZone < nZone; iZone++) {
    switch (config_container[iZone]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
      integration_container[iZone][FLOW_SOL]->SetConvergence(false);
      if (config_container[iZone]->GetKind_Solver() == RANS) integration_container[iZone][TURB_SOL]->SetConvergence(false);
      if(config_container[iZone]->GetKind_Trans_Model() == LM) integration_container[iZone][TRANS_SOL]->SetConvergence(false);
      break;

    case WAVE_EQUATION:
      integration_container[iZone][WAVE_SOL]->SetConvergence(false);
      break;

    case HEAT_EQUATION:
      integration_container[iZone][HEAT_SOL]->SetConvergence(false);
      break;

    case POISSON_EQUATION:
      break;

    case FEM_ELASTICITY:
      integration_container[iZone][FEA_SOL]->SetConvergence(false);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      integration_container[iZone][ADJFLOW_SOL]->SetConvergence(false);
      if( (config_container[iZone]->GetKind_Solver() == ADJ_RANS) || (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS) )
        integration_container[iZone][ADJTURB_SOL]->SetConvergence(false);
      break;
    }
  }

}

void CHBDriver::SetHarmonicBalance(unsigned short iZone) {

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  unsigned short iVar, jZone, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());
  if (adjoint) {
    implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  }

  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

  /*--- Retrieve values from the config file ---*/
  su2double *U = new su2double[nVar];
  su2double *U_old = new su2double[nVar];
  su2double *Psi = new su2double[nVar];
  su2double *Psi_old = new su2double[nVar];
  su2double *Source = new su2double[nVar];
  su2double deltaU, deltaPsi;

  /*--- Compute period of oscillation ---*/
  su2double period = config_container[ZONE_0]->GetHarmonicBalance_Period();

  /*--- Non-dimensionalize the input period, if necessary.  */
  period /= config_container[ZONE_0]->GetTime_Ref();

  if (ExtIter == 0)
    ComputeHB_Operator();

  /*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

    /*--- Loop over each node in the volume mesh ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][iMGlevel]->GetnPoint(); iPoint++) {

      for (iVar = 0; iVar < nVar; iVar++) {
        Source[iVar] = 0.0;
      }

      /*--- Step across the columns ---*/
      for (jZone = 0; jZone < nZone; jZone++) {

        /*--- Retrieve solution at this node in current zone ---*/
        for (iVar = 0; iVar < nVar; iVar++) {

          if (!adjoint) {
            U[iVar] = solver_container[jZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetSolution(iVar);
            Source[iVar] += U[iVar]*D[iZone][jZone];

            if (implicit) {
              U_old[iVar] = solver_container[jZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetSolution_Old(iVar);
              deltaU = U[iVar] - U_old[iVar];
              Source[iVar] += deltaU*D[iZone][jZone];
            }

          }

          else {
            Psi[iVar] = solver_container[jZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->GetSolution(iVar);
            Source[iVar] += Psi[iVar]*D[jZone][iZone];

            if (implicit) {
              Psi_old[iVar] = solver_container[jZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->GetSolution_Old(iVar);
              deltaPsi = Psi[iVar] - Psi_old[iVar];
              Source[iVar] += deltaPsi*D[jZone][iZone];
            }
          }
        }

        /*--- Store sources for current row ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (!adjoint) {
            solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iVar]);
          }
          else {
            solver_container[iZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iVar]);
          }
        }

      }
    }
  }

  /*--- Source term for a turbulence model ---*/
  if (config_container[ZONE_0]->GetKind_Solver() == RANS) {

    /*--- Extra variables needed if we have a turbulence model. ---*/
    unsigned short nVar_Turb = solver_container[ZONE_0][MESH_0][TURB_SOL]->GetnVar();
    su2double *U_Turb = new su2double[nVar_Turb];
    su2double *Source_Turb = new su2double[nVar_Turb];

    /*--- Loop over only the finest mesh level (turbulence is always solved
     on the original grid only). ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][MESH_0]->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar_Turb; iVar++) Source_Turb[iVar] = 0.0;
      for (jZone = 0; jZone < nZone; jZone++) {

        /*--- Retrieve solution at this node in current zone ---*/
        for (iVar = 0; iVar < nVar_Turb; iVar++) {
          U_Turb[iVar] = solver_container[jZone][MESH_0][TURB_SOL]->node[iPoint]->GetSolution(iVar);
          Source_Turb[iVar] += U_Turb[iVar]*D[iZone][jZone];
        }
      }

      /*--- Store sources for current iZone ---*/
      for (iVar = 0; iVar < nVar_Turb; iVar++)
        solver_container[iZone][MESH_0][TURB_SOL]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source_Turb[iVar]);
    }

    delete [] U_Turb;
    delete [] Source_Turb;
  }

  delete [] U;
  delete [] U_old;
  delete [] Psi;
  delete [] Psi_old;

}

void CHBDriver::StabilizeHarmonicBalance() {

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  unsigned short i, j, k, iVar, iZone, jZone, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());
  if (adjoint) {
    implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  }

  /*--- Retrieve values from the config file ---*/
  su2double *Source     = new su2double[nZone];
  su2double *Source_old = new su2double[nZone];
  su2double Delta;

  su2double **Pinv     = new su2double*[nZone];
  su2double **P        = new su2double*[nZone];
  for (iZone = 0; iZone < nZone; iZone++) {
    Pinv[iZone]       = new su2double[nZone];
    P[iZone]          = new su2double[nZone];
  }

  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

    /*--- Loop over each node in the volume mesh ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][iMGlevel]->GetnPoint(); iPoint++) {

      /*--- Get time step for current node ---*/
      Delta = solver_container[ZONE_0][iMGlevel][FLOW_SOL]->node[iPoint]->GetDelta_Time();

      /*--- Setup stabilization matrix for this node ---*/
      for (iZone = 0; iZone < nZone; iZone++) {
        for (jZone = 0; jZone < nZone; jZone++) {
          if (jZone == iZone ) {
            Pinv[iZone][jZone] = 1.0 + Delta*D[iZone][jZone];
          }
          else {
            Pinv[iZone][jZone] = Delta*D[iZone][jZone];
          }
        }
      }

      /*--- Invert stabilization matrix Pinv with Gauss elimination---*/

      /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
      su2double **temp = new su2double*[nZone];
      for (i = 0; i < nZone; i++) {
        temp[i] = new su2double[2 * nZone];
      }

      /*---  Copy the desired matrix into the temporary matrix ---*/
      for (i = 0; i < nZone; i++) {
        for (j = 0; j < nZone; j++) {
          temp[i][j] = Pinv[i][j];
          temp[i][nZone + j] = 0;
        }
        temp[i][nZone + i] = 1;
      }

      su2double max_val;
      unsigned short max_idx;

      /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
      for (k = 0; k < nZone - 1; k++) {
        max_idx = k;
        max_val = abs(temp[k][k]);
        /*---  Find the largest value (pivot) in the column  ---*/
        for (j = k; j < nZone; j++) {
          if (abs(temp[j][k]) > max_val) {
            max_idx = j;
            max_val = abs(temp[j][k]);
          }
        }

        /*---  Move the row with the highest value up  ---*/
        for (j = 0; j < (nZone * 2); j++) {
          su2double d = temp[k][j];
          temp[k][j] = temp[max_idx][j];
          temp[max_idx][j] = d;
        }
        /*---  Subtract the moved row from all other rows ---*/
        for (i = k + 1; i < nZone; i++) {
          su2double c = temp[i][k] / temp[k][k];
          for (j = 0; j < (nZone * 2); j++) {
            temp[i][j] = temp[i][j] - temp[k][j] * c;
          }
        }
      }

      /*---  Back-substitution  ---*/
      for (k = nZone - 1; k > 0; k--) {
        if (temp[k][k] != su2double(0.0)) {
          for (int i = k - 1; i > -1; i--) {
            su2double c = temp[i][k] / temp[k][k];
            for (j = 0; j < (nZone * 2); j++) {
              temp[i][j] = temp[i][j] - temp[k][j] * c;
            }
          }
        }
      }

      /*---  Normalize the inverse  ---*/
      for (i = 0; i < nZone; i++) {
        su2double c = temp[i][i];
        for (j = 0; j < nZone; j++) {
          temp[i][j + nZone] = temp[i][j + nZone] / c;
        }
      }

      /*---  Copy the inverse back to the main program flow ---*/
      for (i = 0; i < nZone; i++) {
        for (j = 0; j < nZone; j++) {
          P[i][j] = temp[i][j + nZone];
        }
      }

      /*---  Delete dynamic template  ---*/
      for (iZone = 0; iZone < nZone; iZone++) {
        delete[] temp[iZone];
      }
      delete[] temp;

      /*--- Loop through variables to precondition ---*/
      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Get current source terms (not yet preconditioned) and zero source array to prepare preconditioning ---*/
        for (iZone = 0; iZone < nZone; iZone++) {
          Source_old[iZone] = solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetHarmonicBalance_Source(iVar);
          Source[iZone] = 0;
        }

        /*--- Step through columns ---*/
        for (iZone = 0; iZone < nZone; iZone++) {
          for (jZone = 0; jZone < nZone; jZone++) {
            Source[iZone] += P[iZone][jZone]*Source_old[jZone];
          }

          /*--- Store updated source terms for current node ---*/
          if (!adjoint) {
            solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iZone]);
          }
          else {
            solver_container[iZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iZone]);
          }
        }

      }
    }
  }

  /*--- Deallocate dynamic memory ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    delete [] P[iZone];
    delete [] Pinv[iZone];
  }
  delete [] P;
  delete [] Pinv;
  delete [] Source;
  delete [] Source_old;

}

void CHBDriver::ComputeHB_Operator() {

  const   complex<su2double> J(0.0,1.0);
  unsigned short i, j, k, iZone;

  su2double *Omega_HB       = new su2double[nZone];
  complex<su2double> **E    = new complex<su2double>*[nZone];
  complex<su2double> **Einv = new complex<su2double>*[nZone];
  complex<su2double> **DD   = new complex<su2double>*[nZone];
  for (iZone = 0; iZone < nZone; iZone++) {
    E[iZone]    = new complex<su2double>[nZone];
    Einv[iZone] = new complex<su2double>[nZone];
    DD[iZone]   = new complex<su2double>[nZone];
  }

  /*--- Get simualation period from config file ---*/
  su2double Period = config_container[ZONE_0]->GetHarmonicBalance_Period();

  /*--- Non-dimensionalize the input period, if necessary.      */
  Period /= config_container[ZONE_0]->GetTime_Ref();

  /*--- Build the array containing the selected frequencies to solve ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    Omega_HB[iZone]  = config_container[iZone]->GetOmega_HB()[iZone];
    Omega_HB[iZone] /= config_container[iZone]->GetOmega_Ref();
  }

  /*--- Build the diagonal matrix of the frequencies DD ---*/
  for (i = 0; i < nZone; i++) {
    for (k = 0; k < nZone; k++) {
      if (k == i ) {
        DD[i][k] = J*Omega_HB[k];
      }
    }
  }


  /*--- Build the harmonic balance inverse matrix ---*/
  for (i = 0; i < nZone; i++) {
    for (k = 0; k < nZone; k++) {
      Einv[i][k] = complex<su2double>(cos(Omega_HB[k]*(i*Period/nZone))) + J*complex<su2double>(sin(Omega_HB[k]*(i*Period/nZone)));
    }
  }

  /*---  Invert inverse harmonic balance Einv with Gauss elimination ---*/

  /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
  complex<su2double> **temp = new complex<su2double>*[nZone];
  for (i = 0; i < nZone; i++) {
    temp[i] = new complex<su2double>[2 * nZone];
  }

  /*---  Copy the desired matrix into the temporary matrix ---*/
  for (i = 0; i < nZone; i++) {
    for (j = 0; j < nZone; j++) {
      temp[i][j] = Einv[i][j];
      temp[i][nZone + j] = 0;
    }
    temp[i][nZone + i] = 1;
  }

  su2double max_val;
  unsigned short max_idx;

  /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
  for (k = 0; k < nZone - 1; k++) {
    max_idx = k;
    max_val = abs(temp[k][k]);
    /*---  Find the largest value (pivot) in the column  ---*/
    for (j = k; j < nZone; j++) {
      if (abs(temp[j][k]) > max_val) {
        max_idx = j;
        max_val = abs(temp[j][k]);
      }
    }
    /*---  Move the row with the highest value up  ---*/
    for (j = 0; j < (nZone * 2); j++) {
      complex<su2double> d = temp[k][j];
      temp[k][j] = temp[max_idx][j];
      temp[max_idx][j] = d;
    }
    /*---  Subtract the moved row from all other rows ---*/
    for (i = k + 1; i < nZone; i++) {
      complex<su2double> c = temp[i][k] / temp[k][k];
      for (j = 0; j < (nZone * 2); j++) {
        temp[i][j] = temp[i][j] - temp[k][j] * c;
      }
    }
  }
  /*---  Back-substitution  ---*/
  for (k = nZone - 1; k > 0; k--) {
    if (temp[k][k] != complex<su2double>(0.0)) {
      for (int i = k - 1; i > -1; i--) {
        complex<su2double> c = temp[i][k] / temp[k][k];
        for (j = 0; j < (nZone * 2); j++) {
          temp[i][j] = temp[i][j] - temp[k][j] * c;
        }
      }
    }
  }
  /*---  Normalize the inverse  ---*/
  for (i = 0; i < nZone; i++) {
    complex<su2double> c = temp[i][i];
    for (j = 0; j < nZone; j++) {
      temp[i][j + nZone] = temp[i][j + nZone] / c;
    }
  }
  /*---  Copy the inverse back to the main program flow ---*/
  for (i = 0; i < nZone; i++) {
    for (j = 0; j < nZone; j++) {
      E[i][j] = temp[i][j + nZone];
    }
  }
  /*---  Delete dynamic template  ---*/
  for (i = 0; i < nZone; i++) {
    delete[] temp[i];
  }
  delete[] temp;


  /*---  Temporary matrix for performing product  ---*/
  complex<su2double> **Temp    = new complex<su2double>*[nZone];

  /*---  Temporary complex HB operator  ---*/
  complex<su2double> **Dcpx    = new complex<su2double>*[nZone];

  for (iZone = 0; iZone < nZone; iZone++){
    Temp[iZone]    = new complex<su2double>[nZone];
    Dcpx[iZone]   = new complex<su2double>[nZone];
  }


  /*---  Calculation of the HB operator matrix ---*/
  for (int row = 0; row < nZone; row++) {
    for (int col = 0; col < nZone; col++) {
      for (int inner = 0; inner < nZone; inner++) {
        Temp[row][col] += Einv[row][inner] * DD[inner][col];
      }
    }
  }

  unsigned short row, col, inner;

  for (row = 0; row < nZone; row++) {
    for (col = 0; col < nZone; col++) {
      for (inner = 0; inner < nZone; inner++) {
        Dcpx[row][col] += Temp[row][inner] * E[inner][col];
      }
    }
  }

  /*---  Take just the real part of the HB operator matrix ---*/
  for (i = 0; i < nZone; i++) {
    for (k = 0; k < nZone; k++) {
      D[i][k] = real(Dcpx[i][k]);
    }
  }

  /*--- Deallocate dynamic memory ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    delete [] E[iZone];
    delete [] Einv[iZone];
    delete [] DD[iZone];
    delete [] Temp[iZone];
    delete [] Dcpx[iZone];
  }
  delete [] E;
  delete [] Einv;
  delete [] DD;
  delete [] Temp;
  delete [] Dcpx;
  delete [] Omega_HB;

}

CFSIDriver::CFSIDriver(char* confFile,
                       unsigned short val_nZone,
                       unsigned short val_nDim,
                       SU2_Comm MPICommunicator) : CDriver(confFile,
                                                          val_nZone,
                                                          val_nDim,
                                                          MPICommunicator) { }

CFSIDriver::~CFSIDriver(void) { }

void CFSIDriver::Run() {

  /*--- As of now, we are coding it for just 2 zones. ---*/
  /*--- This will become more general, but we need to modify the configuration for that ---*/
  unsigned short ZONE_FLOW = 0, ZONE_STRUCT = 1;
  unsigned short iZone;

  /*--- Boolean to determine if we are running a static or dynamic case ---*/
  bool stat_fsi = ((config_container[ZONE_FLOW]->GetUnsteady_Simulation() == STEADY) && (config_container[ZONE_STRUCT]->GetDynamic_Analysis() == STATIC));
  bool dyn_fsi = (((config_container[ZONE_FLOW]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[ZONE_FLOW]->GetUnsteady_Simulation() == DT_STEPPING_2ND))
                   && (config_container[ZONE_STRUCT]->GetDynamic_Analysis() == DYNAMIC));

  unsigned long IntIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetIntIter(IntIter);
  unsigned long FSIIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(FSIIter);
  unsigned long nFSIIter = config_container[ZONE_FLOW]->GetnIterFSI();
  unsigned long nIntIter;

  bool StopCalc_Flow = false;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Be careful with whether or not we load the coords and grid velocity
   from the restart files... this needs to be standardized for the different
   solvers, in particular with FSI. ---*/

  // TODO: Test this
  bool update_geo = false;
  //if (config->GetFSI_Simulation()) update_geo = false;

  /*--- If there is a restart, we need to get the old geometry from the fluid field ---*/
  bool restart = (config_container[ZONE_FLOW]->GetRestart() || config_container[ZONE_FLOW]->GetRestart_Flow());
  ExtIter = config_container[ZONE_FLOW]->GetExtIter();

  if (restart && dyn_fsi && (long)ExtIter == config_container[ZONE_FLOW]->GetUnst_RestartIter()) {
    solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Restart_OldGeometry(geometry_container[ZONE_FLOW][MESH_0],config_container[ZONE_FLOW]);
  } else if (restart && stat_fsi){
    solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[ZONE_FLOW], solver_container[ZONE_FLOW], config_container[ZONE_FLOW], 0, update_geo);
  }

  /*-----------------------------------------------------------------*/
  /*---------------- Predict structural displacements ---------------*/
  /*-----------------------------------------------------------------*/

  Predict_Displacements(ZONE_STRUCT, ZONE_FLOW);

  while (FSIIter < nFSIIter) {

    /*-----------------------------------------------------------------*/
    /*------------------- Transfer Displacements ----------------------*/
    /*-----------------------------------------------------------------*/
  if(transfer_container[ZONE_STRUCT][ZONE_FLOW] != NULL)
      Transfer_Displacements(ZONE_STRUCT, ZONE_FLOW);

    /*-----------------------------------------------------------------*/
    /*--------------------- Mesh deformation --------------------------*/
    /*-----------------------------------------------------------------*/

  iteration_container[ZONE_FLOW]->SetGrid_Movement(geometry_container,surface_movement, grid_movement, FFDBox, solver_container,
        config_container, ZONE_FLOW, 0, ExtIter);

    /*-----------------------------------------------------------------*/
    /*-------------------- Fluid subiteration -------------------------*/
    /*-----------------------------------------------------------------*/

  iteration_container[ZONE_FLOW]->Preprocess(output, integration_container, geometry_container,
      solver_container, numerics_container, config_container,
      surface_movement, grid_movement, FFDBox, ZONE_FLOW);

  if ( stat_fsi ) {

    /*--- For steady-state flow simulations, we need to loop over ExtIter for the number of time steps ---*/
    /*--- However, ExtIter is the number of FSI iterations, so nIntIter is used in this case ---*/

    nIntIter = config_container[ZONE_FLOW]->GetUnst_nIntIter();

    for (IntIter = 0; IntIter < nIntIter; IntIter++){

      /*--- Set ExtIter to iExtIter_FLOW; this is a trick to loop on the steady-state flow solver ---*/
      config_container[ZONE_FLOW]->SetExtIter(IntIter);

      iteration_container[ZONE_FLOW]->Iterate(output, integration_container, geometry_container,
          solver_container, numerics_container, config_container,
          surface_movement, grid_movement, FFDBox, ZONE_FLOW);

      /*--- Write the convergence history for the fluid (only screen output) ---*/

      output->SetConvHistory_Body(&ConvHist_file[ZONE_0], geometry_container, solver_container, config_container, integration_container, false, 0.0, ZONE_FLOW);

      /*--- If the convergence criteria is met for the flow, break the loop ---*/
      StopCalc_Flow = integration_container[ZONE_FLOW][FLOW_SOL]->GetConvergence();
      if (StopCalc_Flow) break;

    }

  }
  else if ( dyn_fsi ) {

    /*--- For unsteady flow simulations, we need to loop over nIntIter for the number of time steps ---*/

    nIntIter = config_container[ZONE_FLOW]->GetUnst_nIntIter();

    for (IntIter = 0; IntIter < nIntIter; IntIter++){

      config_container[ZONE_FLOW]->SetIntIter(IntIter);

      iteration_container[ZONE_FLOW]->Iterate(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_FLOW);

      /*--- If convergence was reached in every zone --*/

      if (integration_container[ZONE_FLOW][FLOW_SOL]->GetConvergence() == 1) break;
    }

    /*--- Write the convergence history for the fluid (only screen output) ---*/

     output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_FLOW);

  } else {

    if (rank == MASTER_NODE) cout << "The definition of Fluid and Structural solvers is inconsistent for FSI applications " << endl;

    exit(EXIT_FAILURE);

  }

  /*--- Set the fluid convergence to false (to make sure FSI subiterations converge) ---*/

  integration_container[ZONE_FLOW][FLOW_SOL]->SetConvergence(false);

  /*-----------------------------------------------------------------*/
  /*------------------- Set FEA loads from fluid --------------------*/
  /*-----------------------------------------------------------------*/
  if(transfer_container[ZONE_FLOW][ZONE_STRUCT] != NULL)
      Transfer_Tractions(ZONE_FLOW, ZONE_STRUCT);

    /*-----------------------------------------------------------------*/
    /*------------------ Structural subiteration ----------------------*/
    /*-----------------------------------------------------------------*/

  iteration_container[ZONE_STRUCT]->Iterate(output, integration_container, geometry_container,
                                  solver_container, numerics_container, config_container,
                                  surface_movement, grid_movement, FFDBox, ZONE_STRUCT);

    /*--- Write the convergence history for the structure (only screen output) ---*/

    output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUCT);

    /*--- Set the fluid convergence to false (to make sure FSI subiterations converge) ---*/

    integration_container[ZONE_STRUCT][FEA_SOL]->SetConvergence(false);

    /*-----------------------------------------------------------------*/
    /*----------------- Displacements relaxation ----------------------*/
    /*-----------------------------------------------------------------*/

    Relaxation_Displacements(ZONE_STRUCT, ZONE_FLOW, FSIIter);

    /*-----------------------------------------------------------------*/
    /*-------------------- Check convergence --------------------------*/
    /*-----------------------------------------------------------------*/

    integration_container[ZONE_STRUCT][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT], solver_container[ZONE_STRUCT][MESH_0][FEA_SOL], FSIIter);

    if (integration_container[ZONE_STRUCT][FEA_SOL]->GetConvergence_FSI()) break;

    /*-----------------------------------------------------------------*/
    /*--------------------- Update FSIIter ---------------------------*/
    /*-----------------------------------------------------------------*/

    FSIIter++; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(FSIIter);

  }

  /*-----------------------------------------------------------------*/
  /*------------------ Update coupled solver ------------------------*/
  /*-----------------------------------------------------------------*/

  Update(ZONE_FLOW, ZONE_STRUCT);

  /*-----------------------------------------------------------------*/
  /*-------------------- Update fluid solver ------------------------*/
  /*-----------------------------------------------------------------*/

  iteration_container[ZONE_FLOW]->Update(output, integration_container, geometry_container,
                       solver_container, numerics_container, config_container,
                       surface_movement, grid_movement, FFDBox, ZONE_FLOW);

  /*-----------------------------------------------------------------*/
  /*----------------- Update structural solver ----------------------*/
  /*-----------------------------------------------------------------*/

  iteration_container[ZONE_STRUCT]->Update(output, integration_container, geometry_container,
                         solver_container, numerics_container, config_container,
                         surface_movement, grid_movement, FFDBox, ZONE_STRUCT);


  /*-----------------------------------------------------------------*/
  /*--------------- Update convergence parameter --------------------*/
  /*-----------------------------------------------------------------*/
  integration_container[ZONE_STRUCT][FEA_SOL]->SetConvergence_FSI(false);

}

void CFSIDriver::Predict_Displacements(unsigned short donorZone, unsigned short targetZone) {

#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  solver_container[donorZone][MESH_0][FEA_SOL]->PredictStruct_Displacement(geometry_container[donorZone], config_container[donorZone],
      solver_container[donorZone]);

  /*--- For parallel simulations we need to communicate the predicted solution before updating the fluid mesh ---*/

  solver_container[donorZone][MESH_0][FEA_SOL]->Set_MPI_Solution_Pred(geometry_container[donorZone][MESH_0], config_container[donorZone]);
  

}

void CFSIDriver::Predict_Tractions(unsigned short donorZone, unsigned short targetZone) {

}

void CFSIDriver::Transfer_Displacements(unsigned short donorZone, unsigned short targetZone) {

#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  bool MatchingMesh = config_container[targetZone]->GetMatchingMesh();

  /*--- Select the transfer method and the appropriate mesh properties (matching or nonmatching mesh) ---*/

  switch (config_container[targetZone]->GetKind_TransferMethod()) {
  case BROADCAST_DATA:
    if (MatchingMesh) {
        transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
                                                                                    geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
                                                                                    config_container[donorZone], config_container[targetZone]);
      /*--- Set the volume deformation for the fluid zone ---*/
      //      grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);
        
      }
      else {
        transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
                                                                                       geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
                                                                                       config_container[donorZone], config_container[targetZone]);
      /*--- Set the volume deformation for the fluid zone ---*/
      //      grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);
    }
    break;
  case SCATTER_DATA:
    if (MatchingMesh) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
                                                                         geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
                                                                         config_container[donorZone], config_container[targetZone]);
      /*--- Set the volume deformation for the fluid zone ---*/
      //      grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);
      }
      else {
        cout << "Scatter method not implemented for non-matching meshes. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    break;
  case ALLGATHER_DATA:
    if (MatchingMesh) {
        cout << "Allgather method not yet implemented for matching meshes. Exiting..." << endl;
      exit(EXIT_FAILURE);
      }
      else {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
                                                                           geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
                                                                           config_container[donorZone], config_container[targetZone]);
      /*--- Set the volume deformation for the fluid zone ---*/
      //      grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);
    }
    break;
  case LEGACY_METHOD:
    if (MatchingMesh) {
        solver_container[targetZone][MESH_0][FLOW_SOL]->SetFlow_Displacement(geometry_container[targetZone], grid_movement[targetZone],
          config_container[targetZone], config_container[donorZone],
          geometry_container[donorZone], solver_container[donorZone]);
      }
      else {
        solver_container[targetZone][MESH_0][FLOW_SOL]->SetFlow_Displacement_Int(geometry_container[targetZone], grid_movement[targetZone],
          config_container[targetZone], config_container[donorZone],
          geometry_container[donorZone], solver_container[donorZone]);
    }
    break;
  }

}

void CFSIDriver::Transfer_Tractions(unsigned short donorZone, unsigned short targetZone) {

#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  bool MatchingMesh = config_container[donorZone]->GetMatchingMesh();

  /*--- Load transfer --  This will have to be modified for non-matching meshes ---*/

  unsigned short SolContainer_Position_fea = config_container[targetZone]->GetContainerPosition(RUNTIME_FEA_SYS);

  /*--- FEA equations -- Necessary as the SetFEA_Load routine is as of now contained in the structural solver ---*/
  unsigned long ExtIter = config_container[targetZone]->GetExtIter();
  config_container[targetZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

  /*--- Select the transfer method and the appropriate mesh properties (matching or nonmatching mesh) ---*/

  switch (config_container[donorZone]->GetKind_TransferMethod()) {
  case BROADCAST_DATA:
    if (MatchingMesh) {
        transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FEA_SOL],
                                                                                    geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
                                                                                    config_container[donorZone], config_container[targetZone]);
      }
      else {
        transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FEA_SOL],
                                                                                       geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
                                                                                       config_container[donorZone], config_container[targetZone]);
    }
    break;
  case SCATTER_DATA:
    if (MatchingMesh) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FEA_SOL],
                                                                         geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
                                                                         config_container[donorZone], config_container[targetZone]);
      }
      else {
        cout << "Scatter method not implemented for non-matching meshes. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    break;
  case ALLGATHER_DATA:
    if (MatchingMesh) {
        cout << "Allgather method not yet implemented for matching meshes. Exiting..." << endl;
      exit(EXIT_FAILURE);
      }
      else {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FEA_SOL],
                                                                           geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
                                                                           config_container[donorZone], config_container[targetZone]);
    }
    break;
  case LEGACY_METHOD:
    if (MatchingMesh) {
        solver_container[targetZone][MESH_0][FEA_SOL]->SetFEA_Load(solver_container[donorZone], geometry_container[targetZone], geometry_container[donorZone],
                                                                   config_container[targetZone], config_container[donorZone], numerics_container[targetZone][MESH_0][SolContainer_Position_fea][FEA_TERM]);
      }
      else {
        solver_container[targetZone][MESH_0][FEA_SOL]->SetFEA_Load_Int(solver_container[donorZone], geometry_container[targetZone], geometry_container[donorZone],
                                                                       config_container[targetZone], config_container[donorZone], numerics_container[targetZone][MESH_0][SolContainer_Position_fea][FEA_TERM]);
    }
    break;
  }

}

void CFSIDriver::Relaxation_Displacements(unsigned short donorZone, unsigned short targetZone, unsigned long FSIIter) {

#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*-------------------- Aitken's relaxation ------------------------*/

  /*------------------- Compute the coefficient ---------------------*/

  solver_container[donorZone][MESH_0][FEA_SOL]->ComputeAitken_Coefficient(geometry_container[donorZone], config_container[donorZone],
      solver_container[donorZone], FSIIter);

  /*----------------- Set the relaxation parameter ------------------*/

  solver_container[donorZone][MESH_0][FEA_SOL]->SetAitken_Relaxation(geometry_container[donorZone], config_container[donorZone],
      solver_container[donorZone]);

  /*----------------- Communicate the predicted solution and the old one ------------------*/
  solver_container[donorZone][MESH_0][FEA_SOL]->Set_MPI_Solution_Pred_Old(geometry_container[donorZone][MESH_0], config_container[donorZone]);
  

}

void CFSIDriver::Relaxation_Tractions(unsigned short donorZone, unsigned short targetZone, unsigned long FSIIter) {

}

void CFSIDriver::Update(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT) {

  unsigned long IntIter = 0; // This doesn't affect here but has to go into the function
  ExtIter = config_container[ZONE_FLOW]->GetExtIter();

  /*-----------------------------------------------------------------*/
  /*--------------------- Enforce continuity ------------------------*/
  /*-----------------------------------------------------------------*/

  /*--- Enforces that the geometry of the flow corresponds to the converged, relaxed solution ---*/

  /*-------------------- Transfer the displacements --------------------*/

  Transfer_Displacements(ZONE_STRUCT, ZONE_FLOW);

  /*-------------------- Set the grid movement -------------------------*/

  iteration_container[ZONE_FLOW]->SetGrid_Movement(geometry_container, surface_movement,
                                                   grid_movement, FFDBox, solver_container,
      config_container, ZONE_FLOW, IntIter, ExtIter);

  /*--- TODO: Temporary output of objective function for Flow OFs. Needs to be integrated into the refurbished output ---*/

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE){

  /*--- Choose the filename of the objective function ---*/

    ofstream myfile_res;
    bool of_output = false;
    su2double objective_function = 0.0;

    switch (config_container[ZONE_FLOW]->GetKind_ObjFunc()) {
      case DRAG_COEFFICIENT:
        myfile_res.open("of_drag.opt");
        objective_function = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CD();
        of_output = true;
        break;
      case LIFT_COEFFICIENT:
        myfile_res.open("of_lift.opt");
        objective_function = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CL();
        of_output = true;
      break;
      case EFFICIENCY:
        myfile_res.open("of_efficiency.opt");
        objective_function = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CEff();
        of_output = true;
        break;
      default:
        of_output = false;
        break;
    }

    if (of_output){

        myfile_res.precision(15);
        myfile_res << scientific << objective_function << endl;
        myfile_res.close();

    }

  }


}


CFSIStatDriver::CFSIStatDriver(char* confFile,
                                   unsigned short val_nZone,
                                   unsigned short val_nDim,
                                   SU2_Comm MPICommunicator) : CFSIDriver(confFile,
                                                                           val_nZone,
                                                                           val_nDim,
                                                                           MPICommunicator) { }

CFSIStatDriver::~CFSIStatDriver(void) { }

void CFSIStatDriver::Run() {

  /*--- As of now, we are coding it for just 2 zones. ---*/
  /*--- This will become more general, but we need to modify the configuration for that ---*/
  unsigned short ZONE_FLOW = 0, ZONE_STRUCT = 1;
  unsigned short iZone;

  int rank = MASTER_NODE;

  unsigned long IntIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetIntIter(IntIter);
  unsigned long FSIIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(FSIIter);
  unsigned long nFSIIter = config_container[ZONE_FLOW]->GetnIterFSI();
  bool update_geo = false;  //TODO: check

  bool StopCalc_Flow = false;

  /*--- For steady flow cases, the loop is in nExtIter - For static FSI nExtIter has to be 1, so we need to define an inner loop. ---*/
  /*--- I will use GetUnst_nIntIter() temporarily, as an analogy with the dynamic solver ---*/
  unsigned long nExtIter_FLOW = config_container[ZONE_FLOW]->GetUnst_nIntIter();
  unsigned long iExtIter_FLOW = 0;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  ofstream ConvHist_file;
  if (rank == MASTER_NODE)
  output->SetConvHistory_Header(&ConvHist_file, config_container[ZONE_0], ZONE_0);

   /*--- If there is a restart, we need to get the old geometry from the fluid field ---*/
   bool restart = (config_container[ZONE_FLOW]->GetRestart() || config_container[ZONE_FLOW]->GetRestart_Flow());
//   unsigned long ExtIter = config_container[ZONE_FLOW]->GetExtIter();
   ExtIter = config_container[ZONE_FLOW]->GetExtIter();

   if (restart){
    unsigned short ZONE_FLOW = 0;
    solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[ZONE_FLOW], solver_container[ZONE_FLOW], config_container[ZONE_FLOW], 0, update_geo);
   }

  /*-----------------------------------------------------------------*/
  /*---------------- Predict structural displacements ---------------*/
  /*-----------------------------------------------------------------*/

  Predict_Displacements(ZONE_STRUCT, ZONE_FLOW);


  while (FSIIter < nFSIIter){

    /*-----------------------------------------------------------------*/
    /*------------------- Transfer Displacements ----------------------*/
    /*-----------------------------------------------------------------*/

    Transfer_Displacements(ZONE_STRUCT, ZONE_FLOW);

    /*-----------------------------------------------------------------*/
    /*------------------- Set the Grid movement -----------------------*/
    /*---- No longer done in the preprocess of the flow iteration -----*/
    /*---- as the flag Grid_Movement is set to false in this case -----*/
    /*-----------------------------------------------------------------*/

    iteration_container[ZONE_FLOW]->SetGrid_Movement(geometry_container, surface_movement,
                                                     grid_movement, FFDBox, solver_container,
                                                     config_container, ZONE_FLOW, IntIter, ExtIter);

    /*-----------------------------------------------------------------*/
    /*-------------------- Fluid subiteration -------------------------*/
    /*---- Unsteady flows loop over the ExtIter: this loop needs to ---*/
    /*------ be moved here as the nExtIter is 1 (FSI iterations) ------*/
    /*-----------------------------------------------------------------*/

    for (iExtIter_FLOW = 0; iExtIter_FLOW < nExtIter_FLOW; iExtIter_FLOW++){

      /*--- Set ExtIter to iExtIter_FLOW; this is a trick to loop on the steady-state flow solver ---*/

      config_container[ZONE_FLOW]->SetExtIter(iExtIter_FLOW);

      /*--- For now only preprocess and iterate are necessary ---*/

      iteration_container[ZONE_FLOW]->Preprocess(output, integration_container, geometry_container,
                                             solver_container, numerics_container, config_container,
                                             surface_movement, grid_movement, FFDBox, ZONE_FLOW);

      iteration_container[ZONE_FLOW]->Iterate(output, integration_container, geometry_container,
                                             solver_container, numerics_container, config_container,
                                             surface_movement, grid_movement, FFDBox, ZONE_FLOW);

      /*--- Write the convergence history for the fluid (only screen output) ---*/
      /*--- test what to do for steady-state screen-only output ---*/

      output->SetConvHistory_Body(&ConvHist_file, geometry_container, solver_container, config_container, integration_container, false, 0.0, ZONE_FLOW);

      switch (config_container[ZONE_FLOW]->GetKind_Solver()) {
        case EULER: case NAVIER_STOKES: case RANS:
          StopCalc_Flow = integration_container[ZONE_0][FLOW_SOL]->GetConvergence(); break;
      }

      /*--- If the convergence criteria is met for the flow, break the loop ---*/
      if (StopCalc_Flow) break;

    }

    /*--- Set the fluid convergence to false (to make sure FSI subiterations converge) ---*/

    integration_container[ZONE_FLOW][FLOW_SOL]->SetConvergence(false);

    /*-----------------------------------------------------------------*/
    /*------------------- Set FEA loads from fluid --------------------*/
    /*-----------------------------------------------------------------*/

    Transfer_Tractions(ZONE_FLOW, ZONE_STRUCT);

    /*-----------------------------------------------------------------*/
    /*------------------ Structural subiteration ----------------------*/
    /*-----------------------------------------------------------------*/

    iteration_container[ZONE_STRUCT]->Iterate(output, integration_container, geometry_container,
                                           solver_container, numerics_container, config_container,
                                           surface_movement, grid_movement, FFDBox, ZONE_STRUCT);

    /*--- Write the convergence history for the structure (only screen output) ---*/

    output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUCT);

    /*--- Set the fluid convergence to false (to make sure FSI subiterations converge) ---*/

    integration_container[ZONE_STRUCT][FEA_SOL]->SetConvergence(false);

    /*-----------------------------------------------------------------*/
    /*----------------- Displacements relaxation ----------------------*/
    /*-----------------------------------------------------------------*/

    Relaxation_Displacements(ZONE_STRUCT, ZONE_FLOW, FSIIter);

    /*-----------------------------------------------------------------*/
    /*-------------------- Check convergence --------------------------*/
    /*-----------------------------------------------------------------*/

    integration_container[ZONE_STRUCT][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT],
                                solver_container[ZONE_STRUCT][MESH_0][FEA_SOL], FSIIter);

    if (integration_container[ZONE_STRUCT][FEA_SOL]->GetConvergence_FSI()) break;

    /*-----------------------------------------------------------------*/
    /*--------------------- Update FSIIter ---------------------------*/
    /*-----------------------------------------------------------------*/

    FSIIter++; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(FSIIter);

  }

  /*-----------------------------------------------------------------*/
  /*------------------ Update coupled solver ------------------------*/
  /*-----------------------------------------------------------------*/

  Update(ZONE_FLOW, ZONE_STRUCT);


  /*-----------------------------------------------------------------*/
  /*----------------- Update structural solver ----------------------*/
  /*-----------------------------------------------------------------*/

  /*--- Output the relaxed result, which is the one transferred into the fluid domain (for restart purposes) ---*/
  switch (config_container[ZONE_STRUCT]->GetKind_TimeIntScheme_FEA()) {
    case (NEWMARK_IMPLICIT):
      solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->ImplicitNewmark_Relaxation(geometry_container[ZONE_STRUCT][MESH_0], solver_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT]);
      break;
  }

  /*-----------------------------------------------------------------*/
  /*--------------- Update convergence parameter --------------------*/
  /*-----------------------------------------------------------------*/
  integration_container[ZONE_STRUCT][FEA_SOL]->SetConvergence_FSI(false);


}




CDiscAdjFSIStatDriver::CDiscAdjFSIStatDriver(char* confFile,
                                                   unsigned short val_nZone,
                                                   unsigned short val_nDim,
                                                   SU2_Comm MPICommunicator) : CFSIStatDriver(confFile,
                                                                                               val_nZone,
                                                                                               val_nDim,
                                                                                               MPICommunicator) { 

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  unsigned short iVar;
  unsigned short nVar_Flow = 0, nVar_Struct = 0;
  RecordingState = 0;
  CurrentRecording = 0;

  switch (config_container[ZONE_0]->GetKind_ObjFunc()){
  case DRAG_COEFFICIENT:
  case LIFT_COEFFICIENT:
  case SIDEFORCE_COEFFICIENT:
  case EFFICIENCY:
  case MOMENT_X_COEFFICIENT:
  case MOMENT_Y_COEFFICIENT:
  case MOMENT_Z_COEFFICIENT:
  case EQUIVALENT_AREA:
    Kind_Objective_Function = FLOW_OBJECTIVE_FUNCTION;
    break;
  case REFERENCE_GEOMETRY:
  case REFERENCE_NODE:
    Kind_Objective_Function = FEM_OBJECTIVE_FUNCTION;
    break;
  default:
    Kind_Objective_Function = NO_OBJECTIVE_FUNCTION;
    break;
  }

  direct_iteration = new CIteration*[nZone];

  unsigned short iZone;
  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_Solver()) {
       case DISC_ADJ_RANS: case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES:
         direct_iteration[iZone] = new CFluidIteration(config_container[iZone]);
         nVar_Flow = solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetnVar();
         flow_criteria = config_container[iZone]->GetMinLogResidual_BGS_F();
         flow_criteria_rel = config_container[iZone]->GetOrderMagResidual_BGS_F();
         break;
       case DISC_ADJ_FEM:
         direct_iteration[iZone] = new CFEM_StructuralAnalysis(config_container[iZone]);
         nVar_Struct = solver_container[iZone][MESH_0][ADJFEA_SOL]->GetnVar();
         structure_criteria    = config_container[iZone]->GetMinLogResidual_BGS_S();
         structure_criteria_rel = config_container[iZone]->GetOrderMagResidual_BGS_S();
         break;
    }
  }

  init_res_flow   = new su2double[nVar_Flow];
  init_res_struct = new su2double[nVar_Struct];

  residual_flow   = new su2double[nVar_Flow];
  residual_struct = new su2double[nVar_Struct];

  residual_flow_rel   = new su2double[nVar_Flow];
  residual_struct_rel = new su2double[nVar_Struct];

  for (iVar = 0; iVar < nVar_Flow; iVar++){
    init_res_flow[iVar] = 0.0;
    residual_flow[iVar] = 0.0;
    residual_flow_rel[iVar] = 0.0;
  }
  for (iVar = 0; iVar < nVar_Struct; iVar++){
    init_res_struct[iVar] = 0.0;
    residual_struct[iVar] = 0.0;
    residual_struct_rel[iVar] = 0.0;
  }


  bool write_history = true;

  /*--- Header of the temporary output file ---*/
  if ((write_history) && (rank == MASTER_NODE)){
    ofstream myfile_res;
    myfile_res.open ("history_adjoint_FSI.csv");

    myfile_res << "BGS_Iter\t";

    for (iVar = 0; iVar < nVar_Flow; iVar++){
      myfile_res << "ResFlow[" << iVar << "]\t";
    }

    for (iVar = 0; iVar < nVar_Struct; iVar++){
      myfile_res << "ResFEA[" << iVar << "]\t";
    }


    bool de_effects = config_container[ZONE_0]->GetDE_Effects();
    for (iVar = 0; iVar < config_container[ZONE_0]->GetnElasticityMod(); iVar++)
        myfile_res << "Sens_E_" << iVar << "\t";

    for (iVar = 0; iVar < config_container[ZONE_0]->GetnPoissonRatio(); iVar++)
      myfile_res << "Sens_Nu_" << iVar << "\t";

    if (de_effects){
        for (iVar = 0; iVar < config_container[ZONE_0]->GetnElectric_Field(); iVar++)
          myfile_res << "Sens_EField_" << iVar << "\t";
    }

    myfile_res << endl;

    myfile_res.close();
  }

  // TEST: for implementation of python framework in standalone structural problems
  if ((config_container[ZONE_1]->GetDV_FEA() != NODV_FEA) && (rank == MASTER_NODE)){

    /*--- Header of the temporary output file ---*/
    ofstream myfile_res;

    switch (config_container[ZONE_1]->GetDV_FEA()) {
      case YOUNG_MODULUS:
        myfile_res.open("grad_young.opt");
        break;
      case POISSON_RATIO:
        myfile_res.open("grad_poisson.opt");
        break;
      case DENSITY_VAL:
      case DEAD_WEIGHT:
        myfile_res.open("grad_density.opt");
        break;
      case ELECTRIC_FIELD:
        myfile_res.open("grad_efield.opt");
        break;
      default:
        myfile_res.open("grad.opt");
        break;
    }

    unsigned short iDV;
    unsigned short nDV = solver_container[ZONE_1][MESH_0][ADJFEA_SOL]->GetnDVFEA();

    myfile_res << "INDEX" << "\t" << "GRAD" << endl;

    myfile_res.precision(15);

    for (iDV = 0; iDV < nDV; iDV++){
      myfile_res << iDV;
      myfile_res << "\t";
      myfile_res << scientific << solver_container[ZONE_1][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_DVFEA(iDV);
      myfile_res << endl;
    }

    myfile_res.close();
  }

  /*--- TODO: This is a workaround until the TestCases.py script incorporates new classes for nested loops. ---*/
  config_container[ZONE_0]->SetnExtIter(1);
  config_container[ZONE_1]->SetnExtIter(1);

}

CDiscAdjFSIStatDriver::~CDiscAdjFSIStatDriver(void) {

  delete [] direct_iteration;
  delete [] init_res_flow;
  delete [] init_res_struct;
  delete [] residual_flow;
  delete [] residual_struct;
  delete [] residual_flow_rel;
  delete [] residual_struct_rel;

}


void CDiscAdjFSIStatDriver::Run( ) {

  /*--- As of now, we are coding it for just 2 zones. ---*/
  /*--- This will become more general, but we need to modify the configuration for that ---*/
  unsigned short ZONE_FLOW = 0, ZONE_STRUCT = 1;
  unsigned short iZone;

  unsigned long IntIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetIntIter(IntIter);
  unsigned long FSIIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(FSIIter);

  Preprocess(ZONE_FLOW, ZONE_STRUCT, ALL_VARIABLES);

  switch (Kind_Objective_Function){
  case FLOW_OBJECTIVE_FUNCTION:
    Iterate_Block_FlowOF(ZONE_FLOW, ZONE_STRUCT, ALL_VARIABLES);
    break;
  case FEM_OBJECTIVE_FUNCTION:
    Iterate_Block_StructuralOF(ZONE_FLOW, ZONE_STRUCT, ALL_VARIABLES);
    break;
  }

  Postprocess(ZONE_FLOW, ZONE_STRUCT);

}


void CDiscAdjFSIStatDriver::Preprocess(unsigned short ZONE_FLOW,
                  unsigned short ZONE_STRUCT,
                  unsigned short kind_recording){

  unsigned long IntIter = 0, iPoint;
  config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned short ExtIter = config_container[ZONE_FLOW]->GetExtIter();

  bool dual_time_1st = (config_container[ZONE_FLOW]->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container[ZONE_FLOW]->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  unsigned short iMesh;
  int Direct_Iter_Flow;
  bool update_geo = false;

  /*----------------------------------------------------------------------------*/
  /*------------------------------ FLOW SOLUTION -------------------------------*/
  /*----------------------------------------------------------------------------*/

  /*--- For the unsteady adjoint, load direct solutions from restart files. ---*/

  if (config_container[ZONE_FLOW]->GetUnsteady_Simulation()) {

    Direct_Iter_Flow = SU2_TYPE::Int(config_container[ZONE_FLOW]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 2;

    /*--- For dual-time stepping we want to load the already converged solution at timestep n ---*/

    if (dual_time) {
      Direct_Iter_Flow += 1;
    }

    if (ExtIter == 0){

      if (dual_time_2nd) {

        /*--- Load solution at timestep n-2 ---*/

        iteration_container[ZONE_FLOW]->LoadUnsteady_Solution(geometry_container, solver_container,config_container, ZONE_FLOW, Direct_Iter_Flow-2);

        /*--- Push solution back to correct array ---*/

        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][iMesh]->GetnPoint();iPoint++) {
            solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
            solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
            if (turbulent) {
              solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
              solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
            }
          }
        }
      }
      if (dual_time) {

        /*--- Load solution at timestep n-1 ---*/

        iteration_container[ZONE_FLOW]->LoadUnsteady_Solution(geometry_container, solver_container,config_container, ZONE_FLOW, Direct_Iter_Flow-1);

        /*--- Push solution back to correct array ---*/

        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][iMesh]->GetnPoint();iPoint++) {
            solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
            if (turbulent) {
              solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
            }
          }
        }
      }

      /*--- Load solution timestep n ---*/

      iteration_container[ZONE_FLOW]->LoadUnsteady_Solution(geometry_container, solver_container,config_container, ZONE_FLOW, Direct_Iter_Flow);

    }


    if ((ExtIter > 0) && dual_time){

      /*--- Load solution timestep n - 2 ---*/

      iteration_container[ZONE_FLOW]->LoadUnsteady_Solution(geometry_container, solver_container,config_container, ZONE_FLOW, Direct_Iter_Flow - 2);

      /*--- Temporarily store the loaded solution in the Solution_Old array ---*/

      for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
        for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][iMesh]->GetnPoint();iPoint++) {
           solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->Set_OldSolution();
           if (turbulent){
             solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->Set_OldSolution();
           }
        }
      }

      /*--- Set Solution at timestep n to solution at n-1 ---*/

      for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
        for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][iMesh]->GetnPoint();iPoint++) {
          solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->GetSolution_time_n());
          if (turbulent) {
            solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->GetSolution_time_n());
          }
        }
      }
      if (dual_time_1st){
      /*--- Set Solution at timestep n-1 to the previously loaded solution ---*/
        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][iMesh]->GetnPoint();iPoint++) {
            solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->GetSolution_time_n1());
            if (turbulent) {
              solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n(solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->GetSolution_time_n1());
            }
          }
        }
      }
      if (dual_time_2nd){
        /*--- Set Solution at timestep n-1 to solution at n-2 ---*/
        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][iMesh]->GetnPoint();iPoint++) {
            solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->GetSolution_time_n1());
            if (turbulent) {
              solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n(solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->GetSolution_time_n1());
            }
          }
        }
        /*--- Set Solution at timestep n-2 to the previously loaded solution ---*/
        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][iMesh]->GetnPoint();iPoint++) {
            solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->GetSolution_Old());
            if (turbulent) {
              solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1(solver_container[ZONE_FLOW][iMesh][TURB_SOL]->node[iPoint]->GetSolution_Old());
            }
          }
        }
      }
    }
  }
  else{

    /*--- Load the restart (we need to use the routine in order to get the GEOMETRY, otherwise it's restarted from the base mesh ---*/

    solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[ZONE_FLOW], solver_container[ZONE_FLOW], config_container[ZONE_FLOW], 0, true);

  if (ExtIter == 0 || dual_time) {

    for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
      for (iPoint = 0; iPoint < geometry_container[ZONE_FLOW][iMesh]->GetnPoint(); iPoint++) {
        solver_container[ZONE_FLOW][iMesh][ADJFLOW_SOL]->node[iPoint]->SetSolution_Direct(solver_container[ZONE_FLOW][iMesh][FLOW_SOL]->node[iPoint]->GetSolution());
      }
    }
    if (turbulent && !config_container[ZONE_FLOW]->GetFrozen_Visc_Disc()) {
      for (iPoint = 0; iPoint < geometry_container[ZONE_FLOW][MESH_0]->GetnPoint(); iPoint++) {
        solver_container[ZONE_FLOW][MESH_0][ADJTURB_SOL]->node[iPoint]->SetSolution_Direct(solver_container[ZONE_FLOW][MESH_0][TURB_SOL]->node[iPoint]->GetSolution());
      }
    }
  }

    /*--- Store geometry of the converged solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry_container[ZONE_FLOW][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->node[iPoint]->SetGeometry_Direct(geometry_container[ZONE_FLOW][MESH_0]->node[iPoint]->GetCoord());
    }

  }

  /*----------------------------------------------------------------------------*/
  /*-------------------------- STRUCTURAL SOLUTION -----------------------------*/
  /*----------------------------------------------------------------------------*/

  IntIter = 0;
  config_container[ZONE_STRUCT]->SetIntIter(IntIter);
  ExtIter = config_container[ZONE_STRUCT]->GetExtIter();
  bool dynamic = (config_container[ZONE_STRUCT]->GetDynamic_Analysis() == DYNAMIC);

  int Direct_Iter_FEA;

  /*--- For the dynamic adjoint, load direct solutions from restart files. ---*/

  if (dynamic) {

    Direct_Iter_FEA = SU2_TYPE::Int(config_container[ZONE_STRUCT]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 1;

    /*--- We want to load the already converged solution at timesteps n and n-1 ---*/

    /*--- Load solution at timestep n-1 ---*/

    iteration_container[ZONE_STRUCT]->LoadDynamic_Solution(geometry_container, solver_container,config_container, ZONE_STRUCT, Direct_Iter_FEA-1);

    /*--- Push solution back to correct array ---*/

    for(iPoint=0; iPoint<geometry_container[ZONE_STRUCT][MESH_0]->GetnPoint();iPoint++){
      solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->node[iPoint]->SetSolution_time_n();
    }

    /*--- Push solution back to correct array ---*/

    for(iPoint=0; iPoint<geometry_container[ZONE_STRUCT][MESH_0]->GetnPoint();iPoint++){
      solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Accel_time_n();
    }

    /*--- Push solution back to correct array ---*/

    for(iPoint=0; iPoint<geometry_container[ZONE_STRUCT][MESH_0]->GetnPoint();iPoint++){
      solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Vel_time_n();
    }

    /*--- Load solution timestep n ---*/

    iteration_container[ZONE_STRUCT]->LoadDynamic_Solution(geometry_container, solver_container,config_container, ZONE_STRUCT, Direct_Iter_FEA);

    /*--- Store FEA solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry_container[ZONE_STRUCT][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->node[iPoint]->SetSolution_Direct(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->node[iPoint]->GetSolution());
    }

    for (iPoint = 0; iPoint < geometry_container[ZONE_STRUCT][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->node[iPoint]->SetSolution_Accel_Direct(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Accel());
    }

    for (iPoint = 0; iPoint < geometry_container[ZONE_STRUCT][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->node[iPoint]->SetSolution_Vel_Direct(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel());
    }

  }
  else {

    solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->LoadRestart(geometry_container[ZONE_STRUCT], solver_container[ZONE_STRUCT], config_container[ZONE_STRUCT], 0, update_geo);

    /*--- Store FEA solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry_container[ZONE_STRUCT][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->node[iPoint]->SetSolution_Direct(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->node[iPoint]->GetSolution());
    }

  }

  /*----------------------------------------------------------------------------*/
  /*--------------------- ADJOINT SOLVER PREPROCESSING -------------------------*/
  /*----------------------------------------------------------------------------*/

  solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->Preprocessing(geometry_container[ZONE_FLOW][MESH_0], solver_container[ZONE_FLOW][MESH_0],  config_container[ZONE_FLOW] , MESH_0, 0, RUNTIME_ADJFLOW_SYS, false);

  if (turbulent){
    solver_container[ZONE_FLOW][MESH_0][ADJTURB_SOL]->Preprocessing(geometry_container[ZONE_FLOW][MESH_0], solver_container[ZONE_FLOW][MESH_0],  config_container[ZONE_FLOW] , MESH_0, 0, RUNTIME_ADJTURB_SYS, false);
  }

  solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->Preprocessing(geometry_container[ZONE_STRUCT][MESH_0], solver_container[ZONE_STRUCT][MESH_0],  config_container[ZONE_STRUCT] , MESH_0, 0, RUNTIME_ADJFEA_SYS, false);



}

void CDiscAdjFSIStatDriver::PrintDirect_Residuals(unsigned short ZONE_FLOW,
                                                          unsigned short ZONE_STRUCT,
                                                          unsigned short kind_recording){

  unsigned short ExtIter = config_container[ZONE_FLOW]->GetExtIter();
  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS);
  bool nonlinear_analysis = (config_container[ZONE_STRUCT]->GetGeometricConditions() == LARGE_DEFORMATIONS);   // Nonlinear analysis.
  bool unsteady = config_container[ZONE_FLOW]->GetUnsteady_Simulation() != NONE;
  bool dynamic = (config_container[ZONE_STRUCT]->GetDynamic_Analysis() == DYNAMIC);

  su2double val_OFunction = 0.0;
  string kind_OFunction;

  cout.precision(6);
  cout.setf(ios::scientific, ios::floatfield);

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if ((kind_recording == FLOW_CONS_VARS) || (kind_recording == MESH_COORDS)) {

    /*--- Print residuals in the first iteration ---*/

    if (rank == MASTER_NODE && ((ExtIter == 0) || unsteady )){
      cout << "log10[RMS Density]: "<< log10(solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetRes_RMS(0))
                     <<", Drag: " <<solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CD()
                     <<", Lift: " << solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CL() << "." << endl;

      if (turbulent){
        cout << "log10[RMS k]: " << log10(solver_container[ZONE_FLOW][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
      }
      if (Kind_Objective_Function == FLOW_OBJECTIVE_FUNCTION){
        switch (config_container[ZONE_FLOW]->GetKind_ObjFunc()){
        case DRAG_COEFFICIENT:
          kind_OFunction = "(Drag coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CD();
          break;
        case LIFT_COEFFICIENT:
          kind_OFunction = "(Lift coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CL();
          break;
        case SIDEFORCE_COEFFICIENT:
          kind_OFunction = "(Sideforce coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CSF();
          break;
        case EFFICIENCY:
          kind_OFunction = "(Efficiency): ";
          val_OFunction = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CEff();
          break;
        case MOMENT_X_COEFFICIENT:
          kind_OFunction = "(Moment X coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CMx();
          break;
        case MOMENT_Y_COEFFICIENT:
          kind_OFunction = "(Moment Y coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CMy();
          break;
        case MOMENT_Z_COEFFICIENT:
          kind_OFunction = "(Moment Z coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CMz();
          break;
        case EQUIVALENT_AREA:
          kind_OFunction = "(Equivalent area): ";
          val_OFunction = solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->GetTotal_CEquivArea();
          break;
        default:
          val_OFunction = 0.0;  // If the objective function is computed in a different physical problem
          break;
        }
        cout << "Objective function " << kind_OFunction << val_OFunction << endl;
      }
    }

  }

  if ((kind_recording == FEA_DISP_VARS) || (kind_recording == FLOW_CROSS_TERM) || (kind_recording == GEOMETRY_CROSS_TERM)) {

    if (rank == MASTER_NODE && ((ExtIter == 0) || dynamic )){
      if (nonlinear_analysis){
        cout << "UTOL-A: "   << log10(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_FEM(0))
             << ", RTOL-A: " << log10(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_FEM(1))
             << ", ETOL-A: " << log10(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_FEM(2)) << "." << endl;
      }
      else{
        if (solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetnVar() == 2){
          cout << "log10[RMS Ux]: "   << log10(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_RMS(0))
               << ", log10[RMS Uy]: " << log10(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_RMS(1)) << "." << endl;

        }
        else{
          cout << "log10[RMS Ux]: "   << log10(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_RMS(0))
               << ", log10[RMS Uy]: " << log10(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_RMS(1))
               << ", log10[RMS Uz]: " << log10(solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_RMS(2))<< "." << endl;
        }

      }
      if (Kind_Objective_Function == FEM_OBJECTIVE_FUNCTION){
        switch (config_container[ZONE_STRUCT]->GetKind_ObjFunc()){
        case REFERENCE_GEOMETRY:
          kind_OFunction = "(Reference Geometry): ";
          val_OFunction = solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetTotal_OFRefGeom();
          break;
        case REFERENCE_NODE:
          kind_OFunction = "(Reference Node): ";
          val_OFunction = solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetTotal_OFRefNode();
          break;
        default:
          val_OFunction = 0.0;  // If the objective function is computed in a different physical problem
          break;
        }
        cout << "Objective function " << kind_OFunction << val_OFunction << endl;
      }
    }

  }

}

void CDiscAdjFSIStatDriver::Iterate_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT, unsigned short kind_recording){

  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == MESH_COORDS)) {

    Fluid_Iteration_Direct(ZONE_FLOW, ZONE_STRUCT);


  }

  if ((kind_recording == FEA_DISP_VARS) ||
      (kind_recording == FLOW_CROSS_TERM) ||
      (kind_recording == GEOMETRY_CROSS_TERM)) {

    Structural_Iteration_Direct(ZONE_FLOW, ZONE_STRUCT);

  }


  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {

    Mesh_Deformation_Direct(ZONE_FLOW, ZONE_STRUCT);

  }


}

void CDiscAdjFSIStatDriver::Fluid_Iteration_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT) {

  /*-----------------------------------------------------------------*/
  /*------------------- Set Dependency on Geometry ------------------*/
  /*-----------------------------------------------------------------*/

  geometry_container[ZONE_FLOW][MESH_0]->UpdateGeometry(geometry_container[ZONE_FLOW], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[ZONE_FLOW][MESH_0],solver_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

  /*-----------------------------------------------------------------*/
  /*----------------- Iterate the flow solver -----------------------*/
  /*---- Sets all the cross dependencies for the flow variables -----*/
  /*-----------------------------------------------------------------*/

  config_container[ZONE_FLOW]->SetIntIter(0);

  direct_iteration[ZONE_FLOW]->Iterate(output, integration_container, geometry_container,
      solver_container, numerics_container, config_container,
      surface_movement, grid_movement, FFDBox, ZONE_FLOW);

  /*-----------------------------------------------------------------*/
  /*--------------------- Set MPI Solution --------------------------*/
  /*-----------------------------------------------------------------*/

  solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW]);

}

void CDiscAdjFSIStatDriver::Structural_Iteration_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT) {


  /*-----------------------------------------------------------------*/
  /*---------- Set Dependencies on Geometry and Flow ----------------*/
  /*-----------------------------------------------------------------*/

  solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->Set_MPI_Solution(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT]);

  geometry_container[ZONE_FLOW][MESH_0]->UpdateGeometry(geometry_container[ZONE_FLOW], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[ZONE_FLOW][MESH_0],solver_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

  /*-----------------------------------------------------------------*/
  /*-------------------- Transfer Tractions -------------------------*/
  /*-----------------------------------------------------------------*/

  Transfer_Tractions(ZONE_FLOW, ZONE_STRUCT);

  /*-----------------------------------------------------------------*/
  /*--------------- Iterate the structural solver -------------------*/
  /*-----------------------------------------------------------------*/

  direct_iteration[ZONE_STRUCT]->Iterate(output, integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox, ZONE_STRUCT);

  /*-----------------------------------------------------------------*/
  /*--------------------- Set MPI Solution --------------------------*/
  /*-----------------------------------------------------------------*/

  solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->Set_MPI_Solution(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT]);


}

void CDiscAdjFSIStatDriver::Mesh_Deformation_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT) {

  unsigned long IntIter = config_container[ZONE_STRUCT]->GetIntIter();
  unsigned long ExtIter = config_container[ZONE_STRUCT]->GetExtIter();

  /*-----------------------------------------------------------------*/
  /*--------------------- Set MPI Solution --------------------------*/
  /*-----------------------------------------------------------------*/

  geometry_container[ZONE_FLOW][MESH_0]->UpdateGeometry(geometry_container[ZONE_FLOW], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[ZONE_FLOW][MESH_0],solver_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

  solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->Set_MPI_Solution(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT]);

  /*-----------------------------------------------------------------*/
  /*------------------- Transfer Displacements ----------------------*/
  /*-----------------------------------------------------------------*/

  Transfer_Displacements(ZONE_STRUCT, ZONE_FLOW);

  /*-----------------------------------------------------------------*/
  /*------------------- Set the Grid movement -----------------------*/
  /*---- No longer done in the preprocess of the flow iteration -----*/
  /*---- as the flag Grid_Movement is set to false in this case -----*/
  /*-----------------------------------------------------------------*/

  direct_iteration[ZONE_FLOW]->SetGrid_Movement(geometry_container, surface_movement,
                                                   grid_movement, FFDBox, solver_container,
                                                   config_container, ZONE_FLOW, IntIter, ExtIter);

  geometry_container[ZONE_FLOW][MESH_0]->UpdateGeometry(geometry_container[ZONE_FLOW], config_container[ZONE_FLOW]);
  solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->Set_MPI_Solution(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT]);


}

void CDiscAdjFSIStatDriver::SetRecording(unsigned short ZONE_FLOW,
                                              unsigned short ZONE_STRUCT,
                                              unsigned short kind_recording){

  unsigned long IntIter = config_container[ZONE_0]->GetIntIter();
  bool unsteady = (config_container[ZONE_FLOW]->GetUnsteady_Simulation() != NONE);
  bool dynamic = (config_container[ZONE_STRUCT]->GetDynamic_Analysis() == DYNAMIC);

  string kind_DirectIteration = " ";
  string kind_AdjointIteration = " ";

  if (unsteady || dynamic){
    cout << "DYNAMIC ADJOINT SOLVER NOT IMPLEMENTED FOR FSI APPLICATIONS" << endl;
    exit(EXIT_FAILURE);
  }

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE){
    cout << endl;
    switch (kind_recording){
    case FLOW_CONS_VARS:
      kind_AdjointIteration = "Flow iteration: flow input -> flow output";
      kind_DirectIteration = "flow ";
      break;
    case MESH_COORDS:
      kind_AdjointIteration = "Geometry cross term from flow: geometry input -> flow output";
      kind_DirectIteration = "flow ";
      break;
    case FEA_DISP_VARS:
      kind_AdjointIteration = "Structural iteration: structural input -> structural output";
      kind_DirectIteration = "structural ";
      break;
    case FLOW_CROSS_TERM:
      kind_AdjointIteration = "Flow cross term: flow input -> structural output";
      kind_DirectIteration = "structural ";
      break;
    case GEOMETRY_CROSS_TERM:
      kind_AdjointIteration = "Geometry cross term from structure: geometry input -> structural output";
      kind_DirectIteration = "structural ";
      break;
    case FEM_CROSS_TERM_GEOMETRY:
      kind_AdjointIteration = "Structural cross term from geometry: structural input -> geometry output";
      kind_DirectIteration = "mesh deformation ";
      break;
    }
    cout << kind_AdjointIteration << endl;
    cout << "Direct " << kind_DirectIteration << "iteration to store computational graph." << endl;
    switch (kind_recording){
    case FLOW_CONS_VARS: case MESH_COORDS:
    case FEA_DISP_VARS: case FLOW_CROSS_TERM: case GEOMETRY_CROSS_TERM:
      cout << "Compute residuals to check the convergence of the direct problem." << endl; break;
    case FEM_CROSS_TERM_GEOMETRY:
      cout << "Deform the grid using the converged solution of the direct problem." << endl; break;
    }
  }


  AD::Reset();

  if (CurrentRecording != kind_recording && (CurrentRecording != NONE) ){

    /*--- Clear indices ---*/

    PrepareRecording(ZONE_FLOW, ZONE_STRUCT, ALL_VARIABLES);

    /*--- Clear indices of coupling variables ---*/

    SetDependencies(ZONE_FLOW, ZONE_STRUCT, ALL_VARIABLES);

    /*--- Run one iteration while tape is passive - this clears all indices ---*/
    Iterate_Direct(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  }

  /*--- Prepare for recording ---*/

  PrepareRecording(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Start the recording of all operations ---*/

  AD::StartRecording();

  /*--- Register input variables ---*/

  RegisterInput(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Set dependencies for flow, geometry and structural solvers ---*/

  SetDependencies(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Run a direct iteration ---*/
  Iterate_Direct(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Register objective function and output variables ---*/

  RegisterOutput(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Stop the recording ---*/
  AD::StopRecording();

  /*--- Set the recording status ---*/

  CurrentRecording = kind_recording;

  /* --- Reset the number of the internal iterations---*/

  config_container[ZONE_0]->SetIntIter(IntIter);


}

void CDiscAdjFSIStatDriver::PrepareRecording(unsigned short ZONE_FLOW,
                                                   unsigned short ZONE_STRUCT,
                                                   unsigned short kind_recording){

  unsigned short iMesh;
  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS);

  /*--- Set fluid variables to direct solver values ---*/
  for (iMesh = 0; iMesh <= config_container[ZONE_FLOW]->GetnMGLevels(); iMesh++){
    solver_container[ZONE_FLOW][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW]);
  }
  if (turbulent){
    solver_container[ZONE_FLOW][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW]);
  }

  /*--- Set geometry to the converged values ---*/

  solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->SetMesh_Recording(geometry_container[ZONE_FLOW], grid_movement[ZONE_FLOW], config_container[ZONE_FLOW]);

  /*--- Set structural variables to direct solver values ---*/

  solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->SetRecording(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT]);

}

void CDiscAdjFSIStatDriver::RegisterInput(unsigned short ZONE_FLOW,
                                               unsigned short ZONE_STRUCT,
                                               unsigned short kind_recording){

  /*--- Register flow variables ---*/
  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == FLOW_CROSS_TERM)) {
    iteration_container[ZONE_FLOW]->RegisterInput(solver_container, geometry_container, config_container, ZONE_FLOW, kind_recording);
  }

  /*--- Register geometry variables ---*/
  if ((kind_recording == MESH_COORDS) ||
      (kind_recording == GEOMETRY_CROSS_TERM)) {
    iteration_container[ZONE_FLOW]->RegisterInput(solver_container, geometry_container, config_container, ZONE_FLOW, kind_recording);
  }

  /*--- Register structural variables ---*/
  if ((kind_recording == FEA_DISP_VARS) ||
      (kind_recording == FEM_CROSS_TERM_GEOMETRY)) {
    iteration_container[ZONE_STRUCT]->RegisterInput(solver_container, geometry_container, config_container, ZONE_STRUCT, kind_recording);
  }


}

void CDiscAdjFSIStatDriver::SetDependencies(unsigned short ZONE_FLOW,
                                                  unsigned short ZONE_STRUCT,
                                                  unsigned short kind_recording){

  /*--- Add dependencies for geometrical and turbulent variables ---*/

  iteration_container[ZONE_FLOW]->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_FLOW, kind_recording);

  /*--- Add dependencies for E, Nu, Rho, and Rho_DL variables ---*/

  iteration_container[ZONE_STRUCT]->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_STRUCT, kind_recording);


}

void CDiscAdjFSIStatDriver::RegisterOutput(unsigned short ZONE_FLOW,
                                                 unsigned short ZONE_STRUCT,
                                                 unsigned short kind_recording){

  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS);

  /*--- Register the objective function as output of the iteration ---*/
  /*--- We need to avoid recording it twice for the crossed terms  ---*/

  /*--- Register a flow-type objective function ---*/
  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == MESH_COORDS)) {
    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[ZONE_FLOW]);
  }

  /*--- The FEM_CROSS_TERM_GEOMETRY evaluates the mesh routines:
   *--- They don't throw any dependency on the objective function ---*/

  /*--- Register a structural-type objective function ---*/
  if ((kind_recording == FEA_DISP_VARS) ||
      (kind_recording == FLOW_CROSS_TERM) ||
      (kind_recording == GEOMETRY_CROSS_TERM)){
    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->RegisterObj_Func(config_container[ZONE_STRUCT]);
  }

  /*--- Register the conservative variables of the flow as output of the iteration ---*/
  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == MESH_COORDS)) {

    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[ZONE_FLOW][MESH_0],config_container[ZONE_FLOW]);

    if (turbulent){
      solver_container[ZONE_FLOW][MESH_0][ADJTURB_SOL]->RegisterOutput(geometry_container[ZONE_FLOW][MESH_0],
          config_container[ZONE_FLOW]);
    }
  }

  /*--- Register the displacements of the nodes of the fluid as output of the iteration ---*/
  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {

    geometry_container[ZONE_FLOW][MESH_0]->RegisterOutput_Coordinates(config_container[ZONE_FLOW]);

  }

  /*--- Register the displacements of the structure as output of the iteration ---*/
  if ((kind_recording == FEA_DISP_VARS) ||
      (kind_recording == FLOW_CROSS_TERM) ||
      (kind_recording == GEOMETRY_CROSS_TERM)) {

    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->RegisterOutput(geometry_container[ZONE_STRUCT][MESH_0],config_container[ZONE_STRUCT]);

  }


}


void CDiscAdjFSIStatDriver::Iterate_Block(unsigned short ZONE_FLOW,
                                                unsigned short ZONE_STRUCT,
                                                unsigned short kind_recording){

  unsigned long IntIter=0, nIntIter = 1;
  bool dual_time_1st = (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  bool dynamic = (config_container[ZONE_STRUCT]->GetDynamic_Analysis() == DYNAMIC);

  bool adjoint_convergence = false;

  /*--- Record one direct iteration with kind_recording as input ---*/

  SetRecording(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Print the residuals of the direct subiteration ---*/

  PrintDirect_Residuals(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Run the iteration ---*/

  switch (kind_recording){
  case FLOW_CONS_VARS:
    nIntIter = config_container[ZONE_FLOW]->GetUnst_nIntIter();
    break;
  case FEA_DISP_VARS:
    nIntIter = config_container[ZONE_STRUCT]->GetDyn_nIntIter();
    break;
  case MESH_COORDS:
  case FEM_CROSS_TERM_GEOMETRY:
  case FLOW_CROSS_TERM:
  case GEOMETRY_CROSS_TERM:
    nIntIter = 1;
    break;
  }

  for (unsigned short iZone = 0; iZone < config_container[ZONE_FLOW]->GetnZone(); iZone++)
    config_container[iZone]->SetIntIter(IntIter);

  for(IntIter = 0; IntIter < nIntIter; IntIter++){

    /*--- Set the internal iteration ---*/

    for (unsigned short iZone = 0; iZone < config_container[ZONE_FLOW]->GetnZone(); iZone++)
      config_container[iZone]->SetIntIter(IntIter);

    /*--- Set the adjoint values of the flow and objective function ---*/

    InitializeAdjoint(ZONE_FLOW, ZONE_STRUCT, kind_recording);

    /*--- Run the adjoint computation ---*/

    AD::ComputeAdjoint();

    /*--- Extract the adjoints of the input variables and store them for the next iteration ---*/

    ExtractAdjoint(ZONE_FLOW, ZONE_STRUCT, kind_recording);

    /*--- Clear all adjoints to re-use the stored computational graph in the next iteration ---*/
    AD::ClearAdjoints();

    /*--- Check the convergence of the adjoint block ---*/

    adjoint_convergence = CheckConvergence(IntIter, ZONE_FLOW, ZONE_STRUCT, kind_recording);

    /*--- Write the convergence history (only screen output) ---*/

    ConvergenceHistory(IntIter, nIntIter, ZONE_FLOW, ZONE_STRUCT, kind_recording);

    /*--- Break the loop if converged ---*/

    if (adjoint_convergence) break;


  }

  if (dual_time){
    integration_container[ZONE_FLOW][ADJFLOW_SOL]->SetConvergence(false);
  }
  if (dynamic){
    integration_container[ZONE_FLOW][ADJFLOW_SOL]->SetConvergence(false);
  }

}


void CDiscAdjFSIStatDriver::InitializeAdjoint(unsigned short ZONE_FLOW,
                                                     unsigned short ZONE_STRUCT,
                                                     unsigned short kind_recording){

  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS);

  /*--- Initialize the adjoint of the objective function (typically with 1.0) ---*/
  /*--- We need to avoid setting it twice for the crossed terms  ---*/

  /*--- Register a flow-type objective function ---*/
  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == MESH_COORDS)) {
    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[ZONE_FLOW][MESH_0], config_container[ZONE_FLOW]);
  }

  /*--- Register a structural-type objective function ---*/
  if ((kind_recording == FEA_DISP_VARS) ||
      (kind_recording == FLOW_CROSS_TERM) ||
      (kind_recording == GEOMETRY_CROSS_TERM)){
    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->SetAdj_ObjFunc(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT]);
  }

  /*--- Adjoint of the fluid conservative variables ---*/

  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == MESH_COORDS)) {

    /*--- Initialize the adjoints the conservative variables ---*/
    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[ZONE_FLOW][MESH_0],
                                                                    config_container[ZONE_FLOW]);

    if (turbulent){
      solver_container[ZONE_FLOW][MESH_0][ADJTURB_SOL]->SetAdjoint_Output(geometry_container[ZONE_FLOW][MESH_0],
                                                                      config_container[ZONE_FLOW]);
    }

  }

  /*--- Adjoint of the positions of the mesh ---*/
  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {

    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->SetAdjoint_OutputMesh(geometry_container[ZONE_FLOW][MESH_0],
                                                                            config_container[ZONE_FLOW]);

  }

  /*--- Adjoint of the structural displacements ---*/
  if ((kind_recording == FEA_DISP_VARS) ||
      (kind_recording == FLOW_CROSS_TERM) ||
      (kind_recording == GEOMETRY_CROSS_TERM)) {

    /*--- Initialize the adjoints the conservative variables ---*/

    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->SetAdjoint_Output(geometry_container[ZONE_STRUCT][MESH_0],
                                                                         config_container[ZONE_STRUCT]);

  }

}

void CDiscAdjFSIStatDriver::ExtractAdjoint(unsigned short ZONE_FLOW,
                                                  unsigned short ZONE_STRUCT,
                                                  unsigned short kind_recording){
													  
  /*--- Extract the adjoint of the fluid conservative variables ---*/

  if (kind_recording == FLOW_CONS_VARS) {

    /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry_container[ZONE_FLOW][MESH_0],
                                                      config_container[ZONE_FLOW]);

    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry_container[ZONE_FLOW][MESH_0],
                                                      config_container[ZONE_FLOW]);

    if (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) {
      solver_container[ZONE_FLOW][MESH_0][ADJTURB_SOL]->ExtractAdjoint_Solution(geometry_container[ZONE_FLOW][MESH_0],
                                                        config_container[ZONE_FLOW]);
    }

  }

  /*--- Extract the adjoint of the mesh coordinates ---*/

  if (kind_recording == MESH_COORDS) {

    /*--- Extract the adjoints of the flow geometry and store them for the next iteration ---*/

    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_CrossTerm_Geometry_Flow(geometry_container[ZONE_FLOW][MESH_0],
                                                      config_container[ZONE_FLOW]);

  }

  /*--- Extract the adjoint of the structural displacements ---*/

  if (kind_recording == FEA_DISP_VARS) {

    /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Solution(geometry_container[ZONE_STRUCT][MESH_0],
                                                                               config_container[ZONE_STRUCT]);

    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Variables(geometry_container[ZONE_STRUCT][MESH_0],
                                                                                config_container[ZONE_STRUCT]);

  }

  /*--- Extract the adjoint cross term from the structural problem with respect to the flow variables ---*/
  if (kind_recording == FLOW_CROSS_TERM) {

    /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_CrossTerm(geometry_container[ZONE_FLOW][MESH_0],
                                                      config_container[ZONE_FLOW]);

    if (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) {
      solver_container[ZONE_FLOW][MESH_0][ADJTURB_SOL]->ExtractAdjoint_CrossTerm(geometry_container[ZONE_FLOW][MESH_0],
                                                        config_container[ZONE_FLOW]);
    }

  }

  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {

    /*--- Extract the adjoints of the displacements (input variables) and store them for the next iteration ---*/

    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->ExtractAdjoint_CrossTerm_Geometry(geometry_container[ZONE_STRUCT][MESH_0],
                                                                                config_container[ZONE_STRUCT]);

  }


  if (kind_recording == GEOMETRY_CROSS_TERM) {

    /*--- Extract the adjoints of the geometry input variables and store them for the next iteration ---*/

    solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_CrossTerm_Geometry(geometry_container[ZONE_FLOW][MESH_0],
                                                                                        config_container[ZONE_FLOW]);

  }

}



bool CDiscAdjFSIStatDriver::CheckConvergence(unsigned long IntIter,
                                                   unsigned short ZONE_FLOW,
                                                   unsigned short ZONE_STRUCT,
                                                   unsigned short kind_recording){

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  bool flow_convergence    = false,
        struct_convergence  = false;

  bool adjoint_convergence = false;

  su2double residual_1, residual_2;

  if (kind_recording == FLOW_CONS_VARS) {

      /*--- Set the convergence criteria (only residual possible as of now) ---*/

      residual_1 = log10(solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0));
      residual_2 = log10(solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->GetRes_RMS(1));

      flow_convergence = ((residual_1 < config_container[ZONE_FLOW]->GetMinLogResidual()) &&
                          (residual_2 < config_container[ZONE_FLOW]->GetMinLogResidual()));

  }

  if (kind_recording == FEA_DISP_VARS) {

    /*--- Set the convergence criteria (only residual possible as of now) ---*/

    residual_1 = log10(solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->GetRes_RMS(0));
    residual_2 = log10(solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->GetRes_RMS(1));

    // Temporary, until function is added
    struct_convergence = ((residual_1 < config_container[ZONE_STRUCT]->GetResidual_FEM_UTOL()) &&
                          (residual_2 < config_container[ZONE_STRUCT]->GetResidual_FEM_UTOL()));

  }

  switch (kind_recording){
  case FLOW_CONS_VARS:      adjoint_convergence = flow_convergence; break;
  case MESH_COORDS:  adjoint_convergence = true; break;
  case FEA_DISP_VARS:       adjoint_convergence = struct_convergence; break;
  case FLOW_CROSS_TERM:     adjoint_convergence = true; break;
  case FEM_CROSS_TERM_GEOMETRY:      adjoint_convergence = true; break;
  case GEOMETRY_CROSS_TERM: adjoint_convergence = true; break;
  default:                  adjoint_convergence = false; break;
  }

  /*--- Apply the same convergence criteria to all the processors ---*/

#ifdef HAVE_MPI

  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;

  /*--- Convergence criteria ---*/

  sbuf_conv[0] = adjoint_convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

  /*-- Compute global convergence criteria in the master node --*/

  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }

  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);

  if (sbuf_conv[0] == 1) { adjoint_convergence = true;}
  else { adjoint_convergence = false;}

  delete [] sbuf_conv;
  delete [] rbuf_conv;

#endif

  return adjoint_convergence;

}

void CDiscAdjFSIStatDriver::ConvergenceHistory(unsigned long IntIter,
                                                      unsigned long nIntIter,
                                                      unsigned short ZONE_FLOW,
                                                      unsigned short ZONE_STRUCT,
                                                      unsigned short kind_recording){

  unsigned long BGS_Iter = config_container[ZONE_FLOW]->GetFSIIter();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  ofstream ConvHist_file;
  if (rank == MASTER_NODE)
    output->SetConvHistory_Header(&ConvHist_file, config_container[ZONE_0], ZONE_0);


  if (kind_recording == FLOW_CONS_VARS) {

    if (rank == MASTER_NODE){
      if (IntIter == 0){
        cout << endl;
        cout << " Iter" << "    BGSIter" << "   Res[Psi_Rho]" << "     Res[Psi_E]" << endl;
      }

      if (IntIter % config_container[ZONE_FLOW]->GetWrt_Con_Freq() == 0){
        /*--- Output the flow convergence ---*/
        /*--- This is temporary as it requires several changes in the output structure ---*/
        cout.width(5);     cout << IntIter;
        cout.width(11);    cout << BGS_Iter + 1;
        cout.precision(6); cout.setf(ios::fixed, ios::floatfield);
        cout.width(15);    cout << log10(solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0));
        cout.width(15);    cout << log10(solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->GetRes_RMS(1));
        cout << endl;
      }

    }
  }

  if (kind_recording == FEA_DISP_VARS) {

    if (rank == MASTER_NODE){
      if (IntIter == 0){
        cout << endl;
        cout << " Iter" << "    BGSIter" << "    Res[Ux_bar]" << "     Res[Uy_bar]";
        cout << "       Sens_E" << "       Sens_Nu" << endl;
      }
    }

    /*--- Set the convergence criteria (only residual possible) ---*/
    output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUCT);


    /*--- Output the structural convergence ---*/
    /*--- This is temporary as it requires several changes in the output structure ---*/
//    cout.width(5);  cout << IntIter;
//    cout.width(11); cout << BGS_Iter;
//    cout.width(15); cout.precision(4); cout.setf(ios::scientific, ios::floatfield);
//    cout << solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_RMS(0);
//    cout << solver_container[ZONE_STRUCT][MESH_0][FEA_SOL]->GetRes_RMS(1);

  }


}

void CDiscAdjFSIStatDriver::Iterate_Block_FlowOF(unsigned short ZONE_FLOW,
                                                       unsigned short ZONE_STRUCT,
                                                       unsigned short kind_recording){

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  bool BGS_Converged = false;

  unsigned short iZone;

  unsigned long iFSIIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(iFSIIter);
  unsigned long nFSIIter = config_container[ZONE_FLOW]->GetnIterFSI();

  for (iFSIIter = 0; iFSIIter < nFSIIter; iFSIIter++){

    if (rank == MASTER_NODE){
      cout << endl << "                    ****** BGS ITERATION ";
      cout << iFSIIter;
      cout << " ******" << endl;
    }

    for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(iFSIIter);

    /*--- Iterate fluid (including cross term) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FLOW_CONS_VARS);

    /*--- Compute mesh (it is a cross term dF / dMv ) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, MESH_COORDS);

    /*--- Compute mesh cross term (dM / dSv) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FEM_CROSS_TERM_GEOMETRY);

    /*--- Iterate structure first ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FEA_DISP_VARS);

    /*--- Compute cross term (dS / dFv) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FLOW_CROSS_TERM);

    /*--- Compute cross term (dS / dMv) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, GEOMETRY_CROSS_TERM);


    /*--- Check convergence of the BGS method ---*/
    BGS_Converged = BGSConvergence(iFSIIter, ZONE_FLOW, ZONE_STRUCT);

    if (BGS_Converged) break;

  }


}

void CDiscAdjFSIStatDriver::Iterate_Block_StructuralOF(unsigned short ZONE_FLOW,
                                                              unsigned short ZONE_STRUCT,
                                                              unsigned short kind_recording){

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  bool BGS_Converged = false;

  unsigned short iZone;

  unsigned long iFSIIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(iFSIIter);
  unsigned long nFSIIter = config_container[ZONE_FLOW]->GetnIterFSI();

  ofstream myfile_struc, myfile_flow, myfile_geo;

  for (iFSIIter = 0; iFSIIter < nFSIIter; iFSIIter++){

    if (rank == MASTER_NODE){
      cout << endl << "                    ****** BGS ITERATION ";
      cout << iFSIIter;
      cout << " ******" << endl;
    }

    for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(iFSIIter);

    /*--- Iterate structure first ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FEA_DISP_VARS);

    /*--- Compute cross term (dS / dFv) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FLOW_CROSS_TERM);

    /*--- Compute cross term (dS / dMv) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, GEOMETRY_CROSS_TERM);

    /*--- Iterate fluid (including cross term) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FLOW_CONS_VARS);

    /*--- Compute mesh (it is a cross term dF / dMv ) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, MESH_COORDS);

    /*--- Compute mesh cross term (dM / dSv) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FEM_CROSS_TERM_GEOMETRY);


    /*--- Check convergence of the BGS method ---*/
    BGS_Converged = BGSConvergence(iFSIIter, ZONE_FLOW, ZONE_STRUCT);

    if (BGS_Converged) break;

  }


}

bool CDiscAdjFSIStatDriver::BGSConvergence(unsigned long IntIter,
                                                 unsigned short ZONE_FLOW,
                                                 unsigned short ZONE_STRUCT){

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  unsigned short iMarker;
  unsigned short nVar_Flow = solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->GetnVar(),
                   nVar_Struct = solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->GetnVar();
  unsigned short iRes;

  bool flow_converged_absolute = false,
        flow_converged_relative = false,
        struct_converged_absolute = false,
        struct_converged_relative = false;

  bool Convergence = false;

  /*--- Apply BC's to the structural adjoint - otherwise, clamped nodes have too values that make no sense... ---*/
  for (iMarker = 0; iMarker < config_container[ZONE_STRUCT]->GetnMarker_All(); iMarker++)
  switch (config_container[ZONE_STRUCT]->GetMarker_All_KindBC(iMarker)) {
    case CLAMPED_BOUNDARY:
    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->BC_Clamped_Post(geometry_container[ZONE_STRUCT][MESH_0],
        solver_container[ZONE_STRUCT][MESH_0], numerics_container[ZONE_STRUCT][MESH_0][FEA_SOL][FEA_TERM],
        config_container[ZONE_STRUCT], iMarker);
    break;
  }

  /*--- Compute the residual for the flow and structural zones ---*/

  /*--- Flow ---*/

  solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->ComputeResidual_BGS(geometry_container[ZONE_FLOW][MESH_0],
                                                                        config_container[ZONE_FLOW]);

  /*--- Structure ---*/

  solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->ComputeResidual_BGS(geometry_container[ZONE_STRUCT][MESH_0],
                                                                         config_container[ZONE_STRUCT]);


  /*--- Retrieve residuals ---*/

  /*--- Flow residuals ---*/

  for (iRes = 0; iRes < nVar_Flow; iRes++){
    residual_flow[iRes] = log10(solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->GetRes_BGS(iRes));
    if (IntIter == 0) init_res_flow[iRes] = residual_flow[iRes];
    residual_flow_rel[iRes] = fabs(residual_flow[iRes] - init_res_flow[iRes]);
  }

  /*--- Structure residuals ---*/

  for (iRes = 0; iRes < nVar_Struct; iRes++){
    residual_struct[iRes] = log10(solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->GetRes_BGS(iRes));
    if (IntIter == 0) init_res_struct[iRes] = residual_struct[iRes];
    residual_struct_rel[iRes] = fabs(residual_struct[iRes] - init_res_struct[iRes]);
  }

  /*--- Check convergence ---*/
  flow_converged_absolute = ((residual_flow[0] < flow_criteria) && (residual_flow[nVar_Flow-1] < flow_criteria));
  flow_converged_relative = ((residual_flow_rel[0] > flow_criteria_rel) && (residual_flow_rel[nVar_Flow-1] > flow_criteria_rel));

  struct_converged_absolute = ((residual_struct[0] < structure_criteria) && (residual_struct[nVar_Flow-1] < structure_criteria));
  struct_converged_relative = ((residual_struct_rel[0] > structure_criteria_rel) && (residual_struct_rel[nVar_Flow-1] > structure_criteria_rel));

  Convergence = ((flow_converged_absolute && struct_converged_absolute) ||
                 (flow_converged_absolute && struct_converged_relative) ||
                 (flow_converged_relative && struct_converged_relative) ||
                 (flow_converged_relative && struct_converged_absolute));

  if (rank == MASTER_NODE){

    cout << endl << "-------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << "Convergence summary for BGS subiteration ";
    cout << IntIter << endl;
    cout << endl;
    /*--- TODO: This is a workaround until the TestCases.py script incorporates new classes for nested loops. ---*/
    cout << "   Iter[ID]" << "    [Psi_Rho]" << "      [Psi_E]" << "     [Psi_Ux]" << "     [Psi_Uy]" << endl;
    cout.precision(6); cout.setf(ios::fixed, ios::floatfield);
    cout.width(11); cout << IntIter*1000;
    cout.width(13); cout << residual_flow[0];
    cout.width(13); cout << residual_flow[nVar_Flow-1];
    cout.width(13); cout << residual_struct[0];
    cout.width(13); cout << residual_struct[1];
    cout << endl;
    cout << endl;
    cout.precision(6); cout.setf(ios::fixed, ios::floatfield);
    cout << "                      "; cout << "   Absolute" << "     Criteria" << "     Relative" << "     Criteria" << endl;
    cout << "Flow       [Psi_Rho]: ";
    cout.width(11); cout << residual_flow[0];
    cout.width(13); cout << flow_criteria;
    cout.width(13); cout << residual_flow_rel[0];
    cout.width(13); cout << flow_criteria_rel << endl;
    cout << "             [Psi_E]: ";
    cout.width(11); cout << residual_flow[nVar_Flow-1];
    cout.width(13); cout << flow_criteria;
    cout.width(13); cout << residual_flow_rel[nVar_Flow-1];
    cout.width(13); cout << flow_criteria_rel << endl;
    cout << "Structure   [Psi_Ux]: ";
    cout.width(11); cout << residual_struct[0];
    cout.width(13); cout << structure_criteria;
    cout.width(13); cout << residual_struct_rel[0];
    cout.width(13); cout << structure_criteria_rel << endl;
    cout << "            [Psi_Uy]: ";
    cout.width(11); cout << residual_struct[1];
    cout.width(13); cout << structure_criteria;
    cout.width(13); cout << residual_struct_rel[1];
    cout.width(13); cout << structure_criteria_rel << endl;
    if (geometry_container[ZONE_FLOW][MESH_0]->GetnDim() == 3 ){
      cout << "            [Psi_Uz]: ";
      cout.width(11); cout << residual_struct[2];
      cout.width(13); cout << structure_criteria;
      cout.width(13); cout << residual_struct_rel[2];
      cout.width(13); cout << structure_criteria_rel << endl;        }
    cout << endl;
    cout << "-------------------------------------------------------------------------" << endl;


    bool write_history = true;
    unsigned short iVar;

    /*--- Header of the temporary output file ---*/
    if ((write_history) && (rank == MASTER_NODE)){
      ofstream myfile_res;
      bool de_effects = config_container[ZONE_STRUCT]->GetDE_Effects();

      myfile_res.open ("history_adjoint_FSI.csv", ios::app);

      myfile_res << IntIter << "\t";

      myfile_res.precision(15);

      for (iVar = 0; iVar < nVar_Flow; iVar++){
        myfile_res << fixed << residual_flow[iVar] << "\t";
      }

      for (iVar = 0; iVar < nVar_Struct; iVar++){
        myfile_res << fixed << residual_struct[iVar] << "\t";
      }

      for (iVar = 0; iVar < config_container[ZONE_STRUCT]->GetnElasticityMod(); iVar++)
         myfile_res << scientific << solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_E(iVar) << "\t";
      for (iVar = 0; iVar < config_container[ZONE_STRUCT]->GetnPoissonRatio(); iVar++)
         myfile_res << scientific << solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_Nu(iVar) << "\t";
      if (de_effects){
        for (iVar = 0; iVar < config_container[ZONE_STRUCT]->GetnElectric_Field(); iVar++)
          myfile_res << scientific << solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_EField(0) << "\t";
      }

      myfile_res << endl;

      myfile_res.close();
    }

    // TEST: for implementation of python framework in coupled FSI problems
    if ((config_container[ZONE_1]->GetDV_FEA() != NODV_FEA) && (rank == MASTER_NODE)){

      /*--- Header of the temporary output file ---*/
      ofstream myfile_res;

      switch (config_container[ZONE_1]->GetDV_FEA()) {
        case YOUNG_MODULUS:
          myfile_res.open("grad_young.opt");
          break;
        case POISSON_RATIO:
          myfile_res.open("grad_poisson.opt");
          break;
        case DENSITY_VAL:
        case DEAD_WEIGHT:
          myfile_res.open("grad_density.opt");
          break;
        case ELECTRIC_FIELD:
          myfile_res.open("grad_efield.opt");
          break;
        default:
          myfile_res.open("grad.opt");
          break;
      }

      unsigned short iDV;
      unsigned short nDV = solver_container[ZONE_1][MESH_0][ADJFEA_SOL]->GetnDVFEA();

      myfile_res << "INDEX" << "\t" << "GRAD" << endl;

      myfile_res.precision(15);

      for (iDV = 0; iDV < nDV; iDV++){
        myfile_res << iDV;
        myfile_res << "\t";
        myfile_res << scientific << solver_container[ZONE_1][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_DVFEA(iDV);
        myfile_res << endl;
      }

      myfile_res.close();
    }


  }

  /*--- Apply the same convergence criteria to all the processors ---*/

#ifdef HAVE_MPI

  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;

  /*--- Convergence criteria ---*/

  sbuf_conv[0] = Convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

  /*-- Compute global convergence criteria in the master node --*/

  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }

  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);

  if (sbuf_conv[0] == 1) { Convergence = true;}
  else { Convergence = false;}

  delete [] sbuf_conv;
  delete [] rbuf_conv;

#endif

  /*--- Update the solution for the flow and structural zones ---*/

  /*--- Flow ---*/

  solver_container[ZONE_FLOW][MESH_0][ADJFLOW_SOL]->UpdateSolution_BGS(geometry_container[ZONE_FLOW][MESH_0],
                                                                       config_container[ZONE_FLOW]);

  /*--- Structure ---*/

  solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->UpdateSolution_BGS(geometry_container[ZONE_STRUCT][MESH_0],
                                                                       config_container[ZONE_STRUCT]);

  return Convergence;
}

void CDiscAdjFSIStatDriver::Postprocess(unsigned short ZONE_FLOW,
                                             unsigned short ZONE_STRUCT) {

  unsigned short iMarker;

  /*--- Apply BC's to the structural adjoint after the solution has converged (to avoid unphysical values in clamped nodes) ---*/
  for (iMarker = 0; iMarker < config_container[ZONE_STRUCT]->GetnMarker_All(); iMarker++)
  switch (config_container[ZONE_STRUCT]->GetMarker_All_KindBC(iMarker)) {
    case CLAMPED_BOUNDARY:
    solver_container[ZONE_STRUCT][MESH_0][ADJFEA_SOL]->BC_Clamped_Post(geometry_container[ZONE_STRUCT][MESH_0],
        solver_container[ZONE_STRUCT][MESH_0], numerics_container[ZONE_STRUCT][MESH_0][FEA_SOL][FEA_TERM],
        config_container[ZONE_STRUCT], iMarker);
    break;
  }


}
