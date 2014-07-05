/*!
 * \file SU2_CFD.cpp
 * \brief Main file of Computational Fluid Dynamics Code (SU2_CFD).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

#include "../include/SU2_CFD.hpp"

/*
 * Pseudocode:
 *  Initialize data structures
 *  Parse the input file and load settings
 *  Load in the geometry from the mesh
 *  Set up the solver structure for the particular PDE
 *  Set up the time integration scheme
 *  Set up the numerical methods for spatial integration
 *  Instantiate the output class
 *
 *  for (ExtIter = 0; ExtIter < MaxIter; ExtIter++) {
 *    Make an iteration (problem dependent)
 *    Check Convergence of the flow solution (problem dependent)
 *    Write convergence history
 *  }
 *  Write solution files
 *  Deallocate memory
 */

using namespace std;

int main(int argc, char *argv[]) {
  
  bool StopCalc = false;
  double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  unsigned long ExtIter = 0;
  unsigned short iMesh, iZone, iSol, nZone, nDim;
  ofstream ConvHist_file;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  
#ifdef HAVE_MPI
  /*--- MPI initialization, and buffer setting ---*/
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  /*--- Create pointers to all of the classes that may be used throughout
   the SU2_CFD code. In general, the pointers are instantiated down a
   heirarchy over all zones, multigrid levels, equation sets, and equation
   terms as described in the comments below. ---*/
  
  COutput *output                       = NULL;
  CIntegration ***integration_container = NULL;
  CGeometry ***geometry_container       = NULL;
  CSolver ****solver_container          = NULL;
  CNumerics *****numerics_container     = NULL;
  CConfig **config_container            = NULL;
  CSurfaceMovement **surface_movement   = NULL;
  CVolumetricMovement **grid_movement   = NULL;
  CFreeFormDefBox*** FFDBox             = NULL;
  
  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/
  
  char config_file_name[200];
  if (argc == 2){ strcpy(config_file_name,argv[1]); }
  else{ strcpy(config_file_name, "default.cfg"); }
  
  /*--- Read the name and format of the input mesh file ---*/
  
  CConfig *config = NULL;
  config = new CConfig(config_file_name);
  
  /*--- Get the number of zones and dimensions from the numerical grid
   (required for variables allocation) ---*/
  
  nZone = GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
  nDim  = GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());
  
  /*--- Definition and of the containers for all possible zones. ---*/
  
  solver_container      = new CSolver***[nZone];
  integration_container = new CIntegration**[nZone];
  numerics_container    = new CNumerics****[nZone];
  config_container      = new CConfig*[nZone];
  geometry_container    = new CGeometry **[nZone];
  surface_movement      = new CSurfaceMovement *[nZone];
  grid_movement         = new CVolumetricMovement *[nZone];
  FFDBox                = new CFreeFormDefBox**[nZone];
  
  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone]       = NULL;
    integration_container[iZone]  = NULL;
    numerics_container[iZone]     = NULL;
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
    surface_movement[iZone]       = NULL;
    grid_movement[iZone]          = NULL;
    FFDBox[iZone]                 = NULL;
  }
  
  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
      /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/
    
    config_container[iZone] = new CConfig(config_file_name, SU2_CFD, iZone, nZone, nDim, VERB_HIGH);
    
#ifdef HAVE_MPI
    /*--- Change the name of the input-output files for a parallel computation ---*/
    config_container[iZone]->SetFileNameDomain(rank+1);
#endif
    
    /*--- Perform the non-dimensionalization for the flow equations using the
     specified reference values. ---*/
    
    config_container[iZone]->SetNondimensionalization(nDim, iZone);
    
    /*--- Definition of the geometry class. Within this constructor, the
     mesh file is read and the primal grid is stored (node coords, connectivity,
     & boundary markers. MESH_0 is the index of the finest mesh. ---*/
    
    geometry_container[iZone] = new CGeometry *[config_container[iZone]->GetMGLevels()+1];
    geometry_container[iZone][MESH_0] = new CPhysicalGeometry(config_container[iZone], iZone, nZone);
    
  }
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Geometry Preprocessing ------------------------" << endl;
  
  /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
   based data structure is constructed, i.e. node and cell neighbors are
   identified and linked, face areas and volumes of the dual mesh cells are
   computed, and the multigrid levels are created using an agglomeration procedure. ---*/
  
  Geometrical_Preprocessing(geometry_container, config_container, nZone);
  
#ifdef HAVE_MPI
  /*--- Synchronization point after the geometrical definition subroutine ---*/
MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Preprocessing --------------------------" << endl;
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Computation of wall distances for turbulence modeling ---*/
    
    if ( (config_container[iZone]->GetKind_Solver() == RANS) ||
        (config_container[iZone]->GetKind_Solver() == ADJ_RANS) )
      geometry_container[iZone][MESH_0]->ComputeWall_Distance(config_container[iZone]);
    
    /*--- Computation of positive surface area in the z-plane which is used for
     the calculation of force coefficient (non-dimensionalization). ---*/
    
    geometry_container[iZone][MESH_0]->SetPositive_ZArea(config_container[iZone]);
    
    /*--- Set the near-field, interface and actuator disk boundary conditions, if necessary. ---*/
    
    for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
      geometry_container[iZone][iMesh]->MatchNearField(config_container[iZone]);
      geometry_container[iZone][iMesh]->MatchInterface(config_container[iZone]);
      geometry_container[iZone][iMesh]->MatchActuator_Disk(config_container[iZone]);
    }
    
    /*--- Definition of the solver class: solver_container[#ZONES][#MG_GRIDS][#EQ_SYSTEMS].
     The solver classes are specific to a particular set of governing equations,
     and they contain the subroutines with instructions for computing each spatial
     term of the PDE, i.e. loops over the edges to compute convective and viscous
     fluxes, loops over the nodes to compute source terms, and routines for
     imposing various boundary condition type for the PDE. ---*/
    
    solver_container[iZone] = new CSolver** [config_container[iZone]->GetMGLevels()+1];
    for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++)
      solver_container[iZone][iMesh] = NULL;
    
    for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
      solver_container[iZone][iMesh] = new CSolver* [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        solver_container[iZone][iMesh][iSol] = NULL;
    }
    Solver_Preprocessing(solver_container[iZone], geometry_container[iZone],
                         config_container[iZone], iZone);
    
#ifdef HAVE_MPI
    /*--- Synchronization point after the solution preprocessing subroutine ---*/
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if (rank == MASTER_NODE)
      cout << endl <<"----------------- Integration and Numerics Preprocessing ----------------" << endl;
    
    /*--- Definition of the integration class: integration_container[#ZONES][#EQ_SYSTEMS].
     The integration class orchestrates the execution of the spatial integration
     subroutines contained in the solver class (including multigrid) for computing
     the residual at each node, R(U) and then integrates the equations to a
     steady state or time-accurately. ---*/
    
    integration_container[iZone] = new CIntegration*[MAX_SOLS];
    Integration_Preprocessing(integration_container[iZone], geometry_container[iZone],
                              config_container[iZone], iZone);
    
    if (rank == MASTER_NODE) cout << "Integration Preprocessing." << endl;
    
#ifdef HAVE_MPI
    /*--- Synchronization point after the integration definition subroutine ---*/
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    /*--- Definition of the numerical method class:
     numerics_container[#ZONES][#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
     The numerics class contains the implementation of the numerical methods for
     evaluating convective or viscous fluxes between any two nodes in the edge-based
     data structure (centered, upwind, galerkin), as well as any source terms
     (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/
    
    numerics_container[iZone] = new CNumerics***[config_container[iZone]->GetMGLevels()+1];
    Numerics_Preprocessing(numerics_container[iZone], solver_container[iZone],
                           geometry_container[iZone], config_container[iZone], iZone);
    
    if (rank == MASTER_NODE) cout << "Numerics Preprocessing." << endl;
    
#ifdef HAVE_MPI
    /*--- Synchronization point after the solver definition subroutine ---*/
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    /*--- Instantiate the geometry movement classes for the solution of unsteady
     flows on dynamic meshes, including rigid mesh transformations, dynamically
     deforming meshes, and time-spectral preprocessing. ---*/
    
    if (config_container[iZone]->GetGrid_Movement()) {
      if (rank == MASTER_NODE)
        cout << "Setting dynamic mesh structure." << endl;
      grid_movement[iZone] = new CVolumetricMovement(geometry_container[iZone][MESH_0]);
      FFDBox[iZone] = new CFreeFormDefBox*[MAX_NUMBER_FFD];
      surface_movement[iZone] = new CSurfaceMovement();
      surface_movement[iZone]->CopyBoundary(geometry_container[iZone][MESH_0], config_container[iZone]);
      if (config_container[iZone]->GetUnsteady_Simulation() == TIME_SPECTRAL)
        SetGrid_Movement(geometry_container[iZone], surface_movement[iZone], grid_movement[iZone],
                         FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, 0, 0);
    }
    
  }
  
  /*--- For the time-spectral solver, set the grid node velocities. ---*/
  
  if (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_SPECTRAL)
    SetTimeSpectral_Velocities(geometry_container, config_container, nZone);
  
  /*--- Coupling between zones (limited to two zones at the moment) ---*/
  
  if (nZone == 2) {
    if (rank == MASTER_NODE)
      cout << endl <<"--------------------- Setting Coupling Between Zones --------------------" << endl;
    geometry_container[ZONE_0][MESH_0]->MatchZone(config_container[ZONE_0], geometry_container[ZONE_1][MESH_0],
                                                  config_container[ZONE_1], ZONE_0, nZone);
    geometry_container[ZONE_1][MESH_0]->MatchZone(config_container[ZONE_1], geometry_container[ZONE_0][MESH_0],
                                                  config_container[ZONE_0], ZONE_1, nZone);
  }
  
  /*--- Definition of the output class (one for all zones). The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/
  
  output = new COutput();
  
  /*--- Open the convergence history file ---*/
  
  if (rank == MASTER_NODE)
    output->SetHistory_Header(&ConvHist_file, config_container[ZONE_0]);
  
  /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
  if (config_container[ZONE_0]->GetWrt_Unsteady() && config_container[ZONE_0]->GetRestart())
    ExtIter = config_container[ZONE_0]->GetUnst_RestartIter();
  
  /*--- Main external loop of the solver. Within this loop, each iteration ---*/
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;
  
  /*--- Set up a timer for performance benchmarking (preprocessing time is not included) ---*/
  
#ifndef HAVE_MPI
  StartTime = double(clock())/double(CLOCKS_PER_SEC);
#else
  MPI_Barrier(MPI_COMM_WORLD);
  StartTime = MPI_Wtime();
#endif
  
  while (ExtIter < config_container[ZONE_0]->GetnExtIter()) {
    
    /*--- Set a timer for each iteration. Store the current iteration and
     update  the value of the CFL number (if there is CFL ramping specified)
     in the config class. ---*/
    
    for (iZone = 0; iZone < nZone; iZone++) {
      config_container[iZone]->SetExtIter(ExtIter);
      config_container[iZone]->UpdateCFL(ExtIter);
    }
    
    /*--- Read the target pressure ---*/
    
    if (config_container[ZONE_0]->GetInvDesign_Cp() == YES)
      output->SetCp_InverseDesign(solver_container[ZONE_0][MESH_0][FLOW_SOL],
                                  geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
    
    /*--- Read the target heat flux ---*/

    if (config_container[ZONE_0]->GetInvDesign_HeatFlux() == YES)
      output->SetHeat_InverseDesign(solver_container[ZONE_0][MESH_0][FLOW_SOL],
                                    geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
    
    /*--- Perform a single iteration of the chosen PDE solver. ---*/
    
    switch (config_container[ZONE_0]->GetKind_Solver()) {
        
      case EULER: case NAVIER_STOKES: case RANS:
        MeanFlowIteration(output, integration_container, geometry_container,
                          solver_container, numerics_container, config_container,
                          surface_movement, grid_movement, FFDBox);
        break;
        
      case TNE2_EULER: case TNE2_NAVIER_STOKES:
        TNE2Iteration(output, integration_container,
                      geometry_container, solver_container,
                      numerics_container, config_container,
                      surface_movement, grid_movement, FFDBox);
        break;
        
      case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
        FluidStructureIteration(output, integration_container, geometry_container,
                                solver_container, numerics_container, config_container,
                                surface_movement, grid_movement, FFDBox);
        break;
        
      case WAVE_EQUATION:
        WaveIteration(output, integration_container, geometry_container,
                      solver_container, numerics_container, config_container,
                      surface_movement, grid_movement, FFDBox);
        break;
        
      case HEAT_EQUATION:
        HeatIteration(output, integration_container, geometry_container,
                      solver_container, numerics_container, config_container,
                      surface_movement, grid_movement, FFDBox);
        break;
        
      case POISSON_EQUATION:
        PoissonIteration(output, integration_container, geometry_container,
                         solver_container, numerics_container, config_container,
                         surface_movement, grid_movement, FFDBox);
        break;
        
      case LINEAR_ELASTICITY:
        FEAIteration(output, integration_container, geometry_container,
                     solver_container, numerics_container, config_container,
                     surface_movement, grid_movement, FFDBox);
        break;
        
      case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
        AdjMeanFlowIteration(output, integration_container, geometry_container,
                             solver_container, numerics_container, config_container,
                             surface_movement, grid_movement, FFDBox);
        break;
        
      case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
        AdjTNE2Iteration(output, integration_container, geometry_container,
                         solver_container, numerics_container, config_container,
                         surface_movement, grid_movement, FFDBox);
        break;
    }
    
    
    /*--- Synchronization point after a single solver iteration. Compute the
     wall clock time required. ---*/
    
#ifndef HAVE_MPI
    StopTime = double(clock())/double(CLOCKS_PER_SEC);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    StopTime = MPI_Wtime();
#endif
    
    UsedTime = (StopTime - StartTime);
    
    /*--- For specific applications, evaluate and plot the equivalent area or flow rate. ---*/
    
    if ((config_container[ZONE_0]->GetKind_Solver() == EULER) &&
        (config_container[ZONE_0]->GetEquivArea() == YES)) {
      output->SetEquivalentArea(solver_container[ZONE_0][MESH_0][FLOW_SOL],
                                geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
    }
    
    /*--- Update the convergence history file (serial and parallel computations). ---*/
    
    output->SetConvergence_History(&ConvHist_file, geometry_container, solver_container,
                                   config_container, integration_container, false, UsedTime, ZONE_0);
    
    /*--- Check whether the current simulation has reached the specified
     convergence criteria, and set StopCalc to true, if so. ---*/
    
    switch (config_container[ZONE_0]->GetKind_Solver()) {
      case EULER: case NAVIER_STOKES: case RANS:
        StopCalc = integration_container[ZONE_0][FLOW_SOL]->GetConvergence(); break;
      case TNE2_EULER: case TNE2_NAVIER_STOKES:
        StopCalc = integration_container[ZONE_0][TNE2_SOL]->GetConvergence(); break;
      case WAVE_EQUATION:
        StopCalc = integration_container[ZONE_0][WAVE_SOL]->GetConvergence(); break;
      case HEAT_EQUATION:
        StopCalc = integration_container[ZONE_0][HEAT_SOL]->GetConvergence(); break;
      case LINEAR_ELASTICITY:
        StopCalc = integration_container[ZONE_0][FEA_SOL]->GetConvergence(); break;
      case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
        StopCalc = integration_container[ZONE_0][ADJFLOW_SOL]->GetConvergence(); break;
      case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
        StopCalc = integration_container[ZONE_0][ADJTNE2_SOL]->GetConvergence(); break;
    }
    
    /*--- Solution output. Determine whether a solution needs to be written
     after the current iteration, and if so, execute the output file writing
     routines. ---*/
    
    if ((ExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
        ((ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (ExtIter != 0) &&
         !((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
           (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
        (StopCalc) ||
        (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
          (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
         ((ExtIter == 0) || (ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {
          
          /*--- Low-fidelity simulations (using a coarser multigrid level
           approximation to the solution) require an interpolation back to the
           finest grid. ---*/
          
          if (config_container[ZONE_0]->GetLowFidelitySim()) {
            integration_container[ZONE_0][FLOW_SOL]->SetProlongated_Solution(RUNTIME_FLOW_SYS, solver_container[ZONE_0][MESH_0], solver_container[ZONE_0][MESH_1], geometry_container[ZONE_0][MESH_0], geometry_container[ZONE_0][MESH_1], config_container[ZONE_0]);
            integration_container[ZONE_0][FLOW_SOL]->Smooth_Solution(RUNTIME_FLOW_SYS, solver_container[ZONE_0][MESH_0], geometry_container[ZONE_0][MESH_0], 3, 1.25, config_container[ZONE_0]);
            solver_container[ZONE_0][MESH_0][config_container[ZONE_0]->GetContainerPosition(RUNTIME_FLOW_SYS)]->Set_MPI_Solution(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
            solver_container[ZONE_0][MESH_0][config_container[ZONE_0]->GetContainerPosition(RUNTIME_FLOW_SYS)]->Preprocessing(geometry_container[ZONE_0][MESH_0], solver_container[ZONE_0][MESH_0], config_container[ZONE_0], MESH_0, 0, RUNTIME_FLOW_SYS, false);
          }
          
          /*--- Execute the routine for writing restart, volume solution,
           surface solution, and surface comma-separated value files. ---*/
          
          output->SetResult_Files(solver_container, geometry_container, config_container, ExtIter, nZone);
          
          /*--- Compute the forces at different sections. ---*/
          if (config_container[ZONE_0]->GetPlot_Section_Forces())
            output->SetForceSections(solver_container[ZONE_0][MESH_0][FLOW_SOL],
                                     geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
          
        }
    
    /*--- If the convergence criteria has been met, terminate the simulation. ---*/
    
    if (StopCalc) break;
    
    ExtIter++;
    
  }
  
  /*--- Close the convergence history file. ---*/
  
  if (rank == MASTER_NODE) {
    ConvHist_file.close();
    cout << endl <<"History file closed." << endl;
  }
  
  /*--- Numerics class deallocation ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
      for (iSol = 0; iSol < MAX_SOLS; iSol++) {
        for (unsigned short iTerm = 0; iTerm < MAX_TERMS; iTerm++){
          if (numerics_container[iZone][iMesh][iSol][iTerm] != NULL) {
            delete numerics_container[iZone][iMesh][iSol][iTerm];
          }
        }
        if (numerics_container[iZone][iMesh][iSol] != NULL)
        delete [] numerics_container[iZone][iMesh][iSol];
      }
      if (numerics_container[iZone][iMesh]!= NULL)
      delete numerics_container[iZone][iMesh];
    }
    if (numerics_container[iZone] != NULL)
    delete numerics_container[iZone];
  }
  delete [] numerics_container;
  if (rank == MASTER_NODE) cout <<"Numerics container deallocated." << endl;
  
  /*--- Solver class deallocation ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
        for (iSol = 0; iSol < MAX_SOLS; iSol++) {
          if (solver_container[iZone][iMesh][iSol] != NULL) {
            delete solver_container[iZone][iMesh][iSol];
          }
        }
        delete [] solver_container[iZone][iMesh];
      }
      delete [] solver_container[iZone];
    }
    delete [] solver_container;
    if (rank == MASTER_NODE) cout <<"Solution container, deallocated." << endl;

  /*--- Geometry class deallocation ---*/
  for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (geometry_container[iZone][iMesh]!=NULL){
        //delete geometry_container[iZone][iMesh];
      }
      //delete geometry_container[iZone];
    }
  }
  delete [] geometry_container;
  cout <<"Geometry container deallocated." << endl;
  
  /*--- Integration class deallocation ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    delete integration_container[iZone];
  }
  delete [] integration_container;
  cout <<"Integration container deallocated." << endl;
  
  /*--- Free-form deformation class deallocation ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    delete FFDBox[iZone];
  }
  delete [] FFDBox;
  cout <<"FFD container deallocated." << endl;
  
  delete [] surface_movement;
  cout <<"Surface movement container deallocated." << endl;
  
  /*--- Grid movement class deallocation ---*/
  delete [] grid_movement;
  cout <<"Grid movement container deallocated." << endl;

    /*Deallocate config container*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config_container[iZone]!=NULL){
      //delete config_container[iZone];
    }
  }
  if (config_container!=NULL) delete [] config_container;
  cout <<"Config container deallocated." << endl;
  /*--- Deallocate output container ---*/
  delete output;
  cout <<"Output container deallocated." << endl;

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/
  
#ifndef HAVE_MPI
  StopTime = double(clock())/double(CLOCKS_PER_SEC);
#else
  MPI_Barrier(MPI_COMM_WORLD);
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
  
#ifdef HAVE_MPI
  /*--- Finalize MPI parallelization ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  
  return EXIT_SUCCESS;
}
