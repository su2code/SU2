/*!
 * \file SU2_FSI.cpp
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

#include "../include/SU2_interface.hpp"



//SU2_interface::SU2_interface(char case_filename[200],MPI_Comm SU2_comm_interface){
SU2_interface::SU2_interface(char case_filename[200]){
    
    
    //MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    //MPI_Comm_dup(SU2_comm_interface,&SU2_comm_local);
    
    //MPI_Comm_dup(SU2_comm_interface,&SU2_comm);
    //SU2_comm_local = MPI_COMM_WORLD;
    
   // SU2_comm = SU2_comm_interface;
//    bool StopCalc = false;
//    double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
    ExtIter = 0;
//    unsigned short iMesh, iZone, iSol, nZone, nDim;
//    ofstream ConvHist_file;
    rank = MASTER_NODE;
    size = SINGLE_NODE;
    
    
#ifdef HAVE_MPI
//    /*--- MPI initialization, and buffer setting ---*/
    int *bptr, bl;
//    MPI_Init(&argc,&argv);
    MPI_Buffer_attach( malloc(BUFSIZE), BUFSIZE );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    cout<<"\n Rank : "<<rank<<endl;
    
    /*--- Read the name and format of the input mesh file ---*/
    config = new CConfig(case_filename,SU2_CFD);
    
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
        
        config_container[iZone] = new CConfig(case_filename, SU2_CFD, iZone, nZone, nDim, VERB_HIGH);
        
        
//        /*--- Definition of the geometry class to store the primal grid in the
//         partitioning process. ---*/
//        *geometry_aux = NULL
        
        /*--- All ranks process the grid and call ParMETIS for partitioning ---*/
        
        geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);
        
        /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/
        
        geometry_aux->SetColorGrid_Parallel(config_container[iZone]);
        
        /*--- Allocate the memory of the current domain, and divide the grid
         between the ranks. ---*/
        
        geometry_container[iZone] = new CGeometry *[config_container[iZone]->GetnMGLevels()+1];
        geometry_container[iZone][MESH_0] = new CPhysicalGeometry(geometry_aux, config_container[iZone], 1);
        
        /*--- Deallocate the memory of geometry_aux ---*/
        
        delete geometry_aux;
        
        /*--- Add the Send/Receive boundaries ---*/
        
        geometry_container[iZone][MESH_0]->SetSendReceive(config_container[iZone]);
        
        /*--- Add the Send/Receive boundaries ---*/
        
        geometry_container[iZone][MESH_0]->SetBoundaries(config_container[iZone]);
        
    }
    
    
    if (rank == MASTER_NODE)
        cout << endl <<"------------------------- Geometry Preprocessing ------------------------" << endl;
    
    /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
     based data structure is constructed, i.e. node and cell neighbors are
     identified and linked, face areas and volumes of the dual mesh cells are
     computed, and the multigrid levels are created using an agglomeration procedure. ---*/
    
    Geometrical_Preprocessing(geometry_container, config_container, nZone);
    

    
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
        
        for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
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
        
        solver_container[iZone] = new CSolver** [config_container[iZone]->GetnMGLevels()+1];
        for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++)
            solver_container[iZone][iMesh] = NULL;
        
        for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
            solver_container[iZone][iMesh] = new CSolver* [MAX_SOLS];
            for (iSol = 0; iSol < MAX_SOLS; iSol++)
                solver_container[iZone][iMesh][iSol] = NULL;
        }
        Solver_Preprocessing(solver_container[iZone], geometry_container[iZone],
                             config_container[iZone], iZone);
        

        
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
        

        
        /*--- Definition of the numerical method class:
         numerics_container[#ZONES][#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
         The numerics class contains the implementation of the numerical methods for
         evaluating convective or viscous fluxes between any two nodes in the edge-based
         data structure (centered, upwind, galerkin), as well as any source terms
         (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/
        
        numerics_container[iZone] = new CNumerics***[config_container[iZone]->GetnMGLevels()+1];
        Numerics_Preprocessing(numerics_container[iZone], solver_container[iZone],
                               geometry_container[iZone], config_container[iZone], iZone);
        
        if (rank == MASTER_NODE) cout << "Numerics Preprocessing." << endl;
        

        
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
    output->SetConvHistory_Header(&ConvHist_file, config_container[ZONE_0]);
    
    /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
    if (config_container[ZONE_0]->GetWrt_Unsteady() && config_container[ZONE_0]->GetRestart())
    ExtIter = config_container[ZONE_0]->GetUnst_RestartIter();


    #ifndef HAVE_MPI
        StartTime = double(clock())/double(CLOCKS_PER_SEC);
    #else
        StartTime = MPI_Wtime();
    #endif
}




void SU2_interface::Run_Steady_Cfd(void){
    
    #ifdef HAVE_MPI
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif
    
    
    StopCalc =  false;
    nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
    while (ExtIter < config_container[ZONE_0]->GetnExtIter()) {
        
        /*--- Set a timer for each iteration. Store the current iteration and
         update  the value of the CFL number (if there is CFL ramping specified)
         in the config class. ---*/
        
        config_container[ZONE_0]->SetExtIter(ExtIter);
        
//        IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
//        ExtIter = config_container[ZONE_0]->GetExtIter();
        
//        /*--- Read the target pressure ---*/
//        
//        if (config_container[ZONE_0]->GetInvDesign_Cp() == YES)
//            output->SetCp_InverseDesign(solver_container[ZONE_0][MESH_0][FLOW_SOL],
//                                        geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
//        
//        /*--- Read the target heat flux ---*/
//        
//        if (config_container[ZONE_0]->GetInvDesign_HeatFlux() == YES)
//            output->SetHeat_InverseDesign(solver_container[ZONE_0][MESH_0][FLOW_SOL],
//                                          geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
        
        MeanFlowIteration(output, integration_container, geometry_container,
                                  solver_container, numerics_container, config_container,
                                  surface_movement, grid_movement, FFDBox);
        
        
        /*--- Synchronization point after a single solver iteration. Compute the
         wall clock time required. ---*/
        
#ifndef HAVE_MPI
        StopTime = double(clock())/double(CLOCKS_PER_SEC);
#else
        StopTime = MPI_Wtime();
#endif
        
        UsedTime = (StopTime - StartTime);
        
        /*--- For specific applications, evaluate and plot the equivalent area. ---*/
        
        if (config_container[ZONE_0]->GetEquivArea() == YES) {
            output->SetEquivalentArea(solver_container[ZONE_0][MESH_0][FLOW_SOL],
                                      geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
        }
        
        /*--- Check if there is any change in the runtime parameters ---*/
        
        CConfig *runtime = NULL;
        strcpy(runtime_file_name, "runtime.dat");
        runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
        
        
        /*--- Update the convergence history file (serial and parallel computations). ---*/
        
        output->SetConvHistory_Body(&ConvHist_file, geometry_container, solver_container,
                                    config_container, integration_container, false, UsedTime, ZONE_0);
        
        
        /*--- Evaluate the new CFL number (adaptive). ---*/
        
        if (config_container[ZONE_0]->GetCFL_Adapt() == YES) {
            output->SetCFL_Number(solver_container, config_container, ZONE_0);
        }
        
        
        /*--- Check whether the current simulation has reached the specified
         convergence criteria, and set StopCalc to true, if so. ---*/
        

        StopCalc = integration_container[ZONE_0][FLOW_SOL]->GetConvergence();
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
                    integration_container[ZONE_0][FLOW_SOL]->SetProlongated_Solution(RUNTIME_FLOW_SYS, solver_container[ZONE_0][MESH_0][FLOW_SOL], solver_container[ZONE_0][MESH_1][FLOW_SOL], geometry_container[ZONE_0][MESH_0], geometry_container[ZONE_0][MESH_1], config_container[ZONE_0]);
                    integration_container[ZONE_0][FLOW_SOL]->Smooth_Solution(RUNTIME_FLOW_SYS, solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], 3, 1.25, config_container[ZONE_0]);
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
                
                //          /*--- Compute 1D output. ---*/
                //          if (config->GetWrt_1D_Output())
                //            output->OneDimensionalOutput(solver_container[ZONE_0][MESH_0][FLOW_SOL],
                //                                         geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
                
            }
        
        /*--- If the convergence criteria has been met, terminate the simulation. ---*/
        
        if (StopCalc) break;
        
        ExtIter++;
        
    }
    
    ExtIter_f = ExtIter;
    
    ExtIter = 0;
    
}




void SU2_interface::Write_Output(void){

    /*--- Output some information to the console. ---*/
    if (rank == MASTER_NODE) {
        cout << endl;
        
        /*--- Print out the number of non-physical points and reconstructions ---*/
        if (config_container[ZONE_0]->GetNonphysical_Points() > 0)
            cout << "Warning: there are " << config_container[ZONE_0]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
        if (config_container[ZONE_0]->GetNonphysical_Reconstr() > 0)
            cout << "Warning: " << config_container[ZONE_0]->GetNonphysical_Reconstr() << " reconstructed states for upwinding are non-physical." << endl;
        
        /*--- Close the convergence history file. ---*/
        ConvHist_file.close();
        cout << "History file, closed." << endl;
    }
    
    //  /*--- Deallocate config container ---*/
    //
    //  for (iZone = 0; iZone < nZone; iZone++) {
    //    if (config_container[iZone] != NULL) {
    //      delete config_container[iZone];
    //    }
    //  }
    //  if (config_container != NULL) delete[] config_container;
    
    
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

    
    
    
}


void SU2_interface::Run_Steady_Iteration(unsigned long no_of_iterations){
    
    double Physical_dt, Physical_t;
    unsigned short iMesh, iZone;
    
    bool time_spectral = (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_SPECTRAL);
    unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
    if (time_spectral){
        nZone = config_container[ZONE_0]->GetnTimeInstances();
    }
    unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
    unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
    
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    
    
    
    
    for (unsigned long iters=0;iters<no_of_iterations;iters++){
        
    /*--- Set the initial condition ---*/
    
    for (iZone = 0; iZone < nZone; iZone++)
        solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
    
    /*--- Initial set up for unsteady problems with dynamic meshes. ---*/
    
    for (iZone = 0; iZone < nZone; iZone++) {
        
        /*--- Dynamic mesh update ---*/
        
        if ((config_container[iZone]->GetGrid_Movement()) && (!time_spectral)){
            SetGrid_Movement(geometry_container[iZone], surface_movement[iZone], grid_movement[iZone], FFDBox[iZone], solver_container[iZone],config_container[iZone], iZone, IntIter, ExtIter);
        }
        
        /*--- Apply a Wind Gust ---*/
        
        if (config_container[ZONE_0]->GetWind_Gust()){
            SetWind_GustField(config_container[iZone],geometry_container[iZone],solver_container[iZone]);
        }
    }
    
    for (iZone = 0; iZone < nZone; iZone++) {
        
        /*--- Set the value of the internal iteration ---*/
        
        IntIter = ExtIter;
        if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
            (config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
        
        /*--- Update global parameters ---*/
        
        if (config_container[iZone]->GetKind_Solver() == EULER){
            config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
        }
        if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES){
            config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
        }
        if (config_container[iZone]->GetKind_Solver() == RANS){
            config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
        }
        
        /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
        
        integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                    config_container, RUNTIME_FLOW_SYS, IntIter, iZone);
        
        if (config_container[iZone]->GetKind_Solver() == RANS) {
            
            /*--- Solve the turbulence model ---*/
            
            config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
            integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                         config_container, RUNTIME_TURB_SYS, IntIter, iZone);
            
            /*--- Solve transition model ---*/
            
            if (config_container[iZone]->GetKind_Trans_Model() == LM) {
                config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
                integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                              config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
            }
            
        }
        
        /*--- Compute & store time-spectral source terms across all zones ---*/
        
        if (time_spectral)
            SetTimeSpectral(geometry_container, solver_container, config_container, nZone, (iZone+1)%nZone);
        
    }
    
//    /*--- Update the convergence history file (serial and parallel computations). ---*/
//    
//    output->SetConvergence_History(&ConvHist_file, geometry_container, solver_container,
//                                   config_container, integration_container, false, UsedTime, ZONE_0);
        
        
        
    }
    
    
}



void SU2_interface::Deform_Mesh(void){

    
#ifdef HAVE_MPI
   // MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif
    
    /*--- Surface grid deformation using design variables ---*/
    
    if (rank == MASTER_NODE) cout << endl << "------------------------- Surface grid deformation ----------------------" << endl;
    
    /*--- Surface grid deformation ---*/
    
    if (rank == MASTER_NODE) cout << "Performing the deformation of the surface grid." << endl;
   // cout<<"Zone 0 : "<<ZONE_0<<"MESH 0 : "<<MESH_0<<endl;
    //surface_movement[ZONE_0]->CopyBoundary(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
    //surface_movement[ZONE_0]->SetSurface_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
    
    
    
    /*--- External surface file based ---*/
    
    if(config_container[ZONE_0]->GetDesign_Variable(0) == SURFACE_FILE) {
        
        /*--- Check whether a surface file exists for input ---*/
        ofstream Surface_File;
        string filename = config_container[ZONE_0]->GetMotion_FileName();
        Surface_File.open(filename.c_str(), ios::in);
        
        /*--- A surface file does not exist, so write a new one for the
         markers that are specified as part of the motion. ---*/
        if (Surface_File.fail()) {
            
            if (rank == MASTER_NODE)
            cout << "No surface file found. Writing a new file: " << filename << "." << endl;
            
            Surface_File.open(filename.c_str(), ios::out);
            Surface_File.precision(15);
            unsigned long iMarker, jPoint, GlobalIndex, iVertex; double *Coords;
            for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
                if (config_container[ZONE_0]->GetMarker_All_DV(iMarker) == YES) {
                    for(iVertex = 0; iVertex < geometry_container[ZONE_0][MESH_0]->nVertex[iMarker]; iVertex++) {
                        jPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
                        GlobalIndex = geometry_container[ZONE_0][MESH_0]->node[jPoint]->GetGlobalIndex();
                        Coords = geometry_container[ZONE_0][MESH_0]->node[jPoint]->GetCoord();
                        Surface_File << GlobalIndex << "\t" << Coords[0] << "\t" << Coords[1];
                        if (geometry_container[ZONE_0][MESH_0]->GetnDim() == 2) Surface_File << endl;
                        else Surface_File << "\t" << Coords[2] << endl;
                        
                    }
                }
            }
            Surface_File.close();
            
            cout<<"\n"<<"Got here with reading file"<<endl;
            
            /*--- A surface file exists, so read in the coordinates ---*/
            
        }
        
        else {
            Surface_File.close();
            if (rank == MASTER_NODE) cout << "Updating the surface coordinates from the input file." << endl;
            surface_movement[ZONE_0]->SetExternal_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ZONE_0, 0);
        }
        
    }

    
    
    
    
    
    if (rank == MASTER_NODE) cout << "Performing the deformation of the volumetric grid." << endl;
    
    grid_movement[ZONE_0]->SetVolume_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], false);
    
    
    if (config_container[ZONE_0]->GetVisualize_Deformation()) {
        
//        //output = new COutput();
//        int rank = MASTER_NODE;
//#ifdef HAVE_MPI
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//#endif
        
//        unsigned short iZone;
//        bool new_file = false;
//        unsigned short val_nZone = SINGLE_ZONE;
//        
//        for (iZone = 0; iZone < val_nZone; iZone++) {
//            
//            /*--- Flags identifying the types of files to be written. ---*/
//            
//            bool Wrt_Vol = config_container[iZone]->GetWrt_Vol_Sol();
//            bool Wrt_Srf = config_container[iZone]->GetWrt_Srf_Sol();
//            
//            /*--- Merge the node coordinates and connectivity if necessary. This
//             is only performed if a volume solution file is requested, and it
//             is active by default. ---*/
//            
//            if (Wrt_Vol || Wrt_Srf) {
//                if (rank == MASTER_NODE) cout <<"Merging grid connectivity." << endl;
//                output->MergeConnectivity(config_container[iZone], geometry_container[iZone][MESH_0], iZone);
//            }
//            
//            /*--- Merge coordinates of all grid nodes (excluding ghost points).
//             The grid coordinates are always merged and included first in the
//             restart files. ---*/
//            
//            output->MergeCoordinates(config_container[iZone], geometry_container[iZone][MESH_0]);
//            
//            /*--- Write restart, CGNS, Tecplot or Paraview files using the merged data.
//             This data lives only on the master, and these routines are currently
//             executed by the master proc alone (as if in serial). ---*/
//            
//            if (rank == MASTER_NODE) {
//                
//                if (Wrt_Vol) {
//                    
//                    if (rank == MASTER_NODE) cout <<"Writing volume mesh file." << endl;
//                    
//                    /*--- Write a Tecplot ASCII file ---*/
//                    
//                    output->SetTecplot_MeshASCII(config_container[iZone], geometry_container[iZone][MESH_0], false, new_file);
//                    output->DeallocateConnectivity(config_container[iZone], geometry_container[iZone][MESH_0], false);
//                    
//                }
//                
//                if (Wrt_Srf) {
//                    
//                    if (rank == MASTER_NODE) cout <<"Writing surface mesh file." << endl;
//                    
//                    /*--- Write a Tecplot ASCII file ---*/
//                    
//                    output->SetTecplot_MeshASCII(config_container[iZone], geometry_container[iZone][MESH_0], true, new_file);
//                    output->DeallocateConnectivity(config_container[iZone], geometry_container[iZone][MESH_0], true);
//                    
//                }
//                
//            }
//            
//            /*--- Final broadcast (informing other procs that the base output
//             file was written) & barrier to sync up after master node writes
//             output files. ---*/
//            
//#ifdef HAVE_MPI
//            //MPI_Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
//            MPI_Barrier(MPI_COMM_WORLD);
//#endif
//            
//        }
        

        
    }
    


    
}



void SU2_interface::Write_Surface_Mesh(unsigned long iterations){
    
    unsigned short val_nZone;
    val_nZone = nZone;
    
    for (iZone = 0; iZone < val_nZone; iZone++) {
        
        bool Wrt_Csv = config_container[iZone]->GetWrt_Csv_Sol();
        //bool Wrt_Rst = config_container[iZone]->GetWrt_Restart();
    
    if (Wrt_Csv) output->SetSurfaceCSV_Flow(config_container[iZone], geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][FLOW_SOL], ExtIter_f, iZone);
    
    }
    
    

    
    
}


void SU2_interface::Write_Final_Output(void){
    
    output->SetResult_Files(solver_container, geometry_container, config_container, ExtIter_f, nZone);
    
    
    
}

