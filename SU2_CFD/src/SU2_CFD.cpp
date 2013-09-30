/*!
 * \file SU2_CFD.cpp
 * \brief Main file of Computational Fluid Dynamics Code (SU2_CFD).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.7
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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
 * SU2_CFD.cpp is the main control file for the SU2 tool.
 * Pseudocode:
 *  BEGIN
 *  Initialize data structures
 *  Load in the geometry
 *  Parse the input file and load settings
 *  "Initialize solution and geometry processing" whatever that means
 *  for ExtIter = 0; ExtIter < MaxIter; ExtIter++{
 *      Make an iteration (problem dependent)
 *      Check Convergence of the flow solution (problem dependent)
 *      ??? Not sure ???
 *      if converged{
 *          break
 *      }
 *  }
 *  Write solution files
 *  END
 *
 *  Loading the geometry happens in ______.cpp
 *  Parsing the input file happens in ______.cpp
 *  Iterating the problem is problem dependent, and is located in solution_MATH_PHYSICAL.cpp
 *      MATH is MATH_PROBLEM from the .cfg file
 *      PHYSICAL is PHYSICAL_PROBLEM from the .cfg file
 */


using namespace std;

int main(int argc, char *argv[]) {
	bool StopCalc = false;
	unsigned long StartTime, StopTime, TimeUsed = 0, ExtIter = 0;
    // Add description of iMesh
	unsigned short iMesh, iZone, iSol, nZone, nDim;
	ofstream ConvHist_file;
	int rank = MASTER_NODE;
	
#ifndef NO_MPI
	/*--- MPI initialization, and buffer setting ---*/
	void *buffer, *old_buffer;
	int size, bufsize;
	bufsize = MAX_MPI_BUFFER;
	buffer = new char[bufsize];
	MPI::Init(argc, argv);
	MPI::Attach_buffer(buffer, bufsize);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
#ifdef TIME
    /*--- Set up a timer for parallel performance benchmarking ---*/
    double start, finish, time;
    MPI::COMM_WORLD.Barrier();
    start = MPI::Wtime();
#endif
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
	
	/*--- Definition and of the containers for all possible zones. ---*/
    
    // MAX_ZONES is a hard-coded constant in option_structure.hpp
	solver_container      = new CSolver***[MAX_ZONES];
	integration_container = new CIntegration**[MAX_ZONES];
	numerics_container    = new CNumerics****[MAX_ZONES];
	config_container      = new CConfig*[MAX_ZONES];
	geometry_container    = new CGeometry **[MAX_ZONES];
	surface_movement      = new CSurfaceMovement *[MAX_ZONES];
	grid_movement         = new CVolumetricMovement *[MAX_ZONES];
	FFDBox                = new CFreeFormDefBox**[MAX_ZONES];
    for (iZone = 0; iZone < MAX_ZONES; iZone++) {
        solver_container[iZone]       = NULL;
        integration_container[iZone]  = NULL;
        numerics_container[iZone]     = NULL;
        config_container[iZone]       = NULL;
        geometry_container[iZone]     = NULL;
        surface_movement[iZone]       = NULL;
        grid_movement[iZone]          = NULL;
        FFDBox[iZone]                 = NULL;
    }
	
	/*--- Load in the number of zones and spatial dimensions in the mesh file ---*/
    // If no config file is specified, default.cfg is used
    char config_file_name[200];
    if (argc == 2){
        strcpy(config_file_name,argv[1]);
    }else{
        strcpy(config_file_name, "default.cfg");
    }
	CConfig *config = NULL;
    config = new CConfig(config_file_name);
	nZone = GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
	nDim = GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());
	
    /*--- Loop over all zones to initialize the various classes. In most
     cases, nZone is equal to one. This represents the solution of a partial
     differential equation on a single block, unstructured mesh. ---*/
    
	for (iZone = 0; iZone < nZone; iZone++) {
		
		/*--- Definition of the configuration option class for all zones. In this
         constructor, the input configuration file is parsed and all options are
         read and stored. ---*/
        
		config_container[iZone] = new CConfig(config_file_name, SU2_CFD, iZone, nZone, VERB_HIGH);
		
#ifndef NO_MPI
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
		geometry_container[iZone][MESH_0] = new CPhysicalGeometry(config_container[iZone],
                                                                  iZone+1, nZone);
		
    }
    
    if (rank == MASTER_NODE)
        cout << endl <<"------------------------- Geometry Preprocessing ------------------------" << endl;
    
    /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
     based data structure is constructed, i.e. node and cell neighbors are
     identified and linked, face areas and volumes of the dual mesh cells are
     computed, and the multigrid levels are created using an agglomeration procedure. ---*/
    
    Geometrical_Preprocessing(geometry_container, config_container, nZone);
    
#ifndef NO_MPI
    /*--- Synchronization point after the geometrical definition subroutine ---*/
    MPI::COMM_WORLD.Barrier();
#endif
    
    if (rank == MASTER_NODE)
        cout << endl <<"------------------------- Solver Preprocessing --------------------------" << endl;
    
    for (iZone = 0; iZone < nZone; iZone++) {
        
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
		
#ifndef NO_MPI
		/*--- Synchronization point after the solution preprocessing subroutine ---*/
		MPI::COMM_WORLD.Barrier();
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
        
#ifndef NO_MPI
		/*--- Synchronization point after the integration definition subroutine ---*/
		MPI::COMM_WORLD.Barrier();
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
        
#ifndef NO_MPI
		/*--- Synchronization point after the solver definition subroutine ---*/
		MPI::COMM_WORLD.Barrier();
#endif
        
		/*--- Computation of wall distances for turbulence modeling ---*/
        
		if ( (config_container[iZone]->GetKind_Solver() == RANS)     ||
         (config_container[iZone]->GetKind_Solver() == ADJ_RANS)    )
			geometry_container[iZone][MESH_0]->ComputeWall_Distance(config_container[iZone]);
        
		/*--- Computation of positive surface area in the z-plane which is used for
         the calculation of force coefficient (non-dimensionalization). ---*/
        
		geometry_container[iZone][MESH_0]->SetPositive_ZArea(config_container[iZone]);
        
		/*--- Set the near-field and interface boundary conditions, if necessary. ---*/
        
		for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
			geometry_container[iZone][iMesh]->MatchNearField(config_container[iZone]);
			geometry_container[iZone][iMesh]->MatchInterface(config_container[iZone]);
		}
        
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
                                 FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, 0);
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
	
	while (ExtIter < config_container[ZONE_0]->GetnExtIter()) {
        
		/*--- Set a timer for each iteration. Store the current iteration and
         update  the value of the CFL number (if there is CFL ramping specified)
         in the config class. ---*/
        
		StartTime = clock();
		for (iZone = 0; iZone < nZone; iZone++) {
			config_container[iZone]->SetExtIter(ExtIter);
			config_container[iZone]->UpdateCFL(ExtIter);
		}
        
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
				
			case PLASMA_EULER: case PLASMA_NAVIER_STOKES:
				PlasmaIteration(output, integration_container, geometry_container,
                                solver_container, numerics_container, config_container,
                                surface_movement, grid_movement, FFDBox);
				break;
				
			case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
				FluidStructureIteration(output, integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox);
				break;
				
			case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES:
				AeroacousticIteration(output, integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox);
				break;
				
			case WAVE_EQUATION:
				WaveIteration(output, integration_container, geometry_container,
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
				
        
			case ADJ_PLASMA_EULER: case ADJ_PLASMA_NAVIER_STOKES:
				AdjPlasmaIteration(output, integration_container, geometry_container,
                                   solver_container, numerics_container, config_container,
                                   surface_movement, grid_movement, FFDBox);
				break;
        
        
			case ADJ_AEROACOUSTIC_EULER:
				AdjAeroacousticIteration(output, integration_container, geometry_container,
                                         solver_container, numerics_container, config_container,
                                         surface_movement, grid_movement, FFDBox);
				break;
		}
		
        
        /*--- Synchronization point after a single solver iteration. Compute the
         wall clock time required. ---*/
        
#ifndef NO_MPI
		MPI::COMM_WORLD.Barrier();
#endif
		StopTime = clock(); TimeUsed += (StopTime - StartTime);
        
        /*--- For specific applications, evaluate and plot the equivalent area or flow rate. ---*/
        
        if ((config_container[ZONE_0]->GetKind_Solver() == EULER) &&
            (config_container[ZONE_0]->GetEquivArea() == YES)) {
            output->SetEquivalentArea(solver_container[ZONE_0][MESH_0][FLOW_SOL],
                                      geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
        }
        
		/*--- Update the convergence history file (serial and parallel computations). ---*/
        
        output->SetConvergence_History(&ConvHist_file, geometry_container, solver_container,
                                       config_container, integration_container, false, TimeUsed, ZONE_0);
        
		/*--- Check whether the current simulation has reached the specified
         convergence criteria, and set StopCalc to true, if so. ---*/
        
		switch (config_container[ZONE_0]->GetKind_Solver()) {
			case EULER: case NAVIER_STOKES: case RANS:
				StopCalc = integration_container[ZONE_0][FLOW_SOL]->GetConvergence(); break;
			case PLASMA_EULER: case PLASMA_NAVIER_STOKES:
				StopCalc = integration_container[ZONE_0][PLASMA_SOL]->GetConvergence(); break;
			case WAVE_EQUATION:
				StopCalc = integration_container[ZONE_0][WAVE_SOL]->GetConvergence(); break;
			case LINEAR_ELASTICITY:
				StopCalc = integration_container[ZONE_0][FEA_SOL]->GetConvergence(); break;
			case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
				StopCalc = integration_container[ZONE_0][ADJFLOW_SOL]->GetConvergence(); break;
			case ADJ_PLASMA_EULER: case ADJ_PLASMA_NAVIER_STOKES:
				StopCalc = integration_container[ZONE_0][ADJPLASMA_SOL]->GetConvergence(); break;
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
                    solver_container[ZONE_0][MESH_0][config_container[ZONE_0]->GetContainerPosition(RUNTIME_FLOW_SYS)]->Preprocessing(geometry_container[ZONE_0][MESH_0], solver_container[ZONE_0][MESH_0], config_container[ZONE_0], MESH_0, 0, RUNTIME_FLOW_SYS);
                }
                
                /*--- Execute the routine for writing restart, volume solution,
                 surface solution, and surface comma-separated value files. ---*/
                
                output->SetResult_Files(solver_container, geometry_container, config_container, ExtIter, nZone);
                
            }
        
		/*--- If the convergence criteria has been met, terminate the simulation. ---*/
        
		if (StopCalc) break;
		ExtIter++;
        
	}
    
    /*--- Close the convergence history file. ---*/
    
    if (rank == MASTER_NODE) {
        ConvHist_file.close();
        cout << endl <<"History file, closed." << endl;
    }
    
    /*--- Solver class deallocation ---*/
    //  for (iZone = 0; iZone < nZone; iZone++) {
    //    for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
    //      for (iSol = 0; iSol < MAX_SOLS; iSol++) {
    //        if (solver_container[iZone][iMesh][iSol] != NULL) {
    //          delete solver_container[iZone][iMesh][iSol];
    //        }
    //      }
    //      delete solver_container[iZone][iMesh];
    //    }
    //    delete solver_container[iZone];
    //  }
    //  delete [] solver_container;
    //  if (rank == MASTER_NODE) cout <<"Solution container, deallocated." << endl;
    
    /*--- Geometry class deallocation ---*/
    //  for (iZone = 0; iZone < nZone; iZone++) {
    //    for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
    //      delete geometry_container[iZone][iMesh];
    //    }
    //    delete geometry_container[iZone];
    //  }
    //  delete [] geometry_container;
    //  cout <<"Geometry container, deallocated." << endl;
    
    /*--- Integration class deallocation ---*/
    //  cout <<"Integration container, deallocated." << endl;
    
#ifndef NO_MPI
    /*--- Compute/print the total time for parallel performance benchmarking. ---*/
#ifdef TIME
    MPI::COMM_WORLD.Barrier();
    finish = MPI::Wtime();
    time = finish-start;
    if (rank == MASTER_NODE) {
        cout << "\nCompleted in " << fixed << time << " seconds on "<< size;
        if (size == 1) cout << " core.\n" << endl;
        else cout << " cores.\n" << endl;
    }
#endif
    /*--- Finalize MPI parallelization ---*/
    old_buffer = buffer;
    MPI::Detach_buffer(old_buffer);
    //	delete [] buffer;
    MPI::Finalize();
#endif
    
    /*--- Exit the solver cleanly ---*/
    
    if (rank == MASTER_NODE)
        cout << endl <<"------------------------- Exit Success (SU2_CFD) ------------------------" << endl << endl;
    
    return EXIT_SUCCESS;
}
