/*!
 * \file SU2_EDU.cpp
 * \brief Main file of Computational Fluid Dynamics Code (SU2_EDU).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.9
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

#include "../include/SU2_EDU.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  bool StopCalc = false;
  unsigned long StartTime, StopTime, TimeUsed = 0, ExtIter = 0;
  unsigned short iMesh, iSol, nZone, nDim;
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
   the SU2_EDU code. In general, the pointers are instantiated down a
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

  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/
  
  char config_file_name[200], InputKey[80];
  do {
    cout << "Set inviscid (Inv) or viscous (Visc) simulation: ";
    fgets (InputKey,80,stdin);
    if (strcmp (InputKey,"Inv\n") == 0) { strcpy(config_file_name, "Inviscid_ConfigFile.cfg"); }
    if (strcmp (InputKey,"Visc\n") == 0) { strcpy(config_file_name, "Viscous_ConfigFile.cfg"); }
  } while ((strcmp (InputKey,"Inv\n") != 0) && (strcmp (InputKey,"Visc\n") != 0));
  
  /*--- Definition and of the containers for a zone. ---*/
  nZone = 1; nDim  = 2;

  solver_container      = new CSolver***[nZone];
  integration_container = new CIntegration**[nZone];
  numerics_container    = new CNumerics****[nZone];
  config_container      = new CConfig*[nZone];
  geometry_container    = new CGeometry **[nZone];
  surface_movement      = new CSurfaceMovement *[nZone];
  grid_movement         = new CVolumetricMovement *[nZone];
  
  solver_container[ZONE_0]       = NULL;
  integration_container[ZONE_0]  = NULL;
  numerics_container[ZONE_0]     = NULL;
  config_container[ZONE_0]       = NULL;
  geometry_container[ZONE_0]     = NULL;
  surface_movement[ZONE_0]       = NULL;
  grid_movement[ZONE_0]          = NULL;
  
  /*--- Definition of the configuration option class for all zones. In this
   constructor, the input configuration file is parsed and all options are
   read and stored. ---*/
  
  config_container[ZONE_0] = new CConfig(config_file_name, SU2_EDU, ZONE_0, nZone, VERB_HIGH);
  
#ifndef NO_MPI
  /*--- Change the name of the input-output files for a parallel computation ---*/
  config_container[ZONE_0]->SetFileNameDomain(rank+1);
#endif
  
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  config_container[ZONE_0]->SetNondimensionalization(nDim, ZONE_0);
  
  /*--- Definition of the geometry class. Within this constructor, the
   mesh file is read and the primal grid is stored (node coords, connectivity,
   & boundary markers. MESH_0 is the index of the finest mesh. ---*/
  
  geometry_container[ZONE_0] = new CGeometry *[config_container[ZONE_0]->GetMGLevels()+1];
  geometry_container[ZONE_0][MESH_0] = new CPhysicalGeometry(config_container[ZONE_0],ZONE_1, nZone);
  
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
  
  /*--- Definition of the solver class: solver_container[#ZONES][#MG_GRIDS][#EQ_SYSTEMS].
   The solver classes are specific to a particular set of governing equations,
   and they contain the subroutines with instructions for computing each spatial
   term of the PDE, i.e. loops over the edges to compute convective and viscous
   fluxes, loops over the nodes to compute source terms, and routines for
   imposing various boundary condition type for the PDE. ---*/
  
  solver_container[ZONE_0] = new CSolver** [config_container[ZONE_0]->GetMGLevels()+1];
  for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++)
    solver_container[ZONE_0][iMesh] = NULL;
  
  for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++) {
    solver_container[ZONE_0][iMesh] = new CSolver* [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      solver_container[ZONE_0][iMesh][iSol] = NULL;
  }
  
  Solver_Preprocessing(solver_container[ZONE_0], geometry_container[ZONE_0],
                       config_container[ZONE_0], ZONE_0);
  
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
  
  integration_container[ZONE_0] = new CIntegration*[MAX_SOLS];
  Integration_Preprocessing(integration_container[ZONE_0], geometry_container[ZONE_0],
                            config_container[ZONE_0], ZONE_0);
  
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
  
  numerics_container[ZONE_0] = new CNumerics***[config_container[ZONE_0]->GetMGLevels()+1];
  Numerics_Preprocessing(numerics_container[ZONE_0], solver_container[ZONE_0],
                         geometry_container[ZONE_0], config_container[ZONE_0], ZONE_0);
  
#ifndef NO_MPI
  /*--- Synchronization point after the solver definition subroutine ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
  /*--- Computation of wall distances for turbulence modeling ---*/
  
  if ((config_container[ZONE_0]->GetKind_Solver() == RANS))
    geometry_container[ZONE_0][MESH_0]->ComputeWall_Distance(config_container[ZONE_0]);
  
  /*--- Computation of positive surface area in the z-plane which is used for
   the calculation of force coefficient (non-dimensionalization). ---*/
  
  geometry_container[ZONE_0][MESH_0]->SetPositive_ZArea(config_container[ZONE_0]);
  
  /*--- Surface grid deformation using design variables ---*/
  if (rank == MASTER_NODE) cout << endl << "------------------------- Surface grid deformation ----------------------" << endl;
  
  /*--- Definition and initialization of the surface deformation class ---*/
  surface_movement[ZONE_0]->CopyBoundary(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
  
  /*--- Surface grid deformation ---*/
  if (rank == MASTER_NODE) cout << "Performing the deformation of the surface grid." << endl;
  surface_movement[ZONE_0]->SetAirfoil(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
  
#ifndef NO_MPI
  /*--- MPI syncronization point ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
  /*--- Volumetric grid deformation ---*/
  if (rank == MASTER_NODE) cout << endl << "----------------------- Volumetric grid deformation ---------------------" << endl;
  
  /*--- Definition of the Class for grid movement ---*/
  
  if (rank == MASTER_NODE) cout << "Performing the deformation of the volumetric grid." << endl;
  grid_movement[ZONE_0] = new CVolumetricMovement(geometry_container[ZONE_0][MESH_0]);
  grid_movement[ZONE_0]->SetVolume_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], true);
  grid_movement[ZONE_0]->UpdateMultiGrid(geometry_container[ZONE_0], config_container[ZONE_0]);
  
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
    config_container[ZONE_0]->SetExtIter(ExtIter);
    config_container[ZONE_0]->UpdateCFL(ExtIter);
    
    /*--- Perform a single iteration of the chosen PDE solver. ---*/
    MeanFlowIteration(output, integration_container, geometry_container,
                      solver_container, numerics_container, config_container,
                      surface_movement, grid_movement);
    
    /*--- Synchronization point after a single solver iteration. Compute the
     wall clock time required. ---*/
    
#ifndef NO_MPI
    MPI::COMM_WORLD.Barrier();
#endif
    StopTime = clock(); TimeUsed += (StopTime - StartTime);
    
    /*--- Update the convergence history file (serial and parallel computations). ---*/
    
    output->SetConvergence_History(&ConvHist_file, geometry_container, solver_container,
                                   config_container, integration_container, false, TimeUsed, ZONE_0);
    
    /*--- Check whether the current simulation has reached the specified
     convergence criteria, and set StopCalc to true, if so. ---*/
    
    StopCalc = integration_container[ZONE_0][FLOW_SOL]->GetConvergence();

    /*--- Solution output. Determine whether a solution needs to be written
     after the current iteration, and if so, execute the output file writing
     routines. ---*/
    
    if ((ExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
        ((ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (ExtIter != 0)) ||
        (StopCalc)) {
          
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
  //    for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++) {
  //      for (iSol = 0; iSol < MAX_SOLS; iSol++) {
  //        if (solver_container[ZONE_0][iMesh][iSol] != NULL) {
  //          delete solver_container[ZONE_0][iMesh][iSol];
  //        }
  //      }
  //      delete solver_container[ZONE_0][iMesh];
  //    }
  //    delete solver_container[ZONE_0];
  //  delete [] solver_container;
  //  if (rank == MASTER_NODE) cout <<"Solution container, deallocated." << endl;
  
  /*--- Geometry class deallocation ---*/
  //    for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++) {
  //      delete geometry_container[ZONE_0][iMesh];
  //    }
  //    delete geometry_container[ZONE_0];
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
    cout << endl <<"------------------------- Exit Success (SU2_EDU) ------------------------" << endl << endl;
  
  return EXIT_SUCCESS;
}
