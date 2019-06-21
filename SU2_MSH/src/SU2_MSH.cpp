/*!
 * \file SU2_MSH.cpp
 * \brief Main file of Mesh Adaptation Code (SU2_MSH).
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/SU2_MSH.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	
	/*--- Variable definitions ---*/
  
  unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  char config_file_name[MAX_STRING_SIZE];
  char file_name[MAX_STRING_SIZE];
  int rank, size;
  string str;
  
  /*--- MPI initialization ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Init(&argc,&argv);
  SU2_MPI::Comm MPICommunicator(MPI_COMM_WORLD);
#else
  SU2_Comm MPICommunicator(0);
#endif

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
	
  /*--- Pointer to different structures that will be used throughout the entire code ---*/
  
  CConfig **config_container         = NULL;
  CGeometry **geometry_container     = NULL;
  
  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/
  
  if (argc == 2) { strcpy(config_file_name,argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }
  
	  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation)  ---*/

  CConfig *config = NULL;
  config = new CConfig(config_file_name, SU2_MSH);

  nZone    = CConfig::GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);

  /*--- Definition of the containers per zones ---*/
  
  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry*[nZone];
  
  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
  }
  
  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/
    
    config_container[iZone] = new CConfig(config_file_name, SU2_MSH, iZone, nZone, 0, true);
    config_container[iZone]->SetMPICommunicator(MPICommunicator);
    
    /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/
    
    CGeometry *geometry_aux = NULL;
    
    /*--- All ranks process the grid and call ParMETIS for partitioning ---*/
    
    geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);
    
    /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/
    
    geometry_aux->SetColorGrid_Parallel(config_container[iZone]);
    
    /*--- Build the grid data structures using the ParMETIS coloring. ---*/

    geometry_container[iZone] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);
    
    /*--- Deallocate the memory of geometry_aux ---*/
    
    delete geometry_aux;
    
    /*--- Add the Send/Receive boundaries ---*/
    
    geometry_container[iZone]->SetSendReceive(config_container[iZone]);
    
    /*--- Add the Send/Receive boundaries ---*/
    
    geometry_container[iZone]->SetBoundaries(config_container[iZone]);
    
  }
  
  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/
  
#ifdef HAVE_MPI
  StartTime = MPI_Wtime();
#else
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif
  
	cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
	/*--- Compute elements surrounding points, points surrounding points, and elements surronding elements ---*/
  
	cout << "Setting local point and element connectivity." <<endl;
	geometry_container[ZONE_0]->SetPoint_Connectivity(); geometry_container[ZONE_0]->SetElement_Connectivity();
	
	/*--- Check the orientation before computing geometrical quantities ---*/
  geometry_container[ZONE_0]->SetBoundVolume();
  if (config_container[ZONE_0]->GetReorientElements()) {
		cout << "Check numerical grid orientation." <<endl;
		geometry_container[ZONE_0]->Check_IntElem_Orientation(config_container[ZONE_0]); geometry_container[ZONE_0]->Check_BoundElem_Orientation(config_container[ZONE_0]);
  }
	
	/*--- Create the edge structure ---*/
  
	cout << "Identify faces, edges and vertices." <<endl;
	geometry_container[ZONE_0]->SetFaces(); geometry_container[ZONE_0]->SetEdges(); geometry_container[ZONE_0]->SetVertex(config_container[ZONE_0]); geometry_container[ZONE_0]->SetCoord_CG();
	
	/*--- Create the control volume structures ---*/
  
	cout << "Set control volume structure." << endl;
	geometry_container[ZONE_0]->SetControlVolume(config_container[ZONE_0], ALLOCATE); geometry_container[ZONE_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);

  /*--- Create the point-to-point MPI communication structures. ---*/
  
  geometry_container[ZONE_0]->PreprocessP2PComms(geometry_container[ZONE_0], config_container[ZONE_0]);
	
	if ((config_container[ZONE_0]->GetKind_Adaptation() != NONE) && (config_container[ZONE_0]->GetKind_Adaptation() != PERIODIC)) {
		
		cout << endl <<"--------------------- Start numerical grid adaptation -------------------" << endl;
		
		/*-- Definition of the Class for grid adaptation ---*/
    
		CGridAdaptation *grid_adaptation;
		grid_adaptation = new CGridAdaptation(geometry_container[ZONE_0], config_container[ZONE_0]);
		
		/*--- Read the flow solution and/or the adjoint solution
		 and choose the elements to adapt ---*/
    
		if ((config_container[ZONE_0]->GetKind_Adaptation() != FULL)
				&& (config_container[ZONE_0]->GetKind_Adaptation() != WAKE) && (config_container[ZONE_0]->GetKind_Adaptation() != SMOOTHING) && (config_container[ZONE_0]->GetKind_Adaptation() != SUPERSONIC_SHOCK))
			grid_adaptation->GetFlowSolution(geometry_container[ZONE_0], config_container[ZONE_0]);
		
		switch (config_container[ZONE_0]->GetKind_Adaptation()) {
			case NONE:
				break;
			case SMOOTHING:
				config_container[ZONE_0]->SetSmoothNumGrid(true);
				grid_adaptation->SetNo_Refinement(geometry_container[ZONE_0], 1);
				break;
			case FULL:
				grid_adaptation->SetComplete_Refinement(geometry_container[ZONE_0], 1);
				break;
			case WAKE:
				grid_adaptation->SetWake_Refinement(geometry_container[ZONE_0], 1);
				break;
			case SUPERSONIC_SHOCK:
				grid_adaptation->SetSupShock_Refinement(geometry_container[ZONE_0], config_container[ZONE_0]);
				break;
			case FULL_FLOW:
				grid_adaptation->SetComplete_Refinement(geometry_container[ZONE_0], 1);
				break;
			case FULL_ADJOINT:
				grid_adaptation->GetAdjSolution(geometry_container[ZONE_0], config_container[ZONE_0]);
				grid_adaptation->SetComplete_Refinement(geometry_container[ZONE_0], 1);
				break;
			case GRAD_FLOW:
				grid_adaptation->SetIndicator_Flow(geometry_container[ZONE_0], config_container[ZONE_0], 1);
				break;
			case GRAD_ADJOINT:
				grid_adaptation->GetAdjSolution(geometry_container[ZONE_0], config_container[ZONE_0]);
				grid_adaptation->SetIndicator_Adj(geometry_container[ZONE_0], config_container[ZONE_0], 1);
				break;
			case GRAD_FLOW_ADJ:
				grid_adaptation->GetAdjSolution(geometry_container[ZONE_0], config_container[ZONE_0]);
				grid_adaptation->SetIndicator_FlowAdj(geometry_container[ZONE_0], config_container[ZONE_0]);
				break;
			case COMPUTABLE:
				grid_adaptation->GetAdjSolution(geometry_container[ZONE_0], config_container[ZONE_0]);
				grid_adaptation->GetFlowResidual(geometry_container[ZONE_0], config_container[ZONE_0]);
				grid_adaptation->SetIndicator_Computable(geometry_container[ZONE_0], config_container[ZONE_0]);
				break;
			case REMAINING:
        SU2_MPI::Error("Adaptation method not implemented.", CURRENT_FUNCTION);
				break;
			default :
				cout << "The adaptation is not defined" << endl;
		}
		
		/*--- Perform an homothetic adaptation of the grid ---*/
    
		CPhysicalGeometry *geo_adapt; geo_adapt = new CPhysicalGeometry;
		
		cout << "Homothetic grid adaptation" << endl;
		if (geometry_container[ZONE_0]->GetnDim() == 2) grid_adaptation->SetHomothetic_Adaptation2D(geometry_container[ZONE_0], geo_adapt, config_container[ZONE_0]);
		if (geometry_container[ZONE_0]->GetnDim() == 3) grid_adaptation->SetHomothetic_Adaptation3D(geometry_container[ZONE_0], geo_adapt, config_container[ZONE_0]);
    
		/*--- Smooth the numerical grid coordinates ---*/
    
		if (config_container[ZONE_0]->GetSmoothNumGrid()) {
			cout << "Preprocessing for doing the implicit smoothing." << endl;
			geo_adapt->SetPoint_Connectivity(); geo_adapt->SetElement_Connectivity();
			geo_adapt->SetBoundVolume();
			if (config_container[ZONE_0]->GetReorientElements()) {
				geo_adapt->Check_IntElem_Orientation(config_container[ZONE_0]); geo_adapt->Check_BoundElem_Orientation(config_container[ZONE_0]);
			}
			geo_adapt->SetEdges(); geo_adapt->SetVertex(config_container[ZONE_0]);
			cout << "Implicit smoothing of the numerical grid coordinates." << endl;
			geo_adapt->SetCoord_Smoothing(5, 1.5, config_container[ZONE_0]);
		}
		
		/*--- Original and adapted grid ---*/
    strcpy (file_name, "original_grid.dat");
    geometry_container[ZONE_0]->SetTecPlot(file_name, true);
    strcpy (file_name, "original_surface.dat");
    geometry_container[ZONE_0]->SetBoundTecPlot(file_name, true, config_container[ZONE_0]);
    
		/*--- Write the adapted grid sensor ---*/
    
    strcpy (file_name, "adapted_grid.dat");
    geo_adapt->SetTecPlot(file_name, true);
    strcpy (file_name, "adapted_surface.dat");
    geo_adapt->SetBoundTecPlot(file_name, true, config_container[ZONE_0]);
		
		/*--- Write the new adapted grid, including the modified boundaries surfaces ---*/
    
		geo_adapt->SetMeshFile(config_container[ZONE_0], config_container[ZONE_0]->GetMesh_Out_FileName());
    
    
		/*--- Write the restart file ---*/
    
		if ((config_container[ZONE_0]->GetKind_Adaptation() != SMOOTHING) && (config_container[ZONE_0]->GetKind_Adaptation() != FULL) &&
				(config_container[ZONE_0]->GetKind_Adaptation() != WAKE) &&
				(config_container[ZONE_0]->GetKind_Adaptation() != SUPERSONIC_SHOCK))
			grid_adaptation->SetRestart_FlowSolution(config_container[ZONE_0], geo_adapt, config_container[ZONE_0]->GetRestart_FlowFileName());
		
		if ((config_container[ZONE_0]->GetKind_Adaptation() == GRAD_FLOW_ADJ) || (config_container[ZONE_0]->GetKind_Adaptation() == GRAD_ADJOINT)
				|| (config_container[ZONE_0]->GetKind_Adaptation() == FULL_ADJOINT) || (config_container[ZONE_0]->GetKind_Adaptation() == COMPUTABLE) ||
				(config_container[ZONE_0]->GetKind_Adaptation() == REMAINING))
			grid_adaptation->SetRestart_AdjSolution(config_container[ZONE_0], geo_adapt, config_container[ZONE_0]->GetRestart_AdjFileName());
		
	}
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;
  
  if (geometry_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (geometry_container[iZone] != NULL) {
        delete geometry_container[iZone];
      }
    }
    delete [] geometry_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;
  
  
  if (config_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != NULL) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;
  
  delete config;
  config = NULL;

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/
  
#ifdef HAVE_MPI
  StopTime = MPI_Wtime();
#else
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif
  
  /*--- Compute/print the total time for performance benchmarking. ---*/
  
  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }
  
  /*--- Exit the solver cleanly ---*/
  
	cout << endl <<"------------------------- Exit Success (SU2_MSH) ------------------------" << endl << endl;
  
  /*--- Finalize MPI parallelization ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}

