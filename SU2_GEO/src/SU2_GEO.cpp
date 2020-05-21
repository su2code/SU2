/*!
 * \file SU2_GEO.cpp
 * \brief Main file of the Geometry Definition Code (SU2_GEO).
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../include/SU2_GEO.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  
  unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  char config_file_name[MAX_STRING_SIZE];
  int rank, size;

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
  
  CConfig **config_container          = nullptr;
  CGeometry **geometry_container      = nullptr;
  CGeometryEvaluation *geo_eval       = nullptr;
  
  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/
  
  if (argc == 2) { strcpy(config_file_name,argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation)  ---*/

  CConfig *config = nullptr;
  config = new CConfig(config_file_name, SU2_GEO);

  nZone    = config->GetnZone();

  /*--- Definition of the containers per zones ---*/
  
  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry*[nZone];
  
  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]       = nullptr;
    geometry_container[iZone]     = nullptr;
  }
  
  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/
    
    config_container[iZone] = new CConfig(config_file_name, SU2_GEO, true);
    config_container[iZone]->SetMPICommunicator(MPICommunicator);
    
    /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/
    
    CGeometry *geometry_aux = nullptr;
    
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

  StartTime = SU2_MPI::Wtime();

  /*--- Evaluation of the objective function ---*/
  
  if (rank == MASTER_NODE)
    cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
  
  /*--- Compute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Setting local point connectivity." <<endl;
  geometry_container[ZONE_0]->SetPoint_Connectivity();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
  if (config_container[ZONE_0]->GetReorientElements()) {
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation of the interior elements." <<endl;
    geometry_container[ZONE_0]->Check_IntElem_Orientation(config_container[ZONE_0]);
  }
  
  /*--- Create the edge structure ---*/
  
  if (rank == MASTER_NODE) cout << "Identify edges and vertices." <<endl;
  geometry_container[ZONE_0]->SetEdges(); geometry_container[ZONE_0]->SetVertex(config_container[ZONE_0]);
  
  /*--- Compute center of gravity ---*/
  
  if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
  geometry_container[ZONE_0]->SetCoord_CG();
  
  /*--- Create the dual control volume structures ---*/
  
  if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
  geometry_container[ZONE_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);
  
  /*--- Compute the surface curvature ---*/
  
  if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
  geometry_container[ZONE_0]->ComputeSurf_Curvature(config_container[ZONE_0]);
  
  /*--- Computation of positive surface area in the z-plane ---*/
  
  if (rank == MASTER_NODE) cout << "Setting reference area and span." << endl;
  geometry_container[ZONE_0]->SetPositive_ZArea(config_container[ZONE_0]);
  
  /*--- Create the point-to-point MPI communication structures. ---*/
  geometry_container[ZONE_0]->PreprocessP2PComms(geometry_container[ZONE_0], config_container[ZONE_0]);

  geo_eval = new CGeometryEvaluation(geometry_container[ZONE_0], config_container[ZONE_0], MPICommunicator);

  /*--- Set the number of sections, and allocate the memory ---*/

  if (geometry->GetnDim() == 2) nPlane = 1;
  else nPlane = config->GetnLocationStations();

  geo_eval->SetSectioningVariables(nPlane);

  geo_eval->ComputeGeometry();

  geo_eval->EvaluateObjectiveFunction();

  geo_eval->OutputFunctionFile();

  geo_eval->EverythingGradient();
		
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;

  delete geo_eval;
  geo_eval = nullptr;

  delete config;
  config = nullptr;

  if (geometry_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (geometry_container[iZone] != nullptr) {
        delete geometry_container[iZone];
      }
    }
    delete [] geometry_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;
  
  if (config_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != nullptr) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;
  
  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/
  
  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\n\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }
  
  /*--- Exit the solver cleanly ---*/
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_GEO) ------------------------" << endl << endl;
  
  
  /*--- Finalize MPI parallelization ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}
