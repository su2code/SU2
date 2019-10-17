/*!
 * \file SU2_DEF.cpp
 * \brief Main file of Mesh Deformation Code (SU2_DEF).
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

#include "../include/SU2_DEF.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  
  unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  char config_file_name[MAX_STRING_SIZE];
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
  
  /*--- Pointer to different structures that will be used throughout
   the entire code ---*/
  
  CConfig **config_container          = NULL;
  CGeometry **geometry_container      = NULL;
  CSurfaceMovement **surface_movement = NULL;
  CVolumetricMovement **grid_movement = NULL;
  COutput **output                     = NULL;
  CConfig *driver_config                = NULL;

  /*--- Load in the number of zones and spatial dimensions in the mesh file
   (if no config file is specified, default.cfg is used) ---*/
  
  if (argc == 2) { strcpy(config_file_name, argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation)  ---*/

  CConfig *config = NULL;
  config = new CConfig(config_file_name, SU2_DEF);

  nZone    = config->GetnZone();

  /*--- Definition of the containers per zones ---*/
  
  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry*[nZone];
  surface_movement   = new CSurfaceMovement*[nZone];
  grid_movement      = new CVolumetricMovement*[nZone];
  output             = new COutput*[nZone];
  
  driver_config       = NULL;

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
    surface_movement[iZone]       = NULL;
    grid_movement[iZone]          = NULL;
    output[iZone]                 = NULL;
  }
  
  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_DEF, false);

  /*--- Initialize a char to store the zone filename ---*/
  char zone_file_name[MAX_STRING_SIZE];

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/
    
    if (driver_config->GetnConfigFiles() > 0){
      strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
      config_container[iZone] = new CConfig(driver_config, zone_file_name, SU2_DEF, iZone, nZone, true);
    }
    else{
      config_container[iZone] = new CConfig(driver_config, config_file_name, SU2_DEF, iZone, nZone, true);
    }
    config_container[iZone]->SetMPICommunicator(MPICommunicator);
  }
  
  /*--- Set the multizone part of the problem. ---*/
  if (driver_config->GetMultizone_Problem()){
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config_container[iZone]->SetMultizone(driver_config, config_container);
    }
  }
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
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
  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Computational grid preprocesing ---*/

    if (rank == MASTER_NODE) cout << endl << "----------------------- Preprocessing computations ----------------------" << endl;

    /*--- Compute elements surrounding points, points surrounding points ---*/

    if (rank == MASTER_NODE) cout << "Setting local point connectivity." <<endl;
    geometry_container[iZone]->SetPoint_Connectivity();

    /*--- Check the orientation before computing geometrical quantities ---*/

    geometry_container[iZone]->SetBoundVolume();
    if (config_container[iZone]->GetReorientElements()) {
      if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation of the interior elements." <<endl;
      geometry_container[iZone]->Check_IntElem_Orientation(config_container[iZone]);
      geometry_container[iZone]->Check_BoundElem_Orientation(config_container[iZone]);
    }

    /*--- Create the edge structure ---*/

    if (rank == MASTER_NODE) cout << "Identify edges and vertices." <<endl;
    geometry_container[iZone]->SetEdges(); geometry_container[iZone]->SetVertex(config_container[iZone]);

    if (config_container[iZone]->GetDesign_Variable(0) != NO_DEFORMATION) {
      
      /*--- Compute center of gravity ---*/
      
      if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
      geometry_container[iZone]->SetCoord_CG();
      
      /*--- Create the dual control volume structures ---*/
      
      if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
      geometry_container[iZone]->SetBoundControlVolume(config_container[iZone], ALLOCATE);
      
    }
    /*--- Create the point-to-point MPI communication structures. ---*/
    
    geometry_container[iZone]->PreprocessP2PComms(geometry_container[iZone], config_container[iZone]);
    
    /*--- Allocate the mesh output ---*/
    
    output[iZone] = new CMeshOutput(config_container[iZone], geometry_container[iZone]->GetnDim());
    
    /*--- Preprocess the volume output ---*/
    
    output[iZone]->PreprocessVolumeOutput(config_container[iZone]);
    
    /*--- Preprocess history --- */
    
    output[iZone]->PreprocessHistoryOutput(config_container[iZone], false);
    

  }
  

  /*--- Output original grid for visualization, if requested (surface and volumetric) ---*/
  
  if ((config_container[ZONE_0]->GetVisualize_Volume_Def() ||
       config_container[ZONE_0]->GetVisualize_Surface_Def()) &&
      config_container[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION) {
    
    for (iZone = 0; iZone < nZone; iZone++){
      
//      /*--- Load the data --- */
      
//      output[iZone]->Load_Data(geometry_container[iZone], config_container[iZone], NULL);
      
//      if (config_container[iZone]->GetVisualize_Volume_Def()){
        
//        /*--- If requested, write the volume output for visualization purposes --- */
        
//        output[iZone]->SetVolume_Output(geometry_container[iZone], config_container[iZone], config->GetOutput_FileFormat(), false);
      
//      } 
      
//      if (config_container[iZone]->GetVisualize_Surface_Def()){
        
//        /*--- If requested, write the volume output for visualization purposes --- */
        
//        output[iZone]->SetSurface_Output(geometry_container[iZone], config_container[iZone], config->GetOutput_FileFormat(), false);
        
//      }
      
//      output[iZone]->DeallocateData_Parallel();
      
    }
  }
  
  /*--- Surface grid deformation using design variables ---*/
  
  for (iZone = 0; iZone < nZone; iZone++){
    
    if (config_container[iZone]->GetDesign_Variable(0) != NO_DEFORMATION) {
      
      /*--- Definition of the Class for grid movement ---*/
      grid_movement[iZone] = new CVolumetricMovement(geometry_container[iZone], config_container[iZone]);
      
      /*--- First check for volumetric grid deformation/transformations ---*/
      
      if (config_container[iZone]->GetDesign_Variable(0) == SCALE_GRID) {
        
        if (rank == MASTER_NODE)
          cout << endl << "--------------------- Volumetric grid scaling (ZONE " << iZone <<") ------------------" << endl;
        grid_movement[iZone]->SetVolume_Scaling(geometry_container[iZone], config_container[iZone], false);
        
      } else if (config_container[iZone]->GetDesign_Variable(0) == TRANSLATE_GRID) {
        
        if (rank == MASTER_NODE)
          cout << endl << "------------------- Volumetric grid translation (ZONE " << iZone <<") ----------------" << endl;
        grid_movement[iZone]->SetVolume_Translation(geometry_container[iZone], config_container[iZone], false);
        
      } else if (config_container[iZone]->GetDesign_Variable(0) == ROTATE_GRID) {
        
        if (rank == MASTER_NODE)
          cout << endl << "--------------------- Volumetric grid rotation (ZONE " << iZone <<") -----------------" << endl;
        grid_movement[iZone]->SetVolume_Rotation(geometry_container[iZone], config_container[iZone], false);
        
      } else {
        
        /*--- If no volume-type deformations are requested, then this is a
         surface-based deformation or FFD set up. ---*/
        
        if (rank == MASTER_NODE)
          cout << endl << "--------------------- Surface grid deformation (ZONE " << iZone <<") -----------------" << endl;

        /*--- Definition and initialization of the surface deformation class ---*/

        surface_movement[iZone] = new CSurfaceMovement();

        /*--- Copy coordinates to the surface structure ---*/

        surface_movement[iZone]->CopyBoundary(geometry_container[iZone], config_container[iZone]);

        /*--- Surface grid deformation ---*/

        if (rank == MASTER_NODE) cout << "Performing the deformation of the surface grid." << endl;
        surface_movement[iZone]->SetSurface_Deformation(geometry_container[iZone], config_container[iZone]);

        if (config_container[iZone]->GetDesign_Variable(0) != FFD_SETTING) {

          if (rank == MASTER_NODE)
            cout << endl << "------------------- Volumetric grid deformation (ZONE " << iZone <<") ----------------" << endl;

          if (rank == MASTER_NODE)
            cout << "Performing the deformation of the volumetric grid." << endl;
          grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone], config_container[iZone], false);
          
        }
        
      }
      
    }
    
  }
  
  /*--- Computational grid preprocesing ---*/
  
  if (rank == MASTER_NODE) cout << endl << "----------------------- Write deformed grid files -----------------------" << endl;
  
  /*--- Output deformed grid for visualization, if requested (surface and volumetric), in parallel
   requires to move all the data to the master node---*/

  for (iZone = 0; iZone < nZone; iZone++){
    
    /*--- Load the data --- */
    
    output[iZone]->Load_Data(geometry_container[iZone], config_container[iZone], NULL);
    
    output[iZone]->WriteToFile(config_container[iZone], geometry_container[iZone], MESH, config->GetMesh_Out_FileName());
    
    /*--- Set the file names for the visualization files ---*/
    
    output[iZone]->SetVolume_Filename("volume_deformed");
    output[iZone]->SetSurface_Filename("surface_deformed");
    
    if (config_container[iZone]->GetVisualize_Volume_Def()){
      for (unsigned short iFile = 0; iFile < config_container[iZone]->GetnVolumeOutputFiles(); iFile++){
        unsigned short* FileFormat = config_container[iZone]->GetVolumeOutputFiles();
        if (FileFormat[iFile] != RESTART_ASCII && FileFormat[iFile] != RESTART_BINARY)
          output[iZone]->WriteToFile(config_container[iZone], geometry_container[iZone], FileFormat[iFile]);
      }
    } 
       
  }

  
  if ((config_container[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION) &&
      (config_container[ZONE_0]->GetDesign_Variable(0) != SCALE_GRID)     &&
      (config_container[ZONE_0]->GetDesign_Variable(0) != TRANSLATE_GRID) &&
      (config_container[ZONE_0]->GetDesign_Variable(0) != ROTATE_GRID)) {
  
    /*--- Write the the free-form deformation boxes after deformation. ---*/
    
    if (rank == MASTER_NODE) cout << "Adding any FFD information to the SU2 file." << endl;
    
    surface_movement[ZONE_0]->WriteFFDInfo(surface_movement, geometry_container, config_container);
    
  }
  
  delete config;
  config = NULL;
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
  
  if (surface_movement != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (surface_movement[iZone] != NULL) {
        delete surface_movement[iZone];
      }
    }
    delete [] surface_movement;
  }
  if (rank == MASTER_NODE) cout << "Deleted CSurfaceMovement class." << endl;
  
  if (grid_movement != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (grid_movement[iZone] != NULL) {
        delete grid_movement[iZone];
      }
    }
    delete [] grid_movement;
  }
  if (rank == MASTER_NODE) cout << "Deleted CVolumetricMovement class." << endl;
  
  if (config_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != NULL) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (output != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (output[iZone] != NULL) {
        delete output[iZone];
      }
    }
    delete [] output;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;
  
  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

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
  
  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Exit Success (SU2_DEF) ------------------------" << endl << endl;
  
  /*--- Finalize MPI parallelization ---*/

#ifdef HAVE_MPI
  SU2_MPI::Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}
