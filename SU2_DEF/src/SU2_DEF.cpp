/*!
 * \file SU2_DEF.cpp
 * \brief Main file of Mesh Deformation Code (SU2_DEF).
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#if defined(HAVE_OMP) && defined(HAVE_MPI)
  int provided;
  SU2_MPI::Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  SU2_MPI::Init(&argc, &argv);
#endif
  SU2_MPI::Comm MPICommunicator = SU2_MPI::GetComm();

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  /*--- Pointer to different structures that will be used throughout
   the entire code ---*/

  CConfig **config_container          = nullptr;
  CGeometry **geometry_container      = nullptr;
  CSurfaceMovement **surface_movement = nullptr;
  CVolumetricMovement **grid_movement = nullptr;
  COutput **output                     = nullptr;
  CConfig *driver_config                = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file
   (if no config file is specified, default.cfg is used) ---*/

  if (argc == 2) { strcpy(config_file_name, argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation)  ---*/

  CConfig *config = nullptr;
  config = new CConfig(config_file_name, SU2_COMPONENT::SU2_DEF);

  nZone    = config->GetnZone();

  /*--- Definition of the containers per zones ---*/

  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry*[nZone];
  surface_movement   = new CSurfaceMovement*[nZone];
  grid_movement      = new CVolumetricMovement*[nZone];
  output             = new COutput*[nZone];

  driver_config       = nullptr;

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]       = nullptr;
    geometry_container[iZone]     = nullptr;
    surface_movement[iZone]       = nullptr;
    grid_movement[iZone]          = nullptr;
    output[iZone]                 = nullptr;
  }

  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_DEF, false);

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
      config_container[iZone] = new CConfig(driver_config, zone_file_name, SU2_COMPONENT::SU2_DEF, iZone, nZone, true);
    }
    else{
      config_container[iZone] = new CConfig(driver_config, config_file_name, SU2_COMPONENT::SU2_DEF, iZone, nZone, true);
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
    geometry_container[iZone]->SetEdges();
    geometry_container[iZone]->SetVertex(config_container[iZone]);

    if (config_container[iZone]->GetDesign_Variable(0) != NO_DEFORMATION) {

      /*--- Create the dual control volume structures ---*/

      if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
      geometry_container[iZone]->SetControlVolume(config_container[iZone], ALLOCATE);
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


  /*--- Surface grid deformation using design variables ---*/

  for (iZone = 0; iZone < nZone; iZone++){

    if (config_container[iZone]->GetDesign_Variable(0) != NO_DEFORMATION) {

      /*--- Definition of the Class for grid movement ---*/
      grid_movement[iZone] = new CVolumetricMovement(geometry_container[iZone], config_container[iZone]);

      /*--- Save original coordinates to be reused in convexity checking procedure ---*/
      auto OriginalCoordinates = geometry_container[iZone]->nodes->GetCoord();

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
        auto TotalDeformation = surface_movement[iZone]->SetSurface_Deformation(geometry_container[iZone], config_container[iZone]);

        if (config_container[iZone]->GetDesign_Variable(0) != FFD_SETTING) {

          if (rank == MASTER_NODE)
            cout << endl << "------------------- Volumetric grid deformation (ZONE " << iZone <<") ----------------" << endl;

          if (rank == MASTER_NODE)
            cout << "Performing the deformation of the volumetric grid." << endl;
          grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone], config_container[iZone], false);

          /*--- Get parameters for convexity check ---*/
          bool ConvexityCheck;
          unsigned short ConvexityCheck_MaxIter, ConvexityCheck_MaxDepth;

          tie(ConvexityCheck, ConvexityCheck_MaxIter, ConvexityCheck_MaxDepth) = config_container[iZone]->GetConvexityCheck();

          /*--- Recursively change deformations if there are nonconvex elements. ---*/

          if (ConvexityCheck && geometry_container[iZone]->GetnNonconvexElements() > 0) {
            if (rank == MASTER_NODE) {
              cout << "Nonconvex elements present after deformation. " << endl;
              cout << "Recursively lowering deformation magnitude." << endl;
            }

            /*--- Load initial deformation values ---*/
            auto InitialDeformation = TotalDeformation;

            unsigned short ConvexityCheckIter, RecursionDepth = 0;
            su2double DeformationFactor = 1.0, DeformationDifference = 1.0;
            for (ConvexityCheckIter = 1; ConvexityCheckIter <= ConvexityCheck_MaxIter; ConvexityCheckIter++) {

              /*--- Recursively change deformation magnitude:
              decrease if there are nonconvex elements, increase otherwise ---*/
              DeformationDifference /= 2.0;

              if (geometry_container[iZone]->GetnNonconvexElements() > 0) {
                DeformationFactor -= DeformationDifference;
              } else {
                RecursionDepth += 1;

                if (RecursionDepth == ConvexityCheck_MaxDepth) {
                  if (rank == MASTER_NODE) {
                    cout << "Maximum recursion depth reached." << endl;
                    cout << "Remaining amount of original deformation: ";
                    cout << DeformationFactor*100.0 << " percent. " << endl;
                  }
                  break;
                }

                DeformationFactor += DeformationDifference;
              }

              /*--- Load mesh to start every iteration with an undeformed grid ---*/
              for (auto iPoint = 0ul; iPoint < OriginalCoordinates.rows(); iPoint++) {
                for (auto iDim = 0ul; iDim < OriginalCoordinates.cols(); iDim++) {
                  geometry_container[iZone]->nodes->SetCoord(iPoint, iDim, OriginalCoordinates(iPoint,iDim));
                }
              }

              /*--- Set deformation magnitude as percentage of initial deformation ---*/
              for (auto iDV = 0u; iDV < config->GetnDV(); iDV++) {
                for (auto iDV_Value = 0u; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {
                  config_container[iZone]->SetDV_Value(iDV, iDV_Value, InitialDeformation[iDV][iDV_Value]*DeformationFactor);
                }
              }

              /*--- Surface grid deformation ---*/
              if (rank == MASTER_NODE) cout << "Performing the deformation of the surface grid." << endl;

              TotalDeformation = surface_movement[iZone]->SetSurface_Deformation(geometry_container[iZone], config_container[iZone]);

              if (rank == MASTER_NODE)
                cout << endl << "------------------- Volumetric grid deformation (ZONE " << iZone <<") ----------------" << endl;

              if (rank == MASTER_NODE)
                cout << "Performing the deformation of the volumetric grid." << endl;
              grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone], config_container[iZone], false);

              if (rank == MASTER_NODE) {
                cout << "Number of nonconvex elements for iteration " << ConvexityCheckIter << ": ";
                cout << geometry_container[iZone]->GetnNonconvexElements() << endl;
                cout << "Remaining amount of original deformation: ";
                cout << DeformationFactor*100.0 << " percent. " << endl;
              }

            }

          }

        }

      }

    }

  }

  /*--- Computational grid preprocesing ---*/

  if (rank == MASTER_NODE) cout << endl << "----------------------- Write deformed grid files -----------------------" << endl;

  /*--- Output deformed grid for visualization, if requested (surface and volumetric), in parallel
   requires to move all the data to the master node---*/

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Compute Mesh Quality if requested. Necessary geometry preprocessing re-done beforehand. ---*/

    if (config_container[iZone]->GetWrt_MeshQuality() && !config->GetStructuralProblem()) {

      if (rank == MASTER_NODE) cout << "Recompute geometry properties necessary to evaluate mesh quality statistics.\n";

      geometry_container[iZone]->SetPoint_Connectivity();
      geometry_container[iZone]->SetBoundVolume();
      geometry_container[iZone]->SetEdges();
      geometry_container[iZone]->SetVertex(config_container[iZone]);
      geometry_container[iZone]->SetControlVolume(config_container[iZone], ALLOCATE);
      geometry_container[iZone]->SetBoundControlVolume(config_container[iZone], ALLOCATE);

      if (rank == MASTER_NODE) cout << "Computing mesh quality statistics for the dual control volumes.\n";
      geometry_container[iZone]->ComputeMeshQualityStatistics(config_container[iZone]);
    }// Mesh Quality Output

    /*--- Load the data --- */

    output[iZone]->Load_Data(geometry_container[iZone], config_container[iZone], nullptr);

    output[iZone]->WriteToFile(config_container[iZone], geometry_container[iZone], OUTPUT_TYPE::MESH, config->GetMesh_Out_FileName());

    /*--- Set the file names for the visualization files ---*/

    output[iZone]->SetVolume_Filename("volume_deformed");
    output[iZone]->SetSurface_Filename("surface_deformed");

    for (unsigned short iFile = 0; iFile < config_container[iZone]->GetnVolumeOutputFiles(); iFile++){
      auto FileFormat = config_container[iZone]->GetVolumeOutputFiles();
      if (FileFormat[iFile] != OUTPUT_TYPE::RESTART_ASCII && FileFormat[iFile] != OUTPUT_TYPE::RESTART_BINARY)
        output[iZone]->WriteToFile(config_container[iZone], geometry_container[iZone], FileFormat[iFile]);
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
  config = nullptr;
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;

  if (geometry_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (geometry_container[iZone] != nullptr) {
        delete geometry_container[iZone];
      }
    }
    delete [] geometry_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

  if (surface_movement != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (surface_movement[iZone] != nullptr) {
        delete surface_movement[iZone];
      }
    }
    delete [] surface_movement;
  }
  if (rank == MASTER_NODE) cout << "Deleted CSurfaceMovement class." << endl;

  if (grid_movement != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (grid_movement[iZone] != nullptr) {
        delete grid_movement[iZone];
      }
    }
    delete [] grid_movement;
  }
  if (rank == MASTER_NODE) cout << "Deleted CVolumetricMovement class." << endl;

  if (config_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != nullptr) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (output != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (output[iZone] != nullptr) {
        delete output[iZone];
      }
    }
    delete [] output;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

  StopTime = SU2_MPI::Wtime();

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
  SU2_MPI::Finalize();

  return EXIT_SUCCESS;

}
