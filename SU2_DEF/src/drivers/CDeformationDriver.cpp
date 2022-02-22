/*!
 * \file CDeformationDriver.hpp
 * \brief Main subroutines for driving the mesh deformation.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser/ General Public
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


#include "../../include/drivers/CDeformationDriver.hpp"

#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../SU2_CFD/include/solvers/CMeshSolver.hpp"
#include "../../../SU2_CFD/include/output/CMeshOutput.hpp"
#include "../../../SU2_CFD/include/numerics/elasticity/CFEALinearElasticity.hpp"

using namespace std;

CDeformationDriver::CDeformationDriver(char* confFile, SU2_Comm MPICommunicator) {

  /*--- Initialize Medipack (must also be here so it is initialized from python) ---*/
  #ifdef HAVE_MPI
    #if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
      SU2_MPI::Init_AMPI();
    #endif
  #endif

  SU2_MPI::SetComm(MPICommunicator);

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  /*--- Copy the config filename ---*/
  strcpy(config_file_name, confFile);

  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_DEF);

  nZone = driver_config->GetnZone();

  /*--- Initialize containers --- */

  SetContainers_Null();

  /*--- Preprocessing of the config files. ---*/

  Input_Preprocessing();

  /*--- Set up a timer for performance benchmarking ---*/

  StartTime = SU2_MPI::Wtime();

  /*--- Preprocessing of the geometry for all zones. ---*/

  Geometrical_Preprocessing();

  /*--- Preprocessing of the output for all zones. ---*/

  Output_Preprocessing();

  if (driver_config->GetDeform_Mesh()){

    /*--- Preprocessing of the mesh solver for all zones. ---*/

    Solver_Preprocessing();

    /*--- Preprocessing of the mesh solver for all zones. ---*/

    Numerics_Preprocessing();

  }

  /*--- Preprocessing time is reported now, but not included in the next compute portion. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  UsedTimePreproc    = UsedTime;
  UsedTimeCompute    = 0.0;

}

CDeformationDriver::~CDeformationDriver(void) {

}

void CDeformationDriver::SetContainers_Null() {

    /*--- Create pointers to all of the classes that may be used throughout
        the SU2_DEF code. In general, the pointers are instantiated down a
        hierarchy over all zones as described in the comments below. ---*/
    config_container   = new CConfig*[nZone];
    output_container   = new COutput*[nZone];
    geometry_container = new CGeometry*[nZone];
    surface_movement   = new CSurfaceMovement*[nZone];
    grid_movement      = new CVolumetricMovement*[nZone];

    solver_container   = new CSolver*[nZone];
    numerics_container = new CNumerics**[nZone];

    for (iZone = 0; iZone < nZone; iZone++) {
        config_container[iZone]       = nullptr;
        output_container[iZone]       = nullptr;
        geometry_container[iZone]     = nullptr;
        surface_movement[iZone]       = nullptr;
        grid_movement[iZone]          = nullptr;
        solver_container[iZone]       = nullptr;
        numerics_container[iZone]     = nullptr;
    }
}

void CDeformationDriver::Input_Preprocessing() {

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
        } else {
            config_container[iZone] = new CConfig(driver_config, config_file_name, SU2_COMPONENT::SU2_DEF, iZone, nZone, true);
        }

        config_container[iZone]->SetMPICommunicator(SU2_MPI::GetComm());
    }

    /*--- Set the multizone part of the problem. ---*/

    if (driver_config->GetMultizone_Problem()){
        for (iZone = 0; iZone < nZone; iZone++) {

            /*--- Set the interface markers for multizone ---*/

            config_container[iZone]->SetMultizone(driver_config, config_container);
        }
    }
}

void CDeformationDriver::Geometrical_Preprocessing() {

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

    }
}

void CDeformationDriver::Output_Preprocessing() {

      for (iZone = 0; iZone < nZone; iZone++) {

        /*--- Allocate the mesh output ---*/

        output_container[iZone] = new CMeshOutput(config_container[iZone], geometry_container[iZone]->GetnDim());

        /*--- Preprocess the volume output ---*/

        output_container[iZone]->PreprocessVolumeOutput(config_container[iZone]);

        /*--- Preprocess history --- */

        output_container[iZone]->PreprocessHistoryOutput(config_container[iZone], false);

    }
}

void CDeformationDriver::Solver_Preprocessing() {

  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone] = new CMeshSolver(geometry_container[iZone], config_container[iZone]);
  }

}

void CDeformationDriver::Numerics_Preprocessing() {

  for (iZone = 0; iZone < nZone; iZone++) {
    numerics_container[iZone] = new CNumerics* [omp_get_num_threads() * MAX_TERMS]();

    for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
        const int iTerm = FEA_TERM + thread * MAX_TERMS;
        const int nDim = geometry_container[iZone]->GetnDim();

        numerics_container[iZone][iTerm] = new CFEAMeshElasticity(nDim, nDim, geometry_container[iZone]->GetnElem(), config_container[iZone]);
    }

  }

}

void CDeformationDriver::Run() {

    /* --- Start measuring computation time ---*/

    StartTime = SU2_MPI::Wtime();

    /*--- Surface grid deformation using design variables ---*/

    if (driver_config->GetDeform_Mesh()) {
        Update();
    }
    else {
        Update_Legacy();
    }

    /*--- Synchronization point after a single solver iteration. Compute the
      wall clock time required. ---*/

    StopTime = SU2_MPI::Wtime();

    UsedTimeCompute = StopTime-StartTime;
    if (rank == MASTER_NODE) {
      cout << "\nCompleted in " << fixed << UsedTimeCompute << " seconds on "<< size;

      if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
    }

    /*--- Output the deformed mesh ---*/
    Output();

}

void CDeformationDriver::Update() {

    for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the stiffness of each element mesh into the mesh numerics ---*/

    solver_container[iZone]->SetMesh_Stiffness(numerics_container[iZone], config_container[iZone]);

    /*--- Deform the volume grid around the new boundary locations ---*/

    solver_container[iZone]->DeformMesh(geometry_container[iZone], numerics_container[iZone], config_container[iZone]);

    }
}

void CDeformationDriver::Update_Legacy() {

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
                for (auto iDV = 0u; iDV < driver_config->GetnDV(); iDV++) {
                  for (auto iDV_Value = 0u; iDV_Value < driver_config->GetnDV_Value(iDV); iDV_Value++) {
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

}

void CDeformationDriver::Output() {

    /*--- Output deformed grid for visualization, if requested (surface and volumetric), in parallel
          requires to move all the data to the master node---*/

    if (rank == MASTER_NODE) cout << endl << "----------------------- Write deformed grid files -----------------------" << endl;

    for (iZone = 0; iZone < nZone; iZone++){

      /*--- Compute Mesh Quality if requested. Necessary geometry preprocessing re-done beforehand. ---*/

      if (config_container[iZone]->GetWrt_MeshQuality() && !driver_config->GetStructuralProblem()) {

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

      output_container[iZone]->Load_Data(geometry_container[iZone], config_container[iZone], nullptr);

      output_container[iZone]->WriteToFile(config_container[iZone], geometry_container[iZone], OUTPUT_TYPE::MESH, driver_config->GetMesh_Out_FileName());

      /*--- Set the file names for the visualization files ---*/

      output_container[iZone]->SetVolume_Filename("volume_deformed");
      output_container[iZone]->SetSurface_Filename("surface_deformed");

      for (unsigned short iFile = 0; iFile < config_container[iZone]->GetnVolumeOutputFiles(); iFile++){
        auto FileFormat = config_container[iZone]->GetVolumeOutputFiles();
        if (FileFormat[iFile] != OUTPUT_TYPE::RESTART_ASCII && FileFormat[iFile] != OUTPUT_TYPE::RESTART_BINARY)
          output_container[iZone]->WriteToFile(config_container[iZone], geometry_container[iZone], FileFormat[iFile]);
      }
    }

    if (!driver_config->GetDeform_Mesh()) {
      if ((config_container[ZONE_0]->GetDesign_Variable(0) != NO_DEFORMATION) &&
          (config_container[ZONE_0]->GetDesign_Variable(0) != SCALE_GRID)     &&
          (config_container[ZONE_0]->GetDesign_Variable(0) != TRANSLATE_GRID) &&
          (config_container[ZONE_0]->GetDesign_Variable(0) != ROTATE_GRID)) {

          /*--- Write the the free-form deformation boxes after deformation. ---*/

          if (rank == MASTER_NODE) cout << "Adding any FFD information to the SU2 file." << endl;

          surface_movement[ZONE_0]->WriteFFDInfo(surface_movement, geometry_container, config_container);

    }
  }
}

void CDeformationDriver::Postprocessing() {

    if (rank == MASTER_NODE)
      cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;

    delete driver_config;
    driver_config = nullptr;

    for (iZone = 0; iZone < nZone; iZone++) {
      if (numerics_container[iZone] != nullptr) {
        for (unsigned int iTerm = 0; iTerm < MAX_TERMS*omp_get_max_threads(); iTerm++) {
          delete numerics_container[iZone][iTerm];
        }
        delete [] numerics_container[iZone];
      }
    }
    delete [] numerics_container;
    if (rank == MASTER_NODE) cout << "Deleted CNumerics container." << endl;

    for (iZone = 0; iZone < nZone; iZone++) {
      delete solver_container[iZone];
    }
    delete [] solver_container;
    if (rank == MASTER_NODE) cout << "Deleted CSolver container." << endl;

    if (geometry_container != nullptr) {
      for (iZone = 0; iZone < nZone; iZone++) {
        delete geometry_container[iZone];
      }
      delete [] geometry_container;
    }
    if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

    if (surface_movement != nullptr) {
      for (iZone = 0; iZone < nZone; iZone++) {
        delete surface_movement[iZone];
      }
      delete [] surface_movement;
    }
    if (rank == MASTER_NODE) cout << "Deleted CSurfaceMovement class." << endl;

    if (grid_movement != nullptr) {
      for (iZone = 0; iZone < nZone; iZone++) {
        delete grid_movement[iZone];
      }
      delete [] grid_movement;
    }
    if (rank == MASTER_NODE) cout << "Deleted CVolumetricMovement class." << endl;

    if (config_container != nullptr) {
      for (iZone = 0; iZone < nZone; iZone++) {
        delete config_container[iZone];
      }
      delete [] config_container;
    }
    if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

    if (output_container != nullptr) {
      for (iZone = 0; iZone < nZone; iZone++) {
        delete output_container[iZone];
      }
      delete [] output_container;
    }
    if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

    /*--- Exit the solver cleanly ---*/

    if (rank == MASTER_NODE)
      cout << endl << "------------------------- Exit Success (SU2_DEF) ------------------------" << endl << endl;
}

vector<string> CDeformationDriver::GetAllDeformMeshMarkersTag() const {

  vector<string> interfaceBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_Deform_Mesh();
  interfaceBoundariesTagList.resize(nBoundariesMarker);

  for(iMarker=0; iMarker < nBoundariesMarker; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_Deform_Mesh_TagBound(iMarker);
    interfaceBoundariesTagList[iMarker] = Marker_Tag;
  }

  return interfaceBoundariesTagList;
}

map<string, int> CDeformationDriver::GetAllBoundaryMarkers() const {

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

map<string, string> CDeformationDriver::GetAllBoundaryMarkersType() const {

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

unsigned long CDeformationDriver::GetNumberDimensions() const {
    return geometry_container[ZONE_0]->GetnDim();
}

unsigned long CDeformationDriver::GetNumberElements() const {
    return geometry_container[ZONE_0]->GetnElem();
}

unsigned long CDeformationDriver::GetNumberElementsMarker(unsigned short iMarker) const {
    return geometry_container[ZONE_0]->GetnElem_Bound(iMarker);
}

unsigned long CDeformationDriver::GetNumberVertices() const {
    return geometry_container[ZONE_0]->GetnPoint();
}

unsigned long CDeformationDriver::GetNumberVerticesMarker(unsigned short iMarker) const {
    return geometry_container[ZONE_0]->GetnVertex(iMarker);
}

unsigned long CDeformationDriver::GetNumberHaloVertices() const {
    
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nPoint   = geometry->GetnPoint();
    
    unsigned long nHaloVertices = 0;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        if (!(geometry->nodes->GetDomain(iPoint))) {
            nHaloVertices += 1;
        }
    }
    
    return nHaloVertices;
}

unsigned long CDeformationDriver::GetNumberHaloVerticesMarker(unsigned short iMarker) const {
    
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    
    unsigned long nHaloVertices = 0;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (!(geometry->nodes->GetDomain(iPoint))) {
            nHaloVertices += 1;
        }
    }
    
    return nHaloVertices;
}

vector<unsigned long> CDeformationDriver::GetVertexIDs() const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nPoint   = geometry->GetnPoint();
    
    vector<unsigned long> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(geometry->nodes->GetGlobalIndex(iPoint));
    }
    
    return values;
}

vector<unsigned long> CDeformationDriver::GetVertexIDsMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    
    vector<unsigned long> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        values.push_back(geometry->nodes->GetGlobalIndex(iPoint));
    }
    
    return values;
}

vector<unsigned long> CDeformationDriver::GetElementIDs() const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nElem    = geometry->GetnElem();
    
    vector<unsigned long> values;
    
    for (auto iElem = 0ul; iElem < nElem; iElem++) {
        values.push_back(geometry->elem[iElem]->GetGlobalIndex());
    }
    
    return values;
}

vector<unsigned long> CDeformationDriver::GetElementIDsMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nBound   = geometry->GetnElem_Bound(iMarker);
    
    vector<unsigned long> values;
    
    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        values.push_back(geometry->bound[iMarker][iBound]->GetGlobalIndex());
    }
    
    return values;
}

vector<vector<unsigned long>> CDeformationDriver::GetConnectivity() const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nElem    = geometry->GetnElem();
    
    vector<vector<unsigned long>> values(nElem);
    
    for (auto iElem = 0ul; iElem < nElem; iElem++) {
        unsigned short nNode = geometry->elem[iElem]->GetnNodes();
        
        for (auto iNode = 0u; iNode < nNode; iNode++) {
            values[iElem].push_back(geometry->elem[iElem]->GetNode(iNode));
        }
    }
    
    return values;
}

vector<vector<unsigned long>> CDeformationDriver::GetConnectivityMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nBound   = geometry->GetnElem_Bound(iMarker);
    
    vector<vector<unsigned long>> values(nBound);
    
    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        unsigned short nNode = geometry->bound[iMarker][iBound]->GetnNodes();
        
        for (auto iNode = 0u; iNode < nNode; iNode++) {
            values[iBound].push_back(geometry->bound[iMarker][iBound]->GetNode(iNode));
        }
    }
    
    return values;
}

vector<bool> CDeformationDriver::GetDomain() const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nPoint   = geometry->GetnPoint();
    
    vector<bool> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(geometry->nodes->GetDomain(iPoint));
    }
    
    return values;
}

vector<bool> CDeformationDriver::GetDomainMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    
    vector<bool> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        values.push_back(geometry->nodes->GetDomain(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDeformationDriver::GetCoordinates() const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nPoint   = geometry->GetnPoint();
    const auto nDim     = geometry->GetnDim();
    
    vector<passivedouble> values(nPoint*nDim, 0.0);
    su2double value;

    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = geometry->nodes->GetCoord(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDeformationDriver::GetCoordinatesMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    const auto nDim     = geometry->GetnDim();
    
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = geometry->nodes->GetCoord(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDeformationDriver::SetCoordinates(vector<passivedouble> values) {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nPoint   = geometry->GetnPoint();
    
    if (values.size() != nPoint*nDim) {
        SU2_MPI::Error("Size does not match nPoint * nDim !", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            geometry->nodes->SetCoord(iPoint, iDim, values[iPoint*nDim + iDim]);
        }
    }
}

void CDeformationDriver::SetCoordinatesMarker(unsigned short iMarker, vector<passivedouble> values) {
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    const auto nDim     = geometry->GetnDim();
    
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim !", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            geometry->nodes->SetCoord(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

vector<passivedouble> CDeformationDriver::GetDisplacementsMarker(unsigned short iMarker) const {
    CConfig* config = config_container[ZONE_0];
    
    if (!config->GetDeform_Mesh()) {
        return {};
    }
    
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    const auto nDim     = geometry->GetnDim();
    
    CSolver* solver = solver_container[ZONE_0];
    
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetBound_Disp(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDeformationDriver::SetDisplacementsMarker(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config = config_container[ZONE_0];
    
    if (!config->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined !", CURRENT_FUNCTION);
    }
    
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    const auto nDim     = geometry->GetnDim();
    
    CSolver* solver = solver_container[ZONE_0];
    
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim !", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetBound_Disp(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

vector<passivedouble> CDeformationDriver::GetVelocitiesMarker(unsigned short iMarker) const {
    CConfig* config = config_container[ZONE_0];
    
    if (!config->GetDeform_Mesh()) {
        return {};
    }
    
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    const auto nDim     = geometry->GetnDim();
    
    CSolver* solver = solver_container[ZONE_0];
    
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetBound_Vel(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDeformationDriver::SetVelocitiesMarker(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config = config_container[ZONE_0];
    
    if (!config->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined !", CURRENT_FUNCTION);
    }
    
    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);
    const auto nDim     = geometry->GetnDim();
    
    CSolver* solver = solver_container[ZONE_0];
    
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim !", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetBound_Vel(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}


unsigned long CDeformationDriver::GetVertexGlobalIndex(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint, GlobalIndex;

  iPoint = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetNode();
  GlobalIndex = geometry_container[ZONE_0]->nodes->GetGlobalIndex(iPoint);

  return GlobalIndex;

}

bool CDeformationDriver::IsAHaloNode(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetNode();
  if(geometry_container[ZONE_0]->nodes->GetDomain(iPoint)) return false;
  else return true;

}

vector<passivedouble> CDeformationDriver::GetInitialMeshCoord(unsigned short iMarker, unsigned long iVertex) const {

  vector<su2double> coord(3,0.0);
  vector<passivedouble> coord_passive(3, 0.0);

  auto iPoint = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetNode();

  for (auto iDim = 0 ; iDim < geometry_container[ZONE_0]->GetnDim(); iDim++){
   coord[iDim] = solver_container[ZONE_0]->GetNodes()->GetMesh_Coord(iPoint,iDim);
  }

  coord_passive[0] = SU2_TYPE::GetValue(coord[0]);
  coord_passive[1] = SU2_TYPE::GetValue(coord[1]);
  coord_passive[2] = SU2_TYPE::GetValue(coord[2]);

  return coord_passive;
}

vector<passivedouble> CDeformationDriver::GetVertexNormal(unsigned short iMarker, unsigned long iVertex, bool unitNormal) const {

  int nDim = geometry_container[ZONE_0]->GetnDim();

  su2double *Normal;
  su2double Area;
  vector<su2double> ret_Normal(3, 0.0);
  vector<passivedouble> ret_Normal_passive(3, 0.0);

  Normal = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetNormal();

  if (!unitNormal) {

    ret_Normal_passive[0] = SU2_TYPE::GetValue(Normal[0]);
    ret_Normal_passive[1] = SU2_TYPE::GetValue(Normal[1]);
    if(nDim>2) ret_Normal_passive[2] = SU2_TYPE::GetValue(Normal[2]);

    return ret_Normal_passive;
  }

  Area = GeometryToolbox::Norm(nDim, Normal);

  ret_Normal[0] = Normal[0]/Area;
  ret_Normal[1] = Normal[1]/Area;
  if(nDim>2) ret_Normal[2] = Normal[2]/Area;

  ret_Normal_passive[0] = SU2_TYPE::GetValue(ret_Normal[0]);
  ret_Normal_passive[1] = SU2_TYPE::GetValue(ret_Normal[1]);
  ret_Normal_passive[2] = SU2_TYPE::GetValue(ret_Normal[2]);

  return ret_Normal_passive;
}

void CDeformationDriver::SetVertexCoordX(unsigned short iMarker, unsigned long iVertex, passivedouble newPosX) {

  CGeometry* geometry = geometry_container[ZONE_0];

  auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  geometry->nodes->SetCoord(iPoint, 0, newPosX);

}

void CDeformationDriver::SetVertexCoordY(unsigned short iMarker, unsigned long iVertex, passivedouble newPosY) {

  CGeometry* geometry = geometry_container[ZONE_0];

  auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  geometry->nodes->SetCoord(iPoint, 1, newPosY);

}

void CDeformationDriver::SetVertexCoordZ(unsigned short iMarker, unsigned long iVertex, passivedouble newPosZ) {

  CGeometry* geometry = geometry_container[ZONE_0];

  auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  geometry->nodes->SetCoord(iPoint, 2, newPosZ);

}

vector<passivedouble> CDeformationDriver::GetMeshDisplacementsMarker(unsigned short iMarker) const {
    CConfig* config = config_container[ZONE_0];

    if (!config->GetDeform_Mesh()) {
        return {};
    }

    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);

    CSolver* solver = solver_container[ZONE_0];

    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;

    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetBound_Disp(iPoint, iDim);

            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }

    return values;
}

void CDeformationDriver::SetMeshDisplacementsMarker(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config = config_container[ZONE_0];

    if (!config->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined !", CURRENT_FUNCTION);
    }

    CGeometry* geometry = geometry_container[ZONE_0];
    const auto nVertex  = geometry->GetnVertex(iMarker);

    CSolver* solver = solver_container[ZONE_0];

    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim !", CURRENT_FUNCTION);
    }

    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetBound_Disp(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

void CDeformationDriver::CommunicateMeshDisplacement(void) {

  solver_container[ZONE_0]->InitiateComms(geometry_container[ZONE_0], config_container[ZONE_0], MESH_DISPLACEMENTS);
  solver_container[ZONE_0]->CompleteComms(geometry_container[ZONE_0], config_container[ZONE_0], MESH_DISPLACEMENTS);

}
