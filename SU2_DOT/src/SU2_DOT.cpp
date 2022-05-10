/*!
 * \file SU2_DOT.cpp
 * \brief Main file of the Gradient Projection Code (SU2_DOT).
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


#include "../include/SU2_DOT.hpp"
using namespace std;

int main(int argc, char *argv[]) {

  unsigned short iZone, iInst;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;

  char config_file_name[MAX_STRING_SIZE];

  /*--- OpenMP initialization ---*/

  omp_initialize();

  /*--- MPI initialization, and buffer setting ---*/

#if defined(HAVE_OMP) && defined(HAVE_MPI)
  int provided;
  SU2_MPI::Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  SU2_MPI::Init(&argc, &argv);
#endif
  SU2_MPI::Comm MPICommunicator = SU2_MPI::GetComm();

  const int rank = SU2_MPI::GetRank();
  const int size = SU2_MPI::GetSize();

  /*--- AD initialization ---*/
#ifdef HAVE_OPDI
  AD::getGlobalTape().initialize();
#endif

  /*--- Pointer to different structures that will be used throughout the entire code ---*/

  CConfig **config_container            = nullptr;
  CConfig *driver_config                = nullptr;
  CGeometry ***geometry_container       = nullptr;
  CSurfaceMovement **surface_movement   = nullptr;
  CVolumetricMovement **grid_movement   = nullptr;
  unsigned short *nInst                 = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/

  if (argc == 2) { strcpy(config_file_name,argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation)  ---*/

  CConfig *config = nullptr;
  config = new CConfig(config_file_name, SU2_COMPONENT::SU2_DOT);

  const auto nZone = config->GetnZone();

  /*--- Definition of the containers per zones ---*/

  config_container    = new CConfig*[nZone] ();
  geometry_container  = new CGeometry**[nZone] ();
  surface_movement    = new CSurfaceMovement*[nZone] ();
  grid_movement       = new CVolumetricMovement*[nZone] ();

  nInst               = new unsigned short[nZone];
  driver_config       = nullptr;

  for (iZone = 0; iZone < nZone; iZone++) {
    nInst[iZone] = 1;
  }

  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_DOT, false);

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
      config_container[iZone] = new CConfig(driver_config, zone_file_name, SU2_COMPONENT::SU2_DOT, iZone, nZone, true);
    }
    else{
      config_container[iZone] = new CConfig(driver_config, config_file_name, SU2_COMPONENT::SU2_DOT, iZone, nZone, true);
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

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Determine whether or not the FEM solver is used, which decides the
     type of geometry classes that are instantiated. ---*/
    const bool fem_solver = config_container[iZone]->GetFEMSolver();

    /*--- Read the number of instances for each zone ---*/

    nInst[iZone] = config_container[iZone]->GetnTimeInstances();

    geometry_container[iZone] = new CGeometry*[nInst[iZone]];

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/

      CGeometry *geometry_aux = nullptr;

      /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

      geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

      /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

      if ( fem_solver ) geometry_aux->SetColorFEMGrid_Parallel(config_container[iZone]);
      else              geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

      /*--- Build the grid data structures using the ParMETIS coloring. ---*/

      if( fem_solver ) {
        switch( config_container[iZone]->GetKind_FEM_Flow() ) {
          case DG: {
            geometry_container[iZone][iInst] = new CMeshFEM_DG(geometry_aux, config_container[iZone]);
            break;
          }
        }
      }
      else {
        geometry_container[iZone][iInst] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);
      }

      /*--- Deallocate the memory of geometry_aux ---*/

      delete geometry_aux;

      /*--- Add the Send/Receive boundaries ---*/

      geometry_container[iZone][iInst]->SetSendReceive(config_container[iZone]);

      /*--- Add the Send/Receive boundaries ---*/

      geometry_container[iZone][iInst]->SetBoundaries(config_container[iZone]);

      /*--- Create the vertex structure (required for MPI) ---*/

      if (rank == MASTER_NODE) cout << "Identify vertices." << endl;
      geometry_container[iZone][iInst]->SetVertex(config_container[iZone]);

      /*--- Store the global to local mapping after preprocessing. ---*/

      if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
      geometry_container[iZone][iInst]->SetGlobal_to_Local_Point();

      /* Test for a fem solver, because some more work must be done. */

      if (fem_solver) {

        /*--- Carry out a dynamic cast to CMeshFEM_DG, such that it is not needed to
         define all virtual functions in the base class CGeometry. ---*/
        CMeshFEM_DG *DGMesh = dynamic_cast<CMeshFEM_DG *>(geometry_container[iZone][iInst]);

        /*--- Determine the standard elements for the volume elements. ---*/
        if (rank == MASTER_NODE) cout << "Creating standard volume elements." << endl;
        DGMesh->CreateStandardVolumeElements(config_container[iZone]);

        /*--- Create the face information needed to compute the contour integral
         for the elements in the Discontinuous Galerkin formulation. ---*/
        if (rank == MASTER_NODE) cout << "Creating face information." << endl;
        DGMesh->CreateFaces(config_container[iZone]);
      }
    }
  }

  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/

  StartTime = SU2_MPI::Wtime();

  for (iZone = 0; iZone < nZone; iZone++) {

    if (rank == MASTER_NODE)
      cout << "\n----------------------- Preprocessing computations ----------------------" << endl;

    /*--- Compute elements surrounding points, points surrounding points ---*/

    if (rank == MASTER_NODE) cout << "Setting local point connectivity." << endl;
      geometry_container[iZone][INST_0]->SetPoint_Connectivity();

    /*--- Check the orientation before computing geometrical quantities ---*/

    geometry_container[iZone][INST_0]->SetBoundVolume();
    if (config_container[iZone]->GetReorientElements()) {
      if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation of the elements." << endl;
      geometry_container[iZone][INST_0]->Check_IntElem_Orientation(config_container[iZone]);
      geometry_container[iZone][INST_0]->Check_BoundElem_Orientation(config_container[iZone]);
    }

    /*--- Create the edge structure ---*/

    if (rank == MASTER_NODE) cout << "Identify edges and vertices." << endl;
      geometry_container[iZone][INST_0]->SetEdges(); geometry_container[iZone][INST_0]->SetVertex(config_container[iZone]);

    /*--- Create the dual control volume structures ---*/

    if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
    geometry_container[iZone][INST_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);

    /*--- Store the global to local mapping after preprocessing. ---*/

    if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
    geometry_container[iZone][INST_0]->SetGlobal_to_Local_Point();

    /*--- Create the point-to-point MPI communication structures. ---*/

    geometry_container[iZone][INST_0]->PreprocessP2PComms(geometry_container[iZone][INST_0], config_container[iZone]);

    /*--- Load the surface sensitivities from file. This is done only
     once: if this is an unsteady problem, a time-average of the surface
     sensitivities at each node is taken within this routine. ---*/

    if (!config_container[iZone]->GetDiscrete_Adjoint()) {
      if (rank == MASTER_NODE) cout << "Reading surface sensitivities at each node from file." << endl;
      geometry_container[iZone][INST_0]->SetBoundSensitivity(config_container[iZone]);
    }
    else {
      if (rank == MASTER_NODE) cout << "Reading volume sensitivities at each node from file." << endl;
      grid_movement[iZone] = new CVolumetricMovement(geometry_container[iZone][INST_0], config_container[iZone]);

      /*--- Read in sensitivities from file. ---*/
      if (config_container[ZONE_0]->GetSensitivity_Format() == UNORDERED_ASCII)
        geometry_container[iZone][INST_0]->ReadUnorderedSensitivity(config_container[iZone]);
      else
        geometry_container[iZone][INST_0]->SetSensitivity(config_container[iZone]);

      if (rank == MASTER_NODE)
        cout << "\n---------------------- Mesh sensitivity computation ---------------------" << endl;

      if(config_container[iZone]->GetDiscrete_Adjoint() && config_container[iZone]->GetSmoothGradient() &&
         config_container[iZone]->GetSobMode() == ENUM_SOBOLEV_MODUS::MESH_LEVEL) {
        DerivativeTreatment_MeshSensitivity(geometry_container[iZone][INST_0], config_container[iZone], grid_movement[iZone]);
      } else {
        grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone][INST_0], config_container[iZone], false, true);
      }

    }
  }


  if (config_container[ZONE_0]->GetDiscrete_Adjoint()) {
    if (rank == MASTER_NODE)
      cout << "\n------------------------ Mesh sensitivity Output ------------------------" << endl;
    SetSensitivity_Files(geometry_container, config_container, nZone);
  }

  /*--- Initialize structure to store the gradient ---*/
  su2double** Gradient = new su2double*[config_container[ZONE_0]->GetnDV()];

  for (auto iDV = 0u; iDV  < config_container[ZONE_0]->GetnDV(); iDV++) {
    /*--- Initialize to zero ---*/
    Gradient[iDV] = new su2double[config_container[ZONE_0]->GetnDV_Value(iDV)]();
  }

  ofstream Gradient_file;
  Gradient_file.precision(config->OptionIsSet("OUTPUT_PRECISION") ? config->GetOutput_Precision() : 6);

  /*--- For multizone computations the gradient contributions are summed up and written into one file. ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    if ((config_container[iZone]->GetDesign_Variable(0) != NONE) &&
        (config_container[iZone]->GetDesign_Variable(0) != SURFACE_FILE)) {

      if (rank == MASTER_NODE)
        cout << "\n---------- Start gradient evaluation using sensitivity information ----------" << endl;

      /*--- Definition of the Class for surface deformation ---*/

      surface_movement[iZone] = new CSurfaceMovement();

      /*--- Copy coordinates to the surface structure ---*/

      surface_movement[iZone]->CopyBoundary(geometry_container[iZone][INST_0], config_container[iZone]);

      /*--- If AD mode is enabled we can use it to compute the projection,
       *    otherwise we use finite differences. ---*/

      if (config_container[iZone]->GetAD_Mode()) {
        if (config_container[iZone]->GetSmoothGradient()) {
          DerivativeTreatment_Gradient(geometry_container[iZone][INST_0], config_container[iZone], grid_movement[iZone], surface_movement[iZone] , Gradient);
        } else {
          SetProjection_AD(geometry_container[iZone][INST_0], config_container[iZone], surface_movement[iZone] , Gradient);
        }
      } else {
        SetProjection_FD(geometry_container[iZone][INST_0], config_container[iZone], surface_movement[iZone] , Gradient);
      }
    }
  } // for iZone

  /*--- Write the gradient to a file ---*/

  if (rank == MASTER_NODE)
    Gradient_file.open(config_container[ZONE_0]->GetObjFunc_Grad_FileName().c_str(), ios::out);

  /*--- Print gradients to screen and writes to file ---*/

  OutputGradient(Gradient, config_container[ZONE_0], Gradient_file);

  for (auto iDV = 0u; iDV  < config_container[ZONE_0]->GetnDV(); iDV++){
    delete [] Gradient[iDV];
  }
  delete [] Gradient;

  delete config;
  config = nullptr;

  if (rank == MASTER_NODE)
    cout << "\n------------------------- Solver Postprocessing -------------------------" << endl;

  if (geometry_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (geometry_container[iZone] != nullptr) {
        for (iInst = 0; iInst < nInst[iZone]; iInst++){
          if (geometry_container[iZone][iInst] != nullptr) {
            delete geometry_container[iZone][iInst];
          }
        }
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

  delete config;
  config = nullptr;
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
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << "\n------------------------- Exit Success (SU2_DOT) ------------------------\n" << endl;

  /*--- Finalize AD, if necessary. ---*/
#ifdef HAVE_OPDI
  AD::getGlobalTape().finalize();
#endif

  /*--- Finalize MPI parallelization. ---*/
  SU2_MPI::Finalize();

  /*--- Finalize OpenMP. ---*/
  omp_finalize();

  return EXIT_SUCCESS;

}

void SetProjection_FD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double** Gradient){

  unsigned short iDV, nDV, iFFDBox, nDV_Value, iMarker, iDim;
  unsigned long iVertex, iPoint;
  su2double delta_eps, my_Gradient, localGradient, *Normal, dS, *VarCoord, Sensitivity,
  dalpha[3], deps[3], dalpha_deps;
  bool *UpdatePoint, MoveSurface, Local_MoveSurface;
  CFreeFormDefBox **FFDBox;

  int rank = SU2_MPI::GetRank();

  nDV = config->GetnDV();

  /*--- Boolean controlling points to be updated ---*/

  UpdatePoint = new bool[geometry->GetnPoint()];

  /*--- Definition of the FFD deformation class ---*/

  unsigned short nFFDBox = MAX_NUMBER_FFD;
  FFDBox = new CFreeFormDefBox*[nFFDBox];
  for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) FFDBox[iFFDBox] = nullptr;

  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value = config->GetnDV_Value(iDV);
    if (nDV_Value != 1){
      SU2_MPI::Error("The projection using finite differences currently only supports a fixed direction of movement for FFD points.", CURRENT_FUNCTION);
    }
  }

  /*--- Continuous adjoint gradient computation ---*/

  if (rank == MASTER_NODE)
    cout << "Evaluate functional gradient using Finite Differences." << endl;

  for (iDV = 0; iDV < nDV; iDV++) {

    MoveSurface = true;
    Local_MoveSurface = true;

    /*--- Free Form deformation based ---*/

    if ((config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ||
        (config->GetDesign_Variable(iDV) == FFD_NACELLE) ||
        (config->GetDesign_Variable(iDV) == FFD_GULL) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST) ||
        (config->GetDesign_Variable(iDV) == FFD_ROTATION) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS) ||
        (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK)) {

      /*--- Read the FFD information in the first iteration ---*/

      if (iDV == 0) {

        if (rank == MASTER_NODE)
          cout << "Read the FFD information from mesh file." << endl;

        /*--- Read the FFD information from the grid file ---*/

        surface_movement->ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());

        /*--- If the FFDBox was not defined in the input file ---*/
        if (!surface_movement->GetFFDBoxDefinition()) {
          SU2_MPI::Error("The input grid doesn't have the entire FFD information!", CURRENT_FUNCTION);
        }

        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

          if (rank == MASTER_NODE) cout << "Checking FFD box dimension." << endl;
          surface_movement->CheckFFDDimension(geometry, config, FFDBox[iFFDBox], iFFDBox);

          if (rank == MASTER_NODE) cout << "Check the FFD box intersections with the solid surfaces." << endl;
          surface_movement->CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);

        }

        if (rank == MASTER_NODE)
          cout <<"-------------------------------------------------------------------------" << endl;

      }

      if (rank == MASTER_NODE) {
        cout << endl << "Design variable number "<< iDV <<"." << endl;
        cout << "Performing 3D deformation of the surface." << endl;
      }

      /*--- Apply the control point change ---*/

      MoveSurface = false;

      for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

        /*--- Reset FFD box ---*/

        switch (config->GetDesign_Variable(iDV) ) {
          case FFD_CONTROL_POINT_2D : Local_MoveSurface = surface_movement->SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER_2D :        Local_MoveSurface = surface_movement->SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS_2D :     Local_MoveSurface = surface_movement->SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST_2D :         Local_MoveSurface = surface_movement->SetFFDTwist_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_POINT :    Local_MoveSurface = surface_movement->SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_NACELLE :          Local_MoveSurface = surface_movement->SetFFDNacelle(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_GULL :             Local_MoveSurface = surface_movement->SetFFDGull(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST :            Local_MoveSurface = surface_movement->SetFFDTwist(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_ROTATION :         Local_MoveSurface = surface_movement->SetFFDRotation(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER :           Local_MoveSurface = surface_movement->SetFFDCamber(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS :        Local_MoveSurface = surface_movement->SetFFDThickness(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_SURFACE :  Local_MoveSurface = surface_movement->SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_ANGLE_OF_ATTACK :  Gradient[iDV][0] = config->GetAoA_Sens(); break;
        }

        /*--- Recompute cartesian coordinates using the new control points position ---*/

        if (Local_MoveSurface) {
          MoveSurface = true;
          surface_movement->SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, true);
        }

      }

    }

    /*--- Hicks Henne design variable ---*/

    else if (config->GetDesign_Variable(iDV) == HICKS_HENNE) {
      surface_movement->SetHicksHenne(geometry, config, iDV, true);
    }

    /*--- Surface bump design variable ---*/

    else if (config->GetDesign_Variable(iDV) == SURFACE_BUMP) {
      surface_movement->SetSurface_Bump(geometry, config, iDV, true);
    }

    /*--- Kulfan (CST) design variable ---*/

    else if (config->GetDesign_Variable(iDV) == CST) {
      surface_movement->SetCST(geometry, config, iDV, true);
    }

    /*--- Displacement design variable ---*/

    else if (config->GetDesign_Variable(iDV) == TRANSLATION) {
      surface_movement->SetTranslation(geometry, config, iDV, true);
    }

    /*--- Angle of Attack design variable ---*/

    else if (config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) {
      Gradient[iDV][0] = config->GetAoA_Sens();
    }

    /*--- Scale design variable ---*/

    else if (config->GetDesign_Variable(iDV) == SCALE) {
      surface_movement->SetScale(geometry, config, iDV, true);
    }

    /*--- Rotation design variable ---*/

    else if (config->GetDesign_Variable(iDV) == ROTATION) {
      surface_movement->SetRotation(geometry, config, iDV, true);
    }

    /*--- NACA_4Digits design variable ---*/

    else if (config->GetDesign_Variable(iDV) == NACA_4DIGITS) {
      surface_movement->SetNACA_4Digits(geometry, config);
    }

    /*--- Parabolic design variable ---*/

    else if (config->GetDesign_Variable(iDV) == PARABOLIC) {
      surface_movement->SetParabolic(geometry, config);
    }

    /*--- Design variable not implemented ---*/

    else {
      if (rank == MASTER_NODE)
        cout << "Design Variable not implemented yet" << endl;
    }

    /*--- Load the delta change in the design variable (finite difference step). ---*/

    if ((config->GetDesign_Variable(iDV) != ANGLE_OF_ATTACK) &&
        (config->GetDesign_Variable(iDV) != FFD_ANGLE_OF_ATTACK)) {

      /*--- If the Angle of attack is not involved, reset the value of the gradient ---*/

      my_Gradient = 0.0; Gradient[iDV][0] = 0.0;

      if (MoveSurface) {

        delta_eps = config->GetDV_Value(iDV);

        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
          UpdatePoint[iPoint] = true;

        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if ((iPoint < geometry->GetnPointDomain()) && UpdatePoint[iPoint]) {

                Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
                Sensitivity = geometry->vertex[iMarker][iVertex]->GetAuxVar();

                dS = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  dS += Normal[iDim]*Normal[iDim];
                  deps[iDim] = VarCoord[iDim] / delta_eps;
                }
                dS = sqrt(dS);

                dalpha_deps = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  dalpha[iDim] = Normal[iDim] / dS;
                  dalpha_deps -= dalpha[iDim]*deps[iDim];
                }

                my_Gradient += Sensitivity*dalpha_deps;
                UpdatePoint[iPoint] = false;
              }
            }
          }
        }
      }

      SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      Gradient[iDV][0] += localGradient;
    }
  }

  /*--- Delete memory for parameterization. ---*/

  if (FFDBox != nullptr) {
    for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) {
      if (FFDBox[iFFDBox] != nullptr) {
        delete FFDBox[iFFDBox];
      }
    }
    delete [] FFDBox;
  }

  delete [] UpdatePoint;

}


void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double** Gradient){

  su2double *VarCoord = nullptr, Sensitivity, my_Gradient, localGradient, *Normal, Area = 0.0;
  unsigned short iDV_Value = 0, iMarker, nMarker, iDim, nDim, iDV, nDV;
  unsigned long iVertex, nVertex, iPoint;

  const int rank = SU2_MPI::GetRank();

  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();

  /*--- Discrete adjoint gradient computation ---*/

  if (rank == MASTER_NODE)
    cout  << endl << "Evaluate functional gradient using Algorithmic Differentiation (ZONE " << config->GetiZone() << ")." << endl;

  /*--- Start recording of operations ---*/

  AD::StartRecording();

  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/

  for (iDV = 0; iDV < nDV; iDV++){

    for (iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++){

      config->SetDV_Value(iDV, iDV_Value, 0.0);

      AD::RegisterInput(config->GetDV_Value(iDV, iDV_Value));
    }
  }

  /*--- Call the surface deformation routine ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- Stop the recording --- */

  AD::StopRecording();

  /*--- Create a structure to identify points that have been already visited.
   * We need that to make sure to set the sensitivity of surface points only once
   *  (Markers share points, so we would visit them more than once in the loop over the markers below) ---*/

  vector<bool> visited(geometry->GetnPoint(), false);

  /*--- Initialize the derivatives of the output of the surface deformation routine
   * with the discrete adjoints from the CFD solution ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint      = geometry->vertex[iMarker][iVertex]->GetNode();
        if (!visited[iPoint]){
          VarCoord    = geometry->vertex[iMarker][iVertex]->GetVarCoord();
          Normal      = geometry->vertex[iMarker][iVertex]->GetNormal();

          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++){
            Area += Normal[iDim]*Normal[iDim];
          }
          Area = sqrt(Area);

          for (iDim = 0; iDim < nDim; iDim++){
            if (config->GetDiscrete_Adjoint()){
              Sensitivity = geometry->GetSensitivity(iPoint, iDim);
            } else {
              Sensitivity = -Normal[iDim]*geometry->vertex[iMarker][iVertex]->GetAuxVar()/Area;
            }
            SU2_TYPE::SetDerivative(VarCoord[iDim], SU2_TYPE::GetValue(Sensitivity));
          }
          visited[iPoint] = true;
        }
      }
    }
  }

  /*--- Compute derivatives and extract gradient ---*/

  AD::ComputeAdjoint();

  for (iDV = 0; iDV  < nDV; iDV++){

    for (iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++){

      my_Gradient = SU2_TYPE::GetDerivative(config->GetDV_Value(iDV, iDV_Value));

      SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

      /*--- Angle of Attack design variable (this is different,
       the value comes form the input file) ---*/

      if ((config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) ||
          (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK))  {
        Gradient[iDV][iDV_Value] = config->GetAoA_Sens();
      }

      Gradient[iDV][iDV_Value] += localGradient;
    }
  }

  AD::Reset();

}

void OutputGradient(su2double** Gradient, CConfig* config, ofstream& Gradient_file){

  unsigned short nDV, iDV, iDV_Value, nDV_Value;

  int rank = SU2_MPI::GetRank();

  nDV = config->GetnDV();

  /*--- Loop through all design variables and their gradients ---*/

  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value = config->GetnDV_Value(iDV);
    if (rank == MASTER_NODE){

      /*--- Print the kind of design variable on screen ---*/

      cout << endl << "Design variable (";
      for (std::map<string, ENUM_PARAM>::const_iterator it = Param_Map.begin(); it != Param_Map.end(); ++it ){
        if (it->second == config->GetDesign_Variable(iDV)){
          cout << it->first << ") number "<< iDV << "." << endl;
        }
      }

      /*--- Print the kind of objective function to screen ---*/

      for (std::map<string, ENUM_OBJECTIVE>::const_iterator it = Objective_Map.begin(); it != Objective_Map.end(); ++it ){
        if (it->second == config->GetKind_ObjFunc()){
          cout << it->first << " gradient : ";
          if (iDV == 0) Gradient_file << it->first << " gradient " << endl;
        }
      }

      /*--- Print the gradient to file and screen ---*/

      for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
        cout << Gradient[iDV][iDV_Value];
        if (iDV_Value != nDV_Value-1 ){
          cout << ", ";
        }
        Gradient_file << Gradient[iDV][iDV_Value] << endl;
      }
      cout << endl;
      cout <<"-------------------------------------------------------------------------" << endl;
    }
  }
}


void SetSensitivity_Files(CGeometry ***geometry, CConfig **config, unsigned short val_nZone) {

  unsigned short iMarker,iDim, nDim, nMarker, nVar;
  unsigned long iVertex, iPoint, nPoint, nVertex;
  su2double *Normal, Prod, Sens = 0.0, SensDim, Area;

  unsigned short iZone;

  CSolver *solver  = nullptr;
  COutput *output  = nullptr;

  for (iZone = 0; iZone < val_nZone; iZone++) {

    nPoint = geometry[iZone][INST_0]->GetnPoint();
    nDim   = geometry[iZone][INST_0]->GetnDim();
    nMarker = config[iZone]->GetnMarker_All();
    nVar = nDim + 1;

    /*--- We create a baseline solver to easily merge the sensitivity information ---*/

    vector<string> fieldnames;
    fieldnames.push_back("\"Point\"");
    fieldnames.push_back("\"x\"");
    fieldnames.push_back("\"y\"");
    if (nDim == 3) {
      fieldnames.push_back("\"z\"");
    }
    fieldnames.push_back("\"Sensitivity_x\"");
    fieldnames.push_back("\"Sensitivity_y\"");
    if (nDim == 3) {
      fieldnames.push_back("\"Sensitivity_z\"");
    }
    fieldnames.push_back("\"Surface_Sensitivity\"");

    solver = new CBaselineSolver(geometry[iZone][INST_0], config[iZone], nVar+nDim, fieldnames);

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        solver->GetNodes()->SetSolution(iPoint, iDim,      geometry[iZone][INST_0]->nodes->GetCoord(iPoint, iDim));
        solver->GetNodes()->SetSolution(iPoint, iDim+nDim, geometry[iZone][INST_0]->GetSensitivity(iPoint, iDim));
      }
    }

    /*--- Compute the sensitivity in normal direction ---*/

    for (iMarker = 0; iMarker < nMarker; iMarker++) {

      if(config[iZone]->GetSolid_Wall(iMarker) || (config[iZone]->GetMarker_All_DV(iMarker) == YES )) {

        nVertex = geometry[iZone][INST_0]->GetnVertex(iMarker);

        for (iVertex = 0; iVertex < nVertex; iVertex++) {
          iPoint = geometry[iZone][INST_0]->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry[iZone][INST_0]->vertex[iMarker][iVertex]->GetNormal();
          Prod = 0.0;
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {

            /*--- Retrieve the gradient calculated with discrete adjoint method ---*/

            SensDim = geometry[iZone][INST_0]->GetSensitivity(iPoint, iDim);

            /*--- Calculate scalar product for projection onto the normal vector ---*/

            Prod += Normal[iDim]*SensDim;

            Area += Normal[iDim]*Normal[iDim];
          }

          Area = sqrt(Area);

          /*--- Projection of the gradient onto the normal vector of the surface ---*/

          Sens = Prod/Area;

          solver->GetNodes()->SetSolution(iPoint, 2*nDim, Sens);

        }
      }
    }

    output = new CBaselineOutput(config[iZone], nDim, solver);
    output->PreprocessVolumeOutput(config[iZone]);
    output->PreprocessHistoryOutput(config[iZone], false);

    /*--- Load the data --- */

    output->Load_Data(geometry[iZone][INST_0], config[iZone], &solver);

    /*--- Set the surface filename ---*/

    output->SetSurface_Filename(config[iZone]->GetSurfSens_FileName());

    /*--- Set the volume filename ---*/

    output->SetVolume_Filename(config[iZone]->GetVolSens_FileName());

    /*--- Write to file ---*/

    for (unsigned short iFile = 0; iFile < config[iZone]->GetnVolumeOutputFiles(); iFile++){
      auto FileFormat = config[iZone]->GetVolumeOutputFiles();
      if (FileFormat[iFile] != OUTPUT_TYPE::RESTART_ASCII &&
          FileFormat[iFile] != OUTPUT_TYPE::RESTART_BINARY &&
          FileFormat[iFile] != OUTPUT_TYPE::CSV)
        output->WriteToFile(config[iZone], geometry[iZone][INST_0], FileFormat[iFile]);
    }

    /*--- Free memory ---*/

    delete output;
    delete solver;

  }

}

void DerivativeTreatment_MeshSensitivity(CGeometry *geometry, CConfig *config, CVolumetricMovement *grid_movement) {

  int rank = SU2_MPI::GetRank();

  /*--- Warning if choosen smoothing mode is unsupported.
   * This is the default option if the user has not specified a mode in the config file. ---*/
  if (config->GetSobMode() == ENUM_SOBOLEV_MODUS::NONE) {
    SU2_MPI::Error("Unsupported operation modus for the Sobolev Smoothing Solver.", CURRENT_FUNCTION);
  }

  /*-- Construct the smoothing solver and numerics ---*/
  std::unique_ptr<CSolver> solver(new CGradientSmoothingSolver(geometry, config));
  unsigned dim = (config->GetSmoothOnSurface() ? geometry->GetnDim() - 1 : geometry->GetnDim());
  std::unique_ptr<CNumerics> numerics(new CGradSmoothing(dim, config));

  if (rank == MASTER_NODE)  cout << "Sobolev Smoothing of derivatives is active." << endl;

  /*--- Apply the smoothing procedure on the mesh level. ---*/
  if (config->GetSobMode() == ENUM_SOBOLEV_MODUS::MESH_LEVEL) {
    if (rank == MASTER_NODE) cout << "  working on mesh level" << endl;

    /*--- Work with the surface derivatives. ---*/
    if (config->GetSmoothOnSurface()) {

      /*--- Project to surface and then smooth on surface. ---*/
      grid_movement->SetVolume_Deformation(geometry, config, false, true);

      /*--- Get the sensitivities from the geometry class to work with. ---*/
      solver->ReadSensFromGeometry(geometry);

      /*--- Perform the smoothing procedure on all boundaries marked as DV marker. ---*/
      solver->ApplyGradientSmoothingSurface(geometry, numerics.get(), config);

      /*--- After appling the solver write the results back ---*/
      solver->WriteSensToGeometry(geometry);

    /*--- Work with the volume derivatives. ---*/
    } else {

      /*--- Get the sensitivities from the geometry class to work with. ---*/
      solver->ReadSensFromGeometry(geometry);

      solver->ApplyGradientSmoothingVolume(geometry, numerics.get(), config);

      /*--- After appling the solver write the results back ---*/
      solver->WriteSensToGeometry(geometry);

      /*--- Projection is applied after smoothing. ---*/
      grid_movement->SetVolume_Deformation(geometry, config, false, true);
    }

  }

}

void DerivativeTreatment_Gradient(CGeometry *geometry, CConfig *config, CVolumetricMovement* grid_movement, CSurfaceMovement *surface_movement, su2double** Gradient) {

  int rank = SU2_MPI::GetRank();

  /*--- Error if choosen smoothing mode is unsupported.
   * This is the default option if the user has not specified a mode in the config file. ---*/
  if (config->GetSobMode() == ENUM_SOBOLEV_MODUS::NONE) {
    SU2_MPI::Error("Unsupported operation modus for the Sobolev Smoothing Solver.", CURRENT_FUNCTION);
  }

  /*-- Construct the smoothing solver and numerics ---*/
  std::unique_ptr<CSolver> solver(new CGradientSmoothingSolver(geometry, config));
  unsigned dim = (config->GetSmoothOnSurface() ? geometry->GetnDim() - 1 : geometry->GetnDim());
  std::unique_ptr<CNumerics> numerics(new CGradSmoothing(dim, config));

  if (rank == MASTER_NODE)  cout << "Sobolev Smoothing of derivatives is active." << endl;

  /*--- Get the sensitivities from the geometry class to work with. ---*/
  solver->ReadSensFromGeometry(geometry);

  /*--- Apply the smoothing procedure on the DV level. ---*/
  if (config->GetSobMode() == ENUM_SOBOLEV_MODUS::PARAM_LEVEL_COMPLETE) {
    solver->ApplyGradientSmoothingDV(geometry, numerics.get(), surface_movement, grid_movement, config, Gradient);

    /*--- If smoothing already took place on the mesh level, or none is requested, just do standard projection. ---*/
  } else if (config->GetSobMode() == ENUM_SOBOLEV_MODUS::ONLY_GRAD ||
             config->GetSobMode() == ENUM_SOBOLEV_MODUS::MESH_LEVEL) {
    solver->RecordTapeAndCalculateOriginalGradient(geometry, surface_movement, grid_movement, config, Gradient);
  }

}
