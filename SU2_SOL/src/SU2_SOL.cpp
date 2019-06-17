/*!
 * \file SU2_SOL.cpp
 * \brief Main file for the solution export/conversion code (SU2_SOL).
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

#include "../include/SU2_SOL.hpp"

using namespace std;

int main(int argc, char *argv[]) {

  unsigned short iZone, nZone = SINGLE_ZONE, iInst;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  ofstream ConvHist_file;
  char config_file_name[MAX_STRING_SIZE];
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  bool fem_solver = false;
  bool multizone = false;

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

  COutput *output                 = NULL;
  CGeometry ***geometry_container = NULL;
  CSolver ***solver_container     = NULL;
  CConfig **config_container      = NULL;
  CConfig *driver_config          = NULL;
  unsigned short *nInst           = NULL;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/

  if (argc == 2 || argc == 3) { strcpy(config_file_name,argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }

  CConfig *config = NULL;
  config = new CConfig(config_file_name, SU2_SOL);

  if (config->GetKind_Solver() == MULTIZONE) nZone  = config->GetnConfigFiles();
  else nZone  = CConfig::GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);

  /*--- Definition of the containers per zones ---*/

  solver_container = new CSolver**[nZone];
  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry**[nZone];
  nInst = new unsigned short[nZone];
  driver_config = NULL;

  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone]       = NULL;
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
    nInst[iZone]                  = 1;
  }

  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_SOL, ZONE_0, nZone, 0, false);

  /*--- Initialize a char to store the zone filename ---*/
  char zone_file_name[MAX_STRING_SIZE];

  /*--- Store a boolean for multizone problems ---*/
  multizone = (driver_config->GetKind_Solver() == MULTIZONE);

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    if (multizone){
      strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
      config_container[iZone] = new CConfig(zone_file_name, SU2_SOL, iZone, nZone, 0, true);
    }
    else{
      config_container[iZone] = new CConfig(config_file_name, SU2_SOL, iZone, nZone, 0, true);
    }

    config_container[iZone]->SetMPICommunicator(MPICommunicator);

  }

  /*--- Set the multizone part of the problem. ---*/
  if (driver_config->GetKind_Solver() == MULTIZONE){
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config_container[iZone]->SetMultizone(driver_config, config_container);
    }
  }

  /*--- Read the geometry for each zone ---*/
  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Determine whether or not the FEM solver is used, which decides the
     type of geometry classes that are instantiated. ---*/
    fem_solver = ((config_container[iZone]->GetKind_Solver() == FEM_EULER)          ||
                  (config_container[iZone]->GetKind_Solver() == FEM_NAVIER_STOKES)  ||
                  (config_container[iZone]->GetKind_Solver() == FEM_RANS)           ||
                  (config_container[iZone]->GetKind_Solver() == FEM_LES)            ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_EULER) ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_NS)    ||
                  (config_container[iZone]->GetKind_Solver() == DISC_ADJ_FEM_RANS));

    /*--- Read the number of instances for each zone ---*/

    nInst[iZone] = config_container[iZone]->GetnTimeInstances();

    geometry_container[iZone] = new CGeometry*[nInst[iZone]];
    solver_container[iZone] = new CSolver*[nInst[iZone]];

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      /*--- Allocate solver. ---*/
      solver_container[iZone][iInst] = NULL;

      config_container[iZone]->SetiInst(iInst);

      /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/

      CGeometry *geometry_aux = NULL;

      /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

      geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

      /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

      if ( fem_solver ) geometry_aux->SetColorFEMGrid_Parallel(config_container[iZone]);
      else              geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

      /*--- Allocate the memory of the current domain, and
     divide the grid between the nodes ---*/

      geometry_container[iZone][iInst] = NULL;

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

      if (rank == MASTER_NODE) cout << "Identify vertices." <<endl;
      geometry_container[iZone][iInst]->SetVertex(config_container[iZone]);

      /*--- Store the global to local mapping after preprocessing. ---*/

      if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
      geometry_container[iZone][iInst]->SetGlobal_to_Local_Point();

      /*--- Create the point-to-point MPI communication structures for the fvm solver. ---*/
      
      if (!fem_solver) geometry_container[iZone][iInst]->PreprocessP2PComms(geometry_container[iZone][iInst], config_container[iZone]);
      
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

  /*--- Determine whether the simulation is a FSI simulation ---*/

  bool fsi = config_container[ZONE_0]->GetFSI_Simulation();

  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/

#ifdef HAVE_MPI
  StartTime = MPI_Wtime();
#else
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solution Postprocessing -----------------------" << endl;
  
	/*--- Definition of the output class (one for all the zones) ---*/
	output = new COutput(config_container[ZONE_0]);
  
  /*---  Check whether this is an FSI, fluid unsteady, harmonic balance or structural dynamic simulation and call the
   solution merging routines accordingly.---*/

  if (multizone){

    bool TimeDomain = driver_config->GetTime_Domain();

    if (TimeDomain){

      su2double Physical_dt, Physical_t;
      unsigned long TimeIter = 0;
      bool StopCalc = false;

      /*--- Physical time step ---*/
      Physical_dt = driver_config->GetTime_Step();

      /*--- Check for an unsteady restart. Update TimeIter if necessary. ---*/
      if (driver_config->GetRestart()){
        TimeIter = driver_config->GetRestart_Iter();
      }

      /*--- Instantiate the solvers for each zone. ---*/
      for (iZone = 0; iZone < nZone; iZone++){
        config_container[iZone]->SetiInst(INST_0);
        config_container[iZone]->SetExtIter(TimeIter);
        solver_container[iZone][INST_0] = new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
      }

      /*--- Loop over the whole time domain ---*/
      while (TimeIter < driver_config->GetnTime_Iter()) {

        /*--- Check if the maximum time has been surpassed. ---*/
        Physical_t  = (TimeIter+1)*Physical_dt;
        if (Physical_t >=  driver_config->GetMax_Time())
          StopCalc = true;

        if ((TimeIter+1 == driver_config->GetnTime_Iter()) || // The last time iteration
            (StopCalc) || // We have surpassed the requested time
            ((TimeIter == 0) || (TimeIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)) // The iteration has been requested
          ){

          /*--- Load the restart for all the zones. ---*/
          for (iZone = 0; iZone < nZone; iZone++){

            /*--- Set the current iteration number in the config class. ---*/
            config_container[iZone]->SetExtIter(TimeIter);
            /*--- So far, only enabled for single-instance problems ---*/
            config_container[iZone]->SetiInst(INST_0);
            solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone], config_container[iZone], SU2_TYPE::Int(MESH_0), true);
          }

          if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << TimeIter << ", t = " << Physical_t << " s ." << endl;
          output->SetBaselineResult_Files(solver_container, geometry_container, config_container, TimeIter, nZone);
        }

        TimeIter++;
        if (StopCalc) break;
      }
    }
    else {

      /*--- Steady simulation: merge the solution files for each zone. ---*/
      for (iZone = 0; iZone < nZone; iZone++) {
        config_container[iZone]->SetiInst(INST_0);
        /*--- Definition of the solution class ---*/
        solver_container[iZone][INST_0] = new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
        solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone], config_container[iZone], SU2_TYPE::Int(MESH_0), true);
      }
      output->SetBaselineResult_Files(solver_container, geometry_container, config_container, 0, nZone);
    }

  }
  else if (fsi){

    if (nZone < 2){
      SU2_MPI::Error("For multizone computations, please add the number of zones as a second argument for SU2_SOL.", CURRENT_FUNCTION);
    }

    su2double Physical_dt, Physical_t;
    unsigned long iExtIter = 0, iExtIterFlow = 0, iExtIterFEM = 0;
    bool StopCalc = false;
    bool SolutionInstantiatedFlow = false, SolutionInstantiatedFEM = false;

    /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
    if (config_container[ZONE_0]->GetRestart()){
      iExtIterFlow = config_container[ZONE_0]->GetUnst_RestartIter();
      iExtIterFEM = config_container[ZONE_1]->GetDyn_RestartIter();
      if (iExtIterFlow != iExtIterFEM) {
        SU2_MPI::Error("For multizone computations, please add the number of zones as a second argument for SU2_SOL.", CURRENT_FUNCTION);
      }
      else {
        iExtIter = iExtIterFlow;
      }
    }


    while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

      /*--- Check several conditions in order to merge the correct time step files. ---*/

      Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
      Physical_t  = (iExtIter+1)*Physical_dt;
      if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
        StopCalc = true;

      if (
          ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
           ((iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (iExtIter != 0) &&
            !((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
              (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
           (StopCalc) ||
           (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
             (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
            ((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0))))

          &&

          ((iExtIter+1 == config_container[ZONE_1]->GetnExtIter()) ||
           (StopCalc) ||
           ((config_container[ZONE_1]->GetDynamic_Analysis() == DYNAMIC) &&
            ((iExtIter == 0) || (iExtIter % config_container[ZONE_1]->GetWrt_Sol_Freq_DualTime() == 0))))

          ){

        /*--- Set the current iteration number in the config class. ---*/
        config_container[ZONE_0]->SetExtIter(iExtIter);
        config_container[ZONE_1]->SetExtIter(iExtIter);

        /*--- Read in the restart file for this time step ---*/

        /*--- For the fluid zone (ZONE_0) ---*/
        /*--- Either instantiate the solution class or load a restart file. ---*/
        if (SolutionInstantiatedFlow == false &&
            (iExtIter == 0 || ((config_container[ZONE_0]->GetRestart() && (SU2_TYPE::Int(iExtIter) == config_container[ZONE_0]->GetUnst_RestartIter())) ||
                               iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
                               iExtIter+1 == config_container[ZONE_0]->GetnExtIter()))) {
          solver_container[ZONE_0][INST_0] = new CBaselineSolver(geometry_container[ZONE_0][INST_0], config_container[ZONE_0]);
          SolutionInstantiatedFlow = true;
        }
          solver_container[ZONE_0][INST_0]->LoadRestart_FSI(geometry_container[ZONE_0][INST_0], config_container[ZONE_0], SU2_TYPE::Int(MESH_0));


        /*--- For the structural zone (ZONE_1) ---*/
        /*--- Either instantiate the solution class or load a restart file. ---*/
        /*--- Either instantiate the solution class or load a restart file. ---*/
        if (SolutionInstantiatedFEM == false &&
            (iExtIter == 0 || ((config_container[ZONE_1]->GetRestart() && (SU2_TYPE::Int(iExtIter) == config_container[ZONE_1]->GetDyn_RestartIter())) ||
                               iExtIter % config_container[ZONE_1]->GetWrt_Sol_Freq_DualTime() == 0 ||
                               iExtIter+1 == config_container[ZONE_1]->GetnExtIter()))) {
          solver_container[ZONE_1][INST_0] = new CBaselineSolver(geometry_container[ZONE_1][INST_0], config_container[ZONE_1]);
          SolutionInstantiatedFEM = true;
        }
          solver_container[ZONE_1][INST_0]->LoadRestart_FSI(geometry_container[ZONE_1][INST_0], config_container[ZONE_1], SU2_TYPE::Int(MESH_0));

        if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
        output->SetBaselineResult_Files(solver_container, geometry_container, config_container, iExtIter, nZone);
      }

      iExtIter++;
      if (StopCalc) break;
    }

  } else if (fem_solver) {

    if (config_container[ZONE_0]->GetWrt_Unsteady()) {

      /*--- Unsteady DG simulation: merge all unsteady time steps. First,
       find the frequency and total number of files to write. ---*/

      su2double Physical_dt, Physical_t;
      unsigned long iExtIter = 0;
      bool StopCalc = false;
      bool *SolutionInstantiated = new bool[nZone];

      for (iZone = 0; iZone < nZone; iZone++)
        SolutionInstantiated[iZone] = false;

      /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
      if (config_container[ZONE_0]->GetWrt_Unsteady() && config_container[ZONE_0]->GetRestart())
        iExtIter = config_container[ZONE_0]->GetUnst_RestartIter();

      while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

        /*--- Check several conditions in order to merge the correct time step files. ---*/
        Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
        Physical_t  = (iExtIter+1)*Physical_dt;
        if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
          StopCalc = true;

        if ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
            ((iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (iExtIter != 0) &&
             !(config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_STEPPING)) ||
            (StopCalc) ||
            ((config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_STEPPING) &&
             ((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {

              /*--- Read in the restart file for this time step ---*/
              for (iZone = 0; iZone < nZone; iZone++) {

                /*--- Set the current iteration number in the config class. ---*/
                config_container[iZone]->SetExtIter(iExtIter);

                /*--- Either instantiate the solution class or load a restart file. ---*/
                if (SolutionInstantiated[iZone] == false &&
                    (iExtIter == 0 ||
                     (config_container[ZONE_0]->GetRestart() && ((long)iExtIter == config_container[ZONE_0]->GetUnst_RestartIter() ||
                                                                                  iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
                                                                                  iExtIter+1 == config_container[ZONE_0]->GetnExtIter())))) {

                  solver_container[iZone][INST_0] = new CBaselineSolver_FEM(geometry_container[iZone][INST_0], config_container[iZone]);
                  SolutionInstantiated[iZone] = true;
                }
                solver_container[iZone][INST_0]->LoadRestart(&geometry_container[iZone][INST_0], &solver_container[iZone],
                                                             config_container[iZone], (int)iExtIter, true);
              }

              if (rank == MASTER_NODE)
                cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
              output->SetBaselineResult_Files_FEM(solver_container, geometry_container, config_container, iExtIter, nZone);
            }
        
        iExtIter++;
        if (StopCalc) break;
      }
      
    } else {

    /*--- Steady simulation: merge the single solution file. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Definition of the solution class ---*/

      solver_container[iZone][INST_0] = new CBaselineSolver_FEM(geometry_container[iZone][INST_0], config_container[iZone]);
      solver_container[iZone][INST_0]->LoadRestart(&geometry_container[iZone][INST_0], &solver_container[iZone], config_container[iZone], SU2_TYPE::Int(MESH_0), true);
    }

    output->SetBaselineResult_Files_FEM(solver_container, geometry_container, config_container, 0, nZone);
    }

  }
  else {

    if (config_container[ZONE_0]->GetWrt_Unsteady()) {

      /*--- Unsteady simulation: merge all unsteady time steps. First,
       find the frequency and total number of files to write. ---*/

      su2double Physical_dt, Physical_t;
      unsigned long iExtIter = 0;
      bool StopCalc = false;
      bool *SolutionInstantiated = new bool[nZone];

      for (iZone = 0; iZone < nZone; iZone++)
        SolutionInstantiated[iZone] = false;

      /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
      if (config_container[ZONE_0]->GetWrt_Unsteady() && config_container[ZONE_0]->GetRestart())
        iExtIter = config_container[ZONE_0]->GetUnst_RestartIter();

      while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

        /*--- Check several conditions in order to merge the correct time step files. ---*/
        Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
        Physical_t  = (iExtIter+1)*Physical_dt;
        if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
          StopCalc = true;

        if ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
            ((iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (iExtIter != 0) &&
             !((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
               (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
            (StopCalc) ||
            (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
              (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
             ((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {



              /*--- Read in the restart file for this time step ---*/
              for (iZone = 0; iZone < nZone; iZone++) {

                /*--- Set the current iteration number in the config class. ---*/
                config_container[iZone]->SetExtIter(iExtIter);

                /*--- Either instantiate the solution class or load a restart file. ---*/
                if (SolutionInstantiated[iZone] == false &&
                    (iExtIter == 0 || (config_container[ZONE_0]->GetRestart() && ((long)iExtIter == config_container[ZONE_0]->GetUnst_RestartIter() ||
                                                                                  iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
                                                                                  iExtIter+1 == config_container[ZONE_0]->GetnExtIter())))) {
                  solver_container[iZone][INST_0] = new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
                  SolutionInstantiated[iZone] = true;
                }
                  config_container[iZone]->SetiInst(INST_0);
                  solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone], config_container[iZone], SU2_TYPE::Int(MESH_0), true);
              }

              if (rank == MASTER_NODE)
                cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
              output->SetBaselineResult_Files(solver_container, geometry_container, config_container, iExtIter, nZone);
            }

        iExtIter++;
        if (StopCalc) break;
      }

    }

    else if (config_container[ZONE_0]->GetUnsteady_Simulation() == HARMONIC_BALANCE) {

      /*--- Read in the restart file for this time step ---*/
      for (iZone = 0; iZone < nZone; iZone++) {

        for (iInst = 0; iInst < nInst[iZone]; iInst++){

          config_container[iZone]->SetiInst(iInst);

          /*--- Either instantiate the solution class or load a restart file. ---*/
          solver_container[iZone][iInst] = new CBaselineSolver(geometry_container[iZone][iInst], config_container[iZone]);
          solver_container[iZone][iInst]->LoadRestart(geometry_container[iZone], &solver_container[iZone], config_container[iZone], SU2_TYPE::Int(MESH_0), true);

          /*--- Print progress in solution writing to the screen. ---*/
          if (rank == MASTER_NODE) {
            cout << "Storing the volume solution for time instance " << iInst << "." << endl;
          }

        }

      }

      output->SetBaselineResult_Files(solver_container, geometry_container, config_container, iZone, nZone);
    }

    else if (config_container[ZONE_0]->GetWrt_Dynamic()){

      /*--- Dynamic simulation: merge all unsteady time steps. First,
       find the frequency and total number of files to write. ---*/

      su2double Physical_dt, Physical_t;
      unsigned long iExtIter = 0;
      bool StopCalc = false;
      bool SolutionInstantiated = false;



      /*--- Check for an dynamic restart (structural analysis). Update ExtIter if necessary. ---*/
      if (config_container[ZONE_0]->GetKind_Solver() == FEM_ELASTICITY &&
          config_container[ZONE_0]->GetWrt_Dynamic() && config_container[ZONE_0]->GetRestart())
        iExtIter = config_container[ZONE_0]->GetDyn_RestartIter();

      while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

        /*--- Check several conditions in order to merge the correct time step files. ---*/
        /*--- If the solver is structural, the total and delta_t are obtained from different functions. ---*/

        Physical_dt = config_container[ZONE_0]->GetDelta_DynTime();
        Physical_t  = (iExtIter+1)*Physical_dt;
        if (Physical_t >=  config_container[ZONE_0]->GetTotal_DynTime())
          StopCalc = true;

        if ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
            (StopCalc) ||
            ((config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC) &&
             ((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {

              /*--- Set the current iteration number in the config class. ---*/
              config_container[ZONE_0]->SetExtIter(iExtIter);

              /*--- Read in the restart file for this time step ---*/
              for (iZone = 0; iZone < nZone; iZone++) {

                /*--- Either instantiate the solution class or load a restart file. ---*/
                if (SolutionInstantiated == false &&
                    (iExtIter == 0 || ((config_container[ZONE_0]->GetRestart() && (SU2_TYPE::Int(iExtIter) == config_container[ZONE_0]->GetDyn_RestartIter())) ||
                                       iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
                                       iExtIter+1 == config_container[ZONE_0]->GetnExtIter()))) {
                  solver_container[iZone][INST_0] = new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
                  SolutionInstantiated = true;
                }
                config_container[iZone]->SetiInst(INST_0);
                solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone], config_container[iZone], SU2_TYPE::Int(MESH_0), true);
              }

              if (rank == MASTER_NODE)
                cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
              output->SetBaselineResult_Files(solver_container, geometry_container, config_container, iExtIter, nZone);
            }
        
        iExtIter++;
        if (StopCalc) break;
      }
      
		  }
    
    else {

      /*--- Steady simulation: merge the single solution file. ---*/

      for (iZone = 0; iZone < nZone; iZone++) {
        config_container[iZone]->SetiInst(INST_0);
        /*--- Definition of the solution class ---*/
        solver_container[iZone][INST_0] = new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
        solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone], config_container[iZone], SU2_TYPE::Int(MESH_0), true);
      }

      output->SetBaselineResult_Files(solver_container, geometry_container, config_container, 0, nZone);

		  }
    
  }
  
  delete config;
  config = NULL;

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;
  
  if (geometry_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        if (geometry_container[iZone][iInst] != NULL) {
          delete geometry_container[iZone][iInst];
        }
      }
      if (geometry_container[iZone] != NULL)
        delete geometry_container[iZone];
    }
    delete [] geometry_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;
  
  if (solver_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        if (solver_container[iZone][iInst] != NULL) {
          delete solver_container[iZone][iInst];
        }
      }
      if (solver_container[iZone] != NULL)
        delete solver_container[iZone];
    }
    delete [] solver_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CSolver class." << endl;
  
  if (config_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != NULL) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;
  
  if (output != NULL) delete output;
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
    cout << endl <<"------------------------- Exit Success (SU2_SOL) ------------------------" << endl << endl;
  
  /*--- Finalize MPI parallelization ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Finalize();
#endif
  
  return EXIT_SUCCESS;
}
