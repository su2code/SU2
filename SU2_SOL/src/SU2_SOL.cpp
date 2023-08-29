/*!
 * \file SU2_SOL.cpp
 * \brief Main file for the solution export/conversion code (SU2_SOL).
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

int main(int argc, char* argv[]) {
  unsigned short iZone, iInst;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;

  char config_file_name[MAX_STRING_SIZE];

  /*--- MPI initialization ---*/

  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Comm MPICommunicator = SU2_MPI::GetComm();

  const int rank = SU2_MPI::GetRank();
  const int size = SU2_MPI::GetSize();

  /*--- Pointer to different structures that will be used throughout the entire code ---*/

  COutput** output = nullptr;
  CGeometry*** geometry_container = nullptr;
  CSolver*** solver_container = nullptr;
  CConfig** config_container = nullptr;
  CConfig* driver_config = nullptr;
  unsigned short* nInst = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/

  if (argc == 2 || argc == 3) {
    strcpy(config_file_name, argv[1]);
  } else {
    strcpy(config_file_name, "default.cfg");
  }

  CConfig* config = nullptr;
  config = new CConfig(config_file_name, SU2_COMPONENT::SU2_SOL);

  const auto nZone = config->GetnZone();

  /*--- Definition of the containers per zones ---*/

  solver_container = new CSolver**[nZone]();
  config_container = new CConfig*[nZone]();
  geometry_container = new CGeometry**[nZone]();
  nInst = new unsigned short[nZone];
  driver_config = nullptr;
  output = new COutput*[nZone]();

  for (iZone = 0; iZone < nZone; iZone++) {
    nInst[iZone] = 1;
  }

  /*--- Initialize the configuration of the driver ---*/
  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_SOL, false);

  /*--- Initialize a char to store the zone filename ---*/
  char zone_file_name[MAX_STRING_SIZE];

  /*--- Store a boolean for multizone problems ---*/
  const bool multizone = config->GetMultizone_Problem();

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    if (driver_config->GetnConfigFiles() > 0) {
      strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
      config_container[iZone] = new CConfig(driver_config, zone_file_name, SU2_COMPONENT::SU2_SOL, iZone, nZone, true);
    } else {
      config_container[iZone] =
          new CConfig(driver_config, config_file_name, SU2_COMPONENT::SU2_SOL, iZone, nZone, true);
    }

    config_container[iZone]->SetMPICommunicator(MPICommunicator);
  }

  /*--- Set the multizone part of the problem. ---*/
  if (multizone) {
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config_container[iZone]->SetMultizone(driver_config, config_container);
    }
  }

  /*--- Read the geometry for each zone ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    /*--- Determine whether or not the FEM solver is used, which decides the
     type of geometry classes that are instantiated. ---*/
    const bool fem_solver = config_container[iZone]->GetFEMSolver();

    /*--- Read the number of instances for each zone ---*/

    nInst[iZone] = config_container[iZone]->GetnTimeInstances();

    geometry_container[iZone] = new CGeometry*[nInst[iZone]];
    solver_container[iZone] = new CSolver*[nInst[iZone]];

    for (iInst = 0; iInst < nInst[iZone]; iInst++) {
      /*--- Allocate solver. ---*/
      solver_container[iZone][iInst] = nullptr;

      config_container[iZone]->SetiInst(iInst);

      /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/

      CGeometry* geometry_aux = nullptr;

      /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

      geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

      /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

      if (fem_solver)
        geometry_aux->SetColorFEMGrid_Parallel(config_container[iZone]);
      else
        geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

      /*--- Allocate the memory of the current domain, and divide the grid between the nodes ---*/

      geometry_container[iZone][iInst] = nullptr;

      /*--- Build the grid data structures using the ParMETIS coloring. ---*/

      if (fem_solver) {
        switch (config_container[iZone]->GetKind_FEM_Flow()) {
          case DG: {
            geometry_container[iZone][iInst] = new CMeshFEM_DG(geometry_aux, config_container[iZone]);
            break;
          }
        }
      } else {
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

      /*--- Create the point-to-point MPI communication structures for the fvm solver. ---*/

      if (!fem_solver)
        geometry_container[iZone][iInst]->PreprocessP2PComms(geometry_container[iZone][iInst], config_container[iZone]);

      /* Test for a fem solver, because some more work must be done. */

      if (fem_solver) {
        /*--- Carry out a dynamic cast to CMeshFEM_DG, such that it is not needed to
         define all virtual functions in the base class CGeometry. ---*/
        auto* DGMesh = dynamic_cast<CMeshFEM_DG*>(geometry_container[iZone][iInst]);

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

  const bool fsi = config_container[ZONE_0]->GetFSI_Simulation();
  const bool fem_solver = config_container[ZONE_0]->GetFEMSolver();

  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/

  StartTime = SU2_MPI::Wtime();

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Solution Postprocessing -----------------------" << endl;

  /*---  Check whether this is an FSI, fluid unsteady, harmonic balance or structural dynamic simulation and call the
   solution merging routines accordingly.---*/

  if (multizone) {
    bool TimeDomain = driver_config->GetTime_Domain();

    if (TimeDomain) {
      su2double Physical_dt, Physical_t;
      unsigned long TimeIter = 0;
      bool StopCalc = false;

      /*--- Physical time step ---*/
      Physical_dt = driver_config->GetTime_Step();

      /*--- Check for an unsteady restart. Update TimeIter if necessary. ---*/
      if (driver_config->GetRestart()) {
        TimeIter = driver_config->GetRestart_Iter();
      }

      /*--- Instantiate the solvers for each zone. ---*/
      for (iZone = 0; iZone < nZone; iZone++) {
        config_container[iZone]->SetiInst(INST_0);
        config_container[iZone]->SetTimeIter(TimeIter);
        solver_container[iZone][INST_0] =
            new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);

        output[iZone] = new CBaselineOutput(config_container[iZone], geometry_container[iZone][INST_0]->GetnDim(),
                                            solver_container[iZone][INST_0]);
        output[iZone]->PreprocessVolumeOutput(config_container[iZone]);
        output[iZone]->PreprocessHistoryOutput(config_container[iZone], false);
      }

      /*--- Loop over the whole time domain ---*/
      while (TimeIter < driver_config->GetnTime_Iter()) {
        /*--- Check if the maximum time has been surpassed. ---*/
        Physical_t = (TimeIter + 1) * Physical_dt;
        if (Physical_t >= driver_config->GetMax_Time()) StopCalc = true;

        if ((TimeIter + 1 == driver_config->GetnTime_Iter()) ||  // The last time iteration
            (StopCalc) ||                                        // We have surpassed the requested time
            ((TimeIter == 0) || (TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) ==
                                 0))  // The iteration has been requested
        ) {
          if (rank == MASTER_NODE)
            cout << "Writing the volume solution for time step " << TimeIter << ", t = " << Physical_t << " s ."
                 << endl;

          /*--- Load the restart for all the zones. ---*/
          for (iZone = 0; iZone < nZone; iZone++) {
            /*--- Set the current iteration number in the config class. ---*/
            config_container[iZone]->SetTimeIter(TimeIter);
            /*--- So far, only enabled for single-instance problems ---*/
            config_container[iZone]->SetiInst(INST_0);
            solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone],
                                                         config_container[iZone], TimeIter, true);
          }

          for (iZone = 0; iZone < nZone; iZone++) {
            WriteFiles(config_container[iZone], geometry_container[iZone][INST_0], &solver_container[iZone][INST_0],
                       output[iZone], TimeIter);
          }
        }

        TimeIter++;
        if (StopCalc) break;
      }
    } else {
      /*--- Steady simulation: merge the solution files for each zone. ---*/
      for (iZone = 0; iZone < nZone; iZone++) {
        config_container[iZone]->SetiInst(INST_0);
        /*--- Definition of the solution class ---*/
        solver_container[iZone][INST_0] =
            new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
        solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone],
                                                     config_container[iZone], 0, true);
        output[iZone] = new CBaselineOutput(config_container[iZone], geometry_container[iZone][INST_0]->GetnDim(),
                                            solver_container[iZone][INST_0]);
        output[iZone]->PreprocessVolumeOutput(config_container[iZone]);
        output[iZone]->PreprocessHistoryOutput(config_container[iZone], false);
      }
      for (iZone = 0; iZone < nZone; iZone++) {
        WriteFiles(config_container[iZone], geometry_container[iZone][INST_0], &solver_container[iZone][INST_0],
                   output[iZone], 0);
      }
    }

  } else if (fsi) {
    if (nZone < 2) {
      SU2_MPI::Error("For multizone computations, please add the number of zones as a second argument for SU2_SOL.",
                     CURRENT_FUNCTION);
    }

    su2double Physical_dt, Physical_t;
    unsigned long TimeIter = 0, TimeIterFlow = 0, TimeIterFEM = 0;
    bool StopCalc = false;
    bool SolutionInstantiatedFlow = false, SolutionInstantiatedFEM = false;

    /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
    if (config_container[ZONE_0]->GetRestart()) {
      TimeIterFlow = config_container[ZONE_0]->GetRestart_Iter();
      TimeIterFEM = config_container[ZONE_1]->GetRestart_Iter();
      if (TimeIterFlow != TimeIterFEM) {
        SU2_MPI::Error("For multizone computations, please add the number of zones as a second argument for SU2_SOL.",
                       CURRENT_FUNCTION);
      } else {
        TimeIter = TimeIterFlow;
      }
    }

    while (TimeIter < config_container[ZONE_0]->GetnTime_Iter()) {
      /*--- Check several conditions in order to merge the correct time step files. ---*/

      Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
      Physical_t = (TimeIter + 1) * Physical_dt;
      if (Physical_t >= config_container[ZONE_0]->GetTotal_UnstTime()) StopCalc = true;

      if (((TimeIter + 1 == config_container[ZONE_0]->GetnTime_Iter()) ||
           ((TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0) && (TimeIter != 0) &&
            (config_container[ZONE_0]->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_1ST) &&
            (config_container[ZONE_0]->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_2ND)) ||
           (StopCalc) ||
           (((config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
             (config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) &&
            ((TimeIter == 0) || (TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0))))

          &&

          ((TimeIter + 1 == config_container[ZONE_1]->GetnTime_Iter()) || (StopCalc) ||
           ((config_container[ZONE_1]->GetTime_Domain()) &&
            ((TimeIter == 0) || (TimeIter % config_container[ZONE_1]->GetVolumeOutputFrequency(0) == 0))))

      ) {
        /*--- Set the current iteration number in the config class. ---*/
        config_container[ZONE_0]->SetTimeIter(TimeIter);
        config_container[ZONE_1]->SetTimeIter(TimeIter);

        /*--- Read in the restart file for this time step ---*/

        /*--- For the fluid zone (ZONE_0) ---*/
        /*--- Either instantiate the solution class or load a restart file. ---*/
        if (!SolutionInstantiatedFlow &&
            (TimeIter == 0 ||
             ((config_container[ZONE_0]->GetRestart() &&
               (SU2_TYPE::Int(TimeIter) == SU2_TYPE::Int(config_container[ZONE_0]->GetRestart_Iter()))) ||
              TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0 ||
              TimeIter + 1 == config_container[ZONE_0]->GetnTime_Iter()))) {
          solver_container[ZONE_0][INST_0] =
              new CBaselineSolver(geometry_container[ZONE_0][INST_0], config_container[ZONE_0]);
          output[ZONE_0] = new CBaselineOutput(config_container[ZONE_0], geometry_container[ZONE_0][INST_0]->GetnDim(),
                                               solver_container[ZONE_0][INST_0]);
          output[ZONE_0]->PreprocessVolumeOutput(config_container[ZONE_0]);
          output[ZONE_0]->PreprocessHistoryOutput(config_container[ZONE_0], false);

          SolutionInstantiatedFlow = true;
        }
        solver_container[ZONE_0][INST_0]->LoadRestart_FSI(geometry_container[ZONE_0][INST_0], config_container[ZONE_0],
                                                          TimeIter);

        /*--- For the structural zone (ZONE_1) ---*/
        /*--- Either instantiate the solution class or load a restart file. ---*/
        /*--- Either instantiate the solution class or load a restart file. ---*/
        if (!SolutionInstantiatedFEM &&
            (TimeIter == 0 ||
             ((config_container[ZONE_1]->GetRestart() &&
               (SU2_TYPE::Int(TimeIter) == SU2_TYPE::Int(config_container[ZONE_1]->GetRestart_Iter()))) ||
              TimeIter % config_container[ZONE_1]->GetVolumeOutputFrequency(0) == 0 ||
              TimeIter + 1 == config_container[ZONE_1]->GetnTime_Iter()))) {
          solver_container[ZONE_1][INST_0] =
              new CBaselineSolver(geometry_container[ZONE_1][INST_0], config_container[ZONE_1]);
          output[ZONE_1] = new CBaselineOutput(config_container[ZONE_1], geometry_container[ZONE_1][INST_0]->GetnDim(),
                                               solver_container[ZONE_1][INST_0]);
          output[ZONE_1]->PreprocessVolumeOutput(config_container[ZONE_1]);
          SolutionInstantiatedFEM = true;
        }
        solver_container[ZONE_1][INST_0]->LoadRestart_FSI(geometry_container[ZONE_1][INST_0], config_container[ZONE_1],
                                                          TimeIter);

        if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << TimeIter << "." << endl;
        for (iZone = 0; iZone < nZone; iZone++) {
          WriteFiles(config_container[iZone], geometry_container[iZone][INST_0], &solver_container[iZone][INST_0],
                     output[iZone], TimeIter);
        }
      }

      TimeIter++;
      if (StopCalc) break;
    }

  } else if (fem_solver) {
    if (config->GetTime_Domain()) {
      /*--- Unsteady DG simulation: merge all unsteady time steps. First,
       find the frequency and total number of files to write. ---*/

      su2double Physical_dt, Physical_t;
      unsigned long TimeIter = 0;
      bool StopCalc = false;
      bool* SolutionInstantiated = new bool[nZone];

      for (iZone = 0; iZone < nZone; iZone++) SolutionInstantiated[iZone] = false;

      /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
      if (config_container[ZONE_0]->GetTime_Domain() && config_container[ZONE_0]->GetRestart())
        TimeIter = config_container[ZONE_0]->GetRestart_Iter();

      while (TimeIter < config_container[ZONE_0]->GetnTime_Iter()) {
        /*--- Check several conditions in order to merge the correct time step files. ---*/
        Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
        Physical_t = (TimeIter + 1) * Physical_dt;
        if (Physical_t >= config_container[ZONE_0]->GetTotal_UnstTime()) StopCalc = true;

        if ((TimeIter + 1 == config_container[ZONE_0]->GetnTime_Iter()) ||
            ((TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0) && (TimeIter != 0) &&
             !(config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING)) ||
            (StopCalc) ||
            ((config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING) &&
             ((TimeIter == 0) || (TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0)))) {
          /*--- Read in the restart file for this time step ---*/
          for (iZone = 0; iZone < nZone; iZone++) {
            /*--- Set the current iteration number in the config class. ---*/
            config_container[iZone]->SetTimeIter(TimeIter);

            /*--- Either instantiate the solution class or load a restart file. ---*/
            if (!SolutionInstantiated[iZone] &&
                (TimeIter == 0 || (config_container[ZONE_0]->GetRestart() &&
                                   ((long)TimeIter == SU2_TYPE::Int(config_container[ZONE_0]->GetRestart_Iter()) ||
                                    TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0 ||
                                    TimeIter + 1 == config_container[ZONE_0]->GetnTime_Iter())))) {
              solver_container[iZone][INST_0] =
                  new CBaselineSolver_FEM(geometry_container[iZone][INST_0], config_container[iZone]);
              output[iZone] =
                  new CBaselineOutput(config_container[ZONE_0], geometry_container[ZONE_0][INST_0]->GetnDim(),
                                      solver_container[ZONE_0][INST_0]);
              output[iZone]->PreprocessVolumeOutput(config_container[ZONE_0]);
              output[iZone]->PreprocessHistoryOutput(config_container[ZONE_0], false);
              SolutionInstantiated[iZone] = true;
            }
            solver_container[iZone][INST_0]->LoadRestart(&geometry_container[iZone][INST_0], &solver_container[iZone],
                                                         config_container[iZone], (int)TimeIter, true);
          }

          if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << TimeIter << "." << endl;

          for (iZone = 0; iZone < nZone; iZone++) {
            WriteFiles(config_container[iZone], geometry_container[iZone][INST_0], &solver_container[iZone][INST_0],
                       output[iZone], TimeIter);
          }
        }

        TimeIter++;
        if (StopCalc) break;
      }

    } else {
      /*--- Steady simulation: merge the single solution file. ---*/

      for (iZone = 0; iZone < nZone; iZone++) {
        /*--- Definition of the solution class ---*/

        solver_container[iZone][INST_0] =
            new CBaselineSolver_FEM(geometry_container[iZone][INST_0], config_container[iZone]);
        output[iZone] = new CBaselineOutput(config_container[ZONE_0], geometry_container[ZONE_0][INST_0]->GetnDim(),
                                            solver_container[ZONE_0][INST_0]);
        output[iZone]->PreprocessVolumeOutput(config_container[ZONE_0]);
        output[iZone]->PreprocessHistoryOutput(config_container[ZONE_0], false);
        solver_container[iZone][INST_0]->LoadRestart(&geometry_container[iZone][INST_0], &solver_container[iZone],
                                                     config_container[iZone], 0, true);
      }

      for (iZone = 0; iZone < nZone; iZone++) {
        WriteFiles(config_container[iZone], geometry_container[iZone][INST_0], &solver_container[iZone][INST_0],
                   output[iZone], 0);
      }
    }

  } else {
    if (config_container[ZONE_0]->GetTime_Domain()) {
      /*--- Unsteady simulation: merge all unsteady time steps. First,
       find the frequency and total number of files to write. ---*/

      su2double Physical_dt, Physical_t;
      unsigned long TimeIter = 0;
      bool StopCalc = false;
      bool* SolutionInstantiated = new bool[nZone];

      for (iZone = 0; iZone < nZone; iZone++) SolutionInstantiated[iZone] = false;

      /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
      if (config_container[ZONE_0]->GetTime_Domain() && config_container[ZONE_0]->GetRestart())
        TimeIter = config_container[ZONE_0]->GetRestart_Iter();

      while (TimeIter < config_container[ZONE_0]->GetnTime_Iter()) {
        /*--- Check several conditions in order to merge the correct time step files. ---*/
        Physical_dt = config_container[ZONE_0]->GetTime_Step();
        Physical_t = (TimeIter + 1) * Physical_dt;
        if (Physical_t >= config_container[ZONE_0]->GetMax_Time()) StopCalc = true;

        if ((TimeIter + 1 == config_container[ZONE_0]->GetnTime_Iter()) ||
            ((TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0) && (TimeIter != 0) &&
             (config_container[ZONE_0]->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_1ST) &&
             (config_container[ZONE_0]->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_2ND)) ||
            (StopCalc) ||
            (((config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
              (config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) &&
             ((TimeIter == 0) || (TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0)))) {
          /*--- Read in the restart file for this time step ---*/
          for (iZone = 0; iZone < nZone; iZone++) {
            /*--- Set the current iteration number in the config class. ---*/
            config_container[iZone]->SetTimeIter(TimeIter);

            /*--- Either instantiate the solution class or load a restart file. ---*/
            if (!SolutionInstantiated[iZone] &&
                (TimeIter == 0 || (config_container[ZONE_0]->GetRestart() &&
                                   ((long)TimeIter == SU2_TYPE::Int(config_container[ZONE_0]->GetRestart_Iter()) ||
                                    TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0 ||
                                    TimeIter + 1 == config_container[ZONE_0]->GetnTime_Iter())))) {
              solver_container[iZone][INST_0] =
                  new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
              output[iZone] = new CBaselineOutput(config_container[iZone], geometry_container[iZone][INST_0]->GetnDim(),
                                                  solver_container[iZone][INST_0]);
              output[iZone]->PreprocessVolumeOutput(config_container[iZone]);
              output[iZone]->PreprocessHistoryOutput(config_container[iZone], false);

              SolutionInstantiated[iZone] = true;
            }
            config_container[iZone]->SetiInst(INST_0);
            solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone],
                                                         config_container[iZone], TimeIter, true);
          }

          if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << TimeIter << "." << endl;

          for (iZone = 0; iZone < nZone; iZone++) {
            WriteFiles(config_container[iZone], geometry_container[iZone][INST_0], &solver_container[iZone][INST_0],
                       output[iZone], TimeIter);
          }
        }

        TimeIter++;
        if (StopCalc) break;
      }

    }

    else if (config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE) {
      /*--- Read in the restart file for this time step ---*/
      for (iZone = 0; iZone < nZone; iZone++) {
        for (iInst = 0; iInst < nInst[iZone]; iInst++) {
          config_container[iZone]->SetiInst(iInst);
          config_container[iZone]->SetTimeIter(iInst);

          /*--- Either instantiate the solution class or load a restart file. ---*/
          solver_container[iZone][iInst] =
              new CBaselineSolver(geometry_container[iZone][iInst], config_container[iZone]);
          solver_container[iZone][iInst]->LoadRestart(geometry_container[iZone], &solver_container[iZone],
                                                      config_container[iZone], iInst, true);
          output[iZone] = new CBaselineOutput(config_container[iZone], geometry_container[iZone][iInst]->GetnDim(),
                                              solver_container[iZone][iInst]);
          output[iZone]->PreprocessVolumeOutput(config_container[iZone]);
          output[iZone]->PreprocessHistoryOutput(config_container[iZone], false);

          /*--- Print progress in solution writing to the screen. ---*/
          if (rank == MASTER_NODE) {
            cout << "Storing the volume solution for time instance " << iInst << "." << endl;
          }

          WriteFiles(config_container[iZone], geometry_container[iZone][iInst], &solver_container[iZone][iInst],
                     output[iZone], iInst);
        }
      }

    }

    else if (config_container[ZONE_0]->GetTime_Domain()) {
      /*--- Dynamic simulation: merge all unsteady time steps. First,
       find the frequency and total number of files to write. ---*/

      su2double Physical_dt, Physical_t;
      unsigned long TimeIter = 0;
      bool StopCalc = false;
      bool SolutionInstantiated = false;

      /*--- Check for an dynamic restart (structural analysis). Update ExtIter if necessary. ---*/
      if (config_container[ZONE_0]->GetKind_Solver() == MAIN_SOLVER::FEM_ELASTICITY &&
          config_container[ZONE_0]->GetRestart())
        TimeIter = config_container[ZONE_0]->GetRestart_Iter();

      while (TimeIter < config_container[ZONE_0]->GetnTime_Iter()) {
        /*--- Check several conditions in order to merge the correct time step files. ---*/
        /*--- If the solver is structural, the total and delta_t are obtained from different functions. ---*/

        Physical_dt = config_container[ZONE_0]->GetTime_Step();
        Physical_t = (TimeIter + 1) * Physical_dt;
        if (Physical_t >= config_container[ZONE_0]->GetMax_Time()) StopCalc = true;

        if ((TimeIter + 1 == config_container[ZONE_0]->GetnTime_Iter()) || (StopCalc) ||
            ((config_container[ZONE_0]->GetTime_Domain()) &&
             ((TimeIter == 0) || (TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0)))) {
          /*--- Set the current iteration number in the config class. ---*/
          config_container[ZONE_0]->SetTimeIter(TimeIter);

          /*--- Read in the restart file for this time step ---*/
          for (iZone = 0; iZone < nZone; iZone++) {
            /*--- Either instantiate the solution class or load a restart file. ---*/
            if (!SolutionInstantiated &&
                (TimeIter == 0 ||
                 ((config_container[ZONE_0]->GetRestart() &&
                   (SU2_TYPE::Int(TimeIter) == SU2_TYPE::Int(config_container[ZONE_0]->GetRestart_Iter()))) ||
                  TimeIter % config_container[ZONE_0]->GetVolumeOutputFrequency(0) == 0 ||
                  TimeIter + 1 == config_container[ZONE_0]->GetnTime_Iter()))) {
              solver_container[iZone][INST_0] =
                  new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
              output[iZone] = new CBaselineOutput(config_container[iZone], geometry_container[iZone][INST_0]->GetnDim(),
                                                  solver_container[iZone][INST_0]);
              output[iZone]->PreprocessVolumeOutput(config_container[iZone]);
              output[iZone]->PreprocessHistoryOutput(config_container[iZone], false);

              SolutionInstantiated = true;
            }
            config_container[iZone]->SetiInst(INST_0);
            solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone],
                                                         config_container[iZone], TimeIter, true);
          }

          if (rank == MASTER_NODE) cout << "Writing the volume solution for time step " << TimeIter << "." << endl;
          for (iZone = 0; iZone < nZone; iZone++) {
            WriteFiles(config_container[iZone], geometry_container[iZone][INST_0], &solver_container[iZone][INST_0],
                       output[iZone], TimeIter);
          }
        }

        TimeIter++;
        if (StopCalc) break;
      }

    }

    else {
      /*--- Steady simulation: merge the single solution file. ---*/

      for (iZone = 0; iZone < nZone; iZone++) {
        config_container[iZone]->SetiInst(INST_0);
        /*--- Definition of the solution class ---*/
        solver_container[iZone][INST_0] =
            new CBaselineSolver(geometry_container[iZone][INST_0], config_container[iZone]);
        solver_container[iZone][INST_0]->LoadRestart(geometry_container[iZone], &solver_container[iZone],
                                                     config_container[iZone], 0, true);
        output[iZone] = new CBaselineOutput(config_container[iZone], geometry_container[iZone][INST_0]->GetnDim(),
                                            solver_container[iZone][INST_0]);
        output[iZone]->PreprocessVolumeOutput(config_container[iZone]);
        output[iZone]->PreprocessHistoryOutput(config_container[iZone], false);
      }
      for (iZone = 0; iZone < nZone; iZone++) {
        WriteFiles(config_container[iZone], geometry_container[iZone][INST_0], &solver_container[iZone][INST_0],
                   output[iZone], 0);
      }
    }
  }

  delete config;
  config = nullptr;

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Finalize Solver -------------------------" << endl;

  if (geometry_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (geometry_container[iZone][iInst] != nullptr) {
          delete geometry_container[iZone][iInst];
        }
      }
      if (geometry_container[iZone] != nullptr) delete[] geometry_container[iZone];
    }
    delete[] geometry_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

  if (solver_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        if (solver_container[iZone][iInst] != nullptr) {
          delete solver_container[iZone][iInst];
        }
      }
      if (solver_container[iZone] != nullptr) delete[] solver_container[iZone];
    }
    delete[] solver_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CSolver class." << endl;

  if (config_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != nullptr) {
        delete config_container[iZone];
      }
    }
    delete[] config_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

  if (output != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (output[iZone] != nullptr) {
        delete output[iZone];
      }
    }
    delete[] output;
  }
  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime - StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on " << size;
    if (size == 1)
      cout << " core." << endl;
    else
      cout << " cores." << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << endl << "------------------------- Exit Success (SU2_SOL) ------------------------" << endl << endl;

  /*--- Finalize MPI parallelization ---*/
  SU2_MPI::Finalize();

  return EXIT_SUCCESS;
}

void WriteFiles(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput* output,
                unsigned long TimeIter) {
  /*--- Load history data (volume output might require some values) --- */

  output->SetHistoryOutput(geometry, solver_container, config, TimeIter, 0, 0);

  /*--- Load the data --- */

  output->LoadData(geometry, config, solver_container);

  /*--- Set the filenames ---*/

  output->SetVolumeFilename(config->GetVolume_FileName());

  output->SetSurfaceFilename(config->GetSurfCoeff_FileName());

  for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++) {
    auto FileFormat = config->GetVolumeOutputFiles();
    if (FileFormat[iFile] != OUTPUT_TYPE::RESTART_ASCII && FileFormat[iFile] != OUTPUT_TYPE::RESTART_BINARY &&
        FileFormat[iFile] != OUTPUT_TYPE::CSV)
      output->WriteToFile(config, geometry, FileFormat[iFile]);
  }
}
