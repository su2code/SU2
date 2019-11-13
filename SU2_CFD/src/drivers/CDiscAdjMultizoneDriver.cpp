/*!
 * \file CDiscAdjMultizoneDriver.cpp
 * \brief The main subroutines for driving adjoint multi-zone problems
 * \author O. Burghardt, T. Albring, R. Sanchez
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../include/drivers/CDiscAdjMultizoneDriver.hpp"

CDiscAdjMultizoneDriver::CDiscAdjMultizoneDriver(char* confFile,
                                                 unsigned short val_nZone,
                                                 SU2_Comm MPICommunicator)

                        : CMultizoneDriver(confFile, val_nZone, MPICommunicator) {

  retape = !config_container[ZONE_0]->GetFull_Tape();

  RecordingState = NONE;

  direct_nInst  = new unsigned short[nZone];
  nInnerIter    = new unsigned short[nZone];


  for (iZone = 0; iZone < nZone; iZone++) {

    direct_nInst[iZone] = 1;
    nInnerIter[iZone]   = config_container[iZone]->GetnInner_Iter();
  }

  direct_iteration = new CIteration**[nZone];
  direct_output = new COutput*[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {

    direct_iteration[iZone] = new CIteration*[direct_nInst[iZone]];

    for(iInst = 0; iInst < direct_nInst[iZone]; iInst++) {

      switch (config_container[iZone]->GetKind_Solver()) {

        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          direct_iteration[iZone][iInst] = new CFluidIteration(config_container[iZone]);
          direct_output[iZone] = new CFlowCompOutput(config_container[iZone], nDim);
          break;
        case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
          direct_iteration[iZone][iInst] = new CFluidIteration(config_container[iZone]);
          direct_output[iZone] = new CFlowIncOutput(config_container[iZone], nDim);
          break;
        case DISC_ADJ_HEAT:
          direct_iteration[iZone][iInst] = new CHeatIteration(config_container[iZone]);
          direct_output[iZone] = new CHeatOutput(config_container[iZone], nDim);
          break;
        case DISC_ADJ_FEM:
          SU2_MPI::Error("There is no multizone discrete adjoint functionality for problems including elasticity yet.",
                         CURRENT_FUNCTION);
          direct_iteration[iZone][iInst] = new CFEAIteration(config_container[iZone]);
          direct_output[iZone] = new CElasticityOutput(config_container[iZone], nDim);
          break;

        default:
          SU2_MPI::Error("There is no discrete adjoint functionality for one of the specified solvers yet.",
                         CURRENT_FUNCTION);
      }
    }

    direct_output[iZone]->PreprocessHistoryOutput(config_container[iZone], false);
  }
}

CDiscAdjMultizoneDriver::~CDiscAdjMultizoneDriver(){

  for (iZone = 0; iZone < nZone; iZone++){
    for (iInst = 0; iInst < direct_nInst[iZone]; iInst++){
      delete direct_iteration[iZone][iInst];
    }
    delete [] direct_iteration[iZone];
  }

  delete[] direct_iteration;
  delete[] direct_nInst;

}

void CDiscAdjMultizoneDriver::StartSolver() {

  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;

  if (rank == MASTER_NODE){
    cout << endl << "Simulation Run using the Discrete Adjoint Multizone Driver" << endl;

    if (driver_config->GetTime_Domain())
      SU2_MPI::Error("The discrete adjoint multizone driver is not ready for unsteady computations yet.",
                     CURRENT_FUNCTION);
  }

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the value of the external iteration to TimeIter. -------------------------------------*/
    /*--- TODO: This should be generalised for an homogeneous criteria throughout the code. --------*/
    config_container[iZone]->SetTimeIter(0);
    
  }
  
  /*--- We directly start the (steady-state) discrete adjoint computation. ---*/

  Run();
  
  /*--- Output the solution in files. ---*/

  Output(TimeIter);

}

void CDiscAdjMultizoneDriver::Run() {

  bool            checkSensitivity  = false;
  unsigned short  jZone             = 0,
                  wrt_sol_freq      = config_container[ZONE_0]->GetVolume_Wrt_Freq();
  unsigned long   nIter             = 0,
                  iOuter_Iter       = 0;

  nIter = driver_config->GetnOuter_Iter();

  for (iZone = 0; iZone < nZone; iZone++) {

    iteration_container[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container,
                                                   solver_container, numerics_container, config_container, surface_movement,
                                                   grid_movement, FFDBox, iZone, INST_0);

    /*--- Set BGS_Solution_k to Solution. ---*/

    Set_BGSSolution(iZone);
  }

  /*--- Loop over the number of outer iterations. ---*/

  for (iOuter_Iter = 0; iOuter_Iter < nIter; iOuter_Iter++) {

    for (iZone = 0; iZone < nZone; iZone++) {
      config_container[iZone]->SetOuterIter(iOuter_Iter);
      driver_config->SetOuterIter(iOuter_Iter);
    }


    /*--- For the adjoint iteration we need the derivatives of the iteration function with
     *    respect to the state (and possibly the mesh coordinate) variables.
     *    Since these derivatives do not change in the steady state case we only have to record
     *    if the current recording is different from them.
     *
     *    To set the tape appropriately, the following recording methods are provided:
     *    (1) NONE: All information from a previous recording is removed.
     *    (2) FLOW_CONS_VARS: State variables of all solvers in a zone as input.
     *    (3) MESH_COORDS: Mesh coordinates as input.
     *    (4) COMBINED: Mesh coordinates and state variables as input.
     *
     *    By default, all (state and mesh coordinate variables) will be declared as output,
     *    since it does not change the computational effort. ---*/


    /*--- If we want to set up zone-specific tapes later on,
     *    we just record the objective function section here.
     *    If not, the whole tape of a coupled run will be created. ---*/

    if (retape) {
      SetRecording(NONE, Kind_Tape::FULL_TAPE, ZONE_0);
      SetRecording(FLOW_CONS_VARS, Kind_Tape::OBJECTIVE_FUNCTION_TAPE, ZONE_0);
    }
    else if (RecordingState != FLOW_CONS_VARS) {

      SetRecording(NONE, Kind_Tape::FULL_TAPE, ZONE_0);
      SetRecording(FLOW_CONS_VARS, Kind_Tape::FULL_TAPE, ZONE_0);
    }


    SetExternal_Zero();

    /*-- Start loop over zones. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {

      if (retape) {
        SetRecording(NONE, Kind_Tape::FULL_TAPE, ZONE_0);
        SetRecording(FLOW_CONS_VARS, Kind_Tape::ZONE_SPECIFIC_TAPE, iZone);
      }

      /*--- Evaluate the objective function gradient w.r.t. to solution contributions from iZone.
       *    (We always initialize from Solution_BGS_k and extract to Solution.) ---*/

      AD::ClearAdjoints();

      SetAdj_ObjFunction();

      AD::ComputeAdjoint(OBJECTIVE_FUNCTION, START);

      iteration_container[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container,
                                                  solver_container, numerics_container, config_container,
                                                  surface_movement, grid_movement, FFDBox, iZone, INST_0);

      /*--- Add the objective function contribution to the external contributions for solvers in iZone. ---*/

      Add_Solution_To_ExternalOld(iZone);

      /*--- Inner loop to allow for multiple adjoint updates with respect to solvers in iZone. ---*/

      for (unsigned short iInnerIter = 0; iInnerIter < nInnerIter[iZone]; iInnerIter++) {

        /*--- Evaluate the tape section belonging to solvers in iZone. ---*/

        ComputeAdjoints(iZone);

        /*--- Extracting adjoints for solvers in iZone w.r.t. to outputs in iZone (diagonal part). ---*/

        iteration_container[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container,
                                                    solver_container, numerics_container, config_container,
                                                    surface_movement, grid_movement, FFDBox, iZone, INST_0);

        /*--- Add off-diagonal contribution (including the OF gradient) to Solution for next inner evaluation. ---*/

        Add_ExternalOld_To_Solution(iZone);

        /*--- Print out the convergence data to screen and history file ---*/

        output_container[iZone]->SetHistory_Output(geometry_container[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0],
                                                   config_container[iZone], config_container[iZone]->GetTimeIter(), iOuter_Iter, iInnerIter);
      }

      /*--- Off-diagonal (coupling term) update for next outer evaluation. ---*/

      for (jZone = 0; jZone < nZone; jZone++) {

        if (jZone != iZone) {

          /*--- Extracting adjoints for solvers in jZone w.r.t. to the output of all solvers in iZone,
           *    that is, for the cases iZone != jZone we are evaluating cross derivatives between zones. ---*/

          iteration_container[jZone][INST_0]->Iterate(output_container[jZone], integration_container, geometry_container,
                                                      solver_container, numerics_container, config_container,
                                                      surface_movement, grid_movement, FFDBox, jZone, INST_0);

          /*--- Add the cross derivatives from iZone<-jZone dependencies to the External vector. ---*/

          Add_Solution_To_External(jZone);
        }
      }

      // TODO: Add an option to _already_ update off-diagonals here (i.e., in the zone-loop)

      /*--- Save Solution to Solution_BGS and compute residual from Solution_BGS and Solution_BGS_k. ---*/

      SetResidual_BGS(iZone);

      /*--- Save Solution to Solution_BGS_k for a next outer iteration.
       *    (Solution might be overwritten when entering another zone because of cross derivatives.) ---*/

      Set_BGSSolution(iZone);
        
      /*--- make sure that everything is loaded into the output container ---*/
      
      output_container[iZone]->SetHistory_Output(geometry_container[iZone][INST_0][MESH_0],solver_container[iZone][INST_0][MESH_0], config_container[iZone]);
      
    }

    /*--- Now all coupling terms are summed up, set External_Old to External for next outer iteration. ---*/

    Set_OldExternal();
    
    /*--- Set the multizone output. ---*/

    driver_output->SetMultizoneHistory_Output(output_container, config_container, driver_config,
                                              driver_config->GetTimeIter(), driver_config->GetOuterIter());

    /*--- Check for convergence. ---*/

    StopCalc = driver_output->GetConvergence();

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    /*--- Compute the geometrical sensitivities and write them to file. ---*/

    checkSensitivity = ((iOuter_Iter+1 >= nIter) ||
                        (iOuter_Iter % wrt_sol_freq == 0));

    if (checkSensitivity || StopCalc){

      /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
       *    as argument ensures that all information from a previous recording is removed. ---*/

      SetRecording(NONE, Kind_Tape::FULL_TAPE, ZONE_0);

      /*--- Store the computational graph of one direct iteration with the mesh coordinates as input. ---*/

      SetRecording(MESH_COORDS, Kind_Tape::FULL_TAPE, ZONE_0);

      /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
       *    of the current iteration. The values are passed to the AD tool. ---*/

      for (iZone = 0; iZone < nZone; iZone++) {

        iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container,
                                                              config_container, iZone, INST_0);
      }

      /*--- Initialize the adjoint of the objective function with 1.0. ---*/

      SetAdj_ObjFunction();

      /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

      AD::ComputeAdjoint();

      /*--- Extract the computed sensitivity values. ---*/

      for (iZone = 0; iZone < nZone; iZone++) {

        switch (config_container[iZone]->GetKind_Solver()) {

          case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
            solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],
                                                                                 NULL, config_container[iZone]);
            break;
          case DISC_ADJ_HEAT:
            solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],
                                                                                 NULL, config_container[iZone]);
            break;
          case DISC_ADJ_FEM:
            solver_container[iZone][INST_0][MESH_0][ADJFEA_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],
                                                                                NULL, config_container[iZone]);
            break;

          default:
            cout << "WARNING: Sensitivities not set for one of the specified discrete adjoint solvers!" << endl;
            break;
        }
      }

      /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

      AD::ClearAdjoints();

      for (iZone = 0; iZone < nZone; iZone++) {

        output_container[iZone]->SetResult_Files(geometry_container[iZone][INST_0][MESH_0],
                                                 config_container[iZone],
                                                 solver_container[iZone][INST_0][MESH_0], iOuter_Iter, StopCalc);
      }
    }

    if (StopCalc) break;
  }
}

void CDiscAdjMultizoneDriver::SetRecording(unsigned short kind_recording, Kind_Tape tape_type, unsigned short record_zone) {

  unsigned short iZone, jZone, iSol, UpdateMesh;
  unsigned long ExtIter = 0;
  bool DeformMesh       = false;

  AD::Reset();

  /*--- Prepare for recording by resetting the flow solution to the initial converged solution---*/

  for(iZone = 0; iZone < nZone; iZone++) {

    for (iSol=0; iSol < MAX_SOLS; iSol++){
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL) {
        if (solver_container[iZone][INST_0][MESH_0][iSol]->GetAdjoint()) {
          solver_container[iZone][INST_0][MESH_0][iSol]->SetRecording(geometry_container[iZone][INST_0][MESH_0],
                                                                      config_container[iZone]);
        }
      }
    }
  }

  /*---Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if(kind_recording != NONE) {

    if (rank == MASTER_NODE && kind_recording == FLOW_CONS_VARS) {
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << "Storing computational graph." << endl;
    }

    AD::StartRecording();

    AD::Push_TapePosition();

    for (iZone = 0; iZone < nZone; iZone++) {

      iteration_container[iZone][INST_0]->RegisterInput(solver_container, geometry_container,
                                                        config_container, iZone, INST_0, kind_recording);
    }
  }

  AD::Push_TapePosition();

  for (iZone = 0; iZone < nZone; iZone++) {

    iteration_container[iZone][INST_0]->SetDependencies(solver_container, geometry_container, numerics_container,
                                                        config_container, iZone, INST_0, kind_recording);
  }

  AD::Push_TapePosition();

  /*--- Extract the objective function and store it --- */

  SetObjFunction(kind_recording);

  AD::Push_TapePosition();

  if (tape_type != Kind_Tape::OBJECTIVE_FUNCTION_TAPE) {

    /*--- We do the communication here to not derive wrt updated boundary data. ---*/

    for(iZone = 0; iZone < nZone; iZone++) {

      /*--- In principle, the mesh does not need to be updated ---*/
      UpdateMesh = 0;

      /*--- Transfer from all the remaining zones ---*/
      for (jZone = 0; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone && interface_container[iZone][jZone] != NULL){
          DeformMesh = Transfer_Data(jZone, iZone);
          if (DeformMesh) UpdateMesh+=1;
        }
      }
      /*--- If a mesh update is required due to the transfer of data ---*/
      if (UpdateMesh > 0) DynamicMeshUpdate(iZone, ExtIter);
    }

    AD::Push_TapePosition();

    for(iZone = 0; iZone < nZone; iZone++) {

      AD::Push_TapePosition();

      if (tape_type == Kind_Tape::ZONE_SPECIFIC_TAPE) {
        if (iZone == record_zone) {
          DirectIteration(iZone, kind_recording);
        }
      }
      else {
        DirectIteration(iZone, kind_recording);
      }

      iteration_container[iZone][INST_0]->RegisterOutput(solver_container, geometry_container,
                                                         config_container, output_container[iZone], iZone, INST_0);

      AD::Push_TapePosition();
    }
  }

  if (rank == MASTER_NODE && kind_recording == FLOW_CONS_VARS) {

    if(config_container[record_zone]->GetWrt_AD_Statistics()) {

      AD::PrintStatistics();
    }

    cout << "-------------------------------------------------------------------------" << endl << endl;
  }

  AD::StopRecording();

  RecordingState = kind_recording;
}

void CDiscAdjMultizoneDriver::DirectIteration(unsigned short iZone, unsigned short kind_recording) {

  /*--- Do one iteration of the direct solver ---*/
  direct_iteration[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container, solver_container,
                                              numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

  /*--- Iterate the zone as a block a single time ---*/
  direct_iteration[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container, solver_container,
                                           numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

  Corrector(iZone);

  /*--- Print residuals in the first iteration ---*/

  if (rank == MASTER_NODE && kind_recording == FLOW_CONS_VARS) {

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES:
      case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES:
        cout << " Zone " << iZone << " (flow) - log10[RMS Solution_0]: "
             << log10(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetRes_RMS(0)) << endl;
        break;
      case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS:
        cout << " Zone " << iZone << " (flow) - log10[RMS Solution_0]: "
             << log10(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetRes_RMS(0)) << endl;

        if (!config_container[iZone]->GetFrozen_Visc_Disc()) {

          cout << " Zone " << iZone << " (turb) - log10[RMS k]         : "
               << log10(solver_container[iZone][INST_0][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
        }
        break;
      case DISC_ADJ_HEAT:
        cout << " Zone " << iZone << " (heat) - log10[RMS Solution_0]: "
             << log10(solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetRes_RMS(0)) << endl;
        break;
    }
  }
}

void CDiscAdjMultizoneDriver::SetObjFunction(unsigned short kind_recording) {

  ObjFunc = 0.0;
  su2double Weight_ObjFunc;

  unsigned short iZone, iMarker_Analyze, nMarker_Analyze;


  /*--- Call objective function calculations. ---*/

  for (iZone = 0; iZone < nZone; iZone++){

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
        solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Pressure_Forces(geometry_container[iZone][INST_0][MESH_0],
                                                                           config_container[iZone]);
        solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Momentum_Forces(geometry_container[iZone][INST_0][MESH_0],
                                                                           config_container[iZone]);
        solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Friction_Forces(geometry_container[iZone][INST_0][MESH_0],
                                                                           config_container[iZone]);

        if(config_container[iZone]->GetWeakly_Coupled_Heat()) {

          solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->Heat_Fluxes(geometry_container[iZone][INST_0][MESH_0],
                                                                         solver_container[iZone][INST_0][MESH_0], config_container[iZone]);
        }
        solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Evaluate_ObjFunc(config_container[iZone]);
        break;
      case DISC_ADJ_HEAT:
        solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->Heat_Fluxes(geometry_container[iZone][INST_0][MESH_0],
                                                                       solver_container[iZone][INST_0][MESH_0], config_container[iZone]);
        break;
    }

    direct_output[iZone]->SetHistory_Output(geometry_container[iZone][INST_0][MESH_0],
                                            solver_container[iZone][INST_0][MESH_0], config_container[iZone]);
  }

  /*--- Extract objective function values. ---*/

  for (iZone = 0; iZone < nZone; iZone++){

    nMarker_Analyze = config_container[iZone]->GetnMarker_Analyze();

    for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

      Weight_ObjFunc = config_container[iZone]->GetWeight_ObjFunc(iMarker_Analyze);

      switch (config_container[iZone]->GetKind_Solver()) {

        case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
          // per-surface output to be added soon
          break;
        case HEAT_EQUATION_FVM: case DISC_ADJ_HEAT:
          // per-surface output to be added soon
          break;
        default:
          break;
      }
    }

    /*--- Not-per-surface objective functions (shall not be included above) ---*/

    Weight_ObjFunc = config_container[iZone]->GetWeight_ObjFunc(0);

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:

        switch(config_container[iZone]->GetKind_ObjFunc()) {

          // Aerodynamic coefficients

          case DRAG_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("DRAG")*Weight_ObjFunc;
            break;
          case LIFT_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("LIFT")*Weight_ObjFunc;
            break;
          case SIDEFORCE_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("SIDEFORCE")*Weight_ObjFunc;
            break;
          case EFFICIENCY:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("EFFICIENCY")*Weight_ObjFunc;
            break;
          case MOMENT_X_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("MOMENT-X")*Weight_ObjFunc;
            break;
          case MOMENT_Y_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("MOMENT-Y")*Weight_ObjFunc;
            break;
          case MOMENT_Z_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("MOMENT-Z")*Weight_ObjFunc;
            break;
          case FORCE_X_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("FORCE-X")*Weight_ObjFunc;
            break;
          case FORCE_Y_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("FORCE-Y")*Weight_ObjFunc;
            break;
          case FORCE_Z_COEFFICIENT:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("FORCE-Z")*Weight_ObjFunc;
            break;

          // Other surface-related output values

          case SURFACE_MASSFLOW:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("AVG_MASSFLOW")*Weight_ObjFunc;
            break;
          case SURFACE_MACH:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("AVG_MACH")*Weight_ObjFunc;
            break;
          case SURFACE_UNIFORMITY:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("UNIFORMITY")*Weight_ObjFunc;
            break;
          case SURFACE_SECONDARY:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("SECONDARY_STRENGTH")*Weight_ObjFunc;
            break;
          case SURFACE_MOM_DISTORTION:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("MOMENTUM_DISTORTION")*Weight_ObjFunc;
            break;
          case SURFACE_SECOND_OVER_UNIFORM:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("SECONDARY_OVER_UNIFORMITY")*Weight_ObjFunc;
            break;
          case TOTAL_AVG_TEMPERATURE:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("AVG_TOTALTEMP")*Weight_ObjFunc;
            break;
          case SURFACE_TOTAL_PRESSURE:
            ObjFunc+=direct_output[iZone]->GetHistoryFieldValue("AVG_TOTALPRESS")*Weight_ObjFunc;
            break;

          // Not yet covered by new output structure. Be careful these use MARKER_MONITORING.

          case SURFACE_PRESSURE_DROP:
            ObjFunc+=config_container[iZone]->GetSurface_PressureDrop(0)*Weight_ObjFunc;
            break;
          case SURFACE_STATIC_PRESSURE:
            ObjFunc+=config_container[iZone]->GetSurface_Pressure(0)*Weight_ObjFunc;
            break;
          case TOTAL_HEATFLUX:
            ObjFunc += solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_HeatFlux()*Weight_ObjFunc;
            break;

          default:
            cout << "Objective function not covered for discrete adjoint multiphysics." << endl;
            break;
        }
        break;

        case DISC_ADJ_HEAT:

        switch(config_container[iZone]->GetKind_ObjFunc()) {

          // Not yet covered by new output structure. Be careful these use MARKER_MONITORING.

          case TOTAL_HEATFLUX:
            ObjFunc += solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetTotal_HeatFlux()*Weight_ObjFunc;
            break;
          case TOTAL_AVG_TEMPERATURE:
            ObjFunc += solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetTotal_AvgTemperature()*Weight_ObjFunc;
            break;

          default:
            cout << "Objective function not covered for discrete adjoint multiphysics." << endl;
            break;
        }
        break;

      default:
        break;
    }
  }

  if (rank == MASTER_NODE) {
    AD::RegisterOutput(ObjFunc);
    AD::SetIndex(ObjFunc_Index, ObjFunc);
    if (kind_recording == FLOW_CONS_VARS) {
      cout << " Objective function                   : " << ObjFunc;
      if (driver_config->GetWrt_AD_Statistics()){
        cout << " (" << ObjFunc_Index << ")" << endl;
      }
      cout << endl;
    }
  }
}

void CDiscAdjMultizoneDriver::SetAdj_ObjFunction() {

  bool TimeDomain = config_container[ZONE_0]->GetTime_Marching() != STEADY;
  unsigned long IterAvg_Obj = config_container[ZONE_0]->GetIter_Avg_Objective();

  su2double seeding = 1.0;

  if (TimeDomain){
    if (TimeIter < IterAvg_Obj){
      seeding = 1.0/((su2double)IterAvg_Obj);
    }
    else{
      seeding = 0.0;
    }
  }

  if (rank == MASTER_NODE) {
    AD::SetDerivative(ObjFunc_Index, SU2_TYPE::GetValue(seeding));
  }
}

void CDiscAdjMultizoneDriver::ComputeAdjoints(unsigned short iZone) {

  unsigned short enter_izone = iZone*2+1 + ITERATION_READY;
  unsigned short leave_izone = iZone*2 + ITERATION_READY;

  AD::ClearAdjoints();

  /*--- Initialize the adjoints in iZone ---*/

  iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container,
                                                        config_container, iZone, INST_0);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint(enter_izone, leave_izone);

  AD::ComputeAdjoint(TRANSFER, OBJECTIVE_FUNCTION);
  AD::ComputeAdjoint(DEPENDENCIES, START);
}

void CDiscAdjMultizoneDriver::Add_ExternalOld_To_Solution(unsigned short iZone) {

  unsigned short iSol;

  for (iSol=0; iSol < MAX_SOLS; iSol++){
    if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL) {
      if (solver_container[iZone][INST_0][MESH_0][iSol]->GetAdjoint()) {
        solver_container[iZone][INST_0][MESH_0][iSol]->Add_ExternalOld_To_Solution(geometry_container[iZone][INST_0][MESH_0]);
      }
    }
  }
}

void CDiscAdjMultizoneDriver::SetExternal_Zero(void) {

  unsigned short iZone, iSol;

  for(iZone=0; iZone < nZone; iZone++) {

    for (iSol=0; iSol < MAX_SOLS; iSol++){
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL) {
        if (solver_container[iZone][INST_0][MESH_0][iSol]->GetAdjoint()) {
          solver_container[iZone][INST_0][MESH_0][iSol]->GetNodes()->SetExternalZero();
        }
      }
    }
  }
}

void CDiscAdjMultizoneDriver::Set_OldExternal(void) {

  unsigned short iZone, iSol;

  for(iZone=0; iZone < nZone; iZone++) {

    for (iSol=0; iSol < MAX_SOLS; iSol++){
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL) {
        if (solver_container[iZone][INST_0][MESH_0][iSol]->GetAdjoint()) {
          solver_container[iZone][INST_0][MESH_0][iSol]->GetNodes()->Set_OldExternal();
        }
      }
    }
  }
}

void CDiscAdjMultizoneDriver::Add_Solution_To_External(unsigned short iZone) {

  unsigned short iSol;

  for (iSol=0; iSol < MAX_SOLS; iSol++){
    if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL) {
      if (solver_container[iZone][INST_0][MESH_0][iSol]->GetAdjoint()) {
        solver_container[iZone][INST_0][MESH_0][iSol]->Add_Solution_To_External(geometry_container[iZone][INST_0][MESH_0]);
      }
    }
  }
}

void CDiscAdjMultizoneDriver::Add_Solution_To_ExternalOld(unsigned short iZone) {

  unsigned short iSol;

  for (iSol=0; iSol < MAX_SOLS; iSol++){
    if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL) {
      if (solver_container[iZone][INST_0][MESH_0][iSol]->GetAdjoint()) {
        solver_container[iZone][INST_0][MESH_0][iSol]->Add_Solution_To_ExternalOld(geometry_container[iZone][INST_0][MESH_0]);
      }
    }
  }
}

void CDiscAdjMultizoneDriver::Set_BGSSolution(unsigned short iZone) {

  unsigned short iSol;

  for (iSol=0; iSol < MAX_SOLS; iSol++){
    if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL) {
      if (solver_container[iZone][INST_0][MESH_0][iSol]->GetAdjoint()) {
        solver_container[iZone][INST_0][MESH_0][iSol]->UpdateSolution_BGS(geometry_container[iZone][INST_0][MESH_0],
                                                                          config_container[iZone]);
      }
    }
  }
}

void CDiscAdjMultizoneDriver::SetResidual_BGS(unsigned short iZone) {

  unsigned short iSol;

  for (iSol=0; iSol < MAX_SOLS; iSol++){
    if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL) {
      if (solver_container[iZone][INST_0][MESH_0][iSol]->GetAdjoint()) {
        solver_container[iZone][INST_0][MESH_0][iSol]->ComputeResidual_Multizone(geometry_container[iZone][INST_0][MESH_0],
                                                                                 config_container[iZone]);
      }
    }
  }
}
