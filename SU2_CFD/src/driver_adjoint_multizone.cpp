/*!
 * \file driver_adjoint_multizone.cpp
 * \brief The main subroutines for driving multi-zone problems.
 * \author O. Burghardt, T. Albring, R. Sanchez
 * \version 6.0.1 "Falcon"
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

#include "../include/driver_structure.hpp"
#include "../include/definition_structure.hpp"

CDiscAdjMultizoneDriver::CDiscAdjMultizoneDriver(char* confFile,
                                      unsigned short val_nZone,
                                      unsigned short val_nDim,
                                      bool val_periodic,
                                      SU2_Comm MPICommunicator) : CMultizoneDriver(confFile,
                                                                                  val_nZone,
                                                                                  val_nDim,
                                                                                  val_periodic,
                                                                                  MPICommunicator) {

  RecordingState = NONE;

  nInst = new unsigned short[nZone];

  direct_iteration = new CIteration**[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {

    nInst[iZone]            = 1;
    direct_iteration[iZone] = new CIteration*[nInst[iZone]];

    for(iInst = 0; iInst < nInst[iZone]; iInst++) {

      switch (config_container[iZone]->GetKind_Solver()) {

        case EULER: case NAVIER_STOKES: case RANS:
        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          direct_iteration[iZone][iInst] = new CFluidIteration(config_container[iZone]);
          break;
        case HEAT_EQUATION_FVM: case DISC_ADJ_HEAT:
          direct_iteration[iZone][iInst] = new CHeatIteration(config_container[iZone]);
          break;

        default:
          SU2_MPI::Error("There is no discrete adjoint functionality for one of the specified solvers yet.", CURRENT_FUNCTION);
      }
    }
  }
}

CDiscAdjMultizoneDriver::~CDiscAdjMultizoneDriver(){

  for (iZone = 0; iZone < nZone; iZone++){
    for (iInst = 0; iInst < nInst[iInst]; iInst++){
      delete direct_iteration[iZone][iInst];
    }
    delete direct_iteration[iZone];
  }

  delete [] direct_iteration;
}

void CDiscAdjMultizoneDriver::StartSolver() {

  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;

  if (rank == MASTER_NODE){
    cout << endl <<"Simulation Run using the Discrete Adjoint Multizone Driver" << endl;
    if (driver_config->GetTime_Domain())
      SU2_MPI::Error("The discrete adjoint multizone driver is not ready for unsteady computations yet.", CURRENT_FUNCTION);
  }

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the value of the external iteration to TimeIter. -------------------------------------*/
    /*--- TODO: This should be generalised for an homogeneous criteria throughout the code. --------*/
    config_container[iZone]->SetTimeIter(0);
    
  }
  
  /*--- We directly start the (steady-state) discrete adjoint computation. ---*/
  Run();
  
  /*--- Output the solution in files. ---*/

  Output(0);

}

void CDiscAdjMultizoneDriver::Output(unsigned long TimeIter) {

  
  for (iZone = 0; iZone < nZone; iZone++){
    
    bool output_files = false;
    
    /*--- Determine whether a solution needs to be written
   after the current iteration ---*/
    
    if (
        
        /*--- General if statements to print output statements ---*/
        
        (TimeIter+1 >= config_container[iZone]->GetnTime_Iter()) || (StopCalc) ||
        
        /*--- Unsteady problems ---*/
        
        (((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
          (config_container[iZone]->GetUnsteady_Simulation() == TIME_STEPPING)) &&
         ((TimeIter == 0) || (ExtIter % config_container[iZone]->GetWrt_Sol_Freq_DualTime() == 0))) ||
        
        ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND) &&
         ((TimeIter == 0) || ((TimeIter % config_container[iZone]->GetWrt_Sol_Freq_DualTime() == 0) ||
                              ((TimeIter-1) % config_container[iZone]->GetWrt_Sol_Freq_DualTime() == 0)))) ||
        
        ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND) &&
         ((TimeIter == 0) || ((TimeIter % config_container[iZone]->GetWrt_Sol_Freq_DualTime() == 0)))) ||
        
        ((config_container[iZone]->GetDynamic_Analysis() == DYNAMIC) &&
         ((TimeIter == 0) || (TimeIter % config_container[iZone]->GetWrt_Sol_Freq_DualTime() == 0))) ||
        
        /*--- No inlet profile file found. Print template. ---*/
        
        (config_container[iZone]->GetWrt_InletFile())
        
        ) {
      
      output_files = true;
      
    }
    
    /*--- Determine whether a solution doesn't need to be written
   after the current iteration ---*/
    
    if (config_container[iZone]->GetFixed_CL_Mode()) {
      if (config_container[iZone]->GetnExtIter()-config_container[iZone]->GetIter_dCL_dAlpha() - 1 < ExtIter) output_files = false;
      if (config_container[iZone]->GetnExtIter() - 1 == ExtIter) output_files = true;
    }
    
    /*--- write the solution ---*/
    
    if (output_files) {
      
      /*--- Time the output for performance benchmarking. ---*/
#ifndef HAVE_MPI
      StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
      StopTime = MPI_Wtime();
#endif
      UsedTimeCompute += StopTime-StartTime;
#ifndef HAVE_MPI
      StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
      StartTime = MPI_Wtime();
#endif
      
      if (rank == MASTER_NODE) cout << endl << "-------------------------- File Output Summary --------------------------";
      
      /*--- Execute the routine for writing restart, volume solution,
     surface solution, and surface comma-separated value files. ---*/
      
      output[iZone]->SetResult_Files_Parallel(solver_container, geometry_container, config_container, TimeIter, iZone, nZone);
      
      
      /*--- Execute the routine for writing special output. ---*/
      // output->SetSpecial_Output(solver_container, geometry_container, config_container, TimeIter, nZone);
      
      
      if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl << endl;
      
    }
    /*--- Store output time and restart the timer for the compute phase. ---*/
#ifndef HAVE_MPI
    StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
    StopTime = MPI_Wtime();
#endif
    UsedTimeOutput += StopTime-StartTime;
    OutputCount++;
    BandwidthSum = config_container[ZONE_0]->GetRestart_Bandwidth_Agg();
#ifndef HAVE_MPI
    StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
    StartTime = MPI_Wtime();
#endif

  }

}


void CDiscAdjMultizoneDriver::Run() {

  unsigned short  checkConvergence  = 0,
                  jZone             = 0;
  unsigned long   iOuter_Iter       = 0,
                  OuterIter         = driver_config->GetnOuter_Iter();

  /*--- Loop over the number of outer iterations ---*/

  for (iOuter_Iter = 0; iOuter_Iter < OuterIter; iOuter_Iter++){

    for (iZone = 0; iZone < nZone; iZone++) {

      config_container[iZone]->SetOuterIter(iOuter_Iter);

      iteration_container[iZone][INST_0]->Preprocess(output[iZone], integration_container, geometry_container,
                                                       solver_container, numerics_container, config_container,
                                                       surface_movement, grid_movement, FFDBox, iZone, INST_0);
    }

    /*--- For the adjoint iteration we need the derivatives of the iteration function with
     *    respect to the state (and possibly the mesh coordinate) variables.
     *    Since these derivatives do not change in the steady state case we only have to record
     *    if the current recording is different from them.
     *
     *    To set the tape appropriately, the following recording methods are provided:
     *    (1) NONE: All information from a previous recording is removed.
     *    (2) STATE_VARS: Store computational graph of one direct iteration with state variables
     *        (e.g. conservatives variables for a flow solver) as input.
     *    (3) MESH_COORDS: Mesh coordinates as input.
     *    (4) COMBINED: Mesh coordinates and state variables as input.
     *
     *    By default, all (state and mesh coordinate variables) will be declared as output,
     *    since it does not change the computational effort. ---*/

    if (RecordingState != FLOW_CONS_VARS) {

      SetRecording(NONE);

      SetRecording(FLOW_CONS_VARS);

    }

    Set_OldAdjoints(); SetIter_Zero();

    AD::ClearAdjoints();

    SetAdj_ObjFunction();

    AD::ComputeAdjoint(3,0);

    for (iZone = 0; iZone < nZone; iZone++) {
      
      iteration_container[iZone][INST_0]->Iterate(output[iZone], integration_container, geometry_container,
                                            solver_container, numerics_container, config_container,
                                            surface_movement, grid_movement, FFDBox, iZone, INST_0);
    }

    Add_IterAdjoints();

    for (iZone = 0; iZone < nZone; iZone++) {
      
      config_container[iZone]->SetInnerIter(0);

      ComputeZonewiseAdjoints(iZone);

      for (jZone = 0; jZone < nZone; jZone++) {

          iteration_container[jZone][INST_0]->Iterate(output[iZone], integration_container, geometry_container,
                                            solver_container, numerics_container, config_container,
                                            surface_movement, grid_movement, FFDBox, jZone, INST_0);
      }

      Add_IterAdjoints();
    }

    SetAdjoints_Iter();

    SetResidual_RMS();

    /*--- Check convergence in each zone --*/

    checkConvergence = 0;
    for (iZone = 0; iZone < nZone; iZone++) {

      switch (config_container[iZone]->GetKind_Solver()) {

        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          checkConvergence += (int) integration_container[iZone][INST_0][ADJFLOW_SOL]->GetConvergence();
          break;
        case DISC_ADJ_HEAT:
          checkConvergence += (int) integration_container[iZone][INST_0][ADJHEAT_SOL]->GetConvergence();
          break;
        default:
          checkConvergence += 1;
          break;
      }
    }

    /*--- If convergence was reached in every zone --*/

    if (checkConvergence == nZone) break;

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();


    /*--- Compute the geometrical sensitivities ---*/

    checkConvergence = 0;
    switch (config_container[ZONE_0]->GetKind_Solver()) {

      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
        checkConvergence += (int) integration_container[ZONE_0][INST_0][ADJFLOW_SOL]->GetConvergence();
        break;
      case DISC_ADJ_HEAT:
        checkConvergence += (int) integration_container[ZONE_0][INST_0][ADJHEAT_SOL]->GetConvergence();
        break;
      default:
        checkConvergence += 1;
        break;
    }

    if ((config_container[ZONE_0]->GetOuterIter()+1 >= OuterIter) || checkConvergence == 1 ||
        (config_container[ZONE_0]->GetOuterIter() % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0)){

      /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
       * as argument ensures that all information from a previous recording is removed. ---*/

      SetRecording(NONE);

      /*--- Store the computational graph of one direct iteration with the mesh coordinates as input. ---*/

      SetRecording(MESH_COORDS);

      /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
       *    of the current iteration. The values are passed to the AD tool. ---*/

      Set_OldAdjoints();

      for (iZone = 0; iZone < nZone; iZone++) {

        iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container, config_container, iZone, INST_0);
      }

      /*--- Initialize the adjoint of the objective function with 1.0. ---*/

      SetAdj_ObjFunction();

      /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

      AD::ComputeAdjoint();

      /*--- Extract the computed sensitivity values. ---*/

      for (iZone = 0; iZone < nZone; iZone++) {

        switch (config_container[iZone]->GetKind_Solver()) {

          case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
            solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);
            break;
          case DISC_ADJ_HEAT:
            solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);
            break;
          default:
            cout << "WARNING: Setting sensitivities failed for one of the specified discrete adjoint solvers!" << endl;
            break;
        }
      }

      /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

      AD::ClearAdjoints();
    }
  }
}

void CDiscAdjMultizoneDriver::SetRecording(unsigned short kind_recording) {

  unsigned short jZone, iMesh, UpdateMesh;
  unsigned long ExtIter = 0;
  bool DeformMesh       = false;
  bool heat             = false;

  AD::Reset();

  /*--- Prepare for recording by resetting the flow solution to the initial converged solution---*/

  for(iZone = 0; iZone < nZone; iZone++) {

    heat = config_container[iZone]->GetWeakly_Coupled_Heat();

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES:
        for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
          solver_container[iZone][INST_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[iZone][INST_0][iMesh], config_container[iZone]);
        }
        if (heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetRecording(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        }
        break;
      case DISC_ADJ_RANS:
        for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
          solver_container[iZone][INST_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[iZone][INST_0][iMesh], config_container[iZone]);
        }
        if (!config_container[iZone]->GetFrozen_Visc_Disc()) {
          solver_container[iZone][INST_0][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        }
        if (heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetRecording(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        }
        break;
      case DISC_ADJ_HEAT:
        solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetRecording(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);
        break;
    }
  }

  /*---Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if(kind_recording != NONE) {

    AD::StartRecording();

    AD::Push_TapePosition();

    if (rank == MASTER_NODE && kind_recording == FLOW_CONS_VARS) {
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << "Direct iteration to store computational graph." << endl;
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
      cout << "-------------------------------------------------------------------------" << endl << endl;
    }

    for (iZone = 0; iZone < nZone; iZone++) {

      iteration_container[iZone][INST_0]->RegisterInput(solver_container, geometry_container, config_container, iZone, INST_0, kind_recording);
    }
  }

  AD::Push_TapePosition();

  for(iZone = 0; iZone < nZone; iZone++) {

    iteration_container[iZone][INST_0]->SetDependencies(solver_container, geometry_container, config_container, iZone, INST_0, kind_recording);
  }

  AD::Push_TapePosition();

  /*--- Extract the objective function and store it --- */

  SetObjFunction();

  AD::Push_TapePosition();

  /*--- We do the communication here to not derive wrt updated boundary data. ---*/

  for(iZone = 0; iZone < nZone; iZone++) {

    /*--- In principle, the mesh does not need to be updated ---*/
    UpdateMesh = 0;

    /*--- Transfer from all the remaining zones ---*/
    for (jZone = 0; jZone < nZone; jZone++){
      /*--- The target zone is iZone ---*/
      if (jZone != iZone && transfer_container[iZone][jZone] != NULL){
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

    /*--- Do one iteration of the direct solver ---*/
    direct_iteration[iZone][INST_0]->Preprocess(output[iZone], integration_container, geometry_container, solver_container,
        numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

    /*--- Iterate the zone as a block a single time ---*/
    direct_iteration[iZone][INST_0]->Iterate(output[iZone], integration_container, geometry_container, solver_container,
        numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

//    Corrector(iZone);

    /*--- Print residuals in the first iteration ---*/

    if (rank == MASTER_NODE && kind_recording == FLOW_CONS_VARS) {

      switch (config_container[iZone]->GetKind_Solver()) {

        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES:
          cout << " Zone " << iZone << ": log10[Conservative 0]: " << log10(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetRes_RMS(0)) << endl;
          break;
        case DISC_ADJ_RANS:
          cout << " Zone " << iZone << ": log10[Conservative 0]: " << log10(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetRes_RMS(0)) << endl;
          if (!config_container[iZone]->GetFrozen_Visc_Disc()) {
            cout <<"       log10[RMS k]: " << log10(solver_container[iZone][INST_0][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
          }
          break;
        case DISC_ADJ_HEAT:
          cout << " Zone " << iZone << ": log10[Conservative 0]: " << log10(solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetRes_RMS(0)) << endl;
          break;
      }
    }

    iteration_container[iZone][INST_0]->RegisterOutput(solver_container, geometry_container, config_container, output, iZone, INST_0);

    AD::Push_TapePosition();
  }

  AD::StopRecording();

  AD::PrintStatistics();

  RecordingState = kind_recording;
}

void CDiscAdjMultizoneDriver::SetObjFunction() {

  bool heat = false;

  ObjFunc = 0.0;

  /*--- Repeat objective function calculations, in case they have been included in
   *    the integration step and therefore have to be excluded from this part of the tape. ---*/

  for (iZone = 0; iZone < nZone; iZone++){

    heat = config_container[iZone]->GetWeakly_Coupled_Heat();

    switch (config_container[iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
        solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Pressure_Forces(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Momentum_Forces(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Friction_Forces(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        if(heat) {
          solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->Heat_Fluxes(geometry_container[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0], config_container[iZone]);
        }
        break;
      case HEAT_EQUATION_FVM: case DISC_ADJ_HEAT:
        solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->Heat_Fluxes(geometry_container[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0], config_container[iZone]);
        break;
    }
  }

  for (iZone = 0; iZone < nZone; iZone++){

    switch (config_container[iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
      solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->SetTotal_ComboObj(0.0);
      break;
    }
  }

  /*--- Specific scalar objective functions ---*/

  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:

        if (config_container[ZONE_0]->GetnMarker_Analyze() != 0)
          output->SpecialOutput_AnalyzeSurface(solver_container[iZone][INST_0][MESH_0][FLOW_SOL], geometry_container[iZone][INST_0][MESH_0], config_container[iZone], false);

        if (config_container[ZONE_0]->GetnMarker_Analyze() != 0)
          output->SpecialOutput_Distortion(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], false);

        if (config_container[ZONE_0]->GetnMarker_NearFieldBound() != 0)
          output->SpecialOutput_SonicBoom(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], false);

        if (config_container[ZONE_0]->GetPlot_Section_Forces())
          output->SpecialOutput_SpanLoad(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], false);

        break;
    }
  }

  /*--- Surface based obj. function ---*/

  for (iZone = 0; iZone < nZone; iZone++){

    heat = config_container[iZone]->GetWeakly_Coupled_Heat();

    switch (config_container[iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
        solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Evaluate_ObjFunc(config_container[iZone]);
        ObjFunc += solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_ComboObj();

        if(heat) {
          if (config_container[iZone]->GetKind_ObjFunc() == TOTAL_HEATFLUX) {
            ObjFunc += solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetTotal_HeatFlux();
          }
          else if (config_container[iZone]->GetKind_ObjFunc() == TOTAL_AVG_TEMPERATURE) {
            ObjFunc += solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetTotal_AvgTemperature();
          }
        }
        break;

      case HEAT_EQUATION_FVM: case DISC_ADJ_HEAT:
        if (config_container[iZone]->GetKind_ObjFunc() == TOTAL_HEATFLUX) {
          ObjFunc += solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetTotal_HeatFlux();
        }
        else if (config_container[iZone]->GetKind_ObjFunc() == TOTAL_AVG_TEMPERATURE) {
          ObjFunc += solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->GetTotal_AvgTemperature();
        }
        break;
    }
  }

  if (rank == MASTER_NODE){
    AD::RegisterOutput(ObjFunc);
    AD::Set_AdjIndex(ObjFunc_Index, ObjFunc);
    cout << "Setting objective function value (in recording): " << ObjFunc << " (" << ObjFunc_Index << ")" << endl;
  }
}

void CDiscAdjMultizoneDriver::SetAdj_ObjFunction() {

  bool time_stepping = config_container[ZONE_0]->GetUnsteady_Simulation() != STEADY;
  unsigned long IterAvg_Obj = config_container[ZONE_0]->GetIter_Avg_Objective();
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  su2double seeding = 1.0;

  if (time_stepping){
    if (ExtIter < IterAvg_Obj){
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

void CDiscAdjMultizoneDriver::ComputeZonewiseAdjoints(unsigned short iZone) {

  /*--- Position markers
   * 0: Recording started
   * 1: Copied solution registered and direct solver resetted to copied solution
   * 2: Dependencies are set
   * 3: Objective function is set
   * 4: Data transferred between zones
   ---*/

  unsigned short leave_izone = iZone*2+5;
  unsigned short enter_izone = iZone*2+6;

  AD::ClearAdjoints();

  /*--- Initialize the adjoints in iZone ---*/

  iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container, config_container, iZone, INST_0);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint(enter_izone, leave_izone);
  AD::ComputeAdjoint(4,3);
  AD::ComputeAdjoint(2,0);
}

void CDiscAdjMultizoneDriver::Set_OldAdjoints(void) {

  unsigned short iZone;
  bool weakly_coupled_heat;

  for(iZone=0; iZone < nZone; iZone++) {

    weakly_coupled_heat = config_container[iZone]->GetWeakly_Coupled_Heat();

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Set_OldSolution(geometry_container[iZone][INST_0][MESH_0]);
        break;
      case DISC_ADJ_NAVIER_STOKES:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Set_OldSolution(geometry_container[iZone][INST_0][MESH_0]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Set_OldSolution(geometry_container[iZone][INST_0][MESH_0]);
        }
        break;
      case DISC_ADJ_RANS:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Set_OldSolution(geometry_container[iZone][INST_0][MESH_0]);
        solver_container[iZone][INST_0][MESH_0][ADJTURB_SOL]->Set_OldSolution(geometry_container[iZone][INST_0][MESH_0]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Set_OldSolution(geometry_container[iZone][INST_0][MESH_0]);
        }
        break;
      case DISC_ADJ_HEAT:
        solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Set_OldSolution(geometry_container[iZone][INST_0][MESH_0]);
        break;
    }
  }
}

void CDiscAdjMultizoneDriver::SetIter_Zero(void) {

  unsigned short iZone;
  bool weakly_coupled_heat;

  for(iZone=0; iZone < nZone; iZone++) {

    weakly_coupled_heat = config_container[iZone]->GetWeakly_Coupled_Heat();

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Set_IterSolution_Zero(geometry_container[iZone][INST_0][MESH_0]);
        break;
      case DISC_ADJ_NAVIER_STOKES:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Set_IterSolution_Zero(geometry_container[iZone][INST_0][MESH_0]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Set_IterSolution_Zero(geometry_container[iZone][INST_0][MESH_0]);
        }
        break;
      case DISC_ADJ_RANS:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Set_IterSolution_Zero(geometry_container[iZone][INST_0][MESH_0]);
        solver_container[iZone][INST_0][MESH_0][ADJTURB_SOL]->Set_IterSolution_Zero(geometry_container[iZone][INST_0][MESH_0]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Set_IterSolution_Zero(geometry_container[iZone][INST_0][MESH_0]);
        }
        break;
      case DISC_ADJ_HEAT:
        solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Set_IterSolution_Zero(geometry_container[iZone][INST_0][MESH_0]);
        break;
    }
  }
}

void CDiscAdjMultizoneDriver::Add_IterAdjoints(void) {

  unsigned short iZone;
  bool weakly_coupled_heat;

  for(iZone=0; iZone < nZone; iZone++) {

    weakly_coupled_heat = config_container[iZone]->GetWeakly_Coupled_Heat();

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Add_IterSolution(geometry_container[iZone][INST_0][MESH_0]);
        break;
      case DISC_ADJ_NAVIER_STOKES:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Add_IterSolution(geometry_container[iZone][INST_0][MESH_0]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Add_IterSolution(geometry_container[iZone][INST_0][MESH_0]);
        }
        break;
      case DISC_ADJ_RANS:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->Add_IterSolution(geometry_container[iZone][INST_0][MESH_0]);
        solver_container[iZone][INST_0][MESH_0][ADJTURB_SOL]->Add_IterSolution(geometry_container[iZone][INST_0][MESH_0]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Add_IterSolution(geometry_container[iZone][INST_0][MESH_0]);
        }
        break;
      case DISC_ADJ_HEAT:
        solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->Add_IterSolution(geometry_container[iZone][INST_0][MESH_0]);
        break;
    }
  }
}

void CDiscAdjMultizoneDriver::SetAdjoints_Iter(void) {

  unsigned short iZone;
  bool weakly_coupled_heat;

  for(iZone=0; iZone < nZone; iZone++) {

    weakly_coupled_heat = config_container[iZone]->GetWeakly_Coupled_Heat();

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSolution_Iter(geometry_container[iZone][INST_0][MESH_0]);
        break;
      case DISC_ADJ_NAVIER_STOKES:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSolution_Iter(geometry_container[iZone][INST_0][MESH_0]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetSolution_Iter(geometry_container[iZone][INST_0][MESH_0]);
        }
        break;
      case DISC_ADJ_RANS:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSolution_Iter(geometry_container[iZone][INST_0][MESH_0]);
        solver_container[iZone][INST_0][MESH_0][ADJTURB_SOL]->SetSolution_Iter(geometry_container[iZone][INST_0][MESH_0]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetSolution_Iter(geometry_container[iZone][INST_0][MESH_0]);
        }
        break;
      case DISC_ADJ_HEAT:
        solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetSolution_Iter(geometry_container[iZone][INST_0][MESH_0]);
        break;
    }
  }
}

void CDiscAdjMultizoneDriver::SetResidual_RMS(void) {

  unsigned short iZone;
  bool weakly_coupled_heat;

  for(iZone=0; iZone < nZone; iZone++) {

    weakly_coupled_heat = config_container[iZone]->GetWeakly_Coupled_Heat();

    switch (config_container[iZone]->GetKind_Solver()) {

      case DISC_ADJ_EULER:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetResidual_Solution(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        break;
      case DISC_ADJ_NAVIER_STOKES:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetResidual_Solution(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetResidual_Solution(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        }
        break;
      case DISC_ADJ_RANS:
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetResidual_Solution(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        solver_container[iZone][INST_0][MESH_0][ADJTURB_SOL]->SetResidual_Solution(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        if (weakly_coupled_heat) {
          solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetResidual_Solution(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        }
        break;
      case DISC_ADJ_HEAT:
        solver_container[iZone][INST_0][MESH_0][ADJHEAT_SOL]->SetResidual_Solution(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        break;
    }
  }
}
