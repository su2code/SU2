/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving multi-zone problems.
 * \author R. Sanchez, O. Burghardt
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

#include "../../include/drivers/CMultizoneDriver.hpp"
#include "../../include/definition_structure.hpp"
#include "../../../Common/include/interface_interpolation/CInterpolator.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIteration.hpp"

CMultizoneDriver::CMultizoneDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator) :
                  CDriver(confFile, val_nZone, MPICommunicator, false) {

  /*--- Initialize the counter for TimeIter ---*/
  TimeIter = 0;

  /*--- Structure for multizone convergence ---*/
  init_res     = new su2double*[nZone];
  residual     = new su2double*[nZone];
  residual_rel = new su2double*[nZone];
  nVarZone     = new unsigned short[nZone];

  for (iZone = 0; iZone < nZone; iZone++){
    nVarZone[iZone] = 0;
    /*--- Account for all the solvers ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      if (solver_container[iZone][INST_0][MESH_0][iSol] != nullptr)
        nVarZone[iZone] += solver_container[iZone][INST_0][MESH_0][iSol]->GetnVar();
    }
    init_res[iZone]     = new su2double[nVarZone[iZone]];
    residual[iZone]     = new su2double[nVarZone[iZone]];
    residual_rel[iZone] = new su2double[nVarZone[iZone]];
    /*--- Initialize the residual containers to 0.0 ---*/
    for (unsigned short iVar = 0; iVar < nVarZone[iZone]; iVar++){
      init_res[iZone][iVar]     = 0.0;
      residual[iZone][iVar]     = 0.0;
      residual_rel[iZone][iVar] = 0.0;
    }
  }

  /*----------------------------------------------------*/
  /*------ Determine the properties of the problem -----*/
  /*----------------------------------------------------*/

  bool structural_zone = false;
  bool fluid_zone = false;
  bool heat_zone = false;

  /*--- If there is at least a fluid and a structural zone ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_Solver()) {
    case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
    case MAIN_SOLVER::INC_EULER: case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
      fluid_zone = true;
      break;
    case MAIN_SOLVER::NEMO_EULER: case MAIN_SOLVER::NEMO_NAVIER_STOKES:
      fluid_zone = true;
      break;
    case MAIN_SOLVER::FEM_ELASTICITY:
      structural_zone = true;
      break;
    case MAIN_SOLVER::HEAT_EQUATION:
      heat_zone = true;
      break;
    default:
      break;
    }
  }

  fsi = false; cht = false;
  /*--- If the problem has FSI properties ---*/
  if (fluid_zone && structural_zone) fsi = true;
  /*--- If the problem has CHT properties ---*/
  if (fluid_zone && heat_zone) cht = true;

  /*----------------------------------------------------*/
  /*- Define if a prefixed motion is imposed in a zone -*/
  /*----------------------------------------------------*/

  prefixed_motion = new bool[nZone];
  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_GridMovement()){
      case RIGID_MOTION:
        prefixed_motion[iZone] = true; break;
      default:
        prefixed_motion[iZone] = false; break;
    }
    if (config_container[iZone]->GetSurface_Movement(AEROELASTIC) ||
        config_container[iZone]->GetSurface_Movement(AEROELASTIC_RIGID_MOTION) ||
        config_container[iZone]->GetSurface_Movement(DEFORMING) ||
        config_container[iZone]->GetSurface_Movement(EXTERNAL) ||
        config_container[iZone]->GetSurface_Movement(EXTERNAL_ROTATION)){
      prefixed_motion[iZone] = true;
    }
  }

}

CMultizoneDriver::~CMultizoneDriver(void) {

  for (iZone = 0; iZone < nZone; iZone++){
    delete [] init_res[iZone];
    delete [] residual[iZone];
    delete [] residual_rel[iZone];
  }

  delete [] nVarZone;
  delete [] init_res;
  delete [] residual;
  delete [] residual_rel;

  delete [] prefixed_motion;

}

void CMultizoneDriver::StartSolver() {

  /*--- Find out the minimum of all references times and then set each zone to this (same) value.
        To ensure that all zones run synchronously in time, be it a dimensional or non-dimensionalized one. ---*/

  su2double Time_Ref = std::numeric_limits<su2double>::max();
  for (iZone = 0; iZone < nZone; iZone++)
    Time_Ref = min(Time_Ref, config_container[iZone]->GetTime_Ref());

  for (iZone = 0; iZone < nZone; iZone++) {

    config_container[iZone]->SetTime_Ref(Time_Ref);

    /*--- Recompute some values as the reference time might has changed in iZone ---*/

    config_container[iZone]->SetDelta_UnstTimeND(config_container[iZone]->GetDelta_UnstTime() / Time_Ref);
    config_container[iZone]->SetTotal_UnstTimeND(config_container[iZone]->GetTotal_UnstTime() / Time_Ref);
  }

  StartTime = SU2_MPI::Wtime();

  driver_config->Set_StartTime(StartTime);

  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/

  if (rank == MASTER_NODE){
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;
  }

  if (rank == MASTER_NODE){
    cout << endl <<"Simulation Run using the Multizone Driver" << endl;
    if (driver_config->GetTime_Domain())
      cout << "The simulation will run until time step " << driver_config->GetnTime_Iter() - driver_config->GetRestart_Iter() << "." << endl;
  }

  /*--- Set the initial time iteration to the restart iteration. ---*/
  if (driver_config->GetRestart() && driver_config->GetTime_Domain())
    TimeIter = driver_config->GetRestart_Iter();

  /*--- Run the problem until the number of time iterations required is reached. ---*/
  while ( TimeIter < driver_config->GetnTime_Iter() ) {

    /*--- Perform some preprocessing before starting the time-step simulation. ---*/

    Preprocess(TimeIter);

    /*--- Run a block iteration of the multizone problem. ---*/

    switch (driver_config->GetKind_MZSolver()){
      case ENUM_MULTIZONE::MZ_BLOCK_GAUSS_SEIDEL: Run_GaussSeidel(); break;  // Block Gauss-Seidel iteration
      case ENUM_MULTIZONE::MZ_BLOCK_JACOBI: Run_Jacobi(); break;             // Block-Jacobi iteration
    }

    /*--- Update the solution for dual time stepping strategy ---*/

    Update();

    /*--- Monitor the computations after each iteration. ---*/

    StopCalc = Monitor(TimeIter);

    /*--- Output the solution in files. ---*/

    Output(TimeIter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    TimeIter++;

  }

}

void CMultizoneDriver::Preprocess(unsigned long TimeIter) {

  /*--- Set the current time iteration in the config ---*/
  driver_config->SetTimeIter(TimeIter);

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the value of the external iteration to TimeIter. -------------------------------------*/
    /*--- TODO: This should be generalised for an homogeneous criteria throughout the code. --------*/
    config_container[iZone]->SetTimeIter(TimeIter);


    /*--- Store the current physical time in the config container, as
     this can be used for verification / MMS. This should also be more
     general once the drivers are more stable. ---*/

    if (driver_config->GetTime_Domain()) {
      config_container[iZone]->SetPhysicalTime(static_cast<su2double>(TimeIter)*config_container[iZone]->GetDelta_UnstTimeND());
    }
    else {
      config_container[iZone]->SetPhysicalTime(0.0);
    }

    /*--- Set the initial condition for EULER/N-S/RANS ---------------------------------------------*/
    /*--- For FSI, the initial conditions are set, after the mesh has been moved. --------------------------------------*/
    if (!fsi && config_container[iZone]->GetFluidProblem()) {
      solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone][INST_0],
                                                                             solver_container[iZone][INST_0],
                                                                             config_container[iZone], TimeIter);
    }
    else if (!fsi && config_container[iZone]->GetHeatProblem()) {
      /*--- Set the initial condition for HEAT equation ---------------------------------------------*/
      solver_container[iZone][INST_0][MESH_0][HEAT_SOL]->SetInitialCondition(geometry_container[iZone][INST_0],
                                                                              solver_container[iZone][INST_0],
                                                                              config_container[iZone], TimeIter);
    }
  }

  SU2_MPI::Barrier(SU2_MPI::GetComm());

  /*--- Run a predictor step ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config_container[iZone]->GetPredictor())
      iteration_container[iZone][INST_0]->Predictor(output_container[iZone], integration_container, geometry_container,
                                                    solver_container, numerics_container, config_container, surface_movement,
                                                    grid_movement, FFDBox, iZone, INST_0);
  }

  /*--- Perform a dynamic mesh update if required. ---*/

  DynamicMeshUpdate(TimeIter);

  /*--- Updating zone interface communication patterns for unsteady problems with pre-fixed motion in the config file ---*/
  if (driver_config->GetTime_Domain()) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (unsigned short jZone = 0; jZone < nZone; jZone++){
        if(jZone != iZone && interpolator_container[iZone][jZone] != nullptr && prefixed_motion[iZone])
          interpolator_container[iZone][jZone]->SetTransferCoeff(config_container);
      }
    }
  }

}

void CMultizoneDriver::Run_GaussSeidel() {

  unsigned short UpdateMesh;
  bool DeformMesh = false;

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]->SetOuterIter(0ul);
  }

  /*--- Loop over the number of outer iterations ---*/
  for (auto iOuter_Iter = 0ul; iOuter_Iter < driver_config->GetnOuter_Iter(); iOuter_Iter++){

    /*--- Loop over the number of zones (IZONE) ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- In principle, the mesh does not need to be updated ---*/
      UpdateMesh = 0;

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);
      config_container[iZone]->Set_StartTime(SU2_MPI::Wtime());
      driver_config->SetOuterIter(iOuter_Iter);

      /*--- Transfer from all the remaining zones ---*/
      for (auto jZone = 0u; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone){
          DeformMesh = Transfer_Data(jZone, iZone);
          if (DeformMesh) UpdateMesh+=1;
        }
      }
      /*--- If a mesh update is required due to the transfer of data ---*/
      if (UpdateMesh > 0) DynamicMeshUpdate(iZone, TimeIter);

      /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
      iteration_container[iZone][INST_0]->Solve(output_container[iZone], integration_container, geometry_container,
                                                solver_container, numerics_container, config_container,
                                                surface_movement, grid_movement, FFDBox, iZone, INST_0);

      /*--- A corrector step can help preventing numerical instabilities ---*/
      Corrector(iZone);

    }

    if (OuterConvergence(iOuter_Iter)) break;

  }

}

void CMultizoneDriver::Run_Jacobi() {

  unsigned short UpdateMesh;
  bool DeformMesh = false;

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]->SetOuterIter(0ul);
  }

  /*--- Loop over the number of outer iterations ---*/
  for (auto iOuter_Iter = 0ul; iOuter_Iter < driver_config->GetnOuter_Iter(); iOuter_Iter++){

    /*--- Transfer from all zones ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- In principle, the mesh does not need to be updated ---*/
      UpdateMesh = 0;

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);
      driver_config->SetOuterIter(iOuter_Iter);

      /*--- Transfer from all the remaining zones ---*/
      for (auto jZone = 0u; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone && interface_container[iZone][jZone] != nullptr){
          DeformMesh = Transfer_Data(jZone, iZone);
          if (DeformMesh) UpdateMesh+=1;
        }
      }
      /*--- If a mesh update is required due to the transfer of data ---*/
      if (UpdateMesh > 0) DynamicMeshUpdate(iZone, TimeIter);

    }

      /*--- Loop over the number of zones (IZONE) ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);
      config_container[iZone]->Set_StartTime(SU2_MPI::Wtime());
      driver_config->SetOuterIter(iOuter_Iter);

      /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
      iteration_container[iZone][INST_0]->Solve(output_container[iZone], integration_container, geometry_container,
                                                solver_container, numerics_container, config_container,
                                                surface_movement, grid_movement, FFDBox, iZone, INST_0);

      /*--- A corrector step can help preventing numerical instabilities ---*/
      Corrector(iZone);

    }

    if (OuterConvergence(iOuter_Iter)) break;

  }

}

void CMultizoneDriver::Corrector(unsigned short val_iZone) {

  if (config_container[val_iZone]->GetRelaxation())
    iteration_container[val_iZone][INST_0]->Relaxation(output_container[ZONE_0], integration_container,
                                            geometry_container, solver_container, numerics_container, config_container,
                                            surface_movement, grid_movement, FFDBox, val_iZone, INST_0);
}

bool CMultizoneDriver::OuterConvergence(unsigned long OuterIter) {

  /*--- Update the residual for the all the zones. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Account for all the solvers in this zone. ---*/

    auto solvers = solver_container[iZone][INST_0][MESH_0];

    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      if (solvers[iSol] != nullptr) {
        solvers[iSol]->ComputeResidual_Multizone(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
        solvers[iSol]->GetNodes()->Set_BGSSolution_k();
      }
    }

    /*--- Make sure that everything is loaded into the output container. ---*/

    output_container[iZone]->SetHistory_Output(geometry_container[iZone][INST_0][MESH_0], solvers, config_container[iZone]);

  }

  /*--- Print out the convergence data to screen and history file. ---*/

  driver_output->SetMultizoneHistory_Output(output_container, config_container, driver_config,
                                            driver_config->GetTimeIter(), driver_config->GetOuterIter());

  return driver_output->GetConvergence();

}

void CMultizoneDriver::Update() {

  /*--- For enabling a consistent restart, we need to update the mesh with the interface information that introduces displacements --*/
  /*--- Loop over the number of zones (IZONE) ---*/
  for (iZone = 0; iZone < nZone; iZone++){

    unsigned short UpdateMesh = 0;

    /*--- Transfer from all the remaining zones (JZONE != IZONE)---*/
    for (auto jZone = 0u; jZone < nZone; jZone++){
      /*--- The target zone is iZone ---*/
      if (jZone != iZone){
        UpdateMesh += Transfer_Data(jZone, iZone);
      }
    }
    /*--- If a mesh update is required due to the transfer of data ---*/
    if (UpdateMesh > 0) DynamicMeshUpdate(iZone, TimeIter);

    iteration_container[iZone][INST_0]->Update(output_container[iZone], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, iZone, INST_0);

    /*--- Set the Convergence_FSI boolean to false for the next time step ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS-1; iSol++){
      if (integration_container[iZone][INST_0][iSol] != nullptr){
        integration_container[iZone][INST_0][iSol]->SetConvergence_FSI(false);
      }
    }
  }

}

void CMultizoneDriver::Output(unsigned long TimeIter) {

  /*--- Time the output for performance benchmarking. ---*/

  StopTime = SU2_MPI::Wtime();

  UsedTimeCompute += StopTime-StartTime;

  StartTime = SU2_MPI::Wtime();

  bool wrote_files = false;

  for (iZone = 0; iZone < nZone; iZone++){
    wrote_files = output_container[iZone]->SetResult_Files(geometry_container[iZone][INST_0][MESH_0],
                                                            config_container[iZone],
                                                            solver_container[iZone][INST_0][MESH_0], TimeIter, StopCalc );
  }

  if (wrote_files){

    StopTime = SU2_MPI::Wtime();

    UsedTimeOutput += StopTime-StartTime;
    OutputCount++;
    BandwidthSum = config_container[ZONE_0]->GetRestart_Bandwidth_Agg();

    StartTime = SU2_MPI::Wtime();

    driver_config->Set_StartTime(StartTime);
  }

}

void CMultizoneDriver::DynamicMeshUpdate(unsigned long TimeIter) {

  bool AnyDeformMesh = false;

  for (iZone = 0; iZone < nZone; iZone++) {
    const auto harmonic_balance = (config_container[iZone]->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
    /*--- Dynamic mesh update ---*/
    if ((config_container[iZone]->GetGrid_Movement()) && (!harmonic_balance) && (!fsi)) {
      iteration_container[iZone][INST_0]->SetGrid_Movement(geometry_container[iZone][INST_0],surface_movement[iZone],
                                                           grid_movement[iZone][INST_0], solver_container[iZone][INST_0],
                                                           config_container[iZone], 0, TimeIter);
      AnyDeformMesh = true;
    }
  }
  /*--- Update the wall distances if the mesh was deformed. ---*/
  if (AnyDeformMesh) {
    CGeometry::ComputeWallDistance(config_container, geometry_container);
  }
}

void CMultizoneDriver::DynamicMeshUpdate(unsigned short val_iZone, unsigned long TimeIter) {

  auto iteration = iteration_container[val_iZone][INST_0];

  /*--- Legacy dynamic mesh update - Only if GRID_MOVEMENT = YES ---*/
  if (config_container[val_iZone]->GetGrid_Movement()) {
    iteration->SetGrid_Movement(geometry_container[val_iZone][INST_0],surface_movement[val_iZone],
                                grid_movement[val_iZone][INST_0], solver_container[val_iZone][INST_0],
                                config_container[val_iZone], 0, TimeIter);
  }

  /*--- New solver - all the other routines in SetGrid_Movement should be adapted to this one ---*/
  /*--- Works if DEFORM_MESH = YES ---*/
  iteration->SetMesh_Deformation(geometry_container[val_iZone][INST_0],
                                 solver_container[val_iZone][INST_0][MESH_0],
                                 numerics_container[val_iZone][INST_0][MESH_0],
                                 config_container[val_iZone], RECORDING::CLEAR_INDICES);

  /*--- Update the wall distances if the mesh was deformed. ---*/
  if (config_container[val_iZone]->GetGrid_Movement() ||
      config_container[val_iZone]->GetDeform_Mesh()) {
    CGeometry::ComputeWallDistance(config_container, geometry_container);
  }
}

bool CMultizoneDriver::Transfer_Data(unsigned short donorZone, unsigned short targetZone) {

  bool UpdateMesh = false;

  int donorSolver = -1, targetSolver = -1;

  /*--- Select the transfer method according to the magnitudes being transferred ---*/

  switch(interface_types[donorZone][targetZone]) {

    case SLIDING_INTERFACE:
    {
      donorSolver  = FLOW_SOL;
      targetSolver = FLOW_SOL;

      /*--- Additional transfer for turbulence variables. ---*/
      if (config_container[targetZone]->GetKind_Solver() == MAIN_SOLVER::RANS ||
          config_container[targetZone]->GetKind_Solver() == MAIN_SOLVER::INC_RANS)
      {
        interface_container[donorZone][targetZone]->BroadcastData(
          *interpolator_container[donorZone][targetZone].get(),
          solver_container[donorZone][INST_0][MESH_0][TURB_SOL],
          solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
          geometry_container[donorZone][INST_0][MESH_0],
          geometry_container[targetZone][INST_0][MESH_0],
          config_container[donorZone],
          config_container[targetZone]);
      }
      break;
    }
    case CONJUGATE_HEAT_FS:
    {
      donorSolver  = FLOW_SOL;
      targetSolver = HEAT_SOL;
      break;
    }
    case CONJUGATE_HEAT_WEAKLY_FS:
    {
      donorSolver  = HEAT_SOL;
      targetSolver = HEAT_SOL;
      break;
    }
    case CONJUGATE_HEAT_SF:
    {
      donorSolver  = HEAT_SOL;
      targetSolver = FLOW_SOL;
      break;
    }
    case CONJUGATE_HEAT_WEAKLY_SF:
    {
      donorSolver  = HEAT_SOL;
      targetSolver = HEAT_SOL;
      break;
    }
    case BOUNDARY_DISPLACEMENTS:
    {
      donorSolver  = FEA_SOL;
      targetSolver = MESH_SOL;
      UpdateMesh = true;
      break;
    }
    case FLOW_TRACTION:
    {
      donorSolver  = FLOW_SOL;
      targetSolver = FEA_SOL;
      break;
    }
    case NO_TRANSFER:
    case ZONES_ARE_EQUAL:
    case NO_COMMON_INTERFACE:
      break;
    default:
      if(rank == MASTER_NODE)
        cout << "WARNING: One of the intended interface transfer routines is not "
             << "known to the chosen driver and has not been executed." << endl;
      break;
  }

  if(donorSolver >= 0 && targetSolver >= 0) {
    interface_container[donorZone][targetZone]->BroadcastData(
      *interpolator_container[donorZone][targetZone].get(),
      solver_container[donorZone][INST_0][MESH_0][donorSolver],
      solver_container[targetZone][INST_0][MESH_0][targetSolver],
      geometry_container[donorZone][INST_0][MESH_0],
      geometry_container[targetZone][INST_0][MESH_0],
      config_container[donorZone],
      config_container[targetZone]);
  }

  return UpdateMesh;
}

bool CMultizoneDriver::Monitor(unsigned long TimeIter){

  /*--- Check whether the inner solver has converged --- */

  if (driver_config->GetTime_Domain() == NO){

    const auto OuterIter  = driver_config->GetOuterIter();
    const auto nOuterIter = driver_config->GetnOuter_Iter();

    const auto InnerConvergence = driver_output->GetConvergence();
    const bool MaxIterationsReached = (OuterIter+1 >= nOuterIter);

    if ((MaxIterationsReached || InnerConvergence) && (rank == MASTER_NODE)) {
      cout << endl << "----------------------------- Solver Exit -------------------------------" << endl;
      if (InnerConvergence) cout << "All convergence criteria satisfied." << endl;
      else cout << endl << "Maximum number of iterations reached (OUTER_ITER = " << OuterIter+1 << ") before convergence." << endl;
      driver_output->PrintConvergenceSummary();
      cout << "-------------------------------------------------------------------------" << endl;
    }

    return (MaxIterationsReached || InnerConvergence);
  }
  else { // i.e. unsteady simulation

    /*--- Check whether the outer time integration has reached the final time ---*/
    const auto TimeConvergence = GetTimeConvergence();

    const auto nTimeIter = driver_config->GetnTime_Iter();
    const auto MaxTime   = driver_config->GetMax_Time();
    const auto CurTime   = driver_output->GetHistoryFieldValue("CUR_TIME");

    const bool FinalTimeReached = (CurTime >= MaxTime);
    const bool MaxIterationsReached = (TimeIter+1 >= nTimeIter);

    if ((TimeConvergence || FinalTimeReached || MaxIterationsReached) && (rank == MASTER_NODE)){
      cout << endl << "----------------------------- Solver Exit -------------------------------";
      if (TimeConvergence)  cout << endl << "All windowed time-averaged convergence criteria are fullfilled." << endl;
      if (FinalTimeReached) cout << endl << "Maximum time reached (MAX_TIME = " << MaxTime << "s)." << endl;
      else cout << endl << "Maximum number of time iterations reached (TIME_ITER = " << nTimeIter << ")." << endl;
      cout << "-------------------------------------------------------------------------" << endl;
    }

    return (FinalTimeReached || MaxIterationsReached);
  }

}

bool CMultizoneDriver::GetTimeConvergence() const{
    return output_container[ZONE_0]->GetCauchyCorrectedTimeConvergence(config_container[ZONE_0]);
}
