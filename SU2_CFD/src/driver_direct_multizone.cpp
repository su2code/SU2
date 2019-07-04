/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving multi-zone problems.
 * \author R. Sanchez, O. Burghardt
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

#include "../include/driver_structure.hpp"
#include "../include/definition_structure.hpp"


CMultizoneDriver::CMultizoneDriver(char* confFile,
                       unsigned short val_nZone,
                       unsigned short val_nDim,
                       SU2_Comm MPICommunicator) : CDriver(confFile,
                                                          val_nZone,
                                                          val_nDim,
                                                          MPICommunicator) {

  /*--- Initialize the counter for TimeIter ---*/
  TimeIter = 0;

  /*--- Initialize some useful booleans ---*/
  fsi = false; cht = false;

  /*--- Structure for multizone convergence ---*/
  init_res     = new su2double*[nZone];
  residual     = new su2double*[nZone];
  residual_rel = new su2double*[nZone];
  nVarZone     = new unsigned short[nZone];

  for (iZone = 0; iZone < nZone; iZone++){
    nVarZone[iZone] = 0;
    /*--- Account for all the solvers ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL)
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
    case EULER: case NAVIER_STOKES: case RANS:
      fluid_zone = true;
      break;
    case FEM_ELASTICITY:
      structural_zone = true;
      break;
    case HEAT_EQUATION_FVM:
      heat_zone = true;
      break;
    }
  }

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
      case RIGID_MOTION: case DEFORMING:
      case EXTERNAL: case EXTERNAL_ROTATION:
      case AEROELASTIC: case AEROELASTIC_RIGID_MOTION:
      case ELASTICITY:
        prefixed_motion[iZone] = true; break;
      case FLUID_STRUCTURE: case FLUID_STRUCTURE_STATIC:
      case STEADY_TRANSLATION: case MOVING_WALL: case ROTATING_FRAME:
      case NO_MOVEMENT: case GUST: default:
        prefixed_motion[iZone] = false; break;
    }
  }

}

CMultizoneDriver::~CMultizoneDriver(void) {

  for (iZone = 0; iZone < nZone; iZone++){
    delete [] init_res[iZone];
    delete [] residual[iZone];
    delete [] residual_rel[iZone];
  }

  delete [] init_res;
  delete [] residual;
  delete [] residual_rel;
  delete [] nVarZone;

  delete [] prefixed_motion;

}

void CMultizoneDriver::StartSolver() {

  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;

  if (rank == MASTER_NODE){
    cout << endl <<"Simulation Run using the Multizone Driver" << endl;
    if (driver_config->GetTime_Domain())
      cout << "The simulation will run for " << driver_config->GetnTime_Iter() << " time steps." << endl;
  }

  /*--- Set the initial time iteration to the restart iteration. ---*/
  if (driver_config->GetRestart())
    TimeIter = driver_config->GetRestart_Iter();

  /*--- Run the problem until the number of time iterations required is reached. ---*/
  while ( TimeIter < driver_config->GetnTime_Iter() ) {

    /*--- Perform some preprocessing before starting the time-step simulation. ---*/

    Preprocess(TimeIter);

    /*--- Run a block iteration of the multizone problem. ---*/

    switch (driver_config->GetKind_MZSolver()){
      case MZ_BLOCK_GAUSS_SEIDEL: Run_GaussSeidel(); break;  // Block Gauss-Seidel iteration
      case MZ_BLOCK_JACOBI: Run_Jacobi(); break;             // Block-Jacobi iteration
      default: Run_GaussSeidel(); break;
    }

    /*--- Update the solution for dual time stepping strategy ---*/

    Update();

    /*--- Monitor the computations after each iteration. ---*/

    Monitor(TimeIter);

    /*--- Output the solution in files. ---*/

    Output(TimeIter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    TimeIter++;

  }

}

void CMultizoneDriver::Preprocess(unsigned long TimeIter) {

  bool unsteady = driver_config->GetTime_Domain();

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the value of the external iteration to TimeIter. -------------------------------------*/
    /*--- TODO: This should be generalised for an homogeneous criteria throughout the code. --------*/
    config_container[iZone]->SetExtIter(TimeIter);

    /*--- Store the current physical time in the config container, as
     this can be used for verification / MMS. This should also be more
     general once the drivers are more stable. ---*/
    
    if (unsteady)
      config_container[iZone]->SetPhysicalTime(static_cast<su2double>(TimeIter)*config_container[iZone]->GetDelta_UnstTimeND());
    else
      config_container[iZone]->SetPhysicalTime(0.0);
    
    /*--- Read the target pressure for inverse design. ---------------------------------------------*/
    /*--- TODO: This routine should be taken out of output, and made general for multiple zones. ---*/
    if (config_container[iZone]->GetInvDesign_Cp() == YES)
      output->SetCp_InverseDesign(solver_container[iZone][INST_0][MESH_0][FLOW_SOL],
          geometry_container[iZone][INST_0][MESH_0], config_container[iZone], TimeIter);

    /*--- Read the target heat flux ----------------------------------------------------------------*/
    /*--- TODO: This routine should be taken out of output, and made general for multiple zones. ---*/
    if (config_container[iZone]->GetInvDesign_HeatFlux() == YES)
      output->SetHeatFlux_InverseDesign(solver_container[iZone][INST_0][MESH_0][FLOW_SOL],
          geometry_container[iZone][INST_0][MESH_0], config_container[iZone], TimeIter);

    /*--- Set the initial condition for EULER/N-S/RANS ---------------------------------------------*/
    /*--- For FSI, this is set after the mesh has been moved. --------------------------------------*/
    if ((config_container[iZone]->GetKind_Solver() ==  EULER) ||
        (config_container[iZone]->GetKind_Solver() ==  NAVIER_STOKES) ||
        (config_container[iZone]->GetKind_Solver() ==  RANS) ) {
        if(!fsi) solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone][INST_0], solver_container[iZone][INST_0], config_container[iZone], TimeIter);
    }

  }

#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  /*--- Run a predictor step ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    if (config_container[iZone]->GetPredictor())
      iteration_container[iZone][INST_0]->Predictor(output, integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  /*--- Perform a dynamic mesh update if required. ---*/

  DynamicMeshUpdate(TimeIter);

  /*--- Updating zone interface communication patterns for unsteady problems with pre-fixed motion in the config file ---*/
  if ( unsteady ) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (unsigned short jZone = 0; jZone < nZone; jZone++){
        if(jZone != iZone && interpolator_container[iZone][jZone] != NULL && prefixed_motion[iZone])
          interpolator_container[iZone][jZone]->Set_TransferCoeff(config_container);
      }
    }
  }

}

void CMultizoneDriver::Run_GaussSeidel() {

  unsigned long iOuter_Iter;
  unsigned short jZone, UpdateMesh;
  bool DeformMesh = false;
  unsigned long ExtIter = 0;
  bool Convergence = false;

  unsigned long OuterIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetOuterIter(OuterIter);

  /*--- Loop over the number of outer iterations ---*/
  for (iOuter_Iter = 0; iOuter_Iter < driver_config->GetnOuter_Iter(); iOuter_Iter++){

    /*--- Loop over the number of zones (IZONE) ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- In principle, the mesh does not need to be updated ---*/
      UpdateMesh = 0;

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);

      /*--- Transfer from all the remaining zones ---*/
      for (jZone = 0; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone){
          DeformMesh = Transfer_Data(jZone, iZone);
          if (DeformMesh) UpdateMesh+=1;
        }
      }
      /*--- If a mesh update is required due to the transfer of data ---*/
      if (UpdateMesh > 0) DynamicMeshUpdate(iZone, ExtIter);

      /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
      iteration_container[iZone][INST_0]->Solve(output, integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

      /*--- A corrector step can help preventing numerical instabilities ---*/
      Corrector(iZone);

    }

    /*--- This is temporary. Each zone has to be monitored independently. Right now, fixes CHT output. ---*/
    Monitor(iOuter_Iter);

    Convergence = OuterConvergence(iOuter_Iter);

    if (Convergence) break;

  }

}

void CMultizoneDriver::Run_Jacobi() {

  unsigned long iOuter_Iter;
  unsigned short jZone, UpdateMesh;
  bool DeformMesh = false;
  unsigned long ExtIter = 0;
  bool Convergence = false;

  unsigned long OuterIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetOuterIter(OuterIter);

  /*--- Loop over the number of outer iterations ---*/
  for (iOuter_Iter = 0; iOuter_Iter < driver_config->GetnOuter_Iter(); iOuter_Iter++){

    /*--- Transfer from all zones ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- In principle, the mesh does not need to be updated ---*/
      UpdateMesh = 0;

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);

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

      /*--- Loop over the number of zones (IZONE) ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);

      /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
      iteration_container[iZone][INST_0]->Solve(output, integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

      /*--- A corrector step can help preventing numerical instabilities ---*/
      Corrector(iZone);

    }

    Convergence = OuterConvergence(iOuter_Iter);

    if (Convergence) break;

  }

}

void CMultizoneDriver::Corrector(unsigned short val_iZone) {

    if (config_container[val_iZone]->GetRelaxation())
      iteration_container[val_iZone][INST_0]->Relaxation(output, integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, val_iZone, INST_0);

}

bool CMultizoneDriver::OuterConvergence(unsigned long OuterIter) {

  bool Convergence = false;
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  unsigned short iRes, iVar;
  unsigned short nVarSol;

  /*--- Compute the residual for the all the zones ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    iVar = 0; // Initialize the variable index for each zone

    /*--- Account for all the solvers ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      /*-- If the solver position iSol is enabled --*/
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL){
        nVarSol = solver_container[iZone][INST_0][MESH_0][iSol]->GetnVar();

        /*--- Compute the block residual on each solver ---*/
        solver_container[iZone][INST_0][MESH_0][iSol]->ComputeResidual_Multizone(geometry_container[iZone][INST_0][MESH_0],
            config_container[iZone]);

        /*--- Loop over all the variables in the solver ---*/
        for (iRes = 0; iRes < nVarSol; iRes++){
          /*--- Store the log10 of the residual value ---*/
          residual[iZone][iVar] = log10(solver_container[iZone][INST_0][MESH_0][iSol]->GetRes_BGS(iRes));
          /*--- If it is the first iteration, the init_res is the current residual ---*/
          if (OuterIter == 0) init_res[iZone][iVar] = residual[iZone][iVar];
          /*--- residual_rel checks the difference in order of magnitude between the current and the initial residual ---*/
          residual_rel[iZone][iVar] = fabs(residual[iZone][iRes] - init_res[iZone][iRes]);
          /*--- Move the position of the container ---*/
          iVar++;
        }
      }
    }
  }

  /*--- This still has to be generalised ---*/
  if (fsi){

    /*--- This should be possible to change dinamically in the config ---*/
    if (rank == MASTER_NODE){

      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << endl;
      cout << "Convergence summary for BGS iteration ";
      cout << OuterIter << endl;
      cout << endl;
      /*--- TODO: This is a workaround until the TestCases.py script incorporates new classes for nested loops. ---*/
      cout << "Iter[ID]" << "    BGSRes[Rho]" << "   BGSRes[RhoE]" << "     BGSRes[Ux]" << "     BGSRes[Uy]" << endl;
      cout.precision(6); cout.setf(ios::fixed, ios::floatfield);
      cout.width(8); cout << OuterIter*1000;
      cout.width(15); cout << residual[ZONE_0][0];
      cout.width(15); cout << residual[ZONE_0][3];
      cout.width(15); cout << residual[ZONE_1][0];
      cout.width(15); cout << residual[ZONE_1][1];
      cout << endl;
    }
    for (iZone = 0; iZone < nZone; iZone++){
      if (config_container[iZone]->GetKind_Solver() == FEM_ELASTICITY){
        integration_container[iZone][INST_0][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[iZone][INST_0][MESH_0], config_container[iZone], solver_container[iZone][INST_0][MESH_0][FEA_SOL], OuterIter);
        Convergence = integration_container[iZone][INST_0][FEA_SOL]->GetConvergence_FSI();
      }
    }
  }

  /*--- Update the residual for the all the zones ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    /*--- Accounting for all the solvers ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      /*-- If the solver position iSol is enabled --*/
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL){
        solver_container[iZone][INST_0][MESH_0][iSol]->UpdateSolution_BGS(geometry_container[iZone][INST_0][MESH_0],
            config_container[iZone]);}
    }
  }

  if (rank == MASTER_NODE) cout.setf(ios::scientific, ios::floatfield);

  /*-----------------------------------------------------------------*/
  /*-------------------- Output FSI history -------------------------*/
  /*-----------------------------------------------------------------*/
  if (fsi){
    bool ZONE_FLOW=0, ZONE_FEA=1;
    /*--- This is a hack to test it works. ---*/
    for (iZone = 0; iZone < nZone; iZone++){
      if (config_container[iZone]->GetKind_Solver() == FEM_ELASTICITY) ZONE_FEA = iZone;
      if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) ZONE_FLOW = iZone;
    }
    output->SpecialOutput_FSI(&FSIHist_file, geometry_container, solver_container,
        config_container, integration_container, 0,
        ZONE_FLOW, ZONE_FEA, false);
  }

  return Convergence;

}

void CMultizoneDriver::Update() {

  unsigned short jZone, UpdateMesh;
  bool DeformMesh = false;

  /*--- For enabling a consistent restart, we need to update the mesh with the interface information that introduces displacements --*/
  /*--- Loop over the number of zones (IZONE) ---*/
  for (iZone = 0; iZone < nZone; iZone++){

    UpdateMesh = 0;

      /*--- Transfer from all the remaining zones (JZONE != IZONE)---*/
      for (jZone = 0; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone){
          DeformMesh = Transfer_Data(jZone, iZone);
          if (DeformMesh) UpdateMesh += 1;
        }
      }
    /*--- If a mesh update is required due to the transfer of data ---*/
    if (UpdateMesh > 0) DynamicMeshUpdate(iZone, ExtIter);

    iteration_container[iZone][INST_0]->Update(output, integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, iZone, INST_0);

    /*--- Set the Convergence_FSI boolean to false for the next time step ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL)
        integration_container[iZone][INST_0][iSol]->SetConvergence_FSI(false);
    }
  }

}

void CMultizoneDriver::Output(unsigned long TimeIter) {

  bool output_files = false;

  /*--- Determine whether a solution needs to be written
   after the current iteration ---*/

  if (

      /*--- General if statements to print output statements ---*/

      (TimeIter+1 >= config_container[ZONE_0]->GetnTime_Iter()) || (StopCalc) ||

      /*--- Unsteady problems ---*/

      (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
        (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_STEPPING)) &&
       ((TimeIter == 0) || (ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0))) ||

      ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND) &&
       ((TimeIter == 0) || ((TimeIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0) ||
                           ((TimeIter-1) % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) ||

      ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND) &&
       ((TimeIter == 0) || ((TimeIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) ||

      ((config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC) &&
       ((TimeIter == 0) || (TimeIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0))) ||

      /*--- No inlet profile file found. Print template. ---*/

      (config_container[ZONE_0]->GetWrt_InletFile())

      ) {

    output_files = true;

  }

  /*--- Determine whether a solution doesn't need to be written
   after the current iteration ---*/

  if (config_container[ZONE_0]->GetFixed_CL_Mode()) {
    if (config_container[ZONE_0]->GetnExtIter()-config_container[ZONE_0]->GetIter_dCL_dAlpha() - 1 < ExtIter) output_files = false;
    if (config_container[ZONE_0]->GetnExtIter() - 1 == ExtIter) output_files = true;
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

    output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, TimeIter, nZone);


    /*--- Execute the routine for writing special output. ---*/
    output->SetSpecial_Output(solver_container, geometry_container, config_container, TimeIter, nZone);


    if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl << endl;

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

void CMultizoneDriver::DynamicMeshUpdate(unsigned long ExtIter) {

  bool harmonic_balance;

  for (iZone = 0; iZone < nZone; iZone++) {
   harmonic_balance = (config_container[iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE);
    /*--- Dynamic mesh update ---*/
    if ((config_container[iZone]->GetGrid_Movement()) && (!harmonic_balance) && (!fsi)) {
      iteration_container[iZone][INST_0]->SetGrid_Movement(geometry_container, surface_movement, grid_movement, FFDBox, solver_container, config_container, iZone, INST_0, 0, ExtIter );
    }
  }
}

void CMultizoneDriver::DynamicMeshUpdate(unsigned short val_iZone, unsigned long ExtIter) {

  iteration_container[val_iZone][INST_0]->SetGrid_Movement(geometry_container,surface_movement, grid_movement, FFDBox, solver_container,
        config_container, val_iZone, INST_0, 0, ExtIter);

}

bool CMultizoneDriver::Transfer_Data(unsigned short donorZone, unsigned short targetZone) {

  bool UpdateMesh = false;

  /*--- Select the transfer method according to the magnitudes being transferred ---*/

  if (transfer_types[donorZone][targetZone] == SLIDING_INTERFACE) {
    transfer_container[donorZone][targetZone]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
    if (config_container[targetZone]->GetKind_Solver() == RANS)
      transfer_container[donorZone][targetZone]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][TURB_SOL],solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
                                                                         geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                         config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_FS) {
    transfer_container[donorZone][targetZone]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_FS) {
    transfer_container[donorZone][targetZone]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_SF) {
    transfer_container[donorZone][targetZone]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_SF) {
    transfer_container[donorZone][targetZone]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == STRUCTURAL_DISPLACEMENTS) {
    transfer_container[donorZone][targetZone]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
    UpdateMesh = true;
  }
  else if (transfer_types[donorZone][targetZone] == FLOW_TRACTION) {
    transfer_container[donorZone][targetZone]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FEA_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if ((transfer_types[donorZone][targetZone] == NO_TRANSFER)
           || (transfer_types[donorZone][targetZone] == ZONES_ARE_EQUAL)
           || (transfer_types[donorZone][targetZone] == NO_COMMON_INTERFACE)) { }
  else {
    cout << "WARNING: One of the intended interface transfer routines is not known to the chosen driver and has not been executed." << endl;
  }

  return UpdateMesh;
}
