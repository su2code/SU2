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


CMultizoneDriver::CMultizoneDriver(char* confFile,
                       unsigned short val_nZone,
                       unsigned short val_nDim,
                       bool val_periodic,
                       SU2_Comm MPICommunicator) : CDriver(confFile,
                                                          val_nZone,
                                                          val_nDim,
                                                          val_periodic,
                                                          MPICommunicator) {

  /*--- Initialize some useful booleans ---*/
  fsi = false;

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

}

void CMultizoneDriver::Run() {

  unsigned long iOuter_Iter;
  unsigned short jZone;
  bool UpdateMesh = false;
  unsigned long ExtIter = 0;
  bool Convergence = false;

  unsigned long OuterIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetOuterIter(OuterIter);

  /*--- First, some preprocessing is required. This will determine the kind of problem(s) involved ---*/
  Preprocess();

  /*--- A predictor can be useful for different kinds of problems ---*/
  Predictor();

  /*--- Loop over the number of outer iterations ---*/
  for (iOuter_Iter = 0; iOuter_Iter < driver_config->GetnOuter_Iter(); iOuter_Iter++){

    /*--- Loop over the number of zones (IZONE) ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);

      /*--- Transfer from all the remaining zones ---*/
      for (jZone = 0; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone){
          UpdateMesh = Transfer_Data(jZone, iZone);
          /*--- If a mesh update is required due to the transfer of data ---*/
          if (UpdateMesh) DynamicMeshUpdate(iZone, ExtIter);
        }
      }

//      if (config_container[iZone]->GetnInner_Iter() == 1)
//        /*--- Run one iteration of the solver ---*/
//        iteration_container[iZone][INST_0]->Iterate(output, integration_container, geometry_container, solver_container,
//            numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);
//      else
        /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
        iteration_container[iZone][INST_0]->Iterate_Block(output, integration_container, geometry_container, solver_container,
            numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

      /*--- A relaxation step can help preventing numerical instabilities ---*/
      Relaxation();

    }

    Convergence = OuterConvergence(iOuter_Iter);

    if (Convergence) break;

  }

  /*--- Temporary ---*/
  if (fsi) Update();

}

void CMultizoneDriver::Preprocess() {

  bool time_domain = driver_config->GetTime_Domain();

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
    case HEAT_EQUATION: case HEAT_EQUATION_FVM:
      heat_zone = true;
      break;
    }
  }

  /*--- If the problem has FSI properties ---*/
  if (fluid_zone && structural_zone) fsi = true;

  if (fsi){
    /*--- Loop over the number of zones (IZONE) ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- For the Fluid zones ---*/
      if ((config_container[iZone]->GetKind_Solver()==EULER) ||
          (config_container[iZone]->GetKind_Solver()==NAVIER_STOKES) ||
          (config_container[iZone]->GetKind_Solver()==RANS)){

        /*--- If there is a restart, we need to get the old geometry from the fluid field ---*/
        bool restart = (config_container[iZone]->GetRestart() || config_container[iZone]->GetRestart_Flow());
        ExtIter = config_container[iZone]->GetExtIter();

        /*--- Different treatment required depending of whether it's time domain problem or not ---*/
        if (restart && time_domain && (long)ExtIter == config_container[iZone]->GetUnst_RestartIter()) {
          if (rank == MASTER_NODE) cout << "Restarting geometry structure for previous time steps." << endl;
          solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Restart_OldGeometry(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);
        } else if (restart && !time_domain){
          if (rank == MASTER_NODE) cout << "Restarting geometry structure from converged solution." << endl;
          solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[iZone][INST_0], solver_container[iZone][INST_0], config_container[iZone], 0, false);
        }

      }

    }

  }

}

void CMultizoneDriver::Predictor() {

  for (iZone = 0; iZone < nZone; iZone++){
    if (config_container[iZone]->GetPredictor())
      iteration_container[iZone][INST_0]->Predictor(output, integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

}

void CMultizoneDriver::Relaxation() {

    if (config_container[iZone]->GetRelaxation())
      iteration_container[iZone][INST_0]->Relaxation(output, integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

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
        solver_container[iZone][INST_0][MESH_0][iSol]->ComputeResidual_BGS(geometry_container[iZone][INST_0][MESH_0],
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

    integration_container[ZONE_1][INST_0][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[ZONE_1][INST_0][MESH_0], config_container[ZONE_1], solver_container[ZONE_1][INST_0][MESH_0][FEA_SOL], OuterIter);

    Convergence = integration_container[ZONE_1][INST_0][FEA_SOL]->GetConvergence_FSI();
  }

//  /*--- Update the residual for the all the zones ---*/
//  for (iZone = 0; iZone < nZone; iZone++){
//    /*--- Accounting for all the solvers ---*/
//    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
//      /*-- If the solver position iSol is enabled --*/
//      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL){
//        solver_container[iZone][INST_0][MESH_0][iSol]->UpdateSolution_BGS(geometry_container[iZone][INST_0][MESH_0],
//                                                                          config_container[iZone]);}
//    }
//  }

  if (rank == MASTER_NODE) cout.setf(ios::scientific, ios::floatfield);

  /*-----------------------------------------------------------------*/
  /*-------------------- Output FSI history -------------------------*/
  /*-----------------------------------------------------------------*/
  if (fsi){
  output->SpecialOutput_FSI(&FSIHist_file, geometry_container, solver_container,
                            config_container, integration_container, 0,
                            ZONE_0, ZONE_1, false);
  }

  return Convergence;

}

void CMultizoneDriver::Update() {

  unsigned short jZone;
  bool UpdateMesh;

  /*--- For enabling a consistent restart, we need to update the mesh with the interface information that introduces displacements --*/
  if (fsi){
  /*--- Loop over the number of zones (IZONE) ---*/
  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Transfer from all the remaining zones ---*/
    for (jZone = 0; jZone < nZone; jZone++){
      /*--- The target zone is iZone ---*/
      if (jZone != iZone){
        UpdateMesh = Transfer_Data(jZone, iZone);
        /*--- If a mesh update is required due to the transfer of data ---*/
        if (UpdateMesh) DynamicMeshUpdate(iZone, ExtIter);
      }
    }

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

  bool MatchingMesh = config_container[targetZone]->GetMatchingMesh();
  bool UpdateMesh = false;

  /*--- Select the transfer method and the appropriate mesh properties (matching or nonmatching mesh) ---*/

  switch (config_container[targetZone]->GetKind_TransferMethod()) {

  case BROADCAST_DATA:
      if (MatchingMesh) {
        if (transfer_types[donorZone][targetZone] == SLIDING_INTERFACE) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
          if (config_container[targetZone]->GetKind_Solver() == RANS)
            transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][INST_0][MESH_0][TURB_SOL],solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
                geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_FS) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_FS) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_SF) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_SF) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == STRUCTURAL_DISPLACEMENTS) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
          UpdateMesh = true;
        }
        else if (transfer_types[donorZone][targetZone] == FLOW_TRACTION) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FEA_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if ((transfer_types[donorZone][targetZone] == NO_TRANSFER)
                 || (transfer_types[donorZone][targetZone] == ZONES_ARE_EQUAL)
                 || (transfer_types[donorZone][targetZone] == NO_COMMON_INTERFACE)) { }
        else {
          cout << "WARNING: One of the specified interface transfer routines is not known to the chosen driver and has not been executed." << endl;
        }
      }
      else {
        if (transfer_types[donorZone][targetZone] == SLIDING_INTERFACE) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
          if (config_container[targetZone]->GetKind_Solver() == RANS)
            transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][INST_0][MESH_0][TURB_SOL],solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
                geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_FS) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_FS) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_SF) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_SF) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if (transfer_types[donorZone][targetZone] == STRUCTURAL_DISPLACEMENTS) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
          UpdateMesh = true;
        }
        else if (transfer_types[donorZone][targetZone] == FLOW_TRACTION) {
          transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FEA_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
        }
        else if ((transfer_types[donorZone][targetZone] == NO_TRANSFER)
                 || (transfer_types[donorZone][targetZone] == ZONES_ARE_EQUAL)
                 || (transfer_types[donorZone][targetZone] == NO_COMMON_INTERFACE)) { }
        else {
          cout << "WARNING: One of the intended interface transfer routines is not known to the chosen driver and has not been executed." << endl;
        }
      }
      break;

  case SCATTER_DATA:
    if (MatchingMesh) {
      if (transfer_types[donorZone][targetZone] == SLIDING_INTERFACE) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
        if (config_container[targetZone]->GetKind_Solver() == RANS)
          transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][INST_0][MESH_0][TURB_SOL],solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_FS) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_FS) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_SF) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == STRUCTURAL_DISPLACEMENTS) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
        UpdateMesh = true;
      }
      else if (transfer_types[donorZone][targetZone] == FLOW_TRACTION) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FEA_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_SF) {
        transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if ((transfer_types[donorZone][targetZone] == NO_TRANSFER)
               || (transfer_types[donorZone][targetZone] == ZONES_ARE_EQUAL)
               || (transfer_types[donorZone][targetZone] == NO_COMMON_INTERFACE)) { }
      else {
        cout << "WARNING: One of the specified interface transfer routines is not known to the chosen driver and has not been executed." << endl;
      }
    }
    else {
      SU2_MPI::Error("Scatter method not implemented for non-matching meshes. ", CURRENT_FUNCTION);
    }
    break;

  case ALLGATHER_DATA:
    if (MatchingMesh) {
      SU2_MPI::Error("Allgather method not implemented for matching meshes. ", CURRENT_FUNCTION);
    }
    else {
      if (transfer_types[donorZone][targetZone] == SLIDING_INTERFACE) {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
        if (config_container[targetZone]->GetKind_Solver() == RANS)
          transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][INST_0][MESH_0][TURB_SOL],solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
              geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
              config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_FS) {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_FS) {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_SF) {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_SF) {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if (transfer_types[donorZone][targetZone] == STRUCTURAL_DISPLACEMENTS) {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
        UpdateMesh = true;
      }
      else if (transfer_types[donorZone][targetZone] == FLOW_TRACTION) {
        transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FEA_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone]);
      }
      else if ((transfer_types[donorZone][targetZone] == NO_TRANSFER)
               || (transfer_types[donorZone][targetZone] == ZONES_ARE_EQUAL)
               || (transfer_types[donorZone][targetZone] == NO_COMMON_INTERFACE)) { }
      else {
        cout << "WARNING: One of the specified interface transfer routines is not known to the chosen driver and can't be executed." << endl;
      }
    }
    break;
  }

  return UpdateMesh;
}
