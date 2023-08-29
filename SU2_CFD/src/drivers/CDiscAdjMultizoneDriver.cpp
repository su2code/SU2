/*!
 * \file CDiscAdjMultizoneDriver.cpp
 * \brief The main subroutines for driving adjoint multi-zone problems
 * \author O. Burghardt, P. Gomes, T. Albring, R. Sanchez
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

#include "../../include/drivers/CDiscAdjMultizoneDriver.hpp"
#include "../../include/solvers/CFEASolver.hpp"
#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIterationFactory.hpp"

CDiscAdjMultizoneDriver::CDiscAdjMultizoneDriver(char* confFile,
                                                 unsigned short val_nZone,
                                                 SU2_Comm MPICommunicator)

                        : CMultizoneDriver(confFile, val_nZone, MPICommunicator) {

  direct_nInst.resize(nZone,1);
  nInnerIter.resize(nZone);

  for (iZone = 0; iZone < nZone; iZone++)
    nInnerIter[iZone] = config_container[iZone]->GetnInner_Iter();

  Has_Deformation.resize(nZone) = false;


  FixPtCorrector.resize(nZone);
  LinSolver.resize(nZone);
  AdjRHS.resize(nZone);
  AdjSol.resize(nZone);

  direct_iteration = new CIteration**[nZone];
  direct_output = new COutput*[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Instantiate a direct iteration for each zone. ---*/

    direct_iteration[iZone] = new CIteration*[direct_nInst[iZone]];

    for(iInst = 0; iInst < direct_nInst[iZone]; iInst++) {

      switch (config_container[iZone]->GetKind_Solver()) {

        case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
          direct_iteration[iZone][iInst] = CIterationFactory::CreateIteration(MAIN_SOLVER::EULER, config_container[iZone]);
          break;
        case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
          direct_iteration[iZone][iInst] = CIterationFactory::CreateIteration(MAIN_SOLVER::INC_EULER, config_container[iZone]);
          break;
        case MAIN_SOLVER::DISC_ADJ_HEAT:
          direct_iteration[iZone][iInst] = CIterationFactory::CreateIteration(MAIN_SOLVER::HEAT_EQUATION, config_container[iZone]);
          break;
        case MAIN_SOLVER::DISC_ADJ_FEM:
          direct_iteration[iZone][iInst] = CIterationFactory::CreateIteration(MAIN_SOLVER::FEM_ELASTICITY, config_container[iZone]);
          break;
        default:
          SU2_MPI::Error("There is no discrete adjoint functionality for one of the specified solvers yet.",
                         CURRENT_FUNCTION);
      }
    }

    /*--- Instantiate a direct output to get the results of each direct zone. ---*/

    switch (config_container[iZone]->GetKind_Solver()) {

      case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
        direct_output[iZone] = COutputFactory::CreateOutput(MAIN_SOLVER::EULER, config_container[iZone], nDim);
        break;
      case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
        direct_output[iZone] = COutputFactory::CreateOutput(MAIN_SOLVER::INC_EULER, config_container[iZone], nDim);
        break;
      case MAIN_SOLVER::DISC_ADJ_HEAT:
        direct_output[iZone] = COutputFactory::CreateOutput(MAIN_SOLVER::HEAT_EQUATION, config_container[iZone], nDim);
        break;
      case MAIN_SOLVER::DISC_ADJ_FEM:
        direct_output[iZone] = COutputFactory::CreateOutput(MAIN_SOLVER::FEM_ELASTICITY, config_container[iZone], nDim);
        break;
      default:
        direct_output[iZone] = nullptr;
        break;
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
    delete direct_output[iZone];
  }

  delete[] direct_iteration;
  delete[] direct_output;

}

void CDiscAdjMultizoneDriver::Preprocess(unsigned long TimeIter) {

  const bool time_domain = driver_config->GetTime_Domain();

  for (iZone = 0; iZone < nZone; iZone++) {
    /*--- Set current time iteration ---*/
    config_container[iZone]->SetTimeIter(TimeIter);

    if (time_domain)
      config_container[iZone]->SetPhysicalTime(static_cast<su2double>(TimeIter)*config_container[iZone]->GetDelta_UnstTimeND());
    else
      config_container[iZone]->SetPhysicalTime(0.0);

    /*--- Preprocess the iteration of each zone. ---*/

    iteration_container[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container,
                                                   solver_container, numerics_container, config_container, surface_movement,
                                                   grid_movement, FFDBox, iZone, INST_0);
  }

  if (TimeIter) {
    /*--- Reset cross-terms before new time iterations. ---*/
    for (auto& matOfMat : Cross_Terms)
      for (auto& vecOfMat : matOfMat)
        for (auto& mat : vecOfMat)
          mat = 0.0;

    /*--- Initialize external with dynamic contributions. ---*/
    SetExternalToDualTimeDer();
  }

}

void CDiscAdjMultizoneDriver::StartSolver() {

  const bool time_domain = driver_config->GetTime_Domain();

  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/

  if (rank == MASTER_NODE) {
    cout <<"\n------------------------------ Begin Solver -----------------------------" << endl;

    cout << "\nSimulation Run using the Discrete Adjoint Multizone Driver" << endl;

    if (time_domain)
      cout << "The simulation will run for " << driver_config->GetnTime_Iter() << " time steps." << endl;
  }

  /*--- General setup. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    wrt_sol_freq = min(wrt_sol_freq, config_container[iZone]->GetVolumeOutputFrequency(0));

    /*--- Set BGS_Solution_k to Solution, this is needed to restart
     * correctly as the first OF gradient will overwrite the solution. ---*/

    Set_BGSSolution_k_To_Solution(iZone);

    /*--- Prepare Krylov or quasi-Newton methods. ---*/

    const auto nPoint = geometry_container[iZone][INST_0][MESH_0]->GetnPoint();
    const auto nPointDomain = geometry_container[iZone][INST_0][MESH_0]->GetnPointDomain();
    const auto nVar = GetTotalNumberOfVariables(iZone, true);

    if (config_container[iZone]->GetNewtonKrylov() &&
        config_container[iZone]->GetnQuasiNewtonSamples() >= KrylovMinIters) {
      AdjRHS[iZone].Initialize(nPoint, nPointDomain, nVar, nullptr);
      AdjSol[iZone].Initialize(nPoint, nPointDomain, nVar, nullptr);
      LinSolver[iZone].SetToleranceType(LinearToleranceType::RELATIVE);
    }
    else if (config_container[iZone]->GetnQuasiNewtonSamples() > 1) {
      FixPtCorrector[iZone].resize(config_container[iZone]->GetnQuasiNewtonSamples(), nPoint, nVar, nPointDomain);
    }
  }

  /*--- Size and initialize the matrix of cross-terms. ---*/

  InitializeCrossTerms();

  /*--- Run time iterations. ---*/

  while (TimeIter < driver_config->GetnTime_Iter()) {

    Preprocess(TimeIter);

    Run();

    TimeIter++;
  }

}

bool CDiscAdjMultizoneDriver::Iterate(unsigned short iZone, unsigned long iInnerIter, bool KrylovMode) {

  config_container[iZone]->SetInnerIter(iInnerIter);

  /*--- Evaluate the tape section belonging to solvers in iZone.
   *    Only evaluate TRANSFER terms on the last iteration or after convergence. ---*/

  eval_transfer = (eval_transfer || (iInnerIter == nInnerIter[iZone]-1)) && !KrylovMode;

  ComputeAdjoints(iZone, eval_transfer);

  /*--- Extracting adjoints for solvers in iZone w.r.t. to outputs in iZone (diagonal part). ---*/

  iteration_container[iZone][INST_0]->IterateDiscAdj(geometry_container, solver_container,
                                                     config_container, iZone, INST_0, false);

  /*--- Use QN driver to improve the solution. ---*/

  if (FixPtCorrector[iZone].size()) {
    GetAllSolutions(iZone, true, FixPtCorrector[iZone].FPresult());
    FixPtCorrector[iZone].compute();
    if(iInnerIter) SetAllSolutions(iZone, true, FixPtCorrector[iZone]);
  }

  /*--- This is done explicitly here for multizone cases, only in inner iterations and not when
   *    extracting cross terms so that the adjoint residuals in each zone still make sense. ---*/

  if (!KrylovMode) SetSolutionOldToSolution(iZone);

  /*--- Print out the convergence data to screen and history file. ---*/

  return iteration_container[iZone][INST_0]->Monitor(output_container[iZone], integration_container, geometry_container,
                                                     solver_container, numerics_container, config_container,
                                                     surface_movement, grid_movement, FFDBox, iZone, INST_0);
}

void CDiscAdjMultizoneDriver::KrylovInnerIters(unsigned short iZone) {

  /*--- Use FGMRES to solve the adjoint system, the RHS is -External,
   * the solution are the iZone adjoint variables + External,
   * Recall that External also contains the OF gradient. ---*/

  GetAdjointRHS(iZone, AdjRHS[iZone]);

  AddExternalToSolution(iZone);

  GetAllSolutions(iZone, true, AdjSol[iZone]);

  const bool monitor = config_container[iZone]->GetWrt_ZoneConv();
  const auto product = AdjointProduct(this, iZone);

  /*--- Manipulate the screen output frequency to avoid printing garbage. ---*/
  const auto wrtFreq = config_container[iZone]->GetScreen_Wrt_Freq(2);
  config_container[iZone]->SetScreen_Wrt_Freq(2, nInnerIter[iZone]);
  LinSolver[iZone].SetMonitoringFrequency(wrtFreq);

  Scalar eps = 1.0;
  for (auto totalIter = nInnerIter[iZone]; totalIter >= KrylovMinIters && eps > KrylovTol;) {
    Scalar eps_l = 0.0;
    Scalar tol_l = KrylovTol / eps;
    auto iter = min(totalIter-2ul, config_container[iZone]->GetnQuasiNewtonSamples()-2ul);
    iter = LinSolver[iZone].FGMRES_LinSolver(AdjRHS[iZone], AdjSol[iZone], product, Identity(),
                                             tol_l, iter, eps_l, monitor, config_container[iZone]);
    totalIter -= iter+1;
    eps *= eps_l;
  }

  /*--- Store the solution and restore user settings. ---*/
  SetAllSolutions(iZone, true, AdjSol[iZone]);
  config_container[iZone]->SetScreen_Wrt_Freq(2, wrtFreq);

  /*--- Set the old solution such that iterating gives meaningful residuals. ---*/
  AdjSol[iZone] += AdjRHS[iZone];
  SetAllSolutionsOld(iZone, true, AdjSol[iZone]);

  /*--- Iterate to evaluate cross terms and residuals, this cannot happen within GMRES
   * because the vectors it multiplies by the Jacobian are not the actual solution. ---*/
  eval_transfer = true;
  Iterate(iZone, product.iInnerIter);

  /*--- Set the solution as obtained from GMRES, otherwise it would be GMRES+Iterate once.
   * This is set without the "External" (by adding RHS above) so that it can be added
   * again in the next outer iteration with new contributions from other zones. ---*/
  SetAllSolutions(iZone, true, AdjSol[iZone]);

}

void CDiscAdjMultizoneDriver::Run() {

  const unsigned long nOuterIter = driver_config->GetnOuter_Iter();
  const bool time_domain = driver_config->GetTime_Domain();
  driver_config->Set_StartTime(SU2_MPI::Wtime());

  /*--- If the gradient of the objective function is 0 so are the adjoint variables.
   * Unless in unsteady problems where there are other contributions to the RHS. ---*/

  const auto zeroGrad = EvaluateObjectiveFunctionGradient();

  if (zeroGrad && !time_domain) {
    if (rank == MASTER_NODE) {
      cout << "\nThe gradient of the objective function is numerically 0.";
      cout << "\nThis implies that the adjoint variables are also 0.\n\n";
    }
    EvaluateSensitivities(0, true);
    return;
  }

  /*--- Loop over the number of outer iterations. ---*/

  for (unsigned long iOuterIter = 0, StopCalc = false; !StopCalc; iOuterIter++) {

    driver_config->SetOuterIter(iOuterIter);

    for (iZone = 0; iZone < nZone; iZone++)
      config_container[iZone]->SetOuterIter(iOuterIter);

    /*--- For the adjoint iteration we need the derivatives of the iteration function with
     *    respect to the state (and possibly the mesh coordinate) variables.
     *    Since these derivatives do not change in the steady state case we only have to record
     *    if the current recording is different from them.
     *
     *    To set the tape appropriately, the following recording methods are provided:
     *    (1) CLEAR_INDICES: All information from a previous recording is removed.
     *    (2) SOLUTION_VARIABLES: State variables of all solvers in a zone as input.
     *    (3) MESH_COORDS / MESH_DEFORM: Mesh coordinates as input.
     *    (4) SOLUTION_AND_MESH: Mesh coordinates and state variables as input.
     *
     *    By default, all (state and mesh coordinate variables) will be declared as output,
     *    since it does not change the computational effort, just the memory consumption of the tape. ---*/


    /*--- If we want to set up zone-specific tapes (retape), we do not need to record
     *    here. Otherwise, the whole tape of a coupled run will be created. ---*/

    if (RecordingState != RECORDING::SOLUTION_VARIABLES) {
      SetRecording(RECORDING::CLEAR_INDICES, Kind_Tape::FULL_TAPE, ZONE_0);
      SetRecording(RECORDING::SOLUTION_VARIABLES, Kind_Tape::FULL_TAPE, ZONE_0);
    }

    /*-- Start loop over zones. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {

      config_container[iZone]->Set_StartTime(SU2_MPI::Wtime());

      /*--- Start inner iterations from where we stopped in previous outer iteration. ---*/

      Set_Solution_To_BGSSolution_k(iZone);

      /*--- Inner loop to allow for multiple adjoint updates with respect to solvers in iZone. ---*/

      const bool restart = config_container[iZone]->GetRestart();
      const bool no_restart = (iOuterIter > 0) || !restart;
      eval_transfer = false;

      /*--- Reset QN driver for new inner iterations. ---*/

      if (FixPtCorrector[iZone].size()) {
        FixPtCorrector[iZone].reset();
        if(restart && (iOuterIter==1)) GetAllSolutions(iZone, true, FixPtCorrector[iZone]);
      }

      if (!config_container[iZone]->GetNewtonKrylov() || !no_restart || nInnerIter[iZone] < KrylovMinIters) {

        /*--- Regular fixed-point, possibly with quasi-Newton method. ---*/

        for (unsigned long iInnerIter = 0; iInnerIter < nInnerIter[iZone]; iInnerIter++) {

          /*--- Add off-diagonal contribution (including the OF gradient) to Solution. ---*/

          if (no_restart || (iInnerIter > 0)) {
            AddExternalToSolution(iZone);
          }
          else {
            /*--- If we restarted, Solution already has all contributions,
             *    we run only one inner iter to compute the cross terms. ---*/
            eval_transfer = true;
          }

          const bool converged = Iterate(iZone, iInnerIter);

          if (eval_transfer) break;

          eval_transfer = converged;
        }
      }
      else {
        /*--- Use Krylov method to drive inner iterations of this zone. ---*/

        KrylovInnerIters(iZone);
      }

      /*--- Off-diagonal (coupling term) BGS update. ---*/

      for (unsigned short jZone = 0; jZone < nZone; jZone++) {

        if (jZone != iZone && interface_container[jZone][iZone] != nullptr) {

          /*--- Extracting adjoints for solvers in jZone w.r.t. to the output of all solvers in iZone,
           *    that is, for the cases iZone != jZone we are evaluating cross derivatives between zones. ---*/

          config_container[jZone]->SetInnerIter(0);
          iteration_container[jZone][INST_0]->IterateDiscAdj(geometry_container, solver_container,
                                                             config_container, jZone, INST_0, true);

          /*--- Extract the cross-term performing a relaxed update of it and of the sum (External) for jZone. ---*/

          UpdateCrossTerm(iZone, jZone);
        }
      }

      /*--- Compute residual from Solution and Solution_BGS_k and update the latter. ---*/

      SetResidual_BGS(iZone);
      Set_BGSSolution_k_To_Solution(iZone);
    }

    /*--- Set the multizone output. ---*/

    driver_output->SetMultizoneHistoryOutput(output_container, config_container, driver_config, TimeIter, iOuterIter);

    /*--- Check for convergence. ---*/

    StopCalc = driver_output->GetConvergence() || (iOuterIter == nOuterIter-1) ||
               ((nZone==1) && output_container[ZONE_0]->GetConvergence());

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    /*--- Compute the geometrical sensitivities and write them to file, except for time_domain. ---*/

    if (time_domain) continue;

    if (StopCalc || ((iOuterIter % wrt_sol_freq == 0) && (iOuterIter != 0)))
      EvaluateSensitivities(iOuterIter, StopCalc);
  }

  if (time_domain) {
    EvaluateSensitivities(TimeIter, (TimeIter+1) == driver_config->GetnTime_Iter());
  }

}

bool CDiscAdjMultizoneDriver::EvaluateObjectiveFunctionGradient() {

  /*--- Evaluate the objective function gradient w.r.t. the solutions of all zones. ---*/

  SetRecording(RECORDING::CLEAR_INDICES, Kind_Tape::OBJECTIVE_FUNCTION_TAPE, ZONE_0);
  SetRecording(RECORDING::SOLUTION_VARIABLES, Kind_Tape::OBJECTIVE_FUNCTION_TAPE, ZONE_0);
  RecordingState = RECORDING::CLEAR_INDICES;

  AD::ClearAdjoints();
  SetAdjObjFunction();
  AD::ComputeAdjoint(OBJECTIVE_FUNCTION, START);

  /*--- Initialize External with the objective function gradient. ---*/

  su2double rhs_norm = 0.0;

  for (iZone = 0; iZone < nZone; iZone++) {

    iteration_container[iZone][INST_0]->IterateDiscAdj(geometry_container, solver_container,
                                                       config_container, iZone, INST_0, false);
    AddSolutionToExternal(iZone);

    for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
      auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
      if (solver && solver->GetAdjoint())
        for (unsigned short iVar=0; iVar < solver->GetnVar(); ++iVar)
          rhs_norm += solver->GetRes_RMS(iVar);
    }
  }

  return rhs_norm < EPS;
}

void CDiscAdjMultizoneDriver::EvaluateSensitivities(unsigned long Iter, bool force_writing) {

  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
   *    as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(RECORDING::CLEAR_INDICES, Kind_Tape::FULL_TAPE, ZONE_0);

  /*--- Store the computational graph of one direct iteration with the mesh coordinates as input. ---*/

  SetRecording(RECORDING::MESH_COORDS, Kind_Tape::FULL_TAPE, ZONE_0);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the current iteration. The values are passed to the AD tool. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    Set_Solution_To_BGSSolution_k(iZone);

    AddExternalToSolution(iZone);

    iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container,
                                                          config_container, iZone, INST_0);
  }

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdjObjFunction();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed sensitivity values. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    auto config = config_container[iZone];
    auto solvers = solver_container[iZone][INST_0][MESH_0];
    auto geometry = geometry_container[iZone][INST_0][MESH_0];

    int IDX_SOL = -1;
    if (config->GetFluidProblem()) IDX_SOL = ADJFLOW_SOL;
    else if (config->GetHeatProblem()) IDX_SOL = ADJHEAT_SOL;
    else if (config->GetStructuralProblem()) IDX_SOL = ADJFEA_SOL;
    else {
      if (rank == MASTER_NODE)
        cout << "WARNING: Sensitivities not set for one of the specified discrete adjoint solvers!" << endl;
      continue;
    }

    if (Has_Deformation(iZone))
      solvers[ADJMESH_SOL]->SetSensitivity(geometry, config, solvers[IDX_SOL]);
    else
      solvers[IDX_SOL]->SetSensitivity(geometry, config);

    iteration_container[iZone][INST_0]->Postprocess(output_container[iZone], integration_container, geometry_container,
                                                    solver_container, numerics_container, config_container,
                                                    surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  /*--- Output files (CMultizoneDriver::Output uses StopCalc to force file output). ---*/

  swap(StopCalc, force_writing);
  Output(Iter);
  swap(StopCalc, force_writing);

}

void CDiscAdjMultizoneDriver::SetRecording(RECORDING kind_recording, Kind_Tape tape_type, unsigned short record_zone) {

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution. ---*/

  for(iZone = 0; iZone < nZone; iZone++) {
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++) {
      for (unsigned short iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
        auto solver = solver_container[iZone][INST_0][iMesh][iSol];
        if (solver && solver->GetAdjoint()) {
          solver->SetRecording(geometry_container[iZone][INST_0][iMesh], config_container[iZone]);
        }
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "\n-------------------------------------------------------------------------\n";
    switch(kind_recording) {
    case RECORDING::CLEAR_INDICES:      cout << "Clearing the computational graph." << endl; break;
    case RECORDING::MESH_COORDS:        cout << "Storing computational graph wrt MESH COORDINATES." << endl; break;
    case RECORDING::SOLUTION_VARIABLES: cout << "Storing computational graph wrt CONSERVATIVE VARIABLES." << endl; break;
    default: break;
    }
  }

  /*--- Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if(kind_recording != RECORDING::CLEAR_INDICES) {

    AD::StartRecording();

    AD::Push_TapePosition(); /// START

    for (iZone = 0; iZone < nZone; iZone++) {

      /*--- In multi-physics, MESH_COORDS is an umbrella term for "geometric sensitivities",
       *    if a zone has mesh deformation its recording type needs to change to MESH_DEFORM
       *    as those sensitivities are managed by the adjoint mesh solver instead. ---*/

      RECORDING type_recording = kind_recording;

      if (Has_Deformation(iZone) && (kind_recording == RECORDING::MESH_COORDS)) {
        type_recording = RECORDING::MESH_DEFORM;
      }

      iteration_container[iZone][INST_0]->RegisterInput(solver_container, geometry_container,
                                                        config_container, iZone, INST_0, type_recording);
    }
  }

  AD::Push_TapePosition(); /// REGISTERED

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->SetDependencies(solver_container, geometry_container, numerics_container,
                                                        config_container, iZone, INST_0, kind_recording);
  }

  AD::Push_TapePosition(); /// DEPENDENCIES

  /*--- Extract the objective function and store it.
   *    It is necessary to include data transfer and mesh updates in this section as some functions
   *    computed in one zone depend explicitly on the variables of others through that path. --- */

  if ((tape_type == Kind_Tape::OBJECTIVE_FUNCTION_TAPE) || (kind_recording == RECORDING::MESH_COORDS)) {
    HandleDataTransfer();
    for (iZone = 0; iZone < nZone; iZone++) {
      if (Has_Deformation(iZone)) {
        iteration_container[iZone][INST_0]->SetDependencies(solver_container, geometry_container, numerics_container,
                                                            config_container, iZone, INST_0, kind_recording);
      }
    }
    SetObjFunction(kind_recording);
  }

  AD::Push_TapePosition(); /// OBJECTIVE_FUNCTION

  if (tape_type != Kind_Tape::OBJECTIVE_FUNCTION_TAPE) {

    /*--- We do the communication here to not differentiate wrt updated boundary data.
     *    For recording w.r.t. mesh coordinates the transfer was included before the
     *    objective function, so we do not repeat it here. ---*/

    if (kind_recording != RECORDING::MESH_COORDS) {
      HandleDataTransfer();
    }

    AD::Push_TapePosition(); /// TRANSFER

    for(iZone = 0; iZone < nZone; iZone++) {

      AD::Push_TapePosition(); /// enter_zone

      DirectIteration(iZone, kind_recording);

      iteration_container[iZone][INST_0]->RegisterOutput(solver_container, geometry_container,
                                                         config_container, iZone, INST_0);
      AD::Push_TapePosition(); /// leave_zone
    }
    PrintDirectResidual(kind_recording);
  }

  if (kind_recording != RECORDING::CLEAR_INDICES && driver_config->GetWrt_AD_Statistics()) {
    if (rank == MASTER_NODE) AD::PrintStatistics();
#ifdef CODI_REVERSE_TYPE
    if (size > SINGLE_NODE) {
      su2double myMem = AD::getTape().getTapeValues().getUsedMemorySize(), totMem = 0.0;
      SU2_MPI::Allreduce(&myMem, &totMem, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      if (rank == MASTER_NODE) {
        cout << "MPI\n";
        cout << "-------------------------------------\n";
        cout << "  Total memory used      :  " << totMem << " MB\n";
        cout << "-------------------------------------\n" << endl;
      }
    }
#endif
  }

  AD::StopRecording();

  RecordingState = kind_recording;
}

void CDiscAdjMultizoneDriver::DirectIteration(unsigned short iZone, RECORDING kind_recording) {

  /*--- Do one iteration of the direct solver ---*/
  direct_iteration[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container,
                                              solver_container, numerics_container, config_container,
                                              surface_movement, grid_movement, FFDBox, iZone, INST_0);

  /*--- Iterate the zone as a block a single time ---*/
  direct_iteration[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container,
                                           solver_container, numerics_container, config_container,
                                           surface_movement, grid_movement, FFDBox, iZone, INST_0);

}

void CDiscAdjMultizoneDriver::SetObjFunction(RECORDING kind_recording) {

  ObjFunc = 0.0;

  /*--- Call objective function calculations. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    auto config = config_container[iZone];
    auto solvers = solver_container[iZone][INST_0][MESH_0];
    auto geometry = geometry_container[iZone][INST_0][MESH_0];

    switch (config->GetKind_Solver()) {

      case MAIN_SOLVER::DISC_ADJ_EULER:     case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:     case MAIN_SOLVER::DISC_ADJ_RANS:
      case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:

        solvers[FLOW_SOL]->Pressure_Forces(geometry, config);
        solvers[FLOW_SOL]->Momentum_Forces(geometry, config);
        solvers[FLOW_SOL]->Friction_Forces(geometry, config);

        if(config->GetWeakly_Coupled_Heat()) {
          solvers[HEAT_SOL]->Heat_Fluxes(geometry, solvers, config);
        }

        direct_output[iZone]->SetHistoryOutput(geometry, solvers, config);
        ObjFunc += solvers[FLOW_SOL]->GetTotal_ComboObj();
        break;

      case MAIN_SOLVER::DISC_ADJ_HEAT:
        solvers[HEAT_SOL]->Heat_Fluxes(geometry, solvers, config);
        direct_output[iZone]->SetHistoryOutput(geometry, solvers, config);
        ObjFunc += solvers[HEAT_SOL]->GetTotal_ComboObj();
        break;

      case MAIN_SOLVER::DISC_ADJ_FEM:
        solvers[FEA_SOL]->Postprocessing(geometry, config, numerics_container[iZone][INST_0][MESH_0][FEA_SOL], true);
        direct_output[iZone]->SetHistoryOutput(geometry, solvers, config);
        ObjFunc += solvers[FEA_SOL]->GetTotal_ComboObj();
        break;

      default:
        break;
    }
  }

  if (rank == MASTER_NODE) {
    AD::RegisterOutput(ObjFunc);
    AD::SetIndex(ObjFunc_Index, ObjFunc);
    if (kind_recording == RECORDING::SOLUTION_VARIABLES) {
      cout << " Objective function                   : " << ObjFunc;
      if (driver_config->GetWrt_AD_Statistics()){
        cout << " (" << ObjFunc_Index << ")\n";
      }
      cout << endl;
    }
  }
}

void CDiscAdjMultizoneDriver::SetAdjObjFunction() {
  su2double seeding = 1.0;

  if (config_container[ZONE_0]->GetTime_Domain()) {
    const auto IterAvg_Obj = config_container[ZONE_0]->GetIter_Avg_Objective();
    if (TimeIter < IterAvg_Obj) {
      /*--- Default behavior when no window is chosen is to use Square-Windowing, i.e. the numerator equals 1.0 ---*/
      auto windowEvaluator = CWindowingTools();
      const su2double weight =
          windowEvaluator.GetWndWeight(config_container[ZONE_0]->GetKindWindow(), TimeIter, IterAvg_Obj - 1);
      seeding = weight / IterAvg_Obj;
    }
    else {
      seeding = 0.0;
    }
  }
  if (rank == MASTER_NODE) {
    AD::ResizeAdjoints();
    AD::SetDerivative(ObjFunc_Index, SU2_TYPE::GetValue(seeding));
  }
}

void CDiscAdjMultizoneDriver::ComputeAdjoints(unsigned short iZone, bool eval_transfer) {

#if defined(CODI_INDEX_REUSE)
  if (nZone > 1 && rank == MASTER_NODE) {
    std::cout << "WARNING: AD types that reuse indices do not support multiple zones." << std::endl;
  }
#endif

  AD::ClearAdjoints();

  /*--- Initialize the adjoints in iZone ---*/

  iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container,
                                                        config_container, iZone, INST_0);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  const unsigned short enter_izone = iZone*2+1 + ITERATION_READY;
  const unsigned short leave_izone = iZone*2 + ITERATION_READY;

  AD::ComputeAdjoint(enter_izone, leave_izone);

  /*--- Compute adjoints of transfer and mesh deformation routines, only strictly needed
   *    on the last inner iteration. Structural problems have some minor issue and we
   *    need to evaluate this section on every iteration. ---*/

  if (eval_transfer || config_container[iZone]->GetStructuralProblem())
    AD::ComputeAdjoint(TRANSFER, OBJECTIVE_FUNCTION);

  /*--- Adjoints of dependencies, needed if derivatives of variables
   *    are extracted (e.g. AoA, Mach, etc.) ---*/

  AD::ComputeAdjoint(DEPENDENCIES, START);

}

void CDiscAdjMultizoneDriver::InitializeCrossTerms() {

  if (Cross_Terms.empty()) {
    Cross_Terms.resize(nZone, vector<vector<su2passivematrix> >(nZone));
  }

  for(auto iZone = 0u; iZone < nZone; iZone++) {
    for (auto jZone = 0u; jZone < nZone; jZone++) {
      if (iZone != jZone || interface_container[jZone][iZone] != nullptr) {

        /*--- If jZone contributes to iZone in the primal problem, then
         *    iZone contributes to jZone in the adjoint problem. ---*/

        Cross_Terms[iZone][jZone].resize(MAX_SOLS);

        for (auto iSol = 0u; iSol < MAX_SOLS; iSol++) {

          auto solver = solver_container[jZone][INST_0][MESH_0][iSol];
          if (solver && solver->GetAdjoint()) {

            const auto nPoint = geometry_container[jZone][INST_0][MESH_0]->GetnPoint();
            const auto nVar = solver->GetnVar();
            Cross_Terms[iZone][jZone][iSol].resize(nPoint,nVar) = 0.0;
          }
        }
      }
    }
  }
}

void CDiscAdjMultizoneDriver::HandleDataTransfer() {

  for(iZone = 0; iZone < nZone; iZone++) {

    /*--- In principle, the mesh does not need to be updated ---*/
    bool DeformMesh = false;

    /*--- Transfer from all the remaining zones ---*/
    for (unsigned short jZone = 0; jZone < nZone; jZone++){
      /*--- The target zone is iZone ---*/
      if (jZone != iZone && interface_container[jZone][iZone] != nullptr) {
        DeformMesh |= TransferData(jZone, iZone);
      }
    }

    /*--- If a mesh update is required due to the transfer of data ---*/
    const unsigned long ExtIter = 0;
    if (DeformMesh) DynamicMeshUpdate(iZone, ExtIter);

    Has_Deformation(iZone) = DeformMesh;
  }
}

void CDiscAdjMultizoneDriver::AddSolutionToExternal(unsigned short iZone) {

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
    if (solver && solver->GetAdjoint())
      solver->Add_Solution_To_External();
  }
}

void CDiscAdjMultizoneDriver::SetExternalToDualTimeDer() {

  for (unsigned short iZone = 0; iZone < nZone; iZone++) {
    for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
      auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
      if (solver && solver->GetAdjoint())
        solver->GetNodes()->Set_External_To_DualTimeDer();
    }
  }
}

void CDiscAdjMultizoneDriver::AddExternalToSolution(unsigned short iZone) {

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
    if (solver && solver->GetAdjoint())
      solver->Add_External_To_Solution();
  }
}

void CDiscAdjMultizoneDriver::SetSolutionOldToSolution(unsigned short iZone) {

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
    if (solver && solver->GetAdjoint())
      solver->Set_OldSolution();
  }
}

void CDiscAdjMultizoneDriver::UpdateCrossTerm(unsigned short iZone, unsigned short jZone) {

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    auto solver = solver_container[jZone][INST_0][MESH_0][iSol];
    if (solver && solver->GetAdjoint())
      solver->Update_Cross_Term(config_container[jZone], Cross_Terms[iZone][jZone][iSol]);
  }
}

void CDiscAdjMultizoneDriver::Set_Solution_To_BGSSolution_k(unsigned short iZone) {

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
    if (solver && solver->GetAdjoint())
      solver->GetNodes()->Restore_BGSSolution_k();
  }
}

void CDiscAdjMultizoneDriver::Set_BGSSolution_k_To_Solution(unsigned short iZone) {

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
    if (solver && solver->GetAdjoint())
      solver->GetNodes()->Set_BGSSolution_k();
  }
}

void CDiscAdjMultizoneDriver::SetResidual_BGS(unsigned short iZone) {

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
    if (solver && solver->GetAdjoint())
      solver->ComputeResidual_Multizone(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
  }
}
