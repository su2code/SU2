/*!
 * \file driver_adjoint_singlezone.cpp
 * \brief The main subroutines for driving adjoint single-zone problems.
 * \author R. Sanchez, H. Patel, A. Gastaldi
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/drivers/CDiscAdjSinglezoneDriver.hpp"
#include "../../include/output/tools/CWindowingTools.hpp"
#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIterationFactory.hpp"
#include "../../include/iteration/CTurboIteration.hpp"
#include "../../../Common/include/toolboxes/CQuasiNewtonInvLeastSquares.hpp"

CDiscAdjSinglezoneDriver::CDiscAdjSinglezoneDriver(char* confFile,
                                                   unsigned short val_nZone,
                                                   SU2_Comm MPICommunicator) : CSinglezoneDriver(confFile,
                                                                                                 val_nZone,
                                                                                                 MPICommunicator) {


  /*--- Store the number of internal iterations that will be run by the adjoint solver ---*/
  nAdjoint_Iter = config_container[ZONE_0]->GetnInner_Iter();

  /*--- Store the pointers ---*/
  config      = config_container[ZONE_0];
  iteration   = iteration_container[ZONE_0][INST_0];
  solver      = solver_container[ZONE_0][INST_0][MESH_0];
  numerics    = numerics_container[ZONE_0][INST_0][MESH_0];
  geometry    = geometry_container[ZONE_0][INST_0][MESH_0];
  integration = integration_container[ZONE_0][INST_0];

  /*--- Store the recording state ---*/
  RecordingState = RECORDING::CLEAR_INDICES;

  /*--- Initialize the direct iteration ---*/

  switch (config->GetKind_Solver()) {

  case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
  case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: Euler/Navier-Stokes/RANS equation." << endl;

    if (config->GetBoolTurbomachinery()) {
      direct_iteration = new CTurboIteration(config);
    }
    else { direct_iteration = CIterationFactory::CreateIteration(MAIN_SOLVER::EULER, config); }

    if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {
      direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::EULER, config, nDim);
    }
    else { direct_output =  COutputFactory::CreateOutput(MAIN_SOLVER::INC_EULER, config, nDim); }

    MainVariables = RECORDING::SOLUTION_VARIABLES;
    if (config->GetDeform_Mesh()) {
      SecondaryVariables = RECORDING::MESH_DEFORM;
    } else {
      SecondaryVariables = RECORDING::MESH_COORDS;
    }
    MainSolver = ADJFLOW_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_FEM_EULER : case MAIN_SOLVER::DISC_ADJ_FEM_NS : case MAIN_SOLVER::DISC_ADJ_FEM_RANS :
    if (rank == MASTER_NODE)
      cout << "Direct iteration: Euler/Navier-Stokes/RANS equation." << endl;
    direct_iteration = CIterationFactory::CreateIteration(MAIN_SOLVER::FEM_EULER, config);
    direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::FEM_EULER, config, nDim);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJFLOW_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_FEM:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: elasticity equation." << endl;
    direct_iteration =  CIterationFactory::CreateIteration(MAIN_SOLVER::FEM_ELASTICITY, config);
    direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::FEM_ELASTICITY, config, nDim);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJFEA_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_HEAT:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: heat equation." << endl;
    direct_iteration = CIterationFactory::CreateIteration(MAIN_SOLVER::HEAT_EQUATION, config);
    direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::HEAT_EQUATION, config, nDim);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJHEAT_SOL;
    break;

  default:
    break;

  }

 direct_output->PreprocessHistoryOutput(config, false);

}

CDiscAdjSinglezoneDriver::~CDiscAdjSinglezoneDriver() {

  delete direct_iteration;
  delete direct_output;
  delete PrimalJacobian;
  delete PrimalPreconditioner;
  delete AdjOperator;
  delete AdjPreconditioner;

}

void CDiscAdjSinglezoneDriver::Preprocess(unsigned long TimeIter) {

  /*--- Set the current time iteration in the config and also in the driver
   * because the python interface doesn't offer an explicit way of doing it. ---*/

  this->TimeIter = TimeIter;
  config_container[ZONE_0]->SetTimeIter(TimeIter);

  /*--- Preprocess the adjoint iteration ---*/

  iteration->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
                        solver_container, numerics_container, config_container,
                        surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- For the adjoint iteration we need the derivatives of the iteration function with
   *--- respect to the conservative variables. Since these derivatives do not change in the steady state case
   *--- we only have to record if the current recording is different from the main variables. ---*/

  if (RecordingState != MainVariables){
    MainRecording();
  }

}

void CDiscAdjSinglezoneDriver::Run() {
  if (config->GetKind_DiscreteAdjoint() == ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    RunResidual();
  } else {
    RunFixedPoint();
  }
}

void CDiscAdjSinglezoneDriver::RunFixedPoint() {

  CQuasiNewtonInvLeastSquares<passivedouble> fixPtCorrector;
  if (config->GetnQuasiNewtonSamples() > 1) {
    fixPtCorrector.resize(config->GetnQuasiNewtonSamples(),
                          geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(),
                          GetTotalNumberOfVariables(ZONE_0,true),
                          geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain());

    if (TimeIter != 0) GetAllSolutions(ZONE_0, true, fixPtCorrector);
  }

  for (auto Adjoint_Iter = 0ul; Adjoint_Iter < nAdjoint_Iter; Adjoint_Iter++) {

    /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
     *--- of the previous iteration. The values are passed to the AD tool.
     *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/

    config->SetInnerIter(Adjoint_Iter);

    /*--- Update the state of the tape. ---*/

    UpdateAdjointsFixedPoint();

    /*--- Output files for steady state simulations. ---*/

    if (!config->GetTime_Domain()) {
      iteration->Output(output_container[ZONE_0], geometry_container, solver_container,
                        config_container, Adjoint_Iter, false, ZONE_0, INST_0);
    }

    if (StopCalc) break;

    /*--- Correct the solution with the quasi-Newton approach. ---*/

    if (fixPtCorrector.size()) {
      GetAllSolutions(ZONE_0, true, fixPtCorrector.FPresult());
      SetAllSolutions(ZONE_0, true, fixPtCorrector.compute());
    }
  }
}

void CDiscAdjSinglezoneDriver::RunResidual() {
  if (!KrylovSet) {
    /*--- Initialize the solution, right-hand-side, and system. ---*/
    const auto nVar = GetTotalNumberOfVariables(ZONE_0, true);
    const auto nPoint = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    const auto nPointDomain = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain();

    AdjRHS.Initialize(nPoint, nPointDomain, nVar, nullptr);
    AdjSol.Initialize(nPoint, nPointDomain, nVar, nullptr);

    AdjSolver.SetToleranceType(LinearToleranceType::RELATIVE);
    KrylovSet = true;

    /*--- Initialize the preconditioner using the (transpose) approximate Jacobian from the primal problem. ---*/

    if (config->GetKind_TimeIntScheme() != EULER_IMPLICIT) {
      std::cout << "Cannot build a preconditioner for the discrete-adjoint system (missing primal Jacobian structure) !"
                << std::endl;

      return;
    };

    UpdateJacobians();

    CopiedJacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      for (unsigned long jPoint = 0; jPoint < nPoint; jPoint++) {
        auto value = solver[FLOW_SOL]->Jacobian.GetBlock(iPoint, jPoint);

        CopiedJacobian.SetBlock(iPoint, jPoint, value);
      }
    }

    CopiedJacobian.TransposeInPlace();
    PrimalJacobian = new CSysMatrixVectorProduct<Scalar>(CopiedJacobian, geometry, config);

    const auto kindPreconditioner = static_cast<ENUM_LINEAR_SOLVER_PREC>(config->GetKind_Linear_Solver_Prec());
    PrimalPreconditioner = CPreconditioner<Scalar>::Create(kindPreconditioner, CopiedJacobian, geometry, config);
    PrimalPreconditioner->Build();
  }

  /*--- Use FGMRES to solve the adjoint system, where:
   *      * the RHS is -dObjective/dStates (and any external contributions),
   *      * the solution are the adjoint variables, and
   *      * the system applies the matrix-vector product with dResidual/dStates.  ---*/
  UpdateAdjointsResidual();

  GetAllSolutions(ZONE_0, true, AdjSol);
  GetAllObjectiveStatesSensitivities(AdjRHS);

  /*--- Manipulate the screen output frequency to avoid printing garbage. ---*/
  const bool monitor = true;
  const auto wrtFreq = 1;

  AdjSolver.SetMonitoringFrequency(wrtFreq);

  /*--- Initialize the linear solver iterations ---*/
  AdjOperator = new LinOperator(this);
  AdjPreconditioner = new LinPreconditioner(this);

  Scalar eps = 1.0;

  unsigned long nKrylov_Iter;
  for (nKrylov_Iter = nAdjoint_Iter; nKrylov_Iter >= KrylovMinIters && eps > KrylovSysTol;) {
    std::cout << "Adjoint iteration: " << nKrylov_Iter << " ... " << std::endl;

    auto nIter = min(nKrylov_Iter - 2ul, config_container[iZone]->GetnQuasiNewtonSamples() - 2ul);
    Scalar eps_l = 0.0;
    Scalar tol_l = KrylovSysTol / eps;

    nIter = AdjSolver.FGMRES_LinSolver(AdjRHS, AdjSol, *AdjOperator, *AdjPreconditioner, tol_l, nIter, eps_l, monitor,
                                       config_container[ZONE_0]);
    nKrylov_Iter -= nIter + 1;

    eps *= eps_l;
  }

  /*--- Store the solution and restore user settings. ---*/
  SetAllSolutions(ZONE_0, true, AdjSol);

  UpdateAdjointsResidual();

  /*--- Apply the solution to obtain the total sensitivities (w.r.t. deformed volume coordinates). ---*/
  Postprocess();

  /*--- HACK ! Force output here until proper convergence monitoring can be implemented. ---*/
  const auto inst = config_container[ZONE_0]->GetiInst();

  for (iInst = 0; iInst < nInst[ZONE_0]; ++iInst) {
    config_container[ZONE_0]->SetiInst(iInst);
    output_container[ZONE_0]->SetResultFiles(geometry_container[ZONE_0][iInst][MESH_0], config_container[ZONE_0],
                                              solver_container[ZONE_0][iInst][MESH_0], nAdjoint_Iter - nKrylov_Iter,
                                              true);
  }
  config_container[ZONE_0]->SetiInst(inst);
}

void CDiscAdjSinglezoneDriver::UpdateAdjoints() {
  if (config->GetKind_DiscreteAdjoint() == ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    UpdateAdjointsResidual();
  } else {
    UpdateAdjointsFixedPoint();
  }
}

void CDiscAdjSinglezoneDriver::UpdateAdjointsFixedPoint() {
  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *--- of the previous iteration. The values are passed to the AD tool.
   *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/
  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/
  SetAdjointObjective();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/
  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  iteration->IterateDiscAdj(geometry_container, solver_container,
                            config_container, ZONE_0, INST_0, false);

  /*--- Monitor the pseudo-time ---*/
  StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                solver_container, numerics_container, config_container,
                                surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/
  AD::ClearAdjoints();

}

void CDiscAdjSinglezoneDriver::UpdateAdjointsResidual() {
  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *--- of the previous iteration. The values are passed to the AD tool.
   *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/
  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/
  AD::ComputeAdjoint();

  /*--- Extract the adjoints of the residuals and store them for the next iteration ---*/
  if (config->GetFluidProblem()) {
    solver[ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry, config, ENUM_VARIABLE::RESIDUALS);
    solver[ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry, config, ENUM_VARIABLE::RESIDUALS);
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/
  AD::ClearAdjoints();

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/
  SetAdjointObjective();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/
  AD::ComputeAdjoint();

  /*--- Extract the adjoints of the objective function and store them for the next iteration ---*/
  if (config->GetFluidProblem()) {
    solver[ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry, config, ENUM_VARIABLE::OBJECTIVE);
    solver[ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry, config, ENUM_VARIABLE::OBJECTIVE);
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/
  AD::ClearAdjoints();

  /*--- Initialize the adjoint of the vertex tractions with the corresponding adjoint vector. ---*/
  solver[FLOW_SOL]->SetVertexTractionsAdjoint(geometry, config);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/
  AD::ComputeAdjoint();

  /*--- Extract the adjoints of the vertex tractions and store them for the next iteration ---*/
  if (config->GetFluidProblem()) {
    solver[ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry, config, ENUM_VARIABLE::TRACTIONS);
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/
  AD::ClearAdjoints();
}

void CDiscAdjSinglezoneDriver::SetAdjointObjective() {
  su2double seeding = 1.0;

  if (config->GetTime_Domain()) {
    const auto IterAvg_Obj = config->GetIter_Avg_Objective();
    if (TimeIter < IterAvg_Obj) {
      /*--- Default behavior when no window is chosen is to use Square-Windowing, i.e. the numerator equals 1.0 ---*/
      auto windowEvaluator = CWindowingTools();
      const su2double weight = windowEvaluator.GetWndWeight(config->GetKindWindow(), TimeIter, IterAvg_Obj - 1);
      seeding = weight / IterAvg_Obj;
    } else {
      seeding = 0.0;
    }
  }
  if (rank == MASTER_NODE) {
    SU2_TYPE::SetDerivative(ObjFunc, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc, 0.0);
  }
}

void CDiscAdjSinglezoneDriver::Postprocess() {

  switch(config->GetKind_Solver())
  {
    case MAIN_SOLVER::DISC_ADJ_EULER :     case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES :     case MAIN_SOLVER::DISC_ADJ_RANS :
    case MAIN_SOLVER::DISC_ADJ_INC_EULER : case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES : case MAIN_SOLVER::DISC_ADJ_INC_RANS :
    case MAIN_SOLVER::DISC_ADJ_HEAT :

      /*--- Compute the geometrical sensitivities ---*/
      SecondaryRecording();

      if (config->GetKind_DiscreteAdjoint() == ENUM_DISC_ADJ_TYPE::RESIDUALS) {
        SecondaryRunResidual();
      } else {
        SecondaryRunFixedPoint();
      }

      break;

    case MAIN_SOLVER::DISC_ADJ_FEM :

      /*--- Compute the geometrical sensitivities ---*/
      SecondaryRecording();

      if (config->GetKind_DiscreteAdjoint() == ENUM_DISC_ADJ_TYPE::RESIDUALS) {
        SecondaryRunResidual();
      } else {
        SecondaryRunFixedPoint();
      }

      iteration->Postprocess(output_container[ZONE_0], integration_container, geometry_container,
                             solver_container, numerics_container, config_container,
                             surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);
      break;

    default:
      break;

  }//switch

}

void CDiscAdjSinglezoneDriver::SetRecording(RECORDING kind_recording){

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution. ---*/

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    for (unsigned short iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
      auto solver = solver_container[ZONE_0][INST_0][iMesh][iSol];
      if (solver && solver->GetAdjoint()) {
        SU2_OMP_PARALLEL_(if (solver->GetHasHybridParallel()))
        solver->SetRecording(geometry_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0]);
        END_SU2_OMP_PARALLEL
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "\n-------------------------------------------------------------------------\n";
    switch(kind_recording) {
    case RECORDING::CLEAR_INDICES: cout << "Clearing the computational graph." << endl; break;
    case RECORDING::MESH_COORDS:   cout << "Storing computational graph wrt MESH COORDINATES." << endl; break;
    case RECORDING::SOLUTION_VARIABLES:
      cout << "Direct iteration to store the primal computational graph." << endl;
      cout << "Computing residuals to check the convergence of the direct problem." << endl; break;
    default: break;
    }
  }

  /*---Enable recording and register input of the iteration --- */

  if (kind_recording != RECORDING::CLEAR_INDICES){

    AD::StartRecording();

    iteration->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);
  }

  /*--- Set the dependencies of the iteration ---*/

  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0,
                             INST_0, kind_recording);

  /*--- Do one iteration of the direct solver. ---*/
  if (config->GetKind_DiscreteAdjoint() == ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    DirectRunResidual(kind_recording);
  } else {
    DirectRunFixedPoint(kind_recording);
  }

  /*--- Store the recording state ---*/

  RecordingState = kind_recording;

  /*--- Register Output of the iteration ---*/

  iteration->RegisterOutput(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Extract the tractions and objective function and store them. --- */
  UpdateTractions();
  UpdateObjective();

  if (rank == MASTER_NODE) {
    AD::RegisterOutput(ObjFunc);
  }

  if (kind_recording != RECORDING::CLEAR_INDICES && config_container[ZONE_0]->GetWrt_AD_Statistics()) {
    AD::PrintStatistics(SU2_MPI::GetComm(), rank == MASTER_NODE);
  }

  AD::StopRecording();

}

void CDiscAdjSinglezoneDriver::DirectRunFixedPoint(RECORDING kind_recording) {
  /*--- Mesh movement ---*/

  direct_iteration->SetMesh_Deformation(geometry_container[ZONE_0][INST_0], solver, numerics, config, kind_recording);

  /*--- Zone preprocessing ---*/

  direct_iteration->Preprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Iterate the direct solver ---*/

  direct_iteration->Iterate(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Postprocess the direct solver ---*/

  direct_iteration->Postprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Print the direct residual to screen. ---*/
  PrintDirectResidual(kind_recording);
}

void CDiscAdjSinglezoneDriver::DirectRunResidual(RECORDING kind_recording) {
  /*--- Deform the mesh. ---*/
  DeformGeometry();

  /*--- Pre-process the primal solver state. ---*/
  UpdateTimeIter();
  UpdateFarfield();
  UpdateGeometry();
  UpdateStates();

  /*--- Run the computation of the outputs of interest. ---*/
  UpdateResiduals();
  UpdateTractions();
  UpdateObjective();

  /*--- Print the direct residual to screen. ---*/
  PrintDirectResidual(kind_recording);
}

void CDiscAdjSinglezoneDriver::MainRecording(){
  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with
   *    RECORDING::CLEAR_INDICES as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(RECORDING::CLEAR_INDICES);

  /*--- Store the computational graph of one direct iteration with the solution variables as input. ---*/

  SetRecording(MainVariables);

}

void CDiscAdjSinglezoneDriver::SecondaryRecording(){
  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with
   *    RECORDING::CLEAR_INDICES as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(RECORDING::CLEAR_INDICES);

  /*--- Store the computational graph of one direct iteration with the secondary variables as input. ---*/

  SetRecording(SecondaryVariables);
}

void CDiscAdjSinglezoneDriver::SecondaryRunFixedPoint() {
  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the current iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/
  SetAdjointObjective();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed sensitivity values. ---*/

  if (SecondaryVariables == RECORDING::MESH_COORDS) {
    solver[MainSolver]->SetSensitivity(geometry, config);
  } else {
    solver[ADJMESH_SOL]->SetSensitivity(geometry, config, solver[MainSolver]);
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/
  AD::ClearAdjoints();
}

void CDiscAdjSinglezoneDriver::SecondaryRunResidual() {
  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *--- of the previous iteration. The values are passed to the AD tool.
   *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the adjoints of the residuals and store them for the next iteration ---*/

  if (config->GetFluidProblem()) {
    if (SecondaryVariables == RECORDING::MESH_COORDS) {
      solver[ADJFLOW_SOL]->ExtractAdjoint_Coordinates(geometry, config, nullptr, ENUM_VARIABLE::RESIDUALS);
    } else {
      solver[ADJFLOW_SOL]->ExtractAdjoint_Coordinates(geometry, config, solver[ADJMESH_SOL], ENUM_VARIABLE::RESIDUALS);
    }
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdjointObjective();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the adjoints of the objective function and store them for the next iteration ---*/

  if (config->GetFluidProblem()) {
    if (SecondaryVariables == RECORDING::MESH_COORDS) {
      solver[ADJFLOW_SOL]->ExtractAdjoint_Coordinates(geometry, config, nullptr, ENUM_VARIABLE::OBJECTIVE);
    } else {
      solver[ADJFLOW_SOL]->ExtractAdjoint_Coordinates(geometry, config, solver[ADJMESH_SOL], ENUM_VARIABLE::OBJECTIVE);
    }
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  /*--- Initialize the adjoint of the vertex tractions with the corresponding adjoint vector. ---*/

  solver[FLOW_SOL]->SetVertexTractionsAdjoint(geometry, config);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the adjoints of the vertex tractions and store them for the next iteration ---*/

  if (config->GetFluidProblem()) {
    if (SecondaryVariables == RECORDING::MESH_COORDS) {
      solver[ADJFLOW_SOL]->ExtractAdjoint_Coordinates(geometry, config, nullptr, ENUM_VARIABLE::TRACTIONS);
    } else {
      solver[ADJFLOW_SOL]->ExtractAdjoint_Coordinates(geometry, config, solver[ADJMESH_SOL], ENUM_VARIABLE::TRACTIONS);
    }
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  /*--- Skip the derivation of the mesh solver if it is not defined ---*/
  if (SecondaryVariables == RECORDING::MESH_DEFORM) {
     /*--- Initialize the adjoint of the volume coordinates with the corresponding adjoint vector. ---*/
     SU2_OMP_PARALLEL_(if(solver[ADJMESH_SOL]->GetHasHybridParallel())) {

         /*--- Initialize the adjoints of the volume coordinates ---*/
         solver[ADJMESH_SOL]->SetAdjoint_Output(geometry, config);
     }
     END_SU2_OMP_PARALLEL

     /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/
     AD::ComputeAdjoint();

     /*--- Extract the adjoints of the volume coordinates and store them for the next iteration ---*/
     if (config->GetFluidProblem()) {
         solver[ADJFLOW_SOL]->ExtractAdjoint_Coordinates(geometry, config, solver[ADJMESH_SOL],
         ENUM_VARIABLE::COORDINATES);
     }

     /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/
     AD::ClearAdjoints();
  }

  /*--- Extract the adjoints of the residuals and store them for the next iteration ---*/
  if (config->GetFluidProblem()) {
    solver[ADJFLOW_SOL]->SetSensitivity(geometry, config);
  }
}

void CDiscAdjSinglezoneDriver::UpdateTimeIter() {
  /*--- Update the primal time iteration. ---*/
  solver[FLOW_SOL]->SetTime_Step(geometry, solver, config, MESH_0, config->GetTimeIter());
}

void CDiscAdjSinglezoneDriver::UpdateFarfield() {
  /*--- Update the primal far-field variables. ---*/

  su2double Velocity_Ref = config->GetVelocity_Ref();
  su2double Alpha = config->GetAoA() * PI_NUMBER / 180.0;
  su2double Beta = config->GetAoS() * PI_NUMBER / 180.0;
  su2double Mach = config->GetMach();
  su2double Temperature = config->GetTemperature_FreeStream();
  su2double Gas_Constant = config->GetGas_Constant();
  su2double Gamma = config->GetGamma();
  su2double SoundSpeed = sqrt(Gamma * Gas_Constant * Temperature);

  if (nDim == 2) {
    config->GetVelocity_FreeStreamND()[0] = cos(Alpha) * Mach * SoundSpeed / Velocity_Ref;
    config->GetVelocity_FreeStreamND()[1] = sin(Alpha) * Mach * SoundSpeed / Velocity_Ref;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStreamND()[0] = cos(Alpha) * cos(Beta) * Mach * SoundSpeed / Velocity_Ref;
    config->GetVelocity_FreeStreamND()[1] = sin(Beta) * Mach * SoundSpeed / Velocity_Ref;
    config->GetVelocity_FreeStreamND()[2] = sin(Alpha) * cos(Beta) * Mach * SoundSpeed / Velocity_Ref;
  }
}

void CDiscAdjSinglezoneDriver::UpdateGeometry() {
  /*--- Update the geometry (i.e. dual grid, without multi-grid). ---*/

  geometry->InitiateComms(geometry, config, MPI_QUANTITIES::COORDINATES);
  geometry->CompleteComms(geometry, config, MPI_QUANTITIES::COORDINATES);

  geometry->SetControlVolume(config, UPDATE);
  geometry->SetBoundControlVolume(config, UPDATE);
  geometry->SetMaxLength(config);
}

void CDiscAdjSinglezoneDriver::DeformGeometry() {
  /*--- Deform the geometry. ---*/

  direct_iteration->SetMesh_Deformation(geometry_container[ZONE_0][INST_0], solver, numerics, config,
                                        SecondaryVariables);
}

void CDiscAdjSinglezoneDriver::UpdateObjective() {
  ObjFunc = 0.0;

  /*--- Specific scalar objective functions ---*/
  switch (config->GetKind_Solver()) {
    case MAIN_SOLVER::DISC_ADJ_INC_EULER:       case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:      case MAIN_SOLVER::DISC_ADJ_INC_RANS:
    case MAIN_SOLVER::DISC_ADJ_EULER:           case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:          case MAIN_SOLVER::DISC_ADJ_RANS:
    case MAIN_SOLVER::DISC_ADJ_FEM_EULER:       case MAIN_SOLVER::DISC_ADJ_FEM_NS:                 case MAIN_SOLVER::DISC_ADJ_FEM_RANS:

      /*--- Surface based obj. function ---*/
      direct_output->SetHistoryOutput(geometry, solver, config, config->GetTimeIter(),
                                      config->GetOuterIter(), config->GetInnerIter());
      ObjFunc += solver[FLOW_SOL]->GetTotal_ComboObj();
      break;

    case MAIN_SOLVER::DISC_ADJ_HEAT:
      direct_output->SetHistoryOutput(geometry, solver, config, config->GetTimeIter(),
                                      config->GetOuterIter(), config->GetInnerIter());
      ObjFunc = solver[HEAT_SOL]->GetTotal_ComboObj();
      break;

    case MAIN_SOLVER::DISC_ADJ_FEM:
      solver[FEA_SOL]->Postprocessing(geometry, config, numerics_container[ZONE_0][INST_0][MESH_0][FEA_SOL], true);

      direct_output->SetHistoryOutput(geometry, solver, config, config->GetTimeIter(),
                                      config->GetOuterIter(), config->GetInnerIter());
      ObjFunc = solver[FEA_SOL]->GetTotal_ComboObj();
      break;

    default:
      break;
  }
}

void CDiscAdjSinglezoneDriver::UpdateStates() {
  /*--- Update the flow and turbulent conservative state variables, preparing for other updates. ---*/
  direct_iteration->Preprocess(direct_output, integration_container, geometry_container, solver_container,
                               numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0,
                               INST_0);
}

void CDiscAdjSinglezoneDriver::UpdateResiduals() {
  /*--- Update the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) residuals and objective function (no
   * system solve). ---*/
  integration[FLOW_SOL]->ComputeResiduals(geometry_container, solver_container, numerics_container, config_container,
                                          FLOW_SOL, ZONE_0, INST_0);
}

void CDiscAdjSinglezoneDriver::UpdateTractions() {
  /*--- Update the surface tractions. ---*/
  direct_iteration->Postprocess(direct_output, integration_container, geometry_container, solver_container,
                                numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0,
                                INST_0);
}

void CDiscAdjSinglezoneDriver::UpdateJacobians() {
  /*--- Compute the approximate Jacobian for preconditioning. ---*/
  solver[FLOW_SOL]->PrepareImplicitIteration(geometry, solver, config);
}

void CDiscAdjSinglezoneDriver::ApplyPreconditioner(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) {
  /*--- Use an approximate diagonal preconditioning based on the transpose of the primal Jacobian. ---*/
  (*PrimalPreconditioner)(u, v);

  /*--- Apply a few FGMRES iterations in addition to the above preconditioner. ---*/
  v.SetValZero();

  Scalar KrylovPreEps = KrylovPreTol;
  auto MaxIter = 5;
  solver[FLOW_SOL]->System.FGMRES_LinSolver(u, v, *PrimalJacobian, *PrimalPreconditioner, KrylovPreTol, MaxIter,
                                                    KrylovPreEps, false, config);
}

void CDiscAdjSinglezoneDriver::ApplyOperator(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) {
  /*--- Set the adjoint variables used in the seeding of the tape. ---*/
  SetAllSolutions(ZONE_0, true, u);

  /*--- Evaluate the tape to and extract the partial derivatives. ---*/
  UpdateAdjointsResidual();

  /*--- Extract the partial residual Jacobian-adjoint product. ---*/
  GetAllResidualsStatesSensitivities(v);
}
