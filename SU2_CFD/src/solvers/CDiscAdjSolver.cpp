/*!
 * \file CDiscAdjSolver.cpp
 * \brief Main subroutines for solving the discrete adjoint problem.
 * \author T. Albring
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

#include "../../include/solvers/CDiscAdjSolver.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config, CSolver *direct_sol,
                               unsigned short Kind_Solver, unsigned short iMesh)  : CSolver() {

  /*-- Store some information about direct solver ---*/

  KindDirect_Solver = Kind_Solver;
  direct_solver = direct_sol;

  adjoint = true;

  nVar = direct_solver->GetnVar();
  nDim = geometry->GetnDim();

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), OMP_MAX_SIZE);

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual_RMS.resize(nVar,1.0);
  Residual_Max.resize(nVar,1.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

  if (config->GetMultizone_Residual()){

    Residual_BGS.resize(nVar,1.0);
    Residual_Max_BGS.resize(nVar,1.0);
    Point_Max_BGS.resize(nVar,0);
    Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
  }

  /*--- Sensitivity definition and coefficient in all the markers ---*/

  CSensitivity.resize(nMarker);
  for (auto iMarker = 0ul; iMarker < nMarker; iMarker++) {
    const auto nVertex = geometry->nVertex[iMarker];
    CSensitivity[iMarker].resize(nVertex, 0.0);
  }

  Sens_Geo.resize(config->GetnMarker_Monitoring(), 0.0);

  /*--- Initialize the discrete adjoint solution to zero everywhere. ---*/

  if (nVar > MAXNVAR) {
    SU2_MPI::Error("Oops! The CDiscAdjSolver static array sizes are not large enough.",CURRENT_FUNCTION);
  }

  vector<su2double> Solution(nVar,1e-16);
  nodes = new CDiscAdjVariable(Solution.data(), nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Allocate extra solution variables, if any are in use. ---*/

  const auto nVarExtra = direct_solver->RegisterSolutionExtra(true, config);
  nodes->AllocateAdjointSolutionExtra(nVarExtra);

  switch(KindDirect_Solver){
  case RUNTIME_FLOW_SYS:
    SolverName = "ADJ.FLOW";
    break;
  case RUNTIME_HEAT_SYS:
    SolverName = "ADJ.HEAT";
    break;
  case RUNTIME_TURB_SYS:
    SolverName = "ADJ.TURB";
    break;
  case RUNTIME_SPECIES_SYS:
    SolverName = "ADJ.SPECIES";
    break;
  case RUNTIME_RADIATION_SYS:
    SolverName = "ADJ.RAD";
    break;
  default:
    SolverName = "ADJ.SOL";
    break;
  }
}

CDiscAdjSolver::~CDiscAdjSolver() { delete nodes; }

void CDiscAdjSolver::SetRecording(CGeometry* geometry, CConfig *config){

  const bool time_n1_needed = config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND;
  const bool time_n_needed = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) || time_n1_needed;

  /*--- Reset the solution to the initial (converged) solution ---*/

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    direct_solver->GetNodes()->SetSolution(iPoint, nodes->GetSolution_Direct(iPoint));
  }
  END_SU2_OMP_FOR

  if (time_n_needed) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
      for (auto iVar = 0u; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n(iPoint)[iVar]);
      }
    }
    END_SU2_OMP_FOR
  }
  if (time_n1_needed) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
      for (auto iVar = 0u; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n1(iPoint)[iVar]);
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CDiscAdjSolver::RegisterSolution(CGeometry *geometry, CConfig *config) {

  const bool time_n1_needed = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool time_n_needed  = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) || time_n1_needed;

  /*--- Register solution at all necessary time instances and other variables on the tape ---*/

  /*--- Boolean true indicates that an input is registered ---*/
  direct_solver->GetNodes()->RegisterSolution(true);

  direct_solver->RegisterSolutionExtra(true, config);

  if (time_n_needed)
    direct_solver->GetNodes()->RegisterSolution_time_n();

  if (time_n1_needed)
    direct_solver->GetNodes()->RegisterSolution_time_n1();
}

void CDiscAdjSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset) {

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {

  /*--- Register farfield values as input ---*/

  if((config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS && !config->GetBoolTurbomachinery())) {

    su2double Velocity_Ref = config->GetVelocity_Ref();
    Alpha                  = config->GetAoA()*PI_NUMBER/180.0;
    Beta                   = config->GetAoS()*PI_NUMBER/180.0;
    Mach                   = config->GetMach();
    Pressure               = config->GetPressure_FreeStreamND();
    Temperature            = config->GetTemperature_FreeStreamND();

    su2double SoundSpeed = 0.0;

    if (nDim == 2) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*Mach); }
    if (nDim == 3) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*cos(Beta)*Mach); }

    if (!reset) {
      AD::RegisterInput(Mach);
      AD::RegisterInput(Alpha);
      AD::RegisterInput(Temperature);
      AD::RegisterInput(Pressure);
    }

    /*--- Recompute the free stream velocity ---*/

    if (nDim == 2) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Alpha)*Mach*SoundSpeed/Velocity_Ref;
    }
    if (nDim == 3) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[2] = sin(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
    }

    config->SetTemperature_FreeStreamND(Temperature);
    direct_solver->SetTemperature_Inf(Temperature);
    config->SetPressure_FreeStreamND(Pressure);
    direct_solver->SetPressure_Inf(Pressure);

  }

  if ((config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && config->GetBoolTurbomachinery()){

    BPressure = config->GetPressureOut_BC();
    Temperature = config->GetTotalTemperatureIn_BC();

    if (!reset){
      AD::RegisterInput(BPressure);
      AD::RegisterInput(Temperature);
    }

    config->SetPressureOut_BC(BPressure);
    config->SetTotalTemperatureIn_BC(Temperature);
  }

  /*--- Register incompressible initialization values as input ---*/

  if ((config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) &&
      ((KindDirect_Solver == RUNTIME_FLOW_SYS &&
        (!config->GetBoolTurbomachinery())))) {

    /*--- Access the velocity (or pressure) and temperature at the
     inlet BC and the back pressure at the outlet. Note that we are
     assuming that have internal flow, which will be true for the
     majority of cases. External flows with far-field BCs will report
     zero for these sensitivities. ---*/

    ModVel    = config->GetIncInlet_BC();
    BPressure = config->GetIncPressureOut_BC();
    Temperature = config->GetIncTemperature_BC();

    /*--- Register the variables for AD. ---*/

    if (!reset) {
      AD::RegisterInput(ModVel);
      AD::RegisterInput(BPressure);
      AD::RegisterInput(Temperature);
    }

    /*--- Set the BC values in the config class. ---*/

    config->SetIncInlet_BC(ModVel);
    config->SetIncPressureOut_BC(BPressure);
    config->SetIncTemperature_BC(Temperature);

  }

  /*--- Register incompressible radiation values as input ---*/

  if ((config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) &&
      ((KindDirect_Solver == RUNTIME_RADIATION_SYS &&
        (!config->GetBoolTurbomachinery())))) {

    /*--- Access the nondimensional freestream temperature. ---*/

    TemperatureRad = config->GetTemperature_FreeStreamND();

    /*--- Register the variables for AD. ---*/

    if (!reset) {
      AD::RegisterInput(TemperatureRad);
    }

    /*--- Set the temperature at infinity in the direct solver class. ---*/

    direct_solver->SetTemperature_Inf(TemperatureRad);

  }

  /*--- Here it is possible to register other variables as input that influence the flow solution
   * and thereby also the objective function. The adjoint values (i.e. the derivatives) can be
   * extracted in the ExtractAdjointVariables routine. ---*/

  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CDiscAdjSolver::RegisterOutput(CGeometry *geometry, CConfig *config) {

  /*--- Register variables as output of the solver iteration. Boolean false indicates that an output is registered ---*/

  direct_solver->GetNodes()->RegisterSolution(false);

  direct_solver->RegisterSolutionExtra(false, config);
}

void CDiscAdjSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config, bool CrossTerm) {

  const bool time_n1_needed = config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND;
  const bool time_n_needed = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) || time_n1_needed;

  const su2double relax = (config->GetInnerIter()==0) ? 1.0 : config->GetRelaxation_Factor_Adjoint();

  /*--- Thread-local residual variables. ---*/
  su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
  unsigned long idxMax[MAXNVAR] = {0};

  /*--- Set the old solution and compute residuals. ---*/

  if (!config->GetMultizone_Problem()) nodes->Set_OldSolution();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {

    /*--- Extract the adjoint solution ---*/

    su2double Solution[MAXNVAR] = {0.0};
    direct_solver->GetNodes()->GetAdjointSolution(iPoint,Solution);

    /*--- Relax and store the adjoint solution, compute the residuals. ---*/

    for (auto iVar = 0u; iVar < nVar; iVar++) {
      su2double residual = Solution[iVar]-nodes->GetSolution_Old(iPoint,iVar);
      nodes->AddSolution(iPoint, iVar, relax*residual);

      if (iPoint < nPointDomain) {
        /*--- "Add" residual at (iPoint,iVar) to local residual variables. ---*/
        ResidualReductions_PerThread(iPoint,iVar,residual,resRMS,resMax,idxMax);
      }
    }
  }
  END_SU2_OMP_FOR

  direct_solver->ExtractAdjoint_SolutionExtra(nodes->GetSolutionExtra(), config);

  /*--- Residuals and time_n terms are not needed when evaluating multizone cross terms. ---*/
  if (CrossTerm) return;

  /*--- "Add" residuals from all threads to global residual variables. ---*/
  ResidualReductions_FromAllThreads(geometry, config, resRMS,resMax,idxMax);

  SU2_OMP_MASTER {
    SetIterLinSolver(direct_solver->System.GetIterations());
    SetResLinSolver(direct_solver->System.GetResidual());
  }
  END_SU2_OMP_MASTER

  /*--- Extract and store the adjoint of the primal solution at time n ---*/
  if (time_n_needed) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
      su2double Solution[MAXNVAR] = {0.0};
      direct_solver->GetNodes()->GetAdjointSolution_time_n(iPoint,Solution);
      nodes->Set_Solution_time_n(iPoint,Solution);
    }
    END_SU2_OMP_FOR
  }

  /*--- Extract and store the adjoint of the primal solution at time n-1 ---*/
  if (time_n1_needed) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
      su2double Solution[MAXNVAR] = {0.0};
      direct_solver->GetNodes()->GetAdjointSolution_time_n1(iPoint,Solution);
      nodes->Set_Solution_time_n1(iPoint,Solution);
    }
    END_SU2_OMP_FOR
  }

}

void CDiscAdjSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) {

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {

  /*--- Extract the adjoint values of the farfield values ---*/

  if ((config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && !config->GetBoolTurbomachinery()) {
    su2double Local_Sens_Press, Local_Sens_Temp, Local_Sens_AoA, Local_Sens_Mach;

    Local_Sens_Mach  = SU2_TYPE::GetDerivative(Mach);
    Local_Sens_AoA   = SU2_TYPE::GetDerivative(Alpha);
    Local_Sens_Temp  = SU2_TYPE::GetDerivative(Temperature);
    Local_Sens_Press = SU2_TYPE::GetDerivative(Pressure);

    SU2_MPI::Allreduce(&Local_Sens_Mach,  &Total_Sens_Mach,  1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&Local_Sens_AoA,   &Total_Sens_AoA,   1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&Local_Sens_Temp,  &Total_Sens_Temp,  1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&Local_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  }

  if ((config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && config->GetBoolTurbomachinery()){
    su2double Local_Sens_BPress, Local_Sens_Temperature;

    Local_Sens_BPress = SU2_TYPE::GetDerivative(BPressure);
    Local_Sens_Temperature = SU2_TYPE::GetDerivative(Temperature);

    SU2_MPI::Allreduce(&Local_Sens_BPress,   &Total_Sens_BPress,   1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&Local_Sens_Temperature,   &Total_Sens_Temp,   1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  }

  if ((config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) &&
      (KindDirect_Solver == RUNTIME_FLOW_SYS &&
       (!config->GetBoolTurbomachinery()))) {

    su2double Local_Sens_ModVel, Local_Sens_BPress, Local_Sens_Temp;

    Local_Sens_ModVel = SU2_TYPE::GetDerivative(ModVel);
    Local_Sens_BPress = SU2_TYPE::GetDerivative(BPressure);
    Local_Sens_Temp   = SU2_TYPE::GetDerivative(Temperature);

    SU2_MPI::Allreduce(&Local_Sens_ModVel, &Total_Sens_ModVel, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&Local_Sens_BPress, &Total_Sens_BPress, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&Local_Sens_Temp,   &Total_Sens_Temp,   1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  }

  if ((config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) &&
      (KindDirect_Solver == RUNTIME_RADIATION_SYS &&
       (!config->GetBoolTurbomachinery()))) {

    su2double Local_Sens_Temp_Rad;
    Local_Sens_Temp_Rad   = SU2_TYPE::GetDerivative(TemperatureRad);

    SU2_MPI::Allreduce(&Local_Sens_Temp_Rad, &Total_Sens_Temp_Rad, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

    /*--- Store it in the Total_Sens_Temp container so it's accessible without the need of a new method ---*/
    Total_Sens_Temp = Total_Sens_Temp_Rad;

  }

  /*--- Extract here the adjoint values of everything else that is registered as input in RegisterInput. ---*/

  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CDiscAdjSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config) {

  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST ||
                          config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool multizone = config->GetMultizone_Problem();

  /*--- Local container to manipulate the adjoint solution. ---*/
  su2double Solution[MAXNVAR] = {0.0};

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {

    /*--- Get and store the adjoint solution of a point. ---*/
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      Solution[iVar] = nodes->GetSolution(iPoint,iVar);
    }

    /*--- Add dual time contributions to the adjoint solution. Two terms stored for DT-2nd-order. ---*/
    if (dual_time && !multizone) {
      for (auto iVar = 0u; iVar < nVar; iVar++) {
        Solution[iVar] += nodes->GetDual_Time_Derivative(iPoint,iVar);
      }
    }

    /*--- Set the adjoint values of the primal solution. ---*/

    direct_solver->GetNodes()->SetAdjointSolution(iPoint,Solution);
  }
  END_SU2_OMP_FOR

  direct_solver->SetAdjoint_SolutionExtra(nodes->GetSolutionExtra(), config);
}

void CDiscAdjSolver::SetSensitivity(CGeometry *geometry, CConfig *config, CSolver*) {

  SU2_OMP_PARALLEL {

  const bool time_stepping = (config->GetTime_Marching() != TIME_MARCHING::STEADY);
  const su2double eps = config->GetAdjSharp_LimiterCoeff()*config->GetRefElemLength();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {

    auto Coord = geometry->nodes->GetCoord(iPoint);

    for (auto iDim = 0u; iDim < nDim; iDim++) {

      su2double Sensitivity = geometry->nodes->GetAdjointSolution(iPoint, iDim);
      AD::ResetInput(Coord[iDim]);

      /*--- If sharp edge, set the sensitivity to 0 on that region ---*/

      if (config->GetSens_Remove_Sharp() && geometry->nodes->GetSharpEdge_Distance(iPoint) < eps) {
        Sensitivity = 0.0;
      }
      if (!time_stepping) {
        nodes->SetSensitivity(iPoint,iDim, Sensitivity);
      } else {
        nodes->SetSensitivity(iPoint,iDim, nodes->GetSensitivity(iPoint,iDim) + Sensitivity);
      }
    }
  }
  END_SU2_OMP_FOR

  SetSurface_Sensitivity(geometry, config);

  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config) {

  SU2_OMP_MASTER
  for (auto& x : Sens_Geo) x = 0.0;
  END_SU2_OMP_MASTER

  /*--- Loop over boundary markers to select those for Euler walls and NS walls ---*/

  for (auto iMarker = 0ul; iMarker < nMarker; iMarker++) {

    if (!config->GetSolid_Wall(iMarker)) continue;

    su2double Sens = 0.0;

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0ul; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Projection of the gradient calculated with AD onto the normal vector of the surface ---*/

      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      const auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

      su2double Sens_Vertex = 0.0;
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        Sens_Vertex += Normal[iDim] * nodes->GetSensitivity(iPoint,iDim);
      }
      Sens_Vertex /= GeometryToolbox::Norm(nDim, Normal);

      CSensitivity[iMarker][iVertex] = -Sens_Vertex;
      Sens += pow(Sens_Vertex,2);
    }
    END_SU2_OMP_FOR

    if (config->GetMarker_All_Monitoring(iMarker) == NO) continue;

    /*--- Compute sensitivity for each surface point ---*/

    const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);

    for (size_t iMarker_Mon = 0; iMarker_Mon < Sens_Geo.size(); iMarker_Mon++) {
      if (Marker_Tag == config->GetMarker_Monitoring_TagBound(iMarker_Mon)) {
        atomicAdd(Sens, Sens_Geo[iMarker_Mon]);
        break;
      }
    }
  }

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    auto local = Sens_Geo;
    SU2_MPI::Allreduce(local.data(), Sens_Geo.data(), Sens_Geo.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

    Total_Sens_Geo = 0.0;
    for (auto& x : Sens_Geo) {
      x = sqrt(x);
      Total_Sens_Geo += x;
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CDiscAdjSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                   unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_OMP_MASTER
  config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);
  END_SU2_OMP_MASTER

  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                         (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);

  if (!dual_time) return;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint<geometry->GetnPoint(); iPoint++) {
    const auto solution_n = nodes->GetSolution_time_n(iPoint);
    const auto solution_n1 = nodes->GetSolution_time_n1(iPoint);

    for (auto iVar = 0u; iVar < nVar; iVar++) {
      nodes->SetDual_Time_Derivative(iPoint, iVar, solution_n[iVar]+nodes->GetDual_Time_Derivative_n(iPoint, iVar));
      nodes->SetDual_Time_Derivative_n(iPoint,iVar, solution_n1[iVar]);
    }
  }
  END_SU2_OMP_FOR
}

void CDiscAdjSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  auto filename = config->GetSolution_AdjFileName();
  auto restart_filename = config->GetObjFunc_Extension(filename);
  restart_filename = config->GetFilename(restart_filename, "", val_iter);

  const bool rans = (config->GetKind_Turb_Model() != TURB_MODEL::NONE);

  /*--- Skip coordinates ---*/
  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Skip flow adjoint variables ---*/
  if (KindDirect_Solver == RUNTIME_TURB_SYS) {
    skipVars += nDim + 2;
  }

  if (KindDirect_Solver == RUNTIME_SPECIES_SYS) {
    // Skip the number of Flow Vars and Turb Vars to get to the adjoint species vars
    skipVars += nDim + 2;
    if (rans) skipVars += solver[MESH_0][TURB_SOL]->GetnVar();
  }

  /*--- Skip flow adjoint and turbulent variables ---*/
  if (KindDirect_Solver == RUNTIME_RADIATION_SYS) {
    skipVars += nDim + 2;
    if (rans) skipVars += solver[MESH_0][TURB_SOL]->GetnVar();
  }

  BasicLoadRestart(geometry[MESH_0], config, restart_filename, skipVars);

  /*--- Interpolate solution on coarse grids ---*/

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {
    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][ADJFLOW_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][ADJFLOW_SOL]->GetNodes()->GetSolution());
  }
}
