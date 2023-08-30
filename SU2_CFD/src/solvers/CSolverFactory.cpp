/*!
 * \file CSolverFactory.cpp
 * \brief Main subroutines for CSolverFactoryclass.
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

#include "../../include/solvers/CSolver.hpp"
#include "../../include/solvers/CSolverFactory.hpp"
#include "../../include/solvers/CEulerSolver.hpp"
#include "../../include/solvers/CIncEulerSolver.hpp"
#include "../../include/solvers/CNSSolver.hpp"
#include "../../include/solvers/CIncNSSolver.hpp"
#include "../../include/solvers/CNEMOEulerSolver.hpp"
#include "../../include/solvers/CNEMONSSolver.hpp"
#include "../../include/solvers/CTurbSASolver.hpp"
#include "../../include/solvers/CTurbSSTSolver.hpp"
#include "../../include/solvers/CTransLMSolver.hpp"
#include "../../include/solvers/CAdjEulerSolver.hpp"
#include "../../include/solvers/CAdjNSSolver.hpp"
#include "../../include/solvers/CAdjTurbSolver.hpp"
#include "../../include/solvers/CHeatSolver.hpp"
#include "../../include/solvers/CFEASolver.hpp"
#include "../../include/solvers/CTemplateSolver.hpp"
#include "../../include/solvers/CDiscAdjSolver.hpp"
#include "../../include/solvers/CDiscAdjFEASolver.hpp"
#include "../../include/solvers/CFEM_DG_EulerSolver.hpp"
#include "../../include/solvers/CFEM_DG_NSSolver.hpp"
#include "../../include/solvers/CMeshSolver.hpp"
#include "../../include/solvers/CDiscAdjMeshSolver.hpp"
#include "../../include/solvers/CBaselineSolver.hpp"
#include "../../include/solvers/CBaselineSolver_FEM.hpp"
#include "../../include/solvers/CRadP1Solver.hpp"
#include "../../include/solvers/CSpeciesSolver.hpp"
#include "../../include/solvers/CSpeciesFlameletSolver.hpp"

map<const CSolver*, SolverMetaData> CSolverFactory::allocatedSolvers;

CSolver** CSolverFactory::CreateSolverContainer(MAIN_SOLVER kindMainSolver, CConfig *config, CGeometry *geometry, int iMGLevel){

  CSolver** solver;

  solver = new CSolver*[MAX_SOLS]();

  switch (kindMainSolver) {
    case MAIN_SOLVER::TEMPLATE_SOLVER:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::TEMPLATE, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::INC_EULER:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::INC_EULER, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]  = CreateSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::EULER:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::EULER, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::NEMO_EULER:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::NEMO_EULER, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::INC_NAVIER_STOKES:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]  = CreateSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::NAVIER_STOKES:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::NEMO_NAVIER_STOKES:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::NEMO_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::RANS:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[TURB_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[TRANS_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::TRANSITION, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::INC_RANS:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      solver[TURB_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[TRANS_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::TRANSITION, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]  = CreateSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::HEAT_EQUATION:
      solver[HEAT_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::ADJ_EULER:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::CONT_ADJ_EULER, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::ADJ_NAVIER_STOKES:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::CONT_ADJ_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::ADJ_RANS:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::CONT_ADJ_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[ADJTURB_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::CONT_ADJ_TURB, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_EULER:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      solver[ADJSPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_SPECIES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_RANS:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[ADJTURB_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_TURB, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      solver[ADJSPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_SPECIES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_INC_EULER:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::INC_EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]     = CreateSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      solver[ADJRAD_SOL]  = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_RADIATION, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[ADJHEAT_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_HEAT, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]     = CreateSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      solver[ADJRAD_SOL]  = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_RADIATION, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      solver[ADJSPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_SPECIES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[ADJHEAT_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_HEAT, solver, geometry, config, iMGLevel);
      solver[SPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::SPECIES, solver, geometry, config, iMGLevel);
      solver[ADJSPECIES_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_SPECIES, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[ADJTURB_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_TURB, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]     = CreateSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      solver[ADJRAD_SOL]  = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_RADIATION, solver, geometry, config, iMGLevel);
    break;
    case MAIN_SOLVER::DISC_ADJ_HEAT:
      solver[HEAT_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[ADJHEAT_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_HEAT, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::FEM_ELASTICITY:
      solver[FEA_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::FEA, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_FEM:
      solver[FEA_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::FEA, solver, geometry, config, iMGLevel);
      solver[ADJFEA_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FEA, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::FEM_EULER:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DG_EULER, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::FEM_NAVIER_STOKES: case MAIN_SOLVER::FEM_LES:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DG_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::FEM_RANS:
      SU2_MPI::Error("FEM RANS not available", CURRENT_FUNCTION);
      break;
    case MAIN_SOLVER::DISC_ADJ_FEM_EULER:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::DG_EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_FEM_NS:
      solver[FLOW_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::DG_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
      SU2_MPI::Error("Adjoint FEM RANS not available", CURRENT_FUNCTION);
      break;
     default:
      solver = nullptr;
  }

  solver[MESH_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::MESH, solver, geometry, config, iMGLevel);
  solver[ADJMESH_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_MESH, solver, geometry, config, iMGLevel);

  return solver;

}

CSolver* CSolverFactory::CreateSubSolver(SUB_SOLVER_TYPE kindSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *genericSolver = nullptr;

  TURB_MODEL kindTurbModel = config->GetKind_Turb_Model();
  TURB_TRANS_MODEL kindTransModel = config->GetKind_Trans_Model();

  SolverMetaData metaData;

  metaData.solverType = kindSolver;

  switch (kindSolver) {
    case SUB_SOLVER_TYPE::CONT_ADJ_EULER:
      genericSolver = new CAdjEulerSolver(geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::MULTIGRID;
      break;
    case SUB_SOLVER_TYPE::CONT_ADJ_NAVIER_STOKES:
      genericSolver = new CAdjNSSolver(geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::MULTIGRID;
      break;
    case SUB_SOLVER_TYPE::CONT_ADJ_TURB:
      genericSolver = CreateTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, true);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_TURB:
      genericSolver = CreateTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, true);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::BASELINE:
      genericSolver = new CBaselineSolver(geometry, config);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::BASELINE_FEM:
      genericSolver = new CBaselineSolver_FEM(geometry, config);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_FEA:
      genericSolver = new CDiscAdjFEASolver(geometry, config, solver[FEA_SOL], RUNTIME_FEA_SYS, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_MESH:
      genericSolver = CreateMeshSolver(solver, geometry, config, iMGLevel, true);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_FLOW:
      genericSolver = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::EULER:
    case SUB_SOLVER_TYPE::INC_EULER:
    case SUB_SOLVER_TYPE::NEMO_EULER:
    case SUB_SOLVER_TYPE::NAVIER_STOKES:
    case SUB_SOLVER_TYPE::INC_NAVIER_STOKES:
    case SUB_SOLVER_TYPE::NEMO_NAVIER_STOKES:
      genericSolver = CreateFlowSolver(kindSolver, solver, geometry, config, iMGLevel);
      if (!config->GetNewtonKrylov() || config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint())
        metaData.integrationType = INTEGRATION_TYPE::MULTIGRID;
      else
        metaData.integrationType = INTEGRATION_TYPE::NEWTON;
      break;
    case SUB_SOLVER_TYPE::FEA:
      genericSolver = new CFEASolver(geometry, config);
      metaData.integrationType = INTEGRATION_TYPE::STRUCTURAL;
      break;
    case SUB_SOLVER_TYPE::MESH:
      genericSolver = CreateMeshSolver(solver, geometry, config, iMGLevel, false);
      metaData.integrationType = INTEGRATION_TYPE::NONE;
      break;
    case SUB_SOLVER_TYPE::DG_EULER:
    case SUB_SOLVER_TYPE::DG_NAVIER_STOKES:
      genericSolver = CreateDGSolver(kindSolver, geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::FEM_DG;
      break;
    case SUB_SOLVER_TYPE::HEAT:
      genericSolver = CreateHeatSolver(solver, geometry, config, iMGLevel, false);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_HEAT:
      genericSolver = CreateHeatSolver(solver, geometry, config, iMGLevel, true);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::TRANSITION:
      genericSolver = CreateTransSolver(kindTransModel, solver, geometry, config, iMGLevel, false);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::SPECIES:
      genericSolver = CreateSpeciesSolver(solver, geometry, config, iMGLevel, false);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_SPECIES:
      genericSolver = CreateSpeciesSolver(solver, geometry, config, iMGLevel, true);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::TURB:
    case SUB_SOLVER_TYPE::TURB_SA:
    case SUB_SOLVER_TYPE::TURB_SST:
      genericSolver = CreateTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, false);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::TEMPLATE:
      genericSolver = new CTemplateSolver(geometry, config);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::RADIATION:
      if (config->AddRadiation()){
        genericSolver = new CRadP1Solver(geometry, config);
      }
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_RADIATION:
      if (config->AddRadiation()){
        genericSolver = new CDiscAdjSolver(geometry, config, solver[RAD_SOL], RUNTIME_RADIATION_SYS, iMGLevel);
      }
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    default:
      SU2_MPI::Error("No proper allocation found for requested sub solver", CURRENT_FUNCTION);
      break;
  }

  if (genericSolver != nullptr)
    allocatedSolvers[genericSolver] = metaData;

  return genericSolver;

}

CSolver* CSolverFactory::CreateTurbSolver(TURB_MODEL kindTurbModel, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, int adjoint){

  CSolver *turbSolver = nullptr;

  if (!adjoint){
    switch (TurbModelFamily(kindTurbModel)) {
      case TURB_FAMILY::SA:
        turbSolver = new CTurbSASolver(geometry, config, iMGLevel, solver[FLOW_SOL]->GetFluidModel());
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        turbSolver->Postprocessing(geometry, solver, config, iMGLevel);
        break;
      case TURB_FAMILY::KW:
        turbSolver = new CTurbSSTSolver(geometry, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        turbSolver->Postprocessing(geometry, solver, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        break;
      case TURB_FAMILY::NONE:
        SU2_MPI::Error("Trying to create TurbSolver container but TURB_MODEL=NONE.", CURRENT_FUNCTION);
        break;
    }
  } else {

    if (config->GetDiscrete_Adjoint()){
      if (!config->GetFrozen_Visc_Disc())
        turbSolver = new CDiscAdjSolver(geometry, config, solver[TURB_SOL], RUNTIME_TURB_SYS, iMGLevel);
    } else {
      if (!config->GetFrozen_Visc_Cont())
        turbSolver = new CAdjTurbSolver(geometry, config, iMGLevel);
    }
  }

  return turbSolver;
}

CSolver* CSolverFactory::CreateTransSolver(TURB_TRANS_MODEL kindTransModel, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, int adjoint){

  CSolver *transSolver = nullptr;

  if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
    switch (kindTransModel) {
      case TURB_TRANS_MODEL::LM :
        transSolver = new CTransLMSolver(geometry, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        transSolver->Postprocessing(geometry, solver, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        break;
      case TURB_TRANS_MODEL::NONE:
        break;
    }
  }

  return transSolver;
}

CSolver* CSolverFactory::CreateSpeciesSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint){

  CSolver *speciesSolver = nullptr;

  if (config->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    if (adjoint){
      speciesSolver = new CDiscAdjSolver(geometry, config, solver[SPECIES_SOL], RUNTIME_SPECIES_SYS, iMGLevel);
    } else {
      if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT)
        speciesSolver = new CSpeciesSolver(geometry, config, iMGLevel);
      else if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET)
        speciesSolver = new CSpeciesFlameletSolver(geometry, config, iMGLevel);
    }
  }
  return speciesSolver;
}

CSolver* CSolverFactory::CreateHeatSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint){

  CSolver *heatSolver = nullptr;

  /*--- Only allocate a heat solver if it should run standalone
   * or if the weakly coupled heat solver is enabled and no energy equation is included ---*/

  if ((config->GetWeakly_Coupled_Heat() && !config->GetEnergy_Equation()) || config->GetHeatProblem()) {
    if (adjoint) {
      if (config->GetDiscrete_Adjoint()) {
        heatSolver = new CDiscAdjSolver(geometry, config, solver[HEAT_SOL], RUNTIME_HEAT_SYS, iMGLevel);
      }
      else {
        SU2_MPI::Error("No continuous adjoint heat solver available.", CURRENT_FUNCTION);
      }
    }
    else {
      heatSolver = new CHeatSolver(geometry, config, iMGLevel);
    }
  }

  return heatSolver;
}

CSolver* CSolverFactory::CreateMeshSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint){

  CSolver *meshSolver = nullptr;

  if (config->GetDeform_Mesh() && iMGLevel == MESH_0) {
    if (!adjoint) {
      meshSolver = new CMeshSolver(geometry, config);
    }
    if (adjoint && config->GetDiscrete_Adjoint()) {
      meshSolver = new CDiscAdjMeshSolver(geometry, config, solver[MESH_SOL]);
    }
  }

  return meshSolver;
}

CSolver* CSolverFactory::CreateDGSolver(SUB_SOLVER_TYPE kindDGSolver, CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *DGSolver = nullptr;

  switch (kindDGSolver) {
    case SUB_SOLVER_TYPE::DG_EULER:
      if (config->GetKind_FEM_DG_Shock() == FEM_SHOCK_CAPTURING_DG::PERSSON){
        DGSolver = new CFEM_DG_NSSolver(geometry, config, iMGLevel);
      } else {
        DGSolver = new CFEM_DG_EulerSolver(geometry, config, iMGLevel);
      }
      break;
    case SUB_SOLVER_TYPE::DG_NAVIER_STOKES:
      DGSolver = new CFEM_DG_NSSolver(geometry, config, iMGLevel);
      break;
    default:
      SU2_MPI::Error("Requested DG solver not found", CURRENT_FUNCTION);
      break;
  }

  return DGSolver;
}

CSolver* CSolverFactory::CreateFlowSolver(SUB_SOLVER_TYPE kindFlowSolver, CSolver **solver,  CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *flowSolver = nullptr;

  switch (kindFlowSolver) {
    case SUB_SOLVER_TYPE::EULER:
      flowSolver = new CEulerSolver(geometry, config, iMGLevel);
      flowSolver->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      break;
    case SUB_SOLVER_TYPE::NAVIER_STOKES:
      flowSolver = new CNSSolver(geometry, config, iMGLevel);
      break;
    case SUB_SOLVER_TYPE::INC_EULER:
      flowSolver = new CIncEulerSolver(geometry, config, iMGLevel);
      flowSolver->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      break;
    case SUB_SOLVER_TYPE::INC_NAVIER_STOKES:
      flowSolver = new CIncNSSolver(geometry, config, iMGLevel);
      break;
    case SUB_SOLVER_TYPE::NEMO_EULER:
      flowSolver = new CNEMOEulerSolver(geometry, config, iMGLevel);
      break;
    case SUB_SOLVER_TYPE::NEMO_NAVIER_STOKES:
      flowSolver = new CNEMONSSolver(geometry, config, iMGLevel);
      break;
    default:
      SU2_MPI::Error("Flow solver not found", CURRENT_FUNCTION);
      break;
  }

  return flowSolver;
}
