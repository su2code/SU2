/*!
 * \file CSolverFactory.cpp
 * \brief Main subroutines for CSolverFactoryclass.
 * \author T. Albring
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

map<const CSolver*, SolverMetaData> CSolverFactory::allocatedSolvers;

CSolver** CSolverFactory::createSolverContainer(ENUM_MAIN_SOLVER kindMainSolver, CConfig *config, CGeometry *geometry, int iMGLevel){

  CSolver** solver;

  solver = new CSolver*[MAX_SOLS]();

  switch (kindMainSolver) {
    case TEMPLATE_SOLVER:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::TEMPLATE, solver, geometry, config, iMGLevel);
      break;
    case INC_EULER:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::INC_EULER, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]  = createSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      break;
    case EULER:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::EULER, solver, geometry, config, iMGLevel);
      break;
    case INC_NAVIER_STOKES:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL] = createSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]  = createSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      break;
    case NAVIER_STOKES:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case RANS:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[TURB_SOL] = createSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      break;
    case INC_RANS:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL] = createSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[TURB_SOL] = createSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]  = createSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      break;
    case HEAT_EQUATION:
      solver[HEAT_SOL] = createSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      break;
    case ADJ_EULER:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::CONT_ADJ_EULER, solver, geometry, config, iMGLevel);
      break;
    case ADJ_NAVIER_STOKES:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::CONT_ADJ_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case ADJ_RANS:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::CONT_ADJ_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = createSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[ADJTURB_SOL] = createSubSolver(SUB_SOLVER_TYPE::CONT_ADJ_TURB, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_EULER:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_NAVIER_STOKES:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_RANS:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = createSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[ADJTURB_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_TURB, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_INC_EULER:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::INC_EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]     = createSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      solver[ADJRAD_SOL]  = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_RADIATION, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_INC_NAVIER_STOKES:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL]    = createSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[ADJHEAT_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_HEAT, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]     = createSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      solver[ADJRAD_SOL]  = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_RADIATION, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_INC_RANS:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL]    = createSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[ADJHEAT_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_HEAT, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = createSubSolver(SUB_SOLVER_TYPE::TURB, solver, geometry, config, iMGLevel);
      solver[ADJTURB_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_TURB, solver, geometry, config, iMGLevel);
      solver[RAD_SOL]     = createSubSolver(SUB_SOLVER_TYPE::RADIATION, solver, geometry, config, iMGLevel);
      solver[ADJRAD_SOL]  = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_RADIATION, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_HEAT:
      solver[HEAT_SOL]    = createSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
      solver[ADJHEAT_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_HEAT, solver, geometry, config, iMGLevel);
      break;
    case FEM_ELASTICITY:
      solver[FEA_SOL] = createSubSolver(SUB_SOLVER_TYPE::FEA, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_FEM:
      solver[FEA_SOL]    = createSubSolver(SUB_SOLVER_TYPE::FEA, solver, geometry, config, iMGLevel);
      solver[ADJFEA_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FEA, solver, geometry, config, iMGLevel);
      break;
    case FEM_EULER:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DG_EULER, solver, geometry, config, iMGLevel);
      break;
    case FEM_NAVIER_STOKES: case FEM_LES:
      solver[FLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DG_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case FEM_RANS:
      SU2_MPI::Error("FEM RANS not available", CURRENT_FUNCTION);
      break;
    case DISC_ADJ_FEM_EULER:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::DG_EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_FEM_NS:
      solver[FLOW_SOL]    = createSubSolver(SUB_SOLVER_TYPE::DG_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_FLOW, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_FEM_RANS:
      SU2_MPI::Error("Adjoint FEM RANS not available", CURRENT_FUNCTION);
      break;
     default:
      solver = nullptr;
  }

  solver[MESH_SOL]    = createSubSolver(SUB_SOLVER_TYPE::MESH, solver, geometry, config, iMGLevel);
  solver[ADJMESH_SOL] = createSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_MESH, solver, geometry, config, iMGLevel);

  return solver;

}

CSolver* CSolverFactory::createSubSolver(SUB_SOLVER_TYPE kindSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *genericSolver = nullptr;

  ENUM_TURB_MODEL kindTurbModel = static_cast<ENUM_TURB_MODEL>(config->GetKind_Turb_Model());
  
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
      genericSolver = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, true);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_TURB:
      genericSolver = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, true);
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
      genericSolver = createMeshSolver(solver, geometry, config, iMGLevel, true);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_FLOW:
      genericSolver = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::EULER:
      genericSolver = createFlowSolver(SUB_SOLVER_TYPE::EULER, solver, geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::MULTIGRID;
      break;
    case SUB_SOLVER_TYPE::NAVIER_STOKES:
      genericSolver = createFlowSolver(SUB_SOLVER_TYPE::NAVIER_STOKES, solver, geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::MULTIGRID;
      break;
    case SUB_SOLVER_TYPE::INC_EULER:
      genericSolver = createFlowSolver(SUB_SOLVER_TYPE::INC_EULER, solver, geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::MULTIGRID;
      break;
    case SUB_SOLVER_TYPE::INC_NAVIER_STOKES:
      genericSolver = createFlowSolver(SUB_SOLVER_TYPE::INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::MULTIGRID;
      break;
    case SUB_SOLVER_TYPE::FEA:
      genericSolver = new CFEASolver(geometry, config);
      metaData.integrationType = INTEGRATION_TYPE::STRUCTURAL;
      break;
    case SUB_SOLVER_TYPE::MESH:
      genericSolver = createMeshSolver(solver, geometry, config, iMGLevel, false);
      metaData.integrationType = INTEGRATION_TYPE::NONE;
      break;
    case SUB_SOLVER_TYPE::DG_EULER:
      genericSolver = createDGSolver(SUB_SOLVER_TYPE::DG_EULER, geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::FEM_DG;
      break;
    case SUB_SOLVER_TYPE::DG_NAVIER_STOKES:
      genericSolver = createDGSolver(SUB_SOLVER_TYPE::DG_NAVIER_STOKES, geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::FEM_DG;
      break;
    case SUB_SOLVER_TYPE::HEAT:
      genericSolver = createHeatSolver(solver, geometry, config, iMGLevel, false);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::DISC_ADJ_HEAT:
      genericSolver = createHeatSolver(solver, geometry, config, iMGLevel, true);
      metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
      break;
    case SUB_SOLVER_TYPE::TRANSITION:
      genericSolver = new CTransLMSolver(geometry, config, iMGLevel);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    case SUB_SOLVER_TYPE::TURB: case SUB_SOLVER_TYPE::TURB_SA: case SUB_SOLVER_TYPE::TURB_SST:
      genericSolver = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, false);
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

CSolver* CSolverFactory::createTurbSolver(ENUM_TURB_MODEL kindTurbModel, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, int adjoint){

  CSolver *turbSolver = nullptr;

  if (!adjoint){
    switch (kindTurbModel) {
      case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
        turbSolver = new CTurbSASolver(geometry, config, iMGLevel, solver[FLOW_SOL]->GetFluidModel());
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        turbSolver->Postprocessing(geometry, solver, config, iMGLevel);
        break;
      case SST: case SST_SUST:
        turbSolver = new CTurbSSTSolver(geometry, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        turbSolver->Postprocessing(geometry, solver, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        break;
      default:
        SU2_MPI::Error("Unknown turbulence model", CURRENT_FUNCTION);
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

CSolver* CSolverFactory::createHeatSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint){

  CSolver *heatSolver = nullptr;

  bool standalone = (config->GetKind_Solver() == HEAT_EQUATION) ||
                    (config->GetKind_Solver() == DISC_ADJ_HEAT);

  /*--- Only allocate a heat solver if it should run standalone
   * or if the weakly coupled heat solver is enabled and no energy equation is included ---*/

  if ((config->GetWeakly_Coupled_Heat() && !config->GetEnergy_Equation()) || standalone){
    if (adjoint){
      if (config->GetDiscrete_Adjoint()){
        heatSolver = new CDiscAdjSolver(geometry, config, solver[HEAT_SOL], RUNTIME_HEAT_SYS, iMGLevel);
      } else {
        SU2_MPI::Error("No continuous adjoint heat solver available.", CURRENT_FUNCTION);
      }
    } else {
      heatSolver = new CHeatSolver(geometry, config, iMGLevel);
    }
  }
  return heatSolver;

}

CSolver* CSolverFactory::createMeshSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint){

  CSolver *meshSolver = nullptr;

  if (config->GetDeform_Mesh() && iMGLevel == MESH_0){
    if (!adjoint){
      meshSolver = new CMeshSolver(geometry, config);
    }
    if (adjoint && config->GetDiscrete_Adjoint()){
      meshSolver = new CDiscAdjMeshSolver(geometry, config, solver[MESH_SOL]);
    }
  }
  return meshSolver;

}

CSolver* CSolverFactory::createDGSolver(SUB_SOLVER_TYPE kindDGSolver, CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *DGSolver = nullptr;

  switch (kindDGSolver) {
    case SUB_SOLVER_TYPE::DG_EULER:
      if (config->GetKind_FEM_DG_Shock() == PERSSON){
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

CSolver* CSolverFactory::createFlowSolver(SUB_SOLVER_TYPE kindFlowSolver, CSolver **solver,  CGeometry *geometry, CConfig *config, int iMGLevel){

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
    default:
      SU2_MPI::Error("Flow solver not found", CURRENT_FUNCTION);
      break;
  }

  return flowSolver;

}
