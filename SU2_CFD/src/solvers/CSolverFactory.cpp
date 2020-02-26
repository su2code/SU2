/*!
 * \file CSolverFactory.cpp
 * \brief Main subroutines for CSolverFactoryclass.
 * \author T. Albring
 * \version 7.0.1 "Blackbird"
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

CSolver** CSolverFactory::createSolverContainer(ENUM_SOLVER kindSolver, CConfig *config, CGeometry *geometry, int iMGLevel){

  CSolver** solver;

  solver = new CSolver*[MAX_SOLS]();

  const bool allocDirect  = false;
  const bool allocAdjoint = true;
  const bool standAlone   = true;

  ENUM_TURB_MODEL kindTurbModel = static_cast<ENUM_TURB_MODEL>(config->GetKind_Turb_Model());

  switch (kindSolver) {
    case TEMPLATE_SOLVER:
      solver[FLOW_SOL] = createSolver(TEMPLATE_SOLVER, solver, geometry, config, iMGLevel);
      break;
    case INC_EULER:
      solver[FLOW_SOL] = createFlowSolver(INC_EULER, solver, geometry, config, iMGLevel);
      break;
    case EULER:
      solver[FLOW_SOL] = createFlowSolver(EULER, solver, geometry, config, iMGLevel);
      break;
    case INC_NAVIER_STOKES:
      solver[FLOW_SOL] = createFlowSolver(INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL] = createHeatSolver(solver, geometry,config, iMGLevel, allocDirect, !standAlone);
      break;
    case NAVIER_STOKES:
      solver[FLOW_SOL] = createFlowSolver(NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case RANS:
      solver[FLOW_SOL] = createFlowSolver(RANS, solver, geometry, config, iMGLevel);
      solver[TURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      break;
    case INC_RANS:
      solver[FLOW_SOL] = createFlowSolver(INC_RANS, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL] = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect, !standAlone);
      solver[TURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      break;
    case HEAT_EQUATION:
      solver[HEAT_SOL] = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect, standAlone);
      break;
    case ADJ_EULER:
      solver[FLOW_SOL]    = createFlowSolver(EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(ADJ_EULER, solver, geometry, config, iMGLevel);
      break;
    case ADJ_NAVIER_STOKES:
      solver[FLOW_SOL]    = createFlowSolver(NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(ADJ_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case ADJ_RANS:
      solver[FLOW_SOL]    = createFlowSolver(RANS, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(ADJ_RANS, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJTURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocAdjoint);
      break;
    case DISC_ADJ_EULER:
      solver[FLOW_SOL]    = createFlowSolver(EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(DISC_ADJ_EULER, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_NAVIER_STOKES:
      solver[FLOW_SOL]    = createFlowSolver(NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(DISC_ADJ_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_RANS:
      solver[FLOW_SOL]    = createFlowSolver(RANS, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(DISC_ADJ_RANS, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJTURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocAdjoint);
      break;
    case DISC_ADJ_INC_EULER:
      solver[FLOW_SOL]    = createFlowSolver(INC_EULER, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(DISC_ADJ_INC_EULER, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_INC_NAVIER_STOKES:
      solver[FLOW_SOL]       = createFlowSolver(INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL]    = createSolver(DISC_ADJ_INC_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      solver[HEAT_SOL]       = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect, !standAlone);
      solver[ADJHEAT_SOL]    = createHeatSolver(solver, geometry, config, iMGLevel, allocAdjoint, !standAlone);
      break;
    case DISC_ADJ_INC_RANS:
      solver[FLOW_SOL]    = createFlowSolver(INC_RANS, solver, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(DISC_ADJ_INC_RANS, solver, geometry, config, iMGLevel);
      solver[TURB_SOL]    = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJTURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocAdjoint);
      solver[HEAT_SOL]    = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect, !standAlone);
      solver[ADJHEAT_SOL] = createHeatSolver(solver, geometry, config, iMGLevel, allocAdjoint, !standAlone);
      break;
    case DISC_ADJ_HEAT:
      solver[HEAT_SOL]    = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect, standAlone);
      solver[ADJHEAT_SOL] = createHeatSolver(solver, geometry, config, iMGLevel, allocAdjoint, standAlone);
      break;
    case FEM_ELASTICITY:
      solver[FEA_SOL] = createSolver(FEM_ELASTICITY, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_FEM:
      solver[FEA_SOL]    = createSolver(FEM_ELASTICITY, solver, geometry, config, iMGLevel);
      solver[ADJFEA_SOL] = createSolver(DISC_ADJ_FEM, solver, geometry, config, iMGLevel);
      break;
    case FEM_EULER:
      solver[FLOW_SOL] = createDGSolver(FEM_EULER, geometry, config, iMGLevel);
      break;
    case FEM_NAVIER_STOKES: case FEM_LES:
      solver[FLOW_SOL] = createDGSolver(FEM_NAVIER_STOKES, geometry, config, iMGLevel);
      break;
    case FEM_RANS:
      solver[FLOW_SOL] = createDGSolver(FEM_RANS, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_FEM_EULER:
      solver[FLOW_SOL]    = createDGSolver(FEM_EULER, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(DISC_ADJ_EULER, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_FEM_NS:
      solver[FLOW_SOL]    = createDGSolver(FEM_NAVIER_STOKES, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(DISC_ADJ_NAVIER_STOKES, solver, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_FEM_RANS:
      solver[FLOW_SOL]    = createDGSolver(FEM_RANS, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = createSolver(DISC_ADJ_RANS, solver, geometry, config, iMGLevel);
      break;
     default:
      solver = nullptr;
  }

  solver[MESH_SOL]    = createMeshSolver(solver, geometry, config, iMGLevel, allocDirect);
  solver[ADJMESH_SOL] = createMeshSolver(solver, geometry, config, iMGLevel, allocAdjoint);

  return solver;

}

CSolver* CSolverFactory::createSolver(ENUM_SOLVER kindSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *genericSolver = nullptr;

  switch (kindSolver) {
    case TEMPLATE_SOLVER:
      genericSolver = new CTemplateSolver(geometry, config);
      break;
    case ADJ_EULER:
      genericSolver = new CAdjEulerSolver(geometry, config, iMGLevel);
      break;
    case ADJ_NAVIER_STOKES: case ADJ_RANS:
      genericSolver = new CAdjNSSolver(geometry, config, iMGLevel);
      break;
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
      genericSolver = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      break;
    case FEM_ELASTICITY:
      genericSolver = new CFEASolver(geometry, config);
      break;
    case DISC_ADJ_FEM:
      genericSolver = new CDiscAdjFEASolver(geometry, config, solver[FEA_SOL], RUNTIME_FEA_SYS, iMGLevel);
      break;

    default:
      SU2_MPI::Error("Solver not found", CURRENT_FUNCTION);
  }

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

CSolver* CSolverFactory::createHeatSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint, bool standalone){

  CSolver *heatSolver = nullptr;

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

CSolver* CSolverFactory::createDGSolver(ENUM_SOLVER kindDGSolver, CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *DGSolver = nullptr;

  switch (kindDGSolver) {
    case FEM_EULER:
      if (config->GetKind_FEM_DG_Shock() == PERSSON){
        DGSolver = new CFEM_DG_NSSolver(geometry, config, iMGLevel);
      } else {
        DGSolver = new CFEM_DG_EulerSolver(geometry, config, iMGLevel);
      }
      break;
    case FEM_NAVIER_STOKES: case FEM_LES:
      DGSolver = new CFEM_DG_NSSolver(geometry, config, iMGLevel);
      break;
    case FEM_RANS:
      SU2_MPI::Error("DG RANS solver not yet implemented.", CURRENT_FUNCTION);
      break;
    default:
      SU2_MPI::Error("Requested DG solver not found", CURRENT_FUNCTION);
      break;
  }

  return DGSolver;
}

CSolver* CSolverFactory::createFlowSolver(ENUM_SOLVER kindFlowSolver, CSolver **solver,  CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *flowSolver = nullptr;

  switch (kindFlowSolver) {
    case EULER:
      flowSolver = new CEulerSolver(geometry, config, iMGLevel);
      flowSolver->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      break;
    case NAVIER_STOKES: case RANS:
      flowSolver = new CNSSolver(geometry, config, iMGLevel);
      break;
    case INC_EULER:
      flowSolver = new CIncEulerSolver(geometry, config, iMGLevel);
      flowSolver->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      break;
    case INC_NAVIER_STOKES: case INC_RANS:
      flowSolver = new CIncNSSolver(geometry, config, iMGLevel);
      break;
    default:
      SU2_MPI::Error("Flow solver not found", CURRENT_FUNCTION);
      break;
  }

  return flowSolver;

}
