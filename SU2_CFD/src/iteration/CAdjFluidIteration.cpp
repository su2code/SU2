/*!
 * \file CAdjFluidIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
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

#include "../../include/iteration/CAdjFluidIteration.hpp"
#include "../../include/output/COutput.hpp"

void CAdjFluidIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                    CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                    CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                    CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  unsigned short iMesh;
  bool harmonic_balance = (config[ZONE_0]->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
  bool dynamic_mesh = config[ZONE_0]->GetGrid_Movement();
  unsigned long InnerIter = 0;
  unsigned long TimeIter = config[ZONE_0]->GetTimeIter();

  /*--- For the unsteady adjoint, load a new direct solution from a restart file. ---*/

  if (((dynamic_mesh && TimeIter == 0) || (config[val_iZone]->GetTime_Marching() != TIME_MARCHING::STEADY)) && !harmonic_balance) {
    int Direct_Iter = SU2_TYPE::Int(config[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(TimeIter) - 1;
    if (rank == MASTER_NODE && val_iZone == ZONE_0 && (config[val_iZone]->GetTime_Marching() != TIME_MARCHING::STEADY))
      cout << endl << " Loading flow solution from direct iteration " << Direct_Iter << "." << endl;
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->LoadRestart(
        geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], Direct_Iter, true);
  }

  /*--- Continuous adjoint Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations ---*/

  if ((InnerIter == 0) || (config[val_iZone]->GetTime_Marching() != TIME_MARCHING::STEADY)) {
    const auto kind_solver = config[val_iZone]->GetKind_Solver();
    switch (kind_solver) {
      case MAIN_SOLVER::ADJ_EULER:
      case MAIN_SOLVER::ADJ_NAVIER_STOKES:
      case MAIN_SOLVER::ADJ_RANS:
        config[val_iZone]->SetGlobalParam(kind_solver, RUNTIME_FLOW_SYS);
        break;
      default:
        break;
    }

    /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/

    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "Begin direct solver to store flow data (single iteration)." << endl;

    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "Compute residuals to check the convergence of the direct problem." << endl;

    integration[val_iZone][val_iInst][FLOW_SOL]->MultiGrid_Iteration(geometry, solver, numerics, config,
                                                                     RUNTIME_FLOW_SYS, val_iZone, val_iInst);

    if (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) {
      /*--- Solve the turbulence model ---*/

      config[val_iZone]->SetGlobalParam(MAIN_SOLVER::ADJ_RANS, RUNTIME_TURB_SYS);
      integration[val_iZone][val_iInst][TURB_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                        RUNTIME_TURB_SYS, val_iZone, val_iInst);

      /*--- Solve transition model ---*/

      if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
        config[val_iZone]->SetGlobalParam(MAIN_SOLVER::ADJ_RANS, RUNTIME_TRANS_SYS);
        integration[val_iZone][val_iInst][TRANS_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                           RUNTIME_TRANS_SYS, val_iZone, val_iInst);
      }
    }

    /*--- Output the residual (visualization purpouses to identify if
     the direct solution is converged)---*/
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "log10[Maximum residual]: " << log10(solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_Max(0))
           << ", located at point " << solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetPoint_Max(0) << "." << endl;

    /*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
     note that in the direct Euler problem we are not computing the gradients of the primitive variables ---*/

    if (config[val_iZone]->GetKind_Gradient_Method() == GREEN_GAUSS)
      solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->SetPrimitive_Gradient_GG(geometry[val_iZone][val_iInst][MESH_0],
                                                                               config[val_iZone]);
    if (config[val_iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
      solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->SetPrimitive_Gradient_LS(geometry[val_iZone][val_iInst][MESH_0],
                                                                               config[val_iZone]);

    /*--- Set contribution from cost function for boundary conditions ---*/

    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      /*--- Set the value of the non-dimensional coefficients in the coarse levels, using the fine level solution ---*/

      solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->SetTotal_CD(
          solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CD());
      solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->SetTotal_CL(
          solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CL());
      solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->SetTotal_CT(
          solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CT());
      solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->SetTotal_CQ(
          solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CQ());

      /*--- Compute the adjoint boundary condition on Euler walls ---*/

      solver[val_iZone][val_iInst][iMesh][ADJFLOW_SOL]->SetForceProj_Vector(
          geometry[val_iZone][val_iInst][iMesh], solver[val_iZone][val_iInst][iMesh], config[val_iZone]);
    }

    if (rank == MASTER_NODE && val_iZone == ZONE_0) cout << "End direct solver, begin adjoint problem." << endl;
  }
}
void CAdjFluidIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                 CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                 CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                 CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  const auto kind_solver = config[val_iZone]->GetKind_Solver();
  switch (kind_solver) {
    case MAIN_SOLVER::ADJ_EULER:
    case MAIN_SOLVER::ADJ_NAVIER_STOKES:
    case MAIN_SOLVER::ADJ_RANS:
      config[val_iZone]->SetGlobalParam(kind_solver, RUNTIME_ADJFLOW_SYS);
      break;
    default:
      break;
  }

  /*--- Iteration of the flow adjoint problem ---*/

  integration[val_iZone][val_iInst][ADJFLOW_SOL]->MultiGrid_Iteration(geometry, solver, numerics, config,
                                                                      RUNTIME_ADJFLOW_SYS, val_iZone, val_iInst);

  /*--- Iteration of the turbulence model adjoint ---*/

  if ((config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) && (!config[val_iZone]->GetFrozen_Visc_Cont())) {
    /*--- Adjoint turbulence model solution ---*/

    config[val_iZone]->SetGlobalParam(MAIN_SOLVER::ADJ_RANS, RUNTIME_ADJTURB_SYS);
    integration[val_iZone][val_iInst][ADJTURB_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_ADJTURB_SYS, val_iZone, val_iInst);
  }
}
void CAdjFluidIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  /*--- Dual time stepping strategy ---*/

  if ((config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
      (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) {
    /*--- Update dual time solver ---*/

    for (unsigned short iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      integration[val_iZone][val_iInst][ADJFLOW_SOL]->SetDualTime_Solver(
          geometry[val_iZone][val_iInst][iMesh], solver[val_iZone][val_iInst][iMesh][ADJFLOW_SOL], config[val_iZone],
          iMesh);

      integration[val_iZone][val_iInst][ADJFLOW_SOL]->SetDualTime_Geometry(
          geometry[val_iZone][val_iInst][iMesh], solver[val_iZone][val_iInst][iMesh][MESH_SOL], config[val_iZone],
          iMesh);
    }
  }
}
