/*!
 * \file CFEAIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/iteration/CFEAIteration.hpp"
#include "../../include/output/COutput.hpp"

void CFEAIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                            CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                            CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                            unsigned short val_iInst) {
  bool StopCalc = false;
  unsigned long IntIter = 0;

  const unsigned long TimeIter = config[val_iZone]->GetTimeIter();
  const unsigned long nIncrements = config[val_iZone]->GetNumberIncrements();

  const bool nonlinear = (config[val_iZone]->GetGeometricConditions() == STRUCT_DEFORMATION::LARGE);
  const bool linear = (config[val_iZone]->GetGeometricConditions() == STRUCT_DEFORMATION::SMALL);
  const bool heat = config[val_iZone]->GetWeakly_Coupled_Heat();
  const bool disc_adj_fem = config[val_iZone]->GetDiscrete_Adjoint();

  /*--- Loads applied in steps (not used for discrete adjoint). ---*/
  const bool incremental_load = config[val_iZone]->GetIncrementalLoad() && !disc_adj_fem;

  /*--- We need to restore the current inner iter in discrete adjoint cases because file output depends on it. ---*/
  const auto CurIter = config[val_iZone]->GetInnerIter();

  CGeometry** geo = geometry[val_iZone][val_iInst];
  CIntegration* feaIntegration = integration[val_iZone][val_iInst][FEA_SOL];
  CSolver* feaSolver = solver[val_iZone][val_iInst][MESH_0][FEA_SOL];

  const auto nPoint = geo[MESH_0]->GetnPoint();
  const auto nDim = geo[MESH_0]->GetnDim();

  su2activematrix Coords;
  if (nonlinear && heat && Coords.empty()) {
    /*--- Save the mesh coordinates to solve the heat equations in the deformed configuration
     * and then restore the coordinates for the FEA solver. We could use CoordsOld of CPoint but
     * that would require accounting for a lot of physics-specific logic in a geometry class. ---*/
    Coords = geo[MESH_0]->nodes->GetCoord();
  }

  auto IterateHeat = [&]() {
    /*--- Update the FVM mesh for the heat solver based on the structural deformations. ---*/
    if (nonlinear) {
      SU2_OMP_PARALLEL_(for schedule(static))
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
        for (auto iDim = 0u; iDim < nDim; ++iDim) {
          geo[MESH_0]->nodes->SetCoord(iPoint, iDim, Coords(iPoint, iDim) + feaSolver->GetNodes()->GetSolution(iPoint, iDim));
        }
      }
      END_SU2_OMP_PARALLEL

      CGeometry::UpdateGeometry(geo, config[val_iZone]);
    }

    config[val_iZone]->SetGlobalParam(MAIN_SOLVER::HEAT_EQUATION, RUNTIME_HEAT_SYS);
    integration[val_iZone][val_iInst][HEAT_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                      RUNTIME_HEAT_SYS, val_iZone, val_iInst);

    /*--- Restore the coordinates for the FEA solver. ---*/
    if (nonlinear) {
      SU2_OMP_PARALLEL_(for schedule(static))
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
        for (auto iDim = 0u; iDim < nDim; ++iDim) {
          geo[MESH_0]->nodes->SetCoord(iPoint, iDim, Coords(iPoint, iDim));
        }
      }
      END_SU2_OMP_PARALLEL
    }
  };

  auto Iterate = [&](unsigned long IntIter) {
    config[val_iZone]->SetInnerIter(IntIter);

    /*--- Heat transfer equations. ---*/
    if (heat) IterateHeat();

    /*--- FEA equations. ---*/
    config[val_iZone]->SetGlobalParam(MAIN_SOLVER::FEM_ELASTICITY, RUNTIME_FEA_SYS);
    feaIntegration->Structural_Iteration(geometry, solver, numerics, config, RUNTIME_FEA_SYS, val_iZone, val_iInst);

    if (!disc_adj_fem) {
      StopCalc = Monitor(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement,
                         FFDBox, val_iZone, INST_0);
    }
  };


  if (linear || (nonlinear && !incremental_load)) {
    /*--- THIS IS THE DIRECT APPROACH (NO INCREMENTAL LOAD APPLIED) ---*/
    /*--- Newton-Raphson subiterations ---*/

    for (IntIter = 0; IntIter < config[val_iZone]->GetnInner_Iter(); IntIter++) {
      Iterate(IntIter);

      /*--- Limit to only one iteration for the discrete adjoint recording, restore inner iter (see above) ---*/
      if (disc_adj_fem) {
        config[val_iZone]->SetInnerIter(CurIter);
        break;
      }
      /*--- Linear elasticity without thermal effects only needs one iteration. ---*/
      if (linear && !heat) {
        output->SetConvergence(true);
        break;
      }
      /*--- Normal stopping criteria. ---*/
      if (StopCalc && IntIter > 0) break;
    }
    return;
  }

  /*--- THIS IS THE INCREMENTAL LOAD APPROACH (only makes sense for nonlinear) ---*/

  /*--- Assume the initial load increment as 1.0 ---*/
  feaSolver->SetLoad_Increment(0, 1.0);
  feaSolver->SetForceCoeff(1.0);

  /*--- Run two nonlinear iterations to check if incremental loading can be skipped ---*/
  for (IntIter = 0; IntIter < 2; ++IntIter) {
    Iterate(IntIter);
  }

  /*--- Early return if we already meet the convergence criteria. ---*/
  if (StopCalc) return;

  /*--- Check user-defined criteria to either increment loads or continue with NR iterations. ---*/
  bool meetCriteria = true;
  for (int i = 0; i < 3; ++i) {
    meetCriteria &= (log10(feaSolver->GetRes_FEM(i)) < config[val_iZone]->GetIncLoad_Criteria(i));
  }

  /*--- If the criteria is met, i.e. the load is not too large, continue the regular calculation. ---*/
  if (meetCriteria) {
    for (IntIter = 2; IntIter < config[val_iZone]->GetnInner_Iter(); IntIter++) {
      Iterate(IntIter);
      if (StopCalc) break;
    }
    return;
  }

  /*--- If the criteria is not met, a whole set of subiterations for the different loads must be done.
   * Restore solution to initial. Because we ramp the load from zero, in multizone cases it is not
   * adequate to take "old values" as those will be for maximum loading on the previous outer iteration. ---*/

  feaSolver->SetInitialCondition(geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone],
                                 TimeIter);

  /*--- For the number of increments ---*/
  for (auto iIncrement = 1ul; iIncrement <= nIncrements; iIncrement++) {
    /*--- Set the load increment and the initial condition, and output the
     * parameters of UTOL, RTOL, ETOL for the previous iteration ---*/

    su2double loadIncrement = su2double(iIncrement) / nIncrements;
    feaSolver->SetLoad_Increment(iIncrement, loadIncrement);

    /*--- Set the convergence monitor to false, to force the solver to converge every subiteration. ---*/
    output->SetConvergence(false);

    if (rank == MASTER_NODE) cout << "\nIncremental load: increment " << iIncrement << endl;

    /*--- Newton-Raphson subiterations ---*/

    for (IntIter = 0; IntIter < config[val_iZone]->GetnInner_Iter(); IntIter++) {
      Iterate(IntIter);
      if (StopCalc && IntIter > 0) break;
    }
  }
  /*--- Just to be sure, set default increment settings. ---*/
  feaSolver->SetLoad_Increment(0, 1.0);
}

void CFEAIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                           CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                           CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                           unsigned short val_iInst) {
  CSolver* feaSolver = solver[val_iZone][val_iInst][MESH_0][FEA_SOL];

  /*----------------- Update structural solver ----------------------*/

  if (config[val_iZone]->GetTime_Domain()) {
    integration[val_iZone][val_iInst][FEA_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][MESH_0], feaSolver,
                                                                   config[val_iZone], MESH_0);
  } else if (config[val_iZone]->GetFSI_Simulation() && config[val_iZone]->GetRelaxation()) {
    /*--- For FSI problems with relaxation, output the relaxed result, which is the one transferred into the fluid
     * domain (for consistent restart purposes). ---*/
    if (config[val_iZone]->GetKind_TimeIntScheme_FEA() == STRUCT_TIME_INT::NEWMARK_IMPLICIT) {
      feaSolver->ImplicitNewmark_Relaxation(geometry[val_iZone][val_iInst][MESH_0], config[val_iZone]);
    }
  }
}

void CFEAIteration::Predictor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                              CSolver***** solver, CNumerics****** numerics, CConfig** config,
                              CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                              CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  CSolver* feaSolver = solver[val_iZone][val_iInst][MESH_0][FEA_SOL];

  feaSolver->PredictStruct_Displacement(geometry[val_iZone][val_iInst][MESH_0], config[val_iZone]);
}

void CFEAIteration::Relaxation(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                               CSolver***** solver, CNumerics****** numerics, CConfig** config,
                               CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                               CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  CSolver* feaSolver = solver[val_iZone][val_iInst][MESH_0][FEA_SOL];

  /*-------------------- Aitken's relaxation ------------------------*/

  /*------------------- Compute the coefficient ---------------------*/

  feaSolver->ComputeAitken_Coefficient(geometry[val_iZone][val_iInst][MESH_0], config[val_iZone],
                                       config[val_iZone]->GetOuterIter());

  /*----------------- Set the relaxation parameter ------------------*/

  feaSolver->SetAitken_Relaxation(geometry[val_iZone][val_iInst][MESH_0], config[val_iZone]);
}

bool CFEAIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                            CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                            CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                            unsigned short val_iInst) {
  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime - StartTime;

  output->SetHistoryOutput(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0],
                           config[val_iZone], config[val_iZone]->GetTimeIter(), config[val_iZone]->GetOuterIter(),
                           config[val_iZone]->GetInnerIter());

  return output->GetConvergence();
}

void CFEAIteration::Solve(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                          CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                          CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                          unsigned short val_iInst) {
  /*------------------ Structural subiteration ----------------------*/
  Iterate(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox, val_iZone,
          val_iInst);

  if (multizone && !config[val_iZone]->GetTime_Domain()) {
    Output(output, geometry, solver, config, config[val_iZone]->GetOuterIter(), false, val_iZone, val_iInst);
  }
}
