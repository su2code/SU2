/*!
 * \file CHeatIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/iteration/CHeatIteration.hpp"
#include "../../include/output/COutput.hpp"

void CHeatIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                             CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                             CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                             unsigned short val_iInst) {
  /*--- Update global parameters ---*/

  config[val_iZone]->SetGlobalParam(HEAT_EQUATION, RUNTIME_HEAT_SYS);

  integration[val_iZone][val_iInst][HEAT_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                    RUNTIME_HEAT_SYS, val_iZone, val_iInst);
}

void CHeatIteration::Solve(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                           CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                           CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                           unsigned short val_iInst) {
  /*--- Boolean to determine if we are running a static or dynamic case ---*/
  bool steady = !config[val_iZone]->GetTime_Domain();

  unsigned long Inner_Iter, nInner_Iter = config[val_iZone]->GetnInner_Iter();
  bool StopCalc = false;

  /*--- Synchronization point before a single solver iteration.
        Compute the wall clock time required. ---*/

  StartTime = SU2_MPI::Wtime();

  /*--- For steady-state flow simulations, we need to loop over ExtIter for the number of time steps ---*/
  /*--- However, ExtIter is the number of FSI iterations, so nIntIter is used in this case ---*/

  for (Inner_Iter = 0; Inner_Iter < nInner_Iter; Inner_Iter++) {
    config[val_iZone]->SetInnerIter(Inner_Iter);

    /*--- Run a single iteration of the solver ---*/
    Iterate(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox, val_iZone,
            INST_0);

    /*--- Monitor the pseudo-time ---*/
    StopCalc = Monitor(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox,
                       val_iZone, INST_0);

    /*--- Output files at intermediate iterations if the problem is single zone ---*/

    if (singlezone && steady) {
      Output(output, geometry, solver, config, Inner_Iter, StopCalc, val_iZone, val_iInst);
    }

    /*--- If the iteration has converged, break the loop ---*/
    if (StopCalc) break;
  }

  if (multizone && steady) {
    Output(output, geometry, solver, config, config[val_iZone]->GetOuterIter(), StopCalc, val_iZone, val_iInst);

    /*--- Set the fluid convergence to false (to make sure outer subiterations converge) ---*/

    integration[val_iZone][INST_0][HEAT_SOL]->SetConvergence(false);
  }
}

void CHeatIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                            CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                            CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                            unsigned short val_iInst) {
  unsigned short iMesh;

  /*--- Dual time stepping strategy ---*/

  if ((config[val_iZone]->GetTime_Marching() == DT_STEPPING_1ST) ||
      (config[val_iZone]->GetTime_Marching() == DT_STEPPING_2ND)) {
    /*--- Update dual time solver ---*/

    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      integration[val_iZone][val_iInst][HEAT_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][iMesh],
                                                                      solver[val_iZone][val_iInst][iMesh][HEAT_SOL],
                                                                      config[val_iZone], iMesh);
      integration[val_iZone][val_iInst][HEAT_SOL]->SetConvergence(false);
    }
  }
}
