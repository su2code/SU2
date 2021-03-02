/*!
 * \file CHeatIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 7.1.0 "Blackbird"
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

      integration[val_iZone][val_iInst][HEAT_SOL]->SetDualTime_Geometry(geometry[val_iZone][val_iInst][iMesh],
                                                                        solver[val_iZone][val_iInst][iMesh][MESH_SOL],
                                                                        config[val_iZone], iMesh);

      integration[val_iZone][val_iInst][HEAT_SOL]->SetConvergence(false);
    }
  }
}
