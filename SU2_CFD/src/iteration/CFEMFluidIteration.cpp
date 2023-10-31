/*!
 * \file CFEMFluidIteration.cpp
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

#include "../../include/iteration/CFEMFluidIteration.hpp"
#include "../../include/output/COutput.hpp"

void CFEMFluidIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                    CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                    CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                    CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  unsigned long TimeIter = config[ZONE_0]->GetTimeIter();
  const bool restart = (config[ZONE_0]->GetRestart() || config[ZONE_0]->GetRestart_Flow());

  /*--- Set the initial condition if this is not a restart. ---*/
  if (TimeIter == 0 && !restart)
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->SetInitialCondition(
        geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], TimeIter);
}

void CFEMFluidIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                 CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                 CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                 CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  /*--- Update global parameters ---*/

  const auto kind_solver = config[val_iZone]->GetKind_Solver();
  switch (kind_solver) {
    case MAIN_SOLVER::FEM_EULER:
    case MAIN_SOLVER::DISC_ADJ_FEM_EULER:
    case MAIN_SOLVER::FEM_NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_FEM_NS:
    case MAIN_SOLVER::FEM_RANS:
    case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
    case MAIN_SOLVER::FEM_LES:
      config[val_iZone]->SetGlobalParam(kind_solver, RUNTIME_FLOW_SYS);
      break;
    default:
      break;
  }

  /*--- Solve the Euler, Navier-Stokes, RANS or LES equations (one iteration) ---*/

  integration[val_iZone][val_iInst][FLOW_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                    RUNTIME_FLOW_SYS, val_iZone, val_iInst);
}

void CFEMFluidIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {}

void CFEMFluidIteration::Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                     CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                     CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                     CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {}
