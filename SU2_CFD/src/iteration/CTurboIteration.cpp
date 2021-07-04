/*!
 * \file CTurboIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/iteration/CTurboIteration.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/objectives/turbomachinery/CTurbomachineryPerformance.hpp"

void CTurboIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                 CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                 CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                 CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  /*--- Average quantities at the inflow and outflow boundaries ---*/
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->TurboAverageProcess(
      solver[val_iZone][val_iInst][MESH_0], geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], INFLOW);
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->TurboAverageProcess(
      solver[val_iZone][val_iInst][MESH_0], geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], OUTFLOW);

  if (config[val_iZone]->GetBoolTurbomachinery()) {
    InitTurboPerformance(geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetFluidModel());
  }
}

void CTurboIteration::Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                  CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                  CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                  CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  /*--- Average quantities at the inflow and outflow boundaries ---*/
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->TurboAverageProcess(
      solver[val_iZone][val_iInst][MESH_0], geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], INFLOW);
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->TurboAverageProcess(
      solver[val_iZone][val_iInst][MESH_0], geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], OUTFLOW);

  /*--- Gather Inflow and Outflow quantities on the Master Node to compute performance ---*/
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GatherInOutAverageValues(config[val_iZone],
                                                                           geometry[val_iZone][val_iInst][MESH_0]);
}

void CTurboIteration::InitTurboPerformance(CGeometry *geometry, CConfig *config, CFluidModel *fluid){
  TurbomachineryPerformance = std::make_shared<CTurbomachineryPerformance>(*config, *geometry, *fluid);
  TurbomachineryStagePerformance = new CTurbomachineryStagePerformance(*fluid);
}