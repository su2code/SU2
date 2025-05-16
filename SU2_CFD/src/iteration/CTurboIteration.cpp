/*!
 * \file CTurboIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../include/output/CTurboOutput.hpp"

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
    InitTurboPerformance(geometry[val_iZone][INST_0][MESH_0], config,
                         solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetFluidModel());
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

void CTurboIteration::InitTurboPerformance(CGeometry* geometry, CConfig** config, CFluidModel* fluid) {
  TurbomachineryPerformance = std::make_shared<CTurboOutput>(config, *geometry, *fluid);
  TurbomachineryStagePerformance = std::make_shared<CTurbomachineryStagePerformance>(*fluid);
}

void CTurboIteration::ComputeTurboPerformance(CSolver***** solver, CGeometry**** geometry_container, CConfig** config_container) {
  unsigned short nDim = geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
  unsigned short nBladesRow = config_container[ZONE_0]->GetnMarker_Turbomachinery();
  unsigned short nZone = config_container[ZONE_0]->GetnZone();
  unsigned short iBlade=0, iSpan;
  vector<su2double> TurboPrimitiveIn, TurboPrimitiveOut;
  std::vector<std::vector<CTurbomachineryCombinedPrimitiveStates>> bladesPrimitives;

  if (rank == MASTER_NODE) {
      for (iBlade = 0; iBlade < nBladesRow; iBlade++){
      /* Blade Primitive initialized per blade */
      std::vector<CTurbomachineryCombinedPrimitiveStates> bladePrimitives;
      auto nSpan = config_container[iBlade]->GetnSpanWiseSections();
      for (iSpan = 0; iSpan < nSpan + 1; iSpan++) {
        TurboPrimitiveIn= solver[iBlade][INST_0][MESH_0][FLOW_SOL]->GetTurboPrimitive(iBlade, iSpan, true);
        TurboPrimitiveOut= solver[iBlade][INST_0][MESH_0][FLOW_SOL]->GetTurboPrimitive(iBlade, iSpan, false);
        auto spanInletPrimitive = CTurbomachineryPrimitiveState(TurboPrimitiveIn, nDim, geometry_container[iBlade][INST_0][MESH_0]->GetTangGridVelIn(iBlade, iSpan));
        auto spanOutletPrimitive = CTurbomachineryPrimitiveState(TurboPrimitiveOut, nDim, geometry_container[iBlade][INST_0][MESH_0]->GetTangGridVelOut(iBlade, iSpan));
        auto spanCombinedPrimitive = CTurbomachineryCombinedPrimitiveStates(spanInletPrimitive, spanOutletPrimitive);
        bladePrimitives.push_back(spanCombinedPrimitive);
      }
      bladesPrimitives.push_back(bladePrimitives);
    }
    TurbomachineryPerformance->ComputeTurbomachineryPerformance(bladesPrimitives);
    
    auto nSpan = config_container[ZONE_0]->GetnSpanWiseSections();
    auto InState = TurbomachineryPerformance->GetBladesPerformances().at(ZONE_0).at(nSpan)->GetInletState();
    nSpan = config_container[nZone-1]->GetnSpanWiseSections();
    auto OutState =  TurbomachineryPerformance->GetBladesPerformances().at(nZone-1).at(nSpan)->GetOutletState();
    
    TurbomachineryStagePerformance->ComputePerformanceStage(InState, OutState, config_container[nZone-1]);

    /*--- Set turbomachinery objective function value in each zone ---*/
    for (auto iBlade = 0u; iBlade < nBladesRow; iBlade++) {
      // Should we set in ZONE_0 or nZone-1?
      auto iBladePerf = TurbomachineryPerformance->GetBladesPerformances().at(iBlade).at(nSpan);
      InState = iBladePerf->GetInletState();
      OutState = iBladePerf->GetOutletState();
      solver[iBlade][INST_0][MESH_0][FLOW_SOL]->SetTurboObjectiveFunction(ENUM_OBJECTIVE::ENTROPY_GENERATION, iBlade, TurbomachineryStagePerformance->GetNormEntropyGen()*100);
      solver[iBlade][INST_0][MESH_0][FLOW_SOL]->SetTurboObjectiveFunction(ENUM_OBJECTIVE::TOTAL_PRESSURE_LOSS, iBlade, iBladePerf->GetTotalPressureLoss());
      solver[iBlade][INST_0][MESH_0][FLOW_SOL]->SetTurboObjectiveFunction(ENUM_OBJECTIVE::KINETIC_ENERGY_LOSS, iBlade, iBladePerf->GetKineticEnergyLoss());
    }
    /*--- Set global turbomachinery objective function (evaluated in final zone as dependent on values from all previous zones ) ---*/
    solver[nBladesRow-1][INST_0][MESH_0][FLOW_SOL]->SetTurboObjectiveFunction(ENUM_OBJECTIVE::ENTROPY_GENERATION, nBladesRow, TurbomachineryStagePerformance->GetNormEntropyGen()*100);
    solver[nBladesRow-1][INST_0][MESH_0][FLOW_SOL]->SetTurboObjectiveFunction(ENUM_OBJECTIVE::TOTAL_PRESSURE_LOSS, nBladesRow, TurbomachineryStagePerformance->GetTotalPressureLoss());
    solver[nBladesRow-1][INST_0][MESH_0][FLOW_SOL]->SetTurboObjectiveFunction(ENUM_OBJECTIVE::KINETIC_ENERGY_LOSS, nBladesRow, TurbomachineryStagePerformance->GetKineticEnergyLoss());
  }
}
