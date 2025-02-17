/*!
 * \file CTurboIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

void CTurboIteration::TurboRamp(CGeometry**** geometry_container, CConfig** config_container, unsigned long iter, unsigned short iZone, TURBO_RAMP_TYPE ramp_flag) {
  /*--- Generic function for handling turbomachinery ramps ---*/
  // Grid updates (i.e. rotation/translation) handled seperately to boundary (i.e. pressure/mass flow) updates
  auto* config = config_container[iZone];
  auto* geometry = geometry_container[iZone][INST_0][ZONE_0];

  std::string msg = "\nUpdated rotating frame grid velocities for zone ";
  if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetKind_GridMovement() == ENUM_GRIDMOVEMENT::STEADY_TRANSLATION) 
    msg = "\nUpdated translation frame grid velocities for zone ";

  if (config_container[ZONE_0]->GetMultizone_Problem())
    iter = config_container[ZONE_0]->GetOuterIter();

  /*-- Update grid velocities (ROTATING_FRAME, STEADY_TRANSLATION)*/
  if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetGrid_Movement()) {
    const auto ini_vel = config->GetRamp_Coeff(ramp_flag, TURBO_RAMP_COEFF::INITIAL_VALUE);
    const auto unsigned rampFreq = SU2_TYPE::Int(config->GetRamp_Coeff(ramp_flag, TURBO_RAMP_COEFF::UPDATE_FREQ));
    const auto unsigned finalRamp_Iter = SU2_TYPE::Int(config->GetRamp_Coeff(ramp_flag, TURBO_RAMP_COEFF::FINAL_ITER));
    
    // Two options needed as if finalRamp_Iter % rampFreq != 0 final value is not set correctly 
    if((iter % rampFreq == 0 && iter < finalRamp_Iter) || (iter == finalRamp_Iter)){   
      const auto final_vel =  config->GetFinalValue(ramp_flag);
      if(fabs(final_vel) > 0.0) {
        const auto vel = ini_vel + iter * (final_vel - ini_vel)/finalRamp_Iter;
        config->SetRate(vel);
        if (rank == MASTER_NODE && iter > 0) cout << msg << iZone << ".\n";
        geometry->SetVelocity(config, true);
        geometry->SetShroudVelocity(config);
      }
      // Update average turbo values
      geometry->SetAvgTurboValue(config, iZone, INFLOW, false);
      geometry->SetAvgTurboValue(config, iZone, OUTFLOW, false);
      geometry->GatherInOutAverageValues(config, false);

      if (iZone < nZone - 1) {
        geometry_container[nZone-1][INST_0][MESH_0]->SetAvgTurboGeoValues(config ,geometry_container[iZone][INST_0][MESH_0], iZone);
      }
    }
  }

  // Boundary ramps (pressure/mass flow)
  if (ramp_flag == TURBO_RAMP_TYPE::BOUNDARY){
    const auto outVal_ini = config->GetRamp_Coeff(ramp_flag, TURBO_RAMP_COEFF::INITIAL_VALUE);
    const long unsigned rampFreq = SU2_TYPE::Int(config->GetRamp_Coeff(ramp_flag, TURBO_RAMP_COEFF::UPDATE_FREQ));
    const long unsigned finalRamp_Iter = SU2_TYPE::Int(config->GetRamp_Coeff(ramp_flag, TURBO_RAMP_COEFF::FINAL_ITER));
    const auto outVal_final = config->GetFinalValue(ramp_flag);

    if ((iter % rampFreq == 0 && iter < finalRamp_Iter) || (iter == finalRamp_Iter)) {
      const su2double outVal = outVal_ini + iter * (outVal_final - outVal_ini) / finalRamp_Iter;
      if (rank == MASTER_NODE) config->SetMonitorValue(outVal);

      for (auto iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        const auto KindBC = config->GetMarker_All_KindBC(iMarker);
        const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        unsigned short KindBCOption;
        switch (KindBC) {
          case RIEMANN_BOUNDARY:
            KindBCOption = config->GetKind_Data_Riemann(Marker_Tag);
            if (KindBCOption == STATIC_PRESSURE || KindBCOption == RADIAL_EQUILIBRIUM) {
              SU2_MPI::Error("Outlet pressure ramp only implemented for NRBC", CURRENT_FUNCTION);
            }
            break;
          case GILES_BOUNDARY:
            KindBCOption = config->GetKind_Data_Giles(Marker_Tag);
            if (KindBCOption == STATIC_PRESSURE || KindBCOption == STATIC_PRESSURE_1D ||
                KindBCOption == RADIAL_EQUILIBRIUM || KindBCOption == MASS_FLOW_OUTLET) {
              config->SetGiles_Var1(outVal, Marker_Tag);
            }
            break;
        }
      }
    }
  }
}