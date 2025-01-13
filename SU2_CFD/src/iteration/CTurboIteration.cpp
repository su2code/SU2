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

void CTurboIteration::TurboRamp(CGeometry**** geometry_container, CConfig** config_container, unsigned long iter, unsigned short iZone, unsigned short ramp_flag) {

  auto* config = config_container[iZone];
  auto* geometry = geometry_container[iZone][INST_0][ZONE_0];

  auto GetRamp_Coeff = [&](CConfig* &config, unsigned short int x) { 
    if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampRotatingFrame()) return config->GetRampRotatingFrame_Coeff(x);
    else if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampTranslationFrame()) return config->GetRampTranslationFrame_Coeff(x);
    else if (ramp_flag == TURBO_RAMP_TYPE::BOUNDARY && config->GetRampOutletPressure()) return config->GetRampOutletPressure_Coeff(x);
    else if (ramp_flag == TURBO_RAMP_TYPE::BOUNDARY && config->GetRampOutletMassFlow()) return config->GetRampOutletMassFlow_Coeff(x);
    else return config->GetRampRotatingFrame_Coeff(x); //No ramp specified, does nothing
  };

  auto SetRate = [&](CConfig* &config, su2double val) { 
    if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampRotatingFrame()) config->SetRotation_Rate(2, val);
    else if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampTranslationFrame()) config->SetTranslation_Rate(1, val);
  };

  auto SetVelocity = [&](CGeometry* &geometry, CConfig* &config, bool print) { 
    if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampRotatingFrame()) geometry->SetRotationalVelocity(config, print);
    else if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampTranslationFrame()) geometry->SetTranslationalVelocity(config, print);
  };

  auto GetFinalValue = [&](CConfig* &config) { 
    if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampRotatingFrame()) return config->GetFinalRotation_Rate_Z();
    else if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampTranslationFrame()) return config->GetFinalTranslation_Rate_Y();
    else if (ramp_flag == TURBO_RAMP_TYPE::BOUNDARY && config->GetRampOutletPressure()) return config->GetFinalOutletPressure();
    else if (ramp_flag == TURBO_RAMP_TYPE::BOUNDARY && config->GetRampOutletMassFlow()) return config->GetFinalOutletMassFlow();
    else return config->GetFinalRotation_Rate_Z();
  };

  auto msg = "\nUpdated rotating frame grid velocities for zone ";
  if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetRampTranslationFrame()) msg = "\nUpdated translation frame grid velocities for zone ";

  auto SetMonitorValue = [&](CConfig* &config, su2double val) { 
    if (ramp_flag == TURBO_RAMP_TYPE::BOUNDARY && config->GetRampOutletPressure()) config->SetMonitorOutletPressure(val);
    else if (ramp_flag == TURBO_RAMP_TYPE::BOUNDARY && config->GetRampOutletMassFlow()) config->SetMonitorOutletMassFlow(val);
  };

  if (config_container[ZONE_0]->GetMultizone_Problem())
    iter = config_container[ZONE_0]->GetOuterIter();

  /*-- Update grid velocities (ROTATING_FRAME, STEADY_TRANSLATION)*/
  if (ramp_flag == TURBO_RAMP_TYPE::GRID && config->GetGrid_Movement()) {
    const long unsigned rampFreq = SU2_TYPE::Int(GetRamp_Coeff(config, 1));
    const long unsigned finalRamp_Iter = SU2_TYPE::Int(GetRamp_Coeff(config, 2));
    const auto ini_vel = GetRamp_Coeff(config, 0);
    const bool print = (config->GetComm_Level() == COMM_FULL);

    if(iter % rampFreq == 0 && iter <= finalRamp_Iter){   
      const auto final_vel =  GetFinalValue(config);
      if(fabs(final_vel) > 0.0) {
        const auto vel = ini_vel + iter * (final_vel - ini_vel)/finalRamp_Iter;
        SetRate(config, vel);
        if (rank == MASTER_NODE && iter > 0) cout << msg << iZone << ".\n";
        SetVelocity(geometry, config, print);
        geometry->SetShroudVelocity(config);
      }
      geometry->SetAvgTurboValue(config, iZone, INFLOW, false);
      geometry->SetAvgTurboValue(config, iZone, OUTFLOW, false);
      geometry->GatherInOutAverageValues(config, false);

      if (iZone < nZone - 1) {
        geometry_container[nZone-1][INST_0][MESH_0]->SetAvgTurboGeoValues(config ,geometry_container[iZone][INST_0][MESH_0], iZone);
      }
    }
  }

  if (ramp_flag == TURBO_RAMP_TYPE::BOUNDARY){
    const long unsigned rampFreq = SU2_TYPE::Int(GetRamp_Coeff(config, 1));
    const long unsigned finalRamp_Iter = SU2_TYPE::Int(GetRamp_Coeff(config, 2));
    const auto outVal_ini = GetRamp_Coeff(config, 0);
    const auto outVal_final = GetFinalValue(config);

    if (iter % rampFreq == 0 && iter <= finalRamp_Iter) {
      const su2double outVal = outVal_ini + iter * (outVal_final - outVal_ini) / finalRamp_Iter;
      if (rank == MASTER_NODE) SetMonitorValue(config, outVal);

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