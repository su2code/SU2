/*!
 * \file output_adj_flow_inc.cpp
 * \brief Main subroutines for flow discrete adjoint output
 * \author R. Sanchez
 * \version 7.0.5 "Blackbird"
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


#include "../../include/output/CAdjFlowIncOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CAdjFlowIncOutput::CAdjFlowIncOutput(CConfig *config, unsigned short nDim) :
  COutput(config, nDim, false, false, moduleManagerPtr(new CModuleManager<Modules, Modifiers>(config, nDim))) {

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    requestedHistoryFields.emplace_back("SENSITIVITY");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }

  if (nRequestedScreenFields == 0){
    if (config->GetTime_Domain()) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_ADJ_PRESSURE");
    requestedScreenFields.emplace_back("RMS_ADJ_VELOCITY_X");
    requestedScreenFields.emplace_back("SENS_GEO");
    requestedScreenFields.emplace_back("SENS_AOA");
    nRequestedScreenFields = requestedScreenFields.size();
  }

  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("SENSITIVITY");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  if (find(requestedVolumeFields.begin(), requestedVolumeFields.end(), string("SENSITIVITY")) == requestedVolumeFields.end()) {
    requestedVolumeFields.emplace_back("SENSITIVITY");
    nRequestedVolumeFields ++;
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Adj. Incomp. Fluid)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetAdj_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfAdjCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_AdjFileName();

  /*--- Add the obj. function extension --- */

  restartFilename = config->GetObjFunc_Extension(restartFilename);

}

CAdjFlowIncOutput::~CAdjFlowIncOutput(void) {}


void CAdjFlowIncOutputModule::DefineHistoryFields(CHistoryOutFieldManager& historyFields){

  historyFields.AddField("RMS_ADJ_PRESSURE",    "rms[A_P]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Pressure.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ADJ_VELOCITY_X", "rms[A_U]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Velocity x-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ADJ_VELOCITY_Y", "rms[A_V]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Velocity y-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ADJ_VELOCITY_Z", "rms[A_W]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Velocity z-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ADJ_TEMPERATURE", "rms[A_T]", ScreenOutputFormat::FIXED, "RMS_RES", " Root-mean square residual of the adjoint temperature.", FieldType::RESIDUAL);
  if (rad_model != NONE)
    historyFields.AddField("RMS_ADJ_RAD_ENERGY", "rms[A_P1]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the P1 radiative energy.",FieldType::RESIDUAL);

  historyFields.AddField("MAX_ADJ_PRESSURE",    "max[A_Rho]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Pressure.", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ADJ_VELOCITY_X", "max[A_RhoU]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Velocity x-component", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ADJ_VELOCITY_Y", "max[A_RhoV]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Velocity y-component", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ADJ_VELOCITY_Z", "max[A_RhoW]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Velocity z-component", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ADJ_TEMPERATURE", "max[A_T]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the temperature.", FieldType::RESIDUAL);

  historyFields.AddField("BGS_ADJ_PRESSURE",    "bgs[A_Rho]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Pressure.", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ADJ_VELOCITY_X", "bgs[A_RhoU]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Velocity x-component", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ADJ_VELOCITY_Y", "bgs[A_RhoV]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Velocity y-component", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ADJ_VELOCITY_Z", "bgs[A_RhoW]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Velocity z-component", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ADJ_TEMPERATURE", "bgs[A_T]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint temperature.", FieldType::RESIDUAL);
  if (rad_model != NONE)
    historyFields.AddField("BGS_ADJ_RAD_ENERGY", "bgs[A_P1]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual  of the P1 radiative energy.",FieldType::RESIDUAL);

  historyFields.AddField("SENS_GEO",   "Sens_Geo",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.", FieldType::DEFAULT);
  historyFields.AddField("SENS_PRESS", "Sens_Press", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the far-field pressure.", FieldType::DEFAULT);
  historyFields.AddField("SENS_TEMP",  "Sens_Temp",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", " Sensitivity of the objective function with respect to the far-field temperature.", FieldType::DEFAULT);
  historyFields.AddField("SENS_VEL_IN", "Sens_Vin", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", " Sensitivity of the objective function with respect to the inlet velocity.", FieldType::DEFAULT);
  historyFields.AddField("SENS_PRESS_OUT",  "Sens_Pout",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the outlet pressure.", FieldType::DEFAULT);

}

void CAdjFlowIncOutputModule::LoadHistoryData(CHistoryOutFieldManager &historyFields, const SolverData &solverData, const IterationInfo &iterationInfo){

  const auto* adjflow_solver = solverData.solver[ADJFLOW_SOL];
  const auto* adjheat_solver = solverData.solver[ADJHEAT_SOL];
  const auto* adjrad_solver = solverData.solver[ADJRAD_SOL];

  const auto* config = solverData.config;

  historyFields.SetFieldValue("RMS_ADJ_PRESSURE", log10(adjflow_solver->GetRes_RMS(0)));
  historyFields.SetFieldValue("RMS_ADJ_VELOCITY_X", log10(adjflow_solver->GetRes_RMS(1)));
  historyFields.SetFieldValue("RMS_ADJ_VELOCITY_Y", log10(adjflow_solver->GetRes_RMS(2)));
  if (nDim == 3) {
    historyFields.SetFieldValue("RMS_ADJ_VELOCITY_Z", log10(adjflow_solver->GetRes_RMS(3)));
  }
  if (weakly_coupled_heat){
    historyFields.SetFieldValue("RMS_ADJ_TEMPERATURE",         log10(adjheat_solver->GetRes_RMS(0)));
  }
  if (heat){
    if (nDim == 3) historyFields.SetFieldValue("RMS_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_RMS(4)));
    else           historyFields.SetFieldValue("RMS_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_RMS(3)));
  }
  if (rad_model != NONE)
    historyFields.SetFieldValue("RMS_ADJ_RAD_ENERGY", log10(adjrad_solver->GetRes_RMS(0)));


  historyFields.SetFieldValue("MAX_ADJ_PRESSURE", log10(adjflow_solver->GetRes_RMS(0)));
  historyFields.SetFieldValue("MAX_ADJ_VELOCITY_X", log10(adjflow_solver->GetRes_RMS(1)));
  historyFields.SetFieldValue("MAX_ADJ_VELOCITY_Y", log10(adjflow_solver->GetRes_RMS(2)));
  if (nDim == 3) {
    historyFields.SetFieldValue("MAX_ADJ_VELOCITY_Z", log10(adjflow_solver->GetRes_RMS(3)));
  }
  if (weakly_coupled_heat){
    historyFields.SetFieldValue("MAX_ADJ_TEMPERATURE",         log10(adjheat_solver->GetRes_RMS(0)));
  }
  if (heat){
    if (nDim == 3) historyFields.SetFieldValue("MAX_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_RMS(4)));
    else           historyFields.SetFieldValue("MAX_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_RMS(3)));
  }

  if (config->GetMultizone_Problem()){
    historyFields.SetFieldValue("BGS_ADJ_PRESSURE", log10(adjflow_solver->GetRes_RMS(0)));
    historyFields.SetFieldValue("BGS_ADJ_VELOCITY_X", log10(adjflow_solver->GetRes_RMS(1)));
    historyFields.SetFieldValue("BGS_ADJ_VELOCITY_Y", log10(adjflow_solver->GetRes_RMS(2)));
    if (nDim == 3) {
      historyFields.SetFieldValue("BGS_ADJ_VELOCITY_Z", log10(adjflow_solver->GetRes_RMS(3)));
    }
    if (weakly_coupled_heat){
      historyFields.SetFieldValue("BGS_ADJ_TEMPERATURE",         log10(adjheat_solver->GetRes_RMS(0)));
    }
    if (heat){
      if (nDim == 3) historyFields.SetFieldValue("BGS_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_RMS(4)));
      else           historyFields.SetFieldValue("BGS_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_RMS(3)));
    }
    if (rad_model != NONE)
      historyFields.SetFieldValue("BGS_ADJ_RAD_ENERGY", log10(adjrad_solver->GetRes_RMS(0)));
  }


  historyFields.SetFieldValue("SENS_GEO", adjflow_solver->GetTotal_Sens_Geo());
  historyFields.SetFieldValue("SENS_PRESS", adjflow_solver->GetTotal_Sens_Press());
  historyFields.SetFieldValue("SENS_TEMP", adjflow_solver->GetTotal_Sens_Temp());
  historyFields.SetFieldValue("SENS_VEL_IN", adjflow_solver->GetTotal_Sens_ModVel());
  historyFields.SetFieldValue("SENS_PRESS_OUT", adjflow_solver->GetTotal_Sens_BPress());
}

void CAdjFlowIncOutputModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){

  volumeFields.AddField("ADJ_PRESS",    "Adjoint_Pressure",    "SOLUTION", "Adjoint pressure", FieldType::DEFAULT);
  volumeFields.AddField("ADJ_VEL_X", "Adjoint_Velocity_x", "SOLUTION", "x-component of the adjoint velocity vector", FieldType::DEFAULT);
  volumeFields.AddField("ADJ_VEL_Y", "Adjoint_Velocity_y", "SOLUTION", "y-component of the adjoint velocity vector", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("ADJ_VEL_Z", "Adjoint_Velocity_z", "SOLUTION", "z-component of the adjoint velocity vector", FieldType::DEFAULT);
  volumeFields.AddField("ADJ_TEMP", "Adjoint_Temperature", "SOLUTION",  "Adjoint temperature", FieldType::DEFAULT);

  if (rad_model != NONE)
    volumeFields.AddField("ADJ_P1_EN",  "Adjoint_Energy(P1)", "SOLUTION", "Adjoint radiative energy", FieldType::DEFAULT);

  volumeFields.AddField("RES_ADJ_PRESS",    "Residual_Adjoint_Pressure",    "RESIDUAL", "Residual of the adjoint pressure", FieldType::DEFAULT);
  volumeFields.AddField("RES_ADJ_VEL_X", "Residual_Adjoint_Velocity_x", "RESIDUAL", "Residual of the adjoint x-velocity", FieldType::DEFAULT);
  volumeFields.AddField("RES_ADJ_VEL_Y", "Residual_Adjoint_Velocity_y", "RESIDUAL", "Residual of the adjoint y-velocity", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("RES_ADJ_VEL_Z", "Residual_Adjoint_Velocity_z", "RESIDUAL", "Residual of the adjoint z-velocity", FieldType::DEFAULT);
  volumeFields.AddField("RES_ADJ_TEMP", "Residual_Adjoint_Heat", "RESIDUAL", "Residual of the adjoint temperature", FieldType::DEFAULT);
  if (rad_model != NONE)
    volumeFields.AddField("RES_ADJ_P1_EN",  "Residual_Adjoint_Energy_P1", "RESIDUAL", "Residual of adjoint radiative energy", FieldType::DEFAULT);

  volumeFields.AddField("SENSITIVITY_X", "Sensitivity_x", "SENSITIVITY", "x-component of the sensitivity vector", FieldType::DEFAULT);
  volumeFields.AddField("SENSITIVITY_Y", "Sensitivity_y", "SENSITIVITY", "y-component of the sensitivity vector", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("SENSITIVITY_Z", "Sensitivity_z", "SENSITIVITY", "z-component of the sensitivity vector", FieldType::DEFAULT);
  volumeFields.AddField("SENSITIVITY", "Surface_Sensitivity", "SENSITIVITY", "sensitivity in normal direction", FieldType::SURFACE_INTEGRATE);

}

void CAdjFlowIncOutputModule::LoadVolumeData(CVolumeOutFieldManager &volumeFields, const SolverData &solverData, const IterationInfo &iterationInfo, const PointInfo &pointInfo){

  const auto* Node_AdjFlow = solverData.solver[ADJFLOW_SOL]->GetNodes();
  const auto* Node_AdjHeat = solverData.solver[ADJHEAT_SOL]->GetNodes();
  const auto* Node_AdjRad  = solverData.solver[ADJRAD_SOL]->GetNodes();
  const auto iPoint        = pointInfo.iPoint;

  volumeFields.SetFieldValue("ADJ_PRESS",    Node_AdjFlow->GetSolution(iPoint, 0));
  volumeFields.SetFieldValue("ADJ_VEL_X", Node_AdjFlow->GetSolution(iPoint, 1));
  volumeFields.SetFieldValue("ADJ_VEL_Y", Node_AdjFlow->GetSolution(iPoint, 2));
  if (nDim == 3){
    volumeFields.SetFieldValue("ADJ_VEL_Z", Node_AdjFlow->GetSolution(iPoint, 3));
  }

  if (weakly_coupled_heat){
    volumeFields.SetFieldValue("ADJ_TEMP", Node_AdjHeat->GetSolution(iPoint, 0));
  }
  else {
    if (nDim == 3) volumeFields.SetFieldValue("ADJ_TEMP", Node_AdjFlow->GetSolution(iPoint, 4));
    else           volumeFields.SetFieldValue("ADJ_TEMP", Node_AdjFlow->GetSolution(iPoint, 3));
  }
  if (rad_model != NONE)
    volumeFields.SetFieldValue("ADJ_P1_ENERGY", Node_AdjRad->GetSolution(iPoint, 0));

  volumeFields.SetFieldValue("RES_ADJ_PRESS", Node_AdjFlow->GetSolution(iPoint, 0) - Node_AdjFlow->GetSolution_Old(iPoint, 0));
  volumeFields.SetFieldValue("RES_ADJ_VEL_X", Node_AdjFlow->GetSolution(iPoint, 1) - Node_AdjFlow->GetSolution_Old(iPoint, 1));
  volumeFields.SetFieldValue("RES_ADJ_VEL_Y", Node_AdjFlow->GetSolution(iPoint, 2) - Node_AdjFlow->GetSolution_Old(iPoint, 2));
  if (nDim == 3){
    volumeFields.SetFieldValue("RES_ADJ_VEL_Z", Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
  }

  if (weakly_coupled_heat){
    volumeFields.SetFieldValue("RES_ADJ_TEMP", Node_AdjHeat->GetSolution(iPoint, 0) - Node_AdjHeat->GetSolution_Old(iPoint, 0));
  }
  else {
    if (nDim == 3) volumeFields.SetFieldValue("RES_ADJ_TEMP", Node_AdjFlow->GetSolution(iPoint, 4) - Node_AdjFlow->GetSolution_Old(iPoint, 4));
    else           volumeFields.SetFieldValue("RES_ADJ_TEMP", Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
  }
  if (rad_model != NONE)
    volumeFields.SetFieldValue("RES_ADJ_P1_EN", Node_AdjRad->GetSolution(iPoint, 0) - Node_AdjRad->GetSolution_Old(iPoint, 0));

  volumeFields.SetFieldValue("SENSITIVITY_X", Node_AdjFlow->GetSensitivity(iPoint, 0));
  volumeFields.SetFieldValue("SENSITIVITY_Y", Node_AdjFlow->GetSensitivity(iPoint, 1));
  if (nDim == 3)
    volumeFields.SetFieldValue("SENSITIVITY_Z", Node_AdjFlow->GetSensitivity(iPoint, 2));
}

void CAdjFlowIncOutputModule::LoadSurfaceData(CVolumeOutFieldManager &volumeFields, const SolverData &solverData, const IterationInfo &iterationInfo, const PointInfo &pointInfo){

  const auto SensX = volumeFields.GetFieldValue("SENSITIVITY_X");
  const auto SensY = volumeFields.GetFieldValue("SENSITIVITY_Y");
  const auto SensZ = nDim == 3 ? volumeFields.GetFieldValue("SENSITIVITY_Z") : 0.0;
  const auto NormalX = volumeFields.GetFieldValue("NORMAL_X");
  const auto NormalY = volumeFields.GetFieldValue("NORMAL_Y");
  const auto NormalZ = nDim == 3 ? volumeFields.GetFieldValue("NORMAL_Z") : 0.0;

  volumeFields.SetFieldValue("SENSITIVITY", SensX*NormalX + SensY*NormalY + SensZ*NormalZ);

}


bool CAdjFlowIncOutput::SetInit_Residuals(CConfig *config){

  return (config->GetTime_Marching() != STEADY && (curInnerIter == 0))||
        (config->GetTime_Marching() == STEADY && (curTimeIter < 2));

}

bool CAdjFlowIncOutput::SetUpdate_Averages(CConfig *config){
  return false;

//  return (config->GetUnsteady_Simulation() != STEADY && !dualtime);

}

