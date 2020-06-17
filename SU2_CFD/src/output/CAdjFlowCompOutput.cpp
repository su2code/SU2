/*!
 * \file output_adj_flow_comp.cpp
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


#include "../../include/output/CAdjFlowOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CAdjFlowCompOutput::CAdjFlowCompOutput(CConfig *config, unsigned short nDim) :
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
    requestedScreenFields.emplace_back("RMS_ADJ_DENSITY");
    requestedScreenFields.emplace_back("RMS_ADJ_MOMENTUM_X");
    requestedScreenFields.emplace_back("SENS_MACH");
    requestedScreenFields.emplace_back("SENS_PRESS");
    requestedScreenFields.emplace_back("SENS_TEMP");
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
  ss << "Zone " << config->GetiZone() << " (Adj. Comp. Fluid)";
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

CAdjFlowCompOutput::~CAdjFlowCompOutput(void) {}

void CAdjFlowCompOutputModule::DefineHistoryFields(CHistoryOutFieldManager& historyFields){

  historyFields.AddField("RMS_ADJ_DENSITY",    "rms[A_Rho]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint density.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ADJ_MOMENTUM_X", "rms[A_RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum x-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ADJ_MOMENTUM_Y", "rms[A_RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum y-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ADJ_MOMENTUM_Z", "rms[A_RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum z-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ADJ_ENERGY",     "rms[A_E]",    ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint energy.", FieldType::RESIDUAL);

  historyFields.AddField("MAX_ADJ_DENSITY",    "max[A_Rho]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint density.", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ADJ_MOMENTUM_X", "max[A_RhoU]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum x-component", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ADJ_MOMENTUM_Y", "max[A_RhoV]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum y-component", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ADJ_MOMENTUM_Z", "max[A_RhoW]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum z-component", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ADJ_ENERGY",     "max[A_E]",    ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint energy.", FieldType::RESIDUAL);

  historyFields.AddField("BGS_ADJ_DENSITY",    "bgs[A_Rho]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint density.", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ADJ_MOMENTUM_X", "bgs[A_RhoU]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint momentum x-component", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ADJ_MOMENTUM_Y", "bgs[A_RhoV]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint momentum y-component", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ADJ_MOMENTUM_Z", "bgs[A_RhoW]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint momentum z-component", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ADJ_ENERGY",     "bgs[A_E]",    ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint energy.", FieldType::RESIDUAL);

  historyFields.AddField("SENS_GEO",   "Sens_Geo",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.", FieldType::DEFAULT);
  historyFields.AddField("SENS_MACH",  "Sens_Mach",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the Mach number (only of compressible solver).", FieldType::DEFAULT);
  historyFields.AddField("SENS_PRESS", "Sens_Press", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the far-field pressure.", FieldType::DEFAULT);
  historyFields.AddField("SENS_TEMP",  "Sens_Temp",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the far-field temperature.", FieldType::DEFAULT);
  historyFields.AddField("SENS_AOA",   "Sens_AoA",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the angle of attack (only for compressible solver).", FieldType::DEFAULT);

}

void CAdjFlowCompOutputModule::LoadHistoryData(CHistoryOutFieldManager &historyFields, const SolverData &solverData, const IterationInfo &iterationInfo){

  const auto* adjflow_solver = solverData.solver[ADJFLOW_SOL];
  const auto* config = solverData.config;

  historyFields.SetFieldValue("RMS_ADJ_DENSITY", log10(adjflow_solver->GetRes_RMS(0)));
  historyFields.SetFieldValue("RMS_ADJ_MOMENTUM_X", log10(adjflow_solver->GetRes_RMS(1)));
  historyFields.SetFieldValue("RMS_ADJ_MOMENTUM_Y", log10(adjflow_solver->GetRes_RMS(2)));
  if (nDim == 3) {
    historyFields.SetFieldValue("RMS_ADJ_MOMENTUM_Z", log10(adjflow_solver->GetRes_RMS(3)));
    historyFields.SetFieldValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(4)));
  } else {
    historyFields.SetFieldValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(3)));
  }

  historyFields.SetFieldValue("MAX_ADJ_DENSITY", log10(adjflow_solver->GetRes_Max(0)));
  historyFields.SetFieldValue("MAX_ADJ_MOMENTUM_X", log10(adjflow_solver->GetRes_Max(1)));
  historyFields.SetFieldValue("MAX_ADJ_MOMENTUM_Y", log10(adjflow_solver->GetRes_Max(2)));
  if (nDim == 3) {
    historyFields.SetFieldValue("MAX_ADJ_MOMENTUM_Z", log10(adjflow_solver->GetRes_Max(3)));
    historyFields.SetFieldValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(4)));
  } else {
    historyFields.SetFieldValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(3)));
  }

  if (config->GetMultizone_Problem()){
    historyFields.SetFieldValue("BGS_ADJ_DENSITY", log10(adjflow_solver->GetRes_BGS(0)));
    historyFields.SetFieldValue("BGS_ADJ_MOMENTUM_X", log10(adjflow_solver->GetRes_BGS(1)));
    historyFields.SetFieldValue("BGS_ADJ_MOMENTUM_Y", log10(adjflow_solver->GetRes_BGS(2)));
    if (nDim == 3) {
      historyFields.SetFieldValue("BGS_ADJ_MOMENTUM_Z", log10(adjflow_solver->GetRes_BGS(3)));
      historyFields.SetFieldValue("BGS_ADJ_ENERGY", log10(adjflow_solver->GetRes_BGS(4)));
    } else {
      historyFields.SetFieldValue("BGS_ADJ_ENERGY", log10(adjflow_solver->GetRes_BGS(3)));
    }

  }

  historyFields.SetFieldValue("SENS_GEO", adjflow_solver->GetTotal_Sens_Geo());
  historyFields.SetFieldValue("SENS_AOA", adjflow_solver->GetTotal_Sens_AoA() * PI_NUMBER / 180.0);
  historyFields.SetFieldValue("SENS_MACH", adjflow_solver->GetTotal_Sens_Mach());
  historyFields.SetFieldValue("SENS_PRESS", adjflow_solver->GetTotal_Sens_Press());
  historyFields.SetFieldValue("SENS_TEMP", adjflow_solver->GetTotal_Sens_Temp());
}

void CAdjFlowCompOutputModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){

  volumeFields.AddField("ADJ_DENSITY",    "Adjoint_Density",    "SOLUTION", "Adjoint density", FieldType::DEFAULT);
  volumeFields.AddField("ADJ_MOM_X", "Adjoint_Momentum_x", "SOLUTION", "x-component of the adjoint momentum vector", FieldType::DEFAULT);
  volumeFields.AddField("ADJ_MOM_Y", "Adjoint_Momentum_y", "SOLUTION", "y-component of the adjoint momentum vector", FieldType::DEFAULT);
  volumeFields.AddField("ADJ_MOM_Z", "Adjoint_Momentum_z", "SOLUTION", "z-component of the adjoint momentum vector", FieldType::DEFAULT);
  volumeFields.AddField("ADJ_ENERGY", "Adjoint_Energy", "SOLUTION", "Adjoint energy", FieldType::DEFAULT);

  volumeFields.AddField("RES_ADJ_DENSITY",    "Residual_Adjoint_Density",    "RESIDUAL", "Residual of the adjoint density", FieldType::DEFAULT);
  volumeFields.AddField("RES_ADJ_MOM_X", "Residual_Adjoint_Momentum_x", "RESIDUAL", "Residual of the adjoint x-momentum", FieldType::DEFAULT);
  volumeFields.AddField("RES_ADJ_MOM_Y", "Residual_Adjoint_Momentum_y", "RESIDUAL", "Residual of the adjoint y-momentum", FieldType::DEFAULT);
  volumeFields.AddField("RES_ADJ_MOM_Z", "Residual_Adjoint_Momentum_z", "RESIDUAL", "Residual of the adjoint z-momentum", FieldType::DEFAULT);
  volumeFields.AddField("RES_ADJ_ENERGY", "Residual_Adjoint_Energy", "RESIDUAL", "Residual of the adjoint energy", FieldType::DEFAULT);

  volumeFields.AddField("SENSITIVITY_X", "Sensitivity_x", "SENSITIVITY", "x-component of the sensitivity vector", FieldType::DEFAULT);
  volumeFields.AddField("SENSITIVITY_Y", "Sensitivity_y", "SENSITIVITY", "y-component of the sensitivity vector", FieldType::DEFAULT);
  volumeFields.AddField("SENSITIVITY_Z", "Sensitivity_z", "SENSITIVITY", "z-component of the sensitivity vector", FieldType::DEFAULT);
  volumeFields.AddField("SENSITIVITY", "Surface_Sensitivity", "SENSITIVITY", "sensitivity in normal direction", FieldType::SURFACE_INTEGRATE);

}

void CAdjFlowCompOutputModule::LoadVolumeData(CVolumeOutFieldManager &volumeFields, const SolverData &solverData, const IterationInfo &iterationInfo, const PointInfo &pointInfo){

  const auto* Node_AdjFlow = solverData.solver[ADJFLOW_SOL]->GetNodes();
  const auto iPoint        = pointInfo.iPoint;

  volumeFields.SetFieldValue("ADJ_DENSITY",    Node_AdjFlow->GetSolution(iPoint, 0));
  volumeFields.SetFieldValue("ADJ_MOM_X", Node_AdjFlow->GetSolution(iPoint, 1));
  volumeFields.SetFieldValue("ADJ_MOM_Y", Node_AdjFlow->GetSolution(iPoint, 2));
  if (nDim == 3){
    volumeFields.SetFieldValue("ADJ_MOM_Z", Node_AdjFlow->GetSolution(iPoint, 3));
    volumeFields.SetFieldValue("ADJ_ENERGY",      Node_AdjFlow->GetSolution(iPoint, 4));
  } else {
    volumeFields.SetFieldValue("ADJ_ENERGY",     Node_AdjFlow->GetSolution(iPoint, 3));
  }

  // Residuals
  volumeFields.SetFieldValue("RES_ADJ_DENSITY",    Node_AdjFlow->GetSolution(iPoint, 0) - Node_AdjFlow->GetSolution_Old(iPoint, 0));
  volumeFields.SetFieldValue("RES_ADJ_MOM_X", Node_AdjFlow->GetSolution(iPoint, 1) - Node_AdjFlow->GetSolution_Old(iPoint, 1));
  volumeFields.SetFieldValue("RES_ADJ_MOM_Y", Node_AdjFlow->GetSolution(iPoint, 2) - Node_AdjFlow->GetSolution_Old(iPoint, 2));
  if (nDim == 3){
    volumeFields.SetFieldValue("RES_ADJ_MOM_Z", Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
    volumeFields.SetFieldValue("RES_ADJ_ENERGY",  Node_AdjFlow->GetSolution(iPoint, 4) - Node_AdjFlow->GetSolution_Old(iPoint, 4));
  } else {
    volumeFields.SetFieldValue("RES_ADJ_ENERGY", Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
  }

  volumeFields.SetFieldValue("SENSITIVITY_X", Node_AdjFlow->GetSensitivity(iPoint, 0));
  volumeFields.SetFieldValue("SENSITIVITY_Y", Node_AdjFlow->GetSensitivity(iPoint, 1));
  if (nDim == 3)
    volumeFields.SetFieldValue("SENSITIVITY_Z", Node_AdjFlow->GetSensitivity(iPoint, 2));
}

void CAdjFlowCompOutputModule::LoadSurfaceData(CVolumeOutFieldManager &volumeFields, const SolverData &solverData, const IterationInfo &iterationInfo, const PointInfo &pointInfo){

  const auto SensX = volumeFields.GetFieldValue("SENSITIVITY_X");
  const auto SensY = volumeFields.GetFieldValue("SENSITIVITY_Y");
  const auto SensZ = nDim == 3 ? volumeFields.GetFieldValue("SENSITIVITY_Z") : 0.0;
  const auto NormalX = volumeFields.GetFieldValue("NORMAL_X");
  const auto NormalY = volumeFields.GetFieldValue("NORMAL_Y");
  const auto NormalZ = nDim == 3 ? volumeFields.GetFieldValue("NORMAL_Z") : 0.0;

  volumeFields.SetFieldValue("SENSITIVITY", SensX*NormalX + SensY*NormalY + SensZ*NormalZ);

}

bool CAdjFlowCompOutput::SetInit_Residuals(CConfig *config){

  return ((config->GetTime_Marching() != STEADY) && (curInnerIter == 0)) ||
         ((config->GetTime_Marching() == STEADY) && (curInnerIter < 2));

}

bool CAdjFlowCompOutput::SetUpdate_Averages(CConfig *config){
  return false;

//  return (config->GetUnsteady_Simulation() != STEADY && !dualtime);

}

