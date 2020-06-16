/*!
 * \file output_flow_inc.cpp
 * \brief Main subroutines for incompressible flow output
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


#include "../../include/output/CFlowIncOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CFlowIncOutput::CFlowIncOutput(CConfig *config, unsigned short nDim) :
  CFlowOutput(config, nDim, false, true, moduleManagerPtr(new CModuleManager<Modules>(config, nDim))){

  turb_model = config->GetKind_Turb_Model();

  heat = config->GetEnergy_Equation();

  weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }

  if (nRequestedScreenFields == 0){
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_PRESSURE");
    requestedScreenFields.emplace_back("RMS_VELOCITY_X");
    requestedScreenFields.emplace_back("RMS_VELOCITY_Y");
    nRequestedScreenFields = requestedScreenFields.size();
  }

  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("PRIMITIVE");
    if (config->GetGrid_Movement()) requestedVolumeFields.emplace_back("GRID_VELOCITY");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Incomp. Fluid)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetVolume_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_PRESSURE");


}

CFlowIncOutput::~CFlowIncOutput(void) {}

void CFlowIncOutputModule::DefineHistoryFields(CHistoryOutFieldManager &historyFields){

  historyFields.AddField("RMS_PRESSURE",    "rms[Rho]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the pressure.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_VELOCITY_X", "rms[RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the velocity x-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_VELOCITY_Y", "rms[RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the velocity y-component.", FieldType::RESIDUAL);
  if (nDim == 3)
    historyFields.AddField("RMS_VELOCITY_Z", "rms[RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the velocity z-component.", FieldType::RESIDUAL);
  if (heat) historyFields.AddField("RMS_TEMPERATURE", "rms[T]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the temperature.", FieldType::RESIDUAL);

  historyFields.AddField("MAX_PRESSURE",    "max[Rho]",  ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the pressure.", FieldType::RESIDUAL);
  historyFields.AddField("MAX_VELOCITY_X", "max[RhoU]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the velocity x-component.", FieldType::RESIDUAL);
  historyFields.AddField("MAX_VELOCITY_Y", "max[RhoV]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the velocity y-component.", FieldType::RESIDUAL);
  if (nDim == 3)
    historyFields.AddField("MAX_VELOCITY_Z", "max[RhoW]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the velocity z-component.", FieldType::RESIDUAL);
  if (heat) historyFields.AddField("MAX_TEMPERATURE", "rms[T]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the temperature.", FieldType::RESIDUAL);

  historyFields.AddField("BGS_PRESSURE",    "bgs[Rho]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the pressure.", FieldType::RESIDUAL);
  historyFields.AddField("BGS_VELOCITY_X", "bgs[RhoU]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the velocity x-component.", FieldType::RESIDUAL);
  historyFields.AddField("BGS_VELOCITY_Y", "bgs[RhoV]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the velocity y-component.", FieldType::RESIDUAL);
  if (nDim == 3)
    historyFields.AddField("BGS_VELOCITY_Z", "bgs[RhoW]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the velocity z-component.", FieldType::RESIDUAL);
  if (heat) historyFields.AddField("BGS_TEMPERATURE", "rms[T]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the temperature.", FieldType::RESIDUAL);


}

void CFlowIncOutputModule::LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                           const IterationInfo&){

  const auto* config = solverData.config;
  auto* flow_solver = solverData.solver[FLOW_SOL];

  historyFields.SetFieldValue("RMS_PRESSURE", log10(flow_solver->GetRes_RMS(0)));
  historyFields.SetFieldValue("RMS_VELOCITY_X", log10(flow_solver->GetRes_RMS(1)));
  historyFields.SetFieldValue("RMS_VELOCITY_Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 3){
    historyFields.SetFieldValue("RMS_VELOCITY_Z", log10(flow_solver->GetRes_RMS(3)));
    if (heat) historyFields.SetFieldValue("RMS_TEMPERATURE", log10(flow_solver->GetRes_RMS(4)));
  } else {
    if (heat) historyFields.SetFieldValue("RMS_TEMPERATURE", log10(flow_solver->GetRes_RMS(3)));
  }

  historyFields.SetFieldValue("MAX_PRESSURE", log10(flow_solver->GetRes_Max(0)));
  historyFields.SetFieldValue("MAX_VELOCITY_X", log10(flow_solver->GetRes_Max(1)));
  historyFields.SetFieldValue("MAX_VELOCITY_Y", log10(flow_solver->GetRes_Max(2)));
  if (nDim == 3){
    historyFields.SetFieldValue("MAX_VELOCITY_Z", log10(flow_solver->GetRes_Max(3)));
    if (heat) historyFields.SetFieldValue("MAX_TEMPERATURE", log10(flow_solver->GetRes_Max(4)));
  } else {
    if (heat) historyFields.SetFieldValue("MAX_TEMPERATURE", log10(flow_solver->GetRes_Max(3)));
  }

  if (config->GetMultizone_Problem()){
    historyFields.SetFieldValue("BGS_PRESSURE", log10(flow_solver->GetRes_BGS(0)));
    historyFields.SetFieldValue("BGS_VELOCITY_X", log10(flow_solver->GetRes_BGS(1)));
    historyFields.SetFieldValue("BGS_VELOCITY_Y", log10(flow_solver->GetRes_BGS(2)));
    if (nDim == 3){
      historyFields.SetFieldValue("BGS_VELOCITY_Z", log10(flow_solver->GetRes_BGS(3)));
      if (heat) historyFields.SetFieldValue("BGS_TEMPERATURE", log10(flow_solver->GetRes_BGS(4)));
    } else {
      if (heat) historyFields.SetFieldValue("BGS_TEMPERATURE", log10(flow_solver->GetRes_BGS(3)));
    }
  }
}

void CFlowIncOutputModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){

  volumeFields.AddField("PRESSURE",    "Density",    "SOLUTION", "Pressure", FieldType::DEFAULT);
  volumeFields.AddField("VELOCITY_X", "Momentum_x", "SOLUTION", "x-component of the velocity vector", FieldType::DEFAULT);
  volumeFields.AddField("VELOCITY_Y", "Momentum_y", "SOLUTION", "y-component of the velocity vector", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("VELOCITY_Z", "Momentum_z", "SOLUTION", "z-component of the velocity vector", FieldType::DEFAULT);
  if (heat) volumeFields.AddField("TEMPERATURE",     "Temperature",     "SOLUTION", "Temperature", FieldType::DEFAULT);

  volumeFields.AddField("DENSITY",    "Pressure",                "PRIMITIVE", "Density", FieldType::DEFAULT);
  volumeFields.AddField("MACH",        "Mach",                    "PRIMITIVE", "Mach number", FieldType::DEFAULT);
  volumeFields.AddField("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient", FieldType::DEFAULT);
  if (viscous){
    volumeFields.AddField("LAMINAR_VISC", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity", FieldType::DEFAULT);
    if (heat) {
      volumeFields.AddField("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heatflux", FieldType::DEFAULT);
      volumeFields.AddField("ENTHALPY", "Enthalpy", "PRIMITIVE", "Enthalpy", FieldType::DEFAULT);
    }
  }
  volumeFields.AddField("TOTAL_PRESS",    "Total Pressure", "TOTAL_QUANTITIES", "Total pressure", FieldType::DEFAULT);
  volumeFields.AddField("TOTAL_TEMP", "Total Temperature", "TOTAL_QUANTITIES", "Total temperature", FieldType::DEFAULT);


}

void CFlowIncOutputModule::LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                           const IterationInfo&, const PointInfo& pointInfo){
  const auto iPoint  = pointInfo.iPoint;
  const auto iVertex = pointInfo.iVertex;
  const auto iMarker = pointInfo.iMarker;

  const auto* config = solverData.config;
  const auto* geometry = solverData.geometry;
  auto* flow_solver = solverData.solver[FLOW_SOL];

  if (viscous && heat){
    const auto& Grad_Sol = flow_solver->GetNodes()->GetGradient_Primitive(iPoint);
    const su2double* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
    const auto RefHeatFlux = config->GetHeat_Flux_Ref();
    su2double Area = 0.0; for (int iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    su2double UnitNormal[3];
    for (int iDim = 0; iDim < nDim; iDim++) {
      UnitNormal[iDim] = Normal[iDim]/Area;
    }
    su2double GradTemperature = 0.0;
    for (int iDim = 0; iDim < nDim; iDim++)
      GradTemperature -= Grad_Sol[nDim+1][iDim]*UnitNormal[iDim];

    const su2double thermal_conductivity = flow_solver->GetNodes()->GetThermalConductivity(iPoint);
    const auto Heatflux = -thermal_conductivity*GradTemperature*RefHeatFlux;

    volumeFields.SetFieldValue("HEAT_FLUX", Heatflux);
  }

}

void CFlowIncOutputModule::LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                          const IterationInfo&, const PointInfo& pointInfo){

  const auto iPoint  = pointInfo.iPoint;
  const auto* config = solverData.config;
  auto* flow_solver = solverData.solver[FLOW_SOL];
  const auto* Node_Flow = flow_solver->GetNodes();

  volumeFields.SetFieldValue("PRESSURE",      Node_Flow->GetSolution(iPoint, 0)*config->GetPressure_Ref());
  volumeFields.SetFieldValue("VELOCITY_X",   Node_Flow->GetSolution(iPoint, 1)*config->GetVelocity_Ref());
  volumeFields.SetFieldValue("VELOCITY_Y",   Node_Flow->GetSolution(iPoint, 2)*config->GetVelocity_Ref());
  if (nDim == 3){
    volumeFields.SetFieldValue("VELOCITY_Z", Node_Flow->GetSolution(iPoint, 3)*config->GetVelocity_Ref());
    if (heat) volumeFields.SetFieldValue("TEMPERATURE",     Node_Flow->GetSolution(iPoint, 4)*config->GetTemperature_Ref());
  } else {
    if (heat) volumeFields.SetFieldValue("TEMPERATURE",     Node_Flow->GetSolution(iPoint, 3)*config->GetTemperature_Ref());
  }

  su2double VelMag = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    VelMag += pow(flow_solver->GetVelocity_Inf(iDim),2.0);
  }
  const su2double factor = 1.0/(0.5*flow_solver->GetDensity_Inf()*VelMag);
  volumeFields.SetFieldValue("PRESSURE_COEFF", (Node_Flow->GetPressure(iPoint) - config->GetPressure_FreeStreamND())*factor);
  volumeFields.SetFieldValue("DENSITY", Node_Flow->GetDensity(iPoint)*config->GetDensity_Ref());

  if (viscous){
    volumeFields.SetFieldValue("LAMINAR_VISC", Node_Flow->GetLaminarViscosity(iPoint)*config->GetViscosity_Ref());
    if (heat) volumeFields.SetFieldValue("ENTHALPY", Node_Flow->GetSpecificHeatCp(iPoint)*Node_Flow->GetTemperature(iPoint)*config->GetEnergy_Ref());
  }

  volumeFields.SetFieldValue("TOTAL_PRESS",  (Node_Flow->GetPressure(iPoint)+
                       0.5*Node_Flow->GetVelocity2(iPoint)*Node_Flow->GetDensity(iPoint))*config->GetPressure_Ref());

  volumeFields.SetFieldValue("TOTAL_TEMP",  (Node_Flow->GetTemperature(iPoint) +
                       0.5*Node_Flow->GetVelocity2(iPoint)/Node_Flow->GetSpecificHeatCp(iPoint))*config->GetTemperature_Ref());
}


bool CFlowIncOutput::SetInit_Residuals(CConfig *config){

  return (config->GetTime_Marching() != STEADY && (curInnerIter == 0))||
        (config->GetTime_Marching() == STEADY && (curInnerIter < 2));

}

bool CFlowIncOutput::SetUpdate_Averages(CConfig *config){

  return (config->GetTime_Marching() != STEADY && (curInnerIter == config->GetnInner_Iter() - 1 || convergence));

}
