/*!
 * \file output_heat.cpp
 * \brief Main subroutines for the heat solver output
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


#include "../../include/output/CHeatOutput.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CHeatOutput::CHeatOutput(CConfig *config, unsigned short nDim) :
  COutput(config, nDim, false, true, moduleManagerPtr(new CModuleManager<Modules, Modifiers>(config, nDim))) {

  multiZone = config->GetMultizone_Problem();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_TEMPERATURE");
    nRequestedScreenFields = requestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Solid Heat)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetVolume_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName();

}

void CHeatOutputModule::DefineHistoryFields(CHistoryOutFieldManager &historyFields){

  historyFields.AddField("RMS_TEMPERATURE", "rms[T]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the temperature", FieldType::RESIDUAL);
  historyFields.AddField("MAX_TEMPERATURE", "max[T]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the temperature", FieldType::RESIDUAL);
  historyFields.AddField("BGS_TEMPERATURE", "bgs[T]", ScreenOutputFormat::FIXED, "BGS_RES", "Block-Gauss seidel residual of the temperature", FieldType::RESIDUAL);

}

void CHeatOutputModule::LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                        const IterationInfo&){

  const auto* config = solverData.config;
  const auto* heat_solver = solverData.solver[HEAT_SOL];

  historyFields.SetFieldValue("AVG_TEMPERATURE", heat_solver->GetTotal_AvgTemperature());
  historyFields.SetFieldValue("RMS_TEMPERATURE", log10(heat_solver->GetRes_RMS(0)));
  historyFields.SetFieldValue("MAX_TEMPERATURE", log10(heat_solver->GetRes_Max(0)));
  if (config->GetMultizone_Problem())
    historyFields.SetFieldValue("BGS_TEMPERATURE", log10(heat_solver->GetRes_BGS(0)));

}

void CHeatOutputModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){

  // SOLUTION
  volumeFields.AddField("TEMPERATURE", "Temperature", "SOLUTION", "Temperature", FieldType::DEFAULT);

  // Primitives
  volumeFields.AddField("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heatflux", FieldType::DEFAULT);

}

void CHeatOutputModule::LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                       const IterationInfo&, const PointInfo& pointInfo){

  auto* heat_solver = solverData.solver[HEAT_SOL];
  auto iPoint = pointInfo.iPoint;

  volumeFields.SetFieldValue("TEMPERATURE", heat_solver->GetNodes()->GetSolution(iPoint,0));

}

void CHeatOutputModule::LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                        const IterationInfo&, const PointInfo& pointInfo){

  const auto iPoint  = pointInfo.iPoint;
  const auto iVertex = pointInfo.iVertex;
  const auto iMarker = pointInfo.iMarker;

  const auto* config   = solverData.config;
  const auto* geometry = solverData.geometry;
  auto* heat_solver = solverData.solver[HEAT_SOL];

  const auto& Grad_Sol = heat_solver->GetNodes()->GetGradient(iPoint);
  const su2double* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
  su2double Area = 0.0; for (int iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
  su2double UnitNormal[3];
  for (int iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = Normal[iDim]/Area;
  }

  su2double GradTemperature = 0.0;
  for (int iDim = 0; iDim < nDim; iDim++)
    GradTemperature += Grad_Sol[0][iDim]*UnitNormal[iDim];

  const auto thermal_diffusivity = config->GetThermalDiffusivity_Solid();
  const auto RefHeatFlux = config->GetHeat_Flux_Ref();
  const auto Heatflux = thermal_diffusivity*GradTemperature*RefHeatFlux;

  volumeFields.SetFieldValue("HEAT_FLUX", Heatflux);
}

CHeatOutput::~CHeatOutput(void) {}


