/*!
 * \file output_flow_comp.cpp
 * \brief Main subroutines for compressible flow output
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


#include "../../include/output/CFlowCompOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CFlowCompOutput::CFlowCompOutput(CConfig *config, unsigned short nDim) :
  CFlowOutput(config, nDim, false, true, moduleManagerPtr(new CModuleManager<Modules, Modifiers>(config, nDim))){

  turb_model = config->GetKind_Turb_Model();
  lastInnerIter = curInnerIter;
  gridMovement = config->GetGrid_Movement();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    if (config->GetTime_Domain()) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_DENSITY");
    requestedScreenFields.emplace_back("RMS_MOMENTUM_X");
    requestedScreenFields.emplace_back("RMS_MOMENTUM_Y");
    requestedScreenFields.emplace_back("RMS_ENERGY");
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
  ss << "Zone " << config->GetiZone() << " (Comp. Fluid)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetVolume_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName();

//  if (config->GetFixed_CL_Mode()) {
//    bool found = false;
//    for (unsigned short iField = 0; iField < convFields.size(); iField++)
//      if (convFields[iField] == "LIFT") found = true;
//    if (!found) {
//      if (rank == MASTER_NODE)
//        cout<<"  Fixed CL: Adding LIFT as Convergence Field to ensure convergence to target CL"<<endl;
//      convFields.emplace_back("LIFT");
//      newFunc.resize(convFields.size());
//      oldFunc.resize(convFields.size());
//      cauchySerie.resize(convFields.size(), vector<su2double>(nCauchy_Elems, 0.0));
//    }
//  }
}

CFlowCompOutput::~CFlowCompOutput(void) {}

void CFlowCompOutputModule::LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                            const IterationInfo&){

  const auto* config = solverData.config;
  auto* flow_solver = solverData.solver[FLOW_SOL];

  historyFields.SetFieldValue("RMS_DENSITY", log10(flow_solver->GetRes_RMS(0)));
  historyFields.SetFieldValue("RMS_MOMENTUM_X", log10(flow_solver->GetRes_RMS(1)));
  historyFields.SetFieldValue("RMS_MOMENTUM_Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 2)
    historyFields.SetFieldValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(3)));
  else {
    historyFields.SetFieldValue("RMS_MOMENTUM_Z", log10(flow_solver->GetRes_RMS(3)));
    historyFields.SetFieldValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(4)));
  }
  historyFields.SetFieldValue("MAX_DENSITY", log10(flow_solver->GetRes_Max(0)));
  historyFields.SetFieldValue("MAX_MOMENTUM_X", log10(flow_solver->GetRes_Max(1)));
  historyFields.SetFieldValue("MAX_MOMENTUM_Y", log10(flow_solver->GetRes_Max(2)));
  if (nDim == 2)
    historyFields.SetFieldValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(3)));
  else {
    historyFields.SetFieldValue("MAX_MOMENTUM_Z", log10(flow_solver->GetRes_Max(3)));
    historyFields.SetFieldValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(4)));
  }
  if (config->GetMultizone_Problem()){
    historyFields.SetFieldValue("BGS_DENSITY", log10(flow_solver->GetRes_BGS(0)));
    historyFields.SetFieldValue("BGS_MOMENTUM_X", log10(flow_solver->GetRes_BGS(1)));
    historyFields.SetFieldValue("BGS_MOMENTUM_Y", log10(flow_solver->GetRes_BGS(2)));
    if (nDim == 2)
      historyFields.SetFieldValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(3)));
    else {
      historyFields.SetFieldValue("BGS_MOMENTUM_Z", log10(flow_solver->GetRes_BGS(3)));
      historyFields.SetFieldValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(4)));
    }
  }
}

void CFlowCompOutputModule::DefineHistoryFields(CHistoryOutFieldManager& historyFields){

  historyFields.AddField("RMS_DENSITY",    "rms[Rho]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the density.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_MOMENTUM_X", "rms[RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum x-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_MOMENTUM_Y", "rms[RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum y-component.", FieldType::RESIDUAL);
  if (nDim == 3)
    historyFields.AddField("RMS_MOMENTUM_Z", "rms[RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum z-component.", FieldType::RESIDUAL);
  historyFields.AddField("RMS_ENERGY",     "rms[RhoE]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the energy.", FieldType::RESIDUAL);

  historyFields.AddField("MAX_DENSITY",    "max[Rho]",  ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the density.", FieldType::RESIDUAL);
  historyFields.AddField("MAX_MOMENTUM_X", "max[RhoU]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum x-component.", FieldType::RESIDUAL);
  historyFields.AddField("MAX_MOMENTUM_Y", "max[RhoV]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum y-component.", FieldType::RESIDUAL);
  if (nDim == 3)
    historyFields.AddField("MAX_MOMENTUM_Z", "max[RhoW]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the momentum z-component.", FieldType::RESIDUAL);
  historyFields.AddField("MAX_ENERGY",     "max[RhoE]", ScreenOutputFormat::FIXED, "MAX_RES", "Max residual of the energy.", FieldType::RESIDUAL);

  historyFields.AddField("BGS_DENSITY",    "bgs[Rho]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the density.", FieldType::RESIDUAL);
  historyFields.AddField("BGS_MOMENTUM_X", "bgs[RhoU]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the momentum x-component.", FieldType::RESIDUAL);
  historyFields.AddField("BGS_MOMENTUM_Y", "bgs[RhoV]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the momentum y-component.", FieldType::RESIDUAL);
  if (nDim == 3)
    historyFields.AddField("BGS_MOMENTUM_Z", "bgs[RhoW]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the momentum z-component.", FieldType::RESIDUAL);
  historyFields.AddField("BGS_ENERGY",     "bgs[RhoE]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the energy.", FieldType::RESIDUAL);

}

void CFlowCompOutputModule::DefineVolumeFields(CVolumeOutFieldManager& volumeFields) {
  // Solution variables

  volumeFields.AddField("DENSITY",    "Density",    "SOLUTION", "Density", FieldType::DEFAULT);
  volumeFields.AddField("MOMENTUM_X", "Momentum_x", "SOLUTION", "x-component of the momentum vector", FieldType::DEFAULT);
  volumeFields.AddField("MOMENTUM_Y", "Momentum_y", "SOLUTION", "y-component of the momentum vector", FieldType::DEFAULT);
  if (nDim == 3)
    volumeFields.AddField("MOMENTUM_Z", "Momentum_z", "SOLUTION", "z-component of the momentum vector", FieldType::DEFAULT);
  volumeFields.AddField("ENERGY",     "Energy",     "SOLUTION", "Energy", FieldType::DEFAULT);

  volumeFields.AddField("PRESSURE",    "Pressure",                "PRIMITIVE", "Pressure", FieldType::DEFAULT);
  volumeFields.AddField("TEMPERATURE", "Temperature",             "PRIMITIVE", "Temperature", FieldType::DEFAULT);
  volumeFields.AddField("MACH",        "Mach",                    "PRIMITIVE", "Mach number", FieldType::DEFAULT);
  volumeFields.AddField("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient", FieldType::DEFAULT);
  if (viscous){
    volumeFields.AddField("LAMINAR_VISC", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity", FieldType::DEFAULT);
    volumeFields.AddField("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heatflux", FieldType::DEFAULT);
  }
  volumeFields.AddField("TOTAL_PRESS",    "Total Pressure", "TOTAL_QUANTITIES", "Total pressure", FieldType::DEFAULT);
  volumeFields.AddField("TOTAL_TEMP", "Total Temperature", "TOTAL_QUANTITIES", "Total temperature", FieldType::DEFAULT);



}

void CFlowCompOutputModule::LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                           const IterationInfo&, const PointInfo& pointInfo) {

  const auto iPoint  = pointInfo.iPoint;
  const auto* config = solverData.config;
  auto* flow_solver = solverData.solver[FLOW_SOL];
  const auto* Node_Flow = flow_solver->GetNodes();

  const su2double momentumRef = config->GetDensity_Ref()*config->GetVelocity_Ref();

  volumeFields.SetFieldValue("DENSITY",      Node_Flow->GetSolution(iPoint, 0)*config->GetDensity_Ref());
  volumeFields.SetFieldValue("MOMENTUM_X",   Node_Flow->GetSolution(iPoint, 1)*momentumRef);
  volumeFields.SetFieldValue("MOMENTUM_Y",   Node_Flow->GetSolution(iPoint, 2)*momentumRef);
  if (nDim == 3){
    volumeFields.SetFieldValue("MOMENTUM_Z", Node_Flow->GetSolution(iPoint, 3)*momentumRef);
    volumeFields.SetFieldValue("ENERGY",     Node_Flow->GetSolution(iPoint, 4)*config->GetEnergy_Ref());
  } else {
    volumeFields.SetFieldValue("ENERGY",     Node_Flow->GetSolution(iPoint, 3)*config->GetEnergy_Ref());
  }

  const su2double Gamma = config->GetGamma();
  su2double Mach = sqrt(Node_Flow->GetVelocity2(iPoint))/Node_Flow->GetSoundSpeed(iPoint);

  su2double VelMag = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    VelMag += pow(flow_solver->GetVelocity_Inf(iDim),2.0);
  }
  const su2double factor = 1.0/(0.5*flow_solver->GetDensity_Inf()*VelMag);

  volumeFields.SetFieldValue("PRESSURE_COEFF",
                             (Node_Flow->GetPressure(iPoint) - flow_solver->GetPressure_Inf())*factor);
  volumeFields.SetFieldValue("MACH", Mach);
  volumeFields.SetFieldValue("PRESSURE", Node_Flow->GetPressure(iPoint)*config->GetPressure_Ref());
  volumeFields.SetFieldValue("TEMPERATURE", Node_Flow->GetTemperature(iPoint)*config->GetTemperature_Ref());
  if (viscous)
    volumeFields.SetFieldValue("LAMINAR_VISC", Node_Flow->GetLaminarViscosity(iPoint)*config->GetViscosity_Ref());

  volumeFields.SetFieldValue("TOTAL_PRESS", (Node_Flow->GetPressure(iPoint) *
                                          pow( 1.0 + Mach * Mach * 0.5 *
                                               (Gamma - 1.0), Gamma / (Gamma - 1.0)))*config->GetPressure_Ref());

  volumeFields.SetFieldValue("TOTAL_TEMP", (Node_Flow->GetTemperature(iPoint) *
                                                  (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0)))*config->GetTemperature_Ref());
}

void CFlowCompOutputModule::LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                            const IterationInfo&, const PointInfo& pointInfo){
  const auto iPoint  = pointInfo.iPoint;
  const auto iVertex = pointInfo.iVertex;
  const auto iMarker = pointInfo.iMarker;

  const auto* config = solverData.config;
  const auto* geometry = solverData.geometry;
  auto* flow_solver = solverData.solver[FLOW_SOL];

  if (viscous){
    const auto& Grad_Primitive = flow_solver->GetNodes()->GetGradient_Primitive(iPoint);
    const su2double* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
    const auto Gamma = config->GetGamma();
    const auto Gas_Constant = config->GetGas_Constant();
    const su2double Viscosity = flow_solver->GetNodes()->GetLaminarViscosity(iPoint);
    const auto Prandtl_Lam = config->GetPrandtl_Lam();
    const auto RefHeatFlux = config->GetHeat_Flux_Ref();
    su2double Area = 0.0; for (int iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    su2double UnitNormal[3];
    for (int iDim = 0; iDim < nDim; iDim++) {
      UnitNormal[iDim] = Normal[iDim]/Area;
    }
    su2double GradTemperature = 0.0;
    for (int iDim = 0; iDim < nDim; iDim++)
      GradTemperature -= Grad_Primitive[0][iDim]*UnitNormal[iDim];

    const su2double Cp = Gamma / (Gamma - 1) * Gas_Constant;
    const su2double thermal_conductivity = Cp * Viscosity/Prandtl_Lam;
    const auto Heatflux = -thermal_conductivity*GradTemperature*RefHeatFlux;

    volumeFields.SetFieldValue("HEAT_FLUX", Heatflux);
  }
}


bool CFlowCompOutput::SetInit_Residuals(CConfig *config){

  return (config->GetTime_Marching() != STEADY && (curInnerIter == 0))||
        (config->GetTime_Marching() == STEADY && (curInnerIter < 2));

}

bool CFlowCompOutput::SetUpdate_Averages(CConfig *config){

  return (config->GetTime_Marching() != STEADY && (curInnerIter == config->GetnInner_Iter() - 1 || convergence));

}


void CFlowCompOutput::SetAdditionalScreenOutput(CConfig *config){

  if (config->GetFixed_CL_Mode()){
    SetFixedCLScreenOutput(config);
  }
}

void CFlowCompOutput::SetFixedCLScreenOutput(CConfig *config){
  PrintingToolbox::CTablePrinter FixedCLSummary(&cout);

  if (fabs(modules->GetHistoryFields().GetFieldValue("CL_DRIVER_COMMAND")) > 1e-16){
    FixedCLSummary.AddColumn("Fixed CL Mode", 40);
    FixedCLSummary.AddColumn("Value", 30);
    FixedCLSummary.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    FixedCLSummary.PrintHeader();
    FixedCLSummary << "Current CL" << modules->GetHistoryFields().GetFieldValue("LIFT");
    FixedCLSummary << "Target CL" << config->GetTarget_CL();
    FixedCLSummary << "Previous AOA" << modules->GetHistoryFields().GetFieldValue("PREV_AOA");
    if (config->GetFinite_Difference_Mode()){
      FixedCLSummary << "Changed AoA by (Finite Difference step)" <<  modules->GetHistoryFields().GetFieldValue("CL_DRIVER_COMMAND");
      lastInnerIter = curInnerIter - 1;
    }
    else
      FixedCLSummary << "Changed AoA by" <<  modules->GetHistoryFields().GetFieldValue("CL_DRIVER_COMMAND");
    FixedCLSummary.PrintFooter();
    SetScreen_Header(config);
  }

  else if (config->GetFinite_Difference_Mode() &&
           modules->GetHistoryFields().GetFieldValue("AOA") ==  modules->GetHistoryFields().GetFieldValue("PREV_AOA")){
    FixedCLSummary.AddColumn("Fixed CL Mode (Finite Difference)", 40);
    FixedCLSummary.AddColumn("Value", 30);
    FixedCLSummary.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    FixedCLSummary.PrintHeader();
    FixedCLSummary << "Delta CL / Delta AoA" << config->GetdCL_dAlpha();
    FixedCLSummary << "Delta CD / Delta CL" << config->GetdCD_dCL();
    if (nDim == 3){
      FixedCLSummary << "Delta CMx / Delta CL" << config->GetdCMx_dCL();
      FixedCLSummary << "Delta CMy / Delta CL" << config->GetdCMy_dCL();
    }
    FixedCLSummary << "Delta CMz / Delta CL" << config->GetdCMz_dCL();
    FixedCLSummary.PrintFooter();
    curInnerIter = lastInnerIter;
    WriteMetaData(config);
    curInnerIter = config->GetInnerIter();
  }
}

bool CFlowCompOutput::WriteHistoryFile_Output(CConfig *config) {
  return !config->GetFinite_Difference_Mode() && COutput::WriteHistoryFile_Output(config);
}
