/*!
 * \file CFlowIncOutput.cpp
 * \brief Main subroutines for incompressible flow output
 * \author R. Sanchez
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


#include "../../include/output/CFlowIncOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CFlowIncOutput::CFlowIncOutput(CConfig *config, unsigned short nDim) : CFlowOutput(config, nDim, false) {

  turb_model = config->GetKind_Turb_Model();

  heat = config->GetEnergy_Equation();

  weakly_coupled_heat = config->GetWeakly_Coupled_Heat();
  flamelet = (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET);

  streamwisePeriodic = (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE);
  streamwisePeriodic_temperature = config->GetStreamwise_Periodic_Temperature();

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
    requestedScreenFields.emplace_back("RMS_VELOCITY-X");
    requestedScreenFields.emplace_back("RMS_VELOCITY-Y");
    nRequestedScreenFields = requestedScreenFields.size();
  }

  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("PRIMITIVE");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  if (gridMovement) {
    auto notFound = requestedVolumeFields.end();
    if (find(requestedVolumeFields.begin(), notFound, string("GRID_VELOCITY")) == notFound) {
      requestedVolumeFields.emplace_back("GRID_VELOCITY");
      nRequestedVolumeFields ++;
    }
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

  if (convFields.empty()) convFields.emplace_back("RMS_PRESSURE");

}

void CFlowIncOutput::SetHistoryOutputFields(CConfig *config){

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables.
  /// DESCRIPTION: Root-mean square residual of the pressure.
  AddHistoryOutput("RMS_PRESSURE",   "rms[P]", ScreenOutputFormat::FIXED,   "RMS_RES", "Root-mean square residual of the pressure.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the velocity x-component.
  AddHistoryOutput("RMS_VELOCITY-X", "rms[U]", ScreenOutputFormat::FIXED,   "RMS_RES", "Root-mean square residual of the velocity x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the velocity y-component.
  AddHistoryOutput("RMS_VELOCITY-Y", "rms[V]", ScreenOutputFormat::FIXED,   "RMS_RES", "Root-mean square residual of the velocity y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the velocity z-component.
  if (nDim == 3) AddHistoryOutput("RMS_VELOCITY-Z", "rms[W]", ScreenOutputFormat::FIXED,   "RMS_RES", "Root-mean square residual of the velocity z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the temperature.
  if (heat || weakly_coupled_heat) AddHistoryOutput("RMS_TEMPERATURE", "rms[T]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the temperature.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_ScalarRMS_RES(config);

  /// DESCRIPTION: Root-mean square residual of the radiative energy (P1 model).
  if (config->AddRadiation()) AddHistoryOutput("RMS_RAD_ENERGY", "rms[E_Rad]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the radiative energy.", HistoryFieldType::RESIDUAL);
  /// END_GROUP

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the pressure.
  AddHistoryOutput("MAX_PRESSURE",   "max[P]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum residual of the pressure.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity x-component.
  AddHistoryOutput("MAX_VELOCITY-X", "max[U]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum residual of the velocity x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity y-component.
  AddHistoryOutput("MAX_VELOCITY-Y", "max[V]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum residual of the velocity y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity z-component.
  if (nDim == 3)
    AddHistoryOutput("MAX_VELOCITY-Z", "max[W]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum residual of the velocity z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the temperature.
  if (heat || weakly_coupled_heat)
    AddHistoryOutput("MAX_TEMPERATURE", "max[T]", ScreenOutputFormat::FIXED, "MAX_RES", "Root-mean square residual of the temperature.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_ScalarMAX_RES(config);
  /// END_GROUP

  /// BEGIN_GROUP: BGS_RES, DESCRIPTION: The block-gauss seidel residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the pressure.
  AddHistoryOutput("BGS_PRESSURE",   "bgs[P]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the pressure.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity x-component.
  AddHistoryOutput("BGS_VELOCITY-X", "bgs[U]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the velocity x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity y-component.
  AddHistoryOutput("BGS_VELOCITY-Y", "bgs[V]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the velocity y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity z-component.
  if (nDim == 3)
    AddHistoryOutput("BGS_VELOCITY-Z", "bgs[W]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the velocity z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the temperature.
  if (heat || weakly_coupled_heat)
    AddHistoryOutput("BGS_TEMPERATURE", "bgs[T]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the temperature.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_ScalarBGS_RES(config);

  /// DESCRIPTION: Multizone residual of the radiative energy (P1 model).
  if (config->AddRadiation()) AddHistoryOutput("BGS_RAD_ENERGY", "bgs[E_Rad]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the radiative energy.", HistoryFieldType::RESIDUAL);
  /// END_GROUP

  /// DESCRIPTION: Angle of attack
  AddHistoryOutput("AOA", "AoA", ScreenOutputFormat::SCIENTIFIC,"AOA", "Angle of attack");
  /// DESCRIPTION: Linear solver iterations
  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");
  AddHistoryOutputFieldsScalarLinsol(config);

  AddHistoryOutput("MIN_DELTA_TIME", "Min DT", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current minimum local time step");
  AddHistoryOutput("MAX_DELTA_TIME", "Max DT", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current maximum local time step");

  AddHistoryOutput("MIN_CFL", "Min CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current minimum of the local CFL numbers");
  AddHistoryOutput("MAX_CFL", "Max CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current maximum of the local CFL numbers");
  AddHistoryOutput("AVG_CFL", "Avg CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current average of the local CFL numbers");

  if (config->GetDeform_Mesh()){
    AddHistoryOutput("DEFORM_MIN_VOLUME", "MinVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Minimum volume in the mesh");
    AddHistoryOutput("DEFORM_MAX_VOLUME", "MaxVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Maximum volume in the mesh");
    AddHistoryOutput("DEFORM_ITER", "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", ScreenOutputFormat::FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");
  }

  if(streamwisePeriodic) {
    AddHistoryOutput("STREAMWISE_MASSFLOW", "SWMassflow", ScreenOutputFormat::FIXED, "STREAMWISE_PERIODIC", "Massflow in streamwise periodic flow");
    AddHistoryOutput("STREAMWISE_DP",       "SWDeltaP",   ScreenOutputFormat::FIXED, "STREAMWISE_PERIODIC", "Pressure drop in streamwise periodic flow");
    AddHistoryOutput("STREAMWISE_HEAT",     "SWHeat",     ScreenOutputFormat::FIXED, "STREAMWISE_PERIODIC", "Integrated heat for streamwise periodic flow");
  }

  AddAnalyzeSurfaceOutput(config);

  AddAerodynamicCoefficients(config);

  if (heat || weakly_coupled_heat) {
    AddHeatCoefficients(config);
    AddHistoryOutput("AVG_TEMPERATURE", "Temp", ScreenOutputFormat::SCIENTIFIC, "HEAT", "Average temperature on all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  }

  AddRotatingFrameCoefficients();

}

void CFlowIncOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  CSolver* flow_solver = solver[FLOW_SOL];
  CSolver* heat_solver = solver[HEAT_SOL];
  CSolver* rad_solver  = solver[RAD_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];

  SetHistoryOutputValue("RMS_PRESSURE", log10(flow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_VELOCITY-X", log10(flow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_VELOCITY-Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 3) SetHistoryOutputValue("RMS_VELOCITY-Z", log10(flow_solver->GetRes_RMS(3)));

  if (config->AddRadiation())
    SetHistoryOutputValue("RMS_RAD_ENERGY", log10(rad_solver->GetRes_RMS(0)));

  SetHistoryOutputValue("MAX_PRESSURE", log10(flow_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_VELOCITY-X", log10(flow_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_VELOCITY-Y", log10(flow_solver->GetRes_Max(2)));
  if (nDim == 3) SetHistoryOutputValue("RMS_VELOCITY-Z", log10(flow_solver->GetRes_Max(3)));

  if (multiZone){
    SetHistoryOutputValue("BGS_PRESSURE", log10(flow_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_VELOCITY-X", log10(flow_solver->GetRes_BGS(1)));
    SetHistoryOutputValue("BGS_VELOCITY-Y", log10(flow_solver->GetRes_BGS(2)));
    if (nDim == 3) SetHistoryOutputValue("BGS_VELOCITY-Z", log10(flow_solver->GetRes_BGS(3)));
    if (config->AddRadiation())
      SetHistoryOutputValue("BGS_RAD_ENERGY", log10(rad_solver->GetRes_BGS(0)));
  }

  if (weakly_coupled_heat){
    SetHeatCoefficients(config, heat_solver);
    SetHistoryOutputValue("AVG_TEMPERATURE", heat_solver->GetTotal_AvgTemperature());
    SetHistoryOutputValue("RMS_TEMPERATURE", log10(heat_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("MAX_TEMPERATURE", log10(heat_solver->GetRes_Max(0)));
    if (multiZone) SetHistoryOutputValue("BGS_TEMPERATURE", log10(heat_solver->GetRes_BGS(0)));
  }
  if (heat) {
    SetHeatCoefficients(config, flow_solver);
    SetHistoryOutputValue("AVG_TEMPERATURE", flow_solver->GetTotal_AvgTemperature());
    SetHistoryOutputValue("RMS_TEMPERATURE", log10(flow_solver->GetRes_RMS(nDim + 1)));
    SetHistoryOutputValue("MAX_TEMPERATURE", log10(flow_solver->GetRes_Max(nDim + 1)));
    if (multiZone) {
      SetHistoryOutputValue("BGS_TEMPERATURE", log10(flow_solver->GetRes_BGS(nDim + 1)));
    }
  }

  SetHistoryOutputValue("LINSOL_ITER", flow_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(flow_solver->GetResLinSolver()));

  if (config->GetDeform_Mesh()){
    SetHistoryOutputValue("DEFORM_MIN_VOLUME", mesh_solver->GetMinimum_Volume());
    SetHistoryOutputValue("DEFORM_MAX_VOLUME", mesh_solver->GetMaximum_Volume());
    SetHistoryOutputValue("DEFORM_ITER", mesh_solver->GetIterLinSolver());
    SetHistoryOutputValue("DEFORM_RESIDUAL", log10(mesh_solver->GetResLinSolver()));
  }

  SetHistoryOutputValue("MIN_DELTA_TIME", flow_solver->GetMin_Delta_Time());
  SetHistoryOutputValue("MAX_DELTA_TIME", flow_solver->GetMax_Delta_Time());

  SetHistoryOutputValue("MIN_CFL", flow_solver->GetMin_CFL_Local());
  SetHistoryOutputValue("MAX_CFL", flow_solver->GetMax_CFL_Local());
  SetHistoryOutputValue("AVG_CFL", flow_solver->GetAvg_CFL_Local());

  LoadHistoryDataScalar(config, solver);

  if(streamwisePeriodic) {
    SetHistoryOutputValue("STREAMWISE_MASSFLOW", flow_solver->GetStreamwisePeriodicValues().Streamwise_Periodic_MassFlow);
    SetHistoryOutputValue("STREAMWISE_DP", flow_solver->GetStreamwisePeriodicValues().Streamwise_Periodic_PressureDrop);
    SetHistoryOutputValue("STREAMWISE_HEAT", flow_solver->GetStreamwisePeriodicValues().Streamwise_Periodic_IntegratedHeatFlow);
  }

  /*--- Set the analyse surface history values --- */

  SetAnalyzeSurface(solver, geometry, config, false);

  /*--- Set aeroydnamic coefficients --- */

  SetAerodynamicCoefficients(config, flow_solver);

  /*--- Set rotating frame coefficients --- */

  SetRotatingFrameCoefficients(flow_solver);

  /*--- Keep this as last, since it uses the history values that were set. ---*/

  SetCustomOutputs(solver, geometry, config);

  SetCustomAndComboObjectives(FLOW_SOL, config, solver);

}

void CFlowIncOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddCoordinates();

  // SOLUTION variables
  AddVolumeOutput("PRESSURE",   "Pressure",   "SOLUTION", "Pressure");
  AddVolumeOutput("VELOCITY-X", "Velocity_x", "SOLUTION", "x-component of the velocity vector");
  AddVolumeOutput("VELOCITY-Y", "Velocity_y", "SOLUTION", "y-component of the velocity vector");
  if (nDim == 3)
    AddVolumeOutput("VELOCITY-Z", "Velocity_z", "SOLUTION", "z-component of the velocity vector");
  if (heat || weakly_coupled_heat || flamelet)
    AddVolumeOutput("TEMPERATURE",  "Temperature","SOLUTION", "Temperature");

  SetVolumeOutputFieldsScalarSolution(config);

  // Radiation variables
  if (config->AddRadiation())
    AddVolumeOutput("P1-RAD", "Radiative_Energy(P1)", "SOLUTION", "Radiative Energy");

  // Grid velocity
  if (gridMovement){
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY", "x-component of the grid velocity vector");
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY", "y-component of the grid velocity vector");
    if (nDim == 3 )
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY", "z-component of the grid velocity vector");
  }

  // Primitive variables
  AddVolumeOutput("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient");
  AddVolumeOutput("DENSITY",        "Density",              "PRIMITIVE", "Density");

  if (config->GetKind_Solver() == MAIN_SOLVER::INC_RANS || config->GetKind_Solver() == MAIN_SOLVER::INC_NAVIER_STOKES){
    AddVolumeOutput("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity");
    AddVolumeOutput("HEAT_CAPACITY", "Heat_Capacity", "PRIMITIVE", "Heat capacity");
    AddVolumeOutput("THERMAL_CONDUCTIVITY", "Thermal_Conductivity", "PRIMITIVE", "Thermal conductivity");

    AddVolumeOutput("SKIN_FRICTION-X", "Skin_Friction_Coefficient_x", "PRIMITIVE", "x-component of the skin friction vector");
    AddVolumeOutput("SKIN_FRICTION-Y", "Skin_Friction_Coefficient_y", "PRIMITIVE", "y-component of the skin friction vector");
    if (nDim == 3)
      AddVolumeOutput("SKIN_FRICTION-Z", "Skin_Friction_Coefficient_z", "PRIMITIVE", "z-component of the skin friction vector");

    AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heat-flux");
    AddVolumeOutput("Y_PLUS", "Y_Plus", "PRIMITIVE", "Non-dim. wall distance (Y-Plus)");

  }

  SetVolumeOutputFieldsScalarPrimitive(config);

  if (config->GetSAParsedOptions().bc) {
    AddVolumeOutput("INTERMITTENCY", "gamma_BC", "INTERMITTENCY", "Intermittency");
  }

  //Residuals
  AddVolumeOutput("RES_PRESSURE", "Residual_Pressure", "RESIDUAL", "Residual of the pressure");
  AddVolumeOutput("RES_VELOCITY-X", "Residual_Velocity_x", "RESIDUAL", "Residual of the x-velocity component");
  AddVolumeOutput("RES_VELOCITY-Y", "Residual_Velocity_y", "RESIDUAL", "Residual of the y-velocity component");
  if (nDim == 3)
    AddVolumeOutput("RES_VELOCITY-Z", "Residual_Velocity_z", "RESIDUAL", "Residual of the z-velocity component");
  if (config->GetEnergy_Equation())
    AddVolumeOutput("RES_TEMPERATURE", "Residual_Temperature", "RESIDUAL", "Residual of the temperature");

  SetVolumeOutputFieldsScalarResidual(config);

  if (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE && config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE) {
    AddVolumeOutput("LIMITER_PRESSURE", "Limiter_Pressure", "LIMITER", "Limiter value of the pressure");
    AddVolumeOutput("LIMITER_VELOCITY-X", "Limiter_Velocity_x", "LIMITER", "Limiter value of the x-velocity");
    AddVolumeOutput("LIMITER_VELOCITY-Y", "Limiter_Velocity_y", "LIMITER", "Limiter value of the y-velocity");
    if (nDim == 3)
      AddVolumeOutput("LIMITER_VELOCITY-Z", "Limiter_Velocity_z", "LIMITER", "Limiter value of the z-velocity");
    if (heat || weakly_coupled_heat || flamelet)
      AddVolumeOutput("LIMITER_TEMPERATURE", "Limiter_Temperature", "LIMITER", "Limiter value of the temperature");
  }

  SetVolumeOutputFieldsScalarLimiter(config);

  SetVolumeOutputFieldsScalarSource(config);

  SetVolumeOutputFieldsScalarLookup(config);

  SetVolumeOutputFieldsScalarMisc(config);

  // Streamwise Periodicity
  if(streamwisePeriodic) {
    AddVolumeOutput("RECOVERED_PRESSURE", "Recovered_Pressure", "SOLUTION", "Recovered physical pressure");
    if (heat && streamwisePeriodic_temperature)
      AddVolumeOutput("RECOVERED_TEMPERATURE", "Recovered_Temperature", "SOLUTION", "Recovered physical temperature");
  }

  AddCommonFVMOutputs(config);

  if (config->GetTime_Domain()) {
    SetTimeAveragedFields();
  }
}

void CFlowIncOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
  const CVariable* Node_Heat = nullptr;
  const CVariable* Node_Rad = nullptr;
  auto* Node_Geo = geometry->nodes;

  if (weakly_coupled_heat){
    Node_Heat = solver[HEAT_SOL]->GetNodes();
  }

  LoadCoordinates(Node_Geo->GetCoord(iPoint), iPoint);

  SetVolumeOutputValue("PRESSURE",   iPoint, Node_Flow->GetSolution(iPoint, 0));
  SetVolumeOutputValue("VELOCITY-X", iPoint, Node_Flow->GetSolution(iPoint, 1));
  SetVolumeOutputValue("VELOCITY-Y", iPoint, Node_Flow->GetSolution(iPoint, 2));
  if (nDim == 3)
    SetVolumeOutputValue("VELOCITY-Z", iPoint, Node_Flow->GetSolution(iPoint, 3));

  if (heat || flamelet) SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Flow->GetSolution(iPoint, nDim+1));
  if (weakly_coupled_heat) SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Heat->GetSolution(iPoint, 0));

  // Radiation solver
  if (config->AddRadiation()){
    Node_Rad = solver[RAD_SOL]->GetNodes();
    SetVolumeOutputValue("P1-RAD", iPoint, Node_Rad->GetSolution(iPoint,0));
  }

  if (gridMovement){
    SetVolumeOutputValue("GRID_VELOCITY-X", iPoint, Node_Geo->GetGridVel(iPoint)[0]);
    SetVolumeOutputValue("GRID_VELOCITY-Y", iPoint, Node_Geo->GetGridVel(iPoint)[1]);
    if (nDim == 3)
      SetVolumeOutputValue("GRID_VELOCITY-Z", iPoint, Node_Geo->GetGridVel(iPoint)[2]);
  }

  const su2double factor = solver[FLOW_SOL]->GetReferenceDynamicPressure();
  SetVolumeOutputValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure(iPoint) - solver[FLOW_SOL]->GetPressure_Inf())/factor);
  SetVolumeOutputValue("DENSITY", iPoint, Node_Flow->GetDensity(iPoint));

  if (config->GetKind_Solver() == MAIN_SOLVER::INC_RANS || config->GetKind_Solver() == MAIN_SOLVER::INC_NAVIER_STOKES){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity(iPoint));
    SetVolumeOutputValue("HEAT_CAPACITY", iPoint, Node_Flow->GetSpecificHeatCp(iPoint));
    SetVolumeOutputValue("THERMAL_CONDUCTIVITY", iPoint, Node_Flow->GetThermalConductivity(iPoint));
  }

  SetVolumeOutputValue("RES_PRESSURE", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 0));
  SetVolumeOutputValue("RES_VELOCITY-X", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 1));
  SetVolumeOutputValue("RES_VELOCITY-Y", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 2));
  if (nDim == 3)
    SetVolumeOutputValue("RES_VELOCITY-Z", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 3));
  if (config->GetEnergy_Equation())
    SetVolumeOutputValue("RES_TEMPERATURE", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nDim+1));

  if (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE && config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE) {
    SetVolumeOutputValue("LIMITER_PRESSURE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 0));
    SetVolumeOutputValue("LIMITER_VELOCITY-X", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 1));
    SetVolumeOutputValue("LIMITER_VELOCITY-Y", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 2));
    if (nDim == 3)
      SetVolumeOutputValue("LIMITER_VELOCITY-Z", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 3));
    if (heat || weakly_coupled_heat)
      SetVolumeOutputValue("LIMITER_TEMPERATURE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nDim+1));
  }

  // All turbulence and species outputs.
  LoadVolumeDataScalar(config, solver, geometry, iPoint);

  // Streamwise Periodicity
  if (streamwisePeriodic) {
    SetVolumeOutputValue("RECOVERED_PRESSURE", iPoint, Node_Flow->GetStreamwise_Periodic_RecoveredPressure(iPoint));
    if (heat && streamwisePeriodic_temperature)
      SetVolumeOutputValue("RECOVERED_TEMPERATURE", iPoint, Node_Flow->GetStreamwise_Periodic_RecoveredTemperature(iPoint));
  }

  LoadCommonFVMOutputs(config, geometry, iPoint);

  if (config->GetTime_Domain()) {
    LoadTimeAveragedData(iPoint, Node_Flow);
  }
}

bool CFlowIncOutput::SetInitResiduals(const CConfig *config){

  return (config->GetTime_Marching() != TIME_MARCHING::STEADY && (curInnerIter == 0))||
         (config->GetTime_Marching() == TIME_MARCHING::STEADY && (curInnerIter < 2));

}
