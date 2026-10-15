/*!
 * \file CNEMOCompOutput.cpp
 * \brief Main subroutines for compressible flow output
 * \author W. Maier, R. Sanchez
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/output/CNEMOCompOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CNEMOCompOutput::CNEMOCompOutput(const CConfig *config, unsigned short nDim) : CFlowOutput(config, nDim, false) {

  turb_model = config->GetKind_Turb_Model();
  nSpecies      = config->GetnSpecies();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    RequestCommonHistory(config->GetTime_Domain());
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    if (config->GetTime_Domain()) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
      requestedScreenFields.emplace_back("RMS_DENSITY_" + std::to_string(iSpecies));
    requestedScreenFields.emplace_back("RMS_MOMENTUM-X");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-Y");
    requestedScreenFields.emplace_back("RMS_ENERGY");
    requestedScreenFields.emplace_back("RMS_ENERGY_VE");
    nRequestedScreenFields = requestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("PRIMITIVE");
    requestedVolumeFields.emplace_back("AUXILIARY");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  if (gridMovement) {
    auto notFound = requestedVolumeFields.end();
    if (find(requestedVolumeFields.begin(), notFound, string("GRID_VELOCITY")) == notFound) {
      requestedVolumeFields.emplace_back("GRID_VELOCITY");
      nRequestedVolumeFields++;
    }
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

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_DENSITY_0");

  if (config->GetFixed_CL_Mode()) {
    bool found = false;
    for (unsigned short iField = 0; iField < convFields.size(); iField++)
      if (convFields[iField] == "LIFT") found = true;
    if (!found) {
      if (rank == MASTER_NODE)
        cout<<"  Fixed CL: Adding LIFT as Convergence Field to ensure convergence to target CL"<<endl;
      convFields.emplace_back("LIFT");
      newFunc.resize(convFields.size());
      oldFunc.resize(convFields.size());
      cauchySerie.resize(convFields.size(), vector<su2double>(nCauchy_Elems, 0.0));
    }
  }
}

void CNEMOCompOutput::SetHistoryOutputFields(CConfig *config){

  const auto addResidualFields = [&](const string& prefix, const string& namePrefix, const string& description,
                                     const string& energyVEDescription) {
    const auto group = prefix + "_RES";
    for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
      AddHistoryOutput(prefix + "_DENSITY_" + std::to_string(iSpecies),
                       namePrefix + "[Rho_" + std::to_string(iSpecies) + "]", ScreenOutputFormat::FIXED, group,
                       description + " residual of the species density " + std::to_string(iSpecies) + ".",
                       HistoryFieldType::RESIDUAL);

    AddHistoryOutput(prefix + "_MOMENTUM-X", namePrefix + "[RhoU]", ScreenOutputFormat::FIXED, group,
                     description + " residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
    AddHistoryOutput(prefix + "_MOMENTUM-Y", namePrefix + "[RhoV]", ScreenOutputFormat::FIXED, group,
                     description + " residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
    if (nDim == 3)
      AddHistoryOutput(prefix + "_MOMENTUM-Z", namePrefix + "[RhoW]", ScreenOutputFormat::FIXED, group,
                       description + " residual of the momentum z-component.", HistoryFieldType::RESIDUAL);
    AddHistoryOutput(prefix + "_ENERGY", namePrefix + "[RhoE]", ScreenOutputFormat::FIXED, group,
                     description + " residual of the energy.", HistoryFieldType::RESIDUAL);
    AddHistoryOutput(prefix + "_ENERGY_VE", namePrefix + "[RhoEve]", ScreenOutputFormat::FIXED, group,
                     description + " residual of the " + energyVEDescription + ".", HistoryFieldType::RESIDUAL);
  };

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables.
  addResidualFields("RMS", "rms", "Root-mean square", "energy");
  AddHistoryOutputFields_ScalarRMS_RES(config);
  /// END_GROUP

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables.
  addResidualFields("MAX", "max", "Maximum", "vibrational-electronic energy");
  AddHistoryOutputFields_ScalarMAX_RES(config);
  /// END_GROUP

  /// BEGIN_GROUP: BGS_RES, DESCRIPTION: The block Gauss Seidel residuals of the SOLUTION variables.
  addResidualFields("BGS", "bgs", "BGS", "vibrational-electronic energy");
  AddHistoryOutputFields_ScalarBGS_RES(config);
  /// END_GROUP

  vector<string> Marker_Monitoring;
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
  }

  if (config->GetAeroelastic_Simulation()) {
    /// BEGIN_GROUP: AEROELASTIC, DESCRIPTION: Aeroelastic plunge, pitch
    /// DESCRIPTION: Aeroelastic plunge
    AddHistoryOutputPerSurface("PLUNGE", "plunge", ScreenOutputFormat::FIXED, "AEROELASTIC", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Aeroelastic pitch
    AddHistoryOutputPerSurface("PITCH",  "pitch",  ScreenOutputFormat::FIXED, "AEROELASTIC", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// END_GROUP
  }

  /// DESCRIPTION: Linear solver iterations
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");
  AddHistoryOutputFieldsScalarLinsol(config);

  AddHistoryOutput("MIN_CFL", "Min CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current minimum of the local CFL numbers");
  AddHistoryOutput("MAX_CFL", "Max CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current maximum of the local CFL numbers");
  AddHistoryOutput("AVG_CFL", "Avg CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current average of the local CFL numbers");

  ///   /// BEGIN_GROUP: FIXED_CL, DESCRIPTION: Relevant outputs for the Fixed CL mode
  if (config->GetFixed_CL_Mode()){
    /// DESCRIPTION: Difference between current and target CL
    AddHistoryOutput("DELTA_CL", "Delta_CL", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL", "Difference between Target CL and current CL", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Angle of attack before the most recent update
    AddHistoryOutput("PREV_AOA", "Previous_AOA", ScreenOutputFormat::FIXED, "FIXED_CL", "Angle of Attack at the previous iteration of the Fixed CL driver");
    /// DESCRIPTION: Last change in angle of attack by the Fixed CL driver
    AddHistoryOutput("CHANGE_IN_AOA", "Change_in_AOA", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL", "Last change in Angle of Attack by Fixed CL Driver", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: AOA control command by the CL Driver
    AddHistoryOutput("CL_DRIVER_COMMAND", "CL_Driver_Command", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL", "CL Driver's control command", HistoryFieldType::RESIDUAL);
  }

  if (config->GetDeform_Mesh()){
    AddHistoryOutput("DEFORM_MIN_VOLUME", "MinVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Minimum volume in the mesh");
    AddHistoryOutput("DEFORM_MAX_VOLUME", "MaxVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Maximum volume in the mesh");
    AddHistoryOutput("DEFORM_ITER", "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", ScreenOutputFormat::FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");
  }

  AddAnalyzeSurfaceOutput(config);

  AddAerodynamicCoefficients(config);

  AddHeatCoefficients(config);

  AddRotatingFrameCoefficients();

  AddCpInverseDesignOutput();

}

void CNEMOCompOutput::SetVolumeOutputFields(CConfig *config){

  unsigned short nSpecies = config->GetnSpecies();

  // Grid coordinates
  AddCoordinates();

  // Solution variables
  for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
    AddVolumeOutput("DENSITY_" + std::to_string(iSpecies),  "Density_" + std::to_string(iSpecies),  "SOLUTION", "Density_"  + std::to_string(iSpecies));

  AddVolumeOutput("MOMENTUM-X", "Momentum_x", "SOLUTION", "x-component of the momentum vector");
  AddVolumeOutput("MOMENTUM-Y", "Momentum_y", "SOLUTION", "y-component of the momentum vector");
  if (nDim == 3)
    AddVolumeOutput("MOMENTUM-Z", "Momentum_z", "SOLUTION", "z-component of the momentum vector");
  AddVolumeOutput("ENERGY",       "Energy",     "SOLUTION", "Energy");
  AddVolumeOutput("ENERGY_VE",    "Energy_ve",  "SOLUTION", "Energy_ve");

  SetVolumeOutputFieldsScalarSolution(config);

  //Auxiliary variables for post-processment
  for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
    AddVolumeOutput("MASSFRAC_" + std::to_string(iSpecies),  "MassFrac_" + std::to_string(iSpecies),  "AUXILIARY", "MassFrac_" + std::to_string(iSpecies));

  // Grid velocity
  if (gridMovement){
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY", "x-component of the grid velocity vector");
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY", "y-component of the grid velocity vector");
    if (nDim == 3 )
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY", "z-component of the grid velocity vector");
  }

  // Primitive variables
  AddVolumeOutput("PRESSURE",       "Pressure",       "PRIMITIVE", "Pressure");
  AddVolumeOutput("TEMPERATURE_TR", "Temperature_tr", "PRIMITIVE", "Temperature_tr");
  AddVolumeOutput("TEMPERATURE_VE", "Temperature_ve", "PRIMITIVE", "Temperature_ve");
  AddVolumeOutput("VELOCITY-X", "Velocity_x", "PRIMITIVE", "x-component of the velocity vector");
  AddVolumeOutput("VELOCITY-Y", "Velocity_y", "PRIMITIVE", "y-component of the velocity vector");
  if (nDim == 3)
    AddVolumeOutput("VELOCITY-Z", "Velocity_z", "PRIMITIVE", "z-component of the velocity vector");

  AddVolumeOutput("MACH",        "Mach",                    "PRIMITIVE", "Mach number");
  AddVolumeOutput("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient");

  if (config->GetViscous()) {
    AddVolumeOutput("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity");

    AddVolumeOutput("SKIN_FRICTION-X", "Skin_Friction_Coefficient_x", "PRIMITIVE", "x-component of the skin friction vector");
    AddVolumeOutput("SKIN_FRICTION-Y", "Skin_Friction_Coefficient_y", "PRIMITIVE", "y-component of the skin friction vector");
    if (nDim == 3)
     AddVolumeOutput("SKIN_FRICTION-Z", "Skin_Friction_Coefficient_z", "PRIMITIVE", "z-component of the skin friction vector");

    AddVolumeOutput("THERMAL_CONDUCTIVITY_TR", "Thermal_Conductivity_tr", "PRIMITIVE", "Translational-rotational thermal conductivity");
    AddVolumeOutput("THERMAL_CONDUCTIVITY_VE", "Thermal_Conductivity_ve", "PRIMITIVE", "Vibrational-electronic thermal conductivity");
    AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heat-flux");
    AddVolumeOutput("Y_PLUS", "Y_Plus", "PRIMITIVE", "Non-dim. wall distance (Y-Plus)");

  }

  SetVolumeOutputFieldsScalarPrimitive(config);

  //Residuals
  for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
    AddVolumeOutput("RES_DENSITY_" + std::to_string(iSpecies), "Residual_Density_" + std::to_string(iSpecies), "RESIDUAL", "Residual of species density " + std::to_string(iSpecies));
  AddVolumeOutput("RES_MOMENTUM-X", "Residual_Momentum_x", "RESIDUAL", "Residual of the x-momentum component");
  AddVolumeOutput("RES_MOMENTUM-Y", "Residual_Momentum_y", "RESIDUAL", "Residual of the y-momentum component");
  if (nDim == 3)
    AddVolumeOutput("RES_MOMENTUM-Z", "Residual_Momentum_z", "RESIDUAL", "Residual of the z-momentum component");
  AddVolumeOutput("RES_ENERGY",    "Residual_Energy",    "RESIDUAL", "Residual of the energy");
  AddVolumeOutput("RES_ENERGY_VE", "Residual_Energy_ve", "RESIDUAL", "Residual of the energy_ve");

  SetVolumeOutputFieldsScalarResidual(config);

  if (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE && config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE) {
    // Limiter values (of the primitive variables, which are reconstructed)
    for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
      AddVolumeOutput("LIMITER_DENSITY_" + std::to_string(iSpecies), "Limiter_Density_" + std::to_string(iSpecies), "LIMITER", "Limiter value of the density of species " + std::to_string(iSpecies));
    AddVolumeOutput("LIMITER_TEMPERATURE_TR", "Limiter_Temperature_tr", "LIMITER", "Limiter value of the temperature_tr");
    AddVolumeOutput("LIMITER_TEMPERATURE_VE", "Limiter_Temperature_ve", "LIMITER", "Limiter value of the temperature_ve");
    AddVolumeOutput("LIMITER_VELOCITY-X", "Limiter_Velocity_x", "LIMITER", "Limiter value of the x-velocity");
    AddVolumeOutput("LIMITER_VELOCITY-Y", "Limiter_Velocity_y", "LIMITER", "Limiter value of the y-velocity");
    if (nDim == 3)
      AddVolumeOutput("LIMITER_VELOCITY-Z", "Limiter_Velocity_z", "LIMITER", "Limiter value of the z-velocity");
    AddVolumeOutput("LIMITER_PRESSURE", "Limiter_Pressure", "LIMITER", "Limiter value of the pressure");
  }

  SetVolumeOutputFieldsScalarLimiter(config);

  SetVolumeOutputFieldsScalarSource(config);

  SetVolumeOutputFieldsScalarLookup(config);

  SetVolumeOutputFieldsScalarMisc(config);

  AddCommonFVMOutputs(config);

  if (config->GetTime_Domain()) {
    SetTimeAveragedFields(config);
  }
}

void CNEMOCompOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
  const CVariable* Node_Turb = nullptr;
  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE) Node_Turb = solver[TURB_SOL]->GetNodes();
  auto* Node_Geo = geometry->nodes;
  const auto nSpecies = config->GetnSpecies();

  LoadCoordinates(Node_Geo->GetCoord(iPoint), iPoint);

  for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
    SetVolumeOutputValue("DENSITY_" + std::to_string(iSpecies),   iPoint, Node_Flow->GetSolution(iPoint, iSpecies));

  SetVolumeOutputValue("MOMENTUM-X", iPoint, Node_Flow->GetSolution(iPoint, nSpecies));
  SetVolumeOutputValue("MOMENTUM-Y", iPoint, Node_Flow->GetSolution(iPoint, nSpecies+1));
  if (nDim == 3){
    SetVolumeOutputValue("MOMENTUM-Z", iPoint, Node_Flow->GetSolution(iPoint, nSpecies+2));
    SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(iPoint, nSpecies+3));
    SetVolumeOutputValue("ENERGY_VE",  iPoint, Node_Flow->GetSolution(iPoint, nSpecies+4));
  } else {
    SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(iPoint, nSpecies+2));
    SetVolumeOutputValue("ENERGY_VE",  iPoint, Node_Flow->GetSolution(iPoint, nSpecies+3));
  }

  for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
    SetVolumeOutputValue("MASSFRAC_" + std::to_string(iSpecies),   iPoint, Node_Flow->GetSolution(iPoint, iSpecies)/Node_Flow->GetDensity(iPoint));

  if (gridMovement){
    SetVolumeOutputValue("GRID_VELOCITY-X", iPoint, Node_Geo->GetGridVel(iPoint)[0]);
    SetVolumeOutputValue("GRID_VELOCITY-Y", iPoint, Node_Geo->GetGridVel(iPoint)[1]);
    if (nDim == 3)
      SetVolumeOutputValue("GRID_VELOCITY-Z", iPoint, Node_Geo->GetGridVel(iPoint)[2]);
  }

  SetVolumeOutputValue("PRESSURE", iPoint, Node_Flow->GetPressure(iPoint));
  SetVolumeOutputValue("TEMPERATURE_TR", iPoint, Node_Flow->GetTemperature(iPoint));
  SetVolumeOutputValue("TEMPERATURE_VE", iPoint, Node_Flow->GetTemperature_ve(iPoint));
  SetVolumeOutputValue("VELOCITY-X", iPoint, Node_Flow->GetVelocity(iPoint, 0));
  SetVolumeOutputValue("VELOCITY-Y", iPoint, Node_Flow->GetVelocity(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("VELOCITY-Z", iPoint, Node_Flow->GetVelocity(iPoint, 2));

  SetVolumeOutputValue("MACH", iPoint, sqrt(Node_Flow->GetVelocity2(iPoint))/Node_Flow->GetSoundSpeed(iPoint));

  const su2double factor = solver[FLOW_SOL]->GetReferenceDynamicPressure();
  SetVolumeOutputValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure(iPoint) - solver[FLOW_SOL]->GetPressure_Inf())/factor);

  if (config->GetViscous()){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity(iPoint));
    SetVolumeOutputValue("THERMAL_CONDUCTIVITY_TR", iPoint, Node_Flow->GetThermalConductivity(iPoint));
    SetVolumeOutputValue("THERMAL_CONDUCTIVITY_VE", iPoint, Node_Flow->GetThermalConductivity_ve(iPoint));
  }

  for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
    SetVolumeOutputValue("RES_DENSITY_" + std::to_string(iSpecies), iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, iSpecies));

  SetVolumeOutputValue("RES_MOMENTUM-X", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies));
  SetVolumeOutputValue("RES_MOMENTUM-Y", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies+1));
  if (nDim == 3){
    SetVolumeOutputValue("RES_MOMENTUM-Z", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies+2));
    SetVolumeOutputValue("RES_ENERGY",     iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies+3));
    SetVolumeOutputValue("RES_ENERGY_VE",  iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies+4));
  } else {
    SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies+2));
    SetVolumeOutputValue("RES_ENERGY_VE", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies+3));
  }

  if (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE && config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE) {
    /*--- Primitive variables: species densities, T, Tve, velocity, pressure. ---*/
    for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
      SetVolumeOutputValue("LIMITER_DENSITY_" + std::to_string(iSpecies), iPoint, Node_Flow->GetLimiter_Primitive(iPoint, iSpecies));
    SetVolumeOutputValue("LIMITER_TEMPERATURE_TR", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies));
    SetVolumeOutputValue("LIMITER_TEMPERATURE_VE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies+1));
    SetVolumeOutputValue("LIMITER_VELOCITY-X", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies+2));
    SetVolumeOutputValue("LIMITER_VELOCITY-Y", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies+3));
    if (nDim == 3)
      SetVolumeOutputValue("LIMITER_VELOCITY-Z", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies+4));
    SetVolumeOutputValue("LIMITER_PRESSURE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies+nDim+2));
  }

  LoadVolumeDataScalar(config, solver, geometry, iPoint);

  LoadCommonFVMOutputs(config, geometry, iPoint);

  if (config->GetTime_Domain()) {
    LoadTimeAveragedData(iPoint, Node_Flow, Node_Turb, config);
  }
}

void CNEMOCompOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver)  {

  CSolver* NEMO_solver = solver[FLOW_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];
  const unsigned short nSpecies = config->GetnSpecies();

  const auto loadResiduals = [&](const string& prefix, auto getResidual) {
    for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
      SetHistoryOutputValue(prefix + "_DENSITY_" + std::to_string(iSpecies),
                            log10((NEMO_solver->*getResidual)(iSpecies)));

    SetHistoryOutputValue(prefix + "_MOMENTUM-X", log10((NEMO_solver->*getResidual)(nSpecies)));
    SetHistoryOutputValue(prefix + "_MOMENTUM-Y", log10((NEMO_solver->*getResidual)(nSpecies + 1)));
    if (nDim == 3)
      SetHistoryOutputValue(prefix + "_MOMENTUM-Z", log10((NEMO_solver->*getResidual)(nSpecies + 2)));
    SetHistoryOutputValue(prefix + "_ENERGY", log10((NEMO_solver->*getResidual)(nSpecies + nDim)));
    SetHistoryOutputValue(prefix + "_ENERGY_VE", log10((NEMO_solver->*getResidual)(nSpecies + nDim + 1)));
  };

  loadResiduals("RMS", &CSolver::GetRes_RMS);
  loadResiduals("MAX", &CSolver::GetRes_Max);
  if (multiZone) loadResiduals("BGS", &CSolver::GetRes_BGS);

  SetHistoryOutputValue("MIN_CFL", NEMO_solver->GetMin_CFL_Local());
  SetHistoryOutputValue("MAX_CFL", NEMO_solver->GetMax_CFL_Local());
  SetHistoryOutputValue("AVG_CFL", NEMO_solver->GetAvg_CFL_Local());

  SetHistoryOutputValue("LINSOL_ITER", NEMO_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(NEMO_solver->GetResLinSolver()));

  if (config->GetDeform_Mesh()){
    SetHistoryOutputValue("DEFORM_MIN_VOLUME", mesh_solver->GetMinimum_Volume());
    SetHistoryOutputValue("DEFORM_MAX_VOLUME", mesh_solver->GetMaximum_Volume());
    SetHistoryOutputValue("DEFORM_ITER", mesh_solver->GetIterLinSolver());
    SetHistoryOutputValue("DEFORM_RESIDUAL", log10(mesh_solver->GetResLinSolver()));
  }

  if(config->GetFixed_CL_Mode()){
    SetHistoryOutputValue("DELTA_CL", fabs(NEMO_solver->GetTotal_CL() - config->GetTarget_CL()));
    SetHistoryOutputValue("PREV_AOA", NEMO_solver->GetPrevious_AoA());
    SetHistoryOutputValue("CHANGE_IN_AOA", config->GetAoA()-NEMO_solver->GetPrevious_AoA());
    SetHistoryOutputValue("CL_DRIVER_COMMAND", NEMO_solver->GetAoA_inc());
  }

  LoadHistoryDataScalar(config, solver);

  /*--- Set the analyse surface history values --- */

  SetAnalyzeSurface(solver, geometry, config, false);

  /*--- Set aeroydnamic coefficients --- */

  SetAerodynamicCoefficients(config, NEMO_solver);

  SetHeatCoefficients(config, NEMO_solver);

  SetRotatingFrameCoefficients(NEMO_solver);

  /*--- Set Cp diff fields ---*/

  SetCpInverseDesign(NEMO_solver, geometry, config);

  /*--- Keep this as last, since it uses the history values that were set. ---*/

  SetCustomOutputs(solver, geometry, config);

  SetCustomAndComboObjectives(FLOW_SOL, config, solver);

}

bool CNEMOCompOutput::SetInitResiduals(const CConfig *config){

  return (config->GetTime_Marching() != TIME_MARCHING::STEADY && (curInnerIter == 0))||
         (config->GetTime_Marching() == TIME_MARCHING::STEADY && (curInnerIter < 2));

}

void CNEMOCompOutput::SetAdditionalScreenOutput(const CConfig *config){

  if (config->GetFixed_CL_Mode()){
    SetFixedCLScreenOutput(config);
  }
}

bool CNEMOCompOutput::WriteHistoryFileOutput(const CConfig *config) {
  return !config->GetFinite_Difference_Mode() && COutput::WriteHistoryFileOutput(config);
}
