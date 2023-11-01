/*!
 * \file CFlowCompOutput.cpp
 * \brief Main subroutines for compressible flow output
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

#include "../../include/output/CFlowCompOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CFlowCompOutput::CFlowCompOutput(const CConfig *config, unsigned short nDim) : CFlowOutput(config, nDim, false) {

  turb_model = config->GetKind_Turb_Model();

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
    requestedScreenFields.emplace_back("RMS_MOMENTUM-X");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-Y");
    requestedScreenFields.emplace_back("RMS_ENERGY");
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
  ss << "Zone " << config->GetiZone() << " (Comp. Fluid)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetVolume_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_DENSITY");

  if (config->GetFixed_CL_Mode()) {
    if (std::find(convFields.begin(), convFields.end(), "LIFT") != convFields.end()) {
      if (rank == MASTER_NODE)
        cout<<"  Fixed CL: Adding LIFT as Convergence Field to ensure convergence to target CL"<<endl;
      convFields.emplace_back("LIFT");
      newFunc.resize(convFields.size());
      oldFunc.resize(convFields.size());
      cauchySerie.resize(convFields.size(), vector<su2double>(nCauchy_Elems, 0.0));
    }
  }
}

void CFlowCompOutput::SetHistoryOutputFields(CConfig *config){

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables.
  /// DESCRIPTION: Root-mean square residual of the density.
  AddHistoryOutput("RMS_DENSITY",    "rms[Rho]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum x-component.
  AddHistoryOutput("RMS_MOMENTUM-X", "rms[RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum y-component.
  AddHistoryOutput("RMS_MOMENTUM-Y", "rms[RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("RMS_MOMENTUM-Z", "rms[RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the energy.
  AddHistoryOutput("RMS_ENERGY",     "rms[RhoE]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the energy.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_ScalarRMS_RES(config);
  /// END_GROUP

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the density.
  AddHistoryOutput("MAX_DENSITY",    "max[Rho]",  ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum x-component.
  AddHistoryOutput("MAX_MOMENTUM-X", "max[RhoU]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum y-component.
  AddHistoryOutput("MAX_MOMENTUM-Y", "max[RhoV]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("MAX_MOMENTUM-Z", "max[RhoW]", ScreenOutputFormat::FIXED,"MAX_RES", "Maximum residual of the z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the energy.
  AddHistoryOutput("MAX_ENERGY",     "max[RhoE]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum residual of the energy.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_ScalarMAX_RES(config);
  /// END_GROUP

  /// BEGIN_GROUP: BGS_RES, DESCRIPTION: The block Gauss Seidel residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the density.
  AddHistoryOutput("BGS_DENSITY",    "bgs[Rho]",  ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum x-component.
  AddHistoryOutput("BGS_MOMENTUM-X", "bgs[RhoU]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum y-component.
  AddHistoryOutput("BGS_MOMENTUM-Y", "bgs[RhoV]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the momentum y-component.",  HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("BGS_MOMENTUM-Z", "bgs[RhoW]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the z-component.",  HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the energy.
  AddHistoryOutput("BGS_ENERGY",     "bgs[RhoE]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the energy.",  HistoryFieldType::RESIDUAL);

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
    AddHistoryOutputPerSurface("PITCH", "pitch", ScreenOutputFormat::FIXED, "AEROELASTIC", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// END_GROUP
  }

  /// DESCRIPTION: Linear solver iterations
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");
  AddHistoryOutputFieldsScalarLinsol(config);

  AddHistoryOutput("MIN_DELTA_TIME", "Min DT", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current minimum local time step");
  AddHistoryOutput("MAX_DELTA_TIME", "Max DT", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current maximum local time step");

  AddHistoryOutput("MIN_CFL", "Min CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current minimum of the local CFL numbers");
  AddHistoryOutput("MAX_CFL", "Max CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current maximum of the local CFL numbers");
  AddHistoryOutput("AVG_CFL", "Avg CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current average of the local CFL numbers");

  /// BEGIN_GROUP: FIXED_CL, DESCRIPTION: Relevant outputs for the Fixed CL mode
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
  /// END_GROUP

  if (config->GetDeform_Mesh()){
    AddHistoryOutput("DEFORM_MIN_VOLUME", "MinVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Minimum volume in the mesh");
    AddHistoryOutput("DEFORM_MAX_VOLUME", "MaxVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Maximum volume in the mesh");
    AddHistoryOutput("DEFORM_ITER", "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", ScreenOutputFormat::FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");
  }

  AddAnalyzeSurfaceOutput(config);

  AddAerodynamicCoefficients(config);

  if (config->GetViscous()) {
    AddHistoryOutput("BUFFET", "Buffet", ScreenOutputFormat::SCIENTIFIC, "AERO_COEFF", "Buffet sensor", HistoryFieldType::COEFFICIENT);
  }

  AddHeatCoefficients(config);

  AddRotatingFrameCoefficients();

  AddCpInverseDesignOutput();

  AddNearfieldInverseDesignOutput();

}

void CFlowCompOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddCoordinates();

  // Solution variables
  AddVolumeOutput("DENSITY",    "Density",    "SOLUTION", "Density");
  AddVolumeOutput("MOMENTUM-X", "Momentum_x", "SOLUTION", "x-component of the momentum vector");
  AddVolumeOutput("MOMENTUM-Y", "Momentum_y", "SOLUTION", "y-component of the momentum vector");
  
  if (nDim == 3)
    AddVolumeOutput("MOMENTUM-Z", "Momentum_z", "SOLUTION", "z-component of the momentum vector");
  AddVolumeOutput("ENERGY",     "Energy",     "SOLUTION", "Energy");

  SetVolumeOutputFieldsScalarSolution(config);

  // Grid velocity
  if (gridMovement){
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY", "x-component of the grid velocity vector");
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY", "y-component of the grid velocity vector");
    if (nDim == 3 )
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY", "z-component of the grid velocity vector");
  }

  // Primitive variables
  AddVolumeOutput("PRESSURE",    "Pressure",                "PRIMITIVE", "Pressure");
  AddVolumeOutput("TEMPERATURE", "Temperature",             "PRIMITIVE", "Temperature");
  AddVolumeOutput("MACH",        "Mach",                    "PRIMITIVE", "Mach number");
  AddVolumeOutput("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient");
  AddVolumeOutput("VELOCITY-X", "Velocity_x", "PRIMITIVE", "x-component of the velocity vector");
  AddVolumeOutput("VELOCITY-Y", "Velocity_y", "PRIMITIVE", "y-component of the velocity vector");
  
  if (nDim == 3)
    AddVolumeOutput("VELOCITY-Z", "Velocity_z", "PRIMITIVE", "z-component of the velocity vector");

  // Datadriven fluid model
  if(config->GetKind_FluidModel() == DATADRIVEN_FLUID){
    AddVolumeOutput("EXTRAPOLATION", "Extrapolation", "PRIMITIVE", "Density, energy outside data range");
    AddVolumeOutput("FLUIDMODEL_NEWTONITER", "nIter_Newton", "PRIMITIVE", "Number of iterations evaluated by the Newton solver");
    AddVolumeOutput("ENTROPY", "Entropy", "PRIMITIVE", "Fluid entropy value");
  }

  if (config->GetViscous()) {
    AddVolumeOutput("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity");

    AddVolumeOutput("SKIN_FRICTION-X", "Skin_Friction_Coefficient_x", "PRIMITIVE", "x-component of the skin friction vector");
    AddVolumeOutput("SKIN_FRICTION-Y", "Skin_Friction_Coefficient_y", "PRIMITIVE", "y-component of the skin friction vector");
    if (nDim == 3)
      AddVolumeOutput("SKIN_FRICTION-Z", "Skin_Friction_Coefficient_z", "PRIMITIVE", "z-component of the skin friction vector");

    AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heat-flux");
    AddVolumeOutput("Y_PLUS", "Y_Plus", "PRIMITIVE", "Non-dim. wall distance (Y-Plus)");
  }

  SetVolumeOutputFieldsScalarPrimitive(config);

  //Residuals
  AddVolumeOutput("RES_DENSITY", "Residual_Density", "RESIDUAL", "Residual of the density");
  AddVolumeOutput("RES_MOMENTUM-X", "Residual_Momentum_x", "RESIDUAL", "Residual of the x-momentum component");
  AddVolumeOutput("RES_MOMENTUM-Y", "Residual_Momentum_y", "RESIDUAL", "Residual of the y-momentum component");
  if (nDim == 3)
    AddVolumeOutput("RES_MOMENTUM-Z", "Residual_Momentum_z", "RESIDUAL", "Residual of the z-momentum component");
  AddVolumeOutput("RES_ENERGY", "Residual_Energy", "RESIDUAL", "Residual of the energy");

  SetVolumeOutputFieldsScalarResidual(config);

  if (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE && config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE) {
    AddVolumeOutput("LIMITER_VELOCITY-X", "Limiter_Velocity_x", "LIMITER", "Limiter value of the x-velocity");
    AddVolumeOutput("LIMITER_VELOCITY-Y", "Limiter_Velocity_y", "LIMITER", "Limiter value of the y-velocity");
    if (nDim == 3) {
      AddVolumeOutput("LIMITER_VELOCITY-Z", "Limiter_Velocity_z", "LIMITER", "Limiter value of the z-velocity");
    }
    AddVolumeOutput("LIMITER_PRESSURE", "Limiter_Pressure", "LIMITER", "Limiter value of the pressure");
    AddVolumeOutput("LIMITER_DENSITY", "Limiter_Density", "LIMITER", "Limiter value of the density");
    AddVolumeOutput("LIMITER_ENTHALPY", "Limiter_Enthalpy", "LIMITER", "Limiter value of the enthalpy");
  }

  SetVolumeOutputFieldsScalarLimiter(config);

  SetVolumeOutputFieldsScalarSource(config);

  SetVolumeOutputFieldsScalarLookup(config);

  SetVolumeOutputFieldsScalarMisc(config);

  // Roe Low Dissipation
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS) {
    AddVolumeOutput("ROE_DISSIPATION", "Roe_Dissipation", "ROE_DISSIPATION", "Value of the Roe dissipation");
  }

  AddCommonFVMOutputs(config);

  if (config->GetTime_Domain()) {
    SetTimeAveragedFields();
  }
}

void CFlowCompOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
  auto* Node_Geo  = geometry->nodes;

  LoadCoordinates(Node_Geo->GetCoord(iPoint), iPoint);

  SetVolumeOutputValue("DENSITY",    iPoint, Node_Flow->GetSolution(iPoint, 0));
  SetVolumeOutputValue("MOMENTUM-X", iPoint, Node_Flow->GetSolution(iPoint, 1));
  SetVolumeOutputValue("MOMENTUM-Y", iPoint, Node_Flow->GetSolution(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("MOMENTUM-Z", iPoint, Node_Flow->GetSolution(iPoint, 3));
    SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(iPoint, 4));
  } else {
    SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(iPoint, 3));
  }

  if (gridMovement){
    SetVolumeOutputValue("GRID_VELOCITY-X", iPoint, Node_Geo->GetGridVel(iPoint)[0]);
    SetVolumeOutputValue("GRID_VELOCITY-Y", iPoint, Node_Geo->GetGridVel(iPoint)[1]);
    if (nDim == 3)
      SetVolumeOutputValue("GRID_VELOCITY-Z", iPoint, Node_Geo->GetGridVel(iPoint)[2]);
  }

  SetVolumeOutputValue("PRESSURE", iPoint, Node_Flow->GetPressure(iPoint));
  SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Flow->GetTemperature(iPoint));
  SetVolumeOutputValue("MACH", iPoint, sqrt(Node_Flow->GetVelocity2(iPoint))/Node_Flow->GetSoundSpeed(iPoint));

  const su2double factor = solver[FLOW_SOL]->GetReferenceDynamicPressure();
  SetVolumeOutputValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure(iPoint) - solver[FLOW_SOL]->GetPressure_Inf())/factor);
  SetVolumeOutputValue("VELOCITY-X", iPoint, Node_Flow->GetVelocity(iPoint, 0));
  SetVolumeOutputValue("VELOCITY-Y", iPoint, Node_Flow->GetVelocity(iPoint, 1));
  if (nDim == 3){
    SetVolumeOutputValue("VELOCITY-Z", iPoint, Node_Flow->GetVelocity(iPoint, 2));
  }

  if(config->GetKind_FluidModel() == DATADRIVEN_FLUID){
    SetVolumeOutputValue("EXTRAPOLATION", iPoint, Node_Flow->GetDataExtrapolation(iPoint));
    SetVolumeOutputValue("FLUIDMODEL_NEWTONITER", iPoint, Node_Flow->GetNewtonSolverIterations(iPoint));
    SetVolumeOutputValue("ENTROPY", iPoint, Node_Flow->GetEntropy(iPoint));
  }

  if (config->GetKind_Solver() == MAIN_SOLVER::RANS || config->GetKind_Solver() == MAIN_SOLVER::NAVIER_STOKES){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity(iPoint));
  }

  SetVolumeOutputValue("RES_DENSITY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 0));
  SetVolumeOutputValue("RES_MOMENTUM-X", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 1));
  SetVolumeOutputValue("RES_MOMENTUM-Y", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("RES_MOMENTUM-Z", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 3));
    SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 4));
  } else {
    SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 3));
  }

  if (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE && config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE) {
    SetVolumeOutputValue("LIMITER_VELOCITY-X", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 1));
    SetVolumeOutputValue("LIMITER_VELOCITY-Y", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 2));
    if (nDim == 3){
      SetVolumeOutputValue("LIMITER_VELOCITY-Z", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 3));
    }
    SetVolumeOutputValue("LIMITER_PRESSURE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nDim+1));
    SetVolumeOutputValue("LIMITER_DENSITY", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nDim+2));
    SetVolumeOutputValue("LIMITER_ENTHALPY", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nDim+3));
  }

  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    SetVolumeOutputValue("ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation(iPoint));
  }

  LoadVolumeDataScalar(config, solver, geometry, iPoint);

  LoadCommonFVMOutputs(config, geometry, iPoint);

  if (config->GetTime_Domain()) {
    LoadTimeAveragedData(iPoint, Node_Flow);
  }
}

void CFlowCompOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver)  {

  CSolver* flow_solver = solver[FLOW_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];

  SetHistoryOutputValue("RMS_DENSITY", log10(flow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_MOMENTUM-X", log10(flow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_MOMENTUM-Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 2)
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(3)));
  else {
    SetHistoryOutputValue("RMS_MOMENTUM-Z", log10(flow_solver->GetRes_RMS(3)));
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(4)));
  }
  SetHistoryOutputValue("MAX_DENSITY", log10(flow_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_MOMENTUM-X", log10(flow_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_MOMENTUM-Y", log10(flow_solver->GetRes_Max(2)));
  if (nDim == 2)
    SetHistoryOutputValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(3)));
  else {
    SetHistoryOutputValue("MAX_MOMENTUM-Z", log10(flow_solver->GetRes_Max(3)));
    SetHistoryOutputValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(4)));
  }
  if (multiZone){
    SetHistoryOutputValue("BGS_DENSITY", log10(flow_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_MOMENTUM-X", log10(flow_solver->GetRes_BGS(1)));
    SetHistoryOutputValue("BGS_MOMENTUM-Y", log10(flow_solver->GetRes_BGS(2)));
    if (nDim == 2)
      SetHistoryOutputValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(3)));
    else {
      SetHistoryOutputValue("BGS_MOMENTUM-Z", log10(flow_solver->GetRes_BGS(3)));
      SetHistoryOutputValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(4)));
    }
  }

  SetHistoryOutputValue("MIN_DELTA_TIME", flow_solver->GetMin_Delta_Time());
  SetHistoryOutputValue("MAX_DELTA_TIME", flow_solver->GetMax_Delta_Time());

  SetHistoryOutputValue("MIN_CFL", flow_solver->GetMin_CFL_Local());
  SetHistoryOutputValue("MAX_CFL", flow_solver->GetMax_CFL_Local());
  SetHistoryOutputValue("AVG_CFL", flow_solver->GetAvg_CFL_Local());

  SetHistoryOutputValue("LINSOL_ITER", flow_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(flow_solver->GetResLinSolver()));

  if (config->GetDeform_Mesh()){
    SetHistoryOutputValue("DEFORM_MIN_VOLUME", mesh_solver->GetMinimum_Volume());
    SetHistoryOutputValue("DEFORM_MAX_VOLUME", mesh_solver->GetMaximum_Volume());
    SetHistoryOutputValue("DEFORM_ITER", mesh_solver->GetIterLinSolver());
    SetHistoryOutputValue("DEFORM_RESIDUAL", log10(mesh_solver->GetResLinSolver()));
  }

  if(config->GetFixed_CL_Mode()){
    SetHistoryOutputValue("DELTA_CL", fabs(flow_solver->GetTotal_CL() - config->GetTarget_CL()));
    SetHistoryOutputValue("PREV_AOA", flow_solver->GetPrevious_AoA());
    SetHistoryOutputValue("CHANGE_IN_AOA", config->GetAoA()-flow_solver->GetPrevious_AoA());
    SetHistoryOutputValue("CL_DRIVER_COMMAND", flow_solver->GetAoA_inc());
  }

  LoadHistoryDataScalar(config, solver);

  /*--- Set the analyse surface history values --- */

  SetAnalyzeSurface(solver, geometry, config, false);

  /*--- Set aerodynamic coefficients --- */

  SetAerodynamicCoefficients(config, flow_solver);

  if (config->GetViscous()) {
    SetHistoryOutputValue("BUFFET", flow_solver->GetTotal_Buffet_Metric());
  }

  SetHeatCoefficients(config, flow_solver);

  /*--- Set rotating frame coefficients --- */

  SetRotatingFrameCoefficients(flow_solver);

  /*--- Set Cp diff fields ---*/

  SetCpInverseDesign(flow_solver, geometry, config);

  /*--- Set nearfield diff fields ---*/

  if (config->GetEquivArea()) SetNearfieldInverseDesign(flow_solver, geometry, config);

  /*--- Keep this as last, since it uses the history values that were set. ---*/

  SetCustomOutputs(solver, geometry, config);

  SetCustomAndComboObjectives(FLOW_SOL, config, solver);

}

bool CFlowCompOutput::SetInitResiduals(const CConfig *config){

  return (config->GetTime_Marching() != TIME_MARCHING::STEADY && (curInnerIter == 0))||
         (config->GetTime_Marching() == TIME_MARCHING::STEADY && (curInnerIter < 2));

}

void CFlowCompOutput::SetAdditionalScreenOutput(const CConfig *config){

  if (config->GetFixed_CL_Mode()){
    SetFixedCLScreenOutput(config);
  }
}

bool CFlowCompOutput::WriteHistoryFileOutput(const CConfig *config) {
  return !config->GetFinite_Difference_Mode() && COutput::WriteHistoryFileOutput(config);
}
