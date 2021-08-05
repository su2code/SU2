/*!
 * \file COneShotOutput.cpp
 * \brief Main subroutines for One Shot output
 * \author T. Dick
 * \version 7.1.1 "Blackbird"
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


#include "../../include/output/COneShotOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

COneShotOutput::COneShotOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  turb_model = config->GetKind_Turb_Model();

  cont_adj = config->GetContinuous_Adjoint();

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
    requestedScreenFields.emplace_back("RMS_DENSITY");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-X");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-Y");
    requestedScreenFields.emplace_back("RMS_ENERGY");
    requestedScreenFields.emplace_back("RMS_ADJ_DENSITY");
    requestedScreenFields.emplace_back("RMS_ADJ_MOMENTUM-X");
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
  ss << "Zone " << config->GetiZone() << " (Adj. Comp. Fluid)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFlowFilename = config->GetVolume_FileName();
  volumeFilename = config->GetAdj_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfAdjCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFlowFilename = config->GetRestart_FileName();
  restartFilename = config->GetRestart_AdjFileName();

  /*--- Add the obj. function extension --- */

  restartFilename = config->GetObjFunc_Extension(restartFilename);

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_ADJ_DENSITY");

}

COneShotOutput::~COneShotOutput(void) {}

void COneShotOutput::SetHistoryOutputFields(CConfig *config){

  // For OneShot we are interested in the flow and adjoint solver outputs!
  // First group, residuals (BGS residuals are not included, only relevant for multizone.)
  {
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

  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).
    AddHistoryOutput("RMS_NU_TILDE", "rms[nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
    break;
  case SST: case SST_SUST:
    /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).
    AddHistoryOutput("RMS_TKE", "rms[k]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).
    AddHistoryOutput("RMS_DISSIPATION", "rms[w]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of dissipation (SST model).", HistoryFieldType::RESIDUAL);
    break;
  default: break;
  }
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

  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Maximum residual of nu tilde (SA model).
    AddHistoryOutput("MAX_NU_TILDE",       "max[nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
    break;
  case SST: case SST_SUST:
    /// DESCRIPTION: Maximum residual of kinetic energy (SST model).
    AddHistoryOutput("MAX_TKE", "max[k]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Maximum residual of the dissipation (SST model).
    AddHistoryOutput("MAX_DISSIPATION",    "max[w]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of dissipation (SST model).", HistoryFieldType::RESIDUAL);
    break;
  default: break;
  }
  /// END_GROUP

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables.
  /// DESCRIPTION: Root-mean square residual of the adjoint density.
  AddHistoryOutput("RMS_ADJ_DENSITY",    "rms[A_Rho]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum x-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-X", "rms[A_RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum y-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-Y", "rms[A_RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum z-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-Z", "rms[A_RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint energy.
  AddHistoryOutput("RMS_ADJ_ENERGY",     "rms[A_E]",    ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint energy.", HistoryFieldType::RESIDUAL);
  if ((!config->GetFrozen_Visc_Disc() && !cont_adj) || (!config->GetFrozen_Visc_Cont() && cont_adj)){
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      /// DESCRIPTION: Root-mean square residual of the adjoint nu tilde.
      AddHistoryOutput("RMS_ADJ_NU_TILDE", "rms[A_nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);
      break;
    case SST:
      /// DESCRIPTION: Root-mean square residual of the adjoint kinetic energy.
      AddHistoryOutput("RMS_ADJ_TKE", "rms[A_k]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint kinetic energy.", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: Root-mean square residual of the adjoint dissipation.
      AddHistoryOutput("RMS_ADJ_DISSIPATION",    "rms[A_w]", ScreenOutputFormat::FIXED, "RMS_RES", " Root-mean square residual of the adjoint dissipation.", HistoryFieldType::RESIDUAL);
      break;
    default: break;
    }
  }
  /// END_GROUP

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the adjoint density.
  AddHistoryOutput("MAX_ADJ_DENSITY",    "max[A_Rho]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint momentum x-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-X", "max[A_RhoU]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum x-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint momentum y-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-Y", "max[A_RhoV]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum y-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint momentum z-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-Z", "max[A_RhoW]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum z-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint energy.
  AddHistoryOutput("MAX_ADJ_ENERGY",     "max[A_E]",    ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint energy.", HistoryFieldType::RESIDUAL);
  if (!config->GetFrozen_Visc_Disc()){
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      /// DESCRIPTION: Maximum residual of the adjoint nu tilde.
      AddHistoryOutput("MAX_ADJ_NU_TILDE", "max[A_nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);
      break;
    case SST:
      /// DESCRIPTION: Maximum residual of the adjoint kinetic energy.
      AddHistoryOutput("MAX_ADJ_TKE", "max[A_k]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint kinetic energy.", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: Maximum residual of the adjoint dissipation.
      AddHistoryOutput("MAX_ADJ_DISSIPATION",    "max[A_w]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint dissipation.", HistoryFieldType::RESIDUAL);
      break;
    default: break;
    }
  }
  /// END_GROUP

  }

  // Second group, aerodynamic coefficients
  {
    /// BEGIN_GROUP: AERO_COEFF, DESCRIPTION: Sum of the aerodynamic coefficients and forces on all surfaces (markers) set with MARKER_MONITORING.
    /// DESCRIPTION: Drag coefficient
    AddHistoryOutput("DRAG",       "CD",   ScreenOutputFormat::FIXED, "AERO_COEFF", "Total drag coefficient on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Lift coefficient
    AddHistoryOutput("LIFT",       "CL",   ScreenOutputFormat::FIXED, "AERO_COEFF", "Total lift coefficient on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Sideforce coefficient
    AddHistoryOutput("SIDEFORCE",  "CSF",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total sideforce coefficient on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Moment around the x-axis
    AddHistoryOutput("MOMENT_X",   "CMx",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total momentum x-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Moment around the y-axis
    AddHistoryOutput("MOMENT_Y",   "CMy",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total momentum y-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Moment around the z-axis
    AddHistoryOutput("MOMENT_Z",   "CMz",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total momentum z-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Force in x direction
    AddHistoryOutput("FORCE_X",    "CFx",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total force x-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Force in y direction
    AddHistoryOutput("FORCE_Y",    "CFy",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total force y-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Force in z direction
    AddHistoryOutput("FORCE_Z",    "CFz",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total force z-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Lift-to-drag ratio
    AddHistoryOutput("EFFICIENCY", "CEff", ScreenOutputFormat::FIXED, "AERO_COEFF", "Total lift-to-drag ratio on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
    /// END_GROUP

    /// BEGIN_GROUP: AERO_COEFF_SURF, DESCRIPTION: Aerodynamic coefficients and forces per surface.
    vector<string> Marker_Monitoring;
    for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
      Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
    }
    /// DESCRIPTION: Drag coefficient
    AddHistoryOutputPerSurface("DRAG_ON_SURFACE",       "CD",   ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Lift coefficient
    AddHistoryOutputPerSurface("LIFT_ON_SURFACE",       "CL",   ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Sideforce coefficient
    AddHistoryOutputPerSurface("SIDEFORCE_ON_SURFACE",  "CSF",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Moment around the x-axis
    AddHistoryOutputPerSurface("MOMENT-X_ON_SURFACE",   "CMx",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Moment around the y-axis
    AddHistoryOutputPerSurface("MOMENT-Y_ON_SURFACE",   "CMy",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Moment around the z-axis
    AddHistoryOutputPerSurface("MOMENT-Z_ON_SURFACE",   "CMz",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Force in x direction
    AddHistoryOutputPerSurface("FORCE-X_ON_SURFACE",    "CFx",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Force in y direction
    AddHistoryOutputPerSurface("FORCE-Y_ON_SURFACE",    "CFy",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Force in z direction
    AddHistoryOutputPerSurface("FORCE-Z_ON_SURFACE",    "CFz",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Lift-to-drag ratio
    AddHistoryOutputPerSurface("EFFICIENCY_ON_SURFACE", "CEff", ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// END_GROUP

    /// DESCRIPTION: Angle of attack
    AddHistoryOutput("AOA", "AoA", ScreenOutputFormat::FIXED, "AOA", "Angle of attack");
  }

  // Third group, sensitivities
  {
  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Sensitivities of different geometrical or boundary values.
  /// DESCRIPTION: Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.
  AddHistoryOutput("SENS_GEO",   "Sens_Geo",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sensitivity of the objective function with respect to the angle of attack (only for compressible solver).
  AddHistoryOutput("SENS_AOA",   "Sens_AoA",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the angle of attack (only for compressible solver).", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sensitivity of the objective function with respect to the Mach number (only of compressible solver).
  AddHistoryOutput("SENS_MACH",  "Sens_Mach",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the Mach number (only of compressible solver).", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sensitivity of the objective function with respect to the far-field pressure.
  AddHistoryOutput("SENS_PRESS", "Sens_Press", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the far-field pressure.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sensitivity of the objective function with respect to the far-field temperature.
  AddHistoryOutput("SENS_TEMP",  "Sens_Temp",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the far-field temperature.", HistoryFieldType::COEFFICIENT);
  /// END_GROUP

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

  if (config->GetDeform_Mesh()){
    AddHistoryOutput("DEFORM_ITER", "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", ScreenOutputFormat::FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");
  }
  }

}

void COneShotOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver){

  CSolver* flow_solver = solver[FLOW_SOL];
  CSolver* turb_solver = solver[TURB_SOL];
  CSolver* adjflow_solver = solver[ADJFLOW_SOL];
  CSolver* adjturb_solver = solver[ADJTURB_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];

  /*--- first group: residuals ---*/

  SetHistoryOutputValue("RMS_DENSITY", log10(flow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_MOMENTUM-X", log10(flow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_MOMENTUM-Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 2)
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(3)));
  else {
    SetHistoryOutputValue("RMS_MOMENTUM-Z", log10(flow_solver->GetRes_RMS(3)));
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(4)));
  }
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("RMS_NU_TILDE", log10(turb_solver->GetRes_RMS(0)));
    break;
  case SST: case SST_SUST:
    SetHistoryOutputValue("RMS_TKE", log10(turb_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_DISSIPATION",    log10(turb_solver->GetRes_RMS(1)));
    break;
  default: break;
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
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("MAX_NU_TILDE", log10(turb_solver->GetRes_Max(0)));
    break;
  case SST: case SST_SUST:
    SetHistoryOutputValue("MAX_TKE", log10(turb_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_DISSIPATION",    log10(turb_solver->GetRes_Max(1)));
    break;
  default: break;
  }

  SetHistoryOutputValue("RMS_ADJ_DENSITY", log10(adjflow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_RMS(2)));
  if (geometry->GetnDim() == 3) {
    SetHistoryOutputValue("RMS_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_RMS(3)));
    SetHistoryOutputValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(4)));
  } else {
    SetHistoryOutputValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(3)));
  }
  if ((!config->GetFrozen_Visc_Disc() && !cont_adj) || (!config->GetFrozen_Visc_Cont() && cont_adj)){
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      SetHistoryOutputValue("RMS_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_RMS(0)));
      break;
    case SST:
      SetHistoryOutputValue("RMS_ADJ_TKE", log10(adjturb_solver->GetRes_RMS(0)));
      SetHistoryOutputValue("RMS_ADJ_DISSIPATION",    log10(adjturb_solver->GetRes_RMS(1)));
      break;
    default: break;
    }
  }
  SetHistoryOutputValue("MAX_ADJ_DENSITY", log10(adjflow_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_Max(2)));
  if (geometry->GetnDim() == 3) {
    SetHistoryOutputValue("MAX_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_Max(3)));
    SetHistoryOutputValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(4)));
  } else {
    SetHistoryOutputValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(3)));
  }
  if ((!config->GetFrozen_Visc_Disc() && !cont_adj) || (!config->GetFrozen_Visc_Cont() && cont_adj)){
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      SetHistoryOutputValue("MAX_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_Max(0)));
      break;
    case SST:
      SetHistoryOutputValue("MAX_ADJ_TKE", log10(adjturb_solver->GetRes_Max(0)));
      SetHistoryOutputValue("MAX_ADJ_DISSIPATION", log10(adjturb_solver->GetRes_Max(1)));
      break;
    default: break;
    }
  }

  /*--- second group: residuals ---*/

  SetHistoryOutputValue("DRAG", flow_solver->GetTotal_CD());
  SetHistoryOutputValue("LIFT", flow_solver->GetTotal_CL());
  if (nDim == 3)
    SetHistoryOutputValue("SIDEFORCE", flow_solver->GetTotal_CSF());
  if (nDim == 3){
    SetHistoryOutputValue("MOMENT_X", flow_solver->GetTotal_CMx());
    SetHistoryOutputValue("MOMENT_Y", flow_solver->GetTotal_CMy());
  }
  SetHistoryOutputValue("MOMENT_Z", flow_solver->GetTotal_CMz());
  SetHistoryOutputValue("FORCE_X", flow_solver->GetTotal_CFx());
  SetHistoryOutputValue("FORCE_Y", flow_solver->GetTotal_CFy());
  if (nDim == 3)
    SetHistoryOutputValue("FORCE_Z", flow_solver->GetTotal_CFz());
  SetHistoryOutputValue("EFFICIENCY", flow_solver->GetTotal_CEff());

  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SetHistoryOutputPerSurfaceValue("DRAG_ON_SURFACE", flow_solver->GetSurface_CD(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("LIFT_ON_SURFACE", flow_solver->GetSurface_CL(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetHistoryOutputPerSurfaceValue("SIDEFORCE_ON_SURFACE", flow_solver->GetSurface_CSF(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3){
      SetHistoryOutputPerSurfaceValue("MOMENT-X_ON_SURFACE", flow_solver->GetSurface_CMx(iMarker_Monitoring), iMarker_Monitoring);
      SetHistoryOutputPerSurfaceValue("MOMENT-Y_ON_SURFACE", flow_solver->GetSurface_CMy(iMarker_Monitoring), iMarker_Monitoring);
    }
    SetHistoryOutputPerSurfaceValue("MOMENT-Z_ON_SURFACE", flow_solver->GetSurface_CMz(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("FORCE-X_ON_SURFACE", flow_solver->GetSurface_CFx(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("FORCE-Y_ON_SURFACE", flow_solver->GetSurface_CFy(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetHistoryOutputPerSurfaceValue("FORCE-Z_ON_SURFACE", flow_solver->GetSurface_CFz(iMarker_Monitoring), iMarker_Monitoring);

    SetHistoryOutputPerSurfaceValue("EFFICIENCY_ON_SURFACE", flow_solver->GetSurface_CEff(iMarker_Monitoring), iMarker_Monitoring);
    if (config->GetAeroelastic_Simulation()){
      SetHistoryOutputPerSurfaceValue("PITCH", config->GetAeroelastic_pitch(iMarker_Monitoring), iMarker_Monitoring);
      SetHistoryOutputPerSurfaceValue("PLUNGE", config->GetAeroelastic_plunge(iMarker_Monitoring), iMarker_Monitoring);
    }
  }

  SetHistoryOutputValue("AOA", config->GetAoA());

  /*--- third group: sensitivities ---*/

  SetHistoryOutputValue("SENS_GEO", adjflow_solver->GetTotal_Sens_Geo());
  SetHistoryOutputValue("SENS_AOA", adjflow_solver->GetTotal_Sens_AoA() * PI_NUMBER / 180.0);
  SetHistoryOutputValue("SENS_MACH", adjflow_solver->GetTotal_Sens_Mach());
  SetHistoryOutputValue("SENS_PRESS", adjflow_solver->GetTotal_Sens_Press());
  SetHistoryOutputValue("SENS_TEMP", adjflow_solver->GetTotal_Sens_Temp());

  SetHistoryOutputValue("LINSOL_ITER", adjflow_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(adjflow_solver->GetResLinSolver()));

  if (config->GetDeform_Mesh()) {
    SetHistoryOutputValue("DEFORM_ITER", mesh_solver->System.GetIterations());
    SetHistoryOutputValue("DEFORM_RESIDUAL", log10(mesh_solver->System.GetResidual()));
  }

}

void COneShotOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  /// BEGIN_GROUP: SOLUTION, DESCRIPTION: The SOLUTION variables of the adjoint solver.
  /// DESCRIPTION: Adjoint density.
  AddVolumeOutput("ADJ_DENSITY",    "Adjoint_Density",    "SOLUTION", "Adjoint density");
  /// DESCRIPTION: Adjoint momentum x-component.
  AddVolumeOutput("ADJ_MOMENTUM-X", "Adjoint_Momentum_x", "SOLUTION", "x-component of the adjoint momentum vector");
  /// DESCRIPTION: Adjoint momentum y-component.
  AddVolumeOutput("ADJ_MOMENTUM-Y", "Adjoint_Momentum_y", "SOLUTION", "y-component of the adjoint momentum vector");
  if (nDim == 3)
    /// DESCRIPTION: Adjoint momentum z-component.
    AddVolumeOutput("ADJ_MOMENTUM-Z", "Adjoint_Momentum_z", "SOLUTION", "z-component of the adjoint momentum vector");
  /// DESCRIPTION: Adjoint energy.
  AddVolumeOutput("ADJ_ENERGY", "Adjoint_Energy", "SOLUTION", "Adjoint energy");
  if ((!config->GetFrozen_Visc_Disc() && !cont_adj) || (!config->GetFrozen_Visc_Cont() && cont_adj)){
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      /// DESCRIPTION: Adjoint nu tilde.
      AddVolumeOutput("ADJ_NU_TILDE", "Adjoint_Nu_Tilde", "SOLUTION", "Adjoint Spalart-Allmaras variable");
      break;
    case SST:
      /// DESCRIPTION: Adjoint kinetic energy.
      AddVolumeOutput("ADJ_TKE", "Adjoint_TKE", "SOLUTION", "Adjoint kinetic energy");
      /// DESCRIPTION: Adjoint dissipation.
      AddVolumeOutput("ADJ_DISSIPATION", "Adjoint_Omega", "SOLUTION", "Adjoint rate of dissipation");
      break;
    default: break;
    }
  }
  /// END_GROUP

  /// BEGIN_GROUP: RESIDUAL, DESCRIPTION: Residuals of the SOLUTION variables.
  /// DESCRIPTION: Residual of the adjoint density.
  AddVolumeOutput("RES_ADJ_DENSITY",    "Residual_Adjoint_Density",    "RESIDUAL", "Residual of the adjoint density");
  /// DESCRIPTION: Residual of the adjoint momentum x-component.
  AddVolumeOutput("RES_ADJ_MOMENTUM-X", "Residual_Adjoint_Momentum_x", "RESIDUAL", "Residual of the adjoint x-momentum");
  /// DESCRIPTION: Residual of the adjoint momentum y-component.
  AddVolumeOutput("RES_ADJ_MOMENTUM-Y", "Residual_Adjoint_Momentum_y", "RESIDUAL", "Residual of the adjoint y-momentum");
  if (nDim == 3)
    /// DESCRIPTION: Residual of the adjoint momentum z-component.
    AddVolumeOutput("RES_ADJ_MOMENTUM-Z", "Residual_Adjoint_Momentum_z", "RESIDUAL", "Residual of the adjoint z-momentum");
  /// DESCRIPTION: Residual of the adjoint energy.
  AddVolumeOutput("RES_ADJ_ENERGY", "Residual_Adjoint_Energy", "RESIDUAL", "Residual of the adjoint energy");
  if ((!config->GetFrozen_Visc_Disc() && !cont_adj) || (!config->GetFrozen_Visc_Cont() && cont_adj)){
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      /// DESCRIPTION: Residual of the nu tilde.
      AddVolumeOutput("RES_ADJ_NU_TILDE", "Residual_Adjoint_Nu_Tilde", "RESIDUAL", "Residual of the Spalart-Allmaras variable");
      break;
    case SST:
      /// DESCRIPTION: Residual of the adjoint kinetic energy.
      AddVolumeOutput("RES_ADJ_TKE", "Residual_Adjoint_TKE", "RESIDUAL", "Residual of the turb. kinetic energy");
      /// DESCRIPTION: Residual of the adjoint dissipation.
      AddVolumeOutput("RES_ADJ_DISSIPATION", "Residual_Adjoint_Omega", "RESIDUAL", "Residual of the rate of dissipation");
      break;
    default: break;
    }
  }
  /// END_GROUP

  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Geometrical sensitivities of the current objective function.
  /// DESCRIPTION: Sensitivity x-component.
  AddVolumeOutput("SENSITIVITY-X", "Sensitivity_x", "SENSITIVITY", "x-component of the sensitivity vector");
  /// DESCRIPTION: Sensitivity y-component.
  AddVolumeOutput("SENSITIVITY-Y", "Sensitivity_y", "SENSITIVITY", "y-component of the sensitivity vector");
  if (nDim == 3)
    /// DESCRIPTION: Sensitivity z-component.
    AddVolumeOutput("SENSITIVITY-Z", "Sensitivity_z", "SENSITIVITY", "z-component of the sensitivity vector");
  /// DESCRIPTION: Sensitivity in normal direction.
  AddVolumeOutput("SENSITIVITY", "Surface_Sensitivity", "SENSITIVITY", "sensitivity in normal direction");
  /// END_GROUP
}

void COneShotOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_AdjFlow = solver[ADJFLOW_SOL]->GetNodes();
  CVariable* Node_AdjTurb = nullptr;
  CPoint*    Node_Geo     = geometry->nodes;

  if (config->GetKind_Turb_Model() != NONE &&
      ((!config->GetFrozen_Visc_Disc() && !cont_adj) ||
       (!config->GetFrozen_Visc_Cont() && cont_adj))){
    Node_AdjTurb = solver[ADJTURB_SOL]->GetNodes();
  }

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

  SetVolumeOutputValue("ADJ_DENSITY",    iPoint, Node_AdjFlow->GetSolution(iPoint, 0));
  SetVolumeOutputValue("ADJ_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(iPoint, 1));
  SetVolumeOutputValue("ADJ_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("ADJ_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(iPoint, 3));
    SetVolumeOutputValue("ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, 4));
  } else {
    SetVolumeOutputValue("ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, 3));
  }

  if ((!config->GetFrozen_Visc_Disc() && !cont_adj) || (!config->GetFrozen_Visc_Cont() && cont_adj)){
    // Turbulent
    switch(turb_model){
    case SST:
      SetVolumeOutputValue("ADJ_TKE",         iPoint, Node_AdjTurb->GetSolution(iPoint, 0));
      SetVolumeOutputValue("ADJ_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(iPoint, 1));
      break;
    case SA: case SA_COMP: case SA_E:
    case SA_E_COMP: case SA_NEG:
      SetVolumeOutputValue("ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(iPoint, 0));
      break;
    case NONE:
      break;
    }
  }

  // Residuals
  SetVolumeOutputValue("RES_ADJ_DENSITY",    iPoint, Node_AdjFlow->GetSolution(iPoint, 0) - Node_AdjFlow->GetSolution_Old(iPoint, 0));
  SetVolumeOutputValue("RES_ADJ_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(iPoint, 1) - Node_AdjFlow->GetSolution_Old(iPoint, 1));
  SetVolumeOutputValue("RES_ADJ_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(iPoint, 2) - Node_AdjFlow->GetSolution_Old(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("RES_ADJ_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
    SetVolumeOutputValue("RES_ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, 4) - Node_AdjFlow->GetSolution_Old(iPoint, 4));
  } else {
    SetVolumeOutputValue("RES_ADJ_ENERGY", iPoint, Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
  }

  if ((!config->GetFrozen_Visc_Disc() && !cont_adj) || (!config->GetFrozen_Visc_Cont() && cont_adj)){
    switch(config->GetKind_Turb_Model()){
    case SST:
      SetVolumeOutputValue("RES_ADJ_TKE",         iPoint, Node_AdjTurb->GetSolution(iPoint, 0) - Node_AdjTurb->GetSolution_Old(iPoint, 0));
      SetVolumeOutputValue("RES_ADJ_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(iPoint, 1) - Node_AdjTurb->GetSolution_Old(iPoint, 1));
      break;
    case SA: case SA_COMP: case SA_E:
    case SA_E_COMP: case SA_NEG:
      SetVolumeOutputValue("RES_ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(iPoint, 0) - Node_AdjTurb->GetSolution_Old(iPoint, 0));
      break;
    case NONE:
      break;
    }
  }

  SetVolumeOutputValue("SENSITIVITY-X", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 0));
  SetVolumeOutputValue("SENSITIVITY-Y", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENSITIVITY-Z", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 2));

}

void COneShotOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){

  SetVolumeOutputValue("SENSITIVITY", iPoint, solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex));

}

bool COneShotOutput::SetInit_Residuals(const CConfig *config){

  return ((config->GetTime_Marching() != TIME_MARCHING::STEADY) && (curInnerIter == 0)) ||
         ((config->GetTime_Marching() == TIME_MARCHING::STEADY) && (curInnerIter < 2));

}
