/*!
 * \file CAdjFlowIncOutput.cpp
 * \brief Main subroutines for flow discrete adjoint output
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


#include "../../include/output/CAdjFlowIncOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CAdjFlowIncOutput::CAdjFlowIncOutput(CConfig *config, unsigned short nDim) : CAdjFlowOutput(config, nDim) {

  heat = config->GetEnergy_Equation();

  weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

  rad_model = config->GetKind_RadiationModel();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0) {
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    requestedHistoryFields.emplace_back("SENSITIVITY");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }

  if (nRequestedScreenFields == 0) {
    if (config->GetTime_Domain()) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_ADJ_PRESSURE");
    requestedScreenFields.emplace_back("RMS_ADJ_VELOCITY-X");
    requestedScreenFields.emplace_back("SENS_GEO");
    requestedScreenFields.emplace_back("SENS_AOA");
    nRequestedScreenFields = requestedScreenFields.size();
  }

  if (nRequestedVolumeFields == 0) {
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

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_ADJ_PRESSURE");

}

CAdjFlowIncOutput::~CAdjFlowIncOutput() = default;

void CAdjFlowIncOutput::SetHistoryOutputFields(CConfig *config) {

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables.
  /// DESCRIPTION: Root-mean square residual of the adjoint Pressure.
  AddHistoryOutput("RMS_ADJ_PRESSURE",    "rms[A_P]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Pressure.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint Velocity x-component.
  AddHistoryOutput("RMS_ADJ_VELOCITY-X", "rms[A_U]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Velocity x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint Velocity y-component.
  AddHistoryOutput("RMS_ADJ_VELOCITY-Y", "rms[A_V]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Velocity y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint Velocity z-component.
  if (nDim == 3) {
    AddHistoryOutput("RMS_ADJ_VELOCITY-Z", "rms[A_W]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Velocity z-component.", HistoryFieldType::RESIDUAL);
  }
  /// DESCRIPTION: Maximum residual of the temperature.
  if (heat || weakly_coupled_heat) {
    AddHistoryOutput("RMS_ADJ_TEMPERATURE", "rms[A_T]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint temperature.", HistoryFieldType::RESIDUAL);
  }

  if (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
    AddHistoryOutput("ADJOINT_SOLEXTRA", "Adjoint_SolExtra", ScreenOutputFormat::FIXED, "ADJOINT_SOLEXTRA", "Adjoint value of the first extra Solution.", HistoryFieldType::COEFFICIENT);
  }

  AddHistoryOutputFields_AdjScalarRMS_RES(config);

  if (config->AddRadiation()) {
    /// DESCRIPTION: Root-mean square residual of the adjoint radiative energy tilde.
    AddHistoryOutput("RMS_ADJ_RAD_ENERGY", "rms[A_P1]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the P1 radiative energy.",HistoryFieldType::RESIDUAL);
  }
  /// END_GROUP

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the adjoint Pressure.
  AddHistoryOutput("MAX_ADJ_PRESSURE",    "max[A_Rho]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Pressure.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint Velocity x-component
  AddHistoryOutput("MAX_ADJ_VELOCITY-X", "max[A_RhoU]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Velocity x-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint Velocity y-component
  AddHistoryOutput("MAX_ADJ_VELOCITY-Y", "max[A_RhoV]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Velocity y-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint Velocity z-component
  if (nDim == 3) {
    AddHistoryOutput("MAX_ADJ_VELOCITY-Z", "max[A_RhoW]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Velocity z-component", HistoryFieldType::RESIDUAL);
  }
  /// DESCRIPTION: Maximum residual of the temperature.
  if (heat || weakly_coupled_heat) {
    AddHistoryOutput("MAX_ADJ_TEMPERATURE", "max[A_T]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the temperature.", HistoryFieldType::RESIDUAL);
  }

  AddHistoryOutputFields_AdjScalarMAX_RES(config);
  /// END_GROUP

  /// BEGIN_GROUP: BGS_RES, DESCRIPTION: The block Gauss Seidel residuals of the SOLUTION variables.
  /// DESCRIPTION: BGS residual of the adjoint Pressure.
  AddHistoryOutput("BGS_ADJ_PRESSURE",    "bgs[A_Rho]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Pressure.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: BGS residual of the adjoint Velocity x-component
  AddHistoryOutput("BGS_ADJ_VELOCITY-X", "bgs[A_RhoU]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Velocity x-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: BGS residual of the adjoint Velocity y-component
  AddHistoryOutput("BGS_ADJ_VELOCITY-Y", "bgs[A_RhoV]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Velocity y-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: BGS residual of the adjoint Velocity z-component
  if (nDim == 3) {
    AddHistoryOutput("BGS_ADJ_VELOCITY-Z", "bgs[A_RhoW]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Velocity z-component", HistoryFieldType::RESIDUAL);
  }
  /// DESCRIPTION: BGS residual of the temperature.
  if (heat || weakly_coupled_heat) {
    AddHistoryOutput("BGS_ADJ_TEMPERATURE", "bgs[A_T]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint temperature.", HistoryFieldType::RESIDUAL);
  }

  AddHistoryOutputFields_AdjScalarBGS_RES(config);

  if (config->AddRadiation()) {
    /// DESCRIPTION: Root-mean square residual of the adjoint radiative energy tilde.
    AddHistoryOutput("BGS_ADJ_RAD_ENERGY", "bgs[A_P1]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual  of the P1 radiative energy.",HistoryFieldType::RESIDUAL);
  }
  /// END_GROUP

  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Sensitivities of different geometrical or boundary values.
  /// DESCRIPTION: Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.
  AddHistoryOutput("SENS_GEO",   "Sens_Geo",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sensitivity of the objective function with respect to the far-field pressure.
  AddHistoryOutput("SENS_PRESS", "Sens_Press", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the far-field pressure.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sensitivity of the objective function with respect to the far-field temperature.
  AddHistoryOutput("SENS_TEMP",  "Sens_Temp",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", " Sensitivity of the objective function with respect to the far-field temperature.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sensitivity of the objective function with respect to the inlet velocity.
  AddHistoryOutput("SENS_VEL_IN", "Sens_Vin", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", " Sensitivity of the objective function with respect to the inlet velocity.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sensitivity of the objective function with respect to the outlet pressure.
  AddHistoryOutput("SENS_PRESS_OUT",  "Sens_Pout",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the outlet pressure.", HistoryFieldType::COEFFICIENT);
  /// END_GROUP

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

  AddHistoryOutputFieldsAdjScalarLinsol(config);

  if (config->GetDeform_Mesh()) {
    AddHistoryOutput("DEFORM_ITER", "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", ScreenOutputFormat::FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");
  }

}

void CAdjFlowIncOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  CSolver* adjflow_solver = solver[ADJFLOW_SOL];
  CSolver* adjheat_solver = solver[ADJHEAT_SOL];
  CSolver* adjrad_solver = solver[ADJRAD_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];

  SetHistoryOutputValue("RMS_ADJ_PRESSURE", log10(adjflow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_ADJ_VELOCITY-X", log10(adjflow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_ADJ_VELOCITY-Y", log10(adjflow_solver->GetRes_RMS(2)));
  if (nDim == 3) {
    SetHistoryOutputValue("RMS_ADJ_VELOCITY-Z", log10(adjflow_solver->GetRes_RMS(3)));
  }
  if (weakly_coupled_heat) {
    SetHistoryOutputValue("RMS_ADJ_TEMPERATURE",         log10(adjheat_solver->GetRes_RMS(0)));
  }
  if (heat) {
    if (nDim == 3) SetHistoryOutputValue("RMS_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_RMS(4)));
    else           SetHistoryOutputValue("RMS_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_RMS(3)));
  }

  if (config->AddRadiation()) {
    SetHistoryOutputValue("RMS_ADJ_RAD_ENERGY", log10(adjrad_solver->GetRes_RMS(0)));
  }

  if (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
    SetHistoryOutputValue("ADJOINT_SOLEXTRA", adjflow_solver->GetNodes()->GetSolutionExtra()[0]);
  }

  SetHistoryOutputValue("MAX_ADJ_PRESSURE", log10(adjflow_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_ADJ_VELOCITY-X", log10(adjflow_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_ADJ_VELOCITY-Y", log10(adjflow_solver->GetRes_Max(2)));
  if (nDim == 3) {
    SetHistoryOutputValue("MAX_ADJ_VELOCITY-Z", log10(adjflow_solver->GetRes_Max(3)));
  }
  if (weakly_coupled_heat) {
    SetHistoryOutputValue("MAX_ADJ_TEMPERATURE",         log10(adjheat_solver->GetRes_Max(0)));
  }
  if (heat) {
    if (nDim == 3) SetHistoryOutputValue("MAX_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_Max(4)));
    else           SetHistoryOutputValue("MAX_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_Max(3)));
  }

  if (multiZone) {
    SetHistoryOutputValue("BGS_ADJ_PRESSURE", log10(adjflow_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_ADJ_VELOCITY-X", log10(adjflow_solver->GetRes_BGS(1)));
    SetHistoryOutputValue("BGS_ADJ_VELOCITY-Y", log10(adjflow_solver->GetRes_BGS(2)));
    if (nDim == 3) {
      SetHistoryOutputValue("BGS_ADJ_VELOCITY-Z", log10(adjflow_solver->GetRes_BGS(3)));
    }
    if (weakly_coupled_heat) {
      SetHistoryOutputValue("BGS_ADJ_TEMPERATURE",         log10(adjheat_solver->GetRes_BGS(0)));
    }
    if (heat) {
      if (nDim == 3) SetHistoryOutputValue("BGS_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_BGS(4)));
      else           SetHistoryOutputValue("BGS_ADJ_TEMPERATURE",         log10(adjflow_solver->GetRes_BGS(3)));
    }

    if (config->AddRadiation()) {
      SetHistoryOutputValue("BGS_ADJ_RAD_ENERGY", log10(adjrad_solver->GetRes_BGS(0)));
    }
  }

  SetHistoryOutputValue("SENS_GEO", adjflow_solver->GetTotal_Sens_Geo());
  SetHistoryOutputValue("SENS_PRESS", adjflow_solver->GetTotal_Sens_Press());
  SetHistoryOutputValue("SENS_TEMP", adjflow_solver->GetTotal_Sens_Temp());
  SetHistoryOutputValue("SENS_VEL_IN", adjflow_solver->GetTotal_Sens_ModVel());
  SetHistoryOutputValue("SENS_PRESS_OUT", adjflow_solver->GetTotal_Sens_BPress());

  SetHistoryOutputValue("LINSOL_ITER", adjflow_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(adjflow_solver->GetResLinSolver()));

  if (config->GetDeform_Mesh()) {
    SetHistoryOutputValue("DEFORM_ITER", mesh_solver->System.GetIterations());
    SetHistoryOutputValue("DEFORM_RESIDUAL", log10(mesh_solver->System.GetResidual()));
  }

  LoadHistoryDataAdjScalar(config, solver);

  ComputeSimpleCustomOutputs(config);
}

void CAdjFlowIncOutput::SetVolumeOutputFields(CConfig *config) {


  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3) {
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");
  }
  /// BEGIN_GROUP: SOLUTION, DESCRIPTION: The SOLUTION variables of the adjoint solver.
  /// DESCRIPTION: Adjoint Pressure.
  AddVolumeOutput("ADJ_PRESSURE",    "Adjoint_Pressure",    "SOLUTION", "Adjoint pressure");
  /// DESCRIPTION: Adjoint Velocity x-component.
  AddVolumeOutput("ADJ_VELOCITY-X", "Adjoint_Velocity_x", "SOLUTION", "x-component of the adjoint velocity vector");
  /// DESCRIPTION: Adjoint Velocity y-component.
  AddVolumeOutput("ADJ_VELOCITY-Y", "Adjoint_Velocity_y", "SOLUTION", "y-component of the adjoint velocity vector");
  if (nDim == 3) {
    /// DESCRIPTION: Adjoint Velocity z-component.
    AddVolumeOutput("ADJ_VELOCITY-Z", "Adjoint_Velocity_z", "SOLUTION", "z-component of the adjoint velocity vector");
  }
  AddVolumeOutput("ADJ_TEMPERATURE", "Adjoint_Temperature", "SOLUTION",  "Adjoint temperature");

  SetVolumeOutputFieldsAdjScalarSolution(config);

  if (config->AddRadiation()) {
    AddVolumeOutput("ADJ_P1_ENERGY",  "Adjoint_Energy(P1)", "SOLUTION", "Adjoint radiative energy");
  }
  /// END_GROUP

  // Grid velocity
  if (config->GetDynamic_Grid()) {
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY", "x-component of the grid velocity vector");
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY", "y-component of the grid velocity vector");
    if (nDim == 3 ) {
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY", "z-component of the grid velocity vector");
    }
  }

  /// BEGIN_GROUP: RESIDUAL, DESCRIPTION: Residuals of the SOLUTION variables.
  /// DESCRIPTION: Residual of the adjoint Pressure.
  AddVolumeOutput("RES_ADJ_PRESSURE",    "Residual_Adjoint_Pressure",    "RESIDUAL", "Residual of the adjoint pressure");
  /// DESCRIPTION: Residual of the adjoint Velocity x-component.
  AddVolumeOutput("RES_ADJ_VELOCITY-X", "Residual_Adjoint_Velocity_x", "RESIDUAL", "Residual of the adjoint x-velocity");
  /// DESCRIPTION: Residual of the adjoint Velocity y-component.
  AddVolumeOutput("RES_ADJ_VELOCITY-Y", "Residual_Adjoint_Velocity_y", "RESIDUAL", "Residual of the adjoint y-velocity");
  if (nDim == 3) {
    /// DESCRIPTION: Residual of the adjoint Velocity z-component.
    AddVolumeOutput("RES_ADJ_VELOCITY-Z", "Residual_Adjoint_Velocity_z", "RESIDUAL", "Residual of the adjoint z-velocity");
  }
  /// DESCRIPTION: Residual of the adjoint energy.
  AddVolumeOutput("RES_ADJ_TEMPERATURE", "Residual_Adjoint_Heat", "RESIDUAL", "Residual of the adjoint temperature");

  SetVolumeOutputFieldsAdjScalarResidual(config);

  if (config->AddRadiation()) {
    AddVolumeOutput("RES_P1_ENERGY",  "Residual_Adjoint_Energy_P1", "RESIDUAL", "Residual of adjoint radiative energy");
  }
  /// END_GROUP

  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Geometrical sensitivities of the current objective function.
  /// DESCRIPTION: Sensitivity x-component.
  AddVolumeOutput("SENSITIVITY-X", "Sensitivity_x", "SENSITIVITY", "x-component of the sensitivity vector");
  /// DESCRIPTION: Sensitivity y-component.
  AddVolumeOutput("SENSITIVITY-Y", "Sensitivity_y", "SENSITIVITY", "y-component of the sensitivity vector");
  if (nDim == 3) {
    /// DESCRIPTION: Sensitivity z-component.
    AddVolumeOutput("SENSITIVITY-Z", "Sensitivity_z", "SENSITIVITY", "z-component of the sensitivity vector");
  }
  /// DESCRIPTION: Sensitivity in normal direction.
  AddVolumeOutput("SENSITIVITY", "Surface_Sensitivity", "SENSITIVITY", "sensitivity in normal direction");
  /// END_GROUP

}

void CAdjFlowIncOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint) {

  CVariable* Node_AdjFlow = solver[ADJFLOW_SOL]->GetNodes();
  CVariable* Node_AdjHeat = nullptr;
  CVariable* Node_AdjRad  = nullptr;
  CPoint*    Node_Geo     = geometry->nodes;

  if (weakly_coupled_heat) {
    Node_AdjHeat = solver[ADJHEAT_SOL]->GetNodes();
  }
  if (config->GetKind_RadiationModel() != RADIATION_MODEL::NONE) {
    Node_AdjRad = solver[ADJRAD_SOL]->GetNodes();
  }

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3) {
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));
  }

  SetVolumeOutputValue("ADJ_PRESSURE",   iPoint, Node_AdjFlow->GetSolution(iPoint, 0));
  SetVolumeOutputValue("ADJ_VELOCITY-X", iPoint, Node_AdjFlow->GetSolution(iPoint, 1));
  SetVolumeOutputValue("ADJ_VELOCITY-Y", iPoint, Node_AdjFlow->GetSolution(iPoint, 2));
  if (nDim == 3) {
    SetVolumeOutputValue("ADJ_VELOCITY-Z", iPoint, Node_AdjFlow->GetSolution(iPoint, 3));
  }

  if (weakly_coupled_heat) {
    SetVolumeOutputValue("ADJ_TEMPERATURE", iPoint, Node_AdjHeat->GetSolution(iPoint, 0));
  }
  else {
    if (nDim == 3) SetVolumeOutputValue("ADJ_TEMPERATURE", iPoint, Node_AdjFlow->GetSolution(iPoint, 4));
    else           SetVolumeOutputValue("ADJ_TEMPERATURE", iPoint, Node_AdjFlow->GetSolution(iPoint, 3));
  }

  // Radiation
  if (config->AddRadiation()) {
    SetVolumeOutputValue("ADJ_P1_ENERGY", iPoint, Node_AdjRad->GetSolution(iPoint, 0));
  }

  // Residuals
  SetVolumeOutputValue("RES_ADJ_PRESSURE",   iPoint, Node_AdjFlow->GetSolution(iPoint, 0) - Node_AdjFlow->GetSolution_Old(iPoint, 0));
  SetVolumeOutputValue("RES_ADJ_VELOCITY-X", iPoint, Node_AdjFlow->GetSolution(iPoint, 1) - Node_AdjFlow->GetSolution_Old(iPoint, 1));
  SetVolumeOutputValue("RES_ADJ_VELOCITY-Y", iPoint, Node_AdjFlow->GetSolution(iPoint, 2) - Node_AdjFlow->GetSolution_Old(iPoint, 2));
  if (nDim == 3) {
    SetVolumeOutputValue("RES_ADJ_VELOCITY-Z", iPoint, Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
    SetVolumeOutputValue("RES_ADJ_TEMPERATURE",     iPoint, Node_AdjFlow->GetSolution(iPoint, 4) - Node_AdjFlow->GetSolution_Old(iPoint, 4));
  } else {
    SetVolumeOutputValue("RES_ADJ_TEMPERATURE",     iPoint, Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
  }

  if (config->AddRadiation()) {
    SetVolumeOutputValue("RES_P1_ENERGY", iPoint, Node_AdjRad->GetSolution(iPoint, 0) - Node_AdjRad->GetSolution_Old(iPoint, 0));
  }

  SetVolumeOutputValue("SENSITIVITY-X", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 0));
  SetVolumeOutputValue("SENSITIVITY-Y", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 1));
  if (nDim == 3) {
    SetVolumeOutputValue("SENSITIVITY-Z", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 2));
  }

  LoadVolumeDataAdjScalar(config, solver, iPoint);
}

void CAdjFlowIncOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex) {

  SetVolumeOutputValue("SENSITIVITY", iPoint, solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex));

}


bool CAdjFlowIncOutput::SetInitResiduals(const CConfig *config) {

  return (config->GetTime_Marching() != TIME_MARCHING::STEADY && (curInnerIter == 0))||
         (config->GetTime_Marching() == TIME_MARCHING::STEADY && (curTimeIter < 2));

}
