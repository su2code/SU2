/*!
 * \file CAdjHeatOutput.cpp
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


#include "../../include/output/CAdjHeatOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CAdjHeatOutput::CAdjHeatOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    requestedHistoryFields.emplace_back("SENSITIVITY");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }

  if (nRequestedScreenFields == 0){
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_ADJ_TEMPERATURE");
    requestedScreenFields.emplace_back("SENS_GEO");
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
  ss << "Zone " << config->GetiZone() << " (Adj. Heat)";
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

  if (convFields.empty() ) convFields.emplace_back("RMS_ADJ_TEMPERATURE");

}

CAdjHeatOutput::~CAdjHeatOutput() = default;

void CAdjHeatOutput::SetHistoryOutputFields(CConfig *config){

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the conservative variables.
  /// DESCRIPTION: Root-mean square residual of the adjoint temperature.
  AddHistoryOutput("RMS_ADJ_TEMPERATURE",    "rms[A_T]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint temperature.", HistoryFieldType::RESIDUAL);
  /// END_GROUP

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the conservative variables.
  /// DESCRIPTION: Maximum residual of the adjoint temperature.
  AddHistoryOutput("MAX_ADJ_TEMPERATURE",    "max[A_T]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint temperature.", HistoryFieldType::RESIDUAL);

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The root-mean-square residuals of the conservative variables.
  /// DESCRIPTION: Root-mean-square residual of the adjoint temperature.
  AddHistoryOutput("BGS_ADJ_TEMPERATURE",    "bgs[A_T]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint temperature.", HistoryFieldType::RESIDUAL);

  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Sensitivities of different geometrical or boundary values.
  /// DESCRIPTION: Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.
  AddHistoryOutput("SENS_GEO",   "Sens_Geo",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// END_GROUP

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

  if (config->GetDeform_Mesh()){
    AddHistoryOutput("DEFORM_ITER", "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", ScreenOutputFormat::FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");
  }

}

void CAdjHeatOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  CSolver* adjheat_solver = solver[ADJHEAT_SOL];

  SetHistoryOutputValue("RMS_ADJ_TEMPERATURE", log10(adjheat_solver->GetRes_RMS(0)));

  SetHistoryOutputValue("MAX_ADJ_TEMPERATURE", log10(adjheat_solver->GetRes_Max(0)));

  if (multiZone) {
    SetHistoryOutputValue("BGS_ADJ_TEMPERATURE", log10(adjheat_solver->GetRes_BGS(0)));
  }

  SetHistoryOutputValue("SENS_GEO", adjheat_solver->GetTotal_Sens_Geo());

  SetHistoryOutputValue("LINSOL_ITER", adjheat_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(adjheat_solver->GetResLinSolver()));

  if (config->GetDeform_Mesh()) {
    SetHistoryOutputValue("DEFORM_ITER", solver[MESH_SOL]->System.GetIterations());
    SetHistoryOutputValue("DEFORM_RESIDUAL", log10(solver[MESH_SOL]->System.GetResidual()));
  }

  ComputeSimpleCustomOutputs(config);
}

void CAdjHeatOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");


  /// BEGIN_GROUP: CONSERVATIVE, DESCRIPTION: The conservative variables of the adjoint solver.
  /// DESCRIPTION: Adjoint Pressure.
  AddVolumeOutput("ADJ_TEMPERATURE",    "Adjoint_Temperature",    "SOLUTION" ,"Adjoint Temperature");
  /// END_GROUP


  /// BEGIN_GROUP: RESIDUAL, DESCRIPTION: Residuals of the conservative variables.
  /// DESCRIPTION: Residual of the adjoint Pressure.
  AddVolumeOutput("RES_ADJ_TEMPERATURE",    "Residual_Adjoint_Temperature",    "RESIDUAL", "Residual of the Adjoint Temperature");
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

  if (!config->GetTime_Domain()) return;

  /*--- Sensitivities with respect to initial conditions. ---*/

  AddVolumeOutput("SENS_TEMP_N", "SensitivityTempN", "SENSITIVITY_N", "sensitivity to the previous temperature");
  if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) {
    AddVolumeOutput("SENS_TEMP_N1", "SensitivityTempN1", "SENSITIVITY_N", "sensitivity to the previous-1 temperature");
  }
}

void CAdjHeatOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_AdjHeat = solver[ADJHEAT_SOL]->GetNodes();
  CPoint*    Node_Geo     = geometry->nodes;


  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

  SetVolumeOutputValue("ADJ_TEMPERATURE", iPoint, Node_AdjHeat->GetSolution(iPoint, 0));

  // Residuals
  SetVolumeOutputValue("RES_ADJ_TEMPERATURE", iPoint, Node_AdjHeat->GetSolution(iPoint, 0) - Node_AdjHeat->GetSolution_Old(iPoint, 0));

  SetVolumeOutputValue("SENSITIVITY-X", iPoint, Node_AdjHeat->GetSensitivity(iPoint, 0));
  SetVolumeOutputValue("SENSITIVITY-Y", iPoint, Node_AdjHeat->GetSensitivity(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENSITIVITY-Z", iPoint, Node_AdjHeat->GetSensitivity(iPoint, 2));

  if (!config->GetTime_Domain()) return;

  SetVolumeOutputValue("SENS_TEMP_N", iPoint, Node_AdjHeat->GetSolution_time_n(iPoint, 0));
  if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) {
    SetVolumeOutputValue("SENS_TEMP_N1", iPoint, Node_AdjHeat->GetSolution_time_n1(iPoint, 0));
  }

}

void CAdjHeatOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){

  SetVolumeOutputValue("SENSITIVITY", iPoint, solver[ADJHEAT_SOL]->GetCSensitivity(iMarker, iVertex));

}

