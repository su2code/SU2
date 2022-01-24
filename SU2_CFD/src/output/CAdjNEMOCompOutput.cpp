/*!
 * \file CAdjNEMOCompOutput.cpp
 * \brief Main subroutines for flow discrete adjoint output
 * \author W.Maier, R. Sanchez
 * \version 7.3.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/output/CAdjNEMOCompOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CAdjNEMOCompOutput::CAdjNEMOCompOutput(CConfig *config, unsigned short nDim) : CAdjFlowOutput(config, nDim) {


  unsigned short nSpecies   = config->GetnSpecies();



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
    for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      requestedScreenFields.emplace_back("RMS_ADJ_DENSITY_" + std::to_string(iSpecies));
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

  volumeFilename = config->GetAdj_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfAdjCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_AdjFileName();

  /*--- Add the obj. function extension --- */

  restartFilename = config->GetObjFunc_Extension(restartFilename);

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_ADJ_DENSITY_0");

}

CAdjNEMOCompOutput::~CAdjNEMOCompOutput(void) {}

void CAdjNEMOCompOutput::SetHistoryOutputFields(CConfig *config){

  unsigned short nSpecies = config->GetnSpecies();

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables.
  /// DESCRIPTION: Root-mean square residual of the adjoint density.
  for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    AddHistoryOutput("RMS_ADJ_DENSITY_" + std::to_string(iSpecies), "rms[A_Rho_" + std::to_string(iSpecies) + "]",   ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the species adjoint density " + std::to_string(iSpecies) + ".", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum x-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-X", "rms[A_RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum y-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-Y", "rms[A_RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum z-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-Z", "rms[A_RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint energy.
  AddHistoryOutput("RMS_ADJ_ENERGY",     "rms[A_E]",    ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint energy.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the energy.
  AddHistoryOutput("RMS_ADJ_ENERGY_VE",  "rms[A_Eve]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint ve-energy.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_AdjScalarRMS_RES(config);
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
  /// DESCRIPTION: Maximum residual of the adjoint ve-energy.
  AddHistoryOutput("MAX_ADJ_ENERGY_VE",     "max[A_Eve]",    ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint ve-energy.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_AdjScalarMAX_RES(config);
  /// END_GROUP


  ///  /// BEGIN_GROUP: BGS_RES, DESCRIPTION: The Block Gauss Seidel residuals of the SOLUTION variables.
  /// DESCRIPTION: BGS residual of the adjoint density.
  AddHistoryOutput("BGS_ADJ_DENSITY",    "bgs[A_Rho]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: BGS residual of the adjoint momentum x-component
  AddHistoryOutput("BGS_ADJ_MOMENTUM-X", "bgs[A_RhoU]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint momentum x-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: BGS residual of the adjoint momentum y-component
  AddHistoryOutput("BGS_ADJ_MOMENTUM-Y", "bgs[A_RhoV]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint momentum y-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: BGS residual of the adjoint momentum z-component
  AddHistoryOutput("BGS_ADJ_MOMENTUM-Z", "bgs[A_RhoW]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint momentum z-component", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: BGS residual of the adjoint energy.
  AddHistoryOutput("BGS_ADJ_ENERGY",     "bgs[A_E]",    ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint energy.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: BGS residual of the adjoint ve-energy.
  AddHistoryOutput("BGS_ADJ_ENERGY_VE",     "bgs[A_Eve]",    ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint ve-energy.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_AdjScalarBGS_RES(config);
  /// END_GROUP

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

  AddHistoryOutputFields_AdjScalarLinsol(config);

  if (config->GetDeform_Mesh()){
    AddHistoryOutput("DEFORM_ITER", "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", ScreenOutputFormat::FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");
  }

}

void CAdjNEMOCompOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver){

  CSolver* adjflow_solver = solver[ADJFLOW_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];

  unsigned short nSpecies = config->GetnSpecies();

  for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    SetHistoryOutputValue("RMS_ADJ_DENSITY_" + std::to_string(iSpecies), log10(adjflow_solver->GetRes_RMS(iSpecies)));

  SetHistoryOutputValue("RMS_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_RMS(nSpecies)));
  SetHistoryOutputValue("RMS_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_RMS(nSpecies+1)));
  if (geometry->GetnDim() == 3) {
    SetHistoryOutputValue("RMS_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_RMS(nSpecies+2)));
    SetHistoryOutputValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(nSpecies+3)));
    SetHistoryOutputValue("RMS_ADJ_ENERGY_VE", log10(adjflow_solver->GetRes_RMS(nSpecies+4)));
  } else {
    SetHistoryOutputValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(nSpecies+2)));
    SetHistoryOutputValue("RMS_ADJ_ENERGY_VE", log10(adjflow_solver->GetRes_RMS(nSpecies+3)));

  }

  for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    SetHistoryOutputValue("MAX_ADJ_DENSITY_" + std::to_string(iSpecies), log10(adjflow_solver->GetRes_RMS(iSpecies)));
  SetHistoryOutputValue("MAX_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_Max(nSpecies)));
  SetHistoryOutputValue("MAX_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_Max(nSpecies+1)));
  if (geometry->GetnDim() == 3) {
    SetHistoryOutputValue("MAX_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_Max(nSpecies+2)));
    SetHistoryOutputValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(nSpecies+3)));
    SetHistoryOutputValue("MAX_ADJ_ENERGY_VE", log10(adjflow_solver->GetRes_Max(nSpecies+4)));
  } else {
    SetHistoryOutputValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(nSpecies+2)));
    SetHistoryOutputValue("MAX_ADJ_ENERGY_VE", log10(adjflow_solver->GetRes_Max(nSpecies+3)));
  }

  if (multiZone){
    for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      SetHistoryOutputValue("BGS_ADJ_DENSITY_" + std::to_string(iSpecies), log10(adjflow_solver->GetRes_BGS(iSpecies)));
    SetHistoryOutputValue("BGS_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_BGS(nSpecies)));
    SetHistoryOutputValue("BGS_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_BGS(nSpecies+1)));
    if (geometry->GetnDim() == 3) {
      SetHistoryOutputValue("BGS_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_BGS(nSpecies+2)));
      SetHistoryOutputValue("BGS_ADJ_ENERGY", log10(adjflow_solver->GetRes_BGS(nSpecies+3)));
      SetHistoryOutputValue("BGS_ADJ_ENERGY_VE", log10(adjflow_solver->GetRes_BGS(nSpecies+4)));
    } else {
      SetHistoryOutputValue("BGS_ADJ_ENERGY", log10(adjflow_solver->GetRes_BGS(nSpecies+2)));
      SetHistoryOutputValue("BGS_ADJ_ENERGY_VE", log10(adjflow_solver->GetRes_BGS(nSpecies+3)));    }
  }

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

  LoadHistoryData_AdjScalar(config, solver);
}

void CAdjNEMOCompOutput::SetVolumeOutputFields(CConfig *config){

  unsigned short nSpecies = config->GetnSpecies();

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  /// BEGIN_GROUP: SOLUTION, DESCRIPTION: The SOLUTION variables of the adjoint solver.
  /// DESCRIPTION: Adjoint density.
  for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    AddVolumeOutput("ADJ_DENSITY_" + std::to_string(iSpecies),  "Density_" + std::to_string(iSpecies),  "SOLUTION", "Density_"  + std::to_string(iSpecies));
  /// DESCRIPTION: Adjoint momentum x-component.
  AddVolumeOutput("ADJ_MOMENTUM-X", "Adjoint_Momentum_x", "SOLUTION", "x-component of the adjoint momentum vector");
  /// DESCRIPTION: Adjoint momentum y-component.
  AddVolumeOutput("ADJ_MOMENTUM-Y", "Adjoint_Momentum_y", "SOLUTION", "y-component of the adjoint momentum vector");
  if (nDim == 3)
    /// DESCRIPTION: Adjoint momentum z-component.
    AddVolumeOutput("ADJ_MOMENTUM-Z", "Adjoint_Momentum_z", "SOLUTION", "z-component of the adjoint momentum vector");
  /// DESCRIPTION: Adjoint energy.
  AddVolumeOutput("ADJ_ENERGY", "Adjoint_Energy", "SOLUTION", "Adjoint energy");
  /// DESCRIPTION: Adjoint energy.
  AddVolumeOutput("ADJ_ENERGY_VE", "Adjoint_Energy_VE", "SOLUTION", "Adjoint ve-energy");

  SetVolumeOutputFields_AdjScalarSolution(config);
  /// END_GROUP

  /// BEGIN_GROUP: RESIDUAL, DESCRIPTION: Residuals of the SOLUTION variables.
  /// DESCRIPTION: Residual of the adjoint density.
  for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    AddVolumeOutput("RES_ADJ_DENSITY_" + std::to_string(iSpecies), "Residual_Density_" + std::to_string(iSpecies), "RESIDUAL", "Residual of species density " + std::to_string(iSpecies));
  /// DESCRIPTION: Residual of the adjoint momentum x-component.
  AddVolumeOutput("RES_ADJ_MOMENTUM-X", "Residual_Adjoint_Momentum_x", "RESIDUAL", "Residual of the adjoint x-momentum");
  /// DESCRIPTION: Residual of the adjoint momentum y-component.
  AddVolumeOutput("RES_ADJ_MOMENTUM-Y", "Residual_Adjoint_Momentum_y", "RESIDUAL", "Residual of the adjoint y-momentum");
  if (nDim == 3)
    /// DESCRIPTION: Residual of the adjoint momentum z-component.
    AddVolumeOutput("RES_ADJ_MOMENTUM-Z", "Residual_Adjoint_Momentum_z", "RESIDUAL", "Residual of the adjoint z-momentum");
  /// DESCRIPTION: Residual of the adjoint energy.
  AddVolumeOutput("RES_ADJ_ENERGY", "Residual_Adjoint_Energy", "RESIDUAL", "Residual of the adjoint energy");
  /// DESCRIPTION: Residual of the adjoint ve-energy.
  AddVolumeOutput("RES_ADJ_ENERGY_VE", "Residual_Adjoint_Energy_VE", "RESIDUAL", "Residual of the adjoint ve-energy");

  SetVolumeOutputFields_AdjScalarResidual(config);
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

void CAdjNEMOCompOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_AdjFlow = solver[ADJFLOW_SOL]->GetNodes();
  CPoint*    Node_Geo     = geometry->nodes;
  unsigned short nSpecies = config->GetnSpecies();


  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

  for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    SetVolumeOutputValue("ADJ_DENSITY_" + std::to_string(iSpecies), iPoint, Node_AdjFlow->GetSolution(iPoint, iSpecies));

  SetVolumeOutputValue("ADJ_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies));
  SetVolumeOutputValue("ADJ_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+1));
  if (nDim == 3){
    SetVolumeOutputValue("ADJ_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+2));
    SetVolumeOutputValue("ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+3));
    SetVolumeOutputValue("ADJ_ENERGY_VE",  iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+4));
  } else {
    SetVolumeOutputValue("ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+2));
    SetVolumeOutputValue("ADJ_ENERGY_VE",  iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+3));
  }


  // Residuals
  for (auto iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    SetVolumeOutputValue("RES_ADJ_DENSITY_" + std::to_string(iSpecies), iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, iSpecies));

  SetVolumeOutputValue("RES_ADJ_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies) - Node_AdjFlow->GetSolution_Old(iPoint, nSpecies));
  SetVolumeOutputValue("RES_ADJ_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+1) - Node_AdjFlow->GetSolution_Old(iPoint, nSpecies+1));
  if (nDim == 3){
    SetVolumeOutputValue("RES_ADJ_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+2) - Node_AdjFlow->GetSolution_Old(iPoint, nSpecies+2));
    SetVolumeOutputValue("RES_ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+3) - Node_AdjFlow->GetSolution_Old(iPoint, nSpecies+3));
    SetVolumeOutputValue("RES_ADJ_ENERGY_VE",  iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+4) - Node_AdjFlow->GetSolution_Old(iPoint, nSpecies+4));
  } else {
    SetVolumeOutputValue("RES_ADJ_ENERGY", iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+2) - Node_AdjFlow->GetSolution_Old(iPoint, nSpecies+2));
    SetVolumeOutputValue("RES_ADJ_ENERGY_VE", iPoint, Node_AdjFlow->GetSolution(iPoint, nSpecies+3) - Node_AdjFlow->GetSolution_Old(iPoint, nSpecies+3));
  }


  SetVolumeOutputValue("SENSITIVITY-X", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 0));
  SetVolumeOutputValue("SENSITIVITY-Y", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENSITIVITY-Z", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 2));

  LoadVolumeData_AdjScalar(config, solver, iPoint);
}

void CAdjNEMOCompOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){

  SetVolumeOutputValue("SENSITIVITY", iPoint, solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex));

}


bool CAdjNEMOCompOutput::SetInit_Residuals(const CConfig *config){

  return ((config->GetTime_Marching() != TIME_MARCHING::STEADY) && (curInnerIter == 0)) ||
         ((config->GetTime_Marching() == TIME_MARCHING::STEADY) && (curInnerIter < 2));

}

