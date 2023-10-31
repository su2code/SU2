/*!
 * \file CAdjElasticityOutput.cpp
 * \brief Main subroutines for elasticity discrete adjoint output
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


#include "../../include/output/CAdjElasticityOutput.hpp"
#include <string>

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CAdjElasticityOutput::CAdjElasticityOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  /*--- Initialize number of variables ---*/
  nVar_FEM = nDim;

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RESIDUALS");
    requestedHistoryFields.emplace_back("SENSITIVITY");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }

  if (nRequestedScreenFields == 0){
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("ADJOINT_DISP_X");
    requestedScreenFields.emplace_back("ADJOINT_DISP_Y");
    requestedScreenFields.emplace_back("SENS_E_0");
    requestedScreenFields.emplace_back("SENS_NU_0");
    nRequestedScreenFields = requestedScreenFields.size();
  }

  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("SENSITIVITY");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  if (find(requestedVolumeFields.begin(), requestedVolumeFields.end(), "SENSITIVITY") == requestedVolumeFields.end()) {
    requestedVolumeFields.emplace_back("SENSITIVITY");
    nRequestedVolumeFields++;
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Adj. Structure)";
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

  if (convFields.empty() ) convFields.emplace_back("ADJOINT_DISP_X");

}

CAdjElasticityOutput::~CAdjElasticityOutput() = default;

void CAdjElasticityOutput::SetHistoryOutputFields(CConfig *config){

  /*--- Residuals ---*/
  AddHistoryOutput("ADJOINT_DISP_X", "rms[Ux_adj]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint of the X displacements.", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("ADJOINT_DISP_Y", "rms[Uy_adj]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint of the Y displacements.", HistoryFieldType::RESIDUAL);
  if (nVar_FEM == 3) {
    AddHistoryOutput("ADJOINT_DISP_Z", "rms[Uz_adj]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint of the Z displacements.", HistoryFieldType::RESIDUAL);
  }

  /*--- Sensitivities ---*/
  for (auto iVar = 0u; iVar < config->GetnElasticityMat(); iVar++) {
    const auto iVarS = std::to_string(iVar);
    AddHistoryOutput("SENS_E_" + iVarS, "Sens[E" + iVarS + ']', ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "d Objective / d Elasticity modulus");
    AddHistoryOutput("SENS_NU_" + iVarS, "Sens[Nu" + iVarS + ']', ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "d Objective / d Poisson ratio");
    if (config->GetTime_Domain() && !config->GetPseudoStatic()) {
      AddHistoryOutput("SENS_RHO_" + iVarS, "Sens[Rho" + iVarS + ']', ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "d Objective / d Material density");
    }
    if (config->GetDeadLoad()) {
      AddHistoryOutput("SENS_RHO_DL_" + iVarS, "Sens[RhoDL" + iVarS + ']', ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "d Objective / d Dead load density");
    }
  }
  if (config->GetDE_Effects()) {
    for (auto iVar = 0u; iVar < config->GetnElectric_Field(); iVar++) {
      const auto iVarS = std::to_string(iVar);
      AddHistoryOutput("SENS_EFIELD_" + iVarS, "Sens[EField" + iVarS + ']', ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "d Objective / d Electric field");
    }
  }

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

  if (multiZone) {
    AddHistoryOutput("BGS_ADJ_DISP_X", "bgs[A_Ux]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint X displacement.", HistoryFieldType::RESIDUAL);
    AddHistoryOutput("BGS_ADJ_DISP_Y", "bgs[A_Uy]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Y displacement.", HistoryFieldType::RESIDUAL);
    if (nVar_FEM == 3) {
      AddHistoryOutput("BGS_ADJ_DISP_Z", "bgs[A_Uz]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint Z displacement.", HistoryFieldType::RESIDUAL);
    }
  }
}

inline void CAdjElasticityOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  SetHistoryOutputValue("ADJOINT_DISP_X", log10(solver[ADJFEA_SOL]->GetRes_RMS(0)));
  SetHistoryOutputValue("ADJOINT_DISP_Y", log10(solver[ADJFEA_SOL]->GetRes_RMS(1)));
  if (nVar_FEM == 3) {
    SetHistoryOutputValue("ADJOINT_DISP_Z", log10(solver[ADJFEA_SOL]->GetRes_RMS(2)));
  }

  for (unsigned short iVar = 0; iVar < config->GetnElasticityMat(); iVar++) {
    const auto iVarS = std::to_string(iVar);
    SetHistoryOutputValue("SENS_E_" + iVarS, solver[ADJFEA_SOL]->GetTotal_Sens_E(iVar));
    SetHistoryOutputValue("SENS_NU_" + iVarS, solver[ADJFEA_SOL]->GetTotal_Sens_Nu(iVar));
    if (config->GetTime_Domain() && !config->GetPseudoStatic()) {
      SetHistoryOutputValue("SENS_RHO_" + iVarS, solver[ADJFEA_SOL]->GetTotal_Sens_Rho(iVar));
    }
    if (config->GetDeadLoad()) {
      SetHistoryOutputValue("SENS_RHO_DL_" + iVarS, solver[ADJFEA_SOL]->GetTotal_Sens_Rho_DL(iVar));
    }
  }
  if (config->GetDE_Effects()) {
    for (auto iVar = 0u; iVar < config->GetnElectric_Field(); iVar++) {
      const auto iVarS = std::to_string(iVar);
      SetHistoryOutputValue("SENS_EFIELD_" + iVarS, solver[ADJFEA_SOL]->GetTotal_Sens_EField(iVar));
    }
  }

  SetHistoryOutputValue("LINSOL_ITER", solver[ADJFEA_SOL]->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(solver[ADJFEA_SOL]->GetResLinSolver()));

  if (multiZone) {
    SetHistoryOutputValue("BGS_ADJ_DISP_X", log10(solver[ADJFEA_SOL]->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_ADJ_DISP_Y", log10(solver[ADJFEA_SOL]->GetRes_BGS(1)));
    if (nVar_FEM == 3){
      SetHistoryOutputValue("BGS_ADJ_DISP_Z", log10(solver[ADJFEA_SOL]->GetRes_BGS(2)));
    }
  }

  ComputeSimpleCustomOutputs(config);
}

void CAdjElasticityOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Struc = solver[ADJFEA_SOL]->GetNodes();
  CPoint*    Node_Geo  = geometry->nodes;

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

  SetVolumeOutputValue("ADJOINT-X", iPoint, Node_Struc->GetSolution(iPoint, 0));
  SetVolumeOutputValue("ADJOINT-Y", iPoint, Node_Struc->GetSolution(iPoint, 1));
  if (nVar_FEM == 3)
    SetVolumeOutputValue("ADJOINT-Z", iPoint, Node_Struc->GetSolution(iPoint, 2));

  SetVolumeOutputValue("SENSITIVITY-X", iPoint, Node_Struc->GetSensitivity(iPoint, 0));
  SetVolumeOutputValue("SENSITIVITY-Y", iPoint, Node_Struc->GetSensitivity(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENSITIVITY-Z", iPoint, Node_Struc->GetSensitivity(iPoint, 2));

  if (!config->GetTime_Domain()) return;

  SetVolumeOutputValue("SENS_DISP-X", iPoint, Node_Struc->GetSolution_time_n(iPoint, 0));
  SetVolumeOutputValue("SENS_DISP-Y", iPoint, Node_Struc->GetSolution_time_n(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENS_DISP-Z", iPoint, Node_Struc->GetSolution_time_n(iPoint, 2));

  SetVolumeOutputValue("SENS_VEL-X", iPoint, Node_Struc->GetSolution_time_n(iPoint, nDim));
  SetVolumeOutputValue("SENS_VEL-Y", iPoint, Node_Struc->GetSolution_time_n(iPoint, nDim + 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENS_VEL-Z", iPoint, Node_Struc->GetSolution_time_n(iPoint, 5));

  SetVolumeOutputValue("SENS_ACCEL-X", iPoint, Node_Struc->GetSolution_time_n(iPoint, 2 * nDim));
  SetVolumeOutputValue("SENS_ACCEL-Y", iPoint, Node_Struc->GetSolution_time_n(iPoint, 2 * nDim + 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENS_ACCEL-Z", iPoint, Node_Struc->GetSolution_time_n(iPoint, 8));
}

void CAdjElasticityOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  /// BEGIN_GROUP: SOLUTION, DESCRIPTION: Adjoint variables of the current objective function.
  /// DESCRIPTION: Adjoint x-component.
  AddVolumeOutput("ADJOINT-X", "Adjoint_x", "SOLUTION", "adjoint of displacement in the x direction");
  /// DESCRIPTION: Adjoint y-component.
  AddVolumeOutput("ADJOINT-Y", "Adjoint_y", "SOLUTION", "adjoint of displacement in the y direction");
  if (nVar_FEM == 3)
    /// DESCRIPTION: Adjoint z-component.
    AddVolumeOutput("ADJOINT-Z", "Adjoint_z", "SOLUTION", "adjoint of displacement in the z direction");
  /// END_GROUP

  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Geometrical sensitivities of the current objective function.
  /// DESCRIPTION: Sensitivity x-component.
  AddVolumeOutput("SENSITIVITY-X", "Sensitivity_x", "SENSITIVITY", "geometric sensitivity in the x direction");
  /// DESCRIPTION: Sensitivity y-component.
  AddVolumeOutput("SENSITIVITY-Y", "Sensitivity_y", "SENSITIVITY", "geometric sensitivity  in the y direction");
  if (nDim == 3)
    /// DESCRIPTION: Sensitivity z-component.
    AddVolumeOutput("SENSITIVITY-Z", "Sensitivity_z", "SENSITIVITY", "geometric sensitivity  in the z direction");
  /// END_GROUP

  if (!config->GetTime_Domain()) return;

  /*--- Sensitivities with respect to initial conditions. ---*/

  AddVolumeOutput("SENS_DISP-X", "SensitivityDispN_x", "SENSITIVITY_N", "sensitivity to the previous x displacement");
  AddVolumeOutput("SENS_DISP-Y", "SensitivityDispN_y", "SENSITIVITY_N", "sensitivity to the previous y displacement");
  if (nDim == 3)
    AddVolumeOutput("SENS_DISP-Z", "SensitivityDispN_z", "SENSITIVITY_N", "sensitivity to the previous z displacement");

  AddVolumeOutput("SENS_VEL-X", "SensitivityVelN_x", "SENSITIVITY_N", "sensitivity to the previous x velocity");
  AddVolumeOutput("SENS_VEL-Y", "SensitivityVelN_y", "SENSITIVITY_N", "sensitivity to the previous y velocity");
  if (nDim == 3)
    AddVolumeOutput("SENS_VEL-Z", "SensitivityVelN_z", "SENSITIVITY_N", "sensitivity to the previous z velocity");

  AddVolumeOutput("SENS_ACCEL-X", "SensitivityAccelN_x", "SENSITIVITY_N", "sensitivity to the previous x acceleration");
  AddVolumeOutput("SENS_ACCEL-Y", "SensitivityAccelN_y", "SENSITIVITY_N", "sensitivity to the previous y acceleration");
  if (nDim == 3)
    AddVolumeOutput("SENS_ACCEL-Z", "SensitivityAccelN_z", "SENSITIVITY_N", "sensitivity to the previous z acceleration");

}
