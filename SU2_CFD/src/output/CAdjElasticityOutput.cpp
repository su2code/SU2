/*!
 * \file CAdjElasticityOutput.cpp
 * \brief Main subroutines for elasticity discrete adjoint output
 * \author R. Sanchez
 * \version 7.0.4 "Blackbird"
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


#include "../../include/output/CAdjElasticityOutput.hpp"

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
    requestedScreenFields.emplace_back("SENS_E");
    requestedScreenFields.emplace_back("SENS_NU");
    nRequestedScreenFields = requestedScreenFields.size();
  }

  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("SENSITIVITY");
    nRequestedVolumeFields = requestedVolumeFields.size();
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

CAdjElasticityOutput::~CAdjElasticityOutput(void) {}

void CAdjElasticityOutput::SetHistoryOutputFields(CConfig *config){

  // Residuals
  AddHistoryOutput("ADJOINT_DISP_X", "rms[Ux_adj]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint of the X displacements.", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("ADJOINT_DISP_Y", "rms[Uy_adj]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint of the Y displacements.", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("ADJOINT_DISP_Z", "rms[Uz_adj]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint of the Z displacements.", HistoryFieldType::RESIDUAL);

  //Sensitivities
  AddHistoryOutput("SENS_E", "Sens[E]",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "");
  AddHistoryOutput("SENS_NU","Sens[Nu]", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "");

  AddHistoryOutput("COMBO", "ObjFun", ScreenOutputFormat::SCIENTIFIC, "COMBO", "", HistoryFieldType::COEFFICIENT);

}

inline void CAdjElasticityOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  SetHistoryOutputValue("ADJOINT_DISP_X", log10(solver[ADJFEA_SOL]->GetRes_RMS(0)));
  SetHistoryOutputValue("ADJOINT_DISP_Y", log10(solver[ADJFEA_SOL]->GetRes_RMS(1)));
  if (nVar_FEM == 3){
    SetHistoryOutputValue("ADJOINT_DISP_Z", log10(solver[ADJFEA_SOL]->GetRes_RMS(2)));
  }

  su2double Total_SensE = 0.0; su2double Total_SensNu = 0.0;
  if (config->GetnElasticityMod() == 1) {
    Total_SensE = solver[ADJFEA_SOL]->GetGlobal_Sens_E(0);
    Total_SensNu = solver[ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
  }
  else {
    for (unsigned short iVar = 0; iVar < config->GetnElasticityMod(); iVar++){
      Total_SensE += pow(solver[ADJFEA_SOL]->GetGlobal_Sens_E(iVar),2);
      Total_SensNu += pow(solver[ADJFEA_SOL]->GetGlobal_Sens_Nu(iVar),2);
    }
    Total_SensE = sqrt(Total_SensE);
    Total_SensNu = sqrt(Total_SensNu);
  }
  SetHistoryOutputValue("SENS_E", Total_SensE);
  SetHistoryOutputValue("SENS_NU", Total_SensNu);

  SetHistoryOutputValue("COMBO", solver[FEA_SOL]->GetTotal_ComboObj());

}

void CAdjElasticityOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Struc = solver[ADJFEA_SOL]->GetNodes();
  CPoint*    Node_Geo  = geometry->node[iPoint];

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));

  SetVolumeOutputValue("ADJOINT-X", iPoint, Node_Struc->GetSolution(iPoint, 0));
  SetVolumeOutputValue("ADJOINT-Y", iPoint, Node_Struc->GetSolution(iPoint, 1));
  if (nVar_FEM == 3)
    SetVolumeOutputValue("ADJOINT-Z", iPoint, Node_Struc->GetSolution(iPoint, 2));

  SetVolumeOutputValue("SENSITIVITY-X", iPoint, Node_Struc->GetSensitivity(iPoint, 0));
  SetVolumeOutputValue("SENSITIVITY-Y", iPoint, Node_Struc->GetSensitivity(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENSITIVITY-Z", iPoint, Node_Struc->GetSensitivity(iPoint, 2));
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

}
