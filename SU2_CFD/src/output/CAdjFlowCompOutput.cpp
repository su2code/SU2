/*!
 * \file output_adj_flow_comp.cpp
 * \brief Main subroutines for flow discrete adjoint output
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


#include "../../include/output/CAdjFlowOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CAdjFlowCompOutput::CAdjFlowCompOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

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

  if (convFields.empty() ) convFields.emplace_back("RMS_ADJ_DENSITY");

}

CAdjFlowCompOutput::~CAdjFlowCompOutput(void) {}

void CAdjFlowCompOutput::SetHistoryOutputFields(CConfig *config){

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
  if (!config->GetFrozen_Visc_Disc()){
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      /// DESCRIPTION: BGS residual of the adjoint nu tilde.
      AddHistoryOutput("BGS_ADJ_NU_TILDE", "bgs[A_nu]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);
      break;
    case SST:
      /// DESCRIPTION: BGS residual of the adjoint kinetic energy.
      AddHistoryOutput("BGS_ADJ_TKE", "bgs[A_k]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint kinetic energy.", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: BGS residual of the adjoint dissipation.
      AddHistoryOutput("BGS_ADJ_DISSIPATION",    "bgs[A_w]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint dissipation.", HistoryFieldType::RESIDUAL);
      break;
    default: break;
    }
  }

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

}

void CAdjFlowCompOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver){

  CSolver* adjflow_solver = solver[ADJFLOW_SOL];
  CSolver* adjturb_solver = solver[ADJTURB_SOL];

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
      SetHistoryOutputValue("MAX_ADJ_DISSIPATION",    log10(adjturb_solver->GetRes_Max(1)));
      break;
    default: break;
    }
  }

  if (multiZone){
    SetHistoryOutputValue("BGS_ADJ_DENSITY", log10(adjflow_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_BGS(1)));
    SetHistoryOutputValue("BGS_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_BGS(2)));
    if (geometry->GetnDim() == 3) {
      SetHistoryOutputValue("BGS_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_BGS(3)));
      SetHistoryOutputValue("BGS_ADJ_ENERGY", log10(adjflow_solver->GetRes_BGS(4)));
    } else {
      SetHistoryOutputValue("BGS_ADJ_ENERGY", log10(adjflow_solver->GetRes_BGS(3)));
    }
    if ((!config->GetFrozen_Visc_Disc() && !cont_adj) || (!config->GetFrozen_Visc_Cont() && cont_adj)){
      switch(turb_model){
      case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
        SetHistoryOutputValue("BGS_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_BGS(0)));
        break;
      case SST:
        SetHistoryOutputValue("BGS_ADJ_TKE", log10(adjturb_solver->GetRes_BGS(0)));
        SetHistoryOutputValue("BGS_ADJ_DISSIPATION",    log10(adjturb_solver->GetRes_BGS(1)));
        break;
      default: break;
      }
    }
  }

  SetHistoryOutputValue("SENS_GEO", adjflow_solver->GetTotal_Sens_Geo());
  SetHistoryOutputValue("SENS_AOA", adjflow_solver->GetTotal_Sens_AoA() * PI_NUMBER / 180.0);
  SetHistoryOutputValue("SENS_MACH", adjflow_solver->GetTotal_Sens_Mach());
  SetHistoryOutputValue("SENS_PRESS", adjflow_solver->GetTotal_Sens_Press());
  SetHistoryOutputValue("SENS_TEMP", adjflow_solver->GetTotal_Sens_Temp());

}

void CAdjFlowCompOutput::SetVolumeOutputFields(CConfig *config){

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

void CAdjFlowCompOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_AdjFlow = solver[ADJFLOW_SOL]->GetNodes();
  CVariable* Node_AdjTurb = NULL;
  CPoint*    Node_Geo     = geometry->node[iPoint];

  if (config->GetKind_Turb_Model() != NONE &&
      ((!config->GetFrozen_Visc_Disc() && !cont_adj) ||
       (!config->GetFrozen_Visc_Cont() && cont_adj))){
    Node_AdjTurb = solver[ADJTURB_SOL]->GetNodes();
  }

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));

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

void CAdjFlowCompOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){

  SetVolumeOutputValue("SENSITIVITY", iPoint, solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex));

}


bool CAdjFlowCompOutput::SetInit_Residuals(CConfig *config){

  return ((config->GetTime_Marching() != STEADY) && (curInnerIter == 0)) ||
         ((config->GetTime_Marching() == STEADY) && (curInnerIter < 2));

}

bool CAdjFlowCompOutput::SetUpdate_Averages(CConfig *config){
  return false;

//  return (config->GetUnsteady_Simulation() != STEADY && !dualtime);

}

