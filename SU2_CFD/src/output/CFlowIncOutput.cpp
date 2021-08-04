/*!
 * \file CFlowIncOutput.cpp
 * \brief Main subroutines for incompressible flow output
 * \author R. Sanchez
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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
  scalar_model = config->GetKind_Scalar_Model();
  heat = config->GetEnergy_Equation();

  weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

  gridMovement = config->GetDynamic_Grid();
  streamwisePeriodic             = (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE);
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

  if (convFields.empty() ) convFields.emplace_back("RMS_PRESSURE");


}

CFlowIncOutput::~CFlowIncOutput(void) {}


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
  /// DESCRIPTION: Root-mean square residual of the radiative energy (P1 model).
  if (config->AddRadiation()) AddHistoryOutput("RMS_RAD_ENERGY", "rms[E_Rad]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the radiative energy.", HistoryFieldType::RESIDUAL);

  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).
    AddHistoryOutput("RMS_NU_TILDE",       "rms[nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
    break;
  case SST: case SST_SUST:
    /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).
    AddHistoryOutput("RMS_TKE", "rms[k]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).
    AddHistoryOutput("RMS_DISSIPATION",    "rms[w]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of dissipation (SST model).", HistoryFieldType::RESIDUAL);
    break;
  default: break;
  }
  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      AddHistoryOutput("RMS_PASSIVE_SCALAR", "rms[c]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean squared residual of the passive scalar equation.", HistoryFieldType::RESIDUAL);
      break;
    case PROGRESS_VARIABLE:
      AddHistoryOutput("RMS_PROGRESS_VARIABLE", "rms[PV]"       , ScreenOutputFormat::FIXED  , "RMS_RES", "Root-mean squared residual of the progress variable equation.", HistoryFieldType::RESIDUAL);
      AddHistoryOutput("RMS_ENTHALPY"         , "rms[Enth]"     , ScreenOutputFormat::FIXED  , "RMS_RES", "Root-mean squared residual of the enthalpy equation."         , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("RMS_Y_CO"             , "rms[Y_CO]"     , ScreenOutputFormat::FIXED  , "RMS_RES", "Root-mean squared residual of the CO mass fraction equation." , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("RMS_Y_NOX"            , "rms[Y_NOx]"    , ScreenOutputFormat::FIXED  , "RMS_RES", "Root-mean squared residual of the NOx mass fraction equation.", HistoryFieldType::RESIDUAL);
      AddHistoryOutput("N_TABLE_MISSES"       , "# table misses", ScreenOutputFormat::INTEGER, "RMS_RES", "number of table misses"                              , HistoryFieldType::RESIDUAL);
      break;
    default: break;
  }
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
  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      AddHistoryOutput("MAX_PASSIVE_SCALAR", "max[c]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the passive scalar equation.", HistoryFieldType::RESIDUAL);
      break;
    case PROGRESS_VARIABLE:
      AddHistoryOutput("MAX_PROGRESS_VARIABLE" , "max[PV]"    , ScreenOutputFormat::FIXED , "MAX_RES", "Maximum residual of the progress variable equation." , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("MAX_ENTHALPY"          , "max[Enth]"  , ScreenOutputFormat::FIXED , "MAX_RES", "Maximum residual of the enthalpy equation."          , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("MAX_Y_CO"              , "max[Y_CO]"  , ScreenOutputFormat::FIXED , "MAX_RES", "Maximum residual of the CO mass fraction equation."  , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("MAX_Y_NOX"             , "max[Y_NOx]" , ScreenOutputFormat::FIXED , "MAX_RES", "Maximum residual of the NOx mass fraction equation." , HistoryFieldType::RESIDUAL);
      break;
    default: break;
  }
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
  /// DESCRIPTION: Multizone residual of the radiative energy (P1 model).
  if (config->AddRadiation()) AddHistoryOutput("BGS_RAD_ENERGY", "bgs[E_Rad]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the radiative energy.", HistoryFieldType::RESIDUAL);

  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Maximum residual of nu tilde (SA model).
    AddHistoryOutput("BGS_NU_TILDE",       "bgs[nu]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
    break;
  case SST: case SST_SUST:
    /// DESCRIPTION: Maximum residual of kinetic energy (SST model).
    AddHistoryOutput("BGS_TKE", "bgs[k]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Maximum residual of the dissipation (SST model).
    AddHistoryOutput("BGS_DISSIPATION",    "bgs[w]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of dissipation (SST model).", HistoryFieldType::RESIDUAL);
    break;
  default: break;
  }
  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      AddHistoryOutput("BGS_PASSIVE_SCALAR", "bgs[c]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the passive scalar equation.", HistoryFieldType::RESIDUAL);
      AddHistoryOutput("LINSOL_ITER_SCALAR", "LinSolIter[c]", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear scalar solver.");
      AddHistoryOutput("LINSOL_RESIDUAL_SCALAR", "LinSolRes[c]", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear scalar solver.");
      break;
    case PROGRESS_VARIABLE:
      AddHistoryOutput("BGS_PROGRESS_VARIABLE" , "bgs[PV]"    , ScreenOutputFormat::FIXED , "BGS_RES", "BGS residual of the progress variable equation." , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("BGS_ENTHALPY"          , "bgs[Enth]"  , ScreenOutputFormat::FIXED , "BGS_RES", "BGS residual of the enthalpy equation."          , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("BGS_Y_CO"              , "bgs[Y_CO]"  , ScreenOutputFormat::FIXED , "BGS_RES", "BGS residual of the CO mass fraction equation."  , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("BGS_Y_NOX"             , "bgs[Y_NOx]" , ScreenOutputFormat::FIXED , "BGS_RES", "BGS residual of the NOx mass fraction equation." , HistoryFieldType::RESIDUAL);
      AddHistoryOutput("LINSOL_ITER_SCALAR",     "LinSolIter[c]", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear scalar solver.");
      AddHistoryOutput("LINSOL_RESIDUAL_SCALAR", "LinSolRes[c]", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear scalar solver.");
      break;
    default: break;
  }
  /// END_GROUP

  /// BEGIN_GROUP: ROTATING_FRAME, DESCRIPTION: Coefficients related to a rotating frame of reference.
  /// DESCRIPTION: Merit
  AddHistoryOutput("FIGURE_OF_MERIT", "CMerit", ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "Merit", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: CT
  AddHistoryOutput("THRUST",    "CT",     ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "CT", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: CQ
  AddHistoryOutput("TORQUE",    "CQ",     ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "CQ", HistoryFieldType::COEFFICIENT);
  /// END_GROUP

  /// BEGIN_GROUP: HEAT_COEFF, DESCRIPTION: Heat coefficients on all surfaces set with MARKER_MONITORING.
  /// DESCRIPTION: Total heatflux
  AddHistoryOutput("TOTAL_HEATFLUX", "HF",      ScreenOutputFormat::SCIENTIFIC, "HEAT", "Total heatflux on all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Maximal heatflux
  AddHistoryOutput("MAXIMUM_HEATFLUX", "maxHF", ScreenOutputFormat::SCIENTIFIC, "HEAT", "Total maximum heatflux on all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Temperature
  if (heat || weakly_coupled_heat)
    AddHistoryOutput("AVG_TEMPERATURE", "Temp", ScreenOutputFormat::SCIENTIFIC, "HEAT",  "Total avg. temperature on all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// END_GROUP

  /// DESCRIPTION: Angle of attack
  AddHistoryOutput("AOA",         "AoA",                      ScreenOutputFormat::SCIENTIFIC,"AOA", "Angle of attack");
  /// DESCRIPTION: Linear solver iterations
  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

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
  /*--- Add analyze surface history fields --- */

  AddAnalyzeSurfaceOutput(config);

  /*--- Add aerodynamic coefficients fields --- */

  AddAerodynamicCoefficients(config);

}

void CFlowIncOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  CSolver* flow_solver = solver[FLOW_SOL];
  CSolver* turb_solver = solver[TURB_SOL];
  CSolver* heat_solver = solver[HEAT_SOL];
  CSolver* rad_solver  = solver[RAD_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];
  CSolver* scalar_solver = solver[SCALAR_SOL];

  SetHistoryOutputValue("RMS_PRESSURE", log10(flow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_VELOCITY-X", log10(flow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_VELOCITY-Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 3) SetHistoryOutputValue("RMS_VELOCITY-Z", log10(flow_solver->GetRes_RMS(3)));

  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("RMS_NU_TILDE", log10(turb_solver->GetRes_RMS(0)));
    break;
  case SST: case SST_SUST:
    SetHistoryOutputValue("RMS_TKE", log10(turb_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_DISSIPATION",    log10(turb_solver->GetRes_RMS(1)));
    break;
  }

  if (config->AddRadiation())
    SetHistoryOutputValue("RMS_RAD_ENERGY", log10(rad_solver->GetRes_RMS(0)));


  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      SetHistoryOutputValue("RMS_PASSIVE_SCALAR", log10(scalar_solver->GetRes_RMS(0)));
      SetHistoryOutputValue("LINSOL_ITER_SCALAR", scalar_solver->GetIterLinSolver());
      SetHistoryOutputValue("LINSOL_RESIDUAL_SCALAR", log10(scalar_solver->GetResLinSolver()));      break;
    case PROGRESS_VARIABLE:
      SetHistoryOutputValue("RMS_PROGRESS_VARIABLE", log10(scalar_solver->GetRes_RMS(I_PROG_VAR)));
      SetHistoryOutputValue("RMS_ENTHALPY"         , log10(scalar_solver->GetRes_RMS(I_ENTHALPY)));
      SetHistoryOutputValue("RMS_Y_CO"             , log10(scalar_solver->GetRes_RMS(I_CO)      ));
      SetHistoryOutputValue("RMS_Y_NOX"            , log10(scalar_solver->GetRes_RMS(I_NOX)     ));
      SetHistoryOutputValue("N_TABLE_MISSES"       , scalar_solver->GetNTableMisses());
      SetHistoryOutputValue("LINSOL_ITER_SCALAR", scalar_solver->GetIterLinSolver());
      SetHistoryOutputValue("LINSOL_RESIDUAL_SCALAR", log10(scalar_solver->GetResLinSolver()));      break;
      break;
  }
  
  SetHistoryOutputValue("MAX_PRESSURE", log10(flow_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_VELOCITY-X", log10(flow_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_VELOCITY-Y", log10(flow_solver->GetRes_Max(2)));
  if (nDim == 3) SetHistoryOutputValue("RMS_VELOCITY-Z", log10(flow_solver->GetRes_Max(3)));

  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("MAX_NU_TILDE", log10(turb_solver->GetRes_Max(0)));
    break;
  case SST: case SST_SUST:
    SetHistoryOutputValue("MAX_TKE", log10(turb_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_DISSIPATION",    log10(turb_solver->GetRes_Max(1)));
    break;
  }
  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      SetHistoryOutputValue("MAX_PASSIVE_SCALAR", log10(scalar_solver->GetRes_Max(0)));
      break;
    case PROGRESS_VARIABLE:
      SetHistoryOutputValue("MAX_PROGRESS_VARIABLE", log10(scalar_solver->GetRes_Max(I_PROG_VAR)));
      SetHistoryOutputValue("MAX_ENTHALPY"         , log10(scalar_solver->GetRes_Max(I_ENTHALPY)));
      SetHistoryOutputValue("MAX_Y_CO"             , log10(scalar_solver->GetRes_Max(I_CO)      ));
      SetHistoryOutputValue("MAX_Y_NOX"            , log10(scalar_solver->GetRes_Max(I_NOX)     ));
      break;
  }
  
  if (multiZone){
    SetHistoryOutputValue("BGS_PRESSURE", log10(flow_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_VELOCITY-X", log10(flow_solver->GetRes_BGS(1)));
    SetHistoryOutputValue("BGS_VELOCITY-Y", log10(flow_solver->GetRes_BGS(2)));
    if (nDim == 3) SetHistoryOutputValue("BGS_VELOCITY-Z", log10(flow_solver->GetRes_BGS(3)));

    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      SetHistoryOutputValue("BGS_NU_TILDE", log10(turb_solver->GetRes_BGS(0)));
      break;
    case SST:
      SetHistoryOutputValue("BGS_TKE", log10(turb_solver->GetRes_BGS(0)));
      SetHistoryOutputValue("BGS_DISSIPATION",    log10(turb_solver->GetRes_BGS(1)));
      break;
    }

    if (config->AddRadiation())
      SetHistoryOutputValue("BGS_RAD_ENERGY", log10(rad_solver->GetRes_BGS(0)));

    
    switch(scalar_model){
      case PASSIVE_SCALAR:
        SetHistoryOutputValue("BGS_PASSIVE_SCALAR", log10(scalar_solver->GetRes_BGS(0)));
        break;
      case PROGRESS_VARIABLE:
        SetHistoryOutputValue("BGS_PROGRESS_VARIABLE", log10(scalar_solver->GetRes_BGS(I_PROG_VAR)));
        SetHistoryOutputValue("BGS_ENTHALPY"         , log10(scalar_solver->GetRes_BGS(I_ENTHALPY)));
        SetHistoryOutputValue("BGS_Y_CO"             , log10(scalar_solver->GetRes_BGS(I_CO)      ));
        SetHistoryOutputValue("BGS_Y_NOX"            , log10(scalar_solver->GetRes_BGS(I_NOX)     )); 
        break;
    }
  }

  if (weakly_coupled_heat){
    SetHistoryOutputValue("TOTAL_HEATFLUX",   heat_solver->GetTotal_HeatFlux());
    SetHistoryOutputValue("MAXIMUM_HEATFLUX", heat_solver->GetTotal_MaxHeatFlux());
    SetHistoryOutputValue("AVG_TEMPERATURE",  heat_solver->GetTotal_AvgTemperature());
    SetHistoryOutputValue("RMS_TEMPERATURE",  log10(heat_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("MAX_TEMPERATURE",  log10(heat_solver->GetRes_Max(0)));
    if (multiZone) SetHistoryOutputValue("BGS_TEMPERATURE", log10(heat_solver->GetRes_BGS(0)));
  }
  if (heat){
    SetHistoryOutputValue("TOTAL_HEATFLUX",   flow_solver->GetTotal_HeatFlux());
    SetHistoryOutputValue("MAXIMUM_HEATFLUX", flow_solver->GetTotal_MaxHeatFlux());
    SetHistoryOutputValue("AVG_TEMPERATURE",  flow_solver->GetTotal_AvgTemperature());
    if (nDim == 3) SetHistoryOutputValue("RMS_TEMPERATURE", log10(flow_solver->GetRes_RMS(4)));
    else           SetHistoryOutputValue("RMS_TEMPERATURE", log10(flow_solver->GetRes_RMS(3)));

    if (nDim == 3) SetHistoryOutputValue("MAX_TEMPERATURE", log10(flow_solver->GetRes_Max(4)));
    else           SetHistoryOutputValue("MAX_TEMPERATURE", log10(flow_solver->GetRes_Max(3)));
    if (multiZone){
      if (nDim == 3) SetHistoryOutputValue("BGS_TEMPERATURE", log10(flow_solver->GetRes_BGS(4)));
      else           SetHistoryOutputValue("BGS_TEMPERATURE", log10(flow_solver->GetRes_BGS(3)));
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

  SetRotatingFrameCoefficients(config, flow_solver);

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
  if (heat || weakly_coupled_heat)
    AddVolumeOutput("TEMPERATURE",  "Temperature","SOLUTION", "Temperature");

  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    AddVolumeOutput("TKE", "Turb_Kin_Energy", "SOLUTION", "Turbulent kinetic energy");
    AddVolumeOutput("DISSIPATION", "Omega", "SOLUTION", "Rate of dissipation");
    break;
  case SA: case SA_COMP: case SA_E:
  case SA_E_COMP: case SA_NEG:
    AddVolumeOutput("NU_TILDE", "Nu_Tilde", "SOLUTION", "Spalart–Allmaras variable");
    break;
  case NONE:
    break;
  }
  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      AddVolumeOutput("PASSIVE_SCALAR", "Passive_Scalar", "SOLUTION", "Passive scalar solution");
      AddVolumeOutput("DIFFUSIVITY"   , "Diffusivity"   , "SOLUTION", "Passive scalar diffusivity");
      AddVolumeOutput("SPECIFIC_HEAT_CP"   , "Specific_Heat_Cp"   , "SOLUTION", "Mixture specific heat cp");
      AddVolumeOutput("CONDUCTIVITY"   , "Conductivity"   , "SOLUTION", "Mixture conductivity");
      AddVolumeOutput("MEAN_MOLECULAR_WEIGHT"   , "Mean_Molecular_Weight"   , "SOLUTION", "Mixture molecular weight");
      break;
    case PROGRESS_VARIABLE:
      AddVolumeOutput("PROGRESS_VARIABLE", "Progress_Variable", "SOLUTION", "Progress variable solution");
      AddVolumeOutput("ENTHALPY"         , "Enthalpy"         , "SOLUTION", "Enthalpy solution"         );
      AddVolumeOutput("Y_CO"             , "Y_CO"             , "SOLUTION", "CO Mass fraction solution" );
      AddVolumeOutput("Y_NOX"            , "Y_NOx"            , "SOLUTION", "NOx Mass fraction solution");
      AddVolumeOutput("DIFFUSIVITY_PV"   , "Diffusivity_PV"   , "SOLUTION", "Diffusivity of the progress variable");
      AddVolumeOutput("DIFFUSIVITY_ENTH" , "Diffusivity_Enth" , "SOLUTION", "Diffusivity of the enthalpy");
      AddVolumeOutput("SPECIFIC_HEAT_CP"   , "Specific_Heat_Cp"   , "SOLUTION", "Mixture specific heat cp");

      for (int i_lookup = 0; i_lookup < config->GetNLookups(); ++i_lookup)
        if (config->GetLookupName(i_lookup)!="NULL"){ 
          string strname1="lookup_"+config->GetLookupName(i_lookup);
          AddVolumeOutput(config->GetLookupName(i_lookup),strname1,"LOOKUP",config->GetLookupName(i_lookup));
        }
      AddVolumeOutput("TABLE_MISSES"       , "Table_misses"       , "SOLUTION", "Lookup table misses");
      

      break;
    case NO_SCALAR_MODEL:
      break;
  }

  // Sources
  switch (scalar_model) {
    case PASSIVE_SCALAR:
      break;
    case PROGRESS_VARIABLE:
      AddVolumeOutput("SOURCE_PROGRESS_VARIABLE", "Source_Progress_Variable", "SOURCE", "Source Progress Variable");
      AddVolumeOutput("SOURCE_ENTHALPY"         , "Source_Enthalpy"         , "SOURCE", "Source Enthalpy"         );
      AddVolumeOutput("SOURCE_Y_CO"             , "Source_Y_CO"             , "SOURCE", "Source Y_CO"             );
      AddVolumeOutput("SOURCE_Y_NOX"            , "Source_Y_NOx"            , "SOURCE", "Source Y_NOx"            );
    case NO_SCALAR_MODEL:
      break;
  }

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

  if (config->GetKind_Solver() == INC_RANS || config->GetKind_Solver() == INC_NAVIER_STOKES){
    AddVolumeOutput("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity");

    AddVolumeOutput("SKIN_FRICTION-X", "Skin_Friction_Coefficient_x", "PRIMITIVE", "x-component of the skin friction vector");
    AddVolumeOutput("SKIN_FRICTION-Y", "Skin_Friction_Coefficient_y", "PRIMITIVE", "y-component of the skin friction vector");
    if (nDim == 3)
      AddVolumeOutput("SKIN_FRICTION-Z", "Skin_Friction_Coefficient_z", "PRIMITIVE", "z-component of the skin friction vector");

    AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heat-flux");
    AddVolumeOutput("Y_PLUS", "Y_Plus", "PRIMITIVE", "Non-dim. wall distance (Y-Plus)");

  }

  if (config->GetKind_Solver() == INC_RANS) {
    AddVolumeOutput("EDDY_VISCOSITY", "Eddy_Viscosity", "PRIMITIVE", "Turbulent eddy viscosity");
  }

  if (config->GetKind_Trans_Model() == BC){
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

  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    AddVolumeOutput("RES_TKE", "Residual_TKE", "RESIDUAL", "Residual of turbulent kinetic energy");
    AddVolumeOutput("RES_DISSIPATION", "Residual_Omega", "RESIDUAL", "Residual of the rate of dissipation.");
    break;
  case SA: case SA_COMP: case SA_E:
  case SA_E_COMP: case SA_NEG:
    AddVolumeOutput("RES_NU_TILDE", "Residual_Nu_Tilde", "RESIDUAL", "Residual of the Spalart–Allmaras variable");
    break;
  case NONE:
    break;
  }
  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      AddVolumeOutput("RES_PASSIVE_SCALAR", "Residual_Passive_Scalar", "RESIDUAL", "Residual of passive scalar equation");
      break;
    case PROGRESS_VARIABLE:
      AddVolumeOutput("RES_PROGRESS_VARIABLE", "Residual_Progress_Variable", "RESIDUAL", "Residual of the Progress Variable equation");
      AddVolumeOutput("RES_ENTHALPY"         , "Residual_Enthalpy"         , "RESIDUAL", "Residual of the Enthalpy equation"         );
      AddVolumeOutput("RES_Y_CO"             , "Residual_Y_CO"             , "RESIDUAL", "Residual of the Y_CO equation"             );
      AddVolumeOutput("RES_Y_NOX"            , "Residual_Y_NOx"            , "RESIDUAL", "Residual of the Y_NOx equation"            );
    break;

    case NO_SCALAR_MODEL:
      break;
  }
  
  if (config->GetKind_SlopeLimit_Flow() != NO_LIMITER && config->GetKind_SlopeLimit_Flow() != VAN_ALBADA_EDGE) {
    AddVolumeOutput("LIMITER_PRESSURE", "Limiter_Pressure", "LIMITER", "Limiter value of the pressure");
    AddVolumeOutput("LIMITER_VELOCITY-X", "Limiter_Velocity_x", "LIMITER", "Limiter value of the x-velocity");
    AddVolumeOutput("LIMITER_VELOCITY-Y", "Limiter_Velocity_y", "LIMITER", "Limiter value of the y-velocity");
    if (nDim == 3)
      AddVolumeOutput("LIMITER_VELOCITY-Z", "Limiter_Velocity_z", "LIMITER", "Limiter value of the z-velocity");
    if (heat || weakly_coupled_heat)
      AddVolumeOutput("LIMITER_TEMPERATURE", "Limiter_Temperature", "LIMITER", "Limiter value of the temperature");
  }

  if (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) {
    switch(config->GetKind_Turb_Model()){
    case SST: case SST_SUST:
      AddVolumeOutput("LIMITER_TKE", "Limiter_TKE", "LIMITER", "Limiter value of turb. kinetic energy.");
      AddVolumeOutput("LIMITER_DISSIPATION", "Limiter_Omega", "LIMITER", "Limiter value of dissipation rate.");
      break;
    case SA: case SA_COMP: case SA_E:
    case SA_E_COMP: case SA_NEG:
      AddVolumeOutput("LIMITER_NU_TILDE", "Limiter_Nu_Tilde", "LIMITER", "Limiter value of Spalart–Allmaras variable.");
      break;
    case NONE:
      break;
    }
  }
  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      AddVolumeOutput("LIMITER_PASSIVE_SCALAR", "Limiter_Passive_Scalar", "LIMITER", "Limiter value for the passive scalar");
      break;
    case PROGRESS_VARIABLE:
      AddVolumeOutput("LIMITER_PROGRESS_VARIABLE", "Limiter_Progress_Variable", "LIMITER", "Limiter value for the Progress Variable equation");
      AddVolumeOutput("LIMITER_ENTHALPY"         , "Limiter_Enthalpy"         , "LIMITER", "Limiter value for the Enthalpy equation"         );
      AddVolumeOutput("LIMITER_Y_CO"             , "Limiter_Y_CO"             , "LIMITER", "Limiter value for the Y_CO equation"             );
      AddVolumeOutput("LIMITER_Y_NOX"            , "Limiter_Y_NOx"            , "LIMITER", "Limiter value for the Y_NOx equation"            );
      break;
    case NO_SCALAR_MODEL:
      break;
  }
  
  // Hybrid RANS-LES
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    AddVolumeOutput("DES_LENGTHSCALE", "DES_LengthScale", "DDES", "DES length scale value");
    AddVolumeOutput("WALL_DISTANCE", "Wall_Distance", "DDES", "Wall distance value");
  }

  // Roe Low Dissipation
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    AddVolumeOutput("ROE_DISSIPATION", "Roe_Dissipation", "ROE_DISSIPATION", "Value of the Roe dissipation");
  }

  if(config->GetKind_Solver() == INC_RANS || config->GetKind_Solver() == INC_NAVIER_STOKES){
    if (nDim == 3){
      AddVolumeOutput("VORTICITY_X", "Vorticity_x", "VORTEX_IDENTIFICATION", "x-component of the vorticity vector");
      AddVolumeOutput("VORTICITY_Y", "Vorticity_y", "VORTEX_IDENTIFICATION", "y-component of the vorticity vector");
      AddVolumeOutput("VORTICITY_Z", "Vorticity_z", "VORTEX_IDENTIFICATION", "z-component of the vorticity vector");
    } else {
      AddVolumeOutput("VORTICITY", "Vorticity", "VORTEX_IDENTIFICATION", "Value of the vorticity");
    }
    AddVolumeOutput("Q_CRITERION", "Q_Criterion", "VORTEX_IDENTIFICATION", "Value of the Q-Criterion");
  }

  if(config->GetKind_TimeIntScheme_Flow()==EULER_IMPLICIT){
    AddVolumeOutput("TIMESTEP", "timestep", "TIMESTEP", "local timestep");
  }

  // Mesh quality metrics, computed in CPhysicalGeometry::ComputeMeshQualityStatistics.
  AddVolumeOutput("ORTHOGONALITY", "Orthogonality", "MESH_QUALITY", "Orthogonality Angle (deg.)");
  AddVolumeOutput("ASPECT_RATIO",  "Aspect_Ratio",  "MESH_QUALITY", "CV Face Area Aspect Ratio");
  AddVolumeOutput("VOLUME_RATIO",  "Volume_Ratio",  "MESH_QUALITY", "CV Sub-Volume Ratio");

  // Streamwise Periodicity
  if(streamwisePeriodic) {
    AddVolumeOutput("RECOVERED_PRESSURE", "Recovered_Pressure", "SOLUTION", "Recovered physical pressure");
    if (heat && streamwisePeriodic_temperature)
      AddVolumeOutput("RECOVERED_TEMPERATURE", "Recovered_Temperature", "SOLUTION", "Recovered physical temperature");
  }

  AddCommonFVMOutputs(config);
}

void CFlowIncOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Flow = solver[FLOW_SOL]->GetNodes();
  CVariable* Node_Heat = nullptr;
  CVariable* Node_Turb = nullptr;
  CVariable* Node_Rad = nullptr;
  const auto Node_Geo = geometry->nodes;
  CVariable* Node_Scalar = nullptr;

  if (config->GetKind_Turb_Model() != NONE){
    Node_Turb = solver[TURB_SOL]->GetNodes();
  }
  if (weakly_coupled_heat){
    Node_Heat = solver[HEAT_SOL]->GetNodes();
  }
  if (config->GetKind_Scalar_Model() != NONE){
    Node_Scalar = solver[SCALAR_SOL]->GetNodes();
  }

  LoadCoordinates(Node_Geo->GetCoord(iPoint), iPoint);

  SetVolumeOutputValue("PRESSURE",   iPoint, Node_Flow->GetSolution(iPoint, 0));
  SetVolumeOutputValue("VELOCITY-X", iPoint, Node_Flow->GetSolution(iPoint, 1));
  SetVolumeOutputValue("VELOCITY-Y", iPoint, Node_Flow->GetSolution(iPoint, 2));
  if (nDim == 3)
    SetVolumeOutputValue("VELOCITY-Z", iPoint, Node_Flow->GetSolution(iPoint, 3));

  if (heat) SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Flow->GetSolution(iPoint, nDim+1));
  if (weakly_coupled_heat) SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Heat->GetSolution(iPoint, 0));

  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    SetVolumeOutputValue("TKE", iPoint, Node_Turb->GetSolution(iPoint, 0));
    SetVolumeOutputValue("DISSIPATION", iPoint, Node_Turb->GetSolution(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E:
  case SA_E_COMP: case SA_NEG:
    SetVolumeOutputValue("NU_TILDE", iPoint, Node_Turb->GetSolution(iPoint, 0));
    break;
  case NONE:
    break;
  }
  
  // Solution data
  su2double *scalars;
  unsigned long table_misses;
  switch(scalar_model){
    case PASSIVE_SCALAR:
      SetVolumeOutputValue("PASSIVE_SCALAR", iPoint, Node_Scalar->GetSolution(iPoint, 0));
      SetVolumeOutputValue("DIFFUSIVITY"   , iPoint, Node_Scalar->GetDiffusivity(iPoint, 0));
      // SetVolumeOutputValue("SPECIFIC_HEAT_CP"   , iPoint, Node_Flow->GetSpecificHeatCp(iPoint));
      SetVolumeOutputValue("SPECIFIC_HEAT_CP"   , iPoint, 2224.43 * Node_Scalar->GetSolution(iPoint)[0] + 1009.39 * (1- Node_Scalar->GetSolution(iPoint)[0]));
      SetVolumeOutputValue("CONDUCTIVITY"   , iPoint, Node_Flow->GetThermalConductivity(iPoint));
      // SetVolumeOutputValue("MEAN_MOLECULAR_WEIGHT"   , iPoint, solver[FLOW_SOL]->GetFluidModel()->GetMeanMolecularWeight()); 
      SetVolumeOutputValue("MEAN_MOLECULAR_WEIGHT"   , iPoint, 1/(Node_Scalar->GetSolution(iPoint)[0]/(16.043/1000) + (1-Node_Scalar->GetSolution(iPoint)[0])/(28.965/1000)));
      break;
    case PROGRESS_VARIABLE:
      SetVolumeOutputValue("PROGRESS_VARIABLE", iPoint, Node_Scalar->GetSolution(iPoint, I_PROG_VAR));
      SetVolumeOutputValue("ENTHALPY"         , iPoint, Node_Scalar->GetSolution(iPoint, I_ENTHALPY));
      SetVolumeOutputValue("Y_CO"             , iPoint, Node_Scalar->GetSolution(iPoint, I_CO      ));
      SetVolumeOutputValue("Y_NOX"            , iPoint, Node_Scalar->GetSolution(iPoint, I_NOX     ));
      SetVolumeOutputValue("DIFFUSIVITY_PV"   , iPoint, Node_Scalar->GetDiffusivity(iPoint, I_PROG_VAR));
      SetVolumeOutputValue("DIFFUSIVITY_ENTH" , iPoint, Node_Scalar->GetDiffusivity(iPoint, I_ENTHALPY));
      SetVolumeOutputValue("SPECIFIC_HEAT_CP"   , iPoint, Node_Flow->GetSpecificHeatCp(iPoint));

      // update lookup
      scalars      = Node_Scalar->GetSolution(iPoint);
      table_misses = solver[FLOW_SOL]->GetFluidModel()->SetScalarLookups(scalars);
      for (int i_lookup = 0; i_lookup < config->GetNLookups(); ++i_lookup){
        if (config->GetLookupName(i_lookup)!="NULL")
          SetVolumeOutputValue(config->GetLookupName(i_lookup), iPoint, solver[FLOW_SOL]->GetFluidModel()->GetScalarLookups(i_lookup));
      }
      SetVolumeOutputValue("TABLE_MISSES"       , iPoint, (su2double)table_misses);

      break;
    case NO_SCALAR_MODEL:
      break;
  }

  // Sources
  switch(scalar_model){
    case PASSIVE_SCALAR:
      break;
    case PROGRESS_VARIABLE:
      scalars      = Node_Scalar->GetSolution(iPoint);
      table_misses = solver[FLOW_SOL]->GetFluidModel()->SetScalarSources(scalars);
      SetVolumeOutputValue("SOURCE_PROGRESS_VARIABLE", iPoint, solver[FLOW_SOL]->GetFluidModel()->GetScalarSources(I_PROG_VAR));
      SetVolumeOutputValue("SOURCE_ENTHALPY"         , iPoint, solver[FLOW_SOL]->GetFluidModel()->GetScalarSources(I_ENTHALPY));
      SetVolumeOutputValue("SOURCE_Y_CO"             , iPoint, solver[FLOW_SOL]->GetFluidModel()->GetScalarSources(I_CO      ));
      SetVolumeOutputValue("SOURCE_Y_NOX"            , iPoint, solver[FLOW_SOL]->GetFluidModel()->GetScalarSources(I_NOX     ));
      break;
    case NO_SCALAR_MODEL:
      break;
  }

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

  su2double VelMag = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    VelMag += pow(solver[FLOW_SOL]->GetVelocity_Inf(iDim),2.0);
  }
  su2double factor = 1.0/(0.5*solver[FLOW_SOL]->GetDensity_Inf()*VelMag);
  SetVolumeOutputValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure(iPoint) - config->GetPressure_FreeStreamND())*factor);
  SetVolumeOutputValue("DENSITY", iPoint, Node_Flow->GetDensity(iPoint));

  if (config->GetKind_Solver() == INC_RANS || config->GetKind_Solver() == INC_NAVIER_STOKES){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity(iPoint));
  }

  if (config->GetKind_Solver() == INC_RANS) {
    SetVolumeOutputValue("EDDY_VISCOSITY", iPoint, Node_Flow->GetEddyViscosity(iPoint));
  }

  if (config->GetKind_Trans_Model() == BC){
    SetVolumeOutputValue("INTERMITTENCY", iPoint, Node_Turb->GetGammaBC(iPoint));
  }

  SetVolumeOutputValue("RES_PRESSURE", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 0));
  SetVolumeOutputValue("RES_VELOCITY-X", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 1));
  SetVolumeOutputValue("RES_VELOCITY-Y", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 2));
  if (nDim == 3)
    SetVolumeOutputValue("RES_VELOCITY-Z", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 3));
  if (config->GetEnergy_Equation())
    SetVolumeOutputValue("RES_TEMPERATURE", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nDim+1));

  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    SetVolumeOutputValue("RES_TKE", iPoint, solver[TURB_SOL]->LinSysRes(iPoint, 0));
    SetVolumeOutputValue("RES_DISSIPATION", iPoint, solver[TURB_SOL]->LinSysRes(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E:
  case SA_E_COMP: case SA_NEG:
    SetVolumeOutputValue("RES_NU_TILDE", iPoint, solver[TURB_SOL]->LinSysRes(iPoint, 0));
    break;
  case NONE:
    break;
  }
  switch(scalar_model){
    case PASSIVE_SCALAR:
      SetVolumeOutputValue("RES_PASSIVE_SCALAR", iPoint, solver[SCALAR_SOL]->LinSysRes(iPoint, 0));
      break;
    case PROGRESS_VARIABLE:
      SetVolumeOutputValue("RES_PROGRESS_VARIABLE", iPoint, solver[SCALAR_SOL]->LinSysRes(iPoint, I_PROG_VAR));
      SetVolumeOutputValue("RES_ENTHALPY"         , iPoint, solver[SCALAR_SOL]->LinSysRes(iPoint, I_ENTHALPY));
      SetVolumeOutputValue("RES_Y_CO"             , iPoint, solver[SCALAR_SOL]->LinSysRes(iPoint, I_CO      ));
      SetVolumeOutputValue("RES_Y_NOX"            , iPoint, solver[SCALAR_SOL]->LinSysRes(iPoint, I_NOX     ));
     break;      
    case NO_SCALAR_MODEL:
      break;
  }
  

  if (config->GetKind_SlopeLimit_Flow() != NO_LIMITER && config->GetKind_SlopeLimit_Flow() != VAN_ALBADA_EDGE) {
    SetVolumeOutputValue("LIMITER_PRESSURE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 0));
    SetVolumeOutputValue("LIMITER_VELOCITY-X", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 1));
    SetVolumeOutputValue("LIMITER_VELOCITY-Y", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 2));
    if (nDim == 3)
      SetVolumeOutputValue("LIMITER_VELOCITY-Z", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 3));
    if (heat || weakly_coupled_heat)
      SetVolumeOutputValue("LIMITER_TEMPERATURE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nDim+1));

  }

  if (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) {
    switch(config->GetKind_Turb_Model()){
    case SST: case SST_SUST:
      SetVolumeOutputValue("LIMITER_TKE", iPoint, Node_Turb->GetLimiter(iPoint, 0));
      SetVolumeOutputValue("LIMITER_DISSIPATION", iPoint, Node_Turb->GetLimiter(iPoint, 1));
      break;
    case SA: case SA_COMP: case SA_E:
    case SA_E_COMP: case SA_NEG:
      SetVolumeOutputValue("LIMITER_NU_TILDE", iPoint, Node_Turb->GetLimiter(iPoint, 0));
      break;
    case NONE:
      break;
    }
  }
  
  switch(scalar_model){
    case PASSIVE_SCALAR:
      SetVolumeOutputValue("LIMITER_PASSIVE_SCALAR", iPoint, Node_Scalar->GetLimiter(iPoint, 0));
      break;
    case PROGRESS_VARIABLE:
      SetVolumeOutputValue("LIMITER_PROGRESS_VARIABLE", iPoint, Node_Scalar->GetLimiter(iPoint, I_PROG_VAR));
      SetVolumeOutputValue("LIMITER_ENTHALPY"         , iPoint, Node_Scalar->GetLimiter(iPoint, I_ENTHALPY));
      SetVolumeOutputValue("LIMITER_Y_CO"             , iPoint, Node_Scalar->GetLimiter(iPoint, I_CO      ));
      SetVolumeOutputValue("LIMITER_Y_NOX"            , iPoint, Node_Scalar->GetLimiter(iPoint, I_NOX     ));
      break;          
    case NO_SCALAR_MODEL:
      break;
  }
  
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    SetVolumeOutputValue("DES_LENGTHSCALE", iPoint, Node_Flow->GetDES_LengthScale(iPoint));
    SetVolumeOutputValue("WALL_DISTANCE", iPoint, Node_Geo->GetWall_Distance(iPoint));
  }

  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    SetVolumeOutputValue("ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation(iPoint));
  }

  if(config->GetKind_Solver() == INC_RANS || config->GetKind_Solver() == INC_NAVIER_STOKES){
    if (nDim == 3){
      SetVolumeOutputValue("VORTICITY_X", iPoint, Node_Flow->GetVorticity(iPoint)[0]);
      SetVolumeOutputValue("VORTICITY_Y", iPoint, Node_Flow->GetVorticity(iPoint)[1]);
      SetVolumeOutputValue("VORTICITY_Z", iPoint, Node_Flow->GetVorticity(iPoint)[2]);
    } else {
      SetVolumeOutputValue("VORTICITY", iPoint, Node_Flow->GetVorticity(iPoint)[2]);
    }
    SetVolumeOutputValue("Q_CRITERION", iPoint, GetQ_Criterion(Node_Flow->GetGradient_Primitive(iPoint,1)));
  }

  if(config->GetKind_TimeIntScheme_Flow()==EULER_IMPLICIT){
    SetVolumeOutputValue("TIMESTEP", iPoint, Node_Flow->GetDelta_Time(iPoint));
  }
    
  // Streamwise Periodicity
  if(streamwisePeriodic) {
    SetVolumeOutputValue("RECOVERED_PRESSURE", iPoint, Node_Flow->GetStreamwise_Periodic_RecoveredPressure(iPoint));
    if (heat && streamwisePeriodic_temperature)
      SetVolumeOutputValue("RECOVERED_TEMPERATURE", iPoint, Node_Flow->GetStreamwise_Periodic_RecoveredTemperature(iPoint));
  }

  LoadCommonFVMOutputs(config, geometry, iPoint);
}

void CFlowIncOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){

  if ((config->GetKind_Solver() == INC_NAVIER_STOKES) || (config->GetKind_Solver()  == INC_RANS)) {
    SetVolumeOutputValue("SKIN_FRICTION-X", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0));
    SetVolumeOutputValue("SKIN_FRICTION-Y", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1));
    if (nDim == 3)
      SetVolumeOutputValue("SKIN_FRICTION-Z", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2));
    if (weakly_coupled_heat)
      SetVolumeOutputValue("HEAT_FLUX", iPoint, solver[HEAT_SOL]->GetHeatFlux(iMarker, iVertex));
    else {
      SetVolumeOutputValue("HEAT_FLUX", iPoint, solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex));
    }
    SetVolumeOutputValue("Y_PLUS", iPoint, solver[FLOW_SOL]->GetYPlus(iMarker, iVertex));
  }

}

bool CFlowIncOutput::SetInit_Residuals(const CConfig *config){

  return (config->GetTime_Marching() != TIME_MARCHING::STEADY && (curInnerIter == 0))||
         (config->GetTime_Marching() == TIME_MARCHING::STEADY && (curInnerIter < 2));

}
