/*!
 * \file output_direct_mean_inc.cpp
 * \brief Main subroutines for incompressible flow output
 * \author R. Sanchez
 * \version 6.0.1 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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


#include "../include/output_structure.hpp"

CIncFlowOutput::CIncFlowOutput(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) : COutput(config) {

  nDim = geometry->GetnDim();
  
  turb_model = config->GetKind_Turb_Model();
  
  heat = config->GetEnergy_Equation();
  
  weakly_coupled_heat = config->GetWeakly_Coupled_Heat();
  
  grid_movement = config->GetGrid_Movement(); 
  
  su2double Gas_Constant, Mach2Vel, Mach_Motion;
  unsigned short iDim;
  su2double Gamma = config->GetGamma();
      
  /*--- Set the non-dimensionalization for coefficients. ---*/
  
  RefArea = config->GetRefArea();
  
  if (grid_movement) {
    Gas_Constant = config->GetGas_ConstantND();
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
  }
  RefDensity  = solver[FLOW_SOL]->GetDensity_Inf();
  RefPressure = solver[FLOW_SOL]->GetPressure_Inf();
  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
  
}

CIncFlowOutput::~CIncFlowOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }


}


void CIncFlowOutput::SetHistoryOutputFields(CConfig *config){
  
  // Iteration numbers
  AddOutputField("INT_ITER",   "Int_Iter",  FORMAT_INTEGER, "INT_ITER");
  AddOutputField("EXT_ITER",   "Ext_Iter",  FORMAT_INTEGER, "EXT_ITER");
  
  // Residuals
  AddOutputField("PRESSURE",   "Res[P]", FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("VELOCITY-X", "Res[U]", FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("VELOCITY-Y", "Res[V]", FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("VELOCITY-Z", "Res[W]", FORMAT_FIXED,   "RESIDUALS");
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    AddOutputField("NU_TILDE", "Res[nu]", FORMAT_FIXED, "RESIDUALS");
    break;  
  case SST:
    AddOutputField("KINETIC_ENERGY", "Res[k]", FORMAT_FIXED, "RESIDUALS");
    AddOutputField("DISSIPATION",    "Res[w]", FORMAT_FIXED, "RESIDUALS");
    break;
  }
  
  if (heat || weakly_coupled_heat){
    AddOutputField("HEAT", "Res[T]", FORMAT_FIXED, "RESIDUALS");
    AddOutputField("HEATFLUX", "HF(Total)",      FORMAT_SCIENTIFIC, "HEAT");
    AddOutputField("HEATFLUX_MAX", "HF(Max)",    FORMAT_SCIENTIFIC, "HEAT");
    AddOutputField("TEMPERATURE", "Temp(Total)", FORMAT_SCIENTIFIC, "HEAT");
  }
  
  // Aerodynamic coefficients
  AddOutputField("DRAG",       "CD(Total)",   FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("LIFT",       "CL(Total)",   FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("SIDEFORCE",  "CSF(Total)",  FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("MOMENT-X",   "CMx(Total)",  FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("MOMENT-Y",   "CMy(Total)",  FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("MOMENT-Z",   "CMz(Total)",  FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("FORCE-X",    "CFx(Total)",  FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("FORCE-Y",    "CFy(Total)",  FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("FORCE-Z",    "CFz(Total)",  FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("EFFICIENCY", "CEff(Total)", FORMAT_SCIENTIFIC, "AERO_COEFF");
  
  vector<string> Marker_Monitoring;
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
  }
  
  // Aerodynamic coefficients (per surface)  
  AddOutputPerSurfaceField("DRAG_ON_SURFACE",       "CD",   FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("LIFT_ON_SURFACE",       "CL",   FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("SIDEFORCE_ON_SURFACE",  "CSF",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("MOMENT-X_ON_SURFACE",   "CMx",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("MOMENT-Y_ON_SURFACE",   "CMy",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("MOMENT-Z_ON_SURFACE",   "CMz",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("FORCE-X_ON_SURFACE",    "CFx",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("FORCE-Y_ON_SURFACE",    "CFy",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("FORCE-Z_ON_SURFACE",    "CFz",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  AddOutputPerSurfaceField("EFFICIENCY_ON_SURFACE", "CEff", FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring);
  
  
  // Misc.
  AddOutputField("AOA",         "AoA",                      FORMAT_SCIENTIFIC, "AOA");
  AddOutputField("PHYS_TIME",   "Time(min)",                FORMAT_SCIENTIFIC, "PHYS_TIME");
  AddOutputField("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER,    "LINSOL_ITER");
  
  // Surface output
  vector<string> Marker_Analyze;
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < config->GetnMarker_Analyze(); iMarker_Analyze++){
    Marker_Analyze.push_back(config->GetMarker_Analyze_TagBound(iMarker_Analyze));
  }  
  AddOutputPerSurfaceField("AVG_MASSFLOW",             "Avg_Massflow",              FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("AVG_MACH",                 "Avg_Mach",                  FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("AVG_TEMP",                 "Avg_Temp",                  FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("AVG_PRESS",                "Avg_Press",                 FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("AVG_DENSITY",              "Avg_Density",               FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("AVG_ENTHALPY",             "Avg_Enthalpy",              FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("AVG_NORMALVEL",            "Avg_NormalVel",             FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("UNIFORMITY",               "Uniformity",                FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("SECONDARY_STRENGTH",       "Secondary_Strength",        FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("MOMENTUM_DISTORTION",      "Momentum_Distortion",       FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("SECONDARY_OVER_UNIFORMITY", "Secondary_Over_Uniformity", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("AVG_TOTALTEMP",            "Avg_TotalTemp",             FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("AVG_TOTALPRESS",           "Avg_TotalPress",            FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);
  AddOutputPerSurfaceField("PRESSURE_DROP",            "Pressure_Drop",             FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze);

  
}
inline bool CIncFlowOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { 
  return true;
}

inline bool CIncFlowOutput::WriteScreen_Header(CConfig *config) {  
  bool write_header;
  write_header = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));
  
  return true;
}
inline bool CIncFlowOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {
  return true;
}


inline void CIncFlowOutput::LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) {
  
  SetOutputFieldValue("INT_ITER", config[val_iZone]->GetIntIter());
  SetOutputFieldValue("EXT_ITER", config[val_iZone]->GetExtIter());
  SetOutputFieldValue("PRESSURE", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(0)));
  SetOutputFieldValue("VELOCITY-X", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(1)));
  SetOutputFieldValue("VELOCITY-Y", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(2)));
  if (nDim == 3) SetOutputFieldValue("VELOCITY-Z", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(3)));
 
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetOutputFieldValue("NU_TILDE", log10(solver_container[val_iZone][val_iInst][MESH_0][TURB_SOL]->GetRes_RMS(0)));
    break;  
  case SST:
    SetOutputFieldValue("KINETIC_ENERGY", log10(solver_container[val_iZone][val_iInst][MESH_0][TURB_SOL]->GetRes_RMS(0)));
    SetOutputFieldValue("DISSIPATION",    log10(solver_container[val_iZone][val_iInst][MESH_0][TURB_SOL]->GetRes_RMS(1)));
    break;
  }
  if (weakly_coupled_heat){
    SetOutputFieldValue("HEATFLUX",     solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetTotal_HeatFlux());
    SetOutputFieldValue("HEATFLUX_MAX", solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetTotal_MaxHeatFlux());
    SetOutputFieldValue("TEMPERATURE",  solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetTotal_AvgTemperature());
    SetOutputFieldValue("HEAT",         log10(solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetRes_RMS(0)));
  }
  if (heat){
    SetOutputFieldValue("HEATFLUX",     solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_HeatFlux());
    SetOutputFieldValue("HEATFLUX_MAX", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_MaxHeatFlux());
    SetOutputFieldValue("TEMPERATURE",  solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_AvgTemperature());
    if (nDim == 3) SetOutputFieldValue("HEAT",         log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(4)));
    else           SetOutputFieldValue("HEAT",         log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(3)));

  }
  SetOutputFieldValue("DRAG", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CD());
  SetOutputFieldValue("LIFT", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CL());
  if (nDim == 3)
    SetOutputFieldValue("SIDEFORCE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CSF());
  SetOutputFieldValue("MOMENT-X", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CMx());
  SetOutputFieldValue("MOMENT-Y", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CMy());
  if (nDim == 3)
    SetOutputFieldValue("MOMENT-Z", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CMz());
  SetOutputFieldValue("FORCE-X", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CFx());
  SetOutputFieldValue("FORCE-Y", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CFy());
  if (nDim == 3)
    SetOutputFieldValue("FORCE-Z", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetTotal_CFz());
  
  SetOutputFieldValue("AOA", config[val_iZone]->GetAoA());
  SetOutputFieldValue("EFFICIENCY", Output_Fields["DRAG"].Value/Output_Fields["LIFT"].Value);
  SetOutputFieldValue("PHYS_TIME", timeused);
  SetOutputFieldValue("LINSOL_ITER", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetIterLinSolver());
  
  
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config[val_iZone]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SetOutputPerSurfaceFieldValue("DRAG_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CD(iMarker_Monitoring), iMarker_Monitoring);
    SetOutputPerSurfaceFieldValue("LIFT_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CL(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetOutputPerSurfaceFieldValue("SIDEFORCE_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CSF(iMarker_Monitoring), iMarker_Monitoring);
    SetOutputPerSurfaceFieldValue("MOMENT-X_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring), iMarker_Monitoring);
    SetOutputPerSurfaceFieldValue("MOMENT-Y_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetOutputPerSurfaceFieldValue("MOMENT-Z_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring), iMarker_Monitoring);
    SetOutputPerSurfaceFieldValue("FORCE-X_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CFx(iMarker_Monitoring), iMarker_Monitoring);
    SetOutputPerSurfaceFieldValue("FORCE-Y_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CFy(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetOutputPerSurfaceFieldValue("FORCE-Z_ON_SURFACE", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetSurface_CFz(iMarker_Monitoring), iMarker_Monitoring);    
  }
  
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < config[val_iZone]->GetnMarker_Analyze(); iMarker_Analyze++){
    
    SetOutputPerSurfaceFieldValue("AVG_MASSFLOW", config[val_iZone]->GetSurface_MassFlow(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("AVG_MACH",     config[val_iZone]->GetSurface_Mach(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("AVG_TEMP",     config[val_iZone]->GetSurface_Temperature(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("AVG_PRESS",    config[val_iZone]->GetSurface_Pressure(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("AVG_DENSITY",  config[val_iZone]->GetSurface_Density(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("AVG_ENTHALPY",  config[val_iZone]->GetSurface_Enthalpy(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("AVG_NORMALVEL",  config[val_iZone]->GetSurface_NormalVelocity(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("UNIFORMITY",  config[val_iZone]->GetSurface_Uniformity(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("SECONDARY_STRENGTH",  config[val_iZone]->GetSurface_SecondaryStrength(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("MOMENTUM_DISTORTION",  config[val_iZone]->GetSurface_MomentumDistortion(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("SECONDARY_OVER_UNIFORMITY",  config[val_iZone]->GetSurface_SecondOverUniform(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("AVG_TOTALTEMP",  config[val_iZone]->GetSurface_TotalTemperature(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("AVG_TOTALPRESS",  config[val_iZone]->GetSurface_TotalPressure(iMarker_Analyze), iMarker_Analyze);
    SetOutputPerSurfaceFieldValue("PRESSURE_DROP",  config[val_iZone]->GetSurface_PressureDrop(iMarker_Analyze), iMarker_Analyze);
    
    
  }

}


void CIncFlowOutput::SetVolumeOutputFields(CConfig *config){
  
  // Grid coordinates
  AddVolumeOutputField("COORD-X", "x", "COORDINATES");
  AddVolumeOutputField("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    AddVolumeOutputField("COORD-Z", "z", "COORDINATES");
  
  // Conservative variables
  AddVolumeOutputField("PRESSURE",   "Pressure",   "CONSERVATIVE");
  AddVolumeOutputField("VELOCITY-X", "Velocity-x", "CONSERVATIVE");
  AddVolumeOutputField("VELOCITY-Y", "Velocity-y", "CONSERVATIVE");
  if (nDim == 3)
    AddVolumeOutputField("VELOCITY-Z", "Velocity-z", "CONSERVATIVE");
  AddVolumeOutputField("TEMPERATURE",  "Temperature","CONSERVATIVE");  
  
  // Turbulent Residuals
  switch(config->GetKind_Turb_Model()){
  case SST:
    AddVolumeOutputField("TKE", "TKE", "CONSERVATIVE");
    AddVolumeOutputField("OMEGA", "Omega", "CONSERVATIVE");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutputField("NU_TILDE", "Nu_Tilde", "CONSERVATIVE");
    break;
  case NONE:
    break;
  }
  
  // Primitive variables
  AddVolumeOutputField("MACH",           "Mach",                 "PRIMITIVE");
  AddVolumeOutputField("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE");
  AddVolumeOutputField("DENSITY",        "Density",              "PRIMITIVE");
  
  if (config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    AddVolumeOutputField("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE");
    
    AddVolumeOutputField("SKIN_FRICTION-X", "Skin_Friction_Coefficient-x", "PRIMITIVE");
    AddVolumeOutputField("SKIN_FRICTION-Y", "Skin_Friction_Coefficient-y", "PRIMITIVE");
    if (nDim == 3)
      AddVolumeOutputField("SKIN_FRICTION-Z", "Skin_Friction_Coefficient-z", "PRIMITIVE");
    
    AddVolumeOutputField("HEAT_FLUX", "Heat_Flux", "PRIMITIVE");
    AddVolumeOutputField("Y_PLUS", "Y_Plus", "PRIMITIVE");
    
  }
  
  if (config->GetKind_Solver() == RANS) {
    AddVolumeOutputField("EDDY_VISCOSITY", "Eddy_Viscosity", "PRIMITIVE");
  }
  
  if (config->GetKind_Trans_Model() == BC){
    AddVolumeOutputField("INTERMITTENCY", "gamma_BC", "INTERMITTENCY");
  }

  //Residuals
  AddVolumeOutputField("RESIDUAL_PRESSURE", "Residual_Density", "RESIDUAL");
  AddVolumeOutputField("RESIDUAL_VELOCITY-X", "Residual_Momentum-x", "RESIDUAL");
  AddVolumeOutputField("RESIDUAL_VELOCITY-Y", "Residual_Momentum-y", "RESIDUAL");
  if (nDim == 3)
    AddVolumeOutputField("RESIDUAL_VELOCITY-Z", "Residual_Momentum-z", "RESIDUAL");
  AddVolumeOutputField("RESIDUAL_TEMPERATURE", "Residual_Energy", "RESIDUAL");
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    AddVolumeOutputField("RESIDUAL_TKE", "Residual_TKE", "RESIDUAL");
    AddVolumeOutputField("RESIDUAL_OMEGA", "Residual_Omega", "RESIDUAL");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutputField("RESIDUAL_NU_TILDE", "Residual_Nu_Tilde", "RESIDUAL");
    break;
  case NONE:
    break;
  }
  
  // Limiter values
  AddVolumeOutputField("LIMITER_PRESSURE", "Limiter_Pressure", "LIMITER");
  AddVolumeOutputField("LIMITER_VELOCITY-X", "Limiter_Velocity-x", "LIMITER");
  AddVolumeOutputField("LIMITER_VELOCITY-Y", "Limiter_Velocity-y", "LIMITER");
  if (nDim == 3)
    AddVolumeOutputField("LIMITER_VELOCITY-Z", "Limiter_Velocity-z", "LIMITER");
  AddVolumeOutputField("LIMITER_TEMPERATURE", "Limiter_Temperature", "LIMITER");
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    AddVolumeOutputField("LIMITER_TKE", "Limiter_TKE", "RESIDUAL");
    AddVolumeOutputField("LIMITER_OMEGA", "Limiter_Omega", "RESIDUAL");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutputField("LIMITER_NU_TILDE", "Limiter_Nu_Tilde", "RESIDUAL");
    break;
  case NONE:
    break;
  }
  
  // Hybrid RANS-LES
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    AddVolumeOutputField("DES_LENGTHSCALE", "DES_LengthScale", "DDES");
    AddVolumeOutputField("WALL_DISTANCE", "Wall_Distance", "DDES");
  }
  
  // Roe Low Dissipation
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    AddVolumeOutputField("ROE_DISSIPATION", "Roe_Dissipation", "ROE_DISSIPATION");
  }
}

void CIncFlowOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Flow = solver[FLOW_SOL]->node[iPoint]; 
  CVariable* Node_Turb;
  
  if (config->GetKind_Turb_Model() != NONE){
    Node_Turb = solver[TURB_SOL]->node[iPoint]; 
  }
  
  CPoint*    Node_Geo  = geometry->node[iPoint];
          
  SetVolumeOutputFieldValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputFieldValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputFieldValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
  
  SetVolumeOutputFieldValue("PRESSURE",    iPoint, Node_Flow->GetSolution(0));
  SetVolumeOutputFieldValue("VELOCITY-X", iPoint, Node_Flow->GetSolution(1));
  SetVolumeOutputFieldValue("VELOCITY-Y", iPoint, Node_Flow->GetSolution(2));
  if (nDim == 3){
    SetVolumeOutputFieldValue("VELOCITY-Z", iPoint, Node_Flow->GetSolution(3));
    SetVolumeOutputFieldValue("TEMPERATURE",     iPoint, Node_Flow->GetSolution(4));
  } else {
    SetVolumeOutputFieldValue("TEMPERATURE",     iPoint, Node_Flow->GetSolution(3));    
  }
  
  // Turbulent Residuals
  switch(config->GetKind_Turb_Model()){
  case SST:
    SetVolumeOutputFieldValue("TKE", iPoint, Node_Turb->GetSolution(0));
    SetVolumeOutputFieldValue("OMEGA", iPoint, Node_Turb->GetSolution(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputFieldValue("NU_TILDE", iPoint, Node_Turb->GetSolution(0));
    break;
  case NONE:
    break;
  }
  
  SetVolumeOutputFieldValue("MACH", iPoint, sqrt(Node_Flow->GetVelocity2())/Node_Flow->GetSoundSpeed());
  SetVolumeOutputFieldValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure() - RefPressure)*factor*RefArea);
  SetVolumeOutputFieldValue("DENSITY", iPoint, Node_Flow->GetDensity());
  
  if (config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    SetVolumeOutputFieldValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity());
  }
  
  if (config->GetKind_Solver() == RANS) {
    SetVolumeOutputFieldValue("EDDY_VISCOSITY", iPoint, Node_Flow->GetEddyViscosity());
  }
  
  if (config->GetKind_Trans_Model() == BC){
    SetVolumeOutputFieldValue("INTERMITTENCY", iPoint, Node_Turb->GetGammaBC());
  }
  
  SetVolumeOutputFieldValue("RESIDUAL_PRESSURE", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 0));
  SetVolumeOutputFieldValue("RESIDUAL_VELOCITY-X", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 1));
  SetVolumeOutputFieldValue("RESIDUAL_VELOCITY-Y", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputFieldValue("RESIDUAL_VELOCITY-Z", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 3));
    SetVolumeOutputFieldValue("RESIDUAL_TEMPERATURE", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 4));
  } else {
    SetVolumeOutputFieldValue("RESIDUAL_TEMPERATURE", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    SetVolumeOutputFieldValue("RESIDUAL_TKE", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 0));
    SetVolumeOutputFieldValue("RESIDUAL_OMEGA", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputFieldValue("RESIDUAL_NU_TILDE", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 0));
    break;
  case NONE:
    break;
  }
  
  SetVolumeOutputFieldValue("LIMITER_PRESSURE", iPoint, Node_Flow->GetLimiter_Primitive(0));
  SetVolumeOutputFieldValue("LIMITER_VELOCITY-X", iPoint, Node_Flow->GetLimiter_Primitive(1));
  SetVolumeOutputFieldValue("LIMITER_VELOCITY-Y", iPoint, Node_Flow->GetLimiter_Primitive(2));
  if (nDim == 3){
    SetVolumeOutputFieldValue("LIMITER_VELOCITY-Z", iPoint, Node_Flow->GetLimiter_Primitive(3));
    SetVolumeOutputFieldValue("LIMITER_TEMPERATURE", iPoint, Node_Flow->GetLimiter_Primitive(4));
  } else {
    SetVolumeOutputFieldValue("LIMITER_TEMPERATURE", iPoint, Node_Flow->GetLimiter_Primitive(3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    SetVolumeOutputFieldValue("LIMITER_TKE", iPoint, Node_Turb->GetLimiter_Primitive(0));
    SetVolumeOutputFieldValue("LIMITER_OMEGA", iPoint, Node_Turb->GetLimiter_Primitive(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputFieldValue("LIMITER_NU_TILDE", iPoint, Node_Turb->GetLimiter_Primitive(0));
    break;
  case NONE:
    break;
  }
  
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    SetVolumeOutputFieldValue("DES_LENGTHSCALE", iPoint, Node_Flow->GetDES_LengthScale());
    SetVolumeOutputFieldValue("WALL_DISTANCE", iPoint, Node_Geo->GetWall_Distance());
  }
  
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    SetVolumeOutputFieldValue("ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation());
  }  
}

void CIncFlowOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){
  
  if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver()  == RANS)) {  
    SetVolumeOutputFieldValue("SKIN_FRICTION-X", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0));
    SetVolumeOutputFieldValue("SKIN_FRICTION-Y", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1));
    if (nDim == 3)
      SetVolumeOutputFieldValue("SKIN_FRICTION-Z", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2));
  
    SetVolumeOutputFieldValue("HEAT_FLUX", iPoint, solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex));
    SetVolumeOutputFieldValue("Y_PLUS", iPoint, solver[FLOW_SOL]->GetYPlus(iMarker, iVertex));
  }
}
