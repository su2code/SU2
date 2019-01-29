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
  
  /*--- Set the default history fields if nothing is set in the config file ---*/
  
  if (nRequestedHistoryFields == 0){
    RequestedHistoryFields.push_back("ITER");
    RequestedHistoryFields.push_back("RMS_RES");
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }
  
  if (nRequestedScreenFields == 0){
    if (multizone) RequestedScreenFields.push_back("OUTER_ITER");
    RequestedScreenFields.push_back("INNER_ITER");
    RequestedScreenFields.push_back("RMS_PRESSURE");
    RequestedScreenFields.push_back("RMS_VELOCITY-X");
    RequestedScreenFields.push_back("RMS_VELOCITY-Y");
    nRequestedScreenFields = RequestedScreenFields.size();
  }
  
  if (nRequestedVolumeFields == 0){
    RequestedVolumeFields.push_back("COORDINATES");
    RequestedVolumeFields.push_back("CONSERVATIVE");
    RequestedVolumeFields.push_back("PRIMITIVE");
    nRequestedVolumeFields = RequestedVolumeFields.size();
  }
  
  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Incomp. Fluid)";
  MultiZoneHeaderString = ss.str();
  
}

CIncFlowOutput::~CIncFlowOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }


}


void CIncFlowOutput::SetHistoryOutputFields(CConfig *config){
  
  /// BEGIN_GROUP: ITERATION, DESCRIPTION: Iteration identifier.
  /// DESCRIPTION: The time iteration index.
  AddHistoryOutput("TIME_ITER",     "Time_Iter",  FORMAT_INTEGER, "ITER"); 
  /// DESCRIPTION: The internal iteration index.
  AddHistoryOutput("OUTER_ITER",   "Outer_Iter",  FORMAT_INTEGER, "ITER"); 
  /// DESCRIPTION: The external iteration index.
  AddHistoryOutput("INNER_ITER",   "Inner_Iter", FORMAT_INTEGER,  "ITER"); 
  /// END_GROUP

  /// DESCRIPTION: Currently used wall-clock time.
  AddHistoryOutput("PHYS_TIME",   "Time(min)", FORMAT_SCIENTIFIC, "PHYS_TIME"); 
  
  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the conservative variables. 
  /// DESCRIPTION: Root-mean square residual of the pressure.
  AddHistoryOutput("RMS_PRESSURE",   "rms[P]", FORMAT_FIXED,   "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the velocity x-component.  
  AddHistoryOutput("RMS_VELOCITY-X", "rms[U]", FORMAT_FIXED,   "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the velocity y-component.  
  AddHistoryOutput("RMS_VELOCITY-Y", "rms[V]", FORMAT_FIXED,   "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the velocity z-component.  
  AddHistoryOutput("RMS_VELOCITY-Z", "rms[W]", FORMAT_FIXED,   "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the temperature.
  AddHistoryOutput("RMS_HEAT", "rms[T]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).  
  AddHistoryOutput("RMS_NU_TILDE",       "rms[nu]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).    
  AddHistoryOutput("RMS_KINETIC_ENERGY", "rms[k]",  FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).    
  AddHistoryOutput("RMS_DISSIPATION",    "rms[w]",  FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);
  /// END_GROUP
  
  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the conservative variables. 
  /// DESCRIPTION: Maximum residual of the pressure.
  AddHistoryOutput("MAX_PRESSURE",   "max[P]", FORMAT_FIXED,   "MAX_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity x-component.   
  AddHistoryOutput("MAX_VELOCITY-X", "max[U]", FORMAT_FIXED,   "MAX_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity y-component.   
  AddHistoryOutput("MAX_VELOCITY-Y", "max[V]", FORMAT_FIXED,   "MAX_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the velocity z-component.   
  AddHistoryOutput("MAX_VELOCITY-Z", "max[W]", FORMAT_FIXED,   "MAX_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the temperature.
  AddHistoryOutput("MAX_HEAT", "max[T]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of nu tilde (SA model).
  AddHistoryOutput("MAX_NU_TILDE",       "max[nu]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of kinetic energy (SST model). 
  AddHistoryOutput("MAX_KINETIC_ENERGY", "max[k]",  FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the dissipation (SST model).   
  AddHistoryOutput("MAX_DISSIPATION",    "max[w]",  FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);  
  /// END_GROUP
  
  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the conservative variables. 
  /// DESCRIPTION: BGS residual of the pressure.
  AddHistoryOutput("BGS_PRESSURE",   "bgs[P]", FORMAT_FIXED,   "BGS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: BGS residual of the momentum x-component.  
  AddHistoryOutput("BGS_VELOCITY-X", "bgs[U]", FORMAT_FIXED,   "BGS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: BGS residual of the momentum x-component.  
  AddHistoryOutput("BGS_VELOCITY-Y", "bgs[V]", FORMAT_FIXED,   "BGS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: BGS residual of the momentum x-component.  
  AddHistoryOutput("BGS_VELOCITY-Z", "bgs[W]", FORMAT_FIXED,   "BGS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the temperature.
  AddHistoryOutput("BGS_HEAT", "bgs[T]", FORMAT_FIXED, "BGS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: BGS residual of nu tilde (SA model).  
  AddHistoryOutput("BGS_NU_TILDE",       "bgs[nu]", FORMAT_FIXED, "BGS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: BGS residual of kinetic energy (SST model).    
  AddHistoryOutput("BGS_KINETIC_ENERGY", "bgs[k]",  FORMAT_FIXED, "BGS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: BGS residual of the dissipation (SST model).    
  AddHistoryOutput("BGS_DISSIPATION",    "bgs[w]",  FORMAT_FIXED, "BGS_RES", TYPE_RESIDUAL);
  /// END_GROUP

  /// BEGIN_GROUP: AERO_COEFF, DESCRIPTION: Sum of the aerodynamic coefficients and forces on all surfaces (markers) set with MARKER_MONITORING.
  /// DESCRIPTION: Drag coefficient 
  AddHistoryOutput("DRAG",       "CD",   FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Lift coefficient 
  AddHistoryOutput("LIFT",       "CL",   FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Sideforce coefficient   
  AddHistoryOutput("SIDEFORCE",  "CSF",  FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the x-axis    
  AddHistoryOutput("MOMENT-X",   "CMx",  FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the y-axis    
  AddHistoryOutput("MOMENT-Y",   "CMy",  FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the z-axis      
  AddHistoryOutput("MOMENT-Z",   "CMz",  FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in x direction    
  AddHistoryOutput("FORCE-X",    "CFx",  FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in y direction    
  AddHistoryOutput("FORCE-Y",    "CFy",  FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in z direction      
  AddHistoryOutput("FORCE-Z",    "CFz",  FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// DESCRIPTION: Lift-to-drag ratio
  AddHistoryOutput("EFFICIENCY", "CEff", FORMAT_SCIENTIFIC, "AERO_COEFF", TYPE_COEFFICIENT);
  /// END_GROUP
 
  /// BEGIN_GROUP: HEAT_COEFF, DESCRIPTION: Heat coefficients on all surfaces set with MARKER_MONITORING.
  /// DESCRIPTION: Total heatflux
  AddHistoryOutput("HEATFLUX", "HF",      FORMAT_SCIENTIFIC, "HEAT", TYPE_COEFFICIENT);
  /// DESCRIPTION: Maximal heatflux  
  AddHistoryOutput("HEATFLUX_MAX", "maxHF",    FORMAT_SCIENTIFIC, "HEAT", TYPE_COEFFICIENT);
  /// DESCRIPTION: Temperature
  AddHistoryOutput("TEMPERATURE", "Temp", FORMAT_SCIENTIFIC, "HEAT", TYPE_COEFFICIENT);
  /// END_GROUP
  
  
  /// BEGIN_GROUP: AERO_COEFF_SURF, DESCRIPTION: Aerodynamic coefficients and forces per surface.
  vector<string> Marker_Monitoring;
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
  }  
  
  /// DESCRIPTION: Drag coefficient   
  AddHistoryOutputPerSurface("DRAG_ON_SURFACE",       "CD",   FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Lift coefficient   
  AddHistoryOutputPerSurface("LIFT_ON_SURFACE",       "CL",   FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Sideforce coefficient     
  AddHistoryOutputPerSurface("SIDEFORCE_ON_SURFACE",  "CSF",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the x-axis      
  AddHistoryOutputPerSurface("MOMENT-X_ON_SURFACE",   "CMx",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the y-axis      
  AddHistoryOutputPerSurface("MOMENT-Y_ON_SURFACE",   "CMy",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the z-axis        
  AddHistoryOutputPerSurface("MOMENT-Z_ON_SURFACE",   "CMz",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in x direction      
  AddHistoryOutputPerSurface("FORCE-X_ON_SURFACE",    "CFx",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in y direction      
  AddHistoryOutputPerSurface("FORCE-Y_ON_SURFACE",    "CFy",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in z direction        
  AddHistoryOutputPerSurface("FORCE-Z_ON_SURFACE",    "CFz",  FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Lift-to-drag ratio  
  AddHistoryOutputPerSurface("EFFICIENCY_ON_SURFACE", "CEff", FORMAT_SCIENTIFIC, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// END_GROUP 
  
  /// DESCRIPTION: Angle of attack  
  AddHistoryOutput("AOA",         "AoA",                      FORMAT_SCIENTIFIC, "AOA");
  /// DESCRIPTION: Linear solver iterations   
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER,    "LINSOL_ITER");
  
  /// BEGIN_GROUP: AERO_COEFF_SURF, DESCRIPTION: Surface values on non-solid markers.
  vector<string> Marker_Analyze;
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < config->GetnMarker_Analyze(); iMarker_Analyze++){
    Marker_Analyze.push_back(config->GetMarker_Analyze_TagBound(iMarker_Analyze));
  }  
  /// DESCRIPTION: Average mass flow    
  AddHistoryOutputPerSurface("AVG_MASSFLOW",             "Avg_Massflow",              FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Mach number      
  AddHistoryOutputPerSurface("AVG_MACH",                 "Avg_Mach",                  FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Temperature        
  AddHistoryOutputPerSurface("AVG_TEMP",                 "Avg_Temp",                  FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Pressure  
  AddHistoryOutputPerSurface("AVG_PRESS",                "Avg_Press",                 FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Density  
  AddHistoryOutputPerSurface("AVG_DENSITY",              "Avg_Density",               FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Enthalpy  
  AddHistoryOutputPerSurface("AVG_ENTHALPY",             "Avg_Enthalpy",              FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average velocity in normal direction of the surface
  AddHistoryOutputPerSurface("AVG_NORMALVEL",            "Avg_NormalVel",             FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Flow uniformity 
  AddHistoryOutputPerSurface("UNIFORMITY",               "Uniformity",                FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Secondary strength
  AddHistoryOutputPerSurface("SECONDARY_STRENGTH",       "Secondary_Strength",        FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Momentum distortion  
  AddHistoryOutputPerSurface("MOMENTUM_DISTORTION",      "Momentum_Distortion",       FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Secondary over uniformity 
  AddHistoryOutputPerSurface("SECONDARY_OVER_UNIFORMITY", "Secondary_Over_Uniformity", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average total temperature  
  AddHistoryOutputPerSurface("AVG_TOTALTEMP",            "Avg_TotalTemp",             FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average total pressure   
  AddHistoryOutputPerSurface("AVG_TOTALPRESS",           "Avg_TotalPress",            FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Pressure drop    
  AddHistoryOutputPerSurface("PRESSURE_DROP",            "Pressure_Drop",             FORMAT_SCIENTIFIC, "SURFACE_OUTPUT", Marker_Analyze, TYPE_COEFFICIENT);
  /// END_GROUP
  
}

bool CIncFlowOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { 
 if (!write_dualtime){
   return true;
 }
 else {
   return false;
 }
}

bool CIncFlowOutput::WriteScreen_Header(CConfig *config) {  
  bool write_header = false;
  if (config->GetUnsteady_Simulation() == STEADY || config->GetUnsteady_Simulation() == TIME_STEPPING) {
    write_header = ((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0) || (config->GetMultizone_Problem() && config->GetInnerIter() == 0);
  } else {
    write_header = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND) && config->GetIntIter() == 0;
  }

  /*--- For multizone problems, print the header only if requested explicitly (default of GetWrt_ZoneConv is false) ---*/
  if(config->GetMultizone_Problem()) write_header = (write_header && config->GetWrt_ZoneConv());

  return write_header;
}

bool CIncFlowOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {
  bool write_output = false;
  
  if (((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) ) 
      && write_dualtime ){
    write_output = (config->GetIntIter() % config->GetWrt_Con_Freq_DualTime() == 0);
  }
  else if (((config->GetUnsteady_Simulation() == STEADY) || (config->GetUnsteady_Simulation() == TIME_STEPPING) )){
    write_output = (config->GetInnerIter() % config->GetWrt_Con_Freq() == 0) ;    
  } 

  /*--- For multizone problems, print the body only if requested explicitly (default of GetWrt_ZoneConv is false) ---*/
  if(config->GetMultizone_Problem()) write_output = (write_output && config->GetWrt_ZoneConv());

  return write_output;
}



inline void CIncFlowOutput::LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) {
  
  CSolver* flow_solver = solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL];
  CSolver* turb_solver = solver_container[val_iZone][val_iInst][MESH_0][TURB_SOL];  
  CSolver* heat_solver = solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL];
  
  SetHistoryOutputValue("INNER_ITER", config[val_iZone]->GetInnerIter());
  SetHistoryOutputValue("OUTER_ITER", config[val_iZone]->GetOuterIter());    
  //SetHistoryOutputValue("EXT_ITER", config[val_iZone]->GetOuterIter());
  
  SetHistoryOutputValue("RMS_PRESSURE", log10(flow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_VELOCITY-X", log10(flow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_VELOCITY-Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 3) SetHistoryOutputValue("RMS_VELOCITY-Z", log10(flow_solver->GetRes_RMS(3)));
 
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("RMS_NU_TILDE", log10(turb_solver->GetRes_RMS(0)));
    break;  
  case SST:
    SetHistoryOutputValue("RMS_KINETIC_ENERGY", log10(turb_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_DISSIPATION",    log10(turb_solver->GetRes_RMS(1)));
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
  case SST:
    SetHistoryOutputValue("MAX_KINETIC_ENERGY", log10(turb_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_DISSIPATION",    log10(turb_solver->GetRes_Max(1)));
    break;
  }
  
  if (config[val_iZone]->GetMultizone_Problem()){
    SetHistoryOutputValue("BGS_PRESSURE", log10(flow_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_VELOCITY-X", log10(flow_solver->GetRes_BGS(1)));
    SetHistoryOutputValue("BGS_VELOCITY-Y", log10(flow_solver->GetRes_BGS(2)));
    if (nDim == 3) SetHistoryOutputValue("BGS_VELOCITY-Z", log10(flow_solver->GetRes_BGS(3)));
    
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      SetHistoryOutputValue("BGS_NU_TILDE", log10(turb_solver->GetRes_BGS(0)));
      break;  
    case SST:
      SetHistoryOutputValue("BGS_KINETIC_ENERGY", log10(turb_solver->GetRes_BGS(0)));
      SetHistoryOutputValue("BGS_DISSIPATION",    log10(turb_solver->GetRes_BGS(1)));
      break;
    }
  }
  
  if (weakly_coupled_heat){
    SetHistoryOutputValue("HEATFLUX",     heat_solver->GetTotal_HeatFlux());
    SetHistoryOutputValue("HEATFLUX_MAX", heat_solver->GetTotal_MaxHeatFlux());
    SetHistoryOutputValue("TEMPERATURE",  heat_solver->GetTotal_AvgTemperature());
    SetHistoryOutputValue("RMS_HEAT",         log10(heat_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("MAX_HEAT",         log10(heat_solver->GetRes_Max(0)));
    if (config[val_iZone]->GetMultizone_Problem()) {SetHistoryOutputValue("BGS_HEAT", log10(heat_solver->GetRes_BGS(0)));}
  }
  if (heat){
    SetHistoryOutputValue("HEATFLUX",     flow_solver->GetTotal_HeatFlux());
    SetHistoryOutputValue("HEATFLUX_MAX", flow_solver->GetTotal_MaxHeatFlux());
    SetHistoryOutputValue("TEMPERATURE",  flow_solver->GetTotal_AvgTemperature());
    if (nDim == 3) SetHistoryOutputValue("RMS_HEAT",         log10(flow_solver->GetRes_RMS(4)));
    else           SetHistoryOutputValue("RMS_HEAT",         log10(flow_solver->GetRes_RMS(3)));
    
    if (nDim == 3) SetHistoryOutputValue("MAX_HEAT",         log10(flow_solver->GetRes_Max(4)));
    else           SetHistoryOutputValue("MAX_HEAT",         log10(flow_solver->GetRes_Max(3)));

  }
  SetHistoryOutputValue("DRAG", flow_solver->GetTotal_CD());
  SetHistoryOutputValue("LIFT", flow_solver->GetTotal_CL());
  if (nDim == 3)
    SetHistoryOutputValue("SIDEFORCE", flow_solver->GetTotal_CSF());
  SetHistoryOutputValue("MOMENT-X", flow_solver->GetTotal_CMx());
  SetHistoryOutputValue("MOMENT-Y", flow_solver->GetTotal_CMy());
  if (nDim == 3)
    SetHistoryOutputValue("MOMENT-Z", flow_solver->GetTotal_CMz());
  SetHistoryOutputValue("FORCE-X", flow_solver->GetTotal_CFx());
  SetHistoryOutputValue("FORCE-Y", flow_solver->GetTotal_CFy());
  if (nDim == 3)
    SetHistoryOutputValue("FORCE-Z", flow_solver->GetTotal_CFz());
  
  SetHistoryOutputValue("AOA", config[val_iZone]->GetAoA());
  SetHistoryOutputValue("EFFICIENCY", HistoryOutput_Map["DRAG"].Value/HistoryOutput_Map["LIFT"].Value);
  SetHistoryOutputValue("PHYS_TIME", timeused);
  SetHistoryOutputValue("LINSOL_ITER", flow_solver->GetIterLinSolver());
  
  
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config[val_iZone]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SetHistoryOutputPerSurfaceValue("DRAG_ON_SURFACE", flow_solver->GetSurface_CD(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("LIFT_ON_SURFACE", flow_solver->GetSurface_CL(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetHistoryOutputPerSurfaceValue("SIDEFORCE_ON_SURFACE", flow_solver->GetSurface_CSF(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("MOMENT-X_ON_SURFACE", flow_solver->GetSurface_CMx(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("MOMENT-Y_ON_SURFACE", flow_solver->GetSurface_CMy(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetHistoryOutputPerSurfaceValue("MOMENT-Z_ON_SURFACE", flow_solver->GetSurface_CMz(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("FORCE-X_ON_SURFACE", flow_solver->GetSurface_CFx(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("FORCE-Y_ON_SURFACE", flow_solver->GetSurface_CFy(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetHistoryOutputPerSurfaceValue("FORCE-Z_ON_SURFACE", flow_solver->GetSurface_CFz(iMarker_Monitoring), iMarker_Monitoring);    
  }
  
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < config[val_iZone]->GetnMarker_Analyze(); iMarker_Analyze++){
    
    SetHistoryOutputPerSurfaceValue("AVG_MASSFLOW", config[val_iZone]->GetSurface_MassFlow(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("AVG_MACH",     config[val_iZone]->GetSurface_Mach(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("AVG_TEMP",     config[val_iZone]->GetSurface_Temperature(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("AVG_PRESS",    config[val_iZone]->GetSurface_Pressure(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("AVG_DENSITY",  config[val_iZone]->GetSurface_Density(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("AVG_ENTHALPY",  config[val_iZone]->GetSurface_Enthalpy(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("AVG_NORMALVEL",  config[val_iZone]->GetSurface_NormalVelocity(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("UNIFORMITY",  config[val_iZone]->GetSurface_Uniformity(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("SECONDARY_STRENGTH",  config[val_iZone]->GetSurface_SecondaryStrength(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("MOMENTUM_DISTORTION",  config[val_iZone]->GetSurface_MomentumDistortion(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("SECONDARY_OVER_UNIFORMITY",  config[val_iZone]->GetSurface_SecondOverUniform(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("AVG_TOTALTEMP",  config[val_iZone]->GetSurface_TotalTemperature(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("AVG_TOTALPRESS",  config[val_iZone]->GetSurface_TotalPressure(iMarker_Analyze), iMarker_Analyze);
    SetHistoryOutputPerSurfaceValue("PRESSURE_DROP",  config[val_iZone]->GetSurface_PressureDrop(iMarker_Analyze), iMarker_Analyze);
    
    
  }

}


void CIncFlowOutput::SetVolumeOutputFields(CConfig *config){
  
  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES");
  
  // Conservative variables
  AddVolumeOutput("PRESSURE",   "Pressure",   "CONSERVATIVE");
  AddVolumeOutput("VELOCITY-X", "Velocity_x", "CONSERVATIVE");
  AddVolumeOutput("VELOCITY-Y", "Velocity_y", "CONSERVATIVE");
  if (nDim == 3)
    AddVolumeOutput("VELOCITY-Z", "Velocity_z", "CONSERVATIVE");
  if (config->GetEnergy_Equation())
    AddVolumeOutput("TEMPERATURE",  "Temperature","CONSERVATIVE");  
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    AddVolumeOutput("TKE", "TKE", "CONSERVATIVE");
    AddVolumeOutput("OMEGA", "Omega", "CONSERVATIVE");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutput("NU_TILDE", "Nu_Tilde", "CONSERVATIVE");
    break;
  case NONE:
    break;
  }
  
  // Grid velocity
  if (config->GetGrid_Movement()){
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY");
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY");
    if (nDim == 3 ) 
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY");
  }
  
  // Primitive variables
  AddVolumeOutput("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE");
  AddVolumeOutput("DENSITY",        "Density",              "PRIMITIVE");
  
  if (config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    AddVolumeOutput("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE");
    
    AddVolumeOutput("SKIN_FRICTION-X", "Skin_Friction_Coefficient_x", "PRIMITIVE");
    AddVolumeOutput("SKIN_FRICTION-Y", "Skin_Friction_Coefficient_y", "PRIMITIVE");
    if (nDim == 3)
      AddVolumeOutput("SKIN_FRICTION-Z", "Skin_Friction_Coefficient_z", "PRIMITIVE");
    
    AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE");
    AddVolumeOutput("Y_PLUS", "Y_Plus", "PRIMITIVE");
    
  }
  
  if (config->GetKind_Solver() == RANS) {
    AddVolumeOutput("EDDY_VISCOSITY", "Eddy_Viscosity", "PRIMITIVE");
  }
  
  if (config->GetKind_Trans_Model() == BC){
    AddVolumeOutput("INTERMITTENCY", "gamma_BC", "INTERMITTENCY");
  }

  //Residuals
  AddVolumeOutput("RESIDUAL_PRESSURE", "Residual_Density", "RESIDUAL");
  AddVolumeOutput("RESIDUAL_VELOCITY-X", "Residual_Momentum_x", "RESIDUAL");
  AddVolumeOutput("RESIDUAL_VELOCITY-Y", "Residual_Momentum_y", "RESIDUAL");
  if (nDim == 3)
    AddVolumeOutput("RESIDUAL_VELOCITY-Z", "Residual_Momentum_z", "RESIDUAL");
  AddVolumeOutput("RESIDUAL_TEMPERATURE", "Residual_Energy", "RESIDUAL");
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    AddVolumeOutput("RESIDUAL_TKE", "Residual_TKE", "RESIDUAL");
    AddVolumeOutput("RESIDUAL_OMEGA", "Residual_Omega", "RESIDUAL");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutput("RESIDUAL_NU_TILDE", "Residual_Nu_Tilde", "RESIDUAL");
    break;
  case NONE:
    break;
  }
  
  // Limiter values
  AddVolumeOutput("LIMITER_PRESSURE", "Limiter_Pressure", "LIMITER");
  AddVolumeOutput("LIMITER_VELOCITY-X", "Limiter_Velocity_x", "LIMITER");
  AddVolumeOutput("LIMITER_VELOCITY-Y", "Limiter_Velocity_y", "LIMITER");
  if (nDim == 3)
    AddVolumeOutput("LIMITER_VELOCITY-Z", "Limiter_Velocity_z", "LIMITER");
  AddVolumeOutput("LIMITER_TEMPERATURE", "Limiter_Temperature", "LIMITER");
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    AddVolumeOutput("LIMITER_TKE", "Limiter_TKE", "RESIDUAL");
    AddVolumeOutput("LIMITER_OMEGA", "Limiter_Omega", "RESIDUAL");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutput("LIMITER_NU_TILDE", "Limiter_Nu_Tilde", "RESIDUAL");
    break;
  case NONE:
    break;
  }
  
  // Hybrid RANS-LES
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    AddVolumeOutput("DES_LENGTHSCALE", "DES_LengthScale", "DDES");
    AddVolumeOutput("WALL_DISTANCE", "Wall_Distance", "DDES");
  }
  
  // Roe Low Dissipation
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    AddVolumeOutput("ROE_DISSIPATION", "Roe_Dissipation", "ROE_DISSIPATION");
  }
  
  if(config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    if (nDim == 3){
      AddVolumeOutput("VORTICITY_X", "Vorticity_x", "VORTEX_IDENTIFICATION");
      AddVolumeOutput("VORTICITY_Y", "Vorticity_y", "VORTEX_IDENTIFICATION");
    }
    AddVolumeOutput("VORTICITY_Z", "Vorticity_z", "VORTEX_IDENTIFICATION");
    AddVolumeOutput("Q_CRITERION", "Q_Criterion", "VORTEX_IDENTIFICATION");  
  }
}

void CIncFlowOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Flow = solver[FLOW_SOL]->node[iPoint]; 
  CVariable* Node_Turb = NULL;
  
  if (config->GetKind_Turb_Model() != NONE){
    Node_Turb = solver[TURB_SOL]->node[iPoint]; 
  }
  
  CPoint*    Node_Geo  = geometry->node[iPoint];
          
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
  
  SetVolumeOutputValue("PRESSURE",    iPoint, Node_Flow->GetSolution(0));
  SetVolumeOutputValue("VELOCITY-X", iPoint, Node_Flow->GetSolution(1));
  SetVolumeOutputValue("VELOCITY-Y", iPoint, Node_Flow->GetSolution(2));
  if (nDim == 3){
    SetVolumeOutputValue("VELOCITY-Z", iPoint, Node_Flow->GetSolution(3));
    SetVolumeOutputValue("TEMPERATURE",     iPoint, Node_Flow->GetSolution(4));
  } else {
    SetVolumeOutputValue("TEMPERATURE",     iPoint, Node_Flow->GetSolution(3));    
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    SetVolumeOutputValue("TKE", iPoint, Node_Turb->GetSolution(0));
    SetVolumeOutputValue("OMEGA", iPoint, Node_Turb->GetSolution(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("NU_TILDE", iPoint, Node_Turb->GetSolution(0));
    break;
  case NONE:
    break;
  }
  
  if (config->GetGrid_Movement()){
    SetVolumeOutputValue("GRID_VELOCITY-X", iPoint, Node_Geo->GetGridVel()[0]);
    SetVolumeOutputValue("GRID_VELOCITY-Y", iPoint, Node_Geo->GetGridVel()[1]);
    if (nDim == 3)
      SetVolumeOutputValue("GRID_VELOCITY-Z", iPoint, Node_Geo->GetGridVel()[2]);
  }
  
  SetVolumeOutputValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure() - RefPressure)*factor*RefArea);
  SetVolumeOutputValue("DENSITY", iPoint, Node_Flow->GetDensity());
  
  if (config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity());
  }
  
  if (config->GetKind_Solver() == RANS) {
    SetVolumeOutputValue("EDDY_VISCOSITY", iPoint, Node_Flow->GetEddyViscosity());
  }
  
  if (config->GetKind_Trans_Model() == BC){
    SetVolumeOutputValue("INTERMITTENCY", iPoint, Node_Turb->GetGammaBC());
  }
  
  SetVolumeOutputValue("RESIDUAL_PRESSURE", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 0));
  SetVolumeOutputValue("RESIDUAL_VELOCITY-X", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 1));
  SetVolumeOutputValue("RESIDUAL_VELOCITY-Y", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("RESIDUAL_VELOCITY-Z", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 3));
    SetVolumeOutputValue("RESIDUAL_TEMPERATURE", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 4));
  } else {
    SetVolumeOutputValue("RESIDUAL_TEMPERATURE", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    SetVolumeOutputValue("RESIDUAL_TKE", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 0));
    SetVolumeOutputValue("RESIDUAL_OMEGA", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("RESIDUAL_NU_TILDE", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 0));
    break;
  case NONE:
    break;
  }
  
  SetVolumeOutputValue("LIMITER_PRESSURE", iPoint, Node_Flow->GetLimiter_Primitive(0));
  SetVolumeOutputValue("LIMITER_VELOCITY-X", iPoint, Node_Flow->GetLimiter_Primitive(1));
  SetVolumeOutputValue("LIMITER_VELOCITY-Y", iPoint, Node_Flow->GetLimiter_Primitive(2));
  if (nDim == 3){
    SetVolumeOutputValue("LIMITER_VELOCITY-Z", iPoint, Node_Flow->GetLimiter_Primitive(3));
    SetVolumeOutputValue("LIMITER_TEMPERATURE", iPoint, Node_Flow->GetLimiter_Primitive(4));
  } else {
    SetVolumeOutputValue("LIMITER_TEMPERATURE", iPoint, Node_Flow->GetLimiter_Primitive(3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    SetVolumeOutputValue("LIMITER_TKE", iPoint, Node_Turb->GetLimiter_Primitive(0));
    SetVolumeOutputValue("LIMITER_OMEGA", iPoint, Node_Turb->GetLimiter_Primitive(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("LIMITER_NU_TILDE", iPoint, Node_Turb->GetLimiter_Primitive(0));
    break;
  case NONE:
    break;
  }
  
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    SetVolumeOutputValue("DES_LENGTHSCALE", iPoint, Node_Flow->GetDES_LengthScale());
    SetVolumeOutputValue("WALL_DISTANCE", iPoint, Node_Geo->GetWall_Distance());
  }
  
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    SetVolumeOutputValue("ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation());
  }  
  
  if(config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    if (nDim == 3){
      SetVolumeOutputValue("VORTICITY_X", iPoint, Node_Flow->GetVorticity()[0]);
      SetVolumeOutputValue("VORTICITY_Y", iPoint, Node_Flow->GetVorticity()[1]);      
    } 
    SetVolumeOutputValue("VORTICITY_Z", iPoint, Node_Flow->GetVorticity()[2]);      
    SetVolumeOutputValue("Q_CRITERION", iPoint, GetQ_Criterion(config, geometry, Node_Flow));      
  }
}

void CIncFlowOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){
  
  if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver()  == RANS)) {  
    SetVolumeOutputValue("SKIN_FRICTION-X", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0));
    SetVolumeOutputValue("SKIN_FRICTION-Y", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1));
    if (nDim == 3)
      SetVolumeOutputValue("SKIN_FRICTION-Z", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2));
  
    SetVolumeOutputValue("HEAT_FLUX", iPoint, solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex));
    SetVolumeOutputValue("Y_PLUS", iPoint, solver[FLOW_SOL]->GetYPlus(iMarker, iVertex));
  }
}

su2double CIncFlowOutput::GetQ_Criterion(CConfig *config, CGeometry *geometry, CVariable* node_flow){
  
  unsigned short iDim, jDim;
  su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double Omega[3][3]    = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double Strain[3][3]   = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  for (iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0 ; jDim < nDim; jDim++) {
      Grad_Vel[iDim][jDim] = node_flow->GetGradient_Primitive(iDim+1, jDim);
      Strain[iDim][jDim]   = 0.5*(Grad_Vel[iDim][jDim] + Grad_Vel[jDim][iDim]);
      Omega[iDim][jDim]    = 0.5*(Grad_Vel[iDim][jDim] - Grad_Vel[jDim][iDim]);
    }
  }
  
  su2double OmegaMag = 0.0, StrainMag = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
      StrainMag += Strain[iDim][jDim]*Strain[iDim][jDim];
      OmegaMag  += Omega[iDim][jDim]*Omega[iDim][jDim];
    }
  }
  StrainMag = sqrt(StrainMag); OmegaMag = sqrt(OmegaMag);
  
  su2double Q = 0.5*(OmegaMag - StrainMag);
  
  return Q;
}

bool CIncFlowOutput::SetInit_Residuals(CConfig *config){
  
  return (config->GetUnsteady_Simulation() != STEADY && (config->GetIntIter() == 0))|| 
        (config->GetUnsteady_Simulation() == STEADY && (config->GetExtIter() < 2)); 
  
}

bool CIncFlowOutput::SetUpdate_Averages(CConfig *config, bool dualtime){
  
  return (config->GetUnsteady_Simulation() != STEADY && !dualtime);
      
}
