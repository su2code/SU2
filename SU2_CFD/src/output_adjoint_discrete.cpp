/*!
 * \file output_adjoint_mean.cpp
 * \brief Main subroutines for flow discrete adjoint output
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

CDiscAdjFlowOutput::CDiscAdjFlowOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {
  
  nDim = geometry->GetnDim();
 
  turb_model = config->GetKind_Turb_Model();
  
  /*--- Set the default history fields if nothing is set in the config file ---*/
  
  if (nRequestedHistoryFields == 0){
    RequestedHistoryFields.push_back("EXT_ITER");
    RequestedHistoryFields.push_back("RMS_RES");
    RequestedHistoryFields.push_back("SENSITIVITIES");
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }
  
  if (nRequestedScreenFields == 0){
    RequestedScreenFields.push_back("EXT_ITER");
    RequestedScreenFields.push_back("RMS_ADJ_DENSITY");
    RequestedScreenFields.push_back("RMS_ADJ_MOMENTUM-X");
    RequestedScreenFields.push_back("SENS_GEO");
    RequestedScreenFields.push_back("SENS_AOA");
    nRequestedScreenFields = RequestedScreenFields.size();
  }
  
  if (nRequestedVolumeFields == 0){
    RequestedVolumeFields.push_back("COORDINATES");
    RequestedVolumeFields.push_back("CONSERVATIVE");    
    RequestedVolumeFields.push_back("SENSITIVITIES");
    nRequestedVolumeFields = RequestedVolumeFields.size();
  }
  
  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Adj. Comp. Fluid)";
  MultiZoneHeaderString = ss.str();
  
}

CDiscAdjFlowOutput::~CDiscAdjFlowOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }

}

void CDiscAdjFlowOutput::SetHistoryOutputFields(CConfig *config){

  /// BEGIN_GROUP: ITERATION, DESCRIPTION: Iteration identifier.
  /// DESCRIPTION: The time iteration index.
  AddHistoryOutput("TIME_ITER",     "Time_Iter",  FORMAT_INTEGER, "ITER"); 
  /// DESCRIPTION: The internal iteration index.
  AddHistoryOutput("OUTER_ITER",   "Outer_Iter",  FORMAT_INTEGER, "ITER"); 
  /// DESCRIPTION: The external iteration index.
  AddHistoryOutput("INNER_ITER",   "Inner_Iter", FORMAT_INTEGER,  "ITER"); 
  /// END_GROUP
  
  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the conservative variables. 
  /// DESCRIPTION: Root-mean square residual of the adjoint density.
  AddHistoryOutput("RMS_ADJ_DENSITY",    "rms[A_Rho]",  FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL); 
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum x-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-X", "rms[A_RhoU]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum y-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-Y", "rms[A_RhoV]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum z-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-Z", "rms[A_RhoW]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint energy.
  AddHistoryOutput("RMS_ADJ_ENERGY",     "rms[A_E]",    FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL); 
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Root-mean square residual of the adjoint nu tilde.
    AddHistoryOutput("RMS_ADJ_NU_TILDE", "rms[A_nu]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);      
    break;  
  case SST:
    /// DESCRIPTION: Root-mean square residual of the adjoint kinetic energy.
    AddHistoryOutput("RMS_ADJ_KINETIC_ENERGY", "rms[A_k]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL); 
    /// DESCRIPTION: Root-mean square residual of the adjoint dissipation.
    AddHistoryOutput("RMS_ADJ_DISSIPATION",    "rms[A_w]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);   
    break;
  default: break;
  }
  /// END_GROUP
  
  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the conservative variables. 
  /// DESCRIPTION: Maximum residual of the adjoint density.
  AddHistoryOutput("MAX_ADJ_DENSITY",    "max[A_Rho]",  FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint momentum x-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-X", "max[A_RhoU]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL); 
  /// DESCRIPTION: Maximum residual of the adjoint momentum y-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-Y", "max[A_RhoV]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL); 
  /// DESCRIPTION: Maximum residual of the adjoint momentum z-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-Z", "max[A_RhoW]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL); 
  /// DESCRIPTION: Maximum residual of the adjoint energy.
  AddHistoryOutput("MAX_ADJ_ENERGY",     "max[A_E]",    FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL); 
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Maximum residual of the adjoint nu tilde.
    AddHistoryOutput("MAX_ADJ_NU_TILDE", "max[A_nu]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);      
    break;  
  case SST:
    /// DESCRIPTION: Maximum residual of the adjoint kinetic energy.
    AddHistoryOutput("MAX_ADJ_KINETIC_ENERGY", "max[A_k]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);   
    /// DESCRIPTION: Maximum residual of the adjoint dissipation.
    AddHistoryOutput("MAX_ADJ_DISSIPATION",    "max[A_w]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL); 
    break;
  default: break;
  }
  /// END_GROUP
  
  /// BEGIN_GROUP: SENSITIVITIES, DESCRIPTION: Sensitivities of different geometrical or boundary values.   
  /// DESCRIPTION: Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.
  AddHistoryOutput("SENS_GEO",   "Sens_Geo",   FORMAT_SCIENTIFIC, "SENSITIVITIES", TYPE_COEFFICIENT); 
  /// DESCRIPTION: Sensitivity of the objective function with respect to the angle of attack (only for compressible solver).
  AddHistoryOutput("SENS_AOA",   "Sens_AoA",   FORMAT_SCIENTIFIC, "SENSITIVITIES", TYPE_COEFFICIENT); 
  /// DESCRIPTION: Sensitivity of the objective function with respect to the Mach number (only of compressible solver).
  AddHistoryOutput("SENS_MACH",  "Sens_Mach",  FORMAT_SCIENTIFIC, "SENSITIVITIES", TYPE_COEFFICIENT); 
  /// DESCRIPTION: Sensitivity of the objective function with respect to the far-field pressure.
  AddHistoryOutput("SENS_PRESS", "Sens_Press", FORMAT_SCIENTIFIC, "SENSITIVITIES", TYPE_COEFFICIENT); 
  /// DESCRIPTION: Sensitivity of the objective function with respect to the far-field temperature.
  AddHistoryOutput("SENS_TEMP",  "Sens_Temp",  FORMAT_SCIENTIFIC, "SENSITIVITIES", TYPE_COEFFICIENT); 
  /// END_GROUP
  
  /// DESCRIPTION: Currently used wall-clock time.
  AddHistoryOutput("PHYS_TIME",   "Time(min)", FORMAT_SCIENTIFIC, "PHYS_TIME"); 
  
}

void CDiscAdjFlowOutput::LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) { 
 
  CSolver* adjflow_solver = solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL];
  CSolver* adjturb_solver = solver_container[val_iZone][val_iInst][MESH_0][ADJTURB_SOL];  

  SetHistoryOutputValue("INT_ITER", config[val_iZone]->GetIntIter());
  SetHistoryOutputValue("EXT_ITER", config[val_iZone]->GetExtIter());
  
  SetHistoryOutputValue("RMS_ADJ_DENSITY", log10(adjflow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_RMS(2)));
  if (geometry[val_iZone][val_iInst][MESH_0]->GetnDim() == 3) {
    SetHistoryOutputValue("RMS_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_RMS(3)));
    SetHistoryOutputValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(4)));
  } else {
    SetHistoryOutputValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(3)));    
  }
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("RMS_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_RMS(0)));
    break;  
  case SST:
    SetHistoryOutputValue("RMS_ADJ_KINETIC_ENERGY", log10(adjturb_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_ADJ_DISSIPATION",    log10(adjturb_solver->GetRes_RMS(1)));
    break;
  default: break;
  }
  SetHistoryOutputValue("MAX_ADJ_DENSITY", log10(adjflow_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_Max(2)));
  if (geometry[val_iZone][val_iInst][MESH_0]->GetnDim() == 3) {
    SetHistoryOutputValue("MAX_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_Max(3)));
    SetHistoryOutputValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(4)));
  } else {
    SetHistoryOutputValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(3)));    
  }
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("MAX_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_Max(0)));
    break;  
  case SST:
    SetHistoryOutputValue("MAX_ADJ_KINETIC_ENERGY", log10(adjturb_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_ADJOINT_DISSIPATION",    log10(adjturb_solver->GetRes_Max(1)));
    break;
  default: break;
  }
  SetHistoryOutputValue("SENS_GEO", adjflow_solver->GetTotal_Sens_Geo());
  SetHistoryOutputValue("SENS_AOA", adjflow_solver->GetTotal_Sens_AoA());
  SetHistoryOutputValue("SENS_MACH", adjflow_solver->GetTotal_Sens_Mach());
  SetHistoryOutputValue("SENS_PRESS", adjflow_solver->GetTotal_Sens_Press());
  SetHistoryOutputValue("SENS_TEMP", adjflow_solver->GetTotal_Sens_Temp());
  SetHistoryOutputValue("PHYS_TIME", timeused);

}

void CDiscAdjFlowOutput::SetVolumeOutputFields(CConfig *config){
  
  /// BEGIN_GROUP: COORDINATES, DESCRIPTION: Coordinates of the mesh nodes.
  /// DESCRIPTION: x coordinates of the mesh nodes.
  AddVolumeOutput("COORD-X", "x", "COORDINATES"); 
  /// DESCRIPTION: y coordinates of the mesh nodes.
  AddVolumeOutput("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    /// DESCRIPTION: z coordinates of the mesh nodes.
    AddVolumeOutput("COORD-Z", "z", "COORDINATES");
  /// END_GROUP
  
  /// BEGIN_GROUP: CONSERVATIVE, DESCRIPTION: The conservative variables of the adjoint solver.
  /// DESCRIPTION: Adjoint density.
  AddVolumeOutput("ADJ_DENSITY",    "Adjoint_Density",    "CONSERVATIVE"); 
  /// DESCRIPTION: Adjoint momentum x-component.
  AddVolumeOutput("ADJ_MOMENTUM-X", "Adjoint_Momentum_x", "CONSERVATIVE"); 
  /// DESCRIPTION: Adjoint momentum y-component.
  AddVolumeOutput("ADJ_MOMENTUM-Y", "Adjoint_Momentum_y", "CONSERVATIVE"); 
  if (nDim == 3)
    /// DESCRIPTION: Adjoint momentum z-component.
    AddVolumeOutput("ADJ_MOMENTUM-Z", "Adjoint_Momentum_z", "CONSERVATIVE"); 
  /// DESCRIPTION: Adjoint energy.
  AddVolumeOutput("ADJ_ENERGY", "Adjoint_Energy", "CONSERVATIVE");           
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Adjoint nu tilde.
    AddVolumeOutput("ADJ_NU_TILDE", "Adjoint_Nu_Tilde", "CONSERVATIVE"); 
    break;  
  case SST:
    /// DESCRIPTION: Adjoint kinetic energy.
    AddVolumeOutput("ADJ_KINETIC_ENERGY", "Adjoint_TKE", "CONSERVATIVE"); 
    /// DESCRIPTION: Adjoint dissipation.
    AddVolumeOutput("ADJ_DISSIPATION", "Adjoint_Omega", "CONSERVATIVE");  
    break;
  default: break;
  }
  /// END_GROUP
  
  /// BEGIN_GROUP: GRID_VELOCITY, DESCRIPTION: The grid velocity in case of a moving grid.  
  if (config->GetGrid_Movement()){
    /// DESCRIPTION: Grid velocity x-component.
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY"); 
    /// DESCRIPTION: Grid velocity y-component.
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY"); 
    if (nDim == 3)    
      /// DESCRIPTION: Grid velocity z-component.
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY"); 
  }
  /// END_GROUP
  
  /// BEGIN_GROUP: RESIDUAL, DESCRIPTION: Residuals of the conservative variables. 
  /// DESCRIPTION: Residual of the adjoint density.
  AddVolumeOutput("RES_ADJ_DENSITY",    "Residual_Adjoint_Density",    "RESIDUAL");  
  /// DESCRIPTION: Residual of the adjoint momentum x-component.
  AddVolumeOutput("RES_ADJ_MOMENTUM-X", "Residual_Adjoint_Momentum_x", "RESIDUAL");  
  /// DESCRIPTION: Residual of the adjoint momentum y-component.
  AddVolumeOutput("RES_ADJ_MOMENTUM-Y", "Residual_Adjoint_Momentum_y", "RESIDUAL");  
  if (nDim == 3)
    /// DESCRIPTION: Residual of the adjoint momentum z-component.
    AddVolumeOutput("RES_ADJ_MOMENTUM-Z", "Residual_Adjoint_Momentum_z", "RESIDUAL"); 
  /// DESCRIPTION: Residual of the adjoint energy. 
  AddVolumeOutput("RES_ADJ_ENERGY", "Residual_Adjoint_Energy", "RESIDUAL");            
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Residual of the nu tilde. 
    AddVolumeOutput("RES_ADJ_NU_TILDE", "Residual_Adjoint_Nu_Tilde", "RESIDUAL"); 
    break;  
  case SST:
    /// DESCRIPTION: Residual of the adjoint kinetic energy. 
    AddVolumeOutput("RES_ADJ_TKE", "Residual_Adjoint_TKE", "RESIDUAL");    
    /// DESCRIPTION: Residual of the adjoint dissipation.
    AddVolumeOutput("RES_ADJ_NU_TILDE", "Residual_Adjoint_Omega", "RESIDUAL");    
    break;
  default: break;
  }
  /// END_GROUP
  
  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Geometrical sensitivities of the current objective function.
  /// DESCRIPTION: Sensitivity x-component.
  AddVolumeOutput("SENSITIVITY_X", "Sensitivity_x", "SENSITIVITY"); 
  /// DESCRIPTION: Sensitivity y-component.
  AddVolumeOutput("SENSITIVITY_Y", "Sensitivity_y", "SENSITIVITY");
  if (nDim == 3)
    /// DESCRIPTION: Sensitivity z-component.
    AddVolumeOutput("SENSITIVITY_Z", "Sensitivity_z", "SENSITIVITY");   
  /// DESCRIPTION: Sensitivity in normal direction.
  AddVolumeOutput("SENSITIVITY", "Surface_Sensitivity", "SENSITIVITY"); 
  /// END_GROUP
 
}

void CDiscAdjFlowOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  
  CVariable* Node_AdjFlow = solver[ADJFLOW_SOL]->node[iPoint]; 
  CVariable* Node_AdjTurb = NULL;
  CPoint*    Node_Geo     = geometry->node[iPoint];
  
  if (config->GetKind_Turb_Model() != NONE){
    Node_AdjTurb = solver[ADJTURB_SOL]->node[iPoint]; 
  }
  
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
  
  SetVolumeOutputValue("ADJ_DENSITY",    iPoint, Node_AdjFlow->GetSolution(0));
  SetVolumeOutputValue("ADJ_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(1));
  SetVolumeOutputValue("ADJ_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(2));
  if (nDim == 3){
    SetVolumeOutputValue("ADJ_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(3));
    SetVolumeOutputValue("ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(4));
  } else {
    SetVolumeOutputValue("ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(3));    
  }
  
  // Turbulent 
  switch(turb_model){
  case SST:
    SetVolumeOutputValue("ADJ_KINETIC_ENERGY", iPoint, Node_AdjTurb->GetSolution(0));
    SetVolumeOutputValue("ADJ_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(0));
    break;
  case NONE:
    break;
  }
  
  // Residuals
  SetVolumeOutputValue("RES_ADJ_DENSITY", iPoint, Node_AdjFlow->GetSolution(0) - Node_AdjFlow->GetSolution_Old(0));
  SetVolumeOutputValue("RES_ADJ_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(1) - Node_AdjFlow->GetSolution_Old(1));
  SetVolumeOutputValue("RES_ADJ_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(2) - Node_AdjFlow->GetSolution_Old(2));
  if (nDim == 3){
    SetVolumeOutputValue("RES_ADJ_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(3) - Node_AdjFlow->GetSolution_Old(3));
    SetVolumeOutputValue("RES_ADJ_ENERGY", iPoint, Node_AdjFlow->GetSolution(4) - Node_AdjFlow->GetSolution_Old(4));
  } else {
    SetVolumeOutputValue("RES_ADJ_ENERGY", iPoint, Node_AdjFlow->GetSolution(3) - Node_AdjFlow->GetSolution_Old(3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    SetVolumeOutputValue("RES_ADJ_KINETIC_ENERGY", iPoint, Node_AdjTurb->GetSolution(0) - Node_AdjTurb->GetSolution_Old(0));
    SetVolumeOutputValue("RES_ADJ_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(1) - Node_AdjTurb->GetSolution_Old(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("RES_ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(0) - Node_AdjTurb->GetSolution_Old(0));
    break;
  case NONE:
    break;
  }
  
  SetVolumeOutputValue("SENSITIVITY_X", iPoint, Node_AdjFlow->GetSensitivity(0));
  SetVolumeOutputValue("SENSITIVITY_Y", iPoint, Node_AdjFlow->GetSensitivity(1));
  if (nDim == 3)
    SetVolumeOutputValue("SENSITIVITY_Z", iPoint, Node_AdjFlow->GetSensitivity(2));
  
}

void CDiscAdjFlowOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){
  
  SetVolumeOutputValue("SENSITIVITY", iPoint, solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex));
  
}

bool CDiscAdjFlowOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { 
 if (!write_dualtime){
   return true;
 }
 else {
   return false;
 }
}

bool CDiscAdjFlowOutput::WriteScreen_Header(CConfig *config) {  
  bool write_header = false;
  if (config->GetUnsteady_Simulation() == STEADY || config->GetUnsteady_Simulation() == TIME_STEPPING) {
    write_header = (config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0;
  } else {
    write_header = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND) && config->GetIntIter() == 0;
  }
  return write_header;
}

bool CDiscAdjFlowOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {
  bool write_output = false;
  
  if (((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) ) 
      && write_dualtime ){
    write_output = (config->GetIntIter() % config->GetWrt_Con_Freq_DualTime() == 0);
  }
  else if (((config->GetUnsteady_Simulation() == STEADY) || (config->GetUnsteady_Simulation() == TIME_STEPPING) )){
    write_output = (config->GetExtIter() % config->GetWrt_Con_Freq() == 0) ;    
  } 
  return write_output;
}

bool CDiscAdjFlowOutput::SetInit_Residuals(CConfig *config){
  
  return (config->GetUnsteady_Simulation() != STEADY && (config->GetIntIter() == 0))|| 
        (config->GetUnsteady_Simulation() == STEADY && (config->GetExtIter() < 2)); 
  
}

bool CDiscAdjFlowOutput::SetUpdate_Averages(CConfig *config, bool dualtime){
  
  return (config->GetUnsteady_Simulation() != STEADY && !dualtime);
      
}

