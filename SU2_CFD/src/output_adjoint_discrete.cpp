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
  
}

CDiscAdjFlowOutput::~CDiscAdjFlowOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }

}

void CDiscAdjFlowOutput::SetHistoryOutputFields(CConfig *config){

  // Iteration numbers
  AddHistoryOutput("INT_ITER",   "Int_Iter",  FORMAT_INTEGER, "INT_ITER");
  AddHistoryOutput("EXT_ITER",   "Ext_Iter",  FORMAT_INTEGER, "EXT_ITER");
  
  AddHistoryOutput("ADJOINT_DENSITY",    "Res[A_Rho]",  FORMAT_FIXED, "RESIDUALS");
  AddHistoryOutput("ADJOINT_MOMENTUM-X", "Res[A_RhoU]", FORMAT_FIXED, "RESIDUALS");
  AddHistoryOutput("ADJOINT_MOMENTUM-Y", "Res[A_RhoV]", FORMAT_FIXED, "RESIDUALS");
  AddHistoryOutput("ADJOINT_MOMENTUM-Z", "Res[A_RhoW]", FORMAT_FIXED, "RESIDUALS");
  AddHistoryOutput("ADJOINT_ENERGY",     "Res[A_E]",    FORMAT_FIXED, "RESIDUALS");
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    AddHistoryOutput("ADJOINT_NU_TILDE", "Res[A_nu]", FORMAT_FIXED, "RESIDUALS");
    break;  
  case SST:
    AddHistoryOutput("ADJOINT_KINETIC_ENERGY", "Res[A_k]", FORMAT_FIXED, "RESIDUALS");
    AddHistoryOutput("ADJOINT_DISSIPATION",    "Res[A_w]", FORMAT_FIXED, "RESIDUALS");
    break;
  default: break;
  }
  AddHistoryOutput("SENS_GEO",   "Sens_Geo",   FORMAT_SCIENTIFIC, "SENSITIVITIES");
  AddHistoryOutput("SENS_AOA",   "Sens_AoA",   FORMAT_SCIENTIFIC, "SENSITIVITIES");
  AddHistoryOutput("SENS_MACH",  "Sens_Mach",  FORMAT_SCIENTIFIC, "SENSITIVITIES");
  AddHistoryOutput("SENS_PRESS", "Sens_Press", FORMAT_SCIENTIFIC, "SENSITIVITIES");
  AddHistoryOutput("SENS_TEMP",  "Sens_Temp",  FORMAT_SCIENTIFIC, "SENSITIVITIES");
  
  AddHistoryOutput("PHYS_TIME",   "Time(min)",                FORMAT_SCIENTIFIC, "PHYS_TIME");
  
}

void CDiscAdjFlowOutput::LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) { 

  SetHistoryOutputField("INT_ITER", config[val_iZone]->GetIntIter());
  SetHistoryOutputField("EXT_ITER", config[val_iZone]->GetExtIter());
  
  SetHistoryOutputField("ADJOINT_DENSITY", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0)));
  SetHistoryOutputField("ADJOINT_MOMENTUM-X", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(1)));
  SetHistoryOutputField("ADJOINT_MOMENTUM-Y", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(2)));
  if (geometry[val_iZone][val_iInst][MESH_0]->GetnDim() == 3) {
    SetHistoryOutputField("ADJOINT_MOMENTUM-Z", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(3)));
    SetHistoryOutputField("ADJOINT_ENERGY", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(4)));
  } else {
    SetHistoryOutputField("ADJOINT_ENERGY", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(3)));    
  }
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputField("ADJOINT_NU_TILDE", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJTURB_SOL]->GetRes_RMS(0)));
    break;  
  case SST:
    SetHistoryOutputField("ADJOINT_KINETIC_ENERGY", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJTURB_SOL]->GetRes_RMS(0)));
    SetHistoryOutputField("ADJOINT_DISSIPATION",    log10(solver_container[val_iZone][val_iInst][MESH_0][ADJTURB_SOL]->GetRes_RMS(1)));
    break;
  default: break;
  }
  SetHistoryOutputField("SENS_GEO", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo());
  SetHistoryOutputField("SENS_AOA", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_AoA());
  SetHistoryOutputField("SENS_MACH", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Mach());
  SetHistoryOutputField("SENS_PRESS", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Press());
  SetHistoryOutputField("SENS_TEMP", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Temp());
  SetHistoryOutputField("PHYS_TIME", timeused);

}

void CDiscAdjFlowOutput::SetVolumeOutputFields(CConfig *config){
  
  
  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES");
  
  // Conservative
  AddVolumeOutput("ADJOINT_DENSITY",    "Adjoint_Density",    "CONSERVATIVE");
  AddVolumeOutput("ADJOINT_MOMENTUM-X", "Adjoint_Momentum_x", "CONSERVATIVE");
  AddVolumeOutput("ADJOINT_MOMENTUM-Y", "Adjoint_Momentum_y", "CONSERVATIVE");
  if (nDim == 3)
    AddVolumeOutput("ADJOINT_MOMENTUM-Z", "Adjoint_Momentum_z", "CONSERVATIVE");
  AddVolumeOutput("ADJOINT_ENERGY", "Adjoint_Energy", "CONSERVATIVE");
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    AddVolumeOutput("ADJOINT_NU_TILDE", "Adjoint_Nu_Tilde", "CONSERVATIVE");
    break;  
  case SST:
    AddVolumeOutput("ADJOINT_TKE", "Adjoint_TKE", "CONSERVATIVE");
    AddVolumeOutput("ADJOINT_NU_TILDE", "Adjoint_Omega", "CONSERVATIVE");
    break;
  default: break;
  }
  
  if (config->GetGrid_Movement()){
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY");
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY");
    if (nDim == 3)    
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY");
  }
  
  // Residuals
  AddVolumeOutput("RESIDUAL_ADJOINT_DENSITY",    "Residual_Adjoint_Density",    "RESIDUAL");
  AddVolumeOutput("RESIDUAL_ADJOINT_MOMENTUM-X", "Residual_Adjoint_Momentum_x", "RESIDUAL");
  AddVolumeOutput("RESIDUAL_ADJOINT_MOMENTUM-Y", "Residual_Adjoint_Momentum_y", "RESIDUAL");
  if (nDim == 3)
    AddVolumeOutput("RESIDUAL_ADJOINT_MOMENTUM-Z", "Residual_Adjoint_Momentum_z", "RESIDUAL");
  AddVolumeOutput("RESIDUAL_ADJOINT_ENERGY", "Residual_Adjoint_Energy", "RESIDUAL");
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    AddVolumeOutput("RESIDUAL_ADJOINT_NU_TILDE", "Residual_Adjoint_Nu_Tilde", "RESIDUAL");
    break;  
  case SST:
    AddVolumeOutput("RESIDUAL_ADJOINT_TKE", "Residual_Adjoint_TKE", "RESIDUAL");
    AddVolumeOutput("RESIDUAL_ADJOINT_NU_TILDE", "Residual_Adjoint_Omega", "RESIDUAL");
    break;
  default: break;
  }
  
  // Sensitivity
  AddVolumeOutput("SENSITIVITY_X", "Sensitivity_x", "SENSITIVITY");
  AddVolumeOutput("SENSITIVITY_Y", "Sensitivity_y", "SENSITIVITY");
  if (nDim == 3)
    AddVolumeOutput("SENSITIVITY_Z", "Sensitivity_z", "SENSITIVITY");  
  AddVolumeOutput("SENSITIVITY", "Surface_Sensitivity", "SENSITIVITY");
 
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
  
  SetVolumeOutputValue("ADJOINT_DENSITY",    iPoint, Node_AdjFlow->GetSolution(0));
  SetVolumeOutputValue("ADJOINT_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(1));
  SetVolumeOutputValue("ADJOINT_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(2));
  if (nDim == 3){
    SetVolumeOutputValue("ADJOINT_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(3));
    SetVolumeOutputValue("ADJOINT_ENERGY",     iPoint, Node_AdjFlow->GetSolution(4));
  } else {
    SetVolumeOutputValue("ADJOINT_ENERGY",     iPoint, Node_AdjFlow->GetSolution(3));    
  }
  
  // Turbulent 
  switch(turb_model){
  case SST:
    SetVolumeOutputValue("ADJOINT_ADJOINT_KINETIC_ENERGY", iPoint, Node_AdjTurb->GetSolution(0));
    SetVolumeOutputValue("ADJOINT_ADJOINT_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("ADJOINT_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(0));
    break;
  case NONE:
    break;
  }
  
  // Residuals
  SetVolumeOutputValue("RESIDUAL_ADJOINT_DENSITY", iPoint, Node_AdjFlow->GetSolution(0) - Node_AdjFlow->GetSolution_Old(0));
  SetVolumeOutputValue("RESIDUAL_ADJOINT_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(1) - Node_AdjFlow->GetSolution_Old(1));
  SetVolumeOutputValue("RESIDUAL_ADJOINT_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(2) - Node_AdjFlow->GetSolution_Old(2));
  if (nDim == 3){
    SetVolumeOutputValue("RESIDUAL_ADJOINT_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(3) - Node_AdjFlow->GetSolution_Old(3));
    SetVolumeOutputValue("RESIDUAL_ADJOINT_ENERGY", iPoint, Node_AdjFlow->GetSolution(4) - Node_AdjFlow->GetSolution_Old(4));
  } else {
    SetVolumeOutputValue("RESIDUAL_ADJOINT_ENERGY", iPoint, Node_AdjFlow->GetSolution(3) - Node_AdjFlow->GetSolution_Old(3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST:
    SetVolumeOutputValue("RESIDUAL_ADJOINT_KINETIC_ENERGY", iPoint, Node_AdjTurb->GetSolution(0) - Node_AdjTurb->GetSolution_Old(0));
    SetVolumeOutputValue("RESIDUAL_ADJOINT_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(1) - Node_AdjTurb->GetSolution_Old(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("RESIDUAL_ADJOINT_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(0) - Node_AdjTurb->GetSolution_Old(0));
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
  return true;
}

bool CDiscAdjFlowOutput::WriteScreen_Header(CConfig *config) {return true;  }

bool CDiscAdjFlowOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {return true;  }

