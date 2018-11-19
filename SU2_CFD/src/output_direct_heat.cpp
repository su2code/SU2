/*!
 * \file output_direct_heat.cpp
 * \brief Main subroutines for the heat solver output
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

CHeatOutput::CHeatOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {

  nDim = geometry->GetnDim();


  /*--- Set the default history fields if nothing is set in the config file ---*/
  
  if (nRequestedHistoryFields == 0){
    RequestedHistoryFields.push_back("EXT_ITER");
    RequestedHistoryFields.push_back("RMS_RES");
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    RequestedScreenFields.push_back("EXT_ITER");    
    RequestedScreenFields.push_back("RMS_TEMPERATURE");
    nRequestedScreenFields = RequestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    RequestedVolumeFields.push_back("COORDINATES");
    RequestedVolumeFields.push_back("CONSERVATIVE");
    RequestedVolumeFields.push_back("PRIMITIVE");
    nRequestedVolumeFields = RequestedVolumeFields.size();
  }
}

CHeatOutput::~CHeatOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }

}


bool CHeatOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { 
 if (!write_dualtime){
  return true;
 }
 else {
   return false;
 }
}

bool CHeatOutput::WriteScreen_Header(CConfig *config) {  
  bool write_header = false;
  if (config->GetUnsteady_Simulation() == STEADY || config->GetUnsteady_Simulation() == TIME_STEPPING) {
    write_header = ((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0);
  } else {
    write_header = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND) && config->GetIntIter() == 0;
  }
  return write_header;
}

bool CHeatOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {
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

void CHeatOutput::LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) {
  
  CSolver* heat_solver = solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL];  
  
  SetHistoryOutputValue("EXT_ITER",     config[val_iZone]->GetExtIter());
  SetHistoryOutputValue("INT_ITER",     config[val_iZone]->GetIntIter());
  
  SetHistoryOutputValue("HEATFLUX",     heat_solver->GetTotal_HeatFlux());
  SetHistoryOutputValue("HEATFLUX_MAX", heat_solver->GetTotal_MaxHeatFlux());
  SetHistoryOutputValue("AVG_TEMPERATURE",  heat_solver->GetTotal_AvgTemperature());
  SetHistoryOutputValue("RMS_TEMPERATURE", log10(heat_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("MAX_TEMPERATURE", log10(heat_solver->GetRes_Max(0)));
  
  SetHistoryOutputValue("PHYS_TIME", timeused);
  SetHistoryOutputValue("LINSOL_ITER", heat_solver->GetIterLinSolver());
  
}
  

void CHeatOutput::SetHistoryOutputFields(CConfig *config){
  
  AddHistoryOutput("EXT_ITER", "Ext_Iter", FORMAT_INTEGER, "EXT_ITER");
  AddHistoryOutput("INT_ITER", "Int_Iter", FORMAT_INTEGER, "INT_ITER");
  
  AddHistoryOutput("PHYS_TIME",   "Time(min)",                FORMAT_SCIENTIFIC, "PHYS_TIME");
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL_ITER");
  
  AddHistoryOutput("HEATFLUX", "HF",      FORMAT_SCIENTIFIC, "HEAT", TYPE_COEFFICIENT);
  AddHistoryOutput("HEATFLUX_MAX", "MaxHF",    FORMAT_SCIENTIFIC, "HEAT", TYPE_COEFFICIENT);
  AddHistoryOutput("AVG_TEMPERATURE", "AvgTemp", FORMAT_SCIENTIFIC, "HEAT", TYPE_COEFFICIENT);
  
  AddHistoryOutput("RMS_TEMPERATURE", "rms[T]", FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL);
  AddHistoryOutput("MAX_TEMPERATURE", "max[T]", FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);
}


void CHeatOutput::SetVolumeOutputFields(CConfig *config){
  
  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES");
  
  // Conservative
  AddVolumeOutput("TEMPERATURE", "Temperature", "CONSERVATIVE");

  // Residuals  
  AddVolumeOutput("RESIDUAL_TEMPERATURE", "Residual_Temperature", "RESIDUAL");
  
}


void CHeatOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  
  CVariable* Node_Heat = solver[HEAT_SOL]->node[iPoint]; 
  
  CPoint*    Node_Geo  = geometry->node[iPoint];
  
  // Grid coordinates
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
 
  // Conservative
  SetVolumeOutputValue("TEMPEATURE", iPoint, Node_Heat->GetSolution(0));
  
  // Residuals    
  SetVolumeOutputValue("RESIDUAL_TEMPERATURE", iPoint, solver[HEAT_SOL]->LinSysRes.GetBlock(iPoint, 0));
  
}

