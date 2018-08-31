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

}

CHeatOutput::~CHeatOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }

}


inline bool CHeatOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) {
  return true;
}

inline bool CHeatOutput::WriteScreen_Header(CConfig *config) { 

//  return (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));
  return true;
}

inline bool CHeatOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {
  return true;
}

void CHeatOutput::LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) {
  
  CSolver* heat_solver = solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL];  
  
  SetHistoryOutputField("EXT_ITER",     config[val_iZone]->GetExtIter());
  SetHistoryOutputField("INT_ITER",     config[val_iZone]->GetIntIter());
  
  SetHistoryOutputField("HEATFLUX",     heat_solver->GetTotal_HeatFlux());
  SetHistoryOutputField("HEATFLUX_MAX", heat_solver->GetTotal_MaxHeatFlux());
  SetHistoryOutputField("TEMPERATURE",  heat_solver->GetTotal_AvgTemperature());
  SetHistoryOutputField("HEAT", log10(heat_solver->GetRes_RMS(0)));
  
  SetHistoryOutputField("PHYS_TIME", timeused);
  SetHistoryOutputField("LINSOL_ITER", heat_solver->GetIterLinSolver());
  
}
  

void CHeatOutput::SetHistoryOutputFields(CConfig *config){
  
  AddHistoryOutput("EXT_ITER", "Ext_Iter", FORMAT_INTEGER, "EXT_ITER");
  AddHistoryOutput("INT_ITER", "Int_Iter", FORMAT_INTEGER, "INT_ITER");
  
  AddHistoryOutput("PHYS_TIME",   "Time(min)",                FORMAT_SCIENTIFIC, "PHYS_TIME");
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL_ITER");
  
  AddHistoryOutput("HEATFLUX", "HF(Total)",      FORMAT_SCIENTIFIC, "HEAT");
  AddHistoryOutput("HEATFLUX_MAX", "HF(Max)",    FORMAT_SCIENTIFIC, "HEAT");
  AddHistoryOutput("TEMPERATURE", "Temp(Total)", FORMAT_SCIENTIFIC, "HEAT");
  
  AddHistoryOutput("HEAT", "Res[T]", FORMAT_FIXED, "RESIDUALS");
  
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

