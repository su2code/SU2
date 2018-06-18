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

  char buffer[50];

  // Retrieve the history filename
  string history_filename = config->GetConv_FileName();

  // Append the zone ID
  if(config->GetnZone() > 1){
    history_filename = config->GetMultizone_HistoryFileName(history_filename, val_iZone);
  }
  strcpy (char_histfile, history_filename.data());

  // Append the restart iteration: if dynamic problem and restart
  if (config->GetWrt_Unsteady() && config->GetRestart()) {
    long iExtIter = config->GetDyn_RestartIter();
    if (SU2_TYPE::Int(iExtIter) < 10) SPRINTF (buffer, "_0000%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d", SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d", SU2_TYPE::Int(iExtIter));
    strcat(char_histfile, buffer);
  }

  // Add the correct file extension depending on the file format
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == FIELDVIEW)) SPRINTF (buffer, ".dat");
  else if ((config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
           (config->GetOutput_FileFormat() == FIELDVIEW_BINARY))  SPRINTF (buffer, ".plt");
  else if (config->GetOutput_FileFormat() == PARAVIEW)  SPRINTF (buffer, ".csv");
  strcat(char_histfile, buffer);
  SetOutputFields(config);
  
  // Open the history file using only the master node
  if (rank == MASTER_NODE){
    cout << "History filename: " << char_histfile << endl;
    HistFile.open(char_histfile, ios::out);
    HistFile.precision(15);
    SetHistoryFile_Header(config);
  }

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

inline void CHeatOutput::LoadOutput_Data(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) {
  
  SetOutputFieldValue("EXT_ITER",     config[val_iZone]->GetExtIter());
  SetOutputFieldValue("INT_ITER",     config[val_iZone]->GetIntIter());
  
  SetOutputFieldValue("HEATFLUX",     solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetTotal_HeatFlux());
  SetOutputFieldValue("HEATFLUX_MAX", solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetTotal_MaxHeatFlux());
  SetOutputFieldValue("TEMPERATURE",  solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetTotal_AvgTemperature());
  SetOutputFieldValue("HEAT", log10(solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetRes_RMS(0)));
  
  SetOutputFieldValue("TIME", timeused);
  SetOutputFieldValue("LINSOL_ITER", solver_container[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetIterLinSolver());
  
}
  

inline void CHeatOutput::SetOutputFields(CConfig *config){
  
  
  AddOutputField("EXT_ITER", "Ext_Iter", FORMAT_INTEGER, "EXT_ITER");
  AddOutputField("INT_ITER", "Int_Iter", FORMAT_INTEGER, "INT_ITER");
  
  AddOutputField("PHYS_TIME", "Time(min)", FORMAT_SCIENTIFIC, "PHYS_TIME");
  AddOutputField("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL_ITER");
  
  AddOutputField("HEATFLUX", "HF(Total)",      FORMAT_SCIENTIFIC, "HEAT");
  AddOutputField("HEATFLUX_MAX", "HF(Max)",    FORMAT_SCIENTIFIC, "HEAT");
  AddOutputField("TEMPERATURE", "Temp(Total)", FORMAT_SCIENTIFIC, "HEAT");
  
  AddOutputField("HEAT", "Res[Heat]", FORMAT_FIXED, "RESIDUALS");
  
}
