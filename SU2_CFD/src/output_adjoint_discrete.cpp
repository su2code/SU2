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

CDiscAdjFlowOutput::~CDiscAdjFlowOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }

}

void CDiscAdjFlowOutput::SetOutputFields(CConfig *config){

  // Iteration numbers
  AddOutputField("INT_ITER",   "Int_Iter",  FORMAT_INTEGER, "INT_ITER");
  AddOutputField("EXT_ITER",   "Ext_Iter",  FORMAT_INTEGER, "EXT_ITER");
  
  AddOutputField("ADJOINT_DENSITY",    "Res[A_Rho]",  FORMAT_FIXED, "RESIDUALS");
  AddOutputField("ADJOINT_MOMENTUM-X", "Res[A_RhoU]", FORMAT_FIXED, "RESIDUALS");
  AddOutputField("ADJOINT_MOMENTUM-Y", "Res[A_RhoV]", FORMAT_FIXED, "RESIDUALS");
  AddOutputField("ADJOINT_MOMENTUM-Z", "Res[A_RhoW]", FORMAT_FIXED, "RESIDUALS");
  AddOutputField("ADJOINT_ENERGY",     "Res[A_E]",    FORMAT_FIXED, "RESIDUALS");
  
  AddOutputField("SENS_GEO",   "Sens_Geo",   FORMAT_SCIENTIFIC, "SENSITIVITIES");
  AddOutputField("SENS_AOA",   "Sens_AoA",   FORMAT_SCIENTIFIC, "SENSITIVITIES");
  AddOutputField("SENS_MACH",  "Sens_Mach",  FORMAT_SCIENTIFIC, "SENSITIVITIES");
  AddOutputField("SENS_PRESS", "Sens_Press", FORMAT_SCIENTIFIC, "SENSITIVITIES");
  AddOutputField("SENS_TEMP",  "Sens_Temp",  FORMAT_SCIENTIFIC, "SENSITIVITIES");
  
  AddOutputField("PHYS_TIME",   "Time(min)",                FORMAT_SCIENTIFIC, "PHYS_TIME");
  
}

inline void CDiscAdjFlowOutput::LoadOutput_Data(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) { 

  SetOutputFieldValue("INT_ITER", config[val_iZone]->GetIntIter());
  SetOutputFieldValue("EXT_ITER", config[val_iZone]->GetExtIter());
  
  SetOutputFieldValue("ADJOINT_DENSITY", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0)));
  SetOutputFieldValue("ADJOINT_MOMENTUM-X", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(1)));
  SetOutputFieldValue("ADJOINT_MOMENTUM-Y", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(2)));
  if (geometry[val_iZone][val_iInst][MESH_0]->GetnDim() == 3) {
    SetOutputFieldValue("ADJOINT_MOMENTUM-Z", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(3)));
    SetOutputFieldValue("ADJOINT_ENERGY", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(4)));
  } else {
    SetOutputFieldValue("ADJOINT_ENERGY", log10(solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetRes_RMS(3)));    
  }
  SetOutputFieldValue("SENS_GEO", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo());
  SetOutputFieldValue("SENS_AOA", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_AoA());
  SetOutputFieldValue("SENS_MACH", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Mach());
  SetOutputFieldValue("SENS_PRESS", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Press());
  SetOutputFieldValue("SENS_TEMP", solver_container[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Temp());
  SetOutputFieldValue("PHYS_TIME", timeused);

}

inline bool CDiscAdjFlowOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { return true; }

inline bool CDiscAdjFlowOutput::WriteScreen_Header(CConfig *config) {return true;  }

inline bool CDiscAdjFlowOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {return true;  }

