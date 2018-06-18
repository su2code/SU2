/*!
 * \file output_direct_mean.cpp
 * \brief Main subroutines for compressible flow output
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

CFlowOutput::CFlowOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {

  // Open the history file using only the master node
  if (rank == MASTER_NODE){

    nDim = geometry->GetnDim();

    bool turbulent = config->GetKind_Turb_Model() != NONE;
    
    char buffer[50], char_histfile[200];

    // Retrieve the history filename
    string history_filename = config->GetConv_FileName();

    // Append the zone ID
    if(config->GetnZone() > 1){
      history_filename = config->GetMultizone_HistoryFileName(history_filename, val_iZone);
    }
    strcpy (char_histfile, history_filename.data());

    // Append the restart iteration: if dynamic problem and restart
    if (config->GetWrt_Dynamic() && config->GetRestart()) {
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

    cout << "History filename: " << char_histfile << endl;
    HistFile.open(char_histfile, ios::out);
    HistFile.precision(15);
    SetHistoryFile_Header(config);

    /*--- Initialize number of variables ---*/
    nVar = nDim + 2;
 
    ResRMS = new su2double[nVar];
    ResMax = new su2double[nVar];
    
    LinSolvIter = 0;
    Time_Used = 0.0;

    iExtIter = 0;
    iIntIter = 0;
    
  }

}

CFlowOutput::~CFlowOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
    
    delete [] ResRMS;
    delete [] ResMax;
  }


}



inline void CFlowOutput::SetOutputFields(CConfig *config){
  
  // Iteration numbers
  AddOutputField("INT_ITER",   "Int_Iter",  FORMAT_INTEGER, "INT_ITER");
  AddOutputField("EXT_ITER",   "Ext_Iter",  FORMAT_INTEGER, "EXT_ITER");
  
  // Residuals
  AddOutputField("DENSITY",    "Res[Rho]",  FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("MOMENTUM-X", "Res[RhoU]", FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("MOMENTUM-Y", "Res[RhoV]", FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("MOMENTUM-Z", "Res[RhoW]", FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("ENERGY",     "Res[RhoE]", FORMAT_FIXED,   "RESIDUALS");
  
  // Aerodynamic coefficients
  AddOutputField("DRAG",       "CD(Total)", FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("LIFT",       "CL(Total)", FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("SIDEFORCE",  "CSF(Total)",FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("MOMENT-X",  "CMx(Total)", FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("MOMENT-Y",  "CMy(Total)", FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("MOMENT-Z",  "CMz(Total)", FORMAT_SCIENTIFIC, "AERO_COEFF");
  AddOutputField("FORCE-X",  "CFx(Total)", FORMAT_SCIENTIFIC,  "AERO_COEFF");
  AddOutputField("FORCE-Y",  "CFy(Total)", FORMAT_SCIENTIFIC,  "AERO_COEFF");
  AddOutputField("FORCE-Z",  "CFz(Total)", FORMAT_SCIENTIFIC,  "AERO_COEFF");
  AddOutputField("EFFICIENCY", "CEff(Total)", FORMAT_SCIENTIFIC, "AERO_COEFF");
  
  // Misc.
  AddOutputField("AOA",      "AoA", FORMAT_SCIENTIFIC, "AOA");
  AddOutputField("PHYS_TIME", "Time(min)", FORMAT_SCIENTIFIC, "PHYS_TIME");
  AddOutputField("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL_ITER");
  
  // Surface output
  AddOutputField("AVG_MASSFLOW", "Avg_Massflow", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("AVG_MACH", "Avg_Mach", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("AVG_TEMP", "Avg_Temp", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("AVG_PRESS","Avg_Press", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("AVG_DENSITY","Avg_Density", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("AVG_ENTHALPY", "Avg_Enthalpy", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("AVG_NORMALVEL", "Avg_NormalVel", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("UNIFORMITY", "Uniformity", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("SECONDARY_STRENGTH", "Secondary_Strength", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("MOMENTUM_DISTORTION", "Momentum_Distortion", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("SECONDARY_OVER_UNFORMITY", "Secondary_Over_Uniformity", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("AVG_TOTALTEMP", "Avg_TotalTemp", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("AVG_TOTALPRESS", "Avg_TotalPress", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
  AddOutputField("PRESSURE_DROP", "Pressure_Drop", FORMAT_SCIENTIFIC, "SURFACE_OUTPUT");
 
  
  // Engine output
  AddOutputField("AEROCDRAG", "AeroCDrag", FORMAT_SCIENTIFIC, "ENGINE_OUTPUT");
  AddOutputField("SOLIDCDRAG", "SolidCDrag", FORMAT_SCIENTIFIC, "ENGINE_OUTPUT");
  AddOutputField("RADIAL_DISTORTION", "Radial_Distortion", FORMAT_SCIENTIFIC, "ENGINE_OUTPUT");
  AddOutputField("CIRCUMFERENTIAL_DISTORTION", "Circumferential_Distortion", FORMAT_SCIENTIFIC, "ENGINE_OUTPUT");
  
  // Rotating Frame
  AddOutputField("MERIT", "CMerit", FORMAT_SCIENTIFIC, "ROTATING_FRAME");
  AddOutputField("CT", "CT", FORMAT_SCIENTIFIC, "ROTATING_FRAME");
  AddOutputField("CQ", "CQ", FORMAT_SCIENTIFIC, "ROTATING_FRAME");
  
  //Equivalent area
  AddOutputField("EQUIV_AREA", "CEquiv_Area", FORMAT_SCIENTIFIC, "EQUIVALENT_AREA");
  AddOutputField("NEARFIELD_OF", "CNearFieldOF", FORMAT_SCIENTIFIC, "EQUIVALENT_AREA");
  
  
}

inline void CFlowOutput::LoadOutput_Data(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) {
  unsigned short iVar;
  
  SetOutputFieldValue("INT_ITER", config[val_iZone]->GetIntIter());
  SetOutputFieldValue("EXT_ITER", config[val_iZone]->GetExtIter());
  SetOutputFieldValue("DENSITY", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(0)));
  SetOutputFieldValue("MOMENTUM-X", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(1)));
  SetOutputFieldValue("MOMENTUM-Y", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(2)));
  if (nDim == 2)
    SetOutputFieldValue("ENERGY", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(3)));
  else {
    SetOutputFieldValue("MOMENTUM-Z", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(3)));
    SetOutputFieldValue("ENERGY", log10(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetRes_RMS(4)));
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
  SetOutputFieldValue("TIME", timeused);
  SetOutputFieldValue("LINSOL_ITER", solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetIterLinSolver());
  
//  SetOutputFieldValue("AVG_MASSFLOW", config[val_iZone]->GetSurface_MassFlow(0));
//  SetOutputFieldValue("AVG_MACH",     config[val_iZone]->GetSurface_Mach(0));
//  SetOutputFieldValue("AVG_TEMP",     config[val_iZone]->GetSurface_Temperature(0));
//  SetOutputFieldValue("AVG_PRESS",    config[val_iZone]->GetSurface_Pressure(0));
//  SetOutputFieldValue("AVG_DENSITY",  config[val_iZone]->GetSurface_Density(0));
}

inline bool CFlowOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { return true;}

inline bool CFlowOutput::WriteScreen_Header(CConfig *config) {  
  bool write_header;
  write_header = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));
  
  return true;
}
inline bool CFlowOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {
  return true;
}

