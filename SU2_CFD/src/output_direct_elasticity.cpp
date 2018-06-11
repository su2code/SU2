/*!
 * \file output_direct_elasticity.cpp
 * \brief Main subroutines for FEA output
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

CFEAOutput::CFEAOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {

  // Open the history file using only the master node
  if (rank == MASTER_NODE){

    unsigned short nDim = geometry->GetnDim();

    bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.

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

    cout << "History filename: " << char_histfile << endl;
    HistFile.open(char_histfile, ios::out);
    HistFile.precision(15);
    SetConvHistory_Header(config, val_iZone, INST_0);

    /*--- Initialize number of variables ---*/
    if (linear_analysis) nVar_FEM = nDim;
    if (nonlinear_analysis) nVar_FEM = 3;

    /*--- Allocate memory for the residual ---*/
    residual_fem        = new su2double[nVar_FEM];

    /*--- Initialize  ---*/
    Total_VMStress = 0.0;
    Total_ForceCoeff = 0.0;
    Total_IncLoad = 0.0;
    LinSolvIter = 0.0;
    Time_Used = 0.0;

    iExtIter = 0;
    iIntIter = 0;
  }

}

CFEAOutput::~CFEAOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();

    if (residual_fem != NULL) delete [] residual_fem;
  }

}

void CFEAOutput::SetConvHistory_Header(CConfig *config, unsigned short val_iZone, unsigned short val_iInst) {

  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
      (config->GetOutput_FileFormat() == FIELDVIEW) ||
      (config->GetOutput_FileFormat() == FIELDVIEW_BINARY)) {
    HistFile << "TITLE = \"SU2 Simulation\"" << endl;
    HistFile << "VARIABLES = ";
  }

  /*--- Write the header, case depending ---*/
  SetHistoryFile_Header(config);

  if (config->GetOutput_FileFormat() == TECPLOT ||
      config->GetOutput_FileFormat() == TECPLOT_BINARY ||
      config->GetOutput_FileFormat() == FIELDVIEW ||
      config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
    HistFile << "ZONE T= \"Convergence history\"" << endl;
  }

}

void CFEAOutput::LoadOutput_Data(CGeometry ****geometry,
                                     CSolver *****solver_container,
                                     CConfig **config,
                                     CIntegration ****integration,
                                     bool DualTime_Iteration,
                                     su2double timeused,
                                     unsigned short val_iZone,
                                     unsigned short val_iInst) {

  bool fem = ((config[val_iZone]->GetKind_Solver() == FEM_ELASTICITY) ||          // FEM structural solver.
              (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM));
  bool linear_analysis = (config[val_iZone]->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config[val_iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool fsi = (config[val_iZone]->GetFSI_Simulation());          // FEM structural solver.
  bool discadj_fem = (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM);

  unsigned short iVar;
  unsigned short nDim = geometry[val_iZone][INST_0][MESH_0]->GetnDim();

  iExtIter = config[val_iZone]->GetExtIter();
  iIntIter = config[val_iZone]->GetIntIter();

  Time_Used = timeused;

  /*--- FEM coefficients -- As of now, this is the Von Mises Stress ---*/
  Total_VMStress   = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetTotal_CFEA();
  Total_ForceCoeff = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetForceCoeff();
  Total_IncLoad    = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetLoad_Increment();
  LinSolvIter  = (unsigned long) solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetIterLinSolver();

  /*--- Residuals: ---*/
  /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
  /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/

  if (linear_analysis) {
    for (iVar = 0; iVar < nVar_FEM; iVar++) {
      residual_fem[iVar] = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(iVar);
    }
  }
  else if (nonlinear_analysis) {
    for (iVar = 0; iVar < nVar_FEM; iVar++) {
      residual_fem[iVar] = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(iVar);
    }
  }

}

bool CFEAOutput::WriteScreen_Header(CConfig *config){

  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.

  bool write_header;
  if (nonlinear_analysis) write_header = (iIntIter == 0);
  else write_header = (((iExtIter % (config->GetWrt_Con_Freq()*40)) == 0));

  return write_header;

}

void CFEAOutput::SetScreen_Header(CConfig *config){

  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.

  bool absolute = (config->GetResidual_Criteria_FEM() == RESFEM_ABSOLUTE);

  if (dynamic && nonlinear_analysis) {
    cout << endl << "Simulation time: " << config->GetCurrent_DynTime() << ". Time step: " << config->GetDelta_DynTime() << ".";
  }

  // Insert line break
  cout << endl;
  // Evaluate the requested output
  for (unsigned short iField = 0; iField < config->GetnScreenOutput(); iField++){
    switch (config->GetScreenOutput_Field(iField)){
    case SOUT_INTITER: cout <<  " IntIter"; break;
    case SOUT_EXTITER: cout <<  " ExtIter"; break;
    case SOUT_UTOL:
      if (absolute) cout << "   Res[UTOL-A]";
      else cout << "     Res[UTOL]";
      break;
    case SOUT_RTOL:
      if (absolute) cout << "   Res[RTOL-A]";
      else cout << "     Res[RTOL]";
      break;
    case SOUT_ETOL:
      if (absolute) cout << "   Res[ETOL-A]";
      else cout << "     Res[ETOL]";
      break;
    case SOUT_DISPX: cout << "   Res[Displx]"; break;
    case SOUT_DISPY: cout << "   Res[Disply]"; break;
    case SOUT_DISPZ: cout << "   Res[Displz]"; break;
    case SOUT_VMS: cout << "      VMS(Max)"; break;
    }
  }
  // Insert line break
  cout << endl;

}

bool CFEAOutput::WriteScreen_Output(CConfig *config, bool write_dualtime){

  return true;

}

void CFEAOutput::SetScreen_Output(CConfig *config){

  // Evaluate the requested output
  for (unsigned short iField = 0; iField < config->GetnScreenOutput(); iField++){
    switch (config->GetScreenOutput_Field(iField)){
    case SOUT_INTITER:
      cout.width(8);
      cout << iIntIter;
      break;
    case SOUT_EXTITER:
      cout.width(8);
      cout << iExtIter;
      break;
    case SOUT_UTOL:
      cout.precision(6); cout.setf(ios::fixed, ios::floatfield); cout.width(14);
      cout << log10(residual_fem[0]);
      break;
    case SOUT_RTOL:
      cout.precision(6); cout.setf(ios::fixed, ios::floatfield); cout.width(14);
      cout << log10(residual_fem[1]);
      break;
    case SOUT_ETOL:
      cout.precision(6); cout.setf(ios::fixed, ios::floatfield); cout.width(14);
      cout << log10(residual_fem[2]);
      break;
    case SOUT_DISPX:
      cout.precision(6); cout.setf(ios::fixed, ios::floatfield); cout.width(14);
      cout << log10(residual_fem[0]);
      break;
    case SOUT_DISPY:
      cout.precision(6); cout.setf(ios::fixed, ios::floatfield); cout.width(14);
      cout << log10(residual_fem[1]);
      break;
    case SOUT_DISPZ:
      cout.precision(6); cout.setf(ios::fixed, ios::floatfield); cout.width(14);
      if (nVar_FEM == 3)
        cout << log10(residual_fem[2]);
      else
        cout << "              ";
      break;
    case SOUT_VMS:
      cout.precision(4); cout.setf(ios::scientific, ios::floatfield); cout.width(14);
      cout << Total_VMStress;
      break;
    }
  }
  cout << endl;
  cout.unsetf(ios::fixed);

}

void CFEAOutput::SetHistoryFile_Header(CConfig *config){

  // This buffer should be long enough
  char fem_header[1000]="";

  // Evaluate the requested output
  for (unsigned short iField = 0; iField < config->GetnHistoryOutput(); iField++){
    switch (config->GetHistoryOutput_Field(iField)){
    case HOUT_INTITER:
      SPRINTF (fem_header + strlen(fem_header), "\"Int_Iter\"");
      break;
    case HOUT_EXTITER:
      SPRINTF (fem_header + strlen(fem_header), "\"Ext_Iter\"");
      break;
    case HOUT_PHYSTIME:
      SPRINTF (fem_header + strlen(fem_header), "\"Time(min)\"");
      break;
    case HOUT_RESIDUALS:
      if (nVar_FEM == 2)
        SPRINTF (fem_header + strlen(fem_header), "\"Res_FEM[0]\",\"Res_FEM[1]\"");
      else
        SPRINTF (fem_header + strlen(fem_header), "\"Res_FEM[0]\",\"Res_FEM[1]\",\"Res_FEM[2]\"");
      break;
    case HOUT_LINSOL_ITER:
      SPRINTF (fem_header + strlen(fem_header), "\"Linear_Solver_Iterations\"");
      break;
    case HOUT_LOAD_RAMP:
      SPRINTF (fem_header + strlen(fem_header), "\"Load_Ramp\"");
      break;
    case HOUT_LOAD_INCREMENT:
      SPRINTF (fem_header + strlen(fem_header), "\"Load_Increment\"");
      break;
    case HOUT_VMS:
      SPRINTF (fem_header + strlen(fem_header), "\"VonMises_Stress\"");
      break;
    }
    // Print a comma in all fields but the last one
    if (iField < (config->GetnHistoryOutput() - 1))
      SPRINTF (fem_header + strlen(fem_header), ", ");
  }
  SPRINTF (fem_header + strlen(fem_header), "\n");

  HistFile << fem_header;
  HistFile.flush();

}

bool CFEAOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime){

  return true;

}

void CFEAOutput::SetHistoryFile_Output(CConfig *config){


  unsigned long ExtIter_OffSet = config->GetExtIter_OffSet();
  bool incload              = config->GetIncrementalLoad();

  // This buffer should be long enough
  char fem_output[1000]="";

  // Evaluate the requested output
  for (unsigned short iField = 0; iField < config->GetnHistoryOutput(); iField++){
    switch (config->GetHistoryOutput_Field(iField)){
    case HOUT_INTITER:
      SPRINTF (fem_output + strlen(fem_output), "%8d", SU2_TYPE::Int(iIntIter));
      break;
    case HOUT_EXTITER:
      SPRINTF (fem_output + strlen(fem_output), "%8d", SU2_TYPE::Int(iExtIter+ExtIter_OffSet));
      break;
    case HOUT_PHYSTIME:
      SPRINTF (fem_output + strlen(fem_output), "%12.10f", Time_Used/60.0);
      break;
    case HOUT_RESIDUALS:
      if (nVar_FEM == 2)
        SPRINTF (fem_output + strlen(fem_output), "%14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]));
      else
        SPRINTF (fem_output + strlen(fem_output), "%14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), log10 (residual_fem[2]));
      break;
    case HOUT_LINSOL_ITER:
      SPRINTF (fem_output + strlen(fem_output), "%lu", LinSolvIter);
      break;
    case HOUT_LOAD_RAMP:
      SPRINTF (fem_output + strlen(fem_output), "%14.8e", Total_ForceCoeff);
      break;
    case HOUT_LOAD_INCREMENT:
      SPRINTF (fem_output + strlen(fem_output), "%14.8e", Total_IncLoad);
      break;
    case HOUT_VMS:
      SPRINTF (fem_output + strlen(fem_output), "%14.8e", Total_VMStress);
      break;
    }
    // Print a comma in all fields but the last one
    if (iField < (config->GetnHistoryOutput() - 1))
      SPRINTF (fem_output + strlen(fem_output), ", ");
  }
  SPRINTF (fem_output + strlen(fem_output), "\n");

  HistFile << fem_output;
  HistFile.flush();

}


