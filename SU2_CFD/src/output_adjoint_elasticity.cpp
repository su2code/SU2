/*!
 * \file output_adjoint_mean.cpp
 * \brief Main subroutines for elasticity discrete adjoint output
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

CDiscAdjFEAOutput::CDiscAdjFEAOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {

  char buffer[50];

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

  // Open the history file using only the master node
  if (rank == MASTER_NODE){
    cout << "History filename: " << char_histfile << endl;
    HistFile.open(char_histfile, ios::out);
    HistFile.precision(15);
    SetConvHistory_Header(config, val_iZone, INST_0);
  }

}

CDiscAdjFEAOutput::~CDiscAdjFEAOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }


}

void CDiscAdjFEAOutput::SetConvHistory_Header(CConfig *config, unsigned short val_iZone, unsigned short val_iInst) {

  /*--- Begin of the header ---*/
  char begin[]= "\"Iteration\"";

  /*--- Header for the coefficients ---*/
  char fem_coeff[]= ",\"VM_Stress\",\"Force_Coeff\"";
  char fem_incload[]= ",\"IncLoad\"";

  /*--- Header for the residuals ---*/
  char fem_resid[]= ",\"Res_FEM[0]\",\"Res_FEM[1]\",\"Res_FEM[2]\"";

  /*--- End of the header ---*/
  char endfea[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";

  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
      (config->GetOutput_FileFormat() == FIELDVIEW) ||
      (config->GetOutput_FileFormat() == FIELDVIEW_BINARY)) {
    HistFile << "TITLE = \"SU2 Simulation\"" << endl;
    HistFile << "VARIABLES = ";
  }

  /*--- Write the header, case depending ---*/
  HistFile << begin << fem_coeff;
  HistFile << fem_resid << endfea;

  if (config->GetOutput_FileFormat() == TECPLOT ||
      config->GetOutput_FileFormat() == TECPLOT_BINARY ||
      config->GetOutput_FileFormat() == FIELDVIEW ||
      config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
    HistFile << "ZONE T= \"Convergence history\"" << endl;
  }

}


void CDiscAdjFEAOutput::SetConvHistory_Body(CGeometry ****geometry,
                                     CSolver *****solver_container,
                                     CConfig **config,
                                     CIntegration ****integration,
                                     bool DualTime_Iteration,
                                     su2double timeused,
                                     unsigned short val_iZone,
                                     unsigned short val_iInst) {


  bool fluid_structure      = (config[val_iZone]->GetFSI_Simulation());
  bool fea                  = ((config[val_iZone]->GetKind_Solver()== FEM_ELASTICITY)||(config[val_iZone]->GetKind_Solver()== DISC_ADJ_FEM));
  unsigned long iIntIter    = config[val_iZone]->GetIntIter();
  unsigned long iExtIter    = config[val_iZone]->GetExtIter();
  unsigned short nZone      = config[val_iZone]->GetnZone();
  bool incload              = config[val_iZone]->GetIncrementalLoad();
  bool output_files         = true;

  /*--- Output using only the master node ---*/

  if (rank == MASTER_NODE) {

    unsigned long ExtIter_OffSet = config[val_iZone]->GetExtIter_OffSet();

    su2double timeiter = timeused/su2double(iExtIter+1);

    unsigned short iVar;
    unsigned short nDim = geometry[val_iZone][INST_0][MESH_0]->GetnDim();

    bool fsi = (config[val_iZone]->GetFSI_Simulation());          // FEM structural solver.

    /*------------------------------------------------------------------------------------------------------*/
    /*--- Retrieve residual and extra data -----------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------*/

    /*--- Initialize number of variables ---*/
    unsigned short nVar_FEM = nDim;

    /*--- Allocate memory for the residual ---*/
    su2double *residual_fem          = NULL;
    residual_fem        = new su2double[nVar_FEM];

    /*--- FEM coefficients -- As of now, this is the Von Mises Stress ---*/
    su2double Total_CFEM = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetTotal_CFEA();
    su2double Total_SensE = 0.0, Total_SensNu = 0.0;

    /*--- Residuals: ---*/
    /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
    /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/

    for (iVar = 0; iVar < nVar_FEM; iVar++) {
      residual_fem[iVar] = solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetRes_RMS(iVar);
    }

    bool dynamic = (config[val_iZone]->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.

    /*------------------------------------------------------------------------------------------------------*/
    /*--- Write the history file ---------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------*/

    // Load data to buffers
    char begin_fem[1000], fem_coeff[1000], fem_resid[1000], end_fem[1000];

    SPRINTF (begin_fem, "%12d", SU2_TYPE::Int(iExtIter+ExtIter_OffSet));

    /*--- Initial variables ---*/

    /*--- FEM residual ---*/
    if (nVar_FEM == 2)
      SPRINTF (fem_resid, ", %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]));
    else
      SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), log10 (residual_fem[2]));

    // Write to history file

    HistFile << begin_fem << fem_resid << end_fem;
    HistFile.flush();

    /*------------------------------------------------------------------------------------------------------*/
    /*--- Write the screen header---------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------*/

    bool write_header = (iIntIter == 0);

    if (write_header) {

      cout << endl << " IntIter" << " ExtIter";

      if (nDim == 2) cout << "    Res[Ux_adj]" << "    Res[Uy_adj]" << "       Sens[E]" << "      Sens[Nu]"<<  endl;
      if (nDim == 3) cout << "    Res[Ux_adj]" << "    Res[Uy_adj]" << "    Res[Uz_adj]" << "       Sens[E]" << "      Sens[Nu]"<<  endl;

    }

    /*------------------------------------------------------------------------------------------------------*/
    /*--- Write the screen output---------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------*/

    cout.precision(6);
    cout.setf(ios::fixed, ios::floatfield);

    cout.width(15); cout << log10(residual_fem[0]);
    cout.width(15); cout << log10(residual_fem[1]);
    if (nDim == 3) { cout.width(15); cout << log10(residual_fem[2]); }

    cout.precision(4);
    cout.setf(ios::scientific, ios::floatfield);


    if (config[val_iZone]->GetnElasticityMod() == 1){
      cout.width(14); cout << solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_E(0);
      cout.width(14); cout << solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
    }
    else{
      Total_SensE = 0.0; Total_SensNu = 0.0;
      for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnElasticityMod(); iVar++){
          Total_SensE += solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_E(0)
              *solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_E(0);
          Total_SensNu += solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_Nu(0)
              *solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
      }
      Total_SensE = sqrt(Total_SensE);
      Total_SensNu = sqrt(Total_SensNu);
      cout.width(14); cout << Total_SensE;
      cout.width(14); cout << Total_SensNu;
    }

    cout << endl;

    cout.unsetf(ios::fixed);


    delete [] residual_fem;

  }
}

inline bool CDiscAdjFEAOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { }

inline bool CDiscAdjFEAOutput::WriteScreen_Header(CConfig *config) { }

inline bool CDiscAdjFEAOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) { }

inline void CDiscAdjFEAOutput::LoadOutput_Data(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) { }

