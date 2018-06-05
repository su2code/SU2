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

CDiscAdjFEAOutput::CDiscAdjFEAOutput(CConfig *config, unsigned short val_iZone) : COutput(config) {

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

  bool incload = config->GetIncrementalLoad();

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
  if (incload) HistFile << fem_incload;
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

    bool fem = ((config[val_iZone]->GetKind_Solver() == FEM_ELASTICITY) ||          // FEM structural solver.
                (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM));
    bool linear_analysis = (config[val_iZone]->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config[val_iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
    bool fsi = (config[val_iZone]->GetFSI_Simulation());          // FEM structural solver.
    bool discadj_fem = (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM);

    /*------------------------------------------------------------------------------------------------------*/
    /*--- Retrieve residual and extra data -----------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------*/

    /*--- Initialize number of variables ---*/
    unsigned short nVar_FEM = 0;
    if (linear_analysis) nVar_FEM = nDim;
    if (nonlinear_analysis) nVar_FEM = 3;

    /*--- Allocate memory for the residual ---*/
    su2double *residual_fem          = NULL;
    residual_fem        = new su2double[nVar_FEM];

    /*--- FEM coefficients -- As of now, this is the Von Mises Stress ---*/
    su2double Total_VMStress   = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetTotal_CFEA();
    su2double Total_ForceCoeff = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetForceCoeff();
    su2double Total_IncLoad    = solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetLoad_Increment();
    unsigned long LinSolvIter  = (unsigned long) solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetIterLinSolver();

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

    bool dynamic = (config[val_iZone]->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.

    /*------------------------------------------------------------------------------------------------------*/
    /*--- Write the history file ---------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------*/

    // Load data to buffers
    char begin_fem[1000], fem_coeff[1000], fem_resid[1000], end_fem[1000];

    SPRINTF (begin_fem, "%12d", SU2_TYPE::Int(iExtIter+ExtIter_OffSet));

    /*--- Initial variables ---*/
    if (incload) SPRINTF (fem_coeff, ", %14.8e, %14.8e, %14.8e", Total_VMStress, Total_ForceCoeff, Total_IncLoad);
    else SPRINTF (fem_coeff, ", %14.8e, %14.8e", Total_VMStress, Total_ForceCoeff);

    /*--- FEM residual ---*/
    if (nVar_FEM == 2)
      SPRINTF (fem_resid, ", %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]));
    else
      SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), log10 (residual_fem[2]));

    /*--- Linear solver iterations and time used ---*/
    SPRINTF (end_fem, ", %lu, %12.10f\n", LinSolvIter, timeused/60.0);

    // Write to history file

    HistFile << begin_fem << fem_coeff << fem_resid << end_fem;
    HistFile.flush();

    /*------------------------------------------------------------------------------------------------------*/
    /*--- Write the screen header---------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------*/

    bool write_header;
    if (nonlinear_analysis) write_header = (iIntIter == 0);
    else write_header = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*40)) == 0));

    if (write_header) {

      if (dynamic && nonlinear_analysis) {
        cout << endl << "Simulation time: " << config[val_iZone]->GetCurrent_DynTime() << ". Time step: " << config[val_iZone]->GetDelta_DynTime() << ".";
      }

      if (!nonlinear_analysis) cout << endl << " Iter" << "    Time(s)";
      else cout << endl << " IntIter" << " ExtIter";

      if (linear_analysis) {
        if (nDim == 2) cout << "    Res[Displx]" << "    Res[Disply]" << "      VMS(Max)"<<  endl;
        if (nDim == 3) cout << "    Res[Displx]" << "    Res[Disply]" << "    Res[Displz]" << "      VMS(Max)"<<  endl;
      }
      else if (nonlinear_analysis) {
        switch (config[val_iZone]->GetResidual_Criteria_FEM()) {
        case RESFEM_RELATIVE:
          cout << "     Res[UTOL]" << "     Res[RTOL]" << "     Res[ETOL]"  << "      VMS(Max)"<<  endl;
          break;
        case RESFEM_ABSOLUTE:
          cout << "   Res[UTOL-A]" << "   Res[RTOL-A]" << "   Res[ETOL-A]"  << "      VMS(Max)"<<  endl;
          break;
        default:
          cout << "     Res[UTOL]" << "     Res[RTOL]" << "     Res[ETOL]"  << "      VMS(Max)"<<  endl;
          break;
        }
      }

    }

    /*------------------------------------------------------------------------------------------------------*/
    /*--- Write the screen output---------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------*/

    if (!nonlinear_analysis) {
      cout.width(5); cout << iExtIter;
      cout.width(11); cout << timeiter;

    } else {
      cout.width(8); cout << iIntIter;
      cout.width(8); cout << iExtIter;
    }

    cout.precision(6);
    cout.setf(ios::fixed, ios::floatfield);
    if (linear_analysis) {
      cout.width(14); cout << log10(residual_fem[0]);
      cout.width(14); cout << log10(residual_fem[1]);
      if (nDim == 3) { cout.width(14); cout << log10(residual_fem[2]); }
    }
    else if (nonlinear_analysis) {
      cout.width(14); cout << log10(residual_fem[0]);
      cout.width(14); cout << log10(residual_fem[1]);
      cout.width(14); cout << log10(residual_fem[2]);
    }

    cout.precision(4);
    cout.setf(ios::scientific, ios::floatfield);
    cout.width(14); cout << Total_VMStress;
    cout << endl;

    cout.unsetf(ios::fixed);


    delete [] residual_fem;

  }
}
