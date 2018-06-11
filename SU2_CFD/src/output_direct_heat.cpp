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

  // Open the history file using only the master node
  if (rank == MASTER_NODE){
    cout << "History filename: " << char_histfile << endl;
    HistFile.open(char_histfile, ios::out);
    HistFile.precision(15);
    SetConvHistory_Header(config, val_iZone, INST_0);
  }

}

CHeatOutput::~CHeatOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }

}

void CHeatOutput::SetConvHistory_Header(CConfig *config, unsigned short val_iZone, unsigned short val_iInst) {
  char cstr[200], buffer[50], turb_resid[1000], adj_turb_resid[1000];
  unsigned short iMarker_Monitoring;
  string Monitoring_Tag, monitoring_coeff, aeroelastic_coeff, turbo_coeff;

  bool thermal = false; /* Flag for whether to print heat flux values */
  bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

  if (config->GetKind_Solver() == RANS || config->GetKind_Solver()  == NAVIER_STOKES) {
    thermal = true;
  }

  /*--- Begin of the header ---*/

  char begin[]= "\"Iteration\"";

  /*--- Header for the coefficients ---*/

  char heat_coeff[]= ",\"HeatFlux_Total\",\"HeatFlux_Maximum\",\"Temperature_Total\"";
  char heat_resid[]= ",\"Res_Heat\"";

  /*--- End of the header ---*/

  char end[]= ",\"Linear_Solver_Iterations\",\"CFL_Number\",\"Time(min)\"\n";
  char endfea[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";

  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
      (config->GetOutput_FileFormat() == FIELDVIEW) ||
      (config->GetOutput_FileFormat() == FIELDVIEW_BINARY)) {
    HistFile << "TITLE = \"SU2 Simulation\"" << endl;
    HistFile << "VARIABLES = ";
  }

  /*--- Write the header ---*/

  HistFile << begin << heat_coeff;
  HistFile << heat_resid << end;

  if (config->GetOutput_FileFormat() == TECPLOT ||
      config->GetOutput_FileFormat() == TECPLOT_BINARY ||
      config->GetOutput_FileFormat() == FIELDVIEW ||
      config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
    HistFile << "ZONE T= \"Convergence history\"" << endl;
  }

}


void CHeatOutput::SetConvHistory_Body(CGeometry ****geometry,
                                  CSolver *****solver_container,
                                  CConfig **config,
                                  CIntegration ****integration,
                                  bool DualTime_Iteration,
                                  su2double timeused,
                                  unsigned short val_iZone,
                                  unsigned short val_iInst) {

  unsigned long iIntIter    = config[val_iZone]->GetIntIter();
  unsigned long iExtIter    = config[val_iZone]->GetExtIter();
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned short nZone      = config[val_iZone]->GetnZone();
  unsigned short nInst      = config[val_iZone]->GetnTimeInstances();
  bool energy               = config[val_iZone]->GetEnergy_Equation();
  bool output_files         = true;

  /*--- Output using only the master node ---*/

  if (rank == MASTER_NODE) {

    unsigned long ExtIter_OffSet = config[val_iZone]->GetExtIter_OffSet();
    if (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST ||
        config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)
      ExtIter_OffSet = 0;

    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/

    char begin[1000], direct_coeff[1000], heat_coeff[1000], equivalent_area_coeff[1000], engine_coeff[1000], rotating_frame_coeff[1000], Cp_inverse_design[1000], Heat_inverse_design[1000], surface_coeff[1000], aeroelastic_coeff[1000], monitoring_coeff[10000],
    adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000], turb_resid[1000], trans_resid[1000],
    adj_turb_resid[1000],
    begin_fem[1000], fem_coeff[1000], heat_resid[1000], combo_obj[1000],
    fem_resid[1000], end[1000], end_fem[1000], surface_outputs[1000], d_direct_coeff[1000], turbo_coeff[10000];

    su2double dummy = 0.0, *Coord;
    unsigned short iVar, iMarker_Monitoring;

    unsigned long LinSolvIter = 0, iPointMaxResid;
    su2double timeiter = timeused/su2double(iExtIter+1);

    unsigned short nDim = geometry[val_iZone][val_iInst][FinestMesh]->GetnDim();

    bool heat =  ((config[val_iZone]->GetKind_Solver() == HEAT_EQUATION_FVM) || (config[val_iZone]->GetWeakly_Coupled_Heat()));
    bool weakly_coupled_heat = config[val_iZone]->GetWeakly_Coupled_Heat();

    unsigned short direct_diff = config[val_iZone]->GetDirectDiff();

    long ExtraHeatOutputZone = config[val_iZone]->GetExtraHeatOutputZone() - 1;
    bool extra_heat_output = false;
    su2double Extra_Total_Heat = 0.0;
    //su2double Extra_Total_Temperature = 0.0;
    su2double Extra_Heat_Residual = 0.0;

    if (ExtraHeatOutputZone > -1) {
      if (ExtraHeatOutputZone > nZone) {
        SU2_MPI::Error("Error in output routine: Extra output zone number exceeds total number of zones.", CURRENT_FUNCTION);
      }
      else if ((config[ExtraHeatOutputZone]->GetKind_Solver() != HEAT_EQUATION_FVM)) {
        SU2_MPI::Error("Error in output routine: No heat solver in extra output zone.", CURRENT_FUNCTION);
      }
      else {
        extra_heat_output = true;
      }
    }

    /*--- Initialize variables to store information from all domains (direct solution) ---*/

    su2double  Total_CHeat = 0.0, Total_Heat = 0.0, Total_MaxHeat = 0.0, Avg_TotalTemp = 0.0;

    unsigned short iSpan;

    /*--- Residual arrays ---*/
    su2double *residual_heat         = NULL;

    /*--- Coefficients Monitored arrays ---*/

    /*--- Initialize number of variables ---*/
    unsigned short nVar_Heat = 0;

    /*--- Direct problem variables ---*/
    if (heat) nVar_Heat = 1;

    /*--- Allocate memory for the residual ---*/
    residual_heat       = new su2double[nVar_Heat];

    /*--- Write information from nodes ---*/

    /*--- Heat coefficients  ---*/

    Total_Heat     = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_HeatFlux();
    Total_MaxHeat  = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_MaxHeatFlux();
    Avg_TotalTemp  = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_AvgTemperature();

    /*--- Heat Residuals ---*/

    for (iVar = 0; iVar < nVar_Heat; iVar++) {
      residual_heat[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetRes_RMS(iVar);
    }

    if (extra_heat_output) {
      Extra_Total_Heat      = solver_container[ExtraHeatOutputZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_HeatFlux();
      //Extra_Total_Temperature   = solver_container[ExtraHeatOutputZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_Temperature();
      Extra_Heat_Residual   = log10(solver_container[ExtraHeatOutputZone][val_iInst][FinestMesh][HEAT_SOL]->GetRes_RMS(0));
    }

    /*--- Header frequency ---*/

    bool Unsteady = ((config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
        (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
    bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
    bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));

    bool write_heads;
    if (Unsteady) write_heads = (iIntIter == 0);
    else write_heads = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*40)) == 0));

    bool write_turbo = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*40)) == 0) || (iExtIter == (config[val_iZone]->GetnExtIter() -1)));

    /*--- Analogous for dynamic problems (as of now I separate the problems, it may be worthy to do all together later on ---*/

    if (((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3))) {

      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {

        /*--- Write the begining of the history file ---*/
        SPRINTF(begin, "%12d", SU2_TYPE::Int(iExtIter+ExtIter_OffSet));

        /*--- Write the end of the history file ---*/
        SPRINTF (end, ", %12.10f, %12.10f, %12.10f\n", su2double(LinSolvIter), config[val_iZone]->GetCFL(MESH_0), timeused/60.0);

        /*--- Write the solution and residual of the history file ---*/
        SPRINTF (direct_coeff, ", %14.8e, %14.8e, %14.8e", Total_Heat, Total_MaxHeat, Avg_TotalTemp);
        SPRINTF (heat_resid, ", %14.8e", log10 (residual_heat[0]));

      }
      if ((val_iZone == 0 && val_iInst == 0)){
        /*--- Write the screen header---*/
        if (((write_heads) && !(!DualTime_Iteration && Unsteady))) {

          if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
          else cout << endl << " IntIter" << "  ExtIter";

          cout <<  "      Res[Heat]" << "   HFlux(Total)";

        }
      }

      /*--- Write the solution on the screen ---*/

      if ((val_iZone == 0 && val_iInst == 0)){
        cout.precision(6);
        cout.setf(ios::fixed, ios::floatfield);

        if (!Unsteady) {
          cout.width(5); cout << iExtIter + ExtIter_OffSet;
          cout.width(11); cout << timeiter;

        } else if (Unsteady && DualTime_Iteration) {
          cout.width(8); cout << iIntIter;
          cout.width(8); cout << iExtIter;
        }

      }

      if (!DualTime_Iteration) {
        HistFile << begin << direct_coeff << heat_resid << end;
        HistFile.flush();
      }

      cout.unsetf(ios::fixed);

    }


    delete [] residual_heat;

  }
}

inline void CHeatOutput::SetHistoryFile_Header(CConfig *config) { }

inline bool CHeatOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { }

inline void CHeatOutput::SetHistoryFile_Output(CConfig *config) { }

inline bool CHeatOutput::WriteScreen_Header(CConfig *config) { }

inline void CHeatOutput::SetScreen_Header(CConfig *config) { }

inline bool CHeatOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) { }

inline void CHeatOutput::SetScreen_Output(CConfig *config) { }

inline void CHeatOutput::LoadOutput_Data(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) { }

