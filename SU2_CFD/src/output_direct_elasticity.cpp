/*!
 * \file output_direct_elasticity.cpp
 * \brief Main subroutines for output solver information
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

CFEAOutput::CFEAOutput(CConfig *config) : COutput(config) {

}

CFEAOutput::~CFEAOutput(void) {

}

void CFEAOutput::SetConvHistory_Header(ofstream *ConvHist_file, CConfig *config, unsigned short val_iZone) {
  char cstr[200], buffer[50], turb_resid[1000], adj_turb_resid[1000];
  unsigned short iMarker_Monitoring;
  string Monitoring_Tag, monitoring_coeff, aeroelastic_coeff, turbo_coeff;

  bool rotating_frame = config->GetRotating_Frame();
  bool aeroelastic = config->GetAeroelastic_Simulation();
  bool equiv_area = config->GetEquivArea();
  bool engine        = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool cont_adj = config->GetContinuous_Adjoint();
  bool disc_adj = config->GetDiscrete_Adjoint();
  bool frozen_visc = (cont_adj && config->GetFrozen_Visc_Cont()) ||( disc_adj && config->GetFrozen_Visc_Disc());
  bool inv_design = (config->GetInvDesign_Cp() || config->GetInvDesign_HeatFlux());

  bool output_surface = (config->GetnMarker_Analyze() != 0);
  bool output_comboObj = (config->GetnObj() > 1);
  bool output_per_surface = config->GetWrt_Surface();
  bool turbo = config->GetBoolTurbomachinery();
  unsigned short direct_diff = config->GetDirectDiff();

  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool incload = config->GetIncrementalLoad();

  bool thermal = false; /* Flag for whether to print heat flux values */
  bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

  if (config->GetKind_Solver() == RANS || config->GetKind_Solver()  == NAVIER_STOKES) {
    thermal = true;
  }

  /*--- Write file name with extension ---*/
  string filename = config->GetConv_FileName();
  if(config->GetnZone() > 1){
    filename = config->GetMultizone_HistoryFileName(filename, val_iZone);
  }
  strcpy (cstr, filename.data());

  if (config->GetWrt_Unsteady() && config->GetRestart()) {
    long iExtIter = config->GetUnst_RestartIter();
    if (SU2_TYPE::Int(iExtIter) < 10) SPRINTF (buffer, "_0000%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d", SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d", SU2_TYPE::Int(iExtIter));
    strcat(cstr, buffer);
  }

  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == FIELDVIEW)) SPRINTF (buffer, ".dat");
  else if ((config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
           (config->GetOutput_FileFormat() == FIELDVIEW_BINARY))  SPRINTF (buffer, ".plt");
  else if (config->GetOutput_FileFormat() == PARAVIEW)  SPRINTF (buffer, ".csv");
  strcat(cstr, buffer);

  ConvHist_file->open(cstr, ios::out);
  ConvHist_file->precision(15);

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
    ConvHist_file[0] << "TITLE = \"SU2 Simulation\"" << endl;
    ConvHist_file[0] << "VARIABLES = ";
  }

  /*--- Write the header, case depending ---*/

  ConvHist_file[0] << begin << fem_coeff;
  if (incload) ConvHist_file[0] << fem_incload;
  ConvHist_file[0] << fem_resid << endfea;

  if (config->GetOutput_FileFormat() == TECPLOT ||
      config->GetOutput_FileFormat() == TECPLOT_BINARY ||
      config->GetOutput_FileFormat() == FIELDVIEW ||
      config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
    ConvHist_file[0] << "ZONE T= \"Convergence history\"" << endl;
  }

}


void CFEAOutput::SetConvHistory_Body(ofstream *ConvHist_file,
                                  CGeometry ***geometry,
                                  CSolver ****solver_container,
                                  CConfig **config,
                                  CIntegration ***integration,
                                  bool DualTime_Iteration,
                                  su2double timeused,
                                  unsigned short val_iZone) {


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
    if (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST ||
        config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)
      ExtIter_OffSet = 0;

    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/

    char begin_fem[1000], fem_coeff[1000], fem_resid[1000], end_fem[1000];

    su2double dummy = 0.0, *Coord;
    unsigned short iVar, iMarker_Monitoring;

    unsigned long LinSolvIter = 0, iPointMaxResid;
    su2double timeiter = timeused/su2double(iExtIter+1);

    unsigned short nDim = geometry[val_iZone][MESH_0]->GetnDim();

    bool fem = ((config[val_iZone]->GetKind_Solver() == FEM_ELASTICITY) ||          // FEM structural solver.
                (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM));
    bool linear_analysis = (config[val_iZone]->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config[val_iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
    bool fsi = (config[val_iZone]->GetFSI_Simulation());          // FEM structural solver.
    bool discadj_fem = (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM);

    /*--- Initialize variables to store information from all domains (direct solution) ---*/

    su2double Total_CFEM = 0.0;
    su2double Total_ForceCoeff = 0.0, Total_VMStress = 0.0, Total_IncLoad = 0.0;

    unsigned short iSpan;

    /*--- Residual arrays ---*/
    su2double *residual_fea          = NULL;
    su2double *residual_fem          = NULL;

    /*--- Initialize number of variables ---*/
    unsigned short nVar_FEM = 0;
    if (linear_analysis) nVar_FEM = nDim;
    if (nonlinear_analysis) nVar_FEM = 3;


    /*--- Allocate memory for the residual ---*/
    residual_fem        = new su2double[nVar_FEM];

    /*--- FEM coefficients -- As of now, this is the Von Mises Stress ---*/

    Total_VMStress = solver_container[val_iZone][MESH_0][FEA_SOL]->GetTotal_CFEA();

    Total_ForceCoeff = solver_container[val_iZone][MESH_0][FEA_SOL]->GetForceCoeff();

    Total_IncLoad = solver_container[val_iZone][MESH_0][FEA_SOL]->GetLoad_Increment();

    LinSolvIter = (unsigned long) solver_container[val_iZone][MESH_0][FEA_SOL]->GetIterLinSolver();

    /*--- Residuals: ---*/
    /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
    /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/

    if (linear_analysis) {
      for (iVar = 0; iVar < nVar_FEM; iVar++) {
        residual_fem[iVar] = solver_container[val_iZone][MESH_0][FEA_SOL]->GetRes_RMS(iVar);
      }
    }
    else if (nonlinear_analysis) {
      for (iVar = 0; iVar < nVar_FEM; iVar++) {
        residual_fem[iVar] = solver_container[val_iZone][MESH_0][FEA_SOL]->GetRes_FEM(iVar);
      }
    }

    /*--- Header frequency ---*/

    /*--- Header frequency: analogy for dynamic structural analysis ---*/
    /*--- DualTime_Iteration is a bool we receive, which is true if it comes from FEM_StructuralIteration and false from SU2_CFD ---*/
    /*--- We maintain the name, as it is an input of the function ---*/
    /*--- The function GetWrt_Con_Freq_DualTime should be modified to be able to define different frequencies ---*/
    /*--- dynamic determines if the problem is, or not, time dependent ---*/
    bool dynamic = (config[val_iZone]->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
    bool In_NoDynamic = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_Dynamic_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
    bool In_Dynamic_1 = (!DualTime_Iteration && nonlinear_analysis);
    bool In_Dynamic_2 = (nonlinear_analysis && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_Dynamic_3 = (nonlinear_analysis && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));

    /*--- Analogous for dynamic problems (as of now I separate the problems, it may be worthy to do all together later on ---*/
    bool write_heads_FEM;
    if (nonlinear_analysis) write_heads_FEM = (iIntIter == 0);
    else write_heads_FEM = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*40)) == 0));

    if (fem  && ( (In_NoDynamic || In_Dynamic_0 || In_Dynamic_1) && (In_NoDynamic || In_Dynamic_2 || In_Dynamic_3))){

      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {

        SPRINTF (begin_fem, "%12d", SU2_TYPE::Int(iExtIter+ExtIter_OffSet));

        if (incload) SPRINTF (fem_coeff, ", %14.8e, %14.8e, %14.8e", Total_VMStress, Total_ForceCoeff, Total_IncLoad);
        else SPRINTF (fem_coeff, ", %14.8e, %14.8e", Total_VMStress, Total_ForceCoeff);
        /*--- FEM residual ---*/
        if (nDim == 2) {
          if (linear_analysis) SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), dummy);
          if (nonlinear_analysis) SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), log10 (residual_fem[2]));
        }
        else {
          SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), log10 (residual_fem[2]));
        }
        SPRINTF (end_fem, ", %lu, %12.10f\n", LinSolvIter, timeused/60.0);

      }
        /*--- Write the screen header---*/
      if (fem && ((write_heads_FEM) && !(!DualTime_Iteration && nonlinear_analysis))) {

        if (dynamic) {
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

      /*--- Write the solution on the screen and history file ---*/


      if (!nonlinear_analysis) {
        cout.width(5); cout << iExtIter;
        cout.width(11); cout << timeiter;

      } else {
        cout.width(8); cout << iIntIter;
        cout.width(8); cout << iExtIter;
      }

          if (!DualTime_Iteration) {
            config[val_iZone]->GetHistFile()[0] << begin_fem << fem_coeff << fem_resid << end_fem;
            config[val_iZone]->GetHistFile()[0].flush();

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
          }

      cout.unsetf(ios::fixed);

    }


    delete [] residual_fea;
    delete [] residual_fem;


  }
}

