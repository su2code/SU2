/*!
 * \file output_direct_mean.cpp
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

CFlowOutput::CFlowOutput(CConfig *config, unsigned short val_iZone) : COutput(config) {

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
    SetConvHistory_Header(NULL, config, val_iZone, INST_0);
  }

}

CFlowOutput::~CFlowOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }


}

void CFlowOutput::SetConvHistory_Header(ofstream *ConvHist_file, CConfig *config, unsigned short val_iZone, unsigned short val_iInst) {
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

  /*--- Begin of the header ---*/

  char begin[]= "\"Iteration\"";

  /*--- Header for the coefficients ---*/

  char flow_coeff[]= ",\"CL\",\"CD\",\"CSF\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\",\"AoA\",\"Custom_ObjFunc\"";
  char heat_coeff[]= ",\"HeatFlux_Total\",\"HeatFlux_Maximum\",\"Temperature_Total\"";
  char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldOF\"";
  char engine_coeff[]= ",\"AeroCDrag\",\"SolidCDrag\",\"Radial_Distortion\",\"Circumferential_Distortion\"";
  char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
  char adj_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";
  char adj_inc_coeff[]=",\"Sens_Geo\",\"Sens_Vin\",\"Sens_Pout\",\"Sens_Temp\"";
  char adj_turbo_coeff[]=",\"Sens_Geo\",\"Sens_PressOut\",\"Sens_TotTempIn\"";
  char surface_outputs[]= ",\"Avg_MassFlow\",\"Avg_Mach\",\"Avg_Temp\",\"Avg_Press\",\"Avg_Density\",\"Avg_Enthalpy\",\"Avg_NormalVel\",\"Uniformity\",\"Secondary_Strength\",\"Momentum_Distortion\",\"Secondary_Over_Uniformity\",\"Avg_TotalTemp\",\"Avg_TotalPress\",\"Pressure_Drop\"";
  char Cp_inverse_design[]= ",\"Cp_Diff\"";
  char Heat_inverse_design[]= ",\"HeatFlux_Diff\"";
  char d_flow_coeff[] = ",\"D(CL)\",\"D(CD)\",\"D(CSF)\",\"D(CMx)\",\"D(CMy)\",\"D(CMz)\",\"D(CFx)\",\"D(CFy)\",\"D(CFz)\",\"D(CL/CD)\",\"D(Custom_ObjFunc)\"";
  char d_thermal_coeff[] = ",\"D(HeatFlux_Total)\",\"D(HeatFlux_Maximum)\"";
  char d_engine[] = ",\"D(AeroCDrag)\",\"D(SolidCDrag)\",\"D(Radial_Distortion)\",\"D(Circumferential_Distortion)\"";
  char d_turbo_coeff[] = ",\"D(TotalPressureLoss_0)\",\"D(FlowAngleOut_0)\",\"D(TotalEfficency)\",\"D(TotalStaticEfficiency)\", \"D(EntropyGen)\"";

  /*--- Find the markers being monitored and create a header for them ---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
    monitoring_coeff += ",\"CL_"  + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CD_"  + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CSF_" + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CL/CD_" + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFx_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFy_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFz_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMx_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMy_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMz_"    + Monitoring_Tag + "\"";
    aeroelastic_coeff += ",\"plunge_" + Monitoring_Tag + "\"";
    aeroelastic_coeff += ",\"pitch_"  + Monitoring_Tag + "\"";
  }

  if (turbo){
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_TurboPerformance(); iMarker_Monitoring++) {

      stringstream tag;
      tag << iMarker_Monitoring + 1;

      turbo_coeff += ",\"TotalPressureLoss_" + tag.str() + "\"";
      turbo_coeff += ",\"KineticEnergyLoss_" + tag.str() + "\"";
      turbo_coeff += ",\"EntropyGen_" + tag.str() + "\"";
      turbo_coeff += ",\"EulerianWork_" + tag.str() + "\"";
      turbo_coeff += ",\"PressureRatio_" + tag.str() + "\"";
      turbo_coeff += ",\"FlowAngleIn_" + tag.str() + "\"";
      turbo_coeff += ",\"FlowAngleOut_" + tag.str() + "\"";
      turbo_coeff += ",\"AbsFlowAngleIn_" + tag.str() + "\"";
      turbo_coeff += ",\"AbsFlowAngleOut_" + tag.str() + "\"";
      turbo_coeff += ",\"MassFlowIn_" + tag.str() + "\"";
      turbo_coeff += ",\"MassFlowOut_" + tag.str() + "\"";
      turbo_coeff += ",\"MachIn_" + tag.str() + "\"";
      turbo_coeff += ",\"MachOut_" + tag.str() + "\"";
      // different from zero only in multi-zone computation
      turbo_coeff += ",\"TotalEfficiency_" + tag.str() + "\"";
      turbo_coeff += ",\"TotalStaticEfficiency_" + tag.str() + "\"";

    }
  }

  char combo_obj[] = ",\"ComboObj\"";

  /*--- Header for the residuals ---*/

  char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
  char adj_flow_resid[]= ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";
  switch (config->GetKind_Turb_Model()) {
    case SA:case SA_NEG:case SA_E: case SA_COMP: case SA_E_COMP:
      SPRINTF (turb_resid, ",\"Res_Turb[0]\"");
      break;
    case SST:     SPRINTF (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\""); break;
  }
  switch (config->GetKind_Turb_Model()) {
    case SA:case SA_NEG:case SA_E: case SA_COMP: case SA_E_COMP:
      SPRINTF (adj_turb_resid, ",\"Res_AdjTurb[0]\"");
      break;
    case SST:     SPRINTF (adj_turb_resid, ",\"Res_AdjTurb[0]\",\"Res_AdjTurb[1]\""); break;
  }
  char heat_resid[]= ",\"Res_Heat\"";

  /*--- End of the header ---*/

  char end[]= ",\"Linear_Solver_Iterations\",\"CFL_Number\",\"Time(min)\"\n";

  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
      (config->GetOutput_FileFormat() == FIELDVIEW) ||
      (config->GetOutput_FileFormat() == FIELDVIEW_BINARY)) {
    HistFile << "TITLE = \"SU2 Simulation\"" << endl;
    HistFile << "VARIABLES = ";
  }

  /*--- Write the header, case depending ---*/

  HistFile << begin;
  if (!turbo) HistFile << flow_coeff;
  if (turbo) HistFile << turbo_coeff;
  if (thermal && !turbo) HistFile << heat_coeff;
  if (equiv_area) HistFile << equivalent_area_coeff;
  if (engine || actuator_disk) HistFile << engine_coeff;
  if (inv_design) {
    HistFile << Cp_inverse_design;
    if (thermal && !turbo) HistFile << Heat_inverse_design;
  }
  if (rotating_frame && !turbo) HistFile << rotating_frame_coeff;

  HistFile << flow_resid;
  if (turbulent) HistFile << turb_resid;
  if (weakly_coupled_heat) HistFile << heat_resid;
  if (aeroelastic) HistFile << aeroelastic_coeff;
  if (output_per_surface) HistFile << monitoring_coeff;
  if (output_surface) HistFile << surface_outputs;
  if (direct_diff != NO_DERIVATIVE) {
    if (!turbo) HistFile << d_flow_coeff;
    else        HistFile << d_turbo_coeff;
    if (engine || actuator_disk) HistFile << d_engine;
    if (thermal) HistFile << d_thermal_coeff;
  }
  if (output_comboObj) HistFile << combo_obj;
  HistFile << end;

  if (config->GetOutput_FileFormat() == TECPLOT ||
      config->GetOutput_FileFormat() == TECPLOT_BINARY ||
      config->GetOutput_FileFormat() == FIELDVIEW ||
      config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
    HistFile << "ZONE T= \"Convergence history\"" << endl;
  }

}


void CFlowOutput::SetConvHistory_Body(ofstream *ConvHist_file,
                                     CGeometry ****geometry,
                                     CSolver *****solver_container,
                                     CConfig **config,
                                     CIntegration ****integration,
                                     bool DualTime_Iteration,
                                     su2double timeused,
                                     unsigned short val_iZone,
                                     unsigned short val_iInst) {

  bool output_surface       = (config[val_iZone]->GetnMarker_Analyze() != 0);
  bool output_comboObj      = (config[val_iZone]->GetnObj() > 1);
  bool fluid_structure      = (config[val_iZone]->GetFSI_Simulation());
  bool fea                  = ((config[val_iZone]->GetKind_Solver()== FEM_ELASTICITY)||(config[val_iZone]->GetKind_Solver()== DISC_ADJ_FEM));
  unsigned long iIntIter    = config[val_iZone]->GetIntIter();
  unsigned long iExtIter    = config[val_iZone]->GetExtIter();
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned short nZone      = config[val_iZone]->GetnZone();
  unsigned short nInst      = config[val_iZone]->GetnTimeInstances();
  bool cont_adj             = config[val_iZone]->GetContinuous_Adjoint();
  bool disc_adj             = config[val_iZone]->GetDiscrete_Adjoint();
  bool energy               = config[val_iZone]->GetEnergy_Equation();
  bool incload              = config[val_iZone]->GetIncrementalLoad();
  bool output_files         = true;

  bool compressible = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);

  if (!disc_adj && !cont_adj && !DualTime_Iteration) {

    if ((config[val_iZone]->GetFixed_CL_Mode()) &&
        (config[val_iZone]->GetnExtIter()-config[val_iZone]->GetIter_dCL_dAlpha() - 1 < iExtIter)) {
      output_files = false;
    }

    if (fea || fluid_structure) output_files = false;

    /*--- We need to evaluate some of the objective functions to write the value on the history file ---*/

    if (((iExtIter % (config[val_iZone]->GetWrt_Sol_Freq())) == 0) ||
        ((!config[val_iZone]->GetFixed_CL_Mode()) && (iExtIter == (config[val_iZone]->GetnExtIter()-1))) ||
        /*--- If CL mode we need to compute the complete solution at two very particular iterations ---*/
        ((config[val_iZone]->GetFixed_CL_Mode()) && (iExtIter == (config[val_iZone]->GetnExtIter()-2))) ||
        ((config[val_iZone]->GetFixed_CL_Mode()) && (config[val_iZone]->GetnExtIter()-config[val_iZone]->GetIter_dCL_dAlpha() - 1 == iExtIter))) {


      if ((rank == MASTER_NODE) && output_files) cout << endl << "------------------------ Evaluate Special Output ------------------------";

      switch (config[val_iZone]->GetKind_Solver()) {

        case EULER: case NAVIER_STOKES: case RANS:

          /*--- For specific applications, evaluate and plot the surface. ---*/

          if (config[val_iZone]->GetnMarker_Analyze() != 0) {
            SpecialOutput_AnalyzeSurface(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                         geometry[val_iZone][val_iInst][MESH_0], config[ZONE_0], output_files);
          }

          /*--- For specific applications, evaluate and plot the surface. ---*/

          if (config[val_iZone]->GetnMarker_Analyze() != 0) {
            SpecialOutput_Distortion(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                     geometry[val_iZone][val_iInst][MESH_0], config[ZONE_0], output_files);
          }

          /*--- For specific applications, evaluate and plot the surface. ---*/

          if ((config[val_iZone]->GetnMarker_Analyze() != 0) && compressible) {
            SpecialOutput_Distortion(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                     geometry[val_iZone][val_iInst][MESH_0], config[ZONE_0], output_files);
          }

          /*--- For specific applications, evaluate and plot the equivalent area. ---*/

          if (config[val_iZone]->GetnMarker_NearFieldBound() != 0) {
            SpecialOutput_SonicBoom(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                    geometry[val_iZone][val_iInst][MESH_0], config[ZONE_0], output_files);
          }

          break;
      }

      /*--- Output a file with the forces breakdown. ---*/

      if (config[val_iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
        SpecialOutput_HarmonicBalance(solver_container, geometry, config, val_iInst, nInst, output_files);
      }

      /*--- Compute span-wise values file for turbomachinery. ---*/

      if (config[val_iZone]->GetBoolTurbomachinery()) {
        SpecialOutput_Turbo(solver_container, geometry, config, val_iZone, output_files);
      }

      /*--- Output a file with the forces breakdown. ---*/

      SpecialOutput_ForcesBreakdown(solver_container, geometry, config, val_iZone, output_files);

      if ((rank == MASTER_NODE) && !(fea || fluid_structure)) cout << endl;

      if ((rank == MASTER_NODE) && output_files) cout << "-------------------------------------------------------------------------" << endl << endl;

    }

  }

  /*--- Output using only the master node ---*/

  if (rank == MASTER_NODE) {

    /*-- Compute the total objective if a "combo" objective is used ---*/

    if (output_comboObj) {
      solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->SetTotal_ComboObj(0.0);
      switch (config[val_iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
        solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->Evaluate_ObjFunc(config[val_iZone]);
        break;
      }
    }

    unsigned long ExtIter_OffSet = config[val_iZone]->GetExtIter_OffSet();
    if (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST ||
        config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)
      ExtIter_OffSet = 0;

    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/

    char begin[1000], direct_coeff[1000], heat_coeff[1000], equivalent_area_coeff[1000], engine_coeff[1000], rotating_frame_coeff[1000], Cp_inverse_design[1000], Heat_inverse_design[1000], surface_coeff[1000], aeroelastic_coeff[1000], monitoring_coeff[10000],
    adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000], turb_resid[1000], trans_resid[1000],
    adj_turb_resid[1000], wave_coeff[1000],
    begin_fem[1000], fem_coeff[1000], wave_resid[1000], heat_resid[1000], combo_obj[1000],
    fem_resid[1000], end[1000], end_fem[1000], surface_outputs[1000], d_direct_coeff[1000], turbo_coeff[10000];

    su2double dummy = 0.0, *Coord;
    unsigned short iVar, iMarker_Monitoring;

    unsigned long LinSolvIter = 0, iPointMaxResid;
    su2double timeiter = timeused/su2double(iExtIter+1);

    unsigned short nDim = geometry[val_iZone][val_iInst][FinestMesh]->GetnDim();


    bool rotating_frame = config[val_iZone]->GetRotating_Frame();
    bool aeroelastic = config[val_iZone]->GetAeroelastic_Simulation();
    bool equiv_area = config[val_iZone]->GetEquivArea();
    bool engine        = ((config[val_iZone]->GetnMarker_EngineInflow() != 0) || (config[val_iZone]->GetnMarker_EngineExhaust() != 0));
    bool actuator_disk = ((config[val_iZone]->GetnMarker_ActDiskInlet() != 0) || (config[val_iZone]->GetnMarker_ActDiskOutlet() != 0));
    bool inv_design = (config[val_iZone]->GetInvDesign_Cp() || config[val_iZone]->GetInvDesign_HeatFlux());
    bool transition = (config[val_iZone]->GetKind_Trans_Model() == LM);
    bool thermal = (config[val_iZone]->GetKind_Solver() == RANS || config[val_iZone]->GetKind_Solver()  == NAVIER_STOKES);
    bool turbulent = ((config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS) ||
                      (config[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS));
    bool adjoint =  cont_adj || disc_adj;
    bool frozen_visc = (cont_adj && config[val_iZone]->GetFrozen_Visc_Cont()) ||( disc_adj && config[val_iZone]->GetFrozen_Visc_Disc());
    bool heat = (config[val_iZone]->GetKind_Solver() == HEAT_EQUATION) || (config[val_iZone]->GetKind_Solver() == HEAT_EQUATION_FVM) || (config[val_iZone]->GetWeakly_Coupled_Heat());
    bool weakly_coupled_heat = config[val_iZone]->GetWeakly_Coupled_Heat();
    bool flow = (config[val_iZone]->GetKind_Solver() == EULER) || (config[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
    (config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_EULER) ||
    (config[val_iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS);

    bool fsi = (config[val_iZone]->GetFSI_Simulation());          // FEM structural solver.

    bool turbo = config[val_iZone]->GetBoolTurbomachinery();

    unsigned short nTurboPerf  = config[val_iZone]->GetnMarker_TurboPerformance();

    bool output_per_surface = config[val_iZone]->GetWrt_Surface();

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
      else if ((config[ExtraHeatOutputZone]->GetKind_Solver() != HEAT_EQUATION) && (config[ExtraHeatOutputZone]->GetKind_Solver() != HEAT_EQUATION_FVM)) {
        SU2_MPI::Error("Error in output routine: No heat solver in extra output zone.", CURRENT_FUNCTION);
      }
      else {
        extra_heat_output = true;
      }
    }

    /*--- Initialize variables to store information from all domains (direct solution) ---*/

    su2double Total_CL = 0.0, Total_CD = 0.0, Total_CSF = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
    Total_CT = 0.0, Total_CQ = 0.0, Total_CWave = 0.0, Total_CHeat = 0.0,
    Total_Heat = 0.0, Total_MaxHeat = 0.0, Avg_TotalTemp = 0.0, Total_Custom_ObjFunc = 0.0,
    Total_ComboObj = 0.0, Total_AeroCD = 0.0, Total_SolidCD = 0.0, Total_IDR = 0.0, Total_IDC = 0.0,
    Total_AoA = 0.0;
    su2double Surface_MassFlow = 0.0, Surface_Mach = 0.0, Surface_Temperature = 0.0, Surface_Pressure = 0.0, Surface_Density = 0.0, Surface_Enthalpy = 0.0, Surface_NormalVelocity = 0.0, Surface_TotalTemperature = 0.0, Surface_TotalPressure = 0.0, Surface_Uniformity = 0.0, Surface_SecondaryStrength = 0.0,Surface_MomentumDistortion = 0.0, Surface_SecondOverUniform = 0.0, Surface_PressureDrop = 0.0;

    su2double Total_ForceCoeff = 0.0, Total_VMStress = 0.0, Total_IncLoad = 0.0;
    su2double Total_SensE = 0.0, Total_SensNu = 0.0;

    unsigned short iSpan;

    /*--- Initialize variables to store information from all domains (adjoint solution) ---*/
    su2double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
    su2double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0;

    su2double Total_Sens_BPressure = 0.0;
    su2double Total_Sens_Density = 0.0;
    su2double Total_Sens_ModVel = 0.0;

    /*--- Initialize variables to store information from all domains (direct differentiation) ---*/
    su2double D_Total_CL = 0.0, D_Total_CD = 0.0, D_Total_CSF = 0.0, D_Total_CMx = 0.0, D_Total_CMy = 0.0, D_Total_CMz = 0.0, D_Total_CEff = 0.0, D_Total_CFx = 0.0,
        D_Total_CFy = 0.0, D_Total_CFz = 0.0, D_Total_AeroCD = 0.0, D_Total_SolidCD = 0.0, D_Total_IDR = 0.0, D_Total_IDC = 0.0, D_Total_Custom_ObjFunc = 0.0, D_Total_Heat = 0.0, D_Total_MaxHeat = 0.0,
        D_TotalPressure_Loss = 0.0, D_FlowAngle_Out = 0.0, D_TotalStaticEfficiency = 0.0,
        D_TotalTotalEfficiency = 0.0, D_EntropyGen = 0.0;

    /*--- Residual arrays ---*/
    su2double *residual_flow         = NULL,
    *residual_turbulent    = NULL,
    *residual_transition   = NULL;
    su2double *residual_adjflow      = NULL,
    *residual_adjturbulent = NULL;
    su2double *residual_heat         = NULL;

    /*--- Coefficients Monitored arrays ---*/
    su2double *aeroelastic_plunge = NULL,
    *aeroelastic_pitch  = NULL,
    *Surface_CL         = NULL,
    *Surface_CD         = NULL,
    *Surface_CSF        = NULL,
    *Surface_CEff       = NULL,
    *Surface_CFx        = NULL,
    *Surface_CFy        = NULL,
    *Surface_CFz        = NULL,
    *Surface_CMx        = NULL,
    *Surface_CMy        = NULL,
    *Surface_CMz        = NULL;

    /*--- Initialize number of variables ---*/
    unsigned short nVar_Flow = 0, nVar_Turb = 0,
    nVar_Trans = 0, nVar_Heat = 0,
    nVar_AdjFlow = 0, nVar_AdjTurb = 0;

    /*--- Direct problem variables ---*/
    nVar_Flow = nDim+2;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()) {
        case SA: case SA_NEG: case SA_E: case SA_E_COMP: case SA_COMP: nVar_Turb = 1; break;
        case SST:    nVar_Turb = 2; break;
      }
    }
    if (transition) nVar_Trans = 2;
    if (heat) nVar_Heat = 1;

    /*--- Adjoint problem variables ---*/
    nVar_AdjFlow = nDim+2;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()) {
        case SA: case SA_NEG: case SA_E: case SA_E_COMP: case SA_COMP: nVar_AdjTurb = 1; break;
        case SST:    nVar_AdjTurb = 2; break;
      }
    }

    /*--- Allocate memory for the residual ---*/
    residual_flow       = new su2double[nVar_Flow];
    residual_turbulent  = new su2double[nVar_Turb];
    residual_transition = new su2double[nVar_Trans];
    residual_heat       = new su2double[nVar_Heat];

    residual_adjflow      = new su2double[nVar_AdjFlow];
    residual_adjturbulent = new su2double[nVar_AdjTurb];

    /*--- Allocate memory for the coefficients being monitored ---*/
    aeroelastic_plunge = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    aeroelastic_pitch  = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CL      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CD      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CSF = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CEff       = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFx        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFy        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFz        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];

    /*--- Write information from nodes ---*/

    /*--- Flow solution coefficients ---*/

    Total_CL             = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CL();
    Total_CD             = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CD();
    Total_CSF            = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CSF();
    Total_CEff           = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CEff();
    Total_CMx            = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CMx();
    Total_CMy            = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CMy();
    Total_CMz            = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CMz();
    Total_CFx            = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CFx();
    Total_CFy            = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CFy();
    Total_CFz            = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CFz();
    Total_ComboObj       = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_ComboObj();
    Total_AoA            = config[val_iZone]->GetAoA() - config[val_iZone]->GetAoA_Offset();
    Total_Custom_ObjFunc = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_Custom_ObjFunc();

    if (thermal) {
      Total_Heat     = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_HeatFlux();
      Total_MaxHeat  = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_MaxHeatFlux();
      Avg_TotalTemp  = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_AvgTemperature();

      if(weakly_coupled_heat) {
        Total_Heat     = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_HeatFlux();
        Total_MaxHeat  = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_MaxHeatFlux();
        Avg_TotalTemp  = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_AvgTemperature();
      }
    }

    if (direct_diff != NO_DERIVATIVE) {
      D_Total_CL             = SU2_TYPE::GetDerivative(Total_CL);
      D_Total_CD             = SU2_TYPE::GetDerivative(Total_CD);
      D_Total_CSF            = SU2_TYPE::GetDerivative(Total_CSF);
      D_Total_CEff           = SU2_TYPE::GetDerivative(Total_CEff);
      D_Total_CMx            = SU2_TYPE::GetDerivative(Total_CMx);
      D_Total_CMy            = SU2_TYPE::GetDerivative(Total_CMy);
      D_Total_CMz            = SU2_TYPE::GetDerivative(Total_CMz);
      D_Total_CFx            = SU2_TYPE::GetDerivative(Total_CFx);
      D_Total_CFy            = SU2_TYPE::GetDerivative(Total_CFy);
      D_Total_CFz            = SU2_TYPE::GetDerivative(Total_CFz);
      D_Total_Custom_ObjFunc = SU2_TYPE::GetDerivative(Total_Custom_ObjFunc);

      if (thermal) {
        D_Total_Heat    = SU2_TYPE::GetDerivative(Total_Heat);
        D_Total_MaxHeat = SU2_TYPE::GetDerivative(Total_MaxHeat);
        //Davg Temp
      }

      if (engine || actuator_disk) {
        D_Total_AeroCD   = SU2_TYPE::GetDerivative(Total_AeroCD);
        D_Total_SolidCD  = SU2_TYPE::GetDerivative(Total_SolidCD);
        D_Total_IDR      = SU2_TYPE::GetDerivative(Total_IDR);
        D_Total_IDC      = SU2_TYPE::GetDerivative(Total_IDC);
      }

    }

    if (equiv_area) {
      Total_CEquivArea    = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
      Total_CNearFieldOF  = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();

      Total_CEquivArea    = config[val_iZone]->GetWeightCd()*Total_CD + (1.0-config[val_iZone]->GetWeightCd())*Total_CEquivArea;
      Total_CNearFieldOF  = config[val_iZone]->GetWeightCd()*Total_CD + (1.0-config[val_iZone]->GetWeightCd())*Total_CNearFieldOF;
    }

    if (engine || actuator_disk) {
      Total_AeroCD  = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_AeroCD();
      Total_SolidCD = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_SolidCD();
      Total_IDR     = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_IDR();
      Total_IDC     = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_IDC();
    }

    if (rotating_frame) {
      Total_CT      = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CT();
      Total_CQ      = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CQ();
      Total_CMerit  = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
    }

    if (aeroelastic) {
      /*--- Look over the markers being monitored and get the desired values ---*/
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        aeroelastic_plunge[iMarker_Monitoring] = config[val_iZone]->GetAeroelastic_plunge(iMarker_Monitoring);
        aeroelastic_pitch[iMarker_Monitoring]  = config[val_iZone]->GetAeroelastic_pitch(iMarker_Monitoring);
      }
    }

    if (output_per_surface) {
      /*--- Look over the markers being monitored and get the desired values ---*/
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Surface_CL[iMarker_Monitoring]      = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CL(iMarker_Monitoring);
        Surface_CD[iMarker_Monitoring]      = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CD(iMarker_Monitoring);
        Surface_CSF[iMarker_Monitoring] = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CSF(iMarker_Monitoring);
        Surface_CEff[iMarker_Monitoring]       = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CEff(iMarker_Monitoring);
        Surface_CFx[iMarker_Monitoring]        = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CFx(iMarker_Monitoring);
        Surface_CFy[iMarker_Monitoring]        = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CFy(iMarker_Monitoring);
        Surface_CFz[iMarker_Monitoring]        = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CFz(iMarker_Monitoring);
        Surface_CMx[iMarker_Monitoring]        = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring);
        Surface_CMy[iMarker_Monitoring]        = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring);
        Surface_CMz[iMarker_Monitoring]        = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring);
      }
    }

    if (turbo) {
      /*--- Loop over the nMarker of turboperformance and get the desired values ---*/
      for (iMarker_Monitoring = 0; iMarker_Monitoring < nTurboPerf; iMarker_Monitoring++) {
        for(iSpan=0; iSpan<nSpanWiseSections+1; iSpan++){
          if ((iMarker_Monitoring == 0) && (direct_diff != NO_DERIVATIVE)){
            D_TotalPressure_Loss = SU2_TYPE::GetDerivative(TotalPressureLoss[iMarker_Monitoring][iSpan]);
            D_FlowAngle_Out      = 180.0/PI_NUMBER*SU2_TYPE::GetDerivative(FlowAngleOut[iMarker_Monitoring][iSpan]);
          }
        }
      }
      if (direct_diff != NO_DERIVATIVE){
        D_TotalStaticEfficiency = SU2_TYPE::GetDerivative(TotalStaticEfficiency[nTurboPerf-1][nSpanWiseSections]);
        D_TotalTotalEfficiency  = SU2_TYPE::GetDerivative(TotalTotalEfficiency[nTurboPerf-1][nSpanWiseSections]);
        D_EntropyGen            = SU2_TYPE::GetDerivative(EntropyGen[nTurboPerf-1][nSpanWiseSections]);
      }
    }

    /*--- Get flux-averaged values at the specified surface ---*/

    if (output_surface) {

      unsigned short iMarker_Analyze = 0;
      Surface_MassFlow = config[ZONE_0]->GetSurface_MassFlow(iMarker_Analyze);
      Surface_Mach = config[ZONE_0]->GetSurface_Mach(iMarker_Analyze);
      Surface_Temperature = config[ZONE_0]->GetSurface_Temperature(iMarker_Analyze);
      Surface_Pressure = config[ZONE_0]->GetSurface_Pressure(iMarker_Analyze);
      Surface_Density = config[ZONE_0]->GetSurface_Density(iMarker_Analyze);
      Surface_Enthalpy = config[ZONE_0]->GetSurface_Enthalpy(iMarker_Analyze);
      Surface_NormalVelocity = config[ZONE_0]->GetSurface_NormalVelocity(iMarker_Analyze);
      Surface_Uniformity = config[ZONE_0]->GetSurface_Uniformity(iMarker_Analyze);
      Surface_SecondaryStrength = config[ZONE_0]->GetSurface_SecondaryStrength(iMarker_Analyze);
      Surface_MomentumDistortion = config[ZONE_0]->GetSurface_MomentumDistortion(iMarker_Analyze);
      Surface_SecondOverUniform = config[ZONE_0]->GetSurface_SecondOverUniform(iMarker_Analyze);
      Surface_TotalTemperature = config[ZONE_0]->GetSurface_TotalTemperature(iMarker_Analyze);
      Surface_TotalPressure = config[ZONE_0]->GetSurface_TotalPressure(iMarker_Analyze);
      Surface_PressureDrop = config[ZONE_0]->GetSurface_PressureDrop(iMarker_Analyze);

    }

    /*--- Flow Residuals ---*/

    for (iVar = 0; iVar < nVar_Flow; iVar++)
      residual_flow[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetRes_RMS(iVar);

    /*--- Turbulent residual ---*/

    if (turbulent) {
      for (iVar = 0; iVar < nVar_Turb; iVar++)
        residual_turbulent[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][TURB_SOL]->GetRes_RMS(iVar);
    }

    if (weakly_coupled_heat) {
      for (iVar = 0; iVar < nVar_Heat; iVar++) {
        residual_heat[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetRes_RMS(iVar);
      }

    }

    /*--- Transition residual ---*/

    if (transition) {
      for (iVar = 0; iVar < nVar_Trans; iVar++)
        residual_transition[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][TRANS_SOL]->GetRes_RMS(iVar);
    }

    /*--- Iterations of the linear solver ---*/

    LinSolvIter = (unsigned long) solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetIterLinSolver();

    /*--- Adjoint solver ---*/

    if (adjoint) {

      /*--- Adjoint solution coefficients ---*/

      Total_Sens_Geo       = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
      Total_Sens_Mach      = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
      Total_Sens_AoA       = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0;
      Total_Sens_Press     = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Press();
      Total_Sens_Temp      = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Temp();
      Total_Sens_BPressure = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_BPress();
      Total_Sens_Density   = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Density();
      Total_Sens_ModVel    = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_ModVel();

      /*--- Adjoint flow residuals ---*/

      for (iVar = 0; iVar < nVar_AdjFlow; iVar++) {
        residual_adjflow[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(iVar);
      }

      /*--- Adjoint turbulent residuals ---*/

      if (turbulent) {
        if (!frozen_visc) {
          for (iVar = 0; iVar < nVar_AdjTurb; iVar++)
            residual_adjturbulent[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][ADJTURB_SOL]->GetRes_RMS(iVar);
        }
      }

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

    if ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3)) {

      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {

        /*--- Write the begining of the history file ---*/
        SPRINTF(begin, "%12d", SU2_TYPE::Int(iExtIter+ExtIter_OffSet));

        /*--- Write the end of the history file ---*/
        SPRINTF (end, ", %12.10f, %12.10f, %12.10f\n", su2double(LinSolvIter), config[val_iZone]->GetCFL(MESH_0), timeused/60.0);

        /*--- Write the solution and residual of the history file ---*/

            /*--- Direct coefficients ---*/
            SPRINTF (direct_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e",
                     Total_CL, Total_CD, Total_CSF, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff, Total_AoA, Total_Custom_ObjFunc);
            if (thermal || heat) SPRINTF (heat_coeff, ", %14.8e, %14.8e, %14.8e",  Total_Heat, Total_MaxHeat, Avg_TotalTemp);
            if (equiv_area) SPRINTF (equivalent_area_coeff, ", %14.8e, %14.8e", Total_CEquivArea, Total_CNearFieldOF);
            if (engine || actuator_disk) SPRINTF (engine_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e", Total_AeroCD, Total_SolidCD, Total_IDR, Total_IDC);
            if (rotating_frame) SPRINTF (rotating_frame_coeff, ", %14.8e, %14.8e, %14.8e", Total_CMerit, Total_CT, Total_CQ);
            if (inv_design) {
              SPRINTF (Cp_inverse_design, ", %14.8e", solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CpDiff());
              if (thermal && !turbo) SPRINTF (Heat_inverse_design, ", %14.8e", solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_HeatFluxDiff());
            }

            if (direct_diff != NO_DERIVATIVE) {
              if (!turbo)
                SPRINTF (d_direct_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e",
                         D_Total_CL, D_Total_CD, D_Total_CSF, D_Total_CMx, D_Total_CMy, D_Total_CMz, D_Total_CFx, D_Total_CFy,
                         D_Total_CFz, D_Total_CEff, D_Total_Custom_ObjFunc);
              else
                SPRINTF (d_direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", D_TotalPressure_Loss, D_FlowAngle_Out,
                         D_TotalTotalEfficiency, D_TotalStaticEfficiency, D_EntropyGen);
              if (engine || actuator_disk)
              SPRINTF (d_direct_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e",
                       D_Total_CL, D_Total_CD, D_Total_CSF, D_Total_CMx, D_Total_CMy, D_Total_CMz, D_Total_CFx, D_Total_CFy,
                       D_Total_CFz, D_Total_CEff, D_Total_Custom_ObjFunc, D_Total_AeroCD, D_Total_SolidCD, D_Total_IDR, D_Total_IDC);
              if (thermal)
                SPRINTF (d_direct_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e",
                         D_Total_CL, D_Total_CD, D_Total_CSF, D_Total_CMx, D_Total_CMy, D_Total_CMz, D_Total_CFx, D_Total_CFy,
                         D_Total_CFz, D_Total_CEff, D_Total_Custom_ObjFunc, D_Total_Heat, D_Total_MaxHeat);
            }

            if (aeroelastic) {
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                //Append one by one the surface coeff to aeroelastic coeff. (Think better way do this, maybe use string)
                if (iMarker_Monitoring == 0) {
                  SPRINTF(aeroelastic_coeff, ", %12.10f", aeroelastic_plunge[iMarker_Monitoring]);
                }
                else {
                  SPRINTF(surface_coeff, ", %12.10f", aeroelastic_plunge[iMarker_Monitoring]);
                  strcat(aeroelastic_coeff, surface_coeff);
                }
                SPRINTF(surface_coeff, ", %12.10f", aeroelastic_pitch[iMarker_Monitoring]);
                strcat(aeroelastic_coeff, surface_coeff);
              }
            }

            if (output_per_surface) {
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                //Append one by one the surface coeff to monitoring coeff. (Think better way do this, maybe use string)
                if (iMarker_Monitoring == 0) {
                  SPRINTF(monitoring_coeff, ", %12.10f", Surface_CL[iMarker_Monitoring]);
                }
                else {
                  SPRINTF(surface_coeff, ", %12.10f", Surface_CL[iMarker_Monitoring]);
                  strcat(monitoring_coeff, surface_coeff);
                }
                SPRINTF(surface_coeff, ", %12.10f", Surface_CD[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CSF[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CEff[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CFx[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CFy[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CFz[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CMx[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CMy[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", Surface_CMz[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
              }
            }

            if (turbo){
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_TurboPerformance(); iMarker_Monitoring++){
                if (iMarker_Monitoring == 0){
                  SPRINTF(turbo_coeff, ", %12.10f", TotalPressureLoss[iMarker_Monitoring][nSpanWiseSections]);
                }else{
                  SPRINTF(surface_coeff, ", %12.10f", TotalPressureLoss[iMarker_Monitoring][nSpanWiseSections]);
                  strcat(turbo_coeff, surface_coeff);
                }
                SPRINTF(surface_coeff, ", %12.10f", KineticEnergyLoss[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", EntropyGen[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", EulerianWork[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", PressureRatio[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", 180.0/PI_NUMBER*FlowAngleIn[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", 180.0/PI_NUMBER*FlowAngleOut[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", 180.0/PI_NUMBER*AbsFlowAngleIn[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", 180.0/PI_NUMBER*AbsFlowAngleOut[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", MassFlowIn[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", MassFlowOut[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", sqrt(MachIn[iMarker_Monitoring][nSpanWiseSections][1]*MachIn[iMarker_Monitoring][nSpanWiseSections][1] + MachIn[iMarker_Monitoring][nSpanWiseSections][0]*MachIn[iMarker_Monitoring][nSpanWiseSections][0]));
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", sqrt(MachOut[iMarker_Monitoring][nSpanWiseSections][1]*MachOut[iMarker_Monitoring][nSpanWiseSections][1] + MachOut[iMarker_Monitoring][nSpanWiseSections][0]*MachOut[iMarker_Monitoring][nSpanWiseSections][0]));
                strcat(turbo_coeff, surface_coeff);
                //
                SPRINTF(surface_coeff, ", %12.10f", TotalTotalEfficiency[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);
                SPRINTF(surface_coeff, ", %12.10f", TotalStaticEfficiency[iMarker_Monitoring][nSpanWiseSections]);
                strcat(turbo_coeff, surface_coeff);

              }
            }


            /*--- Flow residual ---*/
            if (nDim == 2) {
              if (compressible) SPRINTF (flow_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy);
              if (incompressible) SPRINTF (flow_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy);
            }
            else {
              if (compressible) SPRINTF (flow_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), log10 (residual_flow[4]) );
              if (incompressible) SPRINTF (flow_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), log10 (residual_flow[4]));
            }

            /*--- Turbulent residual ---*/
            if (turbulent) {
              switch(nVar_Turb) {
                case 1: SPRINTF (turb_resid, ", %12.10f", log10 (residual_turbulent[0])); break;
                case 2: SPRINTF (turb_resid, ", %12.10f, %12.10f", log10(residual_turbulent[0]), log10(residual_turbulent[1])); break;
              }
            }

            /*---- Averaged stagnation pressure at an exit ----*/

            if (output_surface) {
              SPRINTF( surface_outputs, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", Surface_MassFlow, Surface_Mach, Surface_Temperature, Surface_Pressure, Surface_Density, Surface_Enthalpy, Surface_NormalVelocity, Surface_Uniformity, Surface_SecondaryStrength, Surface_MomentumDistortion, Surface_SecondOverUniform, Surface_TotalTemperature, Surface_TotalPressure, Surface_PressureDrop);
            }

            /*--- Transition residual ---*/
            if (transition) {
              SPRINTF (trans_resid, ", %12.10f, %12.10f", log10(residual_transition[0]), log10(residual_transition[1]));
            }

            /*--- Combo objective ---*/
            if (output_comboObj) {
              SPRINTF(combo_obj,", %12.10f", Total_ComboObj);
            }

            if (adjoint) {

              /*--- Adjoint coefficients ---*/
              if (!turbo) {
                if (compressible) {
                  SPRINTF (adjoint_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);
                }
                if (incompressible) {
                  SPRINTF (adjoint_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e", Total_Sens_Geo, Total_Sens_ModVel, Total_Sens_BPressure, Total_Sens_Temp);
                }
              } else
                SPRINTF (adjoint_coeff, ", %14.8e, %14.8e, %14.8e", Total_Sens_Geo, Total_Sens_BPressure, Total_Sens_Temp);

              /*--- Adjoint flow residuals ---*/
              if (nDim == 2) {
                if (compressible) SPRINTF (adj_flow_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, 0.0", log10 (residual_adjflow[0]), log10 (residual_adjflow[1]), log10 (residual_adjflow[2]), log10 (residual_adjflow[3]) );
                if (incompressible) SPRINTF (adj_flow_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, 0.0", log10 (residual_adjflow[0]), log10 (residual_adjflow[1]), log10 (residual_adjflow[2]), log10 (residual_adjflow[3]) );
              }
              else {
                if (compressible) SPRINTF (adj_flow_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_adjflow[0]), log10 (residual_adjflow[1]), log10 (residual_adjflow[2]), log10 (residual_adjflow[3]), log10 (residual_adjflow[4]) );
                if (incompressible) SPRINTF (adj_flow_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_adjflow[0]), log10 (residual_adjflow[1]), log10 (residual_adjflow[2]), log10 (residual_adjflow[3]), log10 (residual_adjflow[4]));
              }

              /*--- Adjoint turbulent residuals ---*/
              if (turbulent)
                if (!frozen_visc) {
                  if (nVar_AdjTurb == 1) {
                  SPRINTF (adj_turb_resid, ", %14.8e", log10 (residual_adjturbulent[0]));
                  } else if (nVar_AdjTurb > 1) {
                    SPRINTF (adj_turb_resid, ", %14.8e, %14.8e", log10 (residual_adjturbulent[0]), log10 (residual_adjturbulent[1]));
                  }
                }

            }

            if (weakly_coupled_heat) {
              SPRINTF (heat_resid, ", %14.8e", log10 (residual_heat[0]));
            }

      }
      if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
        /*--- Write the screen header---*/
        if ((write_heads) && !(!DualTime_Iteration && Unsteady)){

            if (!Unsteady && (config[val_iZone]->GetUnsteady_Simulation() != TIME_STEPPING)) {

                cout << endl << "---------------------- Local Time Stepping Summary ----------------------" << endl;

                for (unsigned short iMesh = FinestMesh; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
                  cout << "MG level: "<< iMesh << " -> Min. DT: " << solver_container[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetMin_Delta_Time()<<
                  ". Max. DT: " << solver_container[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetMax_Delta_Time() <<
                  ". CFL: " << config[val_iZone]->GetCFL(iMesh)  << "." << endl;

                  if (nZone > 1)
                    cout << "CFL in zone 2: " << config[1]->GetCFL(MESH_0) << endl;

                cout << "-------------------------------------------------------------------------" << endl;

                if (direct_diff != NO_DERIVATIVE) {
                  cout << endl << "---------------------- Direct Differentiation Summary -------------------" << endl;
                  cout << "Coefficients are differentiated with respect to ";
                  switch (direct_diff) {
                  case D_MACH:
                    cout << "Mach number." << endl;
                    break;
                  case D_AOA:
                    cout << "AoA." << endl;
                    break;
                  case D_SIDESLIP:
                    cout << "AoS." << endl;
                    break;
                  case D_REYNOLDS:
                    cout << "Reynolds number." << endl;
                    break;
                  case D_TURB2LAM:
                    cout << "Turb/Lam ratio." << endl;
                    break;
                  case D_PRESSURE:
                    cout << "Freestream Pressure." << endl;
                    break;
                  case D_TEMPERATURE:
                    cout << "Freestream Temperature." << endl;
                    break;
                  case D_DENSITY:
                    cout << "Freestream Density." << endl;
                    break;
                  case D_VISCOSITY:
                    cout << "Freestream Viscosity." << endl;
                    break;
                  case D_DESIGN:
                    cout << "Design Variables." << endl;
                    break;
                  default:
                    break;
                  }

                  cout << "    D_CLift(Total)" << "    D_CDrag(Total)" << "      D_CMz(Total)" <<"     D_CEff(Total)" << endl;
                  cout.width(18); cout << D_Total_CL;
                  cout.width(18); cout << D_Total_CD;
                  cout.width(18); cout << D_Total_CMz;
                  cout.width(18); cout << D_Total_CEff;
                  cout << endl << "-------------------------------------------------------------------------" << endl;
                  cout << endl;
                }
                if (turbo && write_turbo && val_iZone== 0){
                  WriteTurboPerfConvHistory(config[val_iZone]);
                }
            }
            else {
              if (flow) {
                if ((config[val_iZone]->GetUnsteady_Simulation() == TIME_STEPPING) && (config[val_iZone]->GetUnst_CFL()== 0.0))
                {
                  cout << endl << "Min DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<< ".Max DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
                } else if ((config[val_iZone]->GetUnsteady_Simulation() == TIME_STEPPING) && (config[val_iZone]->GetUnst_CFL()!= 0.0)) {
                  cout << endl << "Min DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<< ".Max DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ". Time step: " << solver_container[val_iZone][val_iInst][config[val_iZone]->GetFinestMesh()][FLOW_SOL]->GetMin_Delta_Time() << ". CFL: " << config[val_iZone]->GetUnst_CFL()<<".";
                } else {
                  cout << endl << "Min DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<< ".Max DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
                }
              } else {
                cout << endl << "Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
              }
            }


          switch (config[val_iZone]->GetKind_Solver()) {
          case EULER :                  case NAVIER_STOKES:

            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetPoint_Max_Coord(0);

            cout << endl << "----------------------- Residual Evolution Summary ----------------------" << endl;

            cout << "log10[Maximum residual]: " << log10(solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetRes_Max(0)) << "." << endl;

            if (config[val_iZone]->GetSystemMeasurements() == SI) {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
              if (nDim == 3) cout << ", " << Coord[2];
              cout <<   ")." << endl;
            }
            else {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0]*12.0 << ", " << Coord[1]*12.0;
              if (nDim == 3) cout << ", " << Coord[2]*12.0;
              cout <<   ")." << endl;
            }

            /*--- Print out the number of non-physical points and reconstructions ---*/

            if (config[val_iZone]->GetNonphysical_Points() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
            if (config[val_iZone]->GetNonphysical_Reconstr() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Reconstr() << " non-physical states in the upwind reconstruction." << endl;

            cout << "-------------------------------------------------------------------------" << endl;

            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";

            //            if (!fluid_structure) {
              if (incompressible && !weakly_coupled_heat) {
              if (energy) {cout << "   Res[Press]" << "     Res[Temp]" << "   CLift(Total)" << "   CDrag(Total)" << endl;}
              else {cout << "   Res[Press]" << "     Res[Velx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;}
              }
              else if (incompressible && weakly_coupled_heat) cout << "   Res[Press]" << "     Res[Heat]" << "   HFlux(Total)";
            else if (rotating_frame && nDim == 3 && !turbo) cout << "     Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
            else if (aeroelastic) cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch" << endl;
            else if (equiv_area) cout << "     Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;

            else if (turbo){

              if(nZone  < 2){
                /*--- single zone output ---*/
                cout << "     Res[Rho]" << "     Res[RhoE]"  << "  TotPresLoss(%)" << "  Entropy Gen.(%)";
              }
              else{
                /* --- multi-zone output ---*/
                cout << "     Res[Rho]" << "     Res[RhoE]"  << " TTEfficiency(%)" << " Entropy Gen.(%)";
              }
            }

            else if (actuator_disk) cout << "     Res[Rho]" << "     Res[RhoE]" << "      CL(Total)" << "   CD-CT(Total)";
            else if (engine) cout << "     Res[Rho]" << "     Res[RhoE]" << "      CL(Total)" << "   CD-CT(Total)";
            else cout << "     Res[Rho]" << "     Res[RhoE]" << "      CL(Total)" << "      CD(Total)";

            if(extra_heat_output) {
              cout <<  "     Res[Heat]" << "   HFlux(Total)";
            }

            cout << endl;

            break;

          case RANS :

            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetPoint_Max_Coord(0);

            cout << endl << "----------------------- Residual Evolution Summary ----------------------" << endl;

            cout << "log10[Maximum residual]: " << log10(solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetRes_Max(0)) << "." << endl;
            if (config[val_iZone]->GetSystemMeasurements() == SI) {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
              if (nDim == 3) cout << ", " << Coord[2];
              cout <<   ")." << endl;
            }
            else {
              cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0]*12.0 << ", " << Coord[1]*12.0;
              if (nDim == 3) cout << ", " << Coord[2]*12.0;
              cout <<   ")." << endl;
            }
            cout <<"Maximum Omega " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetOmega_Max() << ", maximum Strain Rate " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetStrainMag_Max() << "." << endl;

            /*--- Print out the number of non-physical points and reconstructions ---*/
            if (config[val_iZone]->GetNonphysical_Points() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
            if (config[val_iZone]->GetNonphysical_Reconstr() > 0)
              cout << "There are " << config[val_iZone]->GetNonphysical_Reconstr() << " non-physical states in the upwind reconstruction." << endl;

            cout << "-------------------------------------------------------------------------" << endl;

            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            if (incompressible) cout << "   Res[Press]";
            else cout << "      Res[Rho]";//, cout << "     Res[RhoE]";

            switch (config[val_iZone]->GetKind_Turb_Model()) {
              case SA: case SA_NEG: case SA_E: case SA_E_COMP: case SA_COMP:        cout << "       Res[nu]"; break;
              case SST:       cout << "     Res[kine]" << "     Res[omega]"; break;
            }

            if (weakly_coupled_heat) {
              cout <<  "     Res[Heat]";
            }

            if (transition) { cout << "      Res[Int]" << "       Res[Re]"; }
            else if (rotating_frame && nDim == 3 && !turbo ) cout << "   CThrust(Total)" << "   CTorque(Total)";
            else if (aeroelastic) cout << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch";
            else if (equiv_area) cout << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)";
            else if (turbo){
              if (nZone < 2){
                /*--- single zone output ---*/
                cout << "  TotPresLoss(%)" << "  Entropy Gen.(%)";
              }
              else{
                /*--- multi zone output ---*/
                cout << " TTEfficiency(%)" << " Entropy Gen.(%)";

              }
            }
            else if (weakly_coupled_heat) {
              cout << "   HFlux(Total)";
            }
            else cout << "   CLift(Total)"   << "   CDrag(Total)";

            if(extra_heat_output) {
              cout <<  "      Res[Heat]" << "   HFlux(Total)";
            }

            cout << endl;

            break;

            case HEAT_EQUATION :
              if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
              else cout << endl << " IntIter" << "  ExtIter";

              cout << "      Res[Heat]" << "   CHeat(Total)"<<  endl;
              break;

            case HEAT_EQUATION_FVM :
              if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
              else cout << endl << " IntIter" << "  ExtIter";

              cout <<  "      Res[Heat]" << "   HFlux(Total)";
              break;

            case ADJ_EULER :              case ADJ_NAVIER_STOKES :
            case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:

              /*--- Visualize the maximum residual ---*/
              iPointMaxResid = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0);
              Coord = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetPoint_Max_Coord(0);
              cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)) << "." << endl;
              if (config[val_iZone]->GetSystemMeasurements() == SI) {
                cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
                if (nDim == 3) cout << ", " << Coord[2];
                cout <<   ")." << endl;
              }
              else {
                cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0]*12.0 << ", " << Coord[1]*12.0;
                if (nDim == 3) cout << ", " << Coord[2]*12.0;
                cout <<   ")." << endl;
              }

              /*--- Print out the number of non-physical points and reconstructions ---*/
              if (config[val_iZone]->GetNonphysical_Points() > 0)
                cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;

              if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
              else cout << endl << " IntIter" << "  ExtIter";

              if (incompressible) {
                if (energy) {cout << "   Res[Psi_Press]" << "   Res[Psi_Temp]";}
                else {cout << "   Res[Psi_Press]" << "   Res[Psi_Velx]";}
              }
              else cout << "   Res[Psi_Rho]" << "     Res[Psi_E]";
              if (disc_adj) {
                if (!turbo){
                  if (compressible) {
                    cout << "    Sens_Press" << "      Sens_AoA" << endl;
                  }
                  if (incompressible) {
                    if (energy) {
                      cout << "      Sens_Vin" << "     Sens_Temp" << endl;
                    } else {
                      cout << "      Sens_Vin" << "     Sens_Pout" << endl;
                    }
                  }                } else {
                  cout << " Sens_PressOut" << " Sens_TotTempIn" << endl;
                }
              } else {
                cout << "      Sens_Geo" << "      Sens_AoA" << endl;
              }
              break;

            case ADJ_RANS : case DISC_ADJ_RANS:

              /*--- Visualize the maximum residual ---*/
              iPointMaxResid = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0);
              Coord = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetPoint_Max_Coord(0);
              cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)) << "." << endl;
              if (config[val_iZone]->GetSystemMeasurements() == SI) {
                cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
                if (nDim == 3) cout << ", " << Coord[2];
                cout <<   ")." << endl;
              }
              else {
                cout <<"Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0]*12.0 << ", " << Coord[1]*12.0;
                if (nDim == 3) cout << ", " << Coord[2]*12.0;
                cout <<   ")." << endl;
              }

              /*--- Print out the number of non-physical points and reconstructions ---*/
              if (config[val_iZone]->GetNonphysical_Points() > 0)
                cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;

              if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
              else cout << endl << " IntIter" << "  ExtIter";

              if (incompressible) cout << "     Res[Psi_Press]";
              else cout << "     Res[Psi_Rho]";

              if (!frozen_visc) {
                cout << "      Res[Psi_Turb[0]]";
              }
              else {
                if (incompressible) {if (energy) {cout << "   Res[Psi_Temp]";}
                  else {cout << "   Res[Psi_Velx]";}}
                else cout << "     Res[Psi_E]";
              }
              if (disc_adj) {
                if (!turbo){
                  if (compressible) {
                  cout << "    Sens_Press" << "      Sens_AoA" << endl;
                  }
                  if (incompressible) {
                    cout << "      Sens_Vin" << "     Sens_Pout" << endl;
                  }
                } else {
                  cout << " Sens_PressOut" << " Sens_TotTempIn" << endl;                }
              } else {
                cout << "      Sens_Geo" << "      Sens_AoA" << endl;
              }
              break;

            case DISC_ADJ_FEM :
              cout << endl << " IntIter" << " ExtIter";

              if (nDim == 2) cout << "    Res[Ux_adj]" << "    Res[Uy_adj]" << "       Sens[E]" << "      Sens[Nu]"<<  endl;
              if (nDim == 3) cout << "    Res[Ux_adj]" << "    Res[Uy_adj]" << "    Res[Uz_adj]" << "       Sens[E]" << "      Sens[Nu]"<<  endl;

           break;

          }

        }
      }

      /*--- Write the solution on the screen ---*/

      if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
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

      switch (config[val_iZone]->GetKind_Solver()) {
        case EULER : case NAVIER_STOKES:

          /*--- Write history file ---*/

          if ((!DualTime_Iteration) && (output_files)) {
            if (!turbo) {
              HistFile << begin << direct_coeff;
              if (thermal) HistFile << heat_coeff;
              if (equiv_area) HistFile << equivalent_area_coeff;
              if (engine || actuator_disk) HistFile << engine_coeff;
              if (inv_design) {
                HistFile << Cp_inverse_design;
                if (thermal) HistFile << Heat_inverse_design;
              }
              if (rotating_frame && !turbo) HistFile << rotating_frame_coeff;
              HistFile << flow_resid;
              if (weakly_coupled_heat) HistFile << heat_resid;
            }
            else {
              HistFile << begin << turbo_coeff << flow_resid;
            }

            if (aeroelastic) HistFile << aeroelastic_coeff;
            if (output_per_surface) HistFile << monitoring_coeff;
            if (output_surface) HistFile << surface_outputs;
            if (direct_diff != NO_DERIVATIVE) HistFile << d_direct_coeff;
            if (output_comboObj) HistFile << combo_obj;
            HistFile << end;
            HistFile.flush();
          }

          /*--- Write screen output ---*/
          if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
            if(DualTime_Iteration || !Unsteady) {
              cout.precision(6);
              cout.setf(ios::fixed, ios::floatfield);
              cout.width(13); cout << log10(residual_flow[0]);
              if (!equiv_area) {
                if (compressible) {
                  if (nDim == 2 ) { cout.width(14); cout << log10(residual_flow[3]); }
                  else { cout.width(14); cout << log10(residual_flow[4]); }
                }
                if (incompressible && !weakly_coupled_heat) {
                  if (energy) {cout.width(14); cout << log10(residual_flow[nDim+1]);}
                  else {cout.width(14); cout << log10(residual_flow[1]);}
                }
                if (incompressible && weakly_coupled_heat)  { cout.width(14); cout << log10(residual_heat[0]);}

              }

              if (rotating_frame && nDim == 3 && !turbo ) {
                cout.setf(ios::scientific, ios::floatfield);
                cout.width(15); cout << Total_CT;
                cout.width(15); cout << Total_CQ;
                cout.unsetf(ios_base::floatfield);
              }
              else if (equiv_area) { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CL)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CD)); cout.width(15);
              cout.precision(4);
              cout.setf(ios::scientific, ios::floatfield);
              cout << Total_CNearFieldOF; }
              else if (turbo) {
                cout.setf(ios::scientific, ios::floatfield);

                if (nZone < 2) {
                  cout.width(15); cout << TotalPressureLoss[0][nSpanWiseSections]*100.0;
                  cout.width(15); cout << EntropyGen[0][nSpanWiseSections]*100.0;
                }
                else {
                  cout.width(15); cout << TotalTotalEfficiency[nTurboPerf -1][nSpanWiseSections]*100.0;
                  cout.width(15); cout << EntropyGen[nTurboPerf -1][nSpanWiseSections]*100.0;
                }

                cout.unsetf(ios_base::floatfield);

              }
              else if (weakly_coupled_heat) { cout.width(14); cout << log10(Total_Heat); }
              else { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CL)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CD)); }
              if (aeroelastic) {
                cout.setf(ios::scientific, ios::floatfield);
                cout.width(15); cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
                cout.width(15); cout << aeroelastic_pitch[0];
                cout.unsetf(ios_base::floatfield);
              }

              if (extra_heat_output) { cout.width(15); cout << Extra_Heat_Residual; cout.width(15); cout << Extra_Total_Heat; }
            }
            cout << endl;
          }
          break;

        case RANS :

          /*--- Write history file ---*/

          if ((!DualTime_Iteration) && (output_files)) {

            if (!turbo) {
              HistFile << begin << direct_coeff;
              if (thermal) HistFile << heat_coeff;
              if (equiv_area) HistFile << equivalent_area_coeff;
              if (engine || actuator_disk) HistFile << engine_coeff;
              if (inv_design) {
                HistFile << Cp_inverse_design;
                if (thermal) HistFile << Heat_inverse_design;
              }
              if (rotating_frame && !turbo) HistFile << rotating_frame_coeff;
              HistFile << flow_resid << turb_resid;
              if (weakly_coupled_heat) HistFile << heat_resid;
            }
            else {
              HistFile << begin << turbo_coeff << flow_resid << turb_resid;
            }

            if (aeroelastic) HistFile << aeroelastic_coeff;
            if (output_per_surface) HistFile << monitoring_coeff;
            if (output_surface) HistFile << surface_outputs;
            if (direct_diff != NO_DERIVATIVE) HistFile << d_direct_coeff;
            if (output_comboObj) HistFile << combo_obj;
            HistFile << end;
            HistFile.flush();
          }

          /*--- Write screen output ---*/

          if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
            if(DualTime_Iteration || !Unsteady) {
              cout.precision(6);
              cout.setf(ios::fixed, ios::floatfield);

              if (incompressible) cout.width(13);
              else  cout.width(14);
              cout << log10(residual_flow[0]);
              switch(nVar_Turb) {
              case 1: cout.width(14); cout << log10(residual_turbulent[0]); break;
              case 2: cout.width(14); cout << log10(residual_turbulent[0]);
              cout.width(15); cout << log10(residual_turbulent[1]); break;
              }

              if (weakly_coupled_heat) {
                cout.width(14); cout << log10(residual_heat[0]);
              }

              if (transition) { cout.width(14); cout << log10(residual_transition[0]); cout.width(14); cout << log10(residual_transition[1]); }

              if (rotating_frame && nDim == 3 && !turbo ) {
                cout.setf(ios::scientific, ios::floatfield);
                cout.width(15); cout << Total_CT; cout.width(15);
                cout << Total_CQ;
                cout.unsetf(ios_base::floatfield);
              }
              else if (equiv_area) { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CL)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CD)); cout.width(15);
              cout.precision(4);
              cout.setf(ios::scientific, ios::floatfield);
              cout << Total_CNearFieldOF; }
              else if (turbo) {
                cout.setf(ios::scientific, ios::floatfield);
                if (nZone < 2){
                  /*--- single zone output ---*/
                  cout.width(15); cout << TotalPressureLoss[0][nSpanWiseSections]*100.0;
                  cout.width(15); cout << EntropyGen[0][nSpanWiseSections]*100.0;
                }
                else{
                  /*--- multi zone output ---*/
                  cout.width(15); cout << TotalTotalEfficiency[nTurboPerf - 1][nSpanWiseSections]*100.0;
                  cout.width(15); cout << EntropyGen[nTurboPerf -1][nSpanWiseSections]*100.0;
                  if (direct_diff){
                    cout.width(15); cout << D_EntropyGen;
                  }
                }
                cout.unsetf(ios_base::floatfield);
              }
              else if (weakly_coupled_heat) { cout.width(15); cout << Total_Heat; }
              else { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CL)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CD)); }

              if (aeroelastic) {
                cout.setf(ios::scientific, ios::floatfield);
                cout.width(15); cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
                cout.width(15); cout << aeroelastic_pitch[0];
                cout.unsetf(ios_base::floatfield);
              }

              if (extra_heat_output) { cout.width(15); cout << Extra_Heat_Residual; cout.width(15); cout << Extra_Total_Heat; }
              cout << endl;
            }
          }
          break;

        case HEAT_EQUATION:

          if (!DualTime_Iteration) {
            HistFile << begin << heat_coeff << heat_resid << end;
            HistFile.flush();
          }
          if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
            cout.precision(6);
            cout.setf(ios::fixed, ios::floatfield);
            cout.width(14); cout << log10(residual_heat[0]);
            cout.width(14); cout << Total_CHeat;
            cout << endl;
          }
          break;

        case HEAT_EQUATION_FVM:

          if (!DualTime_Iteration) {
            HistFile << begin << direct_coeff << heat_resid << end;
            HistFile.flush();
          }
          break;

        case ADJ_EULER :              case ADJ_NAVIER_STOKES :
        case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:

          if (!DualTime_Iteration) {
            HistFile << begin << adjoint_coeff << adj_flow_resid << end;
            HistFile.flush();
          }
          if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
            if (DualTime_Iteration || !Unsteady){
              cout.precision(6);
              cout.setf(ios::fixed, ios::floatfield);
              if (compressible) {
                cout.width(15); cout << log10(residual_adjflow[0]);
                cout.width(15); cout << log10(residual_adjflow[nDim+1]);
              }
              if (incompressible) {
                cout.width(17); cout << log10(residual_adjflow[0]);
                if (energy) {cout.width(16); cout << log10(residual_adjflow[nDim+1]);}
                else {cout.width(16); cout << log10(residual_adjflow[1]);}
              }

              if (disc_adj) {
                cout.precision(4);
                cout.setf(ios::scientific, ios::floatfield);
                if (!turbo){
                  if (compressible) {
                  cout.width(14); cout << Total_Sens_Press;
                  cout.width(14); cout << Total_Sens_AoA;
                  }
                  if (incompressible) {
                    cout.width(14); cout << Total_Sens_ModVel;
                    if (energy) {
                      cout.width(14); cout << Total_Sens_Temp;
                    } else {
                      cout.width(14); cout << Total_Sens_BPressure;
                    }
                  }
                } else {
                  cout.width(14); cout << Total_Sens_BPressure;
                  cout.width(15); cout << Total_Sens_Temp;
                }
              }else {
                cout.precision(4);
                cout.setf(ios::scientific, ios::floatfield);
                cout.width(14); cout << Total_Sens_Geo;
                cout.width(14); cout << Total_Sens_AoA;
              }
              cout << endl;
              cout.unsetf(ios_base::floatfield);
            }
          }
          break;

        case ADJ_RANS : case DISC_ADJ_RANS:

          if (!DualTime_Iteration) {
            HistFile << begin << adjoint_coeff << adj_flow_resid;
            if (!frozen_visc)
              HistFile << adj_turb_resid;
            HistFile << end;
            HistFile.flush();
          }
          if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
            if (DualTime_Iteration || !Unsteady){
              cout.precision(6);
              cout.setf(ios::fixed, ios::floatfield);
              cout.width(17); cout << log10(residual_adjflow[0]);
            if (!frozen_visc) {
                cout.width(17); cout << log10(residual_adjturbulent[0]);
              }
              else {
                if (compressible) {
                  if (geometry[val_iZone][val_iInst][FinestMesh]->GetnDim() == 2 ) { cout.width(15); cout << log10(residual_adjflow[3]); }
                  else { cout.width(15); cout << log10(residual_adjflow[4]); }
                }
                if (incompressible) {
                  if (energy) {cout.width(15); cout << log10(residual_adjflow[nDim+1]);}
                  else {cout.width(15); cout << log10(residual_adjflow[1]);}
                }
              }
              if (disc_adj) {
                if (!turbo){
                  if (compressible) {
                  cout.width(14); cout << Total_Sens_Press;
                  cout.width(14); cout << Total_Sens_AoA;
                  }
                  if (incompressible) {
                    cout.width(14); cout << Total_Sens_ModVel;
                    if (energy) {
                      cout.width(14); cout << Total_Sens_Temp;
                    } else {
                      cout.width(14); cout << Total_Sens_BPressure;
                    }                  }
                } else {
                  cout.width(14); cout << Total_Sens_BPressure;
                  cout.width(15); cout << Total_Sens_Temp;
                }
              }else {
                cout.precision(4);
                cout.setf(ios::scientific, ios::floatfield);
                cout.width(14); cout << Total_Sens_Geo;
                cout.width(14); cout << Total_Sens_AoA;
              }
              cout << endl;
              cout.unsetf(ios_base::floatfield);
            }
          }
          break;

      }
      cout.unsetf(ios::fixed);

    }


    delete [] residual_flow;
    delete [] residual_turbulent;
    delete [] residual_transition;
    delete [] residual_heat;

    delete [] residual_adjflow;
    delete [] residual_adjturbulent;

    delete [] Surface_CL;
    delete [] Surface_CD;
    delete [] Surface_CSF;
    delete [] Surface_CEff;
    delete [] Surface_CFx;
    delete [] Surface_CFy;
    delete [] Surface_CFz;
    delete [] Surface_CMx;
    delete [] Surface_CMy;
    delete [] Surface_CMz;
    delete [] aeroelastic_pitch;
    delete [] aeroelastic_plunge;

  }
}

